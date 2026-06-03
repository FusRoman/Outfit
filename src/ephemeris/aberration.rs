//! Stellar aberration corrections for apparent-position computation.
//!
//! This module provides two levels of aberration correction:
//!
//! | Function | Description |
//! |---|---|
//! | [`correct_aberration_first_order`]  | Linear light-travel-time shift |
//! | [`correct_aberration_second_order`] | Two-step Keplerian back-propagation |
//!
//! The correction to apply is selected at call-site via
//! [`AberrationOrder`](super::AberrationOrder).
//!
//! # Physical background
//!
//! When an observer detects a photon, the emitting body has already moved on.
//! The apparent direction corresponds to the body's position at the **retarded
//! epoch** — the moment the photon was emitted, not the moment it was received.
//!
//! Both corrections account for this light-travel delay; they differ in how
//! accurately they approximate the retarded position:
//!
//! - **First order** computes the delay from the instantaneous separation and
//!   subtracts a *linear* displacement along the current velocity.  Accurate to
//!   \\( O(v/c) \\).
//!
//! - **Second order** iterates the delay twice and back-propagates the orbit
//!   along the **Keplerian two-body solution** at each step, capturing the
//!   orbital curvature during the light-travel time.  Necessary for sub-mas
//!   accuracy on close-approach objects or highly curved orbits.
//!
//! # Coordinate conventions
//!
//! All vectors are in the **equatorial mean J2000** frame, positions in AU,
//! velocities in AU/day.

use nalgebra::Vector3;

use crate::{constants::ROT_ECLMJ2000_TO_EQUMJ2000, EquinoctialElements, OutfitError, VLIGHT_AU};

// ---------------------------------------------------------------------------
// Aberration order
// ---------------------------------------------------------------------------

/// Stellar aberration correction order applied during apparent-position
/// computation.
///
/// Controls whether the pipeline uses the fast linear approximation or the
/// more accurate two-step Keplerian back-propagation.
///
/// For the vast majority of targets (main-belt asteroids, typical NEOs) the
/// difference between the two corrections is sub-milliarcsecond and
/// [`AberrationOrder::First`] is sufficient.  [`AberrationOrder::Second`]
/// becomes relevant for close-approach objects (geocentric distance
/// \\( \lesssim 0.01 \\) AU) or highly curved orbits where the linear
/// approximation breaks down.
///
/// See `correct_aberration_first_order` and `correct_aberration_second_order`
/// for the implementation of both corrections.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub enum AberrationOrder {
    /// First-order correction: linear light-travel-time shift.
    ///
    /// $$\mathbf{d}_\text{corr} = \mathbf{d} - \frac{|\mathbf{d}|}{c}\,\mathbf{v}_\text{body}$$
    ///
    /// Fast and self-contained — requires only the instantaneous body velocity.
    #[default]
    First,

    /// Second-order correction: two-step Keplerian back-propagation.
    ///
    /// Iterates the light-travel delay twice and back-propagates the orbit via
    /// the analytic two-body solution at each step.  More accurate for objects
    /// at short geocentric distance or with high orbital curvature.
    Second,
}

impl AberrationOrder {
    /// Dispatch aberration correction to the first-order or second-order function.
    ///
    /// Exists as a dedicated helper to keep `apparent_position_from_state` free of
    /// match boilerplate and to give the dispatch a named entry point for testing.
    ///
    /// # Errors
    ///
    /// Forwards any error from [`correct_aberration_second_order`].
    pub(crate) fn apply_aberration(
        &self,
        topocentric_vec: Vector3<f64>,
        ast_vel_equ: Vector3<f64>,
        elements: &EquinoctialElements,
        obs_time_mjd: f64,
        obs_pos_equ: Vector3<f64>,
    ) -> Result<Vector3<f64>, OutfitError> {
        match self {
            AberrationOrder::First => {
                Ok(correct_aberration_first_order(topocentric_vec, ast_vel_equ))
            }
            AberrationOrder::Second => correct_aberration_second_order(
                topocentric_vec,
                elements,
                obs_time_mjd,
                obs_pos_equ,
            ),
        }
    }
}

// ---------------------------------------------------------------------------
// First-order correction
// ---------------------------------------------------------------------------

/// Apply the **first-order** stellar aberration correction to a topocentric
/// position vector.
///
/// The apparent direction is obtained by subtracting the linear displacement
/// of the body during the light-travel time \\( \Delta t = |\mathbf{d}| / c \\):
///
/// $$\mathbf{d}_\text{corr} = \mathbf{d} - \frac{|\mathbf{d}|}{c}\,\mathbf{v}_\text{body}$$
///
/// where \\( c \\) is the speed of light in AU/day ([`VLIGHT_AU`]).
///
/// This approximation is valid for the vast majority of solar system targets
/// (main-belt asteroids, typical NEOs).  The residual error with respect to
/// the exact retarded position is of order \\( O((v/c)^2) \\), sub-mas for
/// most objects.
///
/// # Arguments
///
/// - `topocentric_vec` – Vector from observer to body \\( \mathbf{d} \\) \[AU\],
///   equatorial mean J2000.
/// - `body_velocity`   – Body heliocentric velocity \\( \mathbf{v} \\) \[AU/day\],
///   equatorial mean J2000.
///
/// # Returns
///
/// Aberration-corrected topocentric direction vector \[AU\].  The magnitude is
/// slightly different from the input; only the direction is used downstream.
#[inline]
pub(crate) fn correct_aberration_first_order(
    topocentric_vec: Vector3<f64>,
    body_velocity: Vector3<f64>,
) -> Vector3<f64> {
    let dt = topocentric_vec.norm() / VLIGHT_AU;
    topocentric_vec - dt * body_velocity
}

// ---------------------------------------------------------------------------
// Second-order correction
// ---------------------------------------------------------------------------

/// Apply the **second-order** stellar aberration correction via two-step
/// Keplerian back-propagation.
///
/// Rather than shifting the current position linearly, this function
/// propagates the orbit *backwards* by the estimated light-travel time, twice
/// in succession, to recover the retarded position with sub-mas accuracy.
///
/// ## Algorithm
///
/// Let \\( \mathbf{d}_0 = \mathbf{r}_\text{body} - \mathbf{r}_\text{obs} \\)
/// be the instantaneous topocentric vector.
///
/// **Pass 1:**
/// $$\Delta t_0 = |\mathbf{d}_0| / c$$
/// $$\mathbf{r}_1 = \text{propagate\_twobody}(t_\text{obs} - \Delta t_0)$$
/// $$\mathbf{d}_1 = \mathbf{r}_1 - \mathbf{r}_\text{obs}$$
///
/// **Pass 2:**
/// $$\Delta t_1 = |\mathbf{d}_1| / c$$
/// $$\mathbf{r}_2 = \text{propagate\_twobody}(t_\text{obs} - \Delta t_1)$$
///
/// **Result:**
/// $$\mathbf{d}_\text{corr} = \mathbf{r}_2 - \mathbf{r}_\text{obs}$$
///
/// The two-body propagator is used for both passes regardless of the main
/// propagator choice in [`EphemerisConfig`](super::EphemerisConfig).
///
/// # Arguments
///
/// - `topocentric_vec` – Instantaneous topocentric vector \\( \mathbf{d}_0 \\)
///   \[AU\], equatorial mean J2000.
/// - `elements`        – Equinoctial orbital elements at their reference epoch.
/// - `obs_time_mjd`    – Observation epoch \[MJD TT\].
/// - `obs_pos_equ`     – Observer heliocentric position \[AU\], equatorial
///   mean J2000.
///
/// # Returns
///
/// Aberration-corrected topocentric direction vector \[AU\].
///
/// # Errors
///
/// Returns [`OutfitError`] if either two-body propagation step fails to
/// converge (e.g. degenerate orbit).
pub(crate) fn correct_aberration_second_order(
    topocentric_vec: Vector3<f64>,
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    obs_pos_equ: Vector3<f64>,
) -> Result<Vector3<f64>, OutfitError> {
    let r1 = retropropagate(elements, obs_time_mjd, topocentric_vec.norm())?;
    let d1 = r1 - obs_pos_equ;

    let r2 = retropropagate(elements, obs_time_mjd, d1.norm())?;

    Ok(r2 - obs_pos_equ)
}

// ---------------------------------------------------------------------------
// Private helper
// ---------------------------------------------------------------------------

/// Back-propagate the orbit by a light-travel time derived from `separation`
/// and return the body's heliocentric position at the retarded epoch, in the
/// **equatorial mean J2000** frame \[AU\].
///
/// Computes the retarded epoch as
/// \\( t_\text{ret} = t_\text{obs} - |\text{separation}| / c \\)
/// and calls [`EquinoctialElements::propagate_twobody`].
///
/// # Errors
///
/// Returns [`OutfitError`] if the Kepler solver does not converge.
fn retropropagate(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    separation: f64,
) -> Result<Vector3<f64>, OutfitError> {
    let dt_light = separation / VLIGHT_AU;
    let t_retarded = obs_time_mjd - dt_light;
    let dt_orbit = t_retarded - elements.reference_epoch;
    let (pos_ecl, _, _) = elements.propagate_twobody(0.0, dt_orbit, false)?;
    Ok(ROT_ECLMJ2000_TO_EQUMJ2000 * pos_ecl)
}
