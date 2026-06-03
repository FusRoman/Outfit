//! Orbit propagation strategies for differential orbit correction.
//!
//! This module exposes two propagator kinds:
//!
//! - `PropagatorKind::TwoBody` – classic Keplerian propagation through the
//!   analytic two-body solution.  This is the default and requires no external
//!   ephemeris data beyond the standard solar GM.
//!
//! - `PropagatorKind::NBody` – numerical N-body propagation via an 8th-order
//!   Dormand-Prince (DOP853) integrator.  The target object is integrated
//!   together with its 6×6 state transition matrix (STM) under the influence of
//!   the Sun and any additional perturbing bodies selected in `NBodyConfig`.
//!   Planetary positions are looked up from the `JPLEphem` file supplied to the
//!   differential corrector at runtime.

pub mod nbody;
pub mod planet_gm;

use nalgebra::Vector3;

use crate::{
    constants::ROT_ECLMJ2000_TO_EQUMJ2000, jpl_ephem::naif::naif_ids::NaifIds, EquinoctialElements,
    JPLEphem, OutfitError,
};

/// Select which propagator is used during differential orbit correction.
#[derive(Debug, Clone, Default)]
pub enum PropagatorKind {
    /// Keplerian (analytic) two-body propagation.  Fast and self-contained.
    #[default]
    TwoBody,

    /// Numerical N-body propagation with user-specified perturbing bodies.
    NBody(NBodyConfig),
}

impl PropagatorKind {
    /// Propagate `elements` to `obs_time_mjd` and return the heliocentric
    /// position and velocity in the **equatorial mean J2000** frame.
    ///
    /// Dispatches to the analytic two-body solver ([`PropagatorKind::TwoBody`])
    /// or to the numerical DOP853 integrator ([`PropagatorKind::NBody`])
    /// based on the variant of `self`.
    ///
    /// # Returns
    ///
    /// `(position [AU], velocity [AU/day])`, both in equatorial mean J2000.
    ///
    /// # Errors
    ///
    /// Propagates any error returned by the underlying propagator.
    pub(crate) fn propagate_to_epoch(
        &self,
        elements: &EquinoctialElements,
        obs_time_mjd: f64,
        jpl: &JPLEphem,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        match self {
            PropagatorKind::TwoBody => propagate_twobody(elements, obs_time_mjd),
            PropagatorKind::NBody(nbody_config) => {
                propagate_nbody(elements, obs_time_mjd, jpl, nbody_config)
            }
        }
    }
}

/// Propagate the orbit analytically using the Keplerian two-body model.
///
/// Computes the time-of-flight \\( \Delta t = t_\text{obs} - t_\text{ref} \\)
/// \[days\], calls [`EquinoctialElements::propagate_twobody`], and rotates the
/// result from ecliptic to equatorial mean J2000 via [`ecl_state_to_equ`].
///
/// # Returns
///
/// `(position [AU], velocity [AU/day])` in equatorial mean J2000.
///
/// # Errors
///
/// Returns [`OutfitError`] if the Kepler solver does not converge.
fn propagate_twobody(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
    let dt = obs_time_mjd - elements.reference_epoch;
    let (pos_ecl, vel_ecl, _) = elements.propagate_twobody(0.0, dt, false)?;
    Ok(ecl_state_to_equ(pos_ecl, vel_ecl))
}

/// Propagate the orbit numerically using the N-body DOP853 integrator.
///
/// Calls [`EquinoctialElements::propagate_nbody`] with the supplied
/// [`NBodyConfig`] (perturbing bodies, integration tolerances) and rotates
/// the result from ecliptic to equatorial mean J2000 via [`ecl_state_to_equ`].
///
/// # Returns
///
/// `(position [AU], velocity [AU/day])` in equatorial mean J2000.
///
/// # Errors
///
/// Returns [`OutfitError`] if the integrator fails or the JPL data is
/// unavailable for the integration span.
fn propagate_nbody(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    jpl: &JPLEphem,
    config: &NBodyConfig,
) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
    let result = elements.propagate_nbody(obs_time_mjd, jpl, config)?;
    Ok(ecl_state_to_equ(result.position, result.velocity))
}

/// Rotate an ecliptic-frame position+velocity pair into **equatorial mean
/// J2000** using the pre-computed rotation matrix
/// [`ROT_ECLMJ2000_TO_EQUMJ2000`](crate::constants::ROT_ECLMJ2000_TO_EQUMJ2000).
///
/// The rotation corresponds to a positive rotation of \\( +\varepsilon \\)
/// around the X-axis, where \\( \varepsilon \\) is the obliquity of the
/// ecliptic at J2000.
///
/// # Returns
///
/// `(pos_equ [AU], vel_equ [AU/day])` in equatorial mean J2000.
#[inline]
fn ecl_state_to_equ(pos_ecl: Vector3<f64>, vel_ecl: Vector3<f64>) -> (Vector3<f64>, Vector3<f64>) {
    (
        ROT_ECLMJ2000_TO_EQUMJ2000 * pos_ecl,
        ROT_ECLMJ2000_TO_EQUMJ2000 * vel_ecl,
    )
}

/// Configuration for the N-body propagator.
#[derive(Debug, Clone)]
pub struct NBodyConfig {
    /// Bodies whose gravitational attraction perturbs the target orbit.
    ///
    /// Defaults to `[Sun]`.  For each body the GM is taken from [`planet_gm`]
    /// and its ephemeris is obtained from the `JPLEphem` file passed to the
    /// integrator at runtime.
    pub perturbing_bodies: Vec<NaifIds>,

    /// Absolute tolerance for the DOP853 step-size control.
    ///
    /// Units: AU.  Defaults to `1e-12`.
    pub abs_tol: f64,

    /// Relative tolerance for the DOP853 step-size control.
    ///
    /// Defaults to `1e-12`.
    pub rel_tol: f64,
}

impl Default for NBodyConfig {
    fn default() -> Self {
        use crate::jpl_ephem::naif::naif_ids::solar_system_bary::SolarSystemBary;
        Self {
            perturbing_bodies: vec![NaifIds::SSB(SolarSystemBary::Sun)],
            abs_tol: 1e-12,
            rel_tol: 1e-12,
        }
    }
}
