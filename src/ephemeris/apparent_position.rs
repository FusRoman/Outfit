//! Standalone apparent-position computation for a solar system body.
//!
//! This module provides [`ApparentPosition`] as the public return type and a set
//! of small, composable helper functions that together implement the pipeline:
//!
//! ```text
//! OrbitalElements → propagate → topocentric geometry → aberration → (RA, Dec)
//! ```
//!
//! # Coordinate conventions
//!
//! - Positions in **AU**, velocities in **AU/day**.
//! - Intermediate frames: **ecliptic mean J2000**.
//! - Final output: **equatorial mean J2000** (RA, Dec in radians).
//! - Time: **MJD TT**.
//!
//! # Pipeline steps
//!
//! | Step | Function / Method | Purpose |
//! |------|-------------------|---------|
//! | 1 | [`obs_time_to_epoch`] | Convert MJD TT scalar → [`hifitime::Epoch`] |
//! | 2 | [`PropagatorKind::propagate_to_epoch`](crate::propagator::PropagatorKind::propagate_to_epoch) | Propagate orbit; rotate to equatorial J2000 |
//! | 3 | [`observer_pv`] | Resolve observer heliocentric position **and velocity** |
//! | 4 | [`assemble_apparent_position`] | Compute topocentric vector, apply aberration, convert to (RA, Dec) |
//!
//! # Aberration model
//!
//! The first-order stellar aberration correction shifts the topocentric
//! line-of-sight vector by
//!
//! $$\mathbf{x}\_\text{corr} = \mathbf{x}\_\text{topo}
//!     - \frac{|\mathbf{x}\_\text{topo}|}{c}\\,\mathbf{v}\_\text{body}$$
//!
//! where $c$ is the speed of light in AU/day and
//! $\mathbf{v}\_\text{body}$ is the body's heliocentric velocity.

use hifitime::{ut1::Ut1Provider, Epoch, TimeScale};
use nalgebra::Vector3;
use photom::{
    coordinates::{cartesian::CartesianCoord, equatorial::EquCoord},
    observer::Observer,
};

use crate::{
    cache::observer_fixed_cache::ObserverFixedCache, constants::ROT_ECLMJ2000_TO_EQUMJ2000,
    conversion::ToNotNan, observer_extension::ResolvedObserver, EquinoctialElements, JPLEphem,
    OutfitError,
};

use super::{AberrationOrder, EphemerisConfig};

// ---------------------------------------------------------------------------
// Public return type
// ---------------------------------------------------------------------------

/// Predicted apparent position of a solar system body together with geometric
/// distances.
///
/// The `coord` field contains the predicted equatorial sky position in the
/// **equatorial mean J2000** frame.  Because this is a *prediction* rather than
/// a measurement, the error fields `ra_error` and `dec_error` inside `coord`
/// are always set to `0.0`.
///
/// Distances are computed from the unaberrated heliocentric state before the
/// aberration correction is applied to the line-of-sight direction.
#[derive(Debug, Clone, PartialEq)]
pub struct ApparentPosition {
    /// Predicted equatorial coordinates (RA ∈ \[0, 2π), Dec ∈ (−π/2, π/2),
    /// both in radians; error fields = 0 because this is a prediction).
    pub coord: EquCoord,
    /// Distance from the **Earth's centre** to the body \[AU\].
    ///
    /// Computed as $|\mathbf{r}\_\text{body} - \mathbf{r}\_\text{Earth}|$.
    pub geocentric_dist: f64,
    /// Distance from the **Sun** to the body \[AU\].
    ///
    /// Computed as $|\mathbf{r}\_\text{body}|$ (heliocentric norm).
    pub heliocentric_dist: f64,
}

// ---------------------------------------------------------------------------
// Intermediate propagated state (shared with geometry module)
// ---------------------------------------------------------------------------

/// Full propagated state at a given epoch, shared between position and geometry
/// computations.
///
/// Holding this intermediate result allows [`compute_with_geometry`] to run a
/// single orbit propagation and observer-position query and hand the results to
/// both [`assemble_apparent_position`] and
/// [`super::geometry::compute_geometry`] without redundant work.
pub(crate) struct PropagatedState {
    /// Body heliocentric position \[AU\], equatorial mean J2000.
    pub ast_pos_equ: Vector3<f64>,
    /// Body heliocentric velocity \[AU/day\], equatorial mean J2000.
    pub ast_vel_equ: Vector3<f64>,
    /// Observer heliocentric position \[AU\], equatorial mean J2000.
    pub obs_pos_equ: Vector3<f64>,
    /// Observer heliocentric velocity \[AU/day\], equatorial mean J2000.
    pub obs_vel_equ: Vector3<f64>,
    /// Earth heliocentric position \[AU\], equatorial mean J2000.
    pub earth_pos_equ: Vector3<f64>,
    /// Observation epoch as MJD TT scalar (carried for second-order aberration).
    pub obs_time_mjd: f64,
}

// ---------------------------------------------------------------------------
// Internal helpers — propagation and observer geometry
// ---------------------------------------------------------------------------

/// Propagate the orbit and resolve observer geometry, producing a
/// [`PropagatedState`].
///
/// This is the shared kernel called by [`compute`] and [`compute_with_geometry`].
///
/// # Arguments
///
/// - `elements`     – Equinoctial orbital elements.
/// - `obs_time_mjd` – Observation epoch \[MJD TT\].
/// - `fixed_cache`  – Pre-built body-fixed observer cache (epoch-invariant).
///   Must have been constructed from the same [`Observer`] that owns this
///   request slot.  Building it once per observer slot and reusing it across
///   all epochs avoids redundant trigonometric conversions.
/// - `observer`     – Observing site, used only to attach to the result.
/// - `jpl`          – JPL planetary ephemeris.
/// - `ut1`          – UT1 time-scale provider.
/// - `config`       – Ephemeris configuration (propagator, aberration).
///
/// # Errors
///
/// Returns [`OutfitError`] if:
/// - Orbit propagation fails.
/// - The JPL ephemeris data is unavailable for the requested epoch.
/// - The observer geometry cannot be resolved.
pub(crate) fn propagate(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    fixed_cache: &ObserverFixedCache,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Result<PropagatedState, OutfitError> {
    let epoch = obs_time_to_epoch(obs_time_mjd);
    let (ast_pos_equ, ast_vel_equ) =
        config
            .propagator
            .propagate_to_epoch(elements, obs_time_mjd, jpl)?;
    let (obs_pos_equ, obs_vel_equ, earth_pos_equ) = observer_pv(fixed_cache, jpl, ut1, &epoch)?;

    Ok(PropagatedState {
        ast_pos_equ,
        ast_vel_equ,
        obs_pos_equ,
        obs_vel_equ,
        earth_pos_equ,
        obs_time_mjd,
    })
}

// ---------------------------------------------------------------------------
// Entry points
// ---------------------------------------------------------------------------

/// Compute the apparent equatorial position of a solar system body.
///
/// This is the main entry point called by
/// [`OrbitalElements::apparent_position`](crate::OrbitalElements::apparent_position).
/// It propagates the orbit, resolves observer geometry, applies the aberration
/// correction and converts to equatorial coordinates.
///
/// # Arguments
///
/// - `elements`    – Equinoctial orbital elements.
/// - `obs_time_mjd`– Observation epoch \[MJD TT\].
/// - `fixed_cache` – Pre-built body-fixed observer cache (epoch-invariant).
/// - `jpl`         – JPL planetary ephemeris.
/// - `ut1`         – UT1 time-scale provider.
/// - `config`      – Ephemeris configuration.
///
/// # Errors
///
/// Returns [`OutfitError`] if propagation or observer geometry fails.
pub(crate) fn compute(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    fixed_cache: &ObserverFixedCache,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Result<ApparentPosition, OutfitError> {
    let state = propagate(elements, obs_time_mjd, fixed_cache, jpl, ut1, config)?;
    assemble_apparent_position(&state, elements, &config.aberration)
}

/// Compute both the apparent position and the body geometry in a single
/// propagation pass.
///
/// Called by
/// [`OrbitalElements::apparent_position_and_geometry`](crate::OrbitalElements::apparent_position_and_geometry).
/// The orbit is propagated and observer geometry resolved exactly once; the
/// resulting [`PropagatedState`] is then handed to both
/// [`assemble_apparent_position`] and
/// [`super::geometry::compute_geometry`].
///
/// # Arguments
///
/// - `elements`    – Equinoctial orbital elements.
/// - `obs_time_mjd`– Observation epoch \[MJD TT\].
/// - `fixed_cache` – Pre-built body-fixed observer cache (epoch-invariant).
/// - `jpl`         – JPL planetary ephemeris.
/// - `ut1`         – UT1 time-scale provider.
/// - `config`      – Ephemeris configuration.
///
/// # Errors
///
/// Returns [`OutfitError`] if propagation or either assembly step fails.
pub(crate) fn compute_with_geometry(
    elements: &EquinoctialElements,
    obs_time_mjd: f64,
    fixed_cache: &ObserverFixedCache,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Result<(ApparentPosition, super::geometry::BodyGeometry), OutfitError> {
    let state = propagate(elements, obs_time_mjd, fixed_cache, jpl, ut1, config)?;
    let position = assemble_apparent_position(&state, elements, &config.aberration)?;
    let geometry = super::geometry::compute_geometry(&state, elements, &config.aberration)?;
    Ok((position, geometry))
}

// ---------------------------------------------------------------------------
// Step 1 – epoch conversion
// ---------------------------------------------------------------------------

/// Convert a scalar MJD TT value to a [`hifitime::Epoch`].
///
/// The time scale is fixed to [`TimeScale::TT`] (Terrestrial Time), which is
/// the time argument used throughout the library for orbit propagation and
/// ephemeris lookups.
#[inline]
fn obs_time_to_epoch(obs_time_mjd: f64) -> Epoch {
    Epoch::from_mjd_in_time_scale(obs_time_mjd, TimeScale::TT)
}

// ---------------------------------------------------------------------------
// Step 3 – observer position and velocity
// ---------------------------------------------------------------------------

/// Compute the observer's heliocentric position **and velocity**, the Earth's
/// heliocentric position, all in the equatorial mean J2000 frame.
///
/// # Pre-condition — `fixed_cache` is epoch-invariant
///
/// `fixed_cache` must have been constructed once per observer slot (before the
/// epoch loop) via [`ObserverFixedCache::try_from`].  Passing it in avoids
/// rebuilding the body-fixed geocentric coordinates (sin/cos of longitude, ρ
/// factors) for every epoch.
///
/// # Returns
///
/// `(obs_pos_equ [AU], obs_vel_equ [AU/day], earth_pos_equ [AU])`.
///
/// # Errors
///
/// Returns [`OutfitError`] if the geocentric position computation
/// (`pvobs`) fails or a NaN conversion fails.
type ObserverPv = (Vector3<f64>, Vector3<f64>, Vector3<f64>);

fn observer_pv(
    fixed_cache: &ObserverFixedCache,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    epoch: &Epoch,
) -> Result<ObserverPv, OutfitError> {
    // Geocentric position in ecliptic J2000 (velocity not needed here — site
    // rotation is accounted for via earth_vel from JPL in the velocity path).
    let (geo_pos_ecl, _) = Observer::pvobs(epoch, ut1, fixed_cache, false)?;

    // Single JPL Chebyshev evaluation for Earth's heliocentric state.
    let (earth_pos_equ_raw, earth_vel_opt) = jpl.earth_ephemeris(epoch, true);
    let earth_vel_equ_raw = earth_vel_opt
        .expect("JPL earth_ephemeris with compute_velocity=true must return a velocity");

    // Rotation matrix ecliptic → equatorial J2000 (static const, evaluated once).
    let rot = ROT_ECLMJ2000_TO_EQUMJ2000.to_notnan()?;

    // obs_pos_equ = earth_pos + ROT_ecl→equ * geo_pos_ecl
    // (replaces Observer::helio_position)
    let obs_pos_equ: Vector3<f64> = (earth_pos_equ_raw.to_notnan()? + rot * geo_pos_ecl)
        .map(|x: ordered_float::NotNan<f64>| x.into_inner());

    // earth_pos_equ: used in assemble_apparent_position for geocentric distance.
    // (replaces earth_heliocentric_position)
    let earth_pos_equ = earth_pos_equ_raw;

    // obs_vel_equ = earth_vel + ROT_ecl→equ * geo_vel_ecl
    // Since pvobs was called with compute_velocity=false, geo_vel_ecl = 0,
    // so obs_vel_equ = Earth's heliocentric velocity only.
    // (replaces Observer::helio_velocity with a zero geocentric contribution)
    let obs_vel_equ = earth_vel_equ_raw.to_notnan()?.map(|x| x.into_inner());

    Ok((obs_pos_equ, obs_vel_equ, earth_pos_equ))
}

// ---------------------------------------------------------------------------
// Step 4 – apparent position assembly
// ---------------------------------------------------------------------------

/// Assemble the [`ApparentPosition`] from a [`PropagatedState`].
///
/// The function performs the following sub-steps:
///
/// 1. **Distances** — compute heliocentric and geocentric distances.
/// 2. **Topocentric vector** — $\mathbf{d} = \mathbf{r}\_\text{body} - \mathbf{r}\_\text{obs}$.
/// 3. **Aberration correction** — first-order or second-order, per `aberration`.
/// 4. **Sky coordinates** — convert the corrected direction to (RA, Dec).
///
/// # Errors
///
/// Returns [`OutfitError`] if the second-order aberration back-propagation
/// fails to converge.  The first-order path is always infallible.
pub(crate) fn assemble_apparent_position(
    state: &PropagatedState,
    elements: &EquinoctialElements,
    aberration: &AberrationOrder,
) -> Result<ApparentPosition, OutfitError> {
    let heliocentric_dist = state.ast_pos_equ.norm();
    let geocentric_dist = (state.ast_pos_equ - state.earth_pos_equ).norm();

    let topocentric_vec = state.ast_pos_equ - state.obs_pos_equ;
    let corrected = aberration.apply_aberration(
        topocentric_vec,
        state.ast_vel_equ,
        elements,
        state.obs_time_mjd,
        state.obs_pos_equ,
    )?;
    let coord = cartesian_to_equcoord(corrected);

    Ok(ApparentPosition {
        coord,
        geocentric_dist,
        heliocentric_dist,
    })
}

/// Convert a Cartesian direction vector to an [`EquCoord`] (RA, Dec in radians;
/// error fields set to `0.0`).
///
/// Delegates to [`CartesianCoord`] → [`EquCoord`] conversion, which computes:
/// $$\alpha = \operatorname{atan2}(y, x) \bmod 2\pi, \quad
///   \delta = \operatorname{atan2}\left(z,\\,\sqrt{x^2+y^2}\right)$$
///
/// The magnitude of `v` is irrelevant; only the direction is used.
#[inline]
pub(crate) fn cartesian_to_equcoord(v: Vector3<f64>) -> EquCoord {
    EquCoord::from(CartesianCoord {
        x: v[0],
        y: v[1],
        z: v[2],
    })
}
