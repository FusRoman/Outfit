//! # Observations: ingestion, representation, and sky-projection utilities
//!
//! This module defines the core types and helpers to **ingest**, **store**, and **use**
//! optical astrometric observations for orbit determination workflows.
//!
//! ## What lives here?
//!
//! - [`Observation`](crate::observations::Observation) — a single astrometric measurement (RA/DEC at an epoch) with:
//!   - the observing site identifier (`u16`),
//!   - precomputed **geocentric** and **heliocentric** site positions at the epoch,
//!   - astrometric uncertainties for RA/DEC.
//!
//! - Parsing & I/O:
//!   - `from_80col` (private) and `extract_80col` (private) — read **80-column MPC** formatted files.
//!   - [`ades_reader`](crate::trajectories::ades_reader) — ADES ingestion utilities (XML/CSV).
//!   - `parquet_reader` (private) — internal helpers to read columnar batches.
//!
//! - Batch/transform helpers:
//!   - [`trajectory_file`](crate::trajectories::trajectory_file) — build batches of observations (RA/DEC/time + σ) and convert to [`Observation`](crate::observations::Observation)s.
//!   - [`observations_ext`](crate::observations::observations_ext) — higher-level operations on collections (triplet selection, RMS windows, metrics).
//!   - [`triplets_iod`](crate::observations::triplets_iod) — construction of observation triplets for **Gauss IOD**.
//!
//! ## Units & reference frames
//!
//! - **Angles**: radians  
//! - **Time**: MJD (TT scale)  
//! - **Positions**: AU, **equatorial mean J2000** (J2000/ICRS-aligned)  
//!
//! These conventions are enforced by [`Observation::new`](crate::observations::Observation::new), which computes and stores both
//! the **geocentric** and **heliocentric** site positions at the observation epoch using the
//! [`Outfit`](crate::outfit::Outfit) environment (UT1 provider, JPL ephemerides, site database).
//!
//! ## Typical workflow
//!
//! 1. **Ingest** observations:
//!    - From MPC 80-col: \[`extract_80col`\] → `Vec<Observation>` + object identifier.
//!    - From ADES: via [`ades_reader`](crate::trajectories::ades_reader) into typed batches, then \[`observation_from_batch`\].
//!
//! 2. **Precompute/Access positions** per observation:
//!    - `get_observer_earth_position()` — geocentric site vector at epoch.
//!    - `get_observer_helio_position()` — heliocentric site vector at epoch.
//!
//! 3. **Project to sky** (prediction / fitting):
//!    - [`Observation::compute_apparent_position`](crate::observations::Observation::compute_apparent_position) — propagate an orbit (equinoctial elements),
//!      apply frame transforms + aberration, and return apparent `(RA, DEC)`.
//!    - [`Observation::ephemeris_error`](crate::observations::Observation::ephemeris_error) — normalized squared residual for a single observation.
//!
//! 4. **Build triplets and run IOD**:
//!    - Use [`observations_ext`](crate::observations::observations_ext) / [`triplets_iod`](crate::observations::triplets_iod) to form high-quality triplets and feed
//!      them to the Gauss solver (see `initial_orbit_determination::gauss`).
//!
//! ## Key types & functions
//!
//! - [`Observation`](crate::observations::Observation) — single measurement with site & precomputed positions.
//! - [`Observation::compute_apparent_position`](crate::observations::Observation::compute_apparent_position) — apparent `(RA, DEC)` from an orbit.
//! - [`Observation::ephemeris_error`](crate::observations::Observation::ephemeris_error) — per-observation χ²-like contribution.
//!
//! ## Example
//!
//! ```rust,no_run
//! use outfit::observations::Observation;
//! use outfit::outfit::Outfit;
//! use outfit::error_models::ErrorModel;
//! use outfit::constants::RADSEC;
//!
//! # use outfit::orbit_type::equinoctial_element::EquinoctialElements;
//!
//! let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//!
//! // Example: build one Observation manually
//! let obs = Observation::new(
//!     &env,
//!     0,
//!     1.234,          // RA [rad]
//!     0.5 * RADSEC,   // σ_RA [rad]
//!     0.567,          // DEC [rad]
//!     0.5 * RADSEC,   // σ_DEC [rad]
//!     60300.0,        // MJD (TT)
//! )?;
//!
//! // Predict apparent position for this observation given an orbit
//! # let eq: EquinoctialElements = unimplemented!();
//! # let (_ra, _dec) = obs.compute_apparent_position(&env, &eq)?;
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## See also
//!
//! - [`initial_orbit_determination::gauss`] — Gauss IOD over observation triplets.
//! - [`observers`] — site database, Earth-fixed coordinates, and transformations.
//! - [`orbit_type::equinoctial_element::EquinoctialElements`] — propagation utilities used here.
//! - [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) and [`correct_aberration`](crate::observations::correct_aberration) — sky-projection helpers.
pub mod observations_ext;
pub mod triplets_generator;
pub mod triplets_iod;

use crate::{
    constants::{Observations, Radian, DPI, MJD, VLIGHT_AU},
    conversion::cartesian_to_radec,
    observers::Observer,
    orbit_type::equinoctial_element::EquinoctialElements,
    outfit::Outfit,
    outfit_errors::OutfitError,
};
use hifitime::Epoch;
use nalgebra::Vector3;
use std::f64::consts::PI;

/// Astrometric observation with site and precomputed observer positions.
///
/// This structure represents a single optical astrometric measurement
/// (right ascension/declination at a given epoch) together with:
/// - the associated observing site identifier,
/// - the observer’s **geocentric** position vector at the epoch, and
/// - the observer’s **heliocentric** position vector at the epoch.
///
/// Units & frames:
/// - Angles are stored in **radians**.
/// - Times are stored as **MJD (TT scale)**.
/// - Position vectors are expressed in **AU**, in the **equatorial mean J2000** frame.
///
/// Fields
/// -----------------
/// * `observer` – Site identifier (`u16`) referencing an [`Observer`] known by the [`Outfit`] state.
/// * `ra` – Right ascension `[rad]`.
/// * `error_ra` – Uncertainty on right ascension `[rad]`.
/// * `dec` – Declination `[rad]`.
/// * `error_dec` – Uncertainty on declination `[rad]`.
/// * `time` – Observation epoch as MJD (TT scale).
/// * `observer_earth_position` – Geocentric position of the observer at `time` (AU, equatorial mean J2000).
/// * `observer_helio_position` – Heliocentric position of the observer at `time` (AU, equatorial mean J2000).
#[derive(Debug, Clone, PartialEq, Copy)]
pub struct Observation {
    pub(crate) observer: u16,
    pub ra: Radian,
    pub error_ra: Radian,
    pub dec: Radian,
    pub error_dec: Radian,
    pub time: MJD,
    pub(crate) observer_earth_position: Vector3<f64>,
    pub(crate) observer_helio_position: Vector3<f64>,
}

impl Observation {
    /// Create a new astrometric observation and precompute observer positions.
    ///
    /// This constructor stores the astrometric angles and time, and computes the observer’s
    /// **geocentric** and **heliocentric** position vectors at the same epoch using the
    /// provided [`Outfit`] environment (UT1 provider, ephemerides, and site metadata).
    ///
    /// Arguments
    /// -----------------
    /// * `state` – Global environment providing ephemerides, UT1 provider and site database.
    /// * `observer` – Site identifier (`u16`) referencing an [`Observer`] known by `state`.
    /// * `ra` – Right ascension `[rad]`.
    /// * `error_ra` – Uncertainty on right ascension `[rad]`.
    /// * `dec` – Declination `[rad]`.
    /// * `error_dec` – Uncertainty on declination `[rad]`.
    /// * `time` – Observation epoch as **MJD (TT scale)**.
    ///
    /// Return
    /// ----------
    /// * A `Result` with the newly created [`Observation`], or an [`OutfitError`] if:
    ///   - the observer cannot be resolved in `state`,
    ///   - the UT1 provider / ephemeris computation fails.
    ///
    /// Remarks
    /// ------------
    /// * `pvobs` computes the geocentric position (and velocity) of the observer from Earth rotation and site coordinates.
    /// * `helio_position` converts the geocentric position to the heliocentric frame using the selected JPL ephemeris.
    /// * Both positions are expressed in **AU**, **equatorial mean J2000**.
    ///
    /// See also
    /// ------------
    /// * [`Observer::pvobs`] – Geocentric position/velocity of the observing site.
    /// * [`Observer::helio_position`] – Heliocentric position of the observing site.
    /// * [`crate::trajectories::batch_reader::ObservationBatch`] – Batch operations on observations.
    pub fn new(
        state: &Outfit,
        observer: u16,
        ra: Radian,
        error_ra: Radian,
        dec: Radian,
        error_dec: Radian,
        time: MJD,
    ) -> Result<Self, OutfitError> {
        // Observation time in TT
        let obs_mjd = Epoch::from_mjd_in_time_scale(time, hifitime::TimeScale::TT);
        let obs = state.get_observer_from_uint16(observer);
        let (geo_obs_pos, _) = obs.pvobs(&obs_mjd, state.get_ut1_provider())?;
        let helio_obs_pos = obs.helio_position(state, &obs_mjd, &geo_obs_pos)?;

        Ok(Observation {
            observer,
            ra,
            error_ra,
            dec,
            error_dec,
            time,
            observer_earth_position: geo_obs_pos,
            observer_helio_position: helio_obs_pos,
        })
    }

    /// Construct an [`Observation`] from precomputed observer positions.
    ///
    /// This constructor is a **fast-path alternative** to [`Observation::new`]:
    /// it skips all ephemeris calls by directly injecting the observer’s
    /// geocentric and heliocentric positions. This is useful in ingestion
    /// pipelines (e.g. Parquet readers) where positions can be cached and
    /// reused for multiple observations sharing the same `(observer, time)`.
    ///
    /// Arguments
    /// -----------------
    /// * `observer` – Packed observer identifier (`u16`).
    /// * `ra`, `error_ra` – Right ascension in radians and its 1-σ uncertainty (radians).
    /// * `dec`, `error_dec` – Declination in radians and its 1-σ uncertainty (radians).
    /// * `time` – Observation epoch in Modified Julian Date (TT scale).
    /// * `observer_earth_position` – Geocentric observer position vector (AU, equatorial mean J2000).
    /// * `observer_helio_position` – Heliocentric observer position vector (AU, equatorial mean J2000).
    ///
    /// Return
    /// ----------
    /// * A fully initialized [`Observation`] where astrometric quantities are set
    ///   and observer positions are trusted to be externally consistent.
    ///
    /// Remarks
    /// ------------
    /// * Use this constructor only when you can guarantee that positions were
    ///   computed consistently with the same environment (`Outfit`, UT1, ephemerides).
    /// * All getter methods behave identically to those of [`Observation::new`].
    ///
    /// See also
    /// ------------
    /// * [`Observation::new`] – Computes positions internally (slower, but self-contained).
    /// * [`Observer::pvobs`] – Routine for geocentric position/velocity.
    /// * [`Observer::helio_position`] – Routine for heliocentric position.
    #[allow(clippy::too_many_arguments)]
    pub fn with_positions(
        observer: u16,
        ra: Radian,
        error_ra: Radian,
        dec: Radian,
        error_dec: Radian,
        time: MJD,
        observer_earth_position: Vector3<f64>,
        observer_helio_position: Vector3<f64>,
    ) -> Self {
        Self {
            observer,
            ra,
            error_ra,
            dec,
            error_dec,
            time,
            observer_earth_position,
            observer_helio_position,
        }
    }

    /// Get the observer heliocentric position at the observation epoch.
    ///
    /// Arguments
    /// -----------------
    /// * *(none)* – Accessor method.
    ///
    /// Return
    /// ----------
    /// * A copy of the `3D` position vector (AU, equatorial mean J2000) of the observer
    ///   at `self.time` (MJD TT).
    ///
    /// See also
    /// ------------
    /// * [`Observation::new`] – Computes and stores this vector at construction.
    /// * [`Observer::helio_position`] – Underlying routine used to compute the value.
    pub fn get_observer_helio_position(&self) -> Vector3<f64> {
        self.observer_helio_position
    }

    /// Get the observer geocentric position at the observation epoch.
    ///
    /// Arguments
    /// -----------------
    /// * *(none)* – Accessor method.
    ///
    /// Return
    /// ----------
    /// * A copy of the `3D` position vector (AU, equatorial mean J2000) of the observer
    ///   relative to the Earth’s center at `self.time` (MJD TT).
    ///
    /// Remarks
    /// ------------
    /// * This vector is computed at construction via [`Observer::pvobs`].
    /// * Units are astronomical units (AU), in the equatorial mean J2000 frame.
    ///
    /// See also
    /// ------------
    /// * [`Observation::new`] – Computes and stores this vector at construction.
    /// * [`Observer::pvobs`] – Underlying routine used to compute the value.
    pub fn get_observer_earth_position(&self) -> Vector3<f64> {
        self.observer_earth_position
    }

    /// Get the observer from the observation
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    ///
    /// Return
    /// ------
    /// * The observer
    pub fn get_observer<'a>(&self, env_state: &'a Outfit) -> &'a Observer {
        env_state.get_observer_from_uint16(self.observer)
    }

    /// Compute the apparent equatorial coordinates (RA, DEC) of a solar system body
    /// as seen by this observation’s site at its epoch.
    ///
    /// Overview
    /// -----------------
    /// This method determines the apparent sky position of a target body,
    /// described by equinoctial orbital elements, as seen from the observing site
    /// corresponding to this [`Observation`].
    ///
    /// The computation steps are:
    /// 1. **Orbit propagation** – Propagate the body’s state from its reference epoch to the observation epoch using a two-body model.
    /// 2. **Reference frame handling** – Retrieve Earth’s barycentric position from the JPL ephemeris and transform to *ecliptic mean J2000*.
    /// 3. **Observer position** – Compute the observer’s heliocentric position (Earth + site geocentric offset).
    /// 4. **Light-time and aberration correction** – Form the observer–object vector and correct for aberration.
    /// 5. **Conversion to equatorial coordinates** – Convert the corrected line-of-sight vector to (RA, DEC).
    ///
    /// Arguments
    /// -----------------
    /// * `state` – Global environment providing ephemerides, UT1 provider, and frame utilities.
    /// * `equinoctial_element` – Orbital elements of the target body.
    ///
    /// Return
    /// ----------
    /// * `Result<(f64, f64), OutfitError>` – The apparent right ascension and declination `[rad]`.
    ///
    /// Units
    /// ----------
    /// * Positions: AU  
    /// * Velocities: AU/day  
    /// * Angles: radians  
    /// * Time: MJD TT  
    ///
    /// Errors
    /// ----------
    /// Returns [`OutfitError`] if:
    /// - Orbit propagation fails,
    /// - Ephemeris data is unavailable,
    /// - Reference-frame transformation fails.
    ///
    /// See also
    /// ------------
    /// * [`EquinoctialElements::solve_two_body_problem`] – Orbit propagation.
    /// * [`Observer::pvobs`] – Computes observer’s geocentric position.
    /// * [`correct_aberration`] – Aberration correction.
    /// * [`cartesian_to_radec`] – Convert Cartesian vectors to (RA, DEC).
    pub fn compute_apparent_position(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<(f64, f64), OutfitError> {
        // Hyperbolic/parabolic orbits (e >= 1) are not yet supported
        if equinoctial_element.eccentricity() >= 1.0 {
            return Err(OutfitError::InvalidOrbit(
                "Eccentricity >= 1 is not yet supported".to_string(),
            ));
        }

        // 1. Propagate asteroid position/velocity in ecliptic J2000
        let (cart_pos_ast, cart_pos_vel, _) = equinoctial_element.solve_two_body_problem(
            0.,
            self.time - equinoctial_element.reference_epoch,
            false,
        )?;

        // 2. Observation time in TT
        let obs_mjd = Epoch::from_mjd_in_time_scale(self.time, hifitime::TimeScale::TT);

        // 3. Earth's barycentric position in ecliptic J2000
        let (earth_position, _) = state.get_jpl_ephem()?.earth_ephemeris(&obs_mjd, false);

        // 4. get rotation from equatorial mean J2000 to ecliptic mean J2000
        let matrix_elc_transform = state.get_rot_equmj2000_to_eclmj2000();

        let earth_pos_eclj2000 = matrix_elc_transform.transpose() * earth_position;
        let cart_pos_ast_eclj2000 = matrix_elc_transform * cart_pos_ast;
        let cart_pos_vel_eclj2000 = matrix_elc_transform * cart_pos_vel;

        // 5. Observer heliocentric position
        let geo_obs_pos = self.observer_earth_position;
        let xobs = geo_obs_pos + earth_pos_eclj2000;
        let obs_on_earth = matrix_elc_transform * xobs;

        // 6. Relative position and aberration correction
        let relative_position = cart_pos_ast_eclj2000 - obs_on_earth;
        let corrected_pos = correct_aberration(relative_position, cart_pos_vel_eclj2000);

        // 7. Convert to RA/DEC
        let (alpha, delta, _) = cartesian_to_radec(corrected_pos);
        Ok((alpha, delta))
    }

    /// Compute the normalized squared astrometric residuals (RA, DEC)
    /// between an observed position and a propagated ephemeris.
    ///
    /// Overview
    /// -----------------
    /// This method compares the actual astrometric measurement stored in `self`
    /// against the expected position of the target body propagated from
    /// equinoctial elements.  
    /// It returns a scalar representing the sum of squared, normalized residuals
    /// in RA and DEC.
    ///
    /// Arguments
    /// -----------------
    /// * `state` – Global environment providing ephemerides and time conversions.
    /// * `equinoctial_element` – Orbital elements of the target body.
    ///
    /// Return
    /// ----------
    /// * `Result<f64, OutfitError>` – Dimensionless scalar value representing the weighted sum
    ///   of squared residuals. Equivalent to a chi² contribution for a single observation (without division by 2).
    ///
    /// Remarks
    /// ----------
    /// * Residuals are normalized by the astrometric uncertainties `error_ra` and `error_dec`.
    /// * RA residuals are multiplied by `cos(dec)` to account for projection effects.
    /// * All angles are in radians.
    ///
    /// Errors
    /// ----------
    /// Returns [`OutfitError`] if propagation or ephemeris lookup fails.
    ///
    /// See also
    /// ------------
    /// * [`compute_apparent_position`](crate::observations::Observation::compute_apparent_position) – Used internally to obtain predicted RA/DEC.
    /// * [`Observer::pvobs`] – Computes observer’s geocentric position.
    /// * [`correct_aberration`] – Applies aberration correction.
    /// * [`cartesian_to_radec`] – Converts 3D vectors to (RA, DEC).
    /// * [`EquinoctialElements::solve_two_body_problem`] – Two-body propagation.
    pub fn ephemeris_error(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<f64, OutfitError> {
        let (alpha, delta) = self.compute_apparent_position(state, equinoctial_element)?;

        // ΔRA with wrapping to [-π, π]
        let mut diff_alpha = (self.ra - alpha) % DPI;
        if diff_alpha > PI {
            diff_alpha -= DPI;
        }

        let diff_delta = self.dec - delta;

        // Weighted RMS
        let rms_ra = (self.dec.cos() * (diff_alpha / self.error_ra)).powi(2);
        let rms_dec = (diff_delta / self.error_dec).powi(2);

        Ok(rms_ra + rms_dec)
    }
}

/// Apply stellar aberration correction to a relative position vector.
///
/// This function computes the apparent position of a target object by applying
/// the first-order correction for stellar aberration due to the observer's velocity.
/// It assumes the classical limit (v ≪ c), using a linear time-delay model.
///
/// Arguments
/// ---------
/// * `xrel`: relative position vector from observer to object \[AU\].
/// * `vrel`: velocity of the observer relative to the barycenter \[AU/day\].
///
/// Returns
/// --------
/// * Corrected position vector (same units and directionality as `xrel`),
///   shifted by the aberration effect.
///
/// Formula
/// -------
/// The corrected position is given by:
/// ```text
/// x_corr = xrel − (‖xrel‖ / c) · vrel
/// ```
/// where `c` is the speed of light in AU/day (`VLIGHT_AU`).
///
/// Remarks
/// -------
/// * This function does **not** normalize the output.
/// * Suitable for use in astrometric modeling or when computing apparent direction
///   of celestial objects as seen from a moving observer.
pub fn correct_aberration(xrel: Vector3<f64>, vrel: Vector3<f64>) -> Vector3<f64> {
    let norm_vector = xrel.norm();
    let dt = norm_vector / VLIGHT_AU;
    xrel - dt * vrel
}

#[cfg(test)]
#[cfg(feature = "jpl-download")]
mod test_observations {

    use crate::unit_test_global::OUTFIT_HORIZON_TEST;

    use super::*;

    mod tests_compute_apparent_position {

        use super::*;
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;
        use approx::assert_relative_eq;

        /// Helper: simple circular equinoctial elements for a 1 AU, zero inclination orbit.
        fn simple_circular_elements(epoch: f64) -> EquinoctialElements {
            EquinoctialElements {
                reference_epoch: epoch,
                semi_major_axis: 1.0,
                eccentricity_sin_lon: 0.0,
                eccentricity_cos_lon: 0.0,
                tan_half_incl_sin_node: 0.0,
                tan_half_incl_cos_node: 0.0,
                mean_longitude: 0.0,
            }
        }

        #[test]
        fn test_compute_apparent_position_nominal() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_obs = 59000.0; // MJD
            let equinoctial = simple_circular_elements(t_obs);

            let obs = Observation::new(state, observer_code, 0.0, 0.0, 0.0, 0.0, t_obs).unwrap();

            let (ra, dec) = obs
                .compute_apparent_position(state, &equinoctial)
                .expect("Computation should succeed");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
            assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
            assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
        }

        #[test]
        fn test_compute_apparent_position_same_epoch() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_epoch = 60000.0;
            let equinoctial = simple_circular_elements(t_epoch);

            let obs = Observation::new(state, observer_code, 0.0, 0.0, 0.0, 0.0, t_epoch).unwrap();

            let (ra1, dec1) = obs.compute_apparent_position(state, &equinoctial).unwrap();
            let (ra2, dec2) = obs.compute_apparent_position(state, &equinoctial).unwrap();

            // The same input should always produce the same result
            assert_relative_eq!(ra1, ra2, epsilon = 1e-14);
            assert_relative_eq!(dec1, dec2, epsilon = 1e-14);
        }

        #[test]
        fn test_apparent_position_for_distant_object() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());
            let t_obs = 59000.0;
            let mut equinoctial = simple_circular_elements(t_obs);

            // Objet far away
            equinoctial.semi_major_axis = 100.0;

            let obs = Observation::new(state, observer_code, 0.0, 0.0, 0.0, 0.0, t_obs).unwrap();

            let (ra, dec) = obs
                .compute_apparent_position(state, &equinoctial)
                .expect("Should compute apparent position for distant object");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
        }

        #[test]
        fn test_compute_apparent_position_propagation_failure() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            // Invalid orbital elements to force failure in solve_two_body_problem
            let equinoctial = EquinoctialElements {
                reference_epoch: 59000.0,
                semi_major_axis: -1.0, // Physically invalid
                eccentricity_sin_lon: 0.0,
                eccentricity_cos_lon: 0.0,
                tan_half_incl_sin_node: 0.0,
                tan_half_incl_cos_node: 0.0,
                mean_longitude: 0.0,
            };

            let obs = Observation::new(state, observer_code, 0.0, 0.0, 0.0, 0.0, 59000.0).unwrap();

            let result = obs.compute_apparent_position(state, &equinoctial);
            assert!(result.is_err(), "Invalid elements should trigger an error");
        }

        mod proptests_apparent_position {
            use std::sync::Arc;

            use super::*;
            use proptest::prelude::*;

            /// Strategy: generates random but reasonable equinoctial elements
            /// for property-based tests.
            fn arb_equinoctial_elements() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0,                  // reference_epoch (MJD)
                    0.5..30.0,                         // semi-major axis (AU)
                    -0.5..0.5,                         // h = e * sin(Ω+ω)
                    -0.5..0.5,                         // k = e * cos(Ω+ω)
                    -0.5..0.5,                         // p = tan(i/2)*sin Ω
                    -0.5..0.5,                         // q = tan(i/2)*cos Ω
                    0.0..(2.0 * std::f64::consts::PI), // mean longitude (rad)
                )
                    .prop_map(|(epoch, a, h, k, p, q, lambda)| {
                        EquinoctialElements {
                            reference_epoch: epoch,
                            semi_major_axis: a,
                            eccentricity_sin_lon: h,
                            eccentricity_cos_lon: k,
                            tan_half_incl_sin_node: p,
                            tan_half_incl_cos_node: q,
                            mean_longitude: lambda,
                        }
                    })
            }

            /// Strategy: generates random observer locations on Earth.
            /// - longitude in [-180, 180] degrees
            /// - latitude in [-90, 90] degrees
            /// - elevation from 0 to 5 km
            fn arb_observer() -> impl Strategy<Value = Observer> {
                (-180.0..180.0, -90.0..90.0, 0.0..5.0).prop_map(|(lon, lat, elev)| {
                    Observer::new(lon, lat, elev, None, None, None).unwrap()
                })
            }

            /// Strategy: generates equinoctial elements with a wide range,
            /// including extreme eccentricities and inclinations.
            fn arb_extreme_equinoctial_elements() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0,
                    0.1..50.0,
                    -0.99..0.99,
                    -0.99..0.99,
                    -1.0..1.0,
                    -1.0..1.0,
                    0.0..(2.0 * std::f64::consts::PI),
                )
                    .prop_map(|(epoch, a, h, k, p, q, lambda)| EquinoctialElements {
                        reference_epoch: epoch,
                        semi_major_axis: a,
                        eccentricity_sin_lon: h,
                        eccentricity_cos_lon: k,
                        tan_half_incl_sin_node: p,
                        tan_half_incl_cos_node: q,
                        mean_longitude: lambda,
                    })
                    // Filter out cases where eccentricity >= 1
                    .prop_filter(
                        "Only bound (elliptical) orbits are supported",
                        |elem: &EquinoctialElements| {
                            let e = (elem.eccentricity_sin_lon.powi(2)
                                + elem.eccentricity_cos_lon.powi(2))
                            .sqrt();
                            e < 0.99
                        },
                    )
            }

            proptest! {
                /// Property test: RA and DEC are always finite and in the expected ranges.
                #[test]
                fn proptest_ra_dec_are_finite_and_in_range(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

                    let obs = Observation::new(
                        state,
                        observer_code,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         obs_time,
                    ).unwrap();

                    let result = obs.compute_apparent_position(state, &equinoctial);

                    if let Ok((ra, dec)) = result {
                        // Invariant: returned values must be finite
                        prop_assert!(ra.is_finite());
                        prop_assert!(dec.is_finite());

                        // RA must be in [0, 2π), DEC must be in [-π/2, π/2]
                        prop_assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
                        prop_assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
                    }
                }

                /// Property test: Calling the function twice with the same inputs must produce exactly
                /// the same output (determinism).
                #[test]
                fn proptest_repeatability(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

                    let obs = Observation::new(
                        state,
                        observer_code,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         obs_time,
                    ).unwrap();

                    let r1 = obs.compute_apparent_position(state, &equinoctial);
                    let r2 = obs.compute_apparent_position(state, &equinoctial);

                    // Invariant: repeated computation with the same input should be identical
                    prop_assert_eq!(r1, r2);
                }

                /// Property test: A very small change in observation time (1e-3 days ≈ 1.4 min)
                /// should not cause huge jumps in the resulting RA/DEC.
                #[test]
                fn proptest_small_time_change_has_small_effect(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

                    let obs = Observation::new(
                        state,
                        observer_code,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         obs_time,
                    ).unwrap();

                    let obs_eps = Observation {
                        time: obs_time + 1e-3, // shift by 1.4 minutes
                        ..obs
                    };

                    let r1 = obs.compute_apparent_position(state, &equinoctial);
                    let r2 = obs_eps.compute_apparent_position(state, &equinoctial);

                    if let (Ok((ra1, dec1)), Ok((ra2, dec2))) = (r1, r2) {
                        let dra = (ra1 - ra2).abs();
                        let ddec = (dec1 - dec2).abs();

                        // Invariant: no catastrophic jumps (> 1 radian) for a small time shift
                        prop_assert!(dra < 1.0);
                        prop_assert!(ddec < 1.0);
                    }
                }
            }

            proptest! {
                /// Property: RA/DEC remain finite and in valid ranges for extreme orbits and random observers.
                #[test]
                fn proptest_ra_dec_valid_for_extreme_orbits_and_observers(
                    equinoctial in arb_extreme_equinoctial_elements(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.add_observer_internal(Arc::new(observer));

                    let obs = Observation::new(
                        state,
                        observer_code,
                         0.0,
                         0.0,
                         0.0,
                         0.0,
                         obs_time,
                    ).unwrap();

                    let result = obs.compute_apparent_position(state, &equinoctial);

                    if let Ok((ra, dec)) = result {
                        // Values must be finite
                        prop_assert!(ra.is_finite());
                        prop_assert!(dec.is_finite());
                        // Angles must be within their valid intervals
                        prop_assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
                        prop_assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
                    }
                }
            }

            #[test]
            fn test_hyperbolic_orbit_returns_error() {
                let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

                let obs =
                    Observation::new(state, observer_code, 0.0, 0.0, 0.0, 0.0, 59000.0).unwrap();

                let equinoctial = EquinoctialElements {
                    reference_epoch: 59000.0,
                    semi_major_axis: 1.0,
                    eccentricity_sin_lon: 0.8,
                    eccentricity_cos_lon: 0.8, // e ≈ 1.13 > 1
                    tan_half_incl_sin_node: 0.0,
                    tan_half_incl_cos_node: 0.0,
                    mean_longitude: 0.0,
                };

                let result = obs.compute_apparent_position(state, &equinoctial);
                assert!(
                    result.is_err(),
                    "Hyperbolic or parabolic orbits should currently return an error"
                );
            }
        }
    }

    #[test]
    fn test_new_observation() {
        let state = &mut OUTFIT_HORIZON_TEST.0.clone();
        let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

        let observation =
            Observation::new(state, observer_code, 1.0, 0.1, 2.0, 0.2, 59000.0).unwrap();
        assert_eq!(
            observation,
            Observation {
                observer: 2,
                ra: 1.0,
                error_ra: 0.1,
                dec: 2.0,
                error_dec: 0.2,
                time: 59000.0,
                observer_earth_position: [
                    -1.4662164592060655e-6,
                    4.2560356749756634e-5,
                    -2.1126391698196086e-6
                ]
                .into(),
                observer_helio_position: [
                    -0.35113019606349866,
                    -0.8726512942340473,
                    -0.37829699890414364
                ]
                .into(),
            }
        );

        let observation_2 = Observation::new(
            state,
            2,
            343.097_375,
            2.777_777_777_777_778E-6,
            -14.784833333333333,
            2.777_777_777_777_778E-5,
            59001.0,
        )
        .unwrap();

        assert_eq!(
            observation_2,
            Observation {
                observer: 2,
                ra: 343.097375,
                error_ra: 2.777777777777778e-6,
                dec: -14.784833333333333,
                error_dec: 2.777777777777778e-5,
                time: 59001.0,
                observer_earth_position: [
                    -2.1521316017998277e-6,
                    4.2531873012231404e-5,
                    -2.0988352183088593e-6
                ]
                .into(),
                observer_helio_position: [
                    -0.33522248840408125,
                    -0.8780465618894304,
                    -0.380635845615707
                ]
                .into(),
            }
        );
    }

    mod tests_ephemeris_error {
        use super::*;
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;
        use approx::assert_relative_eq;

        fn simple_equinoctial(epoch: f64) -> EquinoctialElements {
            EquinoctialElements {
                reference_epoch: epoch,
                semi_major_axis: 1.0,
                eccentricity_sin_lon: 0.0,
                eccentricity_cos_lon: 0.0,
                tan_half_incl_sin_node: 0.0,
                tan_half_incl_cos_node: 0.0,
                mean_longitude: 0.0,
            }
        }

        #[test]
        fn test_ephem_error() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let obs = Observation::new(
                state,
                observer_code,
                1.7899347771316527,
                1.770_024_520_608_546E-6,
                0.778_996_538_107_973_6,
                1.259_582_891_829_317_7E-6,
                57070.262067592594,
            )
            .unwrap();

            let equinoctial_element = EquinoctialElements {
                reference_epoch: 57_049.242_334_573_75,
                semi_major_axis: 1.8017360713154256,
                eccentricity_sin_lon: 0.269_373_680_909_227_2,
                eccentricity_cos_lon: 8.856_415_260_013_56E-2,
                tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
                tan_half_incl_cos_node: 0.10168201109730375,
                mean_longitude: 1.6936970079414786,
            };

            let rms_error = obs.ephemeris_error(&OUTFIT_HORIZON_TEST.0, &equinoctial_element);
            assert_eq!(rms_error.unwrap(), 75.00445641224026);
        }

        /// When the observed RA/DEC exactly match the propagated RA/DEC,
        /// the ephemeris_error must be zero.
        #[test]
        fn test_zero_error_when_positions_match() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            // Compute the propagated position
            let obs = Observation::new(state, observer_code, 0.0, 1e-6, 0.0, 1e-6, t_obs).unwrap();
            let (alpha, delta) = obs.compute_apparent_position(state, &equinoctial).unwrap();

            // New observation with exact same RA/DEC
            let obs_match =
                Observation::new(state, observer_code, alpha, 1e-6, delta, 1e-6, t_obs).unwrap();

            let error = obs_match.ephemeris_error(state, &equinoctial).unwrap();
            assert_relative_eq!(error, 0.0, epsilon = 1e-14);
        }

        /// Error grows if RA is off by a known amount
        #[test]
        fn test_error_increases_with_offset() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let base_obs = Observation {
                observer: observer_code,
                ra: 0.0,
                error_ra: 1e-3,
                dec: 0.0,
                error_dec: 1e-3,
                time: t_obs,
                observer_earth_position: Vector3::zeros(),
                observer_helio_position: Vector3::zeros(),
            };
            let (alpha, delta) = base_obs
                .compute_apparent_position(state, &equinoctial)
                .unwrap();

            // Same dec, but RA offset by 1 milliradian
            let obs_offset = Observation {
                observer: 0,
                ra: alpha + 1e-3,
                error_ra: 1e-3,
                dec: delta,
                error_dec: 1e-3,
                time: t_obs,
                observer_earth_position: Vector3::zeros(),
                observer_helio_position: Vector3::zeros(),
            };

            let err = obs_offset.ephemeris_error(state, &equinoctial).unwrap();
            assert!(err > 0.0);
        }

        /// Check that wrapping of RA (close to 2π) does not affect the error
        #[test]
        fn test_ra_wrapping_invariance() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let base_obs =
                Observation::new(state, observer_code, 0.0, 1e-6, 0.0, 1e-6, t_obs).unwrap();
            let (alpha, delta) = base_obs
                .compute_apparent_position(state, &equinoctial)
                .unwrap();

            // Same position but RA shifted by ±2π
            let obs_wrapped = Observation::new(
                state,
                observer_code,
                alpha + std::f64::consts::TAU,
                1e-6,
                delta,
                1e-6,
                t_obs,
            )
            .unwrap();

            let err1 = obs_wrapped.ephemeris_error(state, &equinoctial).unwrap();
            assert_relative_eq!(err1, 0.0, epsilon = 1e-12);
        }

        /// When RA/DEC uncertainties are very large, error is small even with a mismatch.
        #[test]
        fn test_large_uncertainty_downweights_error() {
            let state = &mut OUTFIT_HORIZON_TEST.0.clone();
            let observer_code = state.uint16_from_mpc_code(&"F51".to_string());

            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let base_obs =
                Observation::new(state, observer_code, 0.0, 1.0, 0.0, 1.0, t_obs).unwrap();
            let (alpha, delta) = base_obs
                .compute_apparent_position(state, &equinoctial)
                .unwrap();

            let obs_large_uncertainty = Observation::new(
                state,
                observer_code,
                alpha + 0.1,
                10.0,
                delta + 0.1,
                10.0,
                t_obs,
            )
            .unwrap();

            let err = obs_large_uncertainty
                .ephemeris_error(state, &equinoctial)
                .unwrap();
            assert!(
                err < 1.0,
                "Large uncertainties should reduce the error contribution"
            );
        }

        mod proptests_ephemeris_error {
            use std::sync::Arc;

            use super::*;
            use proptest::prelude::*;

            fn arb_observer() -> impl Strategy<Value = Observer> {
                (-180.0..180.0, -90.0..90.0, 0.0..5.0).prop_map(|(lon, lat, elev)| {
                    Observer::new(lon, lat, elev, None, None, None).unwrap()
                })
            }

            fn arb_elliptical_equinoctial() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0,
                    0.5..20.0,
                    -0.8..0.8,
                    -0.8..0.8,
                    -0.8..0.8,
                    -0.8..0.8,
                    0.0..std::f64::consts::TAU,
                )
                    .prop_map(|(epoch, a, h, k, p, q, l)| EquinoctialElements {
                        reference_epoch: epoch,
                        semi_major_axis: a,
                        eccentricity_sin_lon: h,
                        eccentricity_cos_lon: k,
                        tan_half_incl_sin_node: p,
                        tan_half_incl_cos_node: q,
                        mean_longitude: l,
                    })
                    .prop_filter("Bound orbits only", |e: &EquinoctialElements| {
                        e.eccentricity() < 1.0
                    })
            }

            proptest! {
                /// Property: error is always non-negative and finite for valid inputs
                #[test]
                fn proptest_error_is_non_negative(
                    equinoctial in arb_elliptical_equinoctial(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.add_observer_internal(Arc::new(observer));
                    let obs = Observation::new(state, observer_code,
                         0.0,
                         1e-3,
                         0.0,
                         1e-3,
                         obs_time,).unwrap();

                    let result = obs.ephemeris_error(state, &equinoctial);
                    if let Ok(val) = result {
                        prop_assert!(val.is_finite());
                        prop_assert!(val >= 0.0);
                    }
                }

                /// Property: If uncertainties are huge, the error must be small
                #[test]
                fn proptest_error_downweights_large_uncertainties(
                    equinoctial in arb_elliptical_equinoctial(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &mut OUTFIT_HORIZON_TEST.0.clone();
                    let observer_code = state.add_observer_internal(Arc::new(observer));
                    let obs = Observation::new(state, observer_code,
                        0.5,
                        100.0,
                        0.5,
                        100.0,
                        obs_time,).unwrap();

                    let result = obs.ephemeris_error(state, &equinoctial);
                    if let Ok(val) = result {
                        prop_assert!(val < 1.0);
                    }
                }
            }
        }
    }
}
