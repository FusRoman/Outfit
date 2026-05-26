//! Apparent-position computation and astrometric residuals for individual observations.
//!
//! This module provides the [`ObservationEphemeris`](crate::observation_ephemeris::ObservationEphemeris) trait, which extends
//! [`photom::observation_dataset::observation::Observation`] with two methods:
//!
//! - [`ObservationEphemeris::compute_apparent_position`](crate::observation_ephemeris::ObservationEphemeris::compute_apparent_position) — propagates an orbit
//!   from its reference epoch to the observation epoch, applies the observer
//!   geometry and aberration correction, and returns the predicted (RA, DEC).
//! - [`ObservationEphemeris::ephemeris_error`](crate::observation_ephemeris::ObservationEphemeris::ephemeris_error) — computes the sum of squared,
//!   normalised astrometric residuals between the measured and predicted
//!   (RA, DEC), suitable as a χ² contribution.
//!
//! A stand-alone helper function [`correct_aberration`](crate::observation_ephemeris::correct_aberration) is also exported for
//! direct use when only the first-order aberration shift is needed.
//!
//! # Coordinate conventions
//!
//! - Positions are in **AU**, velocities in **AU/day**, angles in **radians**.
//! - Intermediate computations use the **ecliptic mean J2000** frame; the final
//!   apparent coordinates are returned in the **equatorial** frame (RA, DEC).

use std::f64::consts::PI;

use hifitime::Epoch;
use nalgebra::Vector3;
use photom::{constants::DPI, observation_dataset::observation::Observation};

use crate::{
    cache::OutfitCache, constants::ROT_EQUMJ2000_TO_ECLMJ2000, EquinoctialElements, JPLEphem,
    OutfitError, VLIGHT_AU,
};

/// Apparent equatorial coordinates of a solar system body together with their
/// partial derivatives with respect to the body's heliocentric position.
///
/// This struct is the return type of
/// [`ObservationEphemeris::compute_apparent_pos_derivative`].
///
/// ## Coordinate conventions
///
/// - `ra` and `dec` are in **radians**, equatorial mean J2000.
/// - `d_ra_d_pos` and `d_dec_d_pos` are in **rad/AU**, expressed in the
///   **ecliptic mean J2000** frame.  This matches the frame in which
///   [`EquinoctialElements::solve_two_body_problem`] returns its Jacobians,
///   so the chain rule can be applied directly:
///
/// ```text
/// ∂α/∂elemᵢ = d_ra_d_pos  · dpos/delemᵢ   (dot product, 3 components)
/// ∂δ/∂elemᵢ = d_dec_d_pos · dpos/delemᵢ
/// ```
#[derive(Debug, Clone, PartialEq)]
pub struct ApparentPositionWithPartials {
    /// Predicted right ascension \[rad\], in [0, 2π).
    pub ra: f64,
    /// Predicted declination \[rad\], in (−π/2, π/2).
    pub dec: f64,
    /// ∂α/∂(asteroid heliocentric position) \[rad/AU\], ecliptic mean J2000.
    pub d_ra_d_pos: Vector3<f64>,
    /// ∂δ/∂(asteroid heliocentric position) \[rad/AU\], ecliptic mean J2000.
    pub d_dec_d_pos: Vector3<f64>,
}

pub trait ObservationEphemeris {
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
    /// * [`ResolvedObserver::pvobs`](crate::observer_extension::ResolvedObserver::pvobs) – Computes observer's geocentric position.
    /// * [`correct_aberration`] – Aberration correction.
    /// * [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) – Convert Cartesian vectors to (RA, DEC).
    fn compute_apparent_position(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<(f64, f64), OutfitError>;

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
    /// * [`compute_apparent_position`](crate::observation_ephemeris::ObservationEphemeris::compute_apparent_position) – Used internally to obtain predicted RA/DEC.
    /// * [`ResolvedObserver::pvobs`](crate::observer_extension::ResolvedObserver::pvobs) – Computes observer's geocentric position.
    /// * [`correct_aberration`] – Applies aberration correction.
    /// * [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) – Converts 3D vectors to (RA, DEC).
    /// * [`EquinoctialElements::solve_two_body_problem`] – Two-body propagation.
    fn ephemeris_error(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<f64, OutfitError>;

    /// Compute the apparent equatorial coordinates (RA, DEC) of a solar system body
    /// together with the partial derivatives of (RA, DEC) with respect to the body's
    /// heliocentric Cartesian position (ecliptic mean J2000).
    ///
    /// Overview
    /// -----------------
    /// This is the astrometric analogue of `oss_dif()` in OrbFit.  It performs the
    /// same five steps as [`compute_apparent_position`](ObservationEphemeris::compute_apparent_position)
    /// and additionally returns the geometric Jacobians needed to assemble the design
    /// matrix G for the differential-correction least-squares loop.
    ///
    /// ## Computation steps
    ///
    /// 1. **Orbit propagation** — same as `compute_apparent_position`.
    /// 2. **Topocentric vector** — `d_ecl = ast_pos − obs_pos` (ecliptic J2000),
    ///    corrected for first-order aberration.
    /// 3. **Frame rotation** — `d_eq = R_ecl→eq · d_ecl` (equatorial J2000).
    /// 4. **Sky coordinates** — `α = atan2(d_eq.y, d_eq.x)`, `δ = asin(d_eq.z / ‖d_eq‖)`.
    /// 5. **Geometric Jacobians in equatorial frame**:
    ///    ```text
    ///    ∂α/∂pos_eq = (−d_eq.y / ρ_xy², d_eq.x / ρ_xy², 0)
    ///    ∂δ/∂pos_eq = (−d_eq.z·d_eq.x / (ρ_xy·ρ²),
    ///                  −d_eq.z·d_eq.y / (ρ_xy·ρ²),
    ///                   ρ_xy / ρ²)
    ///    ```
    ///    where `ρ_xy = √(d_eq.x² + d_eq.y²)` and `ρ = ‖d_eq‖`.
    /// 6. **Back-rotation to ecliptic** — multiply by `R_eq→ecl` so the Jacobians
    ///    are compatible with the state-transition matrix from `solve_two_body_problem`.
    ///
    /// ## Return
    ///
    /// An [`ApparentPositionWithPartials`] containing `(α, δ, ∂α/∂pos_ecl, ∂δ/∂pos_ecl)`.
    ///
    /// ## Errors
    ///
    /// Same as [`compute_apparent_position`](ObservationEphemeris::compute_apparent_position).
    ///
    /// ## Note on the polar singularity
    ///
    /// When the target is at or very near the celestial pole (`ρ_xy → 0`), the
    /// RA partial diverges.  This edge case is not handled — callers should avoid
    /// observations within a few arcseconds of the pole.
    fn compute_apparent_pos_derivative(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<ApparentPositionWithPartials, OutfitError>;
}

/// Compute topocentric (RA, DEC) and their partial derivatives w.r.t. the
/// asteroid's heliocentric position, all in the ecliptic mean J2000 frame.
///
/// This is the pure-geometry core shared by
/// [`ObservationEphemeris::compute_apparent_position`] and
/// [`ObservationEphemeris::compute_apparent_pos_derivative`].  It corresponds
/// to `oss_dif()` in OrbFit (`src/propag/pred_obs.f90`).
///
/// # Arguments
///
/// * `ast_pos_ecl` — Asteroid heliocentric position \[AU\], ecliptic J2000.
/// * `ast_vel_ecl` — Asteroid heliocentric velocity \[AU/day\], ecliptic J2000.
/// * `obs_pos_ecl` — Observer heliocentric position \[AU\], ecliptic J2000.
///
/// # Returns
///
/// `(α, δ, ∂α/∂pos_ecl, ∂δ/∂pos_ecl)` where:
/// - `α` ∈ \[0, 2π), `δ` ∈ (−π/2, π/2) \[rad\]
/// - Jacobians are in \[rad/AU\], expressed in the ecliptic J2000 frame.
fn topocentric_radec_and_partials(
    ast_pos_ecl: Vector3<f64>,
    ast_vel_ecl: Vector3<f64>,
    obs_pos_ecl: Vector3<f64>,
) -> (f64, f64, Vector3<f64>, Vector3<f64>) {
    // 1. Topocentric vector + aberration correction.
    //    The result is in the same frame as ast_pos_ecl / obs_pos_ecl.
    //    Following the original pipeline (compatible with the old cartesion_from_vec path),
    //    the RA/Dec angles are computed directly from the (x,y,z) components of `corrected`
    //    using atan2 — no additional frame rotation is applied here.
    let relative = ast_pos_ecl - obs_pos_ecl;
    let corrected = correct_aberration(relative, ast_vel_ecl);

    let x = corrected[0];
    let y = corrected[1];
    let z = corrected[2];

    let rho = corrected.norm(); // topocentric distance
    let rho_xy = x.hypot(y); // ρ_xy = √(x²+y²), matches original hypot usage
    let rho_xy_sq = rho_xy * rho_xy; // ρ_xy²

    let dec = z.atan2(rho_xy); // matches original EquCoord::from(CartesianCoord) formula
    let ra = y.atan2(x).rem_euclid(std::f64::consts::TAU);

    // 3. Geometric Jacobians w.r.t. ast_pos_ecl.
    //
    //    Full chain rule through correct_aberration:
    //      corrected = relative − (‖relative‖ / c) · vel
    //    ⟹  ∂corrected/∂pos = I − (1/c) · vel ⊗ (relative / ‖relative‖)ᵀ
    //
    //    Applying the chain rule:
    //      ∂α/∂pos = grad_ra − (grad_ra · vel) / (‖relative‖·c) · relative
    //      ∂δ/∂pos = grad_dec − (grad_dec · vel) / (‖relative‖·c) · relative
    //
    //    where the gradients w.r.t. `corrected` are:
    //      grad_ra  = (−y/ρ_xy², x/ρ_xy², 0)
    //      grad_dec = (−z·x/(ρ_xy·ρ²), −z·y/(ρ_xy·ρ²), ρ_xy/ρ²)
    let rho_sq = rho * rho;
    let grad_ra = Vector3::new(-y / rho_xy_sq, x / rho_xy_sq, 0.0);
    let grad_dec = Vector3::new(
        -z * x / (rho_xy * rho_sq),
        -z * y / (rho_xy * rho_sq),
        rho_xy / rho_sq,
    );

    // Aberration correction factor: 1 / (‖relative‖ · c)
    let rel_norm = relative.norm();
    let aberr_factor = 1.0 / (rel_norm * VLIGHT_AU);

    let d_ra_d_pos = grad_ra - grad_ra.dot(&ast_vel_ecl) * aberr_factor * relative;
    let d_dec_d_pos = grad_dec - grad_dec.dot(&ast_vel_ecl) * aberr_factor * relative;

    (ra, dec, d_ra_d_pos, d_dec_d_pos)
}

impl ObservationEphemeris for Observation {
    fn compute_apparent_position(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
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
            self.mjd_tt() - equinoctial_element.reference_epoch,
            false,
        )?;

        // 2. Observation time in TT
        let obs_mjd = Epoch::from_mjd_in_time_scale(self.mjd_tt(), hifitime::TimeScale::TT);

        // 3. Earth's barycentric position in ecliptic J2000
        let (earth_position, _) = jpl.earth_ephemeris(&obs_mjd, false);

        let earth_pos_eclj2000 = ROT_EQUMJ2000_TO_ECLMJ2000.transpose() * earth_position;
        let ast_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_ast;
        let ast_vel_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_vel;

        // 4. Observer heliocentric position in ecliptic J2000
        let geo_obs_pos = cache
            .get_observer_geocentric_position(self.index())
            .map(|x| x.into_inner());
        let xobs = geo_obs_pos + earth_pos_eclj2000;
        let obs_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * xobs;

        // 5. (α, δ) via shared geometry core (partials discarded)
        let (ra, dec, _, _) = topocentric_radec_and_partials(ast_pos_ecl, ast_vel_ecl, obs_pos_ecl);

        Ok((ra, dec))
    }

    fn compute_apparent_pos_derivative(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<ApparentPositionWithPartials, OutfitError> {
        // Hyperbolic/parabolic orbits (e >= 1) are not yet supported
        if equinoctial_element.eccentricity() >= 1.0 {
            return Err(OutfitError::InvalidOrbit(
                "Eccentricity >= 1 is not yet supported".to_string(),
            ));
        }

        // 1. Propagate asteroid position/velocity in ecliptic J2000
        let (cart_pos_ast, cart_pos_vel, _) = equinoctial_element.solve_two_body_problem(
            0.,
            self.mjd_tt() - equinoctial_element.reference_epoch,
            false,
        )?;

        // 2. Observation time in TT
        let obs_mjd = Epoch::from_mjd_in_time_scale(self.mjd_tt(), hifitime::TimeScale::TT);

        // 3. Earth's barycentric position in ecliptic J2000
        let (earth_position, _) = jpl.earth_ephemeris(&obs_mjd, false);

        let earth_pos_eclj2000 = ROT_EQUMJ2000_TO_ECLMJ2000.transpose() * earth_position;
        let ast_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_ast;
        let ast_vel_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_vel;

        // 4. Observer heliocentric position in ecliptic J2000
        let geo_obs_pos = cache
            .get_observer_geocentric_position(self.index())
            .map(|x| x.into_inner());
        let xobs = geo_obs_pos + earth_pos_eclj2000;
        let obs_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * xobs;

        // 5. (α, δ) + partial derivatives via shared geometry core
        let (ra, dec, d_ra_d_pos, d_dec_d_pos) =
            topocentric_radec_and_partials(ast_pos_ecl, ast_vel_ecl, obs_pos_ecl);

        Ok(ApparentPositionWithPartials {
            ra,
            dec,
            d_ra_d_pos,
            d_dec_d_pos,
        })
    }

    fn ephemeris_error(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<f64, OutfitError> {
        let (alpha, delta) = self.compute_apparent_position(cache, jpl, equinoctial_element)?;

        let (self_ra, self_ra_err, self_dec, self_dec_err) = (
            self.equ_coord().ra,
            self.equ_coord().ra_error,
            self.equ_coord().dec,
            self.equ_coord().dec_error,
        );

        // ΔRA with wrapping to [-π, π]
        let mut diff_alpha = (self_ra - alpha) % DPI;
        if diff_alpha > PI {
            diff_alpha -= DPI;
        }

        let diff_delta = self_dec - delta;

        // Weighted RMS
        let rms_ra = (self_dec.cos() * (diff_alpha / self_ra_err)).powi(2);
        let rms_dec = (diff_delta / self_dec_err).powi(2);

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
mod test_observations_ephemeris {
    use super::*;

    mod tests_compute_apparent_position {

        use crate::test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER};

        use super::*;
        use approx::assert_relative_eq;
        use photom::{
            coordinates::equatorial::EquCoord,
            observation_dataset::{observation::ObservationInput, ObsDataset},
            observer::error_model::{ModelCorrection, ObsErrorModel},
            photometry::{Filter, Photometry},
            MJDTT,
        };

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

        fn obsdataset_with_observation_time(t_epoch: MJDTT) -> ObsDataset {
            let observation_input = ObservationInput::new(
                0,
                EquCoord {
                    ra: 0.0,
                    ra_error: 0.0,
                    dec: 0.0,
                    dec_error: 0.0,
                },
                Photometry {
                    magnitude: 0.0,
                    error: 0.0,
                    filter: Filter::Int(0),
                },
                t_epoch,
                Some(photom::observer::dataset::ObserverId::MpcCode(*b"F51")),
            );

            ObsDataset::empty()
                .push_observation(vec![observation_input])
                .unwrap()
                .0
                .with_error_model(ObsErrorModel::FCCT14)
                .apply_model_errors()
        }

        #[test]
        fn test_compute_apparent_position_nominal() {
            let t_obs = 59000.0; // MJD

            let obs_dataset = obsdataset_with_observation_time(t_obs);
            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();

            let equinoctial = simple_circular_elements(t_obs);

            let (ra, dec) = obs_dataset
                .get_observation(0)
                .unwrap()
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .expect("Computation should succeed");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
            assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
            assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
        }

        #[test]
        fn test_compute_apparent_position_same_epoch() {
            let t_epoch = 60000.0;

            let obs_dataset = obsdataset_with_observation_time(t_epoch);
            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();

            let equinoctial = simple_circular_elements(t_epoch);

            let obs = obs_dataset.get_observation(0).unwrap();

            let (ra1, dec1) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();
            let (ra2, dec2) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            // The same input should always produce the same result
            assert_relative_eq!(ra1, ra2, epsilon = 1e-14);
            assert_relative_eq!(dec1, dec2, epsilon = 1e-14);
        }

        #[test]
        fn test_apparent_position_for_distant_object() {
            let t_obs = 59000.0;

            let obs_dataset = obsdataset_with_observation_time(t_obs);
            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();

            let mut equinoctial = simple_circular_elements(t_obs);

            // Objet far away
            equinoctial.semi_major_axis = 100.0;

            let obs = obs_dataset.get_observation(0).unwrap();

            let (ra, dec) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .expect("Should compute apparent position for distant object");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
        }

        #[test]
        fn test_compute_apparent_position_propagation_failure() {
            let obs_dataset = obsdataset_with_observation_time(0.0);
            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();

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

            let obs = obs_dataset.get_observation(0).unwrap();

            let result = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);
            assert!(result.is_err(), "Invalid elements should trigger an error");
        }

        mod proptests_apparent_position {
            use super::*;
            use crate::test_fixture::UT1_PROVIDER;
            use photom::{
                observation_dataset::{observation::ObservationInput, ObsDataset},
                observer::{
                    dataset::ObserverId,
                    error_model::{ModelCorrection, ObsErrorModel},
                    Observer,
                },
                photometry::{Filter, Photometry},
            };
            use proptest::prelude::*;

            fn arb_equinoctial_elements() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0f64,
                    0.5..30.0f64,
                    -0.5..0.5f64,
                    -0.5..0.5f64,
                    -0.5..0.5f64,
                    -0.5..0.5f64,
                    0.0..(2.0 * std::f64::consts::PI),
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

            fn arb_observer() -> impl Strategy<Value = Observer> {
                (-180.0..180.0f64, -90.0..90.0f64, 0.0..5.0f64).prop_map(|(lon, lat, elev)| {
                    Observer::new(lon.to_radians(), lat.to_radians(), elev, None, None, None)
                        .unwrap()
                })
            }

            fn arb_extreme_equinoctial_elements() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0f64,
                    0.1..50.0f64,
                    -0.99..0.99f64,
                    -0.99..0.99f64,
                    -1.0..1.0f64,
                    -1.0..1.0f64,
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

            fn make_obs_dataset_and_cache(
                t_obs: f64,
                observer_id: ObserverId,
            ) -> (ObsDataset, OutfitCache) {
                let observation_input = ObservationInput::new(
                    0,
                    EquCoord {
                        ra: 0.0,
                        ra_error: 0.0,
                        dec: 0.0,
                        dec_error: 0.0,
                    },
                    Photometry {
                        magnitude: 0.0,
                        error: 0.0,
                        filter: Filter::Int(0),
                    },
                    t_obs,
                    Some(observer_id),
                );

                let obs_dataset = ObsDataset::empty()
                    .push_observation(vec![observation_input])
                    .unwrap()
                    .0
                    .with_error_model(ObsErrorModel::FCCT14)
                    .apply_model_errors();

                let cache =
                    OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false)
                        .unwrap();

                (obs_dataset, cache)
            }

            fn make_obs_dataset_and_cache_with_custom_observer(
                t_obs: f64,
                observer: Observer,
            ) -> (ObsDataset, OutfitCache) {
                let (obs_dataset_with_obs, observer_id) =
                    ObsDataset::empty().push_observer(observer);

                let observation_input = ObservationInput::new(
                    0,
                    EquCoord {
                        ra: 0.0,
                        ra_error: 0.0,
                        dec: 0.0,
                        dec_error: 0.0,
                    },
                    Photometry {
                        magnitude: 0.0,
                        error: 0.0,
                        filter: Filter::Int(0),
                    },
                    t_obs,
                    Some(observer_id),
                );

                let obs_dataset = obs_dataset_with_obs
                    .push_observation(vec![observation_input])
                    .unwrap()
                    .0
                    .with_error_model(ObsErrorModel::FCCT14)
                    .apply_model_errors();

                let cache =
                    OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false)
                        .unwrap();

                (obs_dataset, cache)
            }

            proptest! {
                #[test]
                fn proptest_ra_dec_are_finite_and_in_range(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) = make_obs_dataset_and_cache(
                        obs_time,
                        photom::observer::dataset::ObserverId::MpcCode(*b"F51"),
                    );

                    let obs = obs_dataset.get_observation(0).unwrap();
                    let result = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                    if let Ok((ra, dec)) = result {
                        prop_assert!(ra.is_finite());
                        prop_assert!(dec.is_finite());
                        prop_assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
                        prop_assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
                    }
                }

                #[test]
                fn proptest_repeatability(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) = make_obs_dataset_and_cache(
                        obs_time,
                        photom::observer::dataset::ObserverId::MpcCode(*b"F51"),
                    );

                    let obs = obs_dataset.get_observation(0).unwrap();

                    let r1 = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);
                    let r2 = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                    prop_assert_eq!(r1, r2);
                }

                #[test]
                fn proptest_small_time_change_has_small_effect(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) = make_obs_dataset_and_cache(
                        obs_time,
                        photom::observer::dataset::ObserverId::MpcCode(*b"F51"),
                    );
                    let (obs_dataset_eps, cache_eps) = make_obs_dataset_and_cache(
                        obs_time + 1e-3,
                        photom::observer::dataset::ObserverId::MpcCode(*b"F51"),
                    );

                    let obs = obs_dataset.get_observation(0).unwrap();
                    let obs_eps = obs_dataset_eps.get_observation(0).unwrap();

                    let r1 = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);
                    let r2 = obs_eps.compute_apparent_position(&cache_eps, &JPL_EPHEM_HORIZON, &equinoctial);

                    if let (Ok((ra1, dec1)), Ok((ra2, dec2))) = (r1, r2) {
                        let dra = (ra1 - ra2).abs();
                        let ddec = (dec1 - dec2).abs();

                        prop_assert!(dra < 1.0, "RA jump too large: {}", dra);
                        prop_assert!(ddec < 1.0, "DEC jump too large: {}", ddec);
                    }
                }

                #[test]
                fn proptest_ra_dec_valid_for_extreme_orbits_and_observers(
                    equinoctial in arb_extreme_equinoctial_elements(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) =
                        make_obs_dataset_and_cache_with_custom_observer(obs_time, observer);

                    let obs = obs_dataset.get_observation(0).unwrap();
                    let result = obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                    if let Ok((ra, dec)) = result {
                        prop_assert!(ra.is_finite());
                        prop_assert!(dec.is_finite());
                        prop_assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
                        prop_assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
                    }
                }
            }

            #[test]
            fn test_hyperbolic_orbit_returns_error() {
                let t_obs = 59000.0;
                let (obs_dataset, cache) = make_obs_dataset_and_cache(
                    t_obs,
                    photom::observer::dataset::ObserverId::MpcCode(*b"F51"),
                );

                let equinoctial = EquinoctialElements {
                    reference_epoch: t_obs,
                    semi_major_axis: 1.0,
                    eccentricity_sin_lon: 0.8,
                    eccentricity_cos_lon: 0.8, // e ≈ 1.13 > 1
                    tan_half_incl_sin_node: 0.0,
                    tan_half_incl_cos_node: 0.0,
                    mean_longitude: 0.0,
                };

                let obs = obs_dataset.get_observation(0).unwrap();
                let result =
                    obs.compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                assert!(
                    result.is_err(),
                    "Hyperbolic or parabolic orbits should currently return an error"
                );
            }
        }
    }

    mod tests_ephemeris_error {
        use super::*;
        use crate::test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER};
        use approx::assert_relative_eq;
        use photom::{
            coordinates::equatorial::EquCoord,
            observation_dataset::{observation::ObservationInput, ObsDataset},
            observer::{
                error_model::{ModelCorrection, ObsErrorModel},
                mpc::MpcCode,
                Observer,
            },
            photometry::{Filter, Photometry},
        };

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

        fn make_obs_dataset_and_cache_mpc(
            ra: f64,
            ra_error: f64,
            dec: f64,
            dec_error: f64,
            t_obs: f64,
            mpc_code: MpcCode,
            apply_model_errors: bool,
        ) -> (ObsDataset, OutfitCache) {
            let observation_input = ObservationInput::new(
                0,
                EquCoord {
                    ra,
                    ra_error,
                    dec,
                    dec_error,
                },
                Photometry {
                    magnitude: 0.0,
                    error: 0.0,
                    filter: Filter::Int(0),
                },
                t_obs,
                Some(photom::observer::dataset::ObserverId::MpcCode(mpc_code)),
            );

            let obs_dataset = ObsDataset::empty()
                .push_observation(vec![observation_input])
                .unwrap()
                .0
                .with_error_model(ObsErrorModel::FCCT14);

            let obs_dataset = if apply_model_errors {
                obs_dataset.apply_model_errors()
            } else {
                obs_dataset
            };

            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();
            (obs_dataset, cache)
        }

        fn make_obs_dataset_and_cache_custom(
            ra: f64,
            ra_error: f64,
            dec: f64,
            dec_error: f64,
            t_obs: f64,
            observer: Observer,
        ) -> (ObsDataset, OutfitCache) {
            let (dataset_with_obs, observer_id) = ObsDataset::empty().push_observer(observer);

            let observation_input = ObservationInput::new(
                0,
                EquCoord {
                    ra,
                    ra_error,
                    dec,
                    dec_error,
                },
                Photometry {
                    magnitude: 0.0,
                    error: 0.0,
                    filter: Filter::Int(0),
                },
                t_obs,
                Some(observer_id),
            );

            let obs_dataset = dataset_with_obs
                .push_observation(vec![observation_input])
                .unwrap()
                .0
                .with_error_model(ObsErrorModel::FCCT14)
                .apply_model_errors();

            let cache =
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();
            (obs_dataset, cache)
        }

        #[test]
        fn test_ephem_error() {
            let (obs_dataset, cache) = make_obs_dataset_and_cache_mpc(
                1.7899347771316527,
                1.770_024_520_608_546E-6,
                0.778_996_538_107_973_6,
                1.259_582_891_829_317_7E-6,
                57070.262067592594,
                *b"F51",
                false,
            );

            let equinoctial_element = EquinoctialElements {
                reference_epoch: 57_049.242_334_573_75,
                semi_major_axis: 1.8017360713154256,
                eccentricity_sin_lon: 0.269_373_680_909_227_2,
                eccentricity_cos_lon: 8.856_415_260_013_56E-2,
                tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
                tan_half_incl_cos_node: 0.10168201109730375,
                mean_longitude: 1.6936970079414786,
            };

            let obs = obs_dataset.get_observation(0).unwrap();
            let rms_error = obs.ephemeris_error(&cache, &JPL_EPHEM_HORIZON, &equinoctial_element);
            assert_eq!(rms_error.unwrap(), 75.00445641224026);
        }

        #[test]
        fn test_zero_error_when_positions_match() {
            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let (obs_dataset, cache) =
                make_obs_dataset_and_cache_mpc(0.0, 1e-6, 0.0, 1e-6, t_obs, *b"F51", true);

            let obs = obs_dataset.get_observation(0).unwrap();
            let (alpha, delta) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            let (obs_dataset_match, cache_match) =
                make_obs_dataset_and_cache_mpc(alpha, 1e-6, delta, 1e-6, t_obs, *b"F51", true);

            let obs_match = obs_dataset_match.get_observation(0).unwrap();
            let error = obs_match
                .ephemeris_error(&cache_match, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            assert_relative_eq!(error, 0.0, epsilon = 1e-14);
        }

        #[test]
        fn test_error_increases_with_offset() {
            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let (obs_dataset, cache) =
                make_obs_dataset_and_cache_mpc(0.0, 1e-3, 0.0, 1e-3, t_obs, *b"F51", true);

            let obs = obs_dataset.get_observation(0).unwrap();
            let (alpha, delta) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            let (obs_dataset_offset, cache_offset) = make_obs_dataset_and_cache_mpc(
                alpha + 1e-3,
                1e-3,
                delta,
                1e-3,
                t_obs,
                *b"F51",
                true,
            );

            let obs_offset = obs_dataset_offset.get_observation(0).unwrap();
            let err = obs_offset
                .ephemeris_error(&cache_offset, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            assert!(err > 0.0);
        }

        #[test]
        fn test_ra_wrapping_invariance() {
            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let (obs_dataset, cache) =
                make_obs_dataset_and_cache_mpc(0.0, 1e-6, 0.0, 1e-6, t_obs, *b"F51", true);

            let obs = obs_dataset.get_observation(0).unwrap();
            let (alpha, delta) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            let (obs_dataset_wrapped, cache_wrapped) = make_obs_dataset_and_cache_mpc(
                alpha + std::f64::consts::TAU,
                1e-6,
                delta,
                1e-6,
                t_obs,
                *b"F51",
                true,
            );

            let obs_wrapped = obs_dataset_wrapped.get_observation(0).unwrap();
            let err = obs_wrapped
                .ephemeris_error(&cache_wrapped, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            assert_relative_eq!(err, 0.0, epsilon = 1e-12);
        }

        #[test]
        fn test_large_uncertainty_downweights_error() {
            let t_obs = 59000.0;
            let equinoctial = simple_equinoctial(t_obs);

            let (obs_dataset, cache) =
                make_obs_dataset_and_cache_mpc(0.0, 1.0, 0.0, 1.0, t_obs, *b"F51", true);

            let obs = obs_dataset.get_observation(0).unwrap();
            let (alpha, delta) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            let (obs_dataset_large, cache_large) = make_obs_dataset_and_cache_mpc(
                alpha + 0.1,
                10.0,
                delta + 0.1,
                10.0,
                t_obs,
                *b"F51",
                true,
            );

            let obs_large = obs_dataset_large.get_observation(0).unwrap();
            let err = obs_large
                .ephemeris_error(&cache_large, &JPL_EPHEM_HORIZON, &equinoctial)
                .unwrap();

            assert!(
                err < 1.0,
                "Large uncertainties should reduce the error contribution"
            );
        }

        mod proptests_ephemeris_error {
            use super::*;
            use proptest::prelude::*;

            fn arb_observer() -> impl Strategy<Value = Observer> {
                (-180.0..180.0f64, -90.0..90.0f64, 0.0..5000.0f64).prop_map(|(lon, lat, elev)| {
                    Observer::new(lon.to_radians(), lat.to_radians(), elev, None, None, None)
                        .unwrap()
                })
            }

            fn arb_elliptical_equinoctial() -> impl Strategy<Value = EquinoctialElements> {
                (
                    58000.0..62000.0f64,
                    0.5..20.0f64,
                    -0.8..0.8f64,
                    -0.8..0.8f64,
                    -0.8..0.8f64,
                    -0.8..0.8f64,
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
                #[test]
                fn proptest_error_is_non_negative(
                    equinoctial in arb_elliptical_equinoctial(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) =
                        make_obs_dataset_and_cache_custom(0.0, 1e-3, 0.0, 1e-3, obs_time, observer);

                    let obs = obs_dataset.get_observation(0).unwrap();
                    let result = obs.ephemeris_error(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                    if let Ok(val) = result {
                        prop_assert!(val.is_finite());
                        prop_assert!(val >= 0.0);
                    }
                }

                #[test]
                fn proptest_error_downweights_large_uncertainties(
                    equinoctial in arb_elliptical_equinoctial(),
                    observer in arb_observer(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let (obs_dataset, cache) =
                        make_obs_dataset_and_cache_custom(0.5, 100.0, 0.5, 100.0, obs_time, observer);

                    let obs = obs_dataset.get_observation(0).unwrap();
                    let result = obs.ephemeris_error(&cache, &JPL_EPHEM_HORIZON, &equinoctial);

                    if let Ok(val) = result {
                        prop_assert!(val < 1.0);
                    }
                }
            }
        }
    }

    mod tests_compute_apparent_pos_derivative {
        use super::*;
        use crate::test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER};
        use approx::assert_relative_eq;
        use photom::{
            coordinates::equatorial::EquCoord,
            observation_dataset::{observation::ObservationInput, ObsDataset},
            observer::error_model::{ModelCorrection, ObsErrorModel},
            photometry::{Filter, Photometry},
            MJDTT,
        };

        fn make_obs_and_cache(t_obs: MJDTT) -> (ObsDataset, OutfitCache) {
            let input = ObservationInput::new(
                0,
                EquCoord {
                    ra: 0.0,
                    ra_error: 1e-6,
                    dec: 0.0,
                    dec_error: 1e-6,
                },
                Photometry {
                    magnitude: 0.0,
                    error: 0.0,
                    filter: Filter::Int(0),
                },
                t_obs,
                Some(photom::observer::dataset::ObserverId::MpcCode(*b"F51")),
            );
            let ds = ObsDataset::empty()
                .push_observation(vec![input])
                .unwrap()
                .0
                .with_error_model(ObsErrorModel::FCCT14)
                .apply_model_errors();
            let cache = OutfitCache::build(&ds, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();
            (ds, cache)
        }

        fn nominal_elements(epoch: f64) -> EquinoctialElements {
            EquinoctialElements {
                reference_epoch: epoch,
                semi_major_axis: 1.8,
                eccentricity_sin_lon: 0.1,
                eccentricity_cos_lon: 0.05,
                tan_half_incl_sin_node: 0.01,
                tan_half_incl_cos_node: 0.08,
                mean_longitude: 1.5,
            }
        }

        /// (α, δ) from compute_apparent_pos_derivative must match
        /// compute_apparent_position exactly.
        #[test]
        fn test_ra_dec_consistent_with_apparent_position() {
            let t_obs = 59000.0;
            let (ds, cache) = make_obs_and_cache(t_obs);
            let elem = nominal_elements(t_obs);
            let obs = ds.get_observation(0).unwrap();

            let (ra, dec) = obs
                .compute_apparent_position(&cache, &JPL_EPHEM_HORIZON, &elem)
                .unwrap();
            let result = obs
                .compute_apparent_pos_derivative(&cache, &JPL_EPHEM_HORIZON, &elem)
                .unwrap();

            assert_relative_eq!(result.ra, ra, epsilon = 1e-14);
            assert_relative_eq!(result.dec, dec, epsilon = 1e-14);
        }

        /// Validate ∂α/∂pos and ∂δ/∂pos by finite differences on ast_pos_ecl.
        ///
        /// We perturb each of the 3 ecliptic position components by ε and compare
        /// the numerical gradient to the analytical one.
        #[test]
        fn test_partials_match_finite_differences() {
            use crate::constants::ROT_EQUMJ2000_TO_ECLMJ2000;

            let t_obs = 59000.0;
            let elem = nominal_elements(t_obs);

            let (cart_pos_ast, cart_vel_ast, _) =
                elem.solve_two_body_problem(0., 0., false).unwrap();

            // Use the same frame as the impl
            let obs_pos_ecl = Vector3::zeros(); // simplify: observer at origin
            let ast_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_ast;
            let ast_vel_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_vel_ast;

            let (_, _, d_ra_d_pos, d_dec_d_pos) =
                topocentric_radec_and_partials(ast_pos_ecl, ast_vel_ecl, obs_pos_ecl);

            let eps = 1e-6_f64; // 1e-6 AU ≈ 150 km

            for i in 0..3 {
                // Central finite differences for better accuracy
                let mut pos_fwd = ast_pos_ecl;
                pos_fwd[i] += eps;
                let mut pos_bwd = ast_pos_ecl;
                pos_bwd[i] -= eps;

                let (ra_fwd, dec_fwd, _, _) =
                    topocentric_radec_and_partials(pos_fwd, ast_vel_ecl, obs_pos_ecl);
                let (ra_bwd, dec_bwd, _, _) =
                    topocentric_radec_and_partials(pos_bwd, ast_vel_ecl, obs_pos_ecl);

                // Handle RA wrapping: unwrap the difference to [-π, π]
                let mut dra = ra_fwd - ra_bwd;
                if dra > std::f64::consts::PI {
                    dra -= std::f64::consts::TAU;
                } else if dra < -std::f64::consts::PI {
                    dra += std::f64::consts::TAU;
                }
                let num_dra = dra / (2.0 * eps);
                let num_ddec = (dec_fwd - dec_bwd) / (2.0 * eps);

                assert_relative_eq!(d_ra_d_pos[i], num_dra, epsilon = 1e-6, max_relative = 1e-5);
                assert_relative_eq!(
                    d_dec_d_pos[i],
                    num_ddec,
                    epsilon = 1e-6,
                    max_relative = 1e-5
                );
            }
        }

        /// Hyperbolic orbit must return an error (same guard as compute_apparent_position).
        #[test]
        fn test_hyperbolic_orbit_returns_error() {
            let t_obs = 59000.0;
            let (ds, cache) = make_obs_and_cache(t_obs);
            let hyperbolic = EquinoctialElements {
                reference_epoch: t_obs,
                semi_major_axis: 1.0,
                eccentricity_sin_lon: 0.8,
                eccentricity_cos_lon: 0.8, // e ≈ 1.13 > 1
                tan_half_incl_sin_node: 0.0,
                tan_half_incl_cos_node: 0.0,
                mean_longitude: 0.0,
            };
            let obs = ds.get_observation(0).unwrap();
            assert!(obs
                .compute_apparent_pos_derivative(&cache, &JPL_EPHEM_HORIZON, &hyperbolic)
                .is_err());
        }

        /// Non-regression oracle test for `topocentric_radec_and_partials`.
        ///
        /// Reference values were computed once in Rust (IEEE 754 double precision)
        /// and checked in as ground truth.  The test guards against accidental
        /// regressions in the Jacobian formula.
        ///
        /// Inputs
        /// ------
        /// ast_pos_ecl = [-4.15006894757462969e-2,  1.36234947972925746e0,  8.71837102048550472e-1]  AU
        /// ast_vel_ecl = [-1.41847650703540683e-2,  1.14549098612109205e-4,  4.07819544228800175e-4]  AU/day
        /// obs_pos_ecl = [0, 0, 0]
        ///
        /// Expected outputs (oracle)
        /// ------------------------
        /// ra  = 1.60115231410015513e0   rad
        /// dec = 5.69067704034184052e-1  rad
        /// d_ra_d_pos  = [-7.33348893215877928e-1, -2.23189920137348632e-2, -3.23656066883762261e-5]  rad/AU
        /// d_dec_d_pos = [ 1.01082321534895769e-2, -3.32887445954391237e-1,  5.20657513131207450e-1]  rad/AU
        #[test]
        fn oracle_topocentric_radec_and_partials() {
            use crate::constants::ROT_EQUMJ2000_TO_ECLMJ2000;
            use approx::assert_abs_diff_eq;

            let t_obs = 59000.0;
            let elem = nominal_elements(t_obs);
            let (cart_pos_ast, cart_vel_ast, _) =
                elem.solve_two_body_problem(0., 0., false).unwrap();
            let obs_pos_ecl = Vector3::zeros();
            let ast_pos_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_ast;
            let ast_vel_ecl = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_vel_ast;

            let (ra, dec, d_ra, d_dec) =
                topocentric_radec_and_partials(ast_pos_ecl, ast_vel_ecl, obs_pos_ecl);

            // Oracle values (Rust IEEE 754 double precision reference)
            let tol = 1e-10_f64;

            assert_abs_diff_eq!(ra, 1.601_152_314_100_155_1, epsilon = tol);
            assert_abs_diff_eq!(dec, 5.690_677_040_341_84e-1, epsilon = tol);

            let ref_d_ra = [
                -7.333_488_932_158_779e-1,
                -2.231_899_201_373_486_3e-2,
                -3.236_560_668_837_622_6e-5,
            ];
            let ref_d_dec = [
                1.010_823_215_348_957_7e-2,
                -3.328_874_459_543_912_3e-1,
                5.206_575_131_312_075e-1,
            ];

            for i in 0..3 {
                assert_abs_diff_eq!(d_ra[i], ref_d_ra[i], epsilon = tol);
                assert_abs_diff_eq!(d_dec[i], ref_d_dec[i], epsilon = tol);
            }
        }
    }
}
