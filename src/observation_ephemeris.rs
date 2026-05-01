use std::f64::consts::PI;

use hifitime::Epoch;
use nalgebra::Vector3;
use photom::{
    constants::DPI, coordinates::equatorial::EquCoord,
    observation_dataset::observation::Observation,
};

use crate::{
    cache::OutfitCache, constants::ROT_EQUMJ2000_TO_ECLMJ2000, conversion::cartesion_from_vec,
    EquinoctialElements, JPLEphem, OutfitError, VLIGHT_AU,
};

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
    /// * [`Observer::pvobs`] – Computes observer’s geocentric position.
    /// * [`correct_aberration`] – Aberration correction.
    /// * [`cartesian_to_radec`] – Convert Cartesian vectors to (RA, DEC).
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
    /// * [`compute_apparent_position`](crate::observations::Observation::compute_apparent_position) – Used internally to obtain predicted RA/DEC.
    /// * [`Observer::pvobs`] – Computes observer’s geocentric position.
    /// * [`correct_aberration`] – Applies aberration correction.
    /// * [`cartesian_to_radec`] – Converts 3D vectors to (RA, DEC).
    /// * [`EquinoctialElements::solve_two_body_problem`] – Two-body propagation.
    fn ephemeris_error(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        equinoctial_element: &EquinoctialElements,
    ) -> Result<f64, OutfitError>;
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
        let cart_pos_ast_eclj2000 = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_ast;
        let cart_pos_vel_eclj2000 = ROT_EQUMJ2000_TO_ECLMJ2000 * cart_pos_vel;

        // 4. Observer heliocentric position
        let geo_obs_pos = cache
            .get_observer_geocentric_position(self.index())
            .map(|x| x.into_inner());
        let xobs = geo_obs_pos + earth_pos_eclj2000;
        let obs_on_earth = ROT_EQUMJ2000_TO_ECLMJ2000 * xobs;

        // 5. Relative position and aberration correction
        let relative_position = cart_pos_ast_eclj2000 - obs_on_earth;
        let corrected_pos = correct_aberration(relative_position, cart_pos_vel_eclj2000);
        let cartesion_pos = cartesion_from_vec(corrected_pos);

        // 6. Convert to equatorial coordinates
        let equatorial_pos: EquCoord = cartesion_pos.into();

        Ok((equatorial_pos.ra, equatorial_pos.dec))
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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                    OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                    OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();
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
                OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();
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
}
