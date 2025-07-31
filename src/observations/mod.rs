pub mod ades_reader;
pub mod observations_ext;
mod parquet_reader;
pub mod trajectory_ext;
pub(crate) mod triplets_iod;

use crate::{
    constants::{ObjectNumber, Observations, Radian, DPI, MJD, RADH, RADSEC, VLIGHT_AU},
    conversion::{cartesian_to_radec, parse_dec_to_deg, parse_ra_to_deg},
    equinoctial_element::EquinoctialElements,
    observations::trajectory_ext::ObservationBatch,
    observers::Observer,
    outfit::Outfit,
    outfit_errors::OutfitError,
    ref_system::{rotpn, RefEpoch, RefSystem},
    time::frac_date_to_mjd,
};
use camino::Utf8Path;
use hifitime::Epoch;
use nalgebra::{Matrix3, Vector3};
use std::{f64::consts::PI, ops::Range, sync::Arc};
use thiserror::Error;

/// A struct containing the observer and the time, ra, dec, and rms of an observation
///
/// # Fields
///
/// * `observer` - The observer index (u16), unique identifier for the observer
/// * `ra` - The right ascension of the observation in Radians
/// * `error_ra` - The error in right ascension in Radians
/// * `dec` - The declination of the observation in Radians
/// * `error_dec` - The error in declination in Radians
/// * `time` - The time of the observation
#[derive(Debug, Clone, PartialEq)]
pub struct Observation {
    pub(crate) observer: u16,
    pub ra: Radian,
    pub error_ra: Radian,
    pub dec: Radian,
    pub error_dec: Radian,
    pub time: MJD,
}

impl Observation {
    /// Create a new observation
    ///
    /// Arguments
    /// ---------
    /// * `observer`: the observer
    /// * `ra`: the right ascension of the observation
    /// * `dec`: the declination of the observation
    /// * `time`: the time of the observation
    ///
    /// Return
    /// ------
    /// * a new Observation struct
    ///   The coordinates and the associated errors are converted from degrees and arcseconds respectively to radians.
    pub fn new(
        observer: u16,
        ra: Radian,
        error_ra: Radian,
        dec: Radian,
        error_dec: Radian,
        time: MJD,
    ) -> Self {
        Observation {
            observer,
            ra,
            error_ra,
            dec,
            error_dec,
            time,
        }
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
    /// as seen by a specific observer at the observation time.
    ///
    /// # Overview
    ///
    /// This function determines where an object (described by its equinoctial orbital elements)
    /// will appear on the sky at the time of the current observation.
    ///
    /// The computation involves:
    ///
    /// 1. **Orbit propagation**:
    ///    Propagates the object's position and velocity from its reference epoch to the
    ///    observation time using a two-body model.
    ///
    /// 2. **Reference frame handling**:
    ///    Retrieves Earth's barycentric position from a JPL ephemeris,
    ///    and builds a rotation matrix to express positions in the *ecliptic mean J2000* frame.
    ///
    /// 3. **Observer position**:
    ///    Computes the geocentric position of the observer, adds Earth’s heliocentric position
    ///    to obtain the full heliocentric position of the observer, and transforms it into the
    ///    same frame.
    ///
    /// 4. **Light-time and aberration corrections**:
    ///    Forms the vector from the observer to the object and applies a simple aberration
    ///    correction using the relative velocity of the object.
    ///
    /// 5. **Conversion to equatorial coordinates**:
    ///    Converts the corrected line-of-sight vector into apparent right ascension (RA)
    ///    and declination (DEC).
    ///
    /// # Returns
    ///
    /// * `(alpha, delta)` – A tuple containing the **apparent right ascension** and
    ///   **declination**, in **radians**.
    ///
    /// # Units
    ///
    /// * Positions: astronomical units (AU)
    /// * Velocities: AU/day
    /// * Angles: radians
    /// * Time: Modified Julian Date (MJD)
    ///
    /// # Errors
    ///
    /// Returns an `OutfitError` if:
    /// - Propagation of the orbit fails,
    /// - JPL ephemeris data is unavailable,
    /// - Reference frame transformation fails.
    ///
    /// # See also
    ///
    /// * [`solve_two_body_problem`] – Orbit propagation.
    /// * [`pvobs`] – Computes observer geocentric position including Earth rotation.
    /// * [`correct_aberration`] – Corrects apparent direction due to observer motion.
    /// * [`cartesian_to_radec`] – Converts a Cartesian vector to (RA, DEC).
    pub fn compute_apparent_position(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
        observer: &Observer,
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
        let (earth_position, _) = state.get_jpl_ephem()?.earth_ephemeris(&obs_mjd, true);

        // 4. Build rotation from equatorial mean J2000 to ecliptic mean J2000
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let matrix_elc_transform = Matrix3::from(rotpn(&ref_sys1, &ref_sys2)?);

        let earth_pos_eclj2000 = matrix_elc_transform.transpose() * earth_position;
        let cart_pos_ast_eclj2000 = matrix_elc_transform * cart_pos_ast;
        let cart_pos_vel_eclj2000 = matrix_elc_transform * cart_pos_vel;

        // 5. Observer heliocentric position
        let (geo_obs_pos, _) = observer.pvobs(&obs_mjd, state.get_ut1_provider())?;
        let xobs = geo_obs_pos + earth_pos_eclj2000;
        let obs_on_earth = matrix_elc_transform * xobs;

        // 6. Relative position and aberration correction
        let relative_position = cart_pos_ast_eclj2000 - obs_on_earth;
        let corrected_pos = correct_aberration(relative_position, cart_pos_vel_eclj2000);

        // 7. Convert to RA/DEC
        let (alpha, delta, _) = cartesian_to_radec(corrected_pos);
        Ok((alpha, delta))
    }

    /// Compute the normalized squared astrometric residuals (RA/DEC) between an observed position and a propagated ephemeris.
    ///
    /// This function compares the actual observed position (stored in `self`) to the expected astrometric position
    /// derived from the propagation of equinoctial orbital elements using a two-body model, corrected for light aberration
    /// and transformed into the appropriate equatorial reference frame. It returns a scalar value representing the sum of
    /// squared, normalized residuals in RA and DEC.
    ///
    /// # Arguments
    /// ---------------
    /// * `state` - The current application state providing ephemerides and time conversions.
    /// * `equinoctial_element` - The orbital elements of the object used to propagate its position.
    /// * `observer` - The observer's geographic location and observing context.
    ///
    /// # Returns
    /// ----------
    /// * `Result<f64, OutfitError>` - A scalar value representing the weighted sum of squared residuals (dimensionless).
    ///   This value is roughly equivalent to a reduced chi-squared contribution for a single observation (but not divided by 2).
    ///
    /// # Notes
    /// ----------
    /// - The residuals are normalized using the observation's astrometric uncertainties (`error_ra`, `error_dec`).
    /// - Right ascension differences are corrected by `cos(dec)` to account for projection effects.
    /// - The resulting value is **dimensionless**.
    /// - Units: all angles are in **radians**.
    ///
    /// # Errors
    /// ----------
    /// Returns an `OutfitError` if ephemeris data or propagation fails (e.g. missing JPL ephemeris).
    ///
    /// # See also
    /// * [`pvobs`] – Computes geocentric position of the observer including Earth rotation and nutation.
    /// * [`correct_aberration`] – Applies aberration correction to the apparent direction of the body.
    /// * [`cartesian_to_radec`] – Converts 3D Cartesian vectors into equatorial coordinates (RA/DEC).
    /// * [`solve_two_body_problem`] – Computes heliocentric position and velocity from orbital elements.
    pub fn ephemeris_error(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
        observer: &Observer,
    ) -> Result<f64, OutfitError> {
        let (alpha, delta) =
            self.compute_apparent_position(state, equinoctial_element, observer)?;

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
/// * `xrel`: relative position vector from observer to object [AU].
/// * `vrel`: velocity of the observer relative to the barycenter [AU/day].
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
fn correct_aberration(xrel: Vector3<f64>, vrel: Vector3<f64>) -> Vector3<f64> {
    let norm_vector = xrel.norm();
    let dt = norm_vector / VLIGHT_AU;
    xrel - dt * vrel
}

#[derive(Error, Debug, PartialEq)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
    #[error("The line is not a CCD observation")]
    NotCCDObs,
    #[error("Error parsing RA: {0}")]
    InvalidRA(String),
    #[error("Invalid Dec value: {0}")]
    InvalidDec(String),
    #[error("Invalid date: {0}")]
    InvalidDate(String),
}

/// Parse a line from an 80 column file to an Observation
///
/// Arguments
/// ---------
/// * `line`: a string representing a line from an 80 column file
///
/// Return
/// ------
/// * an Observation struct
fn from_80col(env_state: &mut Outfit, line: &str) -> Result<Observation, OutfitError> {
    if line.len() < 80 {
        return Err(OutfitError::Parsing80ColumnFileError(
            ParseObsError::TooShortLine,
        ));
    }

    if line.chars().nth(14) == Some('s') {
        return Err(OutfitError::Parsing80ColumnFileError(
            ParseObsError::NotCCDObs,
        ));
    }

    let (ra, error_ra) = parse_ra_to_deg(line[32..44].trim()).ok_or_else(|| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidRA(
            line[32..44].trim().to_string(),
        ))
    })?;

    let (dec, error_dec) = parse_dec_to_deg(line[44..56].trim()).ok_or_else(|| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidDec(
            line[44..56].trim().to_string(),
        ))
    })?;

    let time = frac_date_to_mjd(line[15..32].trim()).map_err(|_| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidDate(
            line[15..32].trim().to_string(),
        ))
    })?;

    let observer_id = env_state.uint16_from_mpc_code(&line[77..80].trim().into());
    let observer = env_state.get_observer_from_uint16(observer_id);

    let max_rms = |observation_error: f64, observer_error: f64, factor: f64| {
        f64::max(observation_error, observer_error * factor)
    };

    let dec_radians = dec.to_radians();
    let dec_rad_cos = dec_radians.cos();

    let ra_error = max_rms(
        (error_ra * RADH) / dec_rad_cos,
        observer.ra_accuracy.map(|v| v.into_inner()).unwrap_or(0.0),
        RADSEC / dec_rad_cos,
    );

    let dec_error = max_rms(
        error_dec.to_radians(),
        observer.dec_accuracy.map(|v| v.into_inner()).unwrap_or(0.0),
        RADSEC,
    );

    let observation = Observation::new(
        observer_id,
        ra.to_radians(),
        ra_error,
        dec_radians,
        dec_error,
        time,
    );
    Ok(observation)
}

/// Extract the observations and the object number from a 80 column file
///
/// Arguments
/// ---------
/// * `colfile`: a path to an 80 column file
///
/// Return
/// ------
/// * a tuple containing the observations and the object number
pub(crate) fn extract_80col(
    env_state: &mut Outfit,
    colfile: &Utf8Path,
) -> Result<(Observations, ObjectNumber), OutfitError> {
    let file_content = std::fs::read_to_string(colfile)
        .unwrap_or_else(|_| panic!("Could not read file {}", colfile.as_str()));

    let first_line = file_content
        .lines()
        .next()
        .unwrap_or_else(|| panic!("Could not read first line of file {}", colfile.as_str()));

    fn get_object_number(line: &str, range: Range<usize>) -> String {
        line[range].trim_start_matches('0').trim().to_string()
    }

    let mut object_number = get_object_number(first_line, 0..5);
    if object_number.is_empty() {
        object_number = get_object_number(first_line, 5..12);
    }

    Ok((
        file_content
            .lines()
            .filter_map(|line| match from_80col(env_state, line) {
                Ok(obs) => Some(obs),
                Err(OutfitError::Parsing80ColumnFileError(ParseObsError::NotCCDObs)) => None,
                Err(e) => panic!("Error parsing line: {e:?}"),
            })
            .collect(),
        ObjectNumber::String(object_number),
    ))
}

/// Create a vector of Observations from a batch of right ascension, declination, and time values.
///
/// Each observation in the batch is assumed to come from the same observer.
///
/// # Arguments
/// * `env_state`: global Outfit state
/// * `batch`: RA/DEC/time values with their corresponding uncertainties
/// * `observer`: the observer
///
/// # Returns
/// A vector of [`Observation`] corresponding to the given inputs.
pub(crate) fn observation_from_batch(
    env_state: &mut Outfit,
    batch: &ObservationBatch<'_>,
    observer: Arc<Observer>,
) -> Observations {
    let obs_uin16 = env_state.uint16_from_observer(observer);

    batch
        .ra
        .iter()
        .zip(batch.dec.iter())
        .zip(batch.time.iter())
        .map(|((ra, dec), time)| {
            Observation::new(obs_uin16, *ra, batch.error_ra, *dec, batch.error_dec, *time)
        })
        .collect()
}

#[cfg(test)]
mod test_observations {

    use super::*;

    #[cfg(test)]
    #[cfg(feature = "jpl-download")]
    mod tests_compute_apparent_position {
        use super::*;
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;
        use approx::assert_relative_eq;

        /// Helper: construct a standard observer at Haleakalā Observatory.
        fn haleakala_observer() -> Observer {
            Observer::new(
                203.744083, // longitude in degrees east
                20.706944,  // latitude in degrees
                3.05,       // elevation in km
                Some("Haleakala".to_string()),
                None,
                None,
            )
        }

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
            let state = &OUTFIT_HORIZON_TEST.0;
            let observer = haleakala_observer();
            let t_obs = 59000.0; // MJD
            let equinoctial = simple_circular_elements(t_obs);

            let obs = Observation {
                observer: 0,
                ra: 0.0,
                error_ra: 0.0,
                dec: 0.0,
                error_dec: 0.0,
                time: t_obs,
            };

            let (ra, dec) = obs
                .compute_apparent_position(state, &equinoctial, &observer)
                .expect("Computation should succeed");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
            assert!((0.0..2.0 * std::f64::consts::PI).contains(&ra));
            assert!((-std::f64::consts::FRAC_PI_2..std::f64::consts::FRAC_PI_2).contains(&dec));
        }

        #[test]
        fn test_compute_apparent_position_same_epoch() {
            let state = &OUTFIT_HORIZON_TEST.0;
            let observer = haleakala_observer();
            let t_epoch = 60000.0;
            let equinoctial = simple_circular_elements(t_epoch);

            let obs = Observation {
                observer: 0,
                ra: 0.0,
                error_ra: 0.0,
                dec: 0.0,
                error_dec: 0.0,
                time: t_epoch,
            };

            let (ra1, dec1) = obs
                .compute_apparent_position(state, &equinoctial, &observer)
                .unwrap();
            let (ra2, dec2) = obs
                .compute_apparent_position(state, &equinoctial, &observer)
                .unwrap();

            // The same input should always produce the same result
            assert_relative_eq!(ra1, ra2, epsilon = 1e-14);
            assert_relative_eq!(dec1, dec2, epsilon = 1e-14);
        }

        #[test]
        fn test_apparent_position_for_distant_object() {
            let state = &OUTFIT_HORIZON_TEST.0;
            let observer = haleakala_observer();
            let t_obs = 59000.0;
            let mut equinoctial = simple_circular_elements(t_obs);

            // Objet très éloigné
            equinoctial.semi_major_axis = 100.0;

            let obs = Observation {
                observer: 0,
                ra: 0.0,
                error_ra: 0.0,
                dec: 0.0,
                error_dec: 0.0,
                time: t_obs,
            };

            let (ra, dec) = obs
                .compute_apparent_position(state, &equinoctial, &observer)
                .expect("Should compute apparent position for distant object");

            assert!(ra.is_finite());
            assert!(dec.is_finite());
        }

        #[test]
        fn test_compute_apparent_position_propagation_failure() {
            let state = &OUTFIT_HORIZON_TEST.0;
            let observer = haleakala_observer();

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

            let obs = Observation {
                observer: 0,
                ra: 0.0,
                error_ra: 0.0,
                dec: 0.0,
                error_dec: 0.0,
                time: 59000.0,
            };

            let result = obs.compute_apparent_position(state, &equinoctial, &observer);
            assert!(result.is_err(), "Invalid elements should trigger an error");
        }

        mod proptests_apparent_position {
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
                (-180.0..180.0, -90.0..90.0, 0.0..5.0)
                    .prop_map(|(lon, lat, elev)| Observer::new(lon, lat, elev, None, None, None))
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

            /// Returns a fixed observer located at Haleakalā Observatory.
            fn fixed_observer() -> Observer {
                Observer::new(
                    203.744083,
                    20.706944,
                    3.05,
                    Some("Haleakala".to_string()),
                    None,
                    None,
                )
            }

            proptest! {
                /// Property test: RA and DEC are always finite and in the expected ranges.
                #[test]
                fn proptest_ra_dec_are_finite_and_in_range(
                    equinoctial in arb_equinoctial_elements(),
                    obs_time in 58000.0f64..62000.0
                ) {
                    let state = &OUTFIT_HORIZON_TEST.0;
                    let observer = fixed_observer();

                    let obs = Observation {
                        observer: 0,
                        ra: 0.0,
                        error_ra: 0.0,
                        dec: 0.0,
                        error_dec: 0.0,
                        time: obs_time,
                    };

                    let result = obs.compute_apparent_position(state, &equinoctial, &observer);

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
                    let state = &OUTFIT_HORIZON_TEST.0;
                    let observer = fixed_observer();
                    let obs = Observation {
                        observer: 0,
                        ra: 0.0,
                        error_ra: 0.0,
                        dec: 0.0,
                        error_dec: 0.0,
                        time: obs_time,
                    };

                    let r1 = obs.compute_apparent_position(state, &equinoctial, &observer);
                    let r2 = obs.compute_apparent_position(state, &equinoctial, &observer);

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
                    let state = &OUTFIT_HORIZON_TEST.0;
                    let observer = fixed_observer();
                    let obs = Observation {
                        observer: 0,
                        ra: 0.0,
                        error_ra: 0.0,
                        dec: 0.0,
                        error_dec: 0.0,
                        time: obs_time,
                    };

                    let obs_eps = Observation {
                        time: obs_time + 1e-3, // shift by 1.4 minutes
                        ..obs.clone()
                    };

                    let r1 = obs.compute_apparent_position(state, &equinoctial, &observer);
                    let r2 = obs_eps.compute_apparent_position(state, &equinoctial, &observer);

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
                    let state = &OUTFIT_HORIZON_TEST.0;

                    let obs = Observation {
                        observer: 0,
                        ra: 0.0,
                        error_ra: 0.0,
                        dec: 0.0,
                        error_dec: 0.0,
                        time: obs_time,
                    };

                    let result = obs.compute_apparent_position(state, &equinoctial, &observer);

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
                let state = &OUTFIT_HORIZON_TEST.0;
                let observer = Observer::new(0.0, 0.0, 0.0, None, None, None);

                let equinoctial = EquinoctialElements {
                    reference_epoch: 59000.0,
                    semi_major_axis: 1.0,
                    eccentricity_sin_lon: 0.8,
                    eccentricity_cos_lon: 0.8, // e ≈ 1.13 > 1
                    tan_half_incl_sin_node: 0.0,
                    tan_half_incl_cos_node: 0.0,
                    mean_longitude: 0.0,
                };

                let obs = Observation {
                    observer: 0,
                    ra: 0.0,
                    error_ra: 0.0,
                    dec: 0.0,
                    error_dec: 0.0,
                    time: 59000.0,
                };

                let result = obs.compute_apparent_position(state, &equinoctial, &observer);
                assert!(
                    result.is_err(),
                    "Hyperbolic or parabolic orbits should currently return an error"
                );
            }
        }
    }

    #[test]
    fn test_new_observation() {
        let observation = Observation::new(1, 1.0, 0.1, 2.0, 0.2, 59000.0);
        assert_eq!(
            observation,
            Observation {
                observer: 1,
                ra: 1.0,
                error_ra: 0.1,
                dec: 2.0,
                error_dec: 0.2,
                time: 59000.0
            }
        );

        let observation_2 = Observation::new(
            2,
            343.097_375,
            2.777_777_777_777_778E-6,
            -14.784833333333333,
            2.777_777_777_777_778E-5,
            59001.0,
        );

        assert_eq!(
            observation_2,
            Observation {
                observer: 2,
                ra: 343.097375,
                error_ra: 2.777777777777778e-6,
                dec: -14.784833333333333,
                error_dec: 2.777777777777778e-5,
                time: 59001.0
            }
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_ephem_error() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let obs = Observation {
            observer: 0,
            ra: 1.7899347771316527,
            error_ra: 1.770_024_520_608_546E-6,
            dec: 0.778_996_538_107_973_6,
            error_dec: 1.259_582_891_829_317_7E-6,
            time: 57070.262067592594,
        };

        let observer = OUTFIT_HORIZON_TEST
            .0
            .get_observer_from_mpc_code(&"F51".to_string());

        let equinoctial_element = EquinoctialElements {
            reference_epoch: 57_049.242_334_573_75,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.269_373_680_909_227_2,
            eccentricity_cos_lon: 8.856_415_260_013_56E-2,
            tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
            tan_half_incl_cos_node: 0.10168201109730375,
            mean_longitude: 1.6936970079414786,
        };

        let rms_error =
            obs.ephemeris_error(&OUTFIT_HORIZON_TEST.0, &equinoctial_element, &observer);
        assert_eq!(rms_error.unwrap(), 75.00445641224026);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_from_80col_valid_line() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C2009 09 15.23433 22 52 22.62 -14 47 03.2          20.8 Vr~097wG96";

        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(result.is_ok());
        let obs = result.unwrap();

        assert_eq!(
            obs,
            Observation {
                observer: 0,
                ra: 5.988124307160555,
                error_ra: 1.2535340843609459e-6,
                dec: -0.25803335512429054,
                error_dec: 1.0181086985431635e-6,
                time: 55089.23509601851
            }
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_from_80col_too_short_line() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line = "short line";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::TooShortLine
            ))
        ));
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_from_80col_invalid_date() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C20xx 09 15.23433 22 52 22.62 -14 47 03.2          20.8 Vr~097wG96";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::InvalidDate(_)
            ))
        ));
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_from_80col_invalid_ra_dec() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C2009 09 15.23433 XX YY ZZ.ZZ -AA BB CC.C          20.8 Vr~097wG96";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::InvalidRA(_)
            ))
        ));
    }
}
