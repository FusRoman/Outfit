pub mod ades_reader;
pub mod observations_ext;
mod parquet_reader;
pub mod trajectory_ext;
pub(crate) mod triplets_iod;

use crate::{
    constants::{ObjectNumber, Observations, Radian, DPI, MJD, RADH, RADSEC},
    conversion::{parse_dec_to_deg, parse_ra_to_deg},
    equinoctial_element::EquinoctialElements,
    observations::trajectory_ext::ObservationBatch,
    observers::{observer_position::geo_obs_pos, Observer},
    outfit::Outfit,
    outfit_errors::OutfitError,
    ref_system::{cartesian_to_radec, correct_aberration, rotpn},
    time::frac_date_to_mjd,
};
use camino::Utf8Path;
use hifitime::Epoch;
use nalgebra::Matrix3;
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
    /// * [`geo_obs_pos`] – Computes geocentric position of the observer including Earth rotation and nutation.
    /// * [`correct_aberration`] – Applies aberration correction to the apparent direction of the body.
    /// * [`cartesian_to_radec`] – Converts 3D Cartesian vectors into equatorial coordinates (RA/DEC).
    /// * [`solve_two_body_problem`] – Computes heliocentric position and velocity from orbital elements.
    pub(crate) fn ephemeris_error(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
        observer: &Observer,
    ) -> Result<f64, OutfitError> {
        // Propagate heliocentric position and velocity of the asteroid using a two-body model.
        // Output is in the ecliptic mean J2000 reference frame.
        let (cart_pos_ast, cart_pos_vel, _) = equinoctial_element.solve_two_body_problem(
            0.,
            self.time - equinoctial_element.reference_epoch,
            false,
        )?;

        // Convert observation time to TT (Terrestrial Time) for ephemeris consistency
        let obs_mjd = Epoch::from_mjd_in_time_scale(self.time, hifitime::TimeScale::TT);

        // Get Earth's heliocentric position (barycentric if true is set) at observation epoch
        let (earth_position, _) = state.get_jpl_ephem()?.earth_ephemeris(&obs_mjd, true);

        // Construct rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
        rotpn(&mut roteqec, "EQUM", "J2000", 0., "ECLM", "J2000", 0.);
        let matrix_elc_transform = Matrix3::from(roteqec);

        // Transform Earth's position to ecliptic J2000 frame
        let earth_pos_eclj2000 = matrix_elc_transform.transpose() * earth_position;

        // Transform asteroid's propagated position and velocity into ecliptic J2000 frame
        let cart_pos_ast_eclj2000 = matrix_elc_transform * cart_pos_ast;
        let cart_pos_vel_eclj2000 = matrix_elc_transform * cart_pos_vel;

        // Compute the observer's geocentric position (accounts for Earth rotation, nutation, etc.)
        let (geo_obs_pos, _) = geo_obs_pos(observer, &obs_mjd, state.get_ut1_provider());

        // Observer's heliocentric position = Earth position + observer's geocentric offset
        let xobs = geo_obs_pos + earth_pos_eclj2000;

        // Transform observer's position to ecliptic J2000 frame
        let obs_on_earth = matrix_elc_transform * xobs;

        // Compute the vector from observer to asteroid in heliocentric frame
        let relative_position = cart_pos_ast_eclj2000 - obs_on_earth;

        // Apply aberration correction due to observer motion (velocity of asteroid required)
        let corrected_pos = correct_aberration(relative_position, cart_pos_vel_eclj2000);

        // Convert corrected position to equatorial coordinates (RA, DEC)
        let (alpha, delta, _) = cartesian_to_radec(corrected_pos);

        // Compute difference in right ascension (with wrapping modulo 2π)
        let diff_alpha = (self.ra - alpha) % DPI;
        let diff_alpha = if diff_alpha > PI {
            diff_alpha - DPI
        } else {
            diff_alpha
        };

        // Compute difference in declination
        let diff_delta = self.dec - delta;

        // Normalize residuals by astrometric uncertainties.
        // RA is scaled by cos(DEC) to account for spherical projection
        let rms_ra = (self.dec.cos() * (diff_alpha / self.error_ra)).powi(2);
        let rms_dec = (diff_delta / self.error_dec).powi(2);

        // Return the total normalized squared residuals (RA² + DEC²)
        Ok(rms_ra + rms_dec)
    }
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
                Err(e) => panic!("Error parsing line: {:?}", e),
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
            error_ra: 1.7700245206085460E-006,
            dec: 0.77899653810797365,
            error_dec: 1.2595828918293177E-006,
            time: 57070.262067592594,
        };

        let observer = OUTFIT_HORIZON_TEST
            .0
            .get_observer_from_mpc_code(&"F51".to_string());

        let equinoctial_element = EquinoctialElements {
            reference_epoch: 57049.242334573748,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.26937368090922720,
            eccentricity_cos_lon: 8.8564152600135601E-002,
            tan_half_incl_sin_node: 8.0899701663963020E-004,
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
