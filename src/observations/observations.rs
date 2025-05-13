use crate::{
    constants::{Degree, ObjectNumber, Observations, DPI, MJD},
    conversion::{parse_dec_to_deg, parse_ra_to_deg},
    equinoctial_element::EquinoctialElements,
    kepler::principal_angle,
    observers::{
        observer_position::{geo_obs_pos, pvobs},
        observers::Observer,
    },
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
/// * `observer` - The observer
/// * `ra` - The right ascension of the observation
/// * `dec` - The declination of the observation
/// * `time` - The time of the observation
#[derive(Debug)]
pub struct Observation {
    observer: u16,
    pub ra: Degree,
    pub error_ra: Degree,
    pub dec: Degree,
    pub error_dec: Degree,
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
    pub fn new(observer: u16, ra: Degree, dec: Degree, time: MJD) -> Self {
        Observation {
            observer,
            ra,
            error_ra: 0.,
            dec,
            error_dec: 0.,
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

    fn ephemeris_error(
        &self,
        state: &Outfit,
        equinoctial_element: &EquinoctialElements,
        observer: &Observer,
    ) -> Result<f64, OutfitError> {
        let (cart_pos_ast, cart_pos_vel, _) = equinoctial_element.solve_two_body_problem(
            0.,
            self.time - equinoctial_element.reference_epoch,
            false,
        )?;

        let obs_mjd = Epoch::from_mjd_in_time_scale(self.time, hifitime::TimeScale::TT);

        let (earth_position, _) = state
            .get_jpl_ephem()?
            .earth_ephemeris(&obs_mjd, true);

        let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
        rotpn(&mut roteqec, "EQUM", "J2000", 0., "ECLM", "J2000", 0.);
        let matrix_elc_transform = Matrix3::from(roteqec);

        let earth_pos_eclj2000 = matrix_elc_transform.transpose() * earth_position;

        let cart_pos_ast_eclj2000 = matrix_elc_transform * cart_pos_ast;
        let cart_pos_vel_eclj2000 = matrix_elc_transform * cart_pos_vel;

        let (geo_obs_pos, _) = geo_obs_pos(&observer, &obs_mjd, &state.get_ut1_provider());

        let xobs = geo_obs_pos + earth_pos_eclj2000;

        let obs_on_earth = matrix_elc_transform * xobs;

        let relative_position = cart_pos_ast_eclj2000 - obs_on_earth;

        let corrected_pos = correct_aberration(relative_position, cart_pos_vel_eclj2000);
        let (alpha, delta, _) = cartesian_to_radec(corrected_pos);

        // -----
        // floating point error between orbfit anf outfit is around 1e-5 in absolute and 1e-7 in relative
        // this is something to check in the future for reproducing the same results
        // -----
        let diff_alpha = (self.ra - alpha) % DPI;
        let diff_alpha = if diff_alpha > PI {
            diff_alpha - DPI
        } else {
            diff_alpha
        };

        let diff_delta = self.dec - delta;
        let rms_ra = (self.dec.cos() * (diff_alpha / self.error_ra)).powi(2);
        let rms_dec = (diff_delta / self.error_dec).powi(2);
        Ok(rms_ra + rms_dec)
    }
}

#[derive(Error, Debug)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
    #[error("The line is not a CCD observation")]
    NotCCDObs,
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
fn from_80col(env_state: &mut Outfit, line: &str) -> Result<Observation, ParseObsError> {
    if line.len() < 80 {
        return Err(ParseObsError::TooShortLine);
    }

    if line.chars().nth(14) == Some('s') {
        return Err(ParseObsError::NotCCDObs);
    }

    // TODO: replace error_ra and error_dec with the correct values
    let observation = Observation {
        time: frac_date_to_mjd(line[15..32].trim())
            .expect(format!("Error parsing date: {}", line[15..32].trim()).as_str()),
        ra: parse_ra_to_deg(line[32..44].trim())
            .expect(format!("Error parsing RA: {}", line[32..44].trim()).as_str()),
        error_ra: 0.,
        dec: parse_dec_to_deg(line[44..56].trim())
            .expect(format!("Error parsing DEC: {}", line[44..56].trim()).as_str()),
        error_dec: 0.,
        observer: env_state.uint16_from_mpc_code(&line[77..80].trim().into()),
    };

    Ok(observation)
}

/// Extract the observations and the object number from an 80 column file
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
) -> (Observations, ObjectNumber) {
    let file_content = std::fs::read_to_string(colfile)
        .expect(format!("Could not read file {}", colfile.as_str()).as_str());

    let first_line = file_content
        .lines()
        .next()
        .expect(format!("Could not read first line of file {}", colfile.as_str()).as_str());

    fn get_object_number(line: &str, range: Range<usize>) -> String {
        line[range].trim_start_matches('0').trim().to_string()
    }

    let mut object_number = get_object_number(first_line, 0..5);
    if object_number.is_empty() {
        object_number = get_object_number(first_line, 5..12);
    }

    (
        file_content
            .lines()
            .filter_map(|line| match from_80col(env_state, line) {
                Ok(obs) => Some(obs),
                Err(ParseObsError::NotCCDObs) => None,
                Err(e) => panic!("Error parsing line: {:?}", e),
            })
            .collect(),
        ObjectNumber::String(object_number),
    )
}

/// Create a vector of Observations from vectors of right ascension, declination, time, and observer
/// Each observations should have been observed by the same observer.
///
/// Arguments
/// ---------
/// * `ra`: a vector of right ascension
/// * `dec`: a vector of declination
/// * `time`: a vector of time
/// * `observer`: the observer
///
/// Return
/// ------
/// * a vector of Observations
pub(crate) fn observation_from_vec(
    env_state: &mut Outfit,
    ra: &Vec<Degree>,
    dec: &Vec<Degree>,
    time: &Vec<MJD>,
    observer: Arc<Observer>,
) -> Observations {
    let obs_uin16 = env_state.uint16_from_observer(observer);
    // TODO: replace error_ra and error_dec with the correct values
    ra.iter()
        .zip(dec.iter())
        .zip(time.iter())
        .map(|((ra, dec), time)| Observation {
            ra: *ra,
            error_ra: 0.,
            dec: *dec,
            error_dec: 0.,
            time: *time,
            observer: obs_uin16,
        })
        .collect()
}

#[cfg(test)]
mod test_observations {
    use super::*;

    use crate::unit_test_global::OUTFIT_HORIZON_TEST;

    #[test]
    fn test_ephem_error() {
        let obs = Observation {
            observer: 0,
            ra: 1.7899347771316527,
            error_ra: 1.7700245206085460E-006,
            dec: 0.77899653810797365,
            error_dec: 1.2595828918293177E-006,
            time: 57070.262067592594,
        };

        let observer = OUTFIT_HORIZON_TEST.get_observer_from_mpc_code(&"F51".to_string());

        let equinoctial_element = EquinoctialElements {
            reference_epoch: 57049.242334573748,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.26937368090922720,
            eccentricity_cos_lon: 8.8564152600135601E-002,
            tan_half_incl_sin_node: 8.0899701663963020E-004,
            tan_half_incl_cos_node: 0.10168201109730375,
            mean_longitude: 1.6936970079414786,
        };

        let rms_error = obs.ephemeris_error(&OUTFIT_HORIZON_TEST, &equinoctial_element, &observer);
        assert_eq!(rms_error.unwrap(), 75.00445641224026);
    }
}
