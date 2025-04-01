use crate::{
    constants::{Degree, ObjectNumber, Observations, MJD},
    conversion::{parse_dec_to_deg, parse_ra_to_deg},
    observers::observers::Observer,
    outfit::Outfit,
    time::frac_date_to_mjd,
};
use camino::Utf8Path;
use std::{ops::Range, sync::Arc};
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
    pub dec: Degree,
    pub time: MJD,
}

impl Observation {
    pub fn new(observer: u16, ra: Degree, dec: Degree, time: MJD) -> Self {
        Observation {
            observer,
            ra,
            dec,
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

    let observation = Observation {
        time: frac_date_to_mjd(line[15..32].trim())
            .expect(format!("Error parsing date: {}", line[15..32].trim()).as_str()),
        ra: parse_ra_to_deg(line[32..44].trim())
            .expect(format!("Error parsing RA: {}", line[32..44].trim()).as_str()),
        dec: parse_dec_to_deg(line[44..56].trim())
            .expect(format!("Error parsing DEC: {}", line[44..56].trim()).as_str()),
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
    ra.iter()
        .zip(dec.iter())
        .zip(time.iter())
        .map(|((ra, dec), time)| Observation {
            ra: *ra,
            dec: *dec,
            time: *time,
            observer: obs_uin16,
        })
        .collect()
}
