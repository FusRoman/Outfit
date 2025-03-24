use crate::{constants::{Degree, MpcCode, ObjectNumber, Observations, MJD}, observations::observation_wrapper::{ObservationWrapper, ParseObsError}};
use camino::Utf8Path;
use serde::Deserialize;
use std::ops::Range;

/// A struct containing the observer and the time, ra, dec, and rms of an observation
///
/// # Fields
///
/// * `observer` - The observer
/// * `ra` - The right ascension of the observation
/// * `dec` - The declination of the observation
/// * `time` - The time of the observation
#[derive(Debug, Deserialize)]
pub struct Observation {
    pub observer: MpcCode,
    pub ra: Degree,
    pub dec: Degree,
    pub time: MJD,
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
pub(crate) fn extract_80col(colfile: &Utf8Path) -> (Observations, ObjectNumber) {
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
            .filter_map(|line| match line.parse::<ObservationWrapper>() {
                Ok(wrapper) => Some(wrapper.0),
                Err(ParseObsError::NotCCDObs) => None,
                Err(e) => panic!("Error parsing line: {:?}", e),
            })
            .collect(),
        object_number,
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
    ra: &Vec<Degree>,
    dec: &Vec<Degree>,
    time: &Vec<MJD>,
    observer: &str,
) -> Observations {
    ra.iter()
        .zip(dec.iter())
        .zip(time.iter())
        .map(|((ra, dec), time)| Observation {
            ra: *ra,
            dec: *dec,
            time: *time,
            observer: observer.to_string(),
        })
        .collect()
}
