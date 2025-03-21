use crate::{
    constants::{Degree, MpcCode, ObjectNumber, Observations, TrajectorySet, MJD},
    time::frac_date_to_mjd,
};
use camino::Utf8Path;
use serde::Deserialize;
use serde_with::DeserializeFromStr;
use std::{collections::HashMap, ops::Range, str::FromStr};
use thiserror::Error;

/// Parse a right ascension string to degrees
///
/// Arguments
/// ---------
/// * `ra`: a string representing the right ascension in the format HH MM SS.SS
///
/// Return
/// ------
/// * a float representing the right ascension in degrees
fn parse_ra_to_deg(ra: &str) -> Option<Degree> {
    let parts: Vec<&str> = ra.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }
    let h: f64 = parts[0].parse().ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s: f64 = parts[2].parse().ok()?;
    Some((h + m / 60.0 + s / 3600.0) * 15.0)
}

/// Parse a declination string to degrees
///
/// Arguments
/// ---------
/// * `dec`: a string representing the declination in the format HH MM SS.SS
///
/// Return
/// ------
/// * a float representing the declination in degrees
fn parse_dec_to_deg(dec: &str) -> Option<Degree> {
    let parts: Vec<&str> = dec.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }
    let sign = if parts[0].starts_with('-') { -1.0 } else { 1.0 };
    let d: f64 = parts[0]
        .trim_start_matches(['-', '+'].as_ref())
        .parse()
        .ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s: f64 = parts[2].parse().ok()?;
    Some(sign * (d + m / 60.0 + s / 3600.0))
}

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

#[derive(Error, Debug)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
    #[error("The line is not a CCD observation")]
    NotCCDObs,
}

#[derive(Debug, DeserializeFromStr)]
struct ObservationWrapper(Observation);

impl FromStr for ObservationWrapper {
    type Err = ParseObsError;

    /// Parse a line from an 80 column file to an Observation
    ///
    /// Arguments
    /// ---------
    /// * `line`: a string representing a line from an 80 column file
    ///
    /// Return
    /// ------
    /// * an Observation struct
    fn from_str(line: &str) -> Result<Self, Self::Err> {
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
            observer: line[77..80].trim().to_string(),
        };

        Ok(ObservationWrapper(observation))
    }
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
fn extract_80col(colfile: &Utf8Path) -> (Observations, ObjectNumber) {
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
fn observation_from_vec(
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

/// A trait to extend the TrajectorySet struct
/// It provides methods to create a TrajectorySet from an 80 column file and to add observations from an 80 column file
/// It also provides methods to create a TrajectorySet from vectors of right ascension, declination, time, and observer and to add observations from vectors of right ascension, declination, time, and observer
pub trait TrajectoryExt {
    fn from_80col(colfile: &Utf8Path) -> Self;
    fn add_80col(&mut self, colfile: &Utf8Path);

    fn new_from_vec(
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) -> Self;
    fn add_from_vec(
        &mut self,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    );
}

impl TrajectoryExt for TrajectorySet {
    /// Create a TrajectorySet from an object number, right ascension, declination, time, and one observer.
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `dec`: a vector of declination
    /// * `time`: a vector of time
    /// * `observer`: the observer
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations
    fn new_from_vec(
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) -> Self {
        let observations: Observations = observation_from_vec(ra, dec, time, observer);
        let mut traj_set: HashMap<ObjectNumber, Observations> = HashMap::new();
        traj_set.insert(object_number.to_string(), observations);
        traj_set
    }

    /// Add the observations of an object number, right ascension, declination, time, and one observer to a TrajectorySet
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `dec`: a vector of declination
    /// * `time`: a vector of time
    /// * `observer`: the observer
    fn add_from_vec(
        &mut self,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) {
        let observations: Observations = observation_from_vec(ra, dec, time, observer);
        self.insert(object_number.to_string(), observations);
    }

    /// Create a TrajectorySet from an 80 column file
    ///
    /// Arguments
    /// ---------
    /// * `colfile`: a path to an 80 column file
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations from the 80 column file
    fn from_80col(colfile: &Utf8Path) -> Self {
        let mut traj_set: HashMap<ObjectNumber, Observations> = HashMap::new();
        let (observations, object_number) = extract_80col(colfile);
        traj_set.insert(object_number, observations);
        traj_set
    }

    /// Add the observations of a 80 column file to a TrajectorySet
    ///
    /// Arguments
    /// ---------
    /// * `colfile`: a path to an 80 column file
    fn add_80col(&mut self, colfile: &Utf8Path) {
        let (observations, object_number) = extract_80col(colfile);
        self.insert(object_number, observations);
    }
}

#[cfg(test)]
mod observations_test {
    use super::*;

    #[test]
    fn test_ra_to_deg() {
        assert_eq!(parse_ra_to_deg("23 58 57.68"), Some(359.7403333333333));
        assert_eq!(parse_ra_to_deg("04 41 04.77"), Some(70.269875));
        assert_eq!(parse_ra_to_deg("1 2 3.4.5"), None);
        assert_eq!(parse_ra_to_deg("1 2"), None);
    }

    #[test]
    fn test_dec_to_deg() {
        assert_eq!(parse_dec_to_deg("-00 30 14.2"), Some(-0.5039444444444444));
        assert_eq!(parse_dec_to_deg("+13 55 42.7"), Some(13.928527777777777));
        assert_eq!(parse_dec_to_deg("89 15 50.2"), Some(89.26394444444445));
        assert_eq!(parse_dec_to_deg("89 15 50.2.3"), None);
        assert_eq!(parse_dec_to_deg("89 15"), None);
    }
}
