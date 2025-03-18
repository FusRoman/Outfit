use std::{collections::HashMap, hash::Hash, str::FromStr};
use crate::{constants::Degree, observers::observers::Observer};
use camino::Utf8Path;
use serde::Deserialize;
use serde_with::DeserializeFromStr;
use thiserror::Error;

/// A struct containing the observer and the time, ra, dec, and rms of an observation
/// 
/// # Fields
///
/// * `observer` - The observer
/// * `time` - The time of the observation
/// * `ra` - The right ascension of the observation
/// * `dec` - The declination of the observation
/// * `rms_ra` - The RMS of the right ascension
/// * `rms_dec` - The RMS of the declination
#[derive(Debug, Deserialize)]
pub struct Observation {
    observer: Observer,
    time: f64,
    ra: Degree,
    dec: Degree,
    rms_ra: Degree,
    rms_dec: Degree,
}

pub type Observations = Vec<Observation>;

pub type Trajectory = HashMap<String, Observations>;

#[derive(Error, Debug)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
}

// #[derive(Debug, DeserializeFromStr)]
// struct ObservationWrapper(Trajectory);

// impl FromStr for ObservationWrapper {
//     type Err = ParseObsError;

//     fn from_str(line: &str) -> Result<Self, Self::Err> {
//         if line.len() < 80 {
//             return Err(ParseObsError::TooShortLine);
//         }

//         let observation = Observation {
//             designation: line[0..5].trim().to_string(),
//             packed_designation: line[5..12].trim().to_string(),
//             observation_date: line[15..32].trim().to_string(),
//             ra: line[32..44].trim().to_string(),
//             dec: line[44..56].trim().to_string(),
//             magnitude: line.get(65..70).and_then(|s| s.trim().parse().ok()), 
//             observatory_code: line[77..80].trim().to_string(),
//         };

//         Ok(ObservationWrapper(observation))
//     }
// }

trait TrajectoryExt{
    fn from_ades(ades_file: &Utf8Path) -> Self;
    fn from_80col(colfile: &Utf8Path) -> Self;
}

impl TrajectoryExt for Trajectory {
    fn from_ades(ades_file: &Utf8Path) -> Self {
        todo!()
    }

    fn from_80col(colfile: &Utf8Path) -> Self {
        todo!()
    }
}