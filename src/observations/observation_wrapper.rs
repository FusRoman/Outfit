use std::str::FromStr;

use serde_with::DeserializeFromStr;
use thiserror::Error;

use crate::{conversion::{parse_dec_to_deg, parse_ra_to_deg}, time::frac_date_to_mjd};

use super::observations::Observation;

#[derive(Error, Debug)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
    #[error("The line is not a CCD observation")]
    NotCCDObs,
}

#[derive(Debug, DeserializeFromStr)]
pub (crate) struct ObservationWrapper(pub(crate) Observation);

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
            observer: line[77..80].trim().into(),
        };

        Ok(ObservationWrapper(observation))
    }
}