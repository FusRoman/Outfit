use std::str::FromStr;

use camino::Utf8Path;
use hifitime::Epoch;
use nom::error;
use quick_xml::de::from_str;
use serde::{Deserialize, Deserializer};
use smallvec::SmallVec;

use crate::{
    constants::{ArcSec, ObjectNumber, TrajectorySet},
    outfit::Outfit,
};

use super::observations::Observation;

#[derive(Debug, Deserialize)]
struct StructuredAdes {
    #[serde(rename = "obsBlock")]
    obs_blocks: Vec<ObsBlock>,
}

#[derive(Debug, Deserialize)]
struct FlatAdes {
    #[serde(rename = "optical")]
    opticals: Vec<OpticalObs>,
}

#[derive(Debug, Deserialize)]
struct ObsBlock {
    #[serde(rename = "obsContext")]
    obs_context: ObsContext,

    #[serde(rename = "obsData")]
    obs_data: ObsData,
}

#[derive(Debug, Deserialize)]
struct ObsContext {
    observatory: Observatory,
}

#[derive(Debug, Deserialize)]
struct Observatory {
    #[serde(rename = "mpcCode")]
    mpc_code: String,
}

#[derive(Debug, Deserialize)]
struct ObsData {
    #[serde(rename = "optical")]
    opticals: Vec<OpticalObs>,
}

#[derive(Debug, Deserialize)]
struct OpticalObs {
    #[serde(rename = "permID")]
    perm_id: Option<String>,
    #[serde(rename = "provID")]
    prov_id: Option<String>,
    #[serde(rename = "trkSub")]
    trk_sub: Option<String>,

    #[serde(rename = "obsTime", deserialize_with = "deserialize_mjd")]
    obs_time: f64,

    ra: f64,
    dec: f64,

    #[serde(rename = "precRA")]
    prec_ra: Option<f64>,
    #[serde(rename = "precDec")]
    prec_dec: Option<f64>,

    #[serde(rename = "rmsRA")]
    rms_ra: Option<f64>,
    #[serde(rename = "rmsDec")]
    rms_dec: Option<f64>,

    stn: String,
}

impl OpticalObs {
    /// Returns the trajectory ID for the optical observation.
    /// It first checks for a `perm_id`, then a `prov_id`, and finally falls back to `trk_sub`.
    /// If none of these are available, it panics with an error message.
    /// The ID is parsed as a `u32` if possible, otherwise it is returned as a string.
    ///
    /// Return
    /// ------
    /// * An `ObjectNumber` representing the trajectory ID.
    /// * If the ID is a valid `u32`, it is returned as `ObjectNumber::Int(id)`.
    /// * If the ID is not a valid `u32`, it is returned as `ObjectNumber::String(id)`.
    /// * If no ID is found, it panics with an error message.
    fn get_id(&self) -> ObjectNumber {
        let id = self
            .perm_id
            .clone()
            .or_else(|| self.prov_id.clone())
            .unwrap_or_else(|| self.trk_sub.clone().expect("No ID found"));
        if let Ok(id) = id.parse::<u32>() {
            ObjectNumber::Int(id)
        } else {
            ObjectNumber::String(id)
        }
    }

    fn to_observation(
        &self,
        observer_idx: u16,
        error_ra: Option<ArcSec>,
        error_dec: Option<ArcSec>,
    ) -> Observation {
        let error_ra = self.rms_ra.unwrap_or_else(|| {
            self.prec_ra
                .unwrap_or_else(|| error_ra.expect("No error for RA when parsing ADES file"))
        });

        let error_dec = self.rms_dec.unwrap_or_else(|| {
            self.prec_dec
                .unwrap_or_else(|| error_dec.expect("No error for Dec when parsing ADES file"))
        });
        Observation::new(
            observer_idx,
            self.ra,
            error_ra,
            self.dec,
            error_dec,
            self.obs_time,
        )
    }
}

/// Deserialize a date string in the format "YYYY-MM-DDTHH:MM:SS" into a floating-point number
/// representing the Modified Julian Date (MJD).
/// The date string is expected to be in UTC.
///
/// Arguments
/// ---------
/// * `deserializer`: The deserializer to use for the date string.
///
/// Return
/// ------
/// * A `Result` containing the parsed MJD as a `f64` or an error if the parsing fails.
fn deserialize_mjd<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let date_str = String::deserialize(deserializer)?;

    let time = Epoch::from_str(&date_str).map_err(serde::de::Error::custom)?;

    Ok(time.to_mjd_utc_days())
}

/// Parses a `FlatAdes` file and populates the given `Outfit` and `TrajectorySet`.
/// It iterates through the optical observations, extracting the observer's MPC code and
/// creating an `Observation` for each optical observation.
/// The observations are then added to the `TrajectorySet` using the trajectory ID.
/// If a new observatory is found, it is added to the `Outfit` observatory set.
///
/// Arguments
/// ---------
/// * `outfit`: A mutable reference to the `Outfit` instance.
/// * `flat_ades`: A reference to the `FlatAdes` instance.
/// * `trajs`: A mutable reference to the `TrajectorySet` instance.
fn parse_flat_ades(
    outfit: &mut Outfit,
    flat_ades: &FlatAdes,
    trajs: &mut TrajectorySet,
    error_ra: Option<ArcSec>,
    error_dec: Option<ArcSec>,
) {
    for optical in &flat_ades.opticals {
        let traj_id = optical.get_id();
        let observer = outfit.uint16_from_mpc_code(&optical.stn);
        let observation = optical.to_observation(observer, error_ra, error_dec);
        trajs
            .entry(traj_id)
            .or_insert_with(|| SmallVec::with_capacity(10))
            .push(observation);
    }
}

/// Parses a `StructuredAdes` file and populates the given `Outfit` and `TrajectorySet`.
/// It iterates through the observation blocks, extracting the observer's MPC code and
/// creating an `Observation` for each optical observation.
/// The observations are then added to the `TrajectorySet` using the trajectory ID.
/// If a new observatory is found, it is added to the `Outfit` observatory set.
///
/// Arguments
/// ---------
/// * `outfit`: A mutable reference to the `Outfit` instance.
/// * `structured_ades`: A reference to the `StructuredAdes` instance.
/// * `trajs`: A mutable reference to the `TrajectorySet` instance.
fn parse_structured_ades(
    outfit: &mut Outfit,
    structured_ades: &StructuredAdes,
    trajs: &mut TrajectorySet,
    error_ra: Option<ArcSec>,
    error_dec: Option<ArcSec>,
) {
    for obs_block in &structured_ades.obs_blocks {
        let obs_context = &obs_block.obs_context;
        let mpc_code = &obs_context.observatory.mpc_code;
        let observer = outfit.uint16_from_mpc_code(mpc_code);

        for optical in &obs_block.obs_data.opticals {
            let observation = optical.to_observation(observer, error_ra, error_dec);
            let traj_id = optical.get_id();
            trajs
                .entry(traj_id)
                .or_insert_with(|| SmallVec::with_capacity(10))
                .push(observation);
        }
    }
}

/// Parses an ADES file and populates the given `Outfit` and `TrajectorySet`.
/// It first attempts to parse the file as a `FlatAdes`, and if that fails, it tries to parse it as a `StructuredAdes`.
/// If both parsing attempts fail, it panics with an error message.
/// If new observatory are found, they are added to the `Outfit` observatory set.
///
/// Arguments
/// ---------
/// * `outfit`: A mutable reference to the `Outfit` instance.
/// * `ades`: A reference to the ADES file path.
/// * `trajs`: A mutable reference to the `TrajectorySet` instance.
pub(crate) fn parse_ades(
    outfit: &mut Outfit,
    ades: &Utf8Path,
    trajs: &mut TrajectorySet,
    error_ra: Option<ArcSec>,
    error_dec: Option<ArcSec>,
) {
    let xml = std::fs::read_to_string(ades)
        .expect(format!("Failed to read ADES file: {}", ades.to_string()).as_str());

    match from_str::<FlatAdes>(&xml) {
        Ok(flat_ades) => {
            parse_flat_ades(outfit, &flat_ades, trajs, error_ra, error_dec);
        }
        Err(flat_err) => match from_str::<StructuredAdes>(&xml) {
            Ok(structured_ades) => {
                parse_structured_ades(outfit, &structured_ades, trajs, error_ra, error_dec);
            }
            Err(structured_err) => {
                panic!(
                    "Failed to parse ADES file:\n- Flat error: {}\n- Structured error: {}",
                    flat_err, structured_err
                );
            }
        },
    }
}
