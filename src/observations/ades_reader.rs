use std::str::FromStr;

use camino::Utf8Path;
use hifitime::Epoch;
use quick_xml::de::from_str;
use serde::{Deserialize, Deserializer};
use smallvec::SmallVec;

use crate::{
    constants::{ObjectNumber, TrajectorySet},
    outfit::Outfit,
};

use super::observations::Observation;

#[derive(Debug, Deserialize)]
pub struct StructuredAdes {
    #[serde(rename = "obsBlock")]
    pub obs_blocks: Vec<ObsBlock>,
}

#[derive(Debug, Deserialize)]
pub struct FlatAdes {
    #[serde(rename = "optical")]
    pub opticals: Vec<OpticalObs>,
}

#[derive(Debug, Deserialize)]
pub struct ObsBlock {
    #[serde(rename = "obsContext")]
    pub obs_context: ObsContext,

    #[serde(rename = "obsData")]
    pub obs_data: ObsData,
}

#[derive(Debug, Deserialize)]
pub struct ObsContext {
    pub observatory: Observatory,
}

#[derive(Debug, Deserialize)]
pub struct Observatory {
    #[serde(rename = "mpcCode")]
    pub mpc_code: String,
}

#[derive(Debug, Deserialize)]
pub struct ObsData {
    #[serde(rename = "optical")]
    pub opticals: Vec<OpticalObs>,
}

#[derive(Debug, Deserialize)]
pub struct OpticalObs {
    #[serde(rename = "permID")]
    pub perm_id: Option<String>,
    #[serde(rename = "provID")]
    pub prov_id: Option<String>,
    #[serde(rename = "trkSub")]
    pub trk_sub: Option<String>,

    #[serde(rename = "obsTime", deserialize_with = "deserialize_mjd")]
    pub obs_time: f64,

    pub ra: f32,
    pub dec: f32,

    pub stn: String,
}

impl OpticalObs {
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

    fn to_observation(&self, observer_idx: u16) -> Observation {
        Observation::new(observer_idx, self.ra, self.dec, self.obs_time)
    }
}

pub fn deserialize_mjd<'de, D>(deserializer: D) -> Result<f64, D::Error>
where
    D: Deserializer<'de>,
{
    let date_str = String::deserialize(deserializer)?;

    let time = Epoch::from_str(&date_str).map_err(serde::de::Error::custom)?;

    Ok(time.to_mjd_utc_days())
}

fn parse_flat_ades(outfit: &mut Outfit, flat_ades: &FlatAdes, trajs: &mut TrajectorySet) {
    for optical in &flat_ades.opticals {
        let traj_id = optical.get_id();
        let observer = outfit.uint16_from_mpc_code(&optical.stn);
        let observation = optical.to_observation(observer);
        trajs
            .entry(traj_id)
            .or_insert_with(|| SmallVec::with_capacity(10))
            .push(observation);
    }
}

fn parse_structured_ades(
    outfit: &mut Outfit,
    structured_ades: &StructuredAdes,
    trajs: &mut TrajectorySet,
) {
    for obs_block in &structured_ades.obs_blocks {
        let obs_context = &obs_block.obs_context;
        let mpc_code = &obs_context.observatory.mpc_code;
        let observer = outfit.uint16_from_mpc_code(mpc_code);

        for optical in &obs_block.obs_data.opticals {
            let observation = optical.to_observation(observer);
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
pub(crate) fn parse_ades(outfit: &mut Outfit, ades: &Utf8Path, trajs: &mut TrajectorySet) {
    let xml = std::fs::read_to_string(ades)
        .expect(format!("Failed to read ADES file: {}", ades.to_string()).as_str());

        match from_str::<FlatAdes>(&xml) {
            Ok(flat_ades) => {
                parse_flat_ades(outfit, &flat_ades, trajs);
            }
            Err(flat_err) => {
                match from_str::<StructuredAdes>(&xml) {
                    Ok(structured_ades) => {
                        parse_structured_ades(outfit, &structured_ades, trajs);
                    }
                    Err(structured_err) => {
                        panic!(
                            "Failed to parse ADES file:\n- Flat error: {}\n- Structured error: {}",
                            flat_err, structured_err
                        );
                    }
                }
            }
        }
}
