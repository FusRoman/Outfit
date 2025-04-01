use std::{collections::HashMap, str::FromStr};

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
#[serde(rename_all = "camelCase")]
pub struct Ades {
    #[serde(rename = "obsBlock")]
    pub obs_blocks: Option<Vec<ObsBlock>>,

    #[serde(rename = "optical")]
    pub flat_opticals: Option<Vec<OpticalObs>>,
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

pub(crate) fn parse_ades(outfit: &mut Outfit, ades: &Utf8Path, trajs: &mut TrajectorySet) {
    let xml = std::fs::read_to_string(ades).unwrap();
    let ades: Ades = from_str(&xml).unwrap();

    match ades.obs_blocks {
        Some(obs_blocks) => {
            for obs_block in obs_blocks {
                let obs_context = obs_block.obs_context;
                let mpc_code = obs_context.observatory.mpc_code;
                let observer = outfit.uint16_from_mpc_code(&mpc_code);

                for optical in obs_block.obs_data.opticals {
                    let observation = optical.to_observation(observer);
                    let traj_id = optical.get_id();
                    trajs
                        .entry(traj_id)
                        .or_insert_with(|| SmallVec::with_capacity(10))
                        .push(observation);
                }
            }
        }
        None => {}
    }
}
