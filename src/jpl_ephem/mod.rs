use download_jpl_file::{EphemFilePath, EphemFileSource};
use hifitime::Epoch;
use horizon::{horizon_data::HorizonData, horizon_ids::HorizonID};
use naif::{
    naif_data::NaifData,
    naif_ids::{planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds},
};
use nalgebra::Vector3;

use crate::outfit_errors::OutfitError;

pub mod download_jpl_file;
pub mod horizon;
pub mod naif;

#[derive(Debug, Clone)]
pub enum JPLEphem {
    HorizonFile(HorizonData),
    NaifFile(NaifData),
}

impl JPLEphem {
    pub fn new(file_source: &EphemFileSource) -> Result<Self, OutfitError> {
        let file_path = EphemFilePath::get_ephemeris_file(file_source)?;
        match file_path {
            EphemFilePath::JPLHorizon(..) => {
                let horizon_data = HorizonData::read_horizon_file(&file_path);
                Ok(JPLEphem::HorizonFile(horizon_data))
            }
            EphemFilePath::Naif(..) => {
                let naif_data = NaifData::read_naif_file(&file_path);
                Ok(JPLEphem::NaifFile(naif_data))
            }
        }
    }

    pub fn earth_ephemeris(
        &self,
        ephem_time: &Epoch,
        compute_velocity: bool,
    ) -> (Vector3<f64>, Option<Vector3<f64>>) {
        match self {
            JPLEphem::HorizonFile(horizon_data) => {
                let ephem_res = horizon_data
                    .ephemeris(
                        HorizonID::Earth,
                        HorizonID::Sun,
                        ephem_time.to_mjd_tt_days(),
                        compute_velocity,
                        false,
                    )
                    .to_au();
                (ephem_res.position, ephem_res.velocity)
            }
            JPLEphem::NaifFile(naif_data) => {
                let ephem_res = naif_data
                    .ephemeris(
                        NaifIds::PB(PlanetaryBary::EarthMoon),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        ephem_time.to_et_seconds(),
                    )
                    .to_au();
                (ephem_res.position, ephem_res.velocity)
            }
        }
    }
}
