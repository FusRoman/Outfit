use download_jpl_file::{EphemFilePath, EphemFileSource};
use hifitime::Epoch;
use horizon::{horizon_data::HorizonData, horizon_ids::HorizonID};
use naif::{
    naif_data::NaifData,
    naif_ids::{planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds},
};
use nalgebra::Vector3;

pub(super) mod download_jpl_file;
pub mod horizon;
pub mod jpl_request;
pub mod naif;

pub enum JPLEphem {
    HorizonFile(HorizonData),
    NaifFile(NaifData),
}

impl JPLEphem {
    pub fn new(file_source: &EphemFileSource) -> Self {
        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        match file_path {
            EphemFilePath::JPLHorizon(..) => {
                let horizon_data = HorizonData::read_horizon_file(&file_path);
                JPLEphem::HorizonFile(horizon_data)
            }
            EphemFilePath::NAIF(..) => {
                let naif_data = NaifData::read_naif_file(&file_path);
                JPLEphem::NaifFile(naif_data)
            }
        }
    }

    pub fn earth_position_ephemeris(&self, ephem_time: &Epoch) -> Vector3<f64> {
        match self {
            JPLEphem::HorizonFile(horizon_data) => {
                horizon_data
                    .ephemeris(
                        HorizonID::Earth,
                        HorizonID::Sun,
                        ephem_time.to_jde_et_days(),
                        false,
                        false,
                    )
                    .to_au()
                    .position
            }
            JPLEphem::NaifFile(naif_data) => {
                naif_data
                    .ephemeris(
                        NaifIds::PB(PlanetaryBary::EarthMoon),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        ephem_time.to_et_seconds(),
                    )
                    .to_au()
                    .position
            }
        }
    }
}
