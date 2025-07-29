pub mod planet_bary;
pub mod planet_mass;
pub mod satellite_mass;
pub mod solar_system_bary;

pub mod naif_type;

use std::fmt;

use planet_bary::PlanetaryBary;
use planet_mass::PlanetMassCenter;
use satellite_mass::SatelliteMassCenter;
use solar_system_bary::SolarSystemBary;
use thiserror::Error;

#[derive(Debug, Clone, Copy, Error)]
pub enum ErrorId {
    #[error("Invalid Planetary Barycenter ID: {0}")]
    InvalidPlanetBaryId(i32),

    #[error("Invalid Planet Mass Center ID: {0}")]
    InvalidPlanetMassCenterId(i32),
    #[error("Invalid Satellite Mass Center ID: {0}")]
    InvalidSatelliteMassCenterId(i32),

    #[error("Invalid NAIF ID: {0}")]
    InvalidNaifId(i32),
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NaifIds {
    SSB(SolarSystemBary),
    PB(PlanetaryBary),
    PMC(PlanetMassCenter),
    SMC(SatelliteMassCenter),
}

impl NaifIds {
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        match id {
            0 => Ok(NaifIds::SSB(SolarSystemBary::SSB)),
            10 => Ok(NaifIds::SSB(SolarSystemBary::Sun)),
            1..=999 => {
                if let Ok(planet) = PlanetaryBary::from_id(id) {
                    Ok(NaifIds::PB(planet))
                } else if let Ok(planet_mass_center) = PlanetMassCenter::from_id(id) {
                    return Ok(NaifIds::PMC(planet_mass_center));
                } else if let Ok(satellite_mass_center) = SatelliteMassCenter::from_id(id) {
                    return Ok(NaifIds::SMC(satellite_mass_center));
                } else {
                    return Err(ErrorId::InvalidNaifId(id));
                }
            }
            _ => Err(ErrorId::InvalidNaifId(id)),
        }
    }

    pub fn to_id(&self) -> i32 {
        match self {
            NaifIds::SSB(solar_system_bary) => match solar_system_bary {
                SolarSystemBary::SSB => 0,
                SolarSystemBary::Sun => 10,
            },
            NaifIds::PB(planetary_bary) => planetary_bary.to_id(),
            NaifIds::PMC(planet_mass_center) => planet_mass_center.to_id(),
            NaifIds::SMC(satellite_mass_center) => satellite_mass_center.to_id(),
        }
    }
}

impl From<NaifIds> for i32 {
    fn from(naif_id: NaifIds) -> Self {
        match naif_id {
            NaifIds::SSB(solar_system_bary) => solar_system_bary.to_id(),
            NaifIds::PB(planetary_bary) => planetary_bary.to_id(),
            NaifIds::PMC(planet_mass_center) => planet_mass_center.to_id(),
            NaifIds::SMC(satellite_mass_center) => satellite_mass_center.to_id(),
        }
    }
}

impl TryFrom<i32> for NaifIds {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        NaifIds::from_id(id)
    }
}

impl fmt::Display for NaifIds {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            NaifIds::SSB(solar_system_bary) => match solar_system_bary {
                SolarSystemBary::SSB => write!(f, "Solar System Barycenter"),
                SolarSystemBary::Sun => write!(f, "Sun"),
            },
            NaifIds::PB(planetary_bary) => write!(f, "{planetary_bary}"),
            NaifIds::PMC(planet_mass_center) => write!(f, "{planet_mass_center}"),
            NaifIds::SMC(satellite_mass_center) => write!(f, "{satellite_mass_center}"),
        }
    }
}

#[cfg(test)]
mod test_naif_id {
    use super::*;

    #[test]
    fn test_naif_ids() {
        assert_eq!(
            NaifIds::from_id(0).unwrap(),
            NaifIds::SSB(SolarSystemBary::SSB)
        );
        assert_eq!(
            NaifIds::from_id(10).unwrap(),
            NaifIds::SSB(SolarSystemBary::Sun)
        );
        assert_eq!(
            NaifIds::from_id(1).unwrap(),
            NaifIds::PB(PlanetaryBary::Mercury)
        );
        assert_eq!(
            NaifIds::from_id(199).unwrap(),
            NaifIds::PMC(PlanetMassCenter::Mercury)
        );
        assert_eq!(
            NaifIds::from_id(301).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Moon)
        );
        assert_eq!(
            NaifIds::from_id(401).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Phobos)
        );
        assert_eq!(
            NaifIds::from_id(901).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Charon)
        );
        assert!(NaifIds::from_id(1000).is_err());
        assert!(NaifIds::from_id(11).is_err());
    }

    #[test]
    fn test_naif_ids_to_id() {
        assert_eq!(NaifIds::SSB(SolarSystemBary::SSB).to_id(), 0);
        assert_eq!(NaifIds::SSB(SolarSystemBary::Sun).to_id(), 10);
        assert_eq!(NaifIds::PB(PlanetaryBary::Mercury).to_id(), 1);
        assert_eq!(NaifIds::PMC(PlanetMassCenter::Mercury).to_id(), 199);
        assert_eq!(NaifIds::SMC(SatelliteMassCenter::Moon).to_id(), 301);
        assert_eq!(NaifIds::SMC(SatelliteMassCenter::Phobos).to_id(), 401);
        assert_eq!(NaifIds::SMC(SatelliteMassCenter::Charon).to_id(), 901);
    }

    #[test]
    fn test_naif_ids_to_string() {
        assert_eq!(
            NaifIds::SSB(SolarSystemBary::SSB).to_string(),
            "Solar System Barycenter"
        );
        assert_eq!(NaifIds::SSB(SolarSystemBary::Sun).to_string(), "Sun");
        assert_eq!(NaifIds::PB(PlanetaryBary::Mercury).to_string(), "Mercury");
        assert_eq!(
            NaifIds::PMC(PlanetMassCenter::Mercury).to_string(),
            "Mercury"
        );
        assert_eq!(NaifIds::SMC(SatelliteMassCenter::Moon).to_string(), "Moon");
        assert_eq!(
            NaifIds::SMC(SatelliteMassCenter::Phobos).to_string(),
            "Phobos"
        );
        assert_eq!(
            NaifIds::SMC(SatelliteMassCenter::Charon).to_string(),
            "Charon"
        );
    }

    #[test]
    fn test_naif_ids_try_from() {
        assert_eq!(
            NaifIds::try_from(0).unwrap(),
            NaifIds::SSB(SolarSystemBary::SSB)
        );
        assert_eq!(
            NaifIds::try_from(10).unwrap(),
            NaifIds::SSB(SolarSystemBary::Sun)
        );
        assert_eq!(
            NaifIds::try_from(1).unwrap(),
            NaifIds::PB(PlanetaryBary::Mercury)
        );
        assert_eq!(
            NaifIds::try_from(199).unwrap(),
            NaifIds::PMC(PlanetMassCenter::Mercury)
        );
        assert_eq!(
            NaifIds::try_from(301).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Moon)
        );
        assert_eq!(
            NaifIds::try_from(401).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Phobos)
        );
        assert_eq!(
            NaifIds::try_from(901).unwrap(),
            NaifIds::SMC(SatelliteMassCenter::Charon)
        );
        assert!(NaifIds::try_from(1000).is_err());
    }

    #[test]
    fn test_naif_ids_from() {
        assert_eq!(i32::from(NaifIds::SSB(SolarSystemBary::SSB)), 0);
        assert_eq!(i32::from(NaifIds::SSB(SolarSystemBary::Sun)), 10);
        assert_eq!(i32::from(NaifIds::PB(PlanetaryBary::Mercury)), 1);
        assert_eq!(i32::from(NaifIds::PMC(PlanetMassCenter::Mercury)), 199);
        assert_eq!(i32::from(NaifIds::SMC(SatelliteMassCenter::Moon)), 301);
        assert_eq!(i32::from(NaifIds::SMC(SatelliteMassCenter::Phobos)), 401);
        assert_eq!(i32::from(NaifIds::SMC(SatelliteMassCenter::Charon)), 901);
    }
}
