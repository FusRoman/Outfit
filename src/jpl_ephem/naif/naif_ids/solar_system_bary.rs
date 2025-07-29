use std::fmt;

use super::ErrorId;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SolarSystemBary {
    SSB = 0,
    Sun = 10,
}

impl SolarSystemBary {
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        match id {
            0 => Ok(SolarSystemBary::SSB),
            10 => Ok(SolarSystemBary::Sun),
            _ => Err(ErrorId::InvalidNaifId(id)),
        }
    }

    pub fn to_id(&self) -> i32 {
        match self {
            SolarSystemBary::SSB => 0,
            SolarSystemBary::Sun => 10,
        }
    }
}

impl From<SolarSystemBary> for i32 {
    fn from(solar_system_bary: SolarSystemBary) -> Self {
        match solar_system_bary {
            SolarSystemBary::SSB => 0,
            SolarSystemBary::Sun => 10,
        }
    }
}

impl TryFrom<i32> for SolarSystemBary {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        match id {
            0 => Ok(SolarSystemBary::SSB),
            10 => Ok(SolarSystemBary::Sun),
            _ => Err(ErrorId::InvalidNaifId(id)),
        }
    }
}

impl fmt::Display for SolarSystemBary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SolarSystemBary::SSB => "Solar System Barycenter",
            SolarSystemBary::Sun => "Sun",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod test_solarsystem_bary {
    use super::*;

    #[test]
    fn test_solar_system_bary() {
        assert_eq!(SolarSystemBary::from_id(0).unwrap(), SolarSystemBary::SSB);
        assert_eq!(SolarSystemBary::from_id(10).unwrap(), SolarSystemBary::Sun);
        assert!(SolarSystemBary::from_id(1).is_err());
    }

    #[test]
    fn test_solar_system_bary_to_id() {
        assert_eq!(SolarSystemBary::SSB.to_id(), 0);
        assert_eq!(SolarSystemBary::Sun.to_id(), 10);
    }

    #[test]
    fn test_solar_system_bary_to_string() {
        assert_eq!(SolarSystemBary::SSB.to_string(), "Solar System Barycenter");
        assert_eq!(SolarSystemBary::Sun.to_string(), "Sun");
    }

    #[test]
    fn test_solar_system_bary_from_id() {
        assert_eq!(SolarSystemBary::from_id(0).unwrap(), SolarSystemBary::SSB);
        assert_eq!(SolarSystemBary::from_id(10).unwrap(), SolarSystemBary::Sun);
        assert!(SolarSystemBary::from_id(1).is_err());
    }

    #[test]
    fn test_solar_system_bary_try_from() {
        assert_eq!(SolarSystemBary::try_from(0).unwrap(), SolarSystemBary::SSB);
        assert_eq!(SolarSystemBary::try_from(10).unwrap(), SolarSystemBary::Sun);
        assert!(SolarSystemBary::try_from(1).is_err());
    }
}
