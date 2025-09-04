use std::fmt;

use super::ErrorId;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlanetaryBary {
    Mercury = 1,
    Venus = 2,
    EarthMoon = 3,
    Mars = 4,
    Jupiter = 5,
    Saturn = 6,
    Uranus = 7,
    Neptune = 8,
    Pluto = 9,
}

impl PlanetaryBary {
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        match id {
            1 => Ok(PlanetaryBary::Mercury),
            2 => Ok(PlanetaryBary::Venus),
            3 => Ok(PlanetaryBary::EarthMoon),
            4 => Ok(PlanetaryBary::Mars),
            5 => Ok(PlanetaryBary::Jupiter),
            6 => Ok(PlanetaryBary::Saturn),
            7 => Ok(PlanetaryBary::Uranus),
            8 => Ok(PlanetaryBary::Neptune),
            9 => Ok(PlanetaryBary::Pluto),
            _ => Err(ErrorId::InvalidPlanetBaryId(id)),
        }
    }

    pub fn to_id(&self) -> i32 {
        match self {
            PlanetaryBary::Mercury => 1,
            PlanetaryBary::Venus => 2,
            PlanetaryBary::EarthMoon => 3,
            PlanetaryBary::Mars => 4,
            PlanetaryBary::Jupiter => 5,
            PlanetaryBary::Saturn => 6,
            PlanetaryBary::Uranus => 7,
            PlanetaryBary::Neptune => 8,
            PlanetaryBary::Pluto => 9,
        }
    }
}

impl From<PlanetaryBary> for i32 {
    fn from(planet: PlanetaryBary) -> Self {
        planet.to_id()
    }
}

impl TryFrom<i32> for PlanetaryBary {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        PlanetaryBary::from_id(id)
    }
}

impl fmt::Display for PlanetaryBary {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            PlanetaryBary::Mercury => "Mercury",
            PlanetaryBary::Venus => "Venus",
            PlanetaryBary::EarthMoon => "Earth-Moon Barycenter",
            PlanetaryBary::Mars => "Mars",
            PlanetaryBary::Jupiter => "Jupiter",
            PlanetaryBary::Saturn => "Saturn",
            PlanetaryBary::Uranus => "Uranus",
            PlanetaryBary::Neptune => "Neptune",
            PlanetaryBary::Pluto => "Pluto",
        };
        write!(f, "{s}")
    }
}

#[cfg(test)]
mod test_planet_bary {
    use super::*;

    #[test]
    fn test_planet_bary_from_id() {
        assert_eq!(PlanetaryBary::from_id(1).unwrap(), PlanetaryBary::Mercury);
        assert_eq!(PlanetaryBary::from_id(2).unwrap(), PlanetaryBary::Venus);
        assert_eq!(PlanetaryBary::from_id(3).unwrap(), PlanetaryBary::EarthMoon);
        assert_eq!(PlanetaryBary::from_id(4).unwrap(), PlanetaryBary::Mars);
        assert_eq!(PlanetaryBary::from_id(5).unwrap(), PlanetaryBary::Jupiter);
        assert_eq!(PlanetaryBary::from_id(6).unwrap(), PlanetaryBary::Saturn);
        assert_eq!(PlanetaryBary::from_id(7).unwrap(), PlanetaryBary::Uranus);
        assert_eq!(PlanetaryBary::from_id(8).unwrap(), PlanetaryBary::Neptune);
        assert_eq!(PlanetaryBary::from_id(9).unwrap(), PlanetaryBary::Pluto);
        assert!(PlanetaryBary::from_id(100).is_err());
    }

    #[test]
    fn test_planet_bary_to_id() {
        assert_eq!(PlanetaryBary::Mercury.to_id(), 1);
        assert_eq!(PlanetaryBary::Venus.to_id(), 2);
        assert_eq!(PlanetaryBary::EarthMoon.to_id(), 3);
        assert_eq!(PlanetaryBary::Mars.to_id(), 4);
        assert_eq!(PlanetaryBary::Jupiter.to_id(), 5);
        assert_eq!(PlanetaryBary::Saturn.to_id(), 6);
        assert_eq!(PlanetaryBary::Uranus.to_id(), 7);
        assert_eq!(PlanetaryBary::Neptune.to_id(), 8);
        assert_eq!(PlanetaryBary::Pluto.to_id(), 9);
    }

    #[test]
    fn test_from() {
        assert_eq!(i32::from(PlanetaryBary::Mercury), 1);
        assert_eq!(i32::from(PlanetaryBary::Venus), 2);
        assert_eq!(i32::from(PlanetaryBary::EarthMoon), 3);
        assert_eq!(i32::from(PlanetaryBary::Mars), 4);
        assert_eq!(i32::from(PlanetaryBary::Jupiter), 5);
        assert_eq!(i32::from(PlanetaryBary::Saturn), 6);
        assert_eq!(i32::from(PlanetaryBary::Uranus), 7);
        assert_eq!(i32::from(PlanetaryBary::Neptune), 8);
        assert_eq!(i32::from(PlanetaryBary::Pluto), 9);
    }

    #[test]
    fn test_into() {
        let mercury: i32 = PlanetaryBary::Mercury.into();
        assert_eq!(mercury, 1);
        let venus: i32 = PlanetaryBary::Venus.into();
        assert_eq!(venus, 2);
        let earth_moon: i32 = PlanetaryBary::EarthMoon.into();
        assert_eq!(earth_moon, 3);
        let mars: i32 = PlanetaryBary::Mars.into();
        assert_eq!(mars, 4);
        let jupiter: i32 = PlanetaryBary::Jupiter.into();
        assert_eq!(jupiter, 5);
        let saturn: i32 = PlanetaryBary::Saturn.into();
        assert_eq!(saturn, 6);
        let uranus: i32 = PlanetaryBary::Uranus.into();
        assert_eq!(uranus, 7);
        let neptune: i32 = PlanetaryBary::Neptune.into();
        assert_eq!(neptune, 8);
        let pluto: i32 = PlanetaryBary::Pluto.into();
        assert_eq!(pluto, 9);
    }

    #[test]
    fn test_try_into() {
        let mercury: Result<PlanetaryBary, ErrorId> = 1.try_into();
        assert_eq!(mercury.unwrap(), PlanetaryBary::Mercury);
        let venus: Result<PlanetaryBary, ErrorId> = 2.try_into();
        assert_eq!(venus.unwrap(), PlanetaryBary::Venus);
        let earth_moon: Result<PlanetaryBary, ErrorId> = 3.try_into();
        assert_eq!(earth_moon.unwrap(), PlanetaryBary::EarthMoon);
        let mars: Result<PlanetaryBary, ErrorId> = 4.try_into();
        assert_eq!(mars.unwrap(), PlanetaryBary::Mars);
        let jupiter: Result<PlanetaryBary, ErrorId> = 5.try_into();
        assert_eq!(jupiter.unwrap(), PlanetaryBary::Jupiter);
        let saturn: Result<PlanetaryBary, ErrorId> = 6.try_into();
        assert_eq!(saturn.unwrap(), PlanetaryBary::Saturn);
        let uranus: Result<PlanetaryBary, ErrorId> = 7.try_into();
        assert_eq!(uranus.unwrap(), PlanetaryBary::Uranus);
        let neptune: Result<PlanetaryBary, ErrorId> = 8.try_into();
        assert_eq!(neptune.unwrap(), PlanetaryBary::Neptune);
        let pluto: Result<PlanetaryBary, ErrorId> = 9.try_into();
        assert_eq!(pluto.unwrap(), PlanetaryBary::Pluto);

        let invalid: Result<PlanetaryBary, ErrorId> = 100.try_into();
        assert!(invalid.is_err());
    }
}
