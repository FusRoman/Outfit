use super::ErrorId;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PlanetMassCenter {
    Mercury = 199,
    Venus = 299,
    Earth = 399,
    Mars = 499,
    Jupiter = 599,
    Saturn = 699,
    Uranus = 799,
    Neptune = 899,
    Pluto = 999,
}

impl PlanetMassCenter {
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        match id {
            199 => Ok(PlanetMassCenter::Mercury),
            299 => Ok(PlanetMassCenter::Venus),
            399 => Ok(PlanetMassCenter::Earth),
            499 => Ok(PlanetMassCenter::Mars),
            599 => Ok(PlanetMassCenter::Jupiter),
            699 => Ok(PlanetMassCenter::Saturn),
            799 => Ok(PlanetMassCenter::Uranus),
            899 => Ok(PlanetMassCenter::Neptune),
            999 => Ok(PlanetMassCenter::Pluto),
            _ => Err(ErrorId::InvalidPlanetMassCenterId(id)),
        }
    }

    pub fn to_id(&self) -> i32 {
        match self {
            PlanetMassCenter::Mercury => 199,
            PlanetMassCenter::Venus => 299,
            PlanetMassCenter::Earth => 399,
            PlanetMassCenter::Mars => 499,
            PlanetMassCenter::Jupiter => 599,
            PlanetMassCenter::Saturn => 699,
            PlanetMassCenter::Uranus => 799,
            PlanetMassCenter::Neptune => 899,
            PlanetMassCenter::Pluto => 999,
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            PlanetMassCenter::Mercury => "Mercury".to_string(),
            PlanetMassCenter::Venus => "Venus".to_string(),
            PlanetMassCenter::Earth => "Earth".to_string(),
            PlanetMassCenter::Mars => "Mars".to_string(),
            PlanetMassCenter::Jupiter => "Jupiter".to_string(),
            PlanetMassCenter::Saturn => "Saturn".to_string(),
            PlanetMassCenter::Uranus => "Uranus".to_string(),
            PlanetMassCenter::Neptune => "Neptune".to_string(),
            PlanetMassCenter::Pluto => "Pluto".to_string(),
        }
    }
}

impl From<PlanetMassCenter> for i32 {
    fn from(planet_mass_center: PlanetMassCenter) -> Self {
        planet_mass_center.to_id()
    }
}

impl TryFrom<i32> for PlanetMassCenter {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        PlanetMassCenter::from_id(id)
    }
}

impl std::fmt::Display for PlanetMassCenter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

#[cfg(test)]
mod planet_mass_test {
    use super::*;

    #[test]
    fn test_planet_mass_center() {
        assert_eq!(
            PlanetMassCenter::from_id(199).unwrap(),
            PlanetMassCenter::Mercury
        );
        assert_eq!(
            PlanetMassCenter::from_id(299).unwrap(),
            PlanetMassCenter::Venus
        );
        assert_eq!(
            PlanetMassCenter::from_id(399).unwrap(),
            PlanetMassCenter::Earth
        );
        assert_eq!(
            PlanetMassCenter::from_id(499).unwrap(),
            PlanetMassCenter::Mars
        );
        assert_eq!(
            PlanetMassCenter::from_id(599).unwrap(),
            PlanetMassCenter::Jupiter
        );
        assert_eq!(
            PlanetMassCenter::from_id(699).unwrap(),
            PlanetMassCenter::Saturn
        );
        assert_eq!(
            PlanetMassCenter::from_id(799).unwrap(),
            PlanetMassCenter::Uranus
        );
        assert_eq!(
            PlanetMassCenter::from_id(899).unwrap(),
            PlanetMassCenter::Neptune
        );
        assert_eq!(
            PlanetMassCenter::from_id(999).unwrap(),
            PlanetMassCenter::Pluto
        );
        assert!(PlanetMassCenter::from_id(1000).is_err());
    }

    #[test]
    fn test_to_id() {
        assert_eq!(PlanetMassCenter::Mercury.to_id(), 199);
        assert_eq!(PlanetMassCenter::Venus.to_id(), 299);
        assert_eq!(PlanetMassCenter::Earth.to_id(), 399);
        assert_eq!(PlanetMassCenter::Mars.to_id(), 499);
        assert_eq!(PlanetMassCenter::Jupiter.to_id(), 599);
        assert_eq!(PlanetMassCenter::Saturn.to_id(), 699);
        assert_eq!(PlanetMassCenter::Uranus.to_id(), 799);
        assert_eq!(PlanetMassCenter::Neptune.to_id(), 899);
        assert_eq!(PlanetMassCenter::Pluto.to_id(), 999);
    }

    #[test]
    fn test_to_string() {
        assert_eq!(PlanetMassCenter::Mercury.to_string(), "Mercury");
        assert_eq!(PlanetMassCenter::Venus.to_string(), "Venus");
        assert_eq!(PlanetMassCenter::Earth.to_string(), "Earth");
        assert_eq!(PlanetMassCenter::Mars.to_string(), "Mars");
        assert_eq!(PlanetMassCenter::Jupiter.to_string(), "Jupiter");
        assert_eq!(PlanetMassCenter::Saturn.to_string(), "Saturn");
        assert_eq!(PlanetMassCenter::Uranus.to_string(), "Uranus");
        assert_eq!(PlanetMassCenter::Neptune.to_string(), "Neptune");
        assert_eq!(PlanetMassCenter::Pluto.to_string(), "Pluto");
    }

    #[test]
    fn test_from() {
        assert_eq!(i32::from(PlanetMassCenter::Mercury), 199);
        assert_eq!(i32::from(PlanetMassCenter::Venus), 299);
        assert_eq!(i32::from(PlanetMassCenter::Earth), 399);
        assert_eq!(i32::from(PlanetMassCenter::Mars), 499);
        assert_eq!(i32::from(PlanetMassCenter::Jupiter), 599);
        assert_eq!(i32::from(PlanetMassCenter::Saturn), 699);
        assert_eq!(i32::from(PlanetMassCenter::Uranus), 799);
        assert_eq!(i32::from(PlanetMassCenter::Neptune), 899);
        assert_eq!(i32::from(PlanetMassCenter::Pluto), 999);
    }

    #[test]
    fn test_try_from() {
        assert_eq!(
            PlanetMassCenter::try_from(199).unwrap(),
            PlanetMassCenter::Mercury
        );
        assert_eq!(
            PlanetMassCenter::try_from(299).unwrap(),
            PlanetMassCenter::Venus
        );
        assert_eq!(
            PlanetMassCenter::try_from(399).unwrap(),
            PlanetMassCenter::Earth
        );
        assert_eq!(
            PlanetMassCenter::try_from(499).unwrap(),
            PlanetMassCenter::Mars
        );
        assert_eq!(
            PlanetMassCenter::try_from(599).unwrap(),
            PlanetMassCenter::Jupiter
        );
        assert_eq!(
            PlanetMassCenter::try_from(699).unwrap(),
            PlanetMassCenter::Saturn
        );
        assert_eq!(
            PlanetMassCenter::try_from(799).unwrap(),
            PlanetMassCenter::Uranus
        );
        assert_eq!(
            PlanetMassCenter::try_from(899).unwrap(),
            PlanetMassCenter::Neptune
        );
        assert_eq!(
            PlanetMassCenter::try_from(999).unwrap(),
            PlanetMassCenter::Pluto
        );
        assert!(PlanetMassCenter::try_from(1000).is_err());
    }
}
