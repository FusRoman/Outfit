use super::ErrorId;

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SatelliteMassCenter {
    Moon = 301,
    Phobos = 401,
    Deimos = 402,
    Io = 501,
    Europa = 502,
    Ganymede = 503,
    Callisto = 504,
    Titan = 601,
    Rhea = 602,
    Iapetus = 603,
    Miranda = 701,
    Ariel = 702,
    Umbriel = 703,
    Titania = 704,
    Oberon = 705,
    Triton = 801,
    Charon = 901,
}

impl SatelliteMassCenter {
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        match id {
            301 => Ok(SatelliteMassCenter::Moon),
            401 => Ok(SatelliteMassCenter::Phobos),
            402 => Ok(SatelliteMassCenter::Deimos),
            501 => Ok(SatelliteMassCenter::Io),
            502 => Ok(SatelliteMassCenter::Europa),
            503 => Ok(SatelliteMassCenter::Ganymede),
            504 => Ok(SatelliteMassCenter::Callisto),
            601 => Ok(SatelliteMassCenter::Titan),
            602 => Ok(SatelliteMassCenter::Rhea),
            603 => Ok(SatelliteMassCenter::Iapetus),
            701 => Ok(SatelliteMassCenter::Miranda),
            702 => Ok(SatelliteMassCenter::Ariel),
            703 => Ok(SatelliteMassCenter::Umbriel),
            704 => Ok(SatelliteMassCenter::Titania),
            705 => Ok(SatelliteMassCenter::Oberon),
            801 => Ok(SatelliteMassCenter::Triton),
            901 => Ok(SatelliteMassCenter::Charon),
            _ => Err(ErrorId::InvalidSatelliteMassCenterId(id)),
        }
    }

    pub fn to_id(&self) -> i32 {
        match self {
            SatelliteMassCenter::Moon => 301,
            SatelliteMassCenter::Phobos => 401,
            SatelliteMassCenter::Deimos => 402,
            SatelliteMassCenter::Io => 501,
            SatelliteMassCenter::Europa => 502,
            SatelliteMassCenter::Ganymede => 503,
            SatelliteMassCenter::Callisto => 504,
            SatelliteMassCenter::Titan => 601,
            SatelliteMassCenter::Rhea => 602,
            SatelliteMassCenter::Iapetus => 603,
            SatelliteMassCenter::Miranda => 701,
            SatelliteMassCenter::Ariel => 702,
            SatelliteMassCenter::Umbriel => 703,
            SatelliteMassCenter::Titania => 704,
            SatelliteMassCenter::Oberon => 705,
            SatelliteMassCenter::Triton => 801,
            SatelliteMassCenter::Charon => 901,
        }
    }

    pub fn to_string(&self) -> String {
        match self {
            SatelliteMassCenter::Moon => "Moon".to_string(),
            SatelliteMassCenter::Phobos => "Phobos".to_string(),
            SatelliteMassCenter::Deimos => "Deimos".to_string(),
            SatelliteMassCenter::Io => "Io".to_string(),
            SatelliteMassCenter::Europa => "Europa".to_string(),
            SatelliteMassCenter::Ganymede => "Ganymede".to_string(),
            SatelliteMassCenter::Callisto => "Callisto".to_string(),
            SatelliteMassCenter::Titan => "Titan".to_string(),
            SatelliteMassCenter::Rhea => "Rhea".to_string(),
            SatelliteMassCenter::Iapetus => "Iapetus".to_string(),
            SatelliteMassCenter::Miranda => "Miranda".to_string(),
            SatelliteMassCenter::Ariel => "Ariel".to_string(),
            SatelliteMassCenter::Umbriel => "Umbriel".to_string(),
            SatelliteMassCenter::Titania => "Titania".to_string(),
            SatelliteMassCenter::Oberon => "Oberon".to_string(),
            SatelliteMassCenter::Triton => "Triton".to_string(),
            SatelliteMassCenter::Charon => "Charon".to_string(),
        }
    }
}

impl From<SatelliteMassCenter> for i32 {
    fn from(satellite_mass_center: SatelliteMassCenter) -> Self {
        satellite_mass_center.to_id()
    }
}

impl TryFrom<i32> for SatelliteMassCenter {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        SatelliteMassCenter::from_id(id)
    }
}

impl std::fmt::Display for SatelliteMassCenter {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

#[cfg(test)]

mod test_sat_mass {
    use super::*;

    #[test]
    fn test_satellite_mass_center() {
        assert_eq!(
            SatelliteMassCenter::from_id(301).unwrap(),
            SatelliteMassCenter::Moon
        );
        assert_eq!(
            SatelliteMassCenter::from_id(401).unwrap(),
            SatelliteMassCenter::Phobos
        );
        assert_eq!(
            SatelliteMassCenter::from_id(402).unwrap(),
            SatelliteMassCenter::Deimos
        );
        assert_eq!(
            SatelliteMassCenter::from_id(501).unwrap(),
            SatelliteMassCenter::Io
        );
        assert_eq!(
            SatelliteMassCenter::from_id(502).unwrap(),
            SatelliteMassCenter::Europa
        );
        assert_eq!(
            SatelliteMassCenter::from_id(503).unwrap(),
            SatelliteMassCenter::Ganymede
        );
        assert_eq!(
            SatelliteMassCenter::from_id(504).unwrap(),
            SatelliteMassCenter::Callisto
        );
        assert_eq!(
            SatelliteMassCenter::from_id(601).unwrap(),
            SatelliteMassCenter::Titan
        );
        assert_eq!(
            SatelliteMassCenter::from_id(602).unwrap(),
            SatelliteMassCenter::Rhea
        );
        assert_eq!(
            SatelliteMassCenter::from_id(603).unwrap(),
            SatelliteMassCenter::Iapetus
        );
        assert_eq!(
            SatelliteMassCenter::from_id(701).unwrap(),
            SatelliteMassCenter::Miranda
        );
        assert_eq!(
            SatelliteMassCenter::from_id(702).unwrap(),
            SatelliteMassCenter::Ariel
        );
        assert_eq!(
            SatelliteMassCenter::from_id(703).unwrap(),
            SatelliteMassCenter::Umbriel
        );
        assert_eq!(
            SatelliteMassCenter::from_id(704).unwrap(),
            SatelliteMassCenter::Titania
        );
        assert_eq!(
            SatelliteMassCenter::from_id(705).unwrap(),
            SatelliteMassCenter::Oberon
        );
        assert_eq!(
            SatelliteMassCenter::from_id(801).unwrap(),
            SatelliteMassCenter::Triton
        );
        assert_eq!(
            SatelliteMassCenter::from_id(901).unwrap(),
            SatelliteMassCenter::Charon
        );
        assert!(SatelliteMassCenter::from_id(1000).is_err());
    }

    #[test]
    fn test_to_id() {
        assert_eq!(SatelliteMassCenter::Moon.to_id(), 301);
        assert_eq!(SatelliteMassCenter::Phobos.to_id(), 401);
        assert_eq!(SatelliteMassCenter::Deimos.to_id(), 402);
        assert_eq!(SatelliteMassCenter::Io.to_id(), 501);
        assert_eq!(SatelliteMassCenter::Europa.to_id(), 502);
        assert_eq!(SatelliteMassCenter::Ganymede.to_id(), 503);
        assert_eq!(SatelliteMassCenter::Callisto.to_id(), 504);
        assert_eq!(SatelliteMassCenter::Titan.to_id(), 601);
        assert_eq!(SatelliteMassCenter::Rhea.to_id(), 602);
        assert_eq!(SatelliteMassCenter::Iapetus.to_id(), 603);
        assert_eq!(SatelliteMassCenter::Miranda.to_id(), 701);
        assert_eq!(SatelliteMassCenter::Ariel.to_id(), 702);
        assert_eq!(SatelliteMassCenter::Umbriel.to_id(), 703);
        assert_eq!(SatelliteMassCenter::Titania.to_id(), 704);
        assert_eq!(SatelliteMassCenter::Oberon.to_id(), 705);
        assert_eq!(SatelliteMassCenter::Triton.to_id(), 801);
        assert_eq!(SatelliteMassCenter::Charon.to_id(), 901);
    }

    #[test]
    fn test_to_string() {
        assert_eq!(SatelliteMassCenter::Moon.to_string(), "Moon");
        assert_eq!(SatelliteMassCenter::Phobos.to_string(), "Phobos");
        assert_eq!(SatelliteMassCenter::Deimos.to_string(), "Deimos");
        assert_eq!(SatelliteMassCenter::Io.to_string(), "Io");
        assert_eq!(SatelliteMassCenter::Europa.to_string(), "Europa");
        assert_eq!(SatelliteMassCenter::Ganymede.to_string(), "Ganymede");
        assert_eq!(SatelliteMassCenter::Callisto.to_string(), "Callisto");
        assert_eq!(SatelliteMassCenter::Titan.to_string(), "Titan");
        assert_eq!(SatelliteMassCenter::Rhea.to_string(), "Rhea");
        assert_eq!(SatelliteMassCenter::Iapetus.to_string(), "Iapetus");
        assert_eq!(SatelliteMassCenter::Miranda.to_string(), "Miranda");
        assert_eq!(SatelliteMassCenter::Ariel.to_string(), "Ariel");
        assert_eq!(SatelliteMassCenter::Umbriel.to_string(), "Umbriel");
        assert_eq!(SatelliteMassCenter::Titania.to_string(), "Titania");
        assert_eq!(SatelliteMassCenter::Oberon.to_string(), "Oberon");
        assert_eq!(SatelliteMassCenter::Triton.to_string(), "Triton");
        assert_eq!(SatelliteMassCenter::Charon.to_string(), "Charon");
    }

    #[test]
    fn test_from() {
        assert_eq!(i32::from(SatelliteMassCenter::Moon), 301);
        assert_eq!(i32::from(SatelliteMassCenter::Phobos), 401);
        assert_eq!(i32::from(SatelliteMassCenter::Deimos), 402);
        assert_eq!(i32::from(SatelliteMassCenter::Io), 501);
        assert_eq!(i32::from(SatelliteMassCenter::Europa), 502);
        assert_eq!(i32::from(SatelliteMassCenter::Ganymede), 503);
        assert_eq!(i32::from(SatelliteMassCenter::Callisto), 504);
        assert_eq!(i32::from(SatelliteMassCenter::Titan), 601);
        assert_eq!(i32::from(SatelliteMassCenter::Rhea), 602);
        assert_eq!(i32::from(SatelliteMassCenter::Iapetus), 603);
        assert_eq!(i32::from(SatelliteMassCenter::Miranda), 701);
        assert_eq!(i32::from(SatelliteMassCenter::Ariel), 702);
        assert_eq!(i32::from(SatelliteMassCenter::Umbriel), 703);
        assert_eq!(i32::from(SatelliteMassCenter::Titania), 704);
        assert_eq!(i32::from(SatelliteMassCenter::Oberon), 705);
        assert_eq!(i32::from(SatelliteMassCenter::Triton), 801);
        assert_eq!(i32::from(SatelliteMassCenter::Charon), 901);
    }

    #[test]
    fn test_try_from() {
        assert_eq!(
            SatelliteMassCenter::try_from(301).unwrap(),
            SatelliteMassCenter::Moon
        );
        assert_eq!(
            SatelliteMassCenter::try_from(401).unwrap(),
            SatelliteMassCenter::Phobos
        );
        assert_eq!(
            SatelliteMassCenter::try_from(402).unwrap(),
            SatelliteMassCenter::Deimos
        );
        assert_eq!(
            SatelliteMassCenter::try_from(501).unwrap(),
            SatelliteMassCenter::Io
        );
        assert_eq!(
            SatelliteMassCenter::try_from(502).unwrap(),
            SatelliteMassCenter::Europa
        );
        assert_eq!(
            SatelliteMassCenter::try_from(503).unwrap(),
            SatelliteMassCenter::Ganymede
        );
        assert_eq!(
            SatelliteMassCenter::try_from(504).unwrap(),
            SatelliteMassCenter::Callisto
        );
        assert_eq!(
            SatelliteMassCenter::try_from(601).unwrap(),
            SatelliteMassCenter::Titan
        );
        assert_eq!(
            SatelliteMassCenter::try_from(602).unwrap(),
            SatelliteMassCenter::Rhea
        );
        assert_eq!(
            SatelliteMassCenter::try_from(603).unwrap(),
            SatelliteMassCenter::Iapetus
        );
        assert_eq!(
            SatelliteMassCenter::try_from(701).unwrap(),
            SatelliteMassCenter::Miranda
        );
        assert_eq!(
            SatelliteMassCenter::try_from(702).unwrap(),
            SatelliteMassCenter::Ariel
        );
        assert_eq!(
            SatelliteMassCenter::try_from(703).unwrap(),
            SatelliteMassCenter::Umbriel
        );
        assert_eq!(
            SatelliteMassCenter::try_from(704).unwrap(),
            SatelliteMassCenter::Titania
        );
        assert_eq!(
            SatelliteMassCenter::try_from(705).unwrap(),
            SatelliteMassCenter::Oberon
        );
        assert_eq!(
            SatelliteMassCenter::try_from(801).unwrap(),
            SatelliteMassCenter::Triton
        );
        assert_eq!(
            SatelliteMassCenter::try_from(901).unwrap(),
            SatelliteMassCenter::Charon
        );

        assert!(SatelliteMassCenter::try_from(1000).is_err());
    }
}
