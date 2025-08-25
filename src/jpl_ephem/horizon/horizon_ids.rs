/// Identifier for bodies and reference points in JPL Horizons ephemeris files.
///
/// Horizons binary ephemerides index data blocks by an integer ID (`0..=14`).
/// This enum provides a typed mapping between those raw integers and their
/// corresponding physical meaning (planet, barycenter, etc.).
///
/// Conversions
/// -----------
/// * Use [`TryFrom<u8>`] to convert from a raw Horizons integer ID to a
///   `HorizonID` enum value. Invalid values return an error.
/// * Use [`From<HorizonID>`] to recover the integer index (`u8`) for storage
///   or file access.
///
/// Examples
/// --------
/// ```rust, no_run
/// use crate::jpl_ephem::horizon::horizon_ids::HorizonID;
/// use std::convert::TryFrom;
///
/// let id = HorizonID::try_from(4).unwrap();
/// assert_eq!(id, HorizonID::Jupiter);
///
/// let raw: u8 = HorizonID::Mars.into();
/// assert_eq!(raw, 3);
/// ```
///
/// See also
/// --------
/// * [`HorizonRecord`](crate::jpl_ephem::horizon::horizon_records::HorizonRecord) – per-body ephemeris segments keyed by `HorizonID`.
/// * [`HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData) – container holding header and all per-body records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HorizonID {
    Mercury = 0,
    Venus = 1,
    Earth = 2,
    Mars = 3,
    Jupiter = 4,
    Saturn = 5,
    Uranus = 6,
    Neptune = 7,
    Pluto = 8,
    Moon = 9,
    Sun = 10,
    SolarSystemBarycenter = 11,
    EarthMoonBarycenter = 12,
    Nutation = 13,
    Libration = 14,
}

impl TryFrom<u8> for HorizonID {
    type Error = String;

    fn try_from(value: u8) -> Result<Self, Self::Error> {
        match value {
            0 => Ok(HorizonID::Mercury),
            1 => Ok(HorizonID::Venus),
            2 => Ok(HorizonID::Earth),
            3 => Ok(HorizonID::Mars),
            4 => Ok(HorizonID::Jupiter),
            5 => Ok(HorizonID::Saturn),
            6 => Ok(HorizonID::Uranus),
            7 => Ok(HorizonID::Neptune),
            8 => Ok(HorizonID::Pluto),
            9 => Ok(HorizonID::Moon),
            10 => Ok(HorizonID::Sun),
            11 => Ok(HorizonID::SolarSystemBarycenter),
            12 => Ok(HorizonID::EarthMoonBarycenter),
            13 => Ok(HorizonID::Nutation),
            14 => Ok(HorizonID::Libration),
            _ => Err(format!("Invalid Horizon ID: {value}")),
        }
    }
}

impl From<HorizonID> for u8 {
    fn from(id: HorizonID) -> Self {
        id as u8
    }
}
