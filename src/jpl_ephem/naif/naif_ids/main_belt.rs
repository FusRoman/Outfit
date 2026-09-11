use std::fmt;

use super::ErrorId;

/// NAIF encodes a numbered small body as `2_000_000 + <minor-planet number>`.
const NAIF_ID_OFFSET: i32 = 2_000_000;

/// NAIF small-body identifiers occupy `2_000_000..=2_999_999`.
const NAIF_ID_RANGE_END: i32 = 3_000_000;

/// A numbered main-belt asteroid, identified by its official minor-planet
/// number (e.g. `1` for Ceres, `4` for Vesta).
///
/// This wraps the raw number rather than enumerating every known asteroid: the
/// supplementary ephemeris kernel Outfit loads for this body family
/// (`codes_300ast_20100725.bsp`, see
/// [`crate::propagator::planet_gm::known_main_belt_asteroids`]) covers 300
/// bodies, far more than is practical to name one by one.
///
/// Wrapped in [`crate::jpl_ephem::naif::naif_ids::NaifIds::AST`], this id is
/// resolvable only by the ANISE backend with the supplementary kernel loaded
/// — see that enum's "Ephemeris backend support" section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AsteroidNumber(pub u32);

impl AsteroidNumber {
    /// 1 Ceres.
    pub const CERES: Self = Self(1);
    /// 2 Pallas.
    pub const PALLAS: Self = Self(2);
    /// 4 Vesta.
    pub const VESTA: Self = Self(4);

    /// Build an [`AsteroidNumber`] from a raw NAIF small-body integer code.
    ///
    /// Arguments
    /// -----------------
    /// * `id`: the raw NAIF integer identifier (`2_000_000 + number`).
    ///
    /// Return
    /// ----------
    /// * `Ok(AsteroidNumber)` if `id` falls in the NAIF small-body range
    ///   `2_000_000..=2_999_999`.
    /// * `Err(ErrorId)` otherwise.
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        if (NAIF_ID_OFFSET..NAIF_ID_RANGE_END).contains(&id) {
            Ok(Self((id - NAIF_ID_OFFSET) as u32))
        } else {
            Err(ErrorId::InvalidAsteroidId(id))
        }
    }

    /// Convert back to the raw NAIF integer code (`2_000_000 + number`).
    pub fn to_id(self) -> i32 {
        NAIF_ID_OFFSET + self.0 as i32
    }
}

impl From<AsteroidNumber> for i32 {
    fn from(asteroid: AsteroidNumber) -> Self {
        asteroid.to_id()
    }
}

impl TryFrom<i32> for AsteroidNumber {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        AsteroidNumber::from_id(id)
    }
}

impl fmt::Display for AsteroidNumber {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Asteroid {}", self.0)
    }
}

#[cfg(test)]
mod test_main_belt {
    use super::*;

    #[test]
    fn from_id_round_trips_to_id() {
        for number in [1_u32, 2, 4, 1467] {
            let asteroid = AsteroidNumber(number);
            assert_eq!(AsteroidNumber::from_id(asteroid.to_id()).unwrap(), asteroid);
        }
    }

    #[test]
    fn known_bodies_map_to_expected_naif_ids() {
        assert_eq!(AsteroidNumber::CERES.to_id(), 2_000_001);
        assert_eq!(AsteroidNumber::PALLAS.to_id(), 2_000_002);
        assert_eq!(AsteroidNumber::VESTA.to_id(), 2_000_004);
    }

    #[test]
    fn from_id_rejects_ids_outside_the_small_body_range() {
        assert!(AsteroidNumber::from_id(0).is_err());
        assert!(AsteroidNumber::from_id(399).is_err());
        assert!(AsteroidNumber::from_id(1_999_999).is_err());
        assert!(AsteroidNumber::from_id(3_000_000).is_err());
    }

    #[test]
    fn display_uses_the_generic_label() {
        assert_eq!(AsteroidNumber(1).to_string(), "Asteroid 1");
        assert_eq!(AsteroidNumber(1467).to_string(), "Asteroid 1467");
    }

    #[test]
    fn try_from_matches_from_id() {
        let via_from_id = AsteroidNumber::from_id(2_000_004);
        let via_try_from: Result<AsteroidNumber, ErrorId> = 2_000_004.try_into();
        assert_eq!(via_from_id.unwrap(), via_try_from.unwrap());
    }
}
