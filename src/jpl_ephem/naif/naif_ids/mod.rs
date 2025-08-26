//! NAIF identifier model for solar‑system bodies.
//!
//! This module provides a strongly‑typed interface to the NAIF/SPICE
//! integer identifier space used to address solar‑system objects such as the
//! Solar System Barycenter (SSB), planetary barycenters, planetary mass centers,
//! and natural satellites.
//!
//! # Overview
//! NAIF IDs are compact `i32` codes. This module groups them into safe enums
//! with conversions in both directions (`try_from<i32>` ⇄ `NaifIds`, and
//! `From<NaifIds> for i32`). It exposes four categories:
//!
//! * `SolarSystemBary` — e.g. `0` for SSB, `10` for the Sun.
//! * `PlanetaryBary` — planetary barycenters (e.g. `1` Mercury barycenter, …).
//! * `PlanetMassCenter` — mass centers (e.g. `199` Mercury, `399` Earth).
//! * `SatelliteMassCenter` — natural satellites (e.g. `301` Moon, `401` Phobos).
//!
//! The exact mapping tables live in the submodules of this namespace
//! (`planet_bary`, `planet_mass`, `satellite_mass`, …). This `mod.rs`
//! concentrates the *dispatching* logic and the ergonomic `NaifIds` enum.
//!
//! # Typical usage
//! ```rust, no_run
//! use outfit::jpl_ephem::naif::naif_ids::NaifIds;
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Parse an integer NAIF code coming from a file or API:
//! let code = 399_i32; // Earth mass center
//! let id = NaifIds::try_from(code)?;
//! assert!(matches!(id, NaifIds::PMC(_)));
//!
//! // Convert back to the raw i32 code (lossless):
//! let roundtrip: i32 = id.into();
//! assert_eq!(roundtrip, 399);
//!
//! // Human‑readable label (Display impl):
//! assert_eq!(id.to_string(), "Earth");
//! # Ok(()) }
//! ```
//!
//! # Error handling
//! Parsing invalid or out‑of‑range integers returns a descriptive `ErrorId`
//! variant. Errors are granular enough to distinguish whether a number failed
//! to match a planetary barycenter, a mass center, or a satellite mass center.
//!
//! # Design notes
//! * The public `NaifIds` enum is a *sum type* that preserves the category
//!   of the ID, which is often semantically useful downstream (e.g. when
//!   formatting, routing to different kernels, or resolving ephemerides).
//! * Submodules own their respective mapping tables and their own error types;
//!   this module lifts them into a unified error surface (`ErrorId`).
//!
//! # See also
//! -------------
//! * [`planet_bary`] — Planetary barycenter mapping table.
//! * [`planet_mass`] — Planet mass‑center mapping table (e.g. 199, 299, …).
//! * [`satellite_mass`] — Natural satellite mapping table (e.g. 301, 401, …).
//! * [`solar_system_bary`] — SSB (0) and Sun (10) codes.
//! * [`naif_type`] — Additional NAIF‑related type definitions.

/// Planetary barycenter identifiers (e.g., Mercury barycenter = 1).
pub mod planet_bary;

/// Planet mass center identifiers (e.g., Earth = 399, Mercury = 199).
pub mod planet_mass;

/// Natural satellite mass center identifiers (e.g., Moon = 301, Phobos = 401).
pub mod satellite_mass;

/// Solar System barycenter and Sun identifiers (0 and 10).
pub mod solar_system_bary;

/// Shared NAIF type definitions and helper utilities.
pub mod naif_type;

use std::fmt;

use planet_bary::PlanetaryBary;
use planet_mass::PlanetMassCenter;
use satellite_mass::SatelliteMassCenter;
use solar_system_bary::SolarSystemBary;
use thiserror::Error;

/// Error variants produced when converting raw integers into NAIF identifiers.
///
/// These errors are specific enough to signal which *category* the conversion
/// failed to match, and include the offending integer.
///
/// /// See also
/// ------------
/// * [`NaifIds::from_id`] – Fallible constructor from `i32`.
/// * Submodule errors in [`planet_bary`], [`planet_mass`], [`satellite_mass`].
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

/// Discriminated union over the main NAIF ID families.
///
/// This type preserves the semantic category of the code, which is useful
/// for routing, formatting, or applying category‑specific logic downstream.
///
/// Examples
/// ----------
/// ```rust
/// # use outfit::jpl_ephem::naif::naif_ids::{NaifIds, solar_system_bary::SolarSystemBary};
/// assert_eq!(NaifIds::SSB(SolarSystemBary::SSB).to_string(), "Solar System Barycenter");
/// ```
///
/// /// See also
/// ------------
/// * [`NaifIds::from_id`] – Parse a raw integer into a typed ID.
/// * [`NaifIds::to_id`] – Convert back to the raw integer form.
/// * [`TryFrom<i32>`] and [`From<NaifIds> for i32`] – Idiomatic conversions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum NaifIds {
    SSB(SolarSystemBary),
    PB(PlanetaryBary),
    PMC(PlanetMassCenter),
    SMC(SatelliteMassCenter),
}

impl NaifIds {
    /// Build a [`NaifIds`] from a raw NAIF integer code.
    ///
    /// Arguments
    /// -----------------
    /// * `id`: The raw NAIF integer identifier.
    ///
    /// Return
    /// ----------
    /// * `Ok(NaifIds)` if the code matches one of the supported families.
    /// * `Err(ErrorId)` if the code does not map to a known identifier.
    ///
    /// Notes
    /// -------
    /// * Special‑case handling for `0` (SSB) and `10` (Sun).
    /// * For `1..=999`, this function attempts, in order:
    ///   planetary barycenter → planet mass center → satellite mass center.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// # use outfit::jpl_ephem::naif::naif_ids::{NaifIds, solar_system_bary::SolarSystemBary};
    /// let ssb = NaifIds::from_id(0).unwrap();
    /// assert!(matches!(ssb, NaifIds::SSB(SolarSystemBary::SSB)));
    /// assert!(NaifIds::from_id(1000).is_err());
    /// ```
    ///
    /// See also
    /// ------------
    /// * [`TryFrom<i32>`] – Idiomatic conversion from integers.
    /// * [`NaifIds::to_id`] – The inverse of this function.
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

    /// Convert a typed [`NaifIds`] back to its raw integer code.
    ///
    /// Arguments
    /// -----------------
    /// * `&self`: The typed NAIF identifier.
    ///
    /// Return
    /// ----------
    /// * The corresponding `i32` NAIF code.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// # use outfit::jpl_ephem::naif::naif_ids::{NaifIds, solar_system_bary::SolarSystemBary};
    /// let code = NaifIds::SSB(SolarSystemBary::Sun).to_id();
    /// assert_eq!(code, 10);
    /// ```
    ///
    /// /// See also
    /// ------------
    /// * [`From<NaifIds> for i32`] – Idiomatic conversion to integers.
    /// * [`NaifIds::from_id`] – Parsing from integers.
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
    /// Lossless conversion from a typed NAIF ID into its raw integer code.
    ///
    /// See also
    /// ------------
    /// * [`NaifIds::to_id`] – Underlying implementation.
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

    /// Fallible conversion from raw integer to typed NAIF ID.
    ///
    /// Arguments
    /// -----------------
    /// * `id`: The raw NAIF code to parse.
    ///
    /// Return
    /// ----------
    /// * `Result<NaifIds, ErrorId>` with a category‑aware error on failure.
    ///
    /// See also
    /// ------------
    /// * [`NaifIds::from_id`] – Same logic in inherent form.
    fn try_from(id: i32) -> Result<Self, Self::Error> {
        NaifIds::from_id(id)
    }
}

impl fmt::Display for NaifIds {
    /// Human‑readable label for NAIF IDs.
    ///
    /// Planet and satellite names come from the corresponding submodule tables.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// # use outfit::jpl_ephem::naif::naif_ids::{NaifIds, solar_system_bary::SolarSystemBary};
    /// assert_eq!(NaifIds::SSB(SolarSystemBary::SSB).to_string(), "Solar System Barycenter");
    /// ```
    ///
    /// See also
    /// ------------
    /// * [`planet_bary`], [`planet_mass`], [`satellite_mass`].
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
