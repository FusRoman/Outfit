//! # Constants and type definitions for Outfit
//!
//! This module centralizes the **physical constants**, **conversion factors**, and **common type
//! definitions** used throughout the `Outfit` library. It also defines key data structures for
//! organizing observations and trajectories.
//!
//! ## Overview
//!
//! - Astronomical and geophysical constants
//! - Unit conversions (degrees ↔ radians, days ↔ seconds, AU ↔ km)
//! - Core type aliases used across the crate
//! - Identifiers for solar system bodies
//! - Container types for storing observations and trajectories
//!
//! These definitions are used by all main modules, including orbit determination, observers,
//! and ephemerides.

use crate::observations::Observation;
use crate::observers::Observer;
use smallvec::SmallVec;
use std::borrow::Cow;
use std::collections::HashMap;
use std::convert::TryFrom;
use std::sync::Arc;

// -------------------------------------------------------------------------------------------------
// Physical constants and unit conversions
// -------------------------------------------------------------------------------------------------

/// 2π, useful for trigonometric conversions
pub const DPI: f64 = 2. * std::f64::consts::PI;

/// Number of seconds in a Julian day
pub const SECONDS_PER_DAY: f64 = 86_400.0;

/// Astronomical Unit in kilometers (IAU 2012)
pub const AU: f64 = 149_597_870.7;

/// Numerical epsilon used for floating-point comparisons
pub const EPS: f64 = 1e-6;

/// MJD epoch of J2000.0 (2000-01-01 12:00:00 TT)
pub const T2000: f64 = 51544.5;

/// Conversion factor between Julian Date and Modified Julian Date
pub const JDTOMJD: f64 = 2400000.5;

/// Degrees → radians
pub const RADEG: f64 = std::f64::consts::PI / 180.0;

/// Arcseconds → radians
pub const RADSEC: f64 = std::f64::consts::PI / 648000.0;

/// Radians → arcseconds
pub const RAD2ARC: f64 = 648000.0 / std::f64::consts::PI;

/// Hours → radians
pub const RADH: f64 = DPI / 24.0;

/// Earth equatorial radius in meters (GRS1980/WGS84)
pub const EARTH_MAJOR_AXIS: f64 = 6_378_137.0;

/// Earth polar radius in meters (GRS1980/WGS84)
pub const EARTH_MINOR_AXIS: f64 = 6_356_752.3;

/// Earth radius expressed in astronomical units
pub const ERAU: f64 = (EARTH_MAJOR_AXIS / 1000.) / AU;

/// Gaussian gravitational constant k (used in classical orbit dynamics)
pub const GAUSS_GRAV: f64 = 0.01720209895;

/// k², often used in Kepler’s third law
pub const GAUSS_GRAV_SQUARED: f64 = GAUSS_GRAV * GAUSS_GRAV;

/// Speed of light in km/s
pub const VLIGHT: f64 = 2.99792458e5;

/// Speed of light in astronomical units per day
pub const VLIGHT_AU: f64 = VLIGHT / AU * SECONDS_PER_DAY;

// -------------------------------------------------------------------------------------------------
// Type aliases
// -------------------------------------------------------------------------------------------------

/// Angle in degrees
pub type Degree = f64;
/// Angle in arcseconds
pub type ArcSec = f64;
/// Angle in radians
pub type Radian = f64;
/// Distance in kilometers
pub type Kilometer = f64;
/// Distance in meters
pub type Meter = f64;
/// MPC code identifying an observatory (3 characters)
pub type MpcCode = String;

/// Lookup table from MPC code to [`Observer`] metadata
pub type MpcCodeObs = HashMap<MpcCode, Arc<Observer>>;

/// Modified Julian Date (days)
pub type MJD = f64;

// -------------------------------------------------------------------------------------------------
// Identifiers and data containers
// -------------------------------------------------------------------------------------------------

/// Identifier of a solar system object.
///
/// This can be:
/// - An asteroid number (e.g. `Int(1234)`)
/// - A comet number (e.g. `"1234P"`)
/// - A provisional designation (e.g. `"K25D50B"`)
#[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum ObjectNumber {
    /// Integer-based MPC designation (e.g. 1, 433…)
    Int(u32),
    /// String-based designation (provisional, comet, etc.)
    String(String),
}

impl std::fmt::Display for ObjectNumber {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ObjectNumber::Int(n) => write!(f, "{n}"),
            ObjectNumber::String(s) => write!(f, "{s}"),
        }
    }
}

// --- Infallible conversions (enable `.into()` directly) ----------------------

impl From<u32> for ObjectNumber {
    #[inline]
    fn from(n: u32) -> Self {
        ObjectNumber::Int(n)
    }
}

impl From<u16> for ObjectNumber {
    #[inline]
    fn from(n: u16) -> Self {
        ObjectNumber::Int(n as u32)
    }
}

impl From<u8> for ObjectNumber {
    #[inline]
    fn from(n: u8) -> Self {
        ObjectNumber::Int(n as u32)
    }
}

impl From<&u32> for ObjectNumber {
    /// Convenience to allow `(&n).into()` without dereferencing at call sites.
    #[inline]
    fn from(n: &u32) -> Self {
        ObjectNumber::Int(*n)
    }
}

impl From<String> for ObjectNumber {
    #[inline]
    fn from(s: String) -> Self {
        ObjectNumber::String(s)
    }
}

impl From<&String> for ObjectNumber {
    /// Clones the string to build a `String`-backed identifier.
    #[inline]
    fn from(s: &String) -> Self {
        ObjectNumber::String(s.clone())
    }
}

impl From<&str> for ObjectNumber {
    /// Note: this **does not** parse numeric strings into `Int`. Use `FromStr` if you want
    /// `"1234"` to become `ObjectNumber::Int(1234)`.
    #[inline]
    fn from(s: &str) -> Self {
        ObjectNumber::String(s.to_string())
    }
}

impl<'a> From<Cow<'a, str>> for ObjectNumber {
    /// Accept both borrowed and owned `Cow<str>`.
    #[inline]
    fn from(c: Cow<'a, str>) -> Self {
        match c {
            Cow::Borrowed(s) => ObjectNumber::String(s.to_string()),
            Cow::Owned(s) => ObjectNumber::String(s),
        }
    }
}

// --- Fallible conversions (use `.try_into()` to be overflow-safe) ------------

impl TryFrom<usize> for ObjectNumber {
    type Error = std::num::TryFromIntError;

    /// Convert a `usize` into `Int(u32)` if it fits.
    #[inline]
    fn try_from(n: usize) -> Result<Self, Self::Error> {
        Ok(ObjectNumber::Int(u32::try_from(n)?))
    }
}

impl TryFrom<u64> for ObjectNumber {
    type Error = std::num::TryFromIntError;

    /// Convert a `u64` into `Int(u32)` if it fits.
    #[inline]
    fn try_from(n: u64) -> Result<Self, Self::Error> {
        Ok(ObjectNumber::Int(u32::try_from(n)?))
    }
}

impl TryFrom<i64> for ObjectNumber {
    type Error = &'static str;

    /// Convert a non-negative `i64` into `Int(u32)` if it fits.
    #[inline]
    fn try_from(n: i64) -> Result<Self, Self::Error> {
        if n < 0 {
            return Err("negative value is not a valid ObjectNumber::Int");
        }
        let n = u64::try_from(n).map_err(|_| "conversion failed")?;
        let n = u32::try_from(n).map_err(|_| "value exceeds u32 range")?;
        Ok(ObjectNumber::Int(n))
    }
}

// --- Smart parsing from &str via `FromStr` (optional) ------------------------

impl std::str::FromStr for ObjectNumber {
    type Err = std::num::ParseIntError;

    /// Try to parse an `ObjectNumber` from a string.
    ///
    /// Rules
    /// -----
    /// - Pure digits that fit in `u32` → `Int(u32)`.
    /// - Otherwise                         → `String(String)`.
    ///
    /// Note
    /// ----
    /// If the string is *only* digits but **does not** fit in `u32`, this returns the
    /// original `ParseIntError`. If you prefer to always fallback to `String` on
    /// overflow, we can change the policy (but it’s usually better to fail loudly).
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.parse::<u32>() {
            Ok(n) => Ok(ObjectNumber::Int(n)),
            Err(e) => {
                if s.chars().any(|c| !c.is_ascii_digit()) {
                    Ok(ObjectNumber::String(s.to_string()))
                } else {
                    Err(e)
                }
            }
        }
    }
}

/// A small, inline-optimized container for observations of a single object.
pub type Observations = SmallVec<[Observation; 6]>;
