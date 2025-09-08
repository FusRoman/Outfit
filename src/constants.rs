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
use ahash::RandomState;
use smallvec::SmallVec;
use std::collections::HashMap;
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
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ObjectNumber {
    /// Integer-based MPC designation (e.g. 1, 433…)
    Int(u32),
    /// String-based designation (provisional, comet, etc.)
    String(String),
}

/// A small, inline-optimized container for observations of a single object.
pub type Observations = SmallVec<[Observation; 6]>;

/// A full set of trajectories for multiple objects.
///
/// The key is the [`ObjectNumber`] (identifier of an object).
/// The value is the list of [`Observation`] associated with this object.
///
/// Uses [`ahash`](https://docs.rs/ahash) for fast hashing.
pub type TrajectorySet = HashMap<ObjectNumber, Observations, RandomState>;
