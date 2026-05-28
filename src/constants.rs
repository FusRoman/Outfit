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

// -------------------------------------------------------------------------------------------------
// Physical constants and unit conversions
// -------------------------------------------------------------------------------------------------

use std::collections::HashMap;

use ahash::RandomState;
use nalgebra::{Matrix3, Vector3};
use photom::TrajId;

use crate::{GaussResult, OrbitalElements, OutfitError};

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

// Angular velocity of Earth rotation (rad/day) on the z-axis.
pub const EARTH_ROTATION: Vector3<f64> = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

// Hard coded rotation matrices for coordinate transformations between mean equatorial J2000 and mean ecliptic J2000 frames.
// Can be computed using the rotpn function in the ref_system module

/// Rotation matrix from mean equatorial J2000 to mean ecliptic J2000.
///
/// Rotation of $-\varepsilon$ around the X-axis, where $\varepsilon$ is the
/// obliquity of the ecliptic at J2000.
///
/// Equivalent to `rotpn(RefSystem::Equm(RefEpoch::J2000), RefSystem::Eclm(RefEpoch::J2000))`.
pub const ROT_EQUMJ2000_TO_ECLMJ2000: Matrix3<f64> = Matrix3::new(
    1.0e0,
    0.0e0,
    0.0e0,
    0.0e0,
    9.174_820_620_691_818e-1,
    3.977_771_559_319_137e-1,
    0.0e0,
    -3.977_771_559_319_137e-1,
    9.174_820_620_691_818e-1,
);

/// Rotation matrix from mean ecliptic J2000 to mean equatorial J2000.
///
/// Rotation of $+\varepsilon$ around the X-axis (transpose / inverse of
/// [`ROT_EQUMJ2000_TO_ECLMJ2000`]).
///
/// Equivalent to `rotpn(RefSystem::Eclm(RefEpoch::J2000), RefSystem::Equm(RefEpoch::J2000))`.
pub const ROT_ECLMJ2000_TO_EQUMJ2000: Matrix3<f64> = Matrix3::new(
    1.0e0,
    0.0e0,
    0.0e0,
    0.0e0,
    9.174_820_620_691_818e-1,
    -3.977_771_559_319_137e-1,
    0.0e0,
    3.977_771_559_319_137e-1,
    9.174_820_620_691_818e-1,
);

/// Modified Julian Date (Scale Ephemeris Time, ET)
pub type MJDET = f64;

/// Type alias for the RMS of normalized residuals from an IOD fit.
/// This is a single scalar value representing the overall fit quality of the IOD solution.
pub type IODRMS = f64;

/// Type alias for the chi-squared value of a fit, used in differential correction.
pub type Chi2 = f64;

pub enum FitOrbitResult {
    IODGauss((GaussResult, IODRMS)),
    DifferentialCorrection((OrbitalElements, Chi2)),
}

impl FitOrbitResult {
    /// Returns a scalar measure of orbit quality for this fit result.
    /// For IODGauss, this is the RMS of normalized residuals; for DifferentialCorrection, this is the chi-squared value.
    ///
    /// # Returns
    /// - `f64` — a single scalar representing the fit quality.
    pub fn orbit_quality(&self) -> f64 {
        match self {
            FitOrbitResult::IODGauss((_, rms)) => *rms,
            FitOrbitResult::DifferentialCorrection((_, chi2)) => *chi2,
        }
    }

    /// Returns a reference to the orbital elements associated with this fit result.
    /// For `IODGauss`, this extracts the orbital elements from the `GaussResult`; for `DifferentialCorrection`, it returns the orbital elements directly.
    ///
    /// # Returns
    /// - `&OrbitalElements` — a reference to the orbital elements of the fitted orbit.
    pub fn orbital_elements(&self) -> &OrbitalElements {
        match self {
            FitOrbitResult::IODGauss((gauss_result, _)) => gauss_result.get_orbit(),
            FitOrbitResult::DifferentialCorrection((orbital_elements, _)) => orbital_elements,
        }
    }
}

/// Full batch orbit determination results.
///
/// Each entry maps an [`TrajId`] to the outcome of a full
/// Initial Orbit Determination (IOD) attempt on its set of observations.
///
/// Internally, this is implemented as:
///
/// ```ignore
/// HashMap<TrajId, Result<FitOrbitResult, OutfitError>, RandomState>
/// ```
///
/// Return semantics
/// -----------------
/// * `Ok(FitOrbitResult::IODGauss((GaussResult, IODRMS)))` – a successful IOD with its RMS of normalized residuals.
/// * `Ok(FitOrbitResult::DifferentialCorrection((OrbitalElements, Chi2)))` – a successful differential correction with its chi-squared value.
/// * `Err(OutfitError)` – a failure isolated to that object.
///
/// Use RandomState from the ahash crate for efficient hashing of TrajId keys.
pub type FullOrbitResult = HashMap<TrajId, Result<FitOrbitResult, OutfitError>, RandomState>;
