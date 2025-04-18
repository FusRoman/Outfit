use crate::observations::observations::Observation;
use crate::observers::observers::Observer;
use std::collections::HashMap;
use std::sync::Arc;
use ahash::RandomState;
use smallvec::SmallVec;

use crate::initial_orbit_determination::gauss::GaussObs;

pub const EPS: f64 = 1e-6;
pub const T2000: f64 = 51544.5; // J2000 Epoch for MJD
pub const JDTOMJD: f64 = 2400000.5; // Julian Date -> Modified Julian Date conversion factor
pub const RADEG: f64 = std::f64::consts::PI / 180.0; // Degrees -> radians conversion factor
pub const RADSEC: f64 = std::f64::consts::PI / 648000.0; // arcsecond -> radians conversion factor

// Constants
pub const DPI: f64 = 2. * std::f64::consts::PI;
pub const SECONDS_PER_DAY: f64 = 86_400.0;
pub const AU: f64 = 149_597_870.7;

/// Earth ellipsoid constant from (GRS1980/WGS84)
pub const EARTH_MAJOR_AXIS: f64 = 6_378_137.0; // Earth equatorial radius in meter
pub const EARTH_MINOR_AXIS: f64 = 6_356_752.3; // Earth polar radius in meter
pub const ERAU: f64 = (EARTH_MAJOR_AXIS / 1000.) / AU;

/// Gaussian gravitational constant
pub const GAUSS_GRAV: f64 = 0.01720209895;
pub const GAUSS_GRAV_SQUARED: f64 = GAUSS_GRAV * GAUSS_GRAV;

///speed of light (Km/s)
pub const VLIGHT: f64 = 2.99792458e5;
/// speed of light (AU/day)
pub const VLIGHT_AU: f64 = VLIGHT / AU * SECONDS_PER_DAY;

// type def
pub type Degree = f64;
pub type Kilometer = f64;
/// a mpc code observatory made of three characters
pub type MpcCode = String;
pub type MpcCodeObs = HashMap<MpcCode, Arc<Observer>>;
/// Modified Julian Date
pub type MJD = f64; // time in modified julian date

// Type definitions for the trajectories and observations

/// An object number if either an asteroid number (e.g. 1234)
/// or a comet number (e.g. 1234P) or a provisional designation (e.g. K25D50B)
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum ObjectNumber {
    Int(u32),
    String(String)
}

pub type Observations = SmallVec<[Observation; 6]>;

/// A set of trajectories
/// The key is the object number
/// The value is a vector of observations associated to this object
pub type TrajectorySet = HashMap<ObjectNumber, Observations, RandomState>;

pub type Triplets = HashMap<ObjectNumber, Vec<GaussObs>>;