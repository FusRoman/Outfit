use crate::observers::observers::Observer;
use std::collections::HashMap;

pub const EPS: f64 = 1e-6;
pub const T2000: f64 = 51544.5; // J2000 Epoch for MJD
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

/// type def
pub type Degree = f64;
pub type Kilometer = f64;
pub type MpcCodeObs = HashMap<String, Observer>;
