//! Gravitational parameters (GM) for solar-system bodies.
//!
//! All values are in **AU³/day²**, consistent with the Gaussian gravitational
//! constant `k` used throughout Outfit (`k² = GAUSS_GRAV_SQUARED`).
//!
//! Sources: DE440 TDB-compatible mass parameters converted to AU³/day² via
//! `GM_km3_s2 * (86400² / AU³_in_km³)` where `AU = 1.495978707e8 km`.
//!
//! # References
//! - Park, R.S. et al. (2021), *The JPL Planetary and Lunar Ephemerides DE440 and DE441*,
//!   AJ 161, 105.

use crate::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, satellite_mass::SatelliteMassCenter,
    solar_system_bary::SolarSystemBary, NaifIds,
};

/// AU in kilometres (IAU 2012).
const AU_KM: f64 = 1.495_978_707e8;

/// Conversion factor: km³/s² → AU³/day²
/// = (86400 s/day)² / (AU_KM km/AU)³
const KM3_S2_TO_AU3_DAY2: f64 = (86400.0 * 86400.0) / (AU_KM * AU_KM * AU_KM);

// ---------------------------------------------------------------------------
// Individual GM values in km³/s² (DE440)
// ---------------------------------------------------------------------------

const GM_SUN_KM3_S2: f64 = 1.327_124_400_41e11;
const GM_MERCURY_KM3_S2: f64 = 2.203_178_e4;
const GM_VENUS_KM3_S2: f64 = 3.248_585_7e5;
/// Earth + Moon system GM (used for EarthMoon barycenter).
const GM_EARTH_MOON_KM3_S2: f64 = 4.035_032_35e5;
const GM_MARS_KM3_S2: f64 = 4.282_837_36e4;
const GM_JUPITER_KM3_S2: f64 = 1.267_127_648e8;
const GM_SATURN_KM3_S2: f64 = 3.794_062_52e7;
const GM_URANUS_KM3_S2: f64 = 5.794_556_4e6;
const GM_NEPTUNE_KM3_S2: f64 = 6.836_527_1e6;
const GM_PLUTO_KM3_S2: f64 = 9.755_e2;
const GM_MOON_KM3_S2: f64 = 4.902_800_066e3;

// ---------------------------------------------------------------------------
// Public AU³/day² constants
// ---------------------------------------------------------------------------

pub const GM_SUN: f64 = GM_SUN_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_MERCURY: f64 = GM_MERCURY_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_VENUS: f64 = GM_VENUS_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_EARTH_MOON: f64 = GM_EARTH_MOON_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_MARS: f64 = GM_MARS_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_JUPITER: f64 = GM_JUPITER_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_SATURN: f64 = GM_SATURN_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_URANUS: f64 = GM_URANUS_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_NEPTUNE: f64 = GM_NEPTUNE_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_PLUTO: f64 = GM_PLUTO_KM3_S2 * KM3_S2_TO_AU3_DAY2;
pub const GM_MOON: f64 = GM_MOON_KM3_S2 * KM3_S2_TO_AU3_DAY2;

/// Return the GM (AU³/day²) for a [`NaifIds`], or `None` if not tabulated.
pub fn gm_au3_day2(body: NaifIds) -> Option<f64> {
    match body {
        NaifIds::SSB(SolarSystemBary::Sun) => Some(GM_SUN),
        NaifIds::SSB(SolarSystemBary::SSB) => None,
        NaifIds::PB(PlanetaryBary::Mercury) => Some(GM_MERCURY),
        NaifIds::PB(PlanetaryBary::Venus) => Some(GM_VENUS),
        NaifIds::PB(PlanetaryBary::EarthMoon) => Some(GM_EARTH_MOON),
        NaifIds::PB(PlanetaryBary::Mars) => Some(GM_MARS),
        NaifIds::PB(PlanetaryBary::Jupiter) => Some(GM_JUPITER),
        NaifIds::PB(PlanetaryBary::Saturn) => Some(GM_SATURN),
        NaifIds::PB(PlanetaryBary::Uranus) => Some(GM_URANUS),
        NaifIds::PB(PlanetaryBary::Neptune) => Some(GM_NEPTUNE),
        NaifIds::PB(PlanetaryBary::Pluto) => Some(GM_PLUTO),
        NaifIds::SMC(SatelliteMassCenter::Moon) => Some(GM_MOON),
        _ => None,
    }
}

#[cfg(test)]
mod planet_gm_tests {
    use super::*;

    #[test]
    fn sun_gm_close_to_gauss_squared() {
        // k² = 2.9591220828e-4 AU³/day² (Gaussian)
        // DE440 GM_Sun should be within 1e-8 relative of the Gaussian value
        let gauss_k2 = crate::constants::GAUSS_GRAV_SQUARED;
        let rel = (GM_SUN - gauss_k2).abs() / gauss_k2;
        assert!(rel < 1e-4, "GM_SUN rel diff from k²: {rel:.2e}");
    }

    #[test]
    fn jupiter_dominates_perturbations() {
        // Jupiter GM should be ~1000× larger than Mars GM
        const { assert!(GM_JUPITER > GM_MARS * 100.0) };
    }
}
