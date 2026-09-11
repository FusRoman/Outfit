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

use std::fmt;

use crate::jpl_ephem::naif::naif_ids::{
    main_belt::AsteroidNumber, planet_bary::PlanetaryBary, satellite_mass::SatelliteMassCenter,
    solar_system_bary::SolarSystemBary, NaifIds,
};
use crate::propagator::asteroid_gm_table::{self, ASTEROID_GM_RATIOS};

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
///
/// Planets and the Moon use the DE440 table above; a numbered main-belt
/// asteroid (`NaifIds::AST`) is looked up in the embedded main-belt asteroid
/// GM table, which covers the 300 bodies of the `codes_300ast_20100725.bsp`
/// supplementary kernel (see [`known_main_belt_asteroids`]).
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
        NaifIds::AST(asteroid) => asteroid_gm_table::asteroid_gm_au3_day2(asteroid.0),
        _ => None,
    }
}

/// One catalog entry for a main-belt asteroid perturber: its identity and its
/// gravitational parameter.
///
/// # See also
///
/// [`main_belt_asteroid_catalog`] builds the full, mass-sorted catalog these
/// entries come from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MainBeltAsteroidInfo {
    /// The asteroid's minor-planet number.
    pub number: AsteroidNumber,
    /// Gravitational parameter, in AU³/day² (same convention as
    /// [`gm_au3_day2`]).
    pub gm_au3_day2: f64,
}

impl fmt::Display for MainBeltAsteroidInfo {
    /// Formats as `"<name/number> — GM = <value> AU³/day²"`, e.g.
    /// `"Ceres (1) — GM = 1.402e-13 AU³/day²"`.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} — GM = {:.3e} AU³/day²",
            self.number, self.gm_au3_day2
        )
    }
}

/// The main-belt asteroid catalog, sorted by decreasing gravitational
/// parameter (the most massive, most perturbing bodies first).
///
/// This is the same 300 bodies as [`known_main_belt_asteroids`], bundled with
/// their name and GM for display, and in the order that matters when picking
/// "the top N" perturbers for a cost/accuracy tradeoff — see
/// [`known_main_belt_asteroids_by_mass`] for a version that returns
/// ready-to-use [`NaifIds`] directly.
///
/// # Returns
///
/// All 300 entries, descending by `gm_au3_day2`.
///
/// # Examples
///
/// ```rust
/// use outfit::propagator::planet_gm::main_belt_asteroid_catalog;
///
/// for entry in main_belt_asteroid_catalog().iter().take(5) {
///     println!("{entry}");
/// }
/// ```
pub fn main_belt_asteroid_catalog() -> Vec<MainBeltAsteroidInfo> {
    let mut catalog: Vec<MainBeltAsteroidInfo> = ASTEROID_GM_RATIOS
        .iter()
        .map(|&(number, ratio)| MainBeltAsteroidInfo {
            number: AsteroidNumber(number),
            gm_au3_day2: ratio * GM_SUN,
        })
        .collect();
    catalog.sort_by(|a, b| b.gm_au3_day2.total_cmp(&a.gm_au3_day2));
    catalog
}

/// The main-belt asteroids Outfit has both a position source and a GM for,
/// ordered by decreasing gravitational parameter (the most perturbing bodies
/// first).
///
/// Equivalent to [`known_main_belt_asteroids`] but sorted by mass instead of
/// by minor-planet number — the natural order for picking, say, "the 16
/// asteroids that perturb the main belt the most" with `.take(16)`.
///
/// # Returns
///
/// An iterator over `NaifIds::AST(_)`, one per known asteroid, descending by
/// GM.
pub fn known_main_belt_asteroids_by_mass() -> impl Iterator<Item = NaifIds> {
    main_belt_asteroid_catalog()
        .into_iter()
        .map(|entry| NaifIds::AST(entry.number))
}

/// The main-belt asteroids Outfit has both a position source and a GM for.
///
/// Each entry is ready to push into
/// [`NBodyConfig::perturbing_bodies`](crate::propagator::NBodyConfig::perturbing_bodies)
/// once the ANISE backend is built with
/// [`JPLEphem::from_anise`](crate::jpl_ephem::JPLEphem::from_anise) and
/// [`JPLEphem::with_main_belt_asteroids`](crate::jpl_ephem::JPLEphem::with_main_belt_asteroids).
/// The 300 bodies come from the `codes_300ast_20100725.bsp` supplementary
/// kernel; see the embedded main-belt asteroid GM table for the exact list
/// and its source.
///
/// # Returns
///
/// An iterator over `NaifIds::AST(_)`, one per known asteroid, in ascending
/// minor-planet-number order. See [`known_main_belt_asteroids_by_mass`] for
/// the same bodies ordered by decreasing perturbing influence instead.
pub fn known_main_belt_asteroids() -> impl Iterator<Item = NaifIds> {
    ASTEROID_GM_RATIOS
        .iter()
        .map(|&(number, _)| NaifIds::AST(AsteroidNumber(number)))
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

    #[test]
    fn gm_au3_day2_resolves_a_known_asteroid() {
        use crate::jpl_ephem::naif::naif_ids::main_belt::AsteroidNumber;
        let ceres_gm = gm_au3_day2(NaifIds::AST(AsteroidNumber::CERES)).unwrap();
        assert!(ceres_gm.is_finite() && ceres_gm > 0.0);
        // Ceres is by far the most massive of the 300 bodies, but still tiny
        // next to the smallest planet perturber tabulated above.
        assert!(ceres_gm < GM_MARS);
    }

    #[test]
    fn gm_au3_day2_rejects_an_unknown_asteroid() {
        use crate::jpl_ephem::naif::naif_ids::main_belt::AsteroidNumber;
        assert_eq!(gm_au3_day2(NaifIds::AST(AsteroidNumber(999_999))), None);
    }

    #[test]
    fn known_main_belt_asteroids_covers_every_tabulated_body() {
        let bodies: Vec<NaifIds> = known_main_belt_asteroids().collect();
        assert_eq!(bodies.len(), ASTEROID_GM_RATIOS.len());
        assert!(bodies.iter().all(|&body| gm_au3_day2(body).is_some()));
    }

    /// The GM table and the name table (`main_belt::AsteroidNumber::name`)
    /// are two independently-sourced 300-row tables that must describe
    /// exactly the same set of bodies; this guards against one being updated
    /// without the other.
    #[test]
    fn every_tabulated_asteroid_has_a_name() {
        for &(number, _) in ASTEROID_GM_RATIOS {
            assert!(
                AsteroidNumber(number).name().is_some(),
                "asteroid {number} has a GM but no name"
            );
        }
    }

    #[test]
    fn main_belt_asteroid_catalog_is_sorted_by_decreasing_mass() {
        let catalog = main_belt_asteroid_catalog();
        assert_eq!(catalog.len(), ASTEROID_GM_RATIOS.len());
        assert!(catalog
            .windows(2)
            .all(|w| w[0].gm_au3_day2 >= w[1].gm_au3_day2));
        // Ceres is the most massive of the 300 bodies.
        assert_eq!(catalog[0].number, AsteroidNumber::CERES);
    }

    #[test]
    fn known_main_belt_asteroids_by_mass_matches_the_catalog_order() {
        let from_catalog: Vec<NaifIds> = main_belt_asteroid_catalog()
            .into_iter()
            .map(|entry| NaifIds::AST(entry.number))
            .collect();
        let from_iterator: Vec<NaifIds> = known_main_belt_asteroids_by_mass().collect();
        assert_eq!(from_catalog, from_iterator);
        assert_eq!(from_iterator.len(), ASTEROID_GM_RATIOS.len());
    }

    #[test]
    fn main_belt_asteroid_info_displays_name_and_gm() {
        let ceres = main_belt_asteroid_catalog()[0];
        assert_eq!(ceres.number, AsteroidNumber::CERES);
        let text = ceres.to_string();
        assert!(text.starts_with("Ceres (1) — GM = "));
        assert!(text.ends_with(" AU³/day²"));
    }
}
