//! Integration test: main-belt asteroids as N-body perturbers (ANISE backend).
//!
//! Loads the supplementary asteroid kernel through
//! [`JPLEphem::from_anise_with_main_belt_asteroids`] and checks that adding
//! Ceres, Vesta and Pallas to `NBodyConfig::perturbing_bodies` perturbs the
//! propagated orbit away from the Sun-only solution by a small, finite
//! amount — the physical-bounds style used elsewhere for N-body
//! non-regression (see `tests/test_diff_cor.rs`'s `test_diff_cor_nbody`)
//! rather than a tight numeric oracle, since the exact perturbation size
//! depends on the geometry of the day.
#![cfg(feature = "ephem-anise")]

use nalgebra::Vector3;
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::jpl_ephem::naif::naif_ids::main_belt::AsteroidNumber;
use outfit::jpl_ephem::naif::naif_ids::{solar_system_bary::SolarSystemBary, NaifIds};
use outfit::orbit_type::equinoctial_element::EquinoctialElements;
use outfit::propagator::NBodyConfig;
use outfit::JPLEphem;

/// Load DE440 with the main-belt asteroid supplementary kernel through ANISE.
fn load_ephem_with_asteroids() -> JPLEphem {
    let source: EphemFileSource = "naif:DE440"
        .try_into()
        .expect("failed to parse JPL ephemeris source");
    JPLEphem::from_anise_with_main_belt_asteroids(&source)
        .expect("failed to load DE440 with the main-belt asteroid supplementary kernel")
}

/// A main-belt-like test orbit whose reference epoch is `epoch` (MJD TT).
fn test_orbit(epoch: f64) -> EquinoctialElements {
    EquinoctialElements {
        reference_epoch: epoch,
        semi_major_axis: 2.7,
        eccentricity_sin_lon: 0.10,
        eccentricity_cos_lon: 0.08,
        tan_half_incl_sin_node: 0.030,
        tan_half_incl_cos_node: 0.082,
        mean_longitude: 0.5,
    }
}

fn sun_only_config() -> NBodyConfig {
    NBodyConfig {
        perturbing_bodies: vec![NaifIds::SSB(SolarSystemBary::Sun)],
        ..NBodyConfig::default()
    }
}

fn sun_and_asteroids_config() -> NBodyConfig {
    NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::AST(AsteroidNumber::CERES),
            NaifIds::AST(AsteroidNumber::VESTA),
            NaifIds::AST(AsteroidNumber::PALLAS),
        ],
        ..NBodyConfig::default()
    }
}

/// Maximum absolute component error between two vectors.
fn max_abs_diff(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    (a - b).abs().max()
}

#[test]
fn asteroid_perturbers_move_the_propagated_orbit_by_a_small_finite_amount() {
    let jpl = load_ephem_with_asteroids();
    let t0 = 60_000.0_f64;
    let elements = test_orbit(t0);

    let sun_only = sun_only_config();
    let sun_and_asteroids = sun_and_asteroids_config();

    for span in [90.0_f64, -90.0] {
        let baseline = elements
            .propagate_nbody(t0 + span, &jpl, &sun_only, None)
            .expect("sun-only n-body propagation failed");
        let perturbed = elements
            .propagate_nbody(t0 + span, &jpl, &sun_and_asteroids, None)
            .expect("n-body propagation with asteroid perturbers failed");

        assert!(
            baseline.position.iter().all(|c| c.is_finite())
                && perturbed.position.iter().all(|c| c.is_finite()),
            "span {span} d: propagated position must be finite"
        );

        let pos_shift = max_abs_diff(baseline.position, perturbed.position);
        // The three asteroids are far less massive than any planet, so their
        // effect over a ~90-day arc is small; it must still be strictly
        // positive (the perturbers are actually applied) and stay well below
        // a planet-sized perturbation (> 1e-3 AU would indicate something is
        // off, e.g. a misplaced body or a unit error).
        assert!(
            pos_shift > 0.0 && pos_shift.is_finite(),
            "span {span} d: asteroid perturbers had no measurable effect"
        );
        assert!(
            pos_shift < 1e-3,
            "span {span} d: asteroid-perturber shift {pos_shift:e} AU is implausibly large"
        );
    }
}

#[test]
fn known_main_belt_asteroids_are_all_usable_as_perturbers() {
    use outfit::propagator::planet_gm::known_main_belt_asteroids;

    let jpl = load_ephem_with_asteroids();
    let t0 = 60_000.0_f64;
    let elements = test_orbit(t0);

    // Build one config from a subset of the full catalog (mirrors how a user
    // would pick "most" of the main-belt perturbers) and check the
    // propagation succeeds end to end.
    let config = NBodyConfig {
        perturbing_bodies: std::iter::once(NaifIds::SSB(SolarSystemBary::Sun))
            .chain(known_main_belt_asteroids().take(16))
            .collect(),
        ..NBodyConfig::default()
    };

    let result = elements
        .propagate_nbody(t0 + 30.0, &jpl, &config, None)
        .expect("n-body propagation with 16 asteroid perturbers failed");
    assert!(result.position.iter().all(|c| c.is_finite()));
}
