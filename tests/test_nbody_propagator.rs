//! Non-regression guardrail for the N-body integrator.
//!
//! With the Sun as the only perturbing body the N-body DOP853 propagation
//! reduces to a numerically solved two-body problem, so it must reproduce the
//! analytic Keplerian solution to near machine precision.  This pins the
//! integrator, the augmented-state layout, the frame conversions and the new
//! time-dependent perturber path: the Sun is special-cased to an exactly-zero
//! heliocentric position, so switching to the Chebyshev interpolation must not
//! move this result.

mod common;

use nalgebra::Vector3;
use outfit::jpl_ephem::naif::naif_ids::{solar_system_bary::SolarSystemBary, NaifIds};
use outfit::orbit_type::equinoctial_element::EquinoctialElements;
use outfit::propagator::perturber_ephemeris::PerturberEphemerisSet;
use outfit::propagator::NBodyConfig;
use outfit::JPLEphem;

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

fn load_ephem() -> JPLEphem {
    common::load_ephem()
}

/// Maximum absolute component error between two vectors.
fn max_abs_diff(a: Vector3<f64>, b: Vector3<f64>) -> f64 {
    (a - b).abs().max()
}

#[test]
fn sun_only_nbody_matches_analytic_two_body() {
    let jpl = load_ephem();
    let t0 = 60_000.0_f64;
    let elements = test_orbit(t0);
    let config = sun_only_config();

    // Forward and backward spans, including a long backward arc like the
    // ground-truth propagation in the fitting pipeline.
    for span in [30.0_f64, 90.0, 200.0, -90.0, -200.0] {
        let (pos_2b, vel_2b, _) = elements
            .propagate_twobody(0.0, span, false)
            .expect("two-body propagation failed");

        let nbody = elements
            .propagate_nbody(t0 + span, &jpl, &config, None)
            .expect("n-body propagation failed");

        let pos_err = max_abs_diff(pos_2b, nbody.position);
        let vel_err = max_abs_diff(vel_2b, nbody.velocity);

        assert!(
            pos_err < 1e-9,
            "span {span} d: position mismatch {pos_err:e} AU exceeds 1e-9"
        );
        assert!(
            vel_err < 1e-11,
            "span {span} d: velocity mismatch {vel_err:e} AU/day exceeds 1e-11"
        );
    }
}

#[test]
fn precomputed_table_matches_on_the_fly_build() {
    let jpl = load_ephem();
    let t0 = 60_000.0_f64;
    let elements = test_orbit(t0);
    let config = NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::PB(outfit::jpl_ephem::naif::naif_ids::planet_bary::PlanetaryBary::Jupiter),
        ],
        ..NBodyConfig::default()
    };

    // Table that comfortably covers the propagation span.
    let covering = PerturberEphemerisSet::build(&config, &jpl, t0 - 260.0, t0 + 260.0).unwrap();
    // Table that does NOT cover it: propagate_nbody must fall back to an
    // on-the-fly build rather than extrapolate.
    let too_small = PerturberEphemerisSet::build(&config, &jpl, t0 - 5.0, t0 + 5.0).unwrap();

    for span in [150.0_f64, -220.0] {
        let with_table = elements
            .propagate_nbody(t0 + span, &jpl, &config, Some(&covering))
            .unwrap();
        let without_table = elements
            .propagate_nbody(t0 + span, &jpl, &config, None)
            .unwrap();
        let with_small = elements
            .propagate_nbody(t0 + span, &jpl, &config, Some(&too_small))
            .unwrap();

        // Same physics whether the table is supplied or rebuilt internally.
        assert!(
            max_abs_diff(with_table.position, without_table.position) < 1e-10,
            "span {span}: supplied vs on-the-fly table disagree"
        );
        assert!(
            max_abs_diff(with_small.position, without_table.position) < 1e-10,
            "span {span}: non-covering table was not rebuilt"
        );
    }
}
