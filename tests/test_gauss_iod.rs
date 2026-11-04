//! Exact-oracle non-regression tests for the Gauss IOD pipeline. The reference
//! orbits below reproduce the in-house ephemeris reader's arithmetic bit-for-bit,
//! so the whole file is gated on that backend; the ANISE backend is exercised by
//! `test_ephemeris`, `test_diff_cor_nbody` and the parity tests in the library.
#![cfg(feature = "ephem-builtin")]

mod common;

use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::initial_orbit_determination::IODParams;
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::orbit_type::keplerian_element::KeplerianElements;
use outfit::orbit_type::OrbitalElements;
use outfit::FitIOD;
use outfit::{FullOrbitResult, JPLEphem};
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::common::approx_equal;

struct ExpectedResult {
    orbit: OrbitalElements,
    rms: f64,
}

fn expected_results() -> Vec<ExpectedResult> {
    vec![
        ExpectedResult {
            orbit: OrbitalElements::Keplerian {
                elements: KeplerianElements {
                    reference_epoch: 57049.2684537375,
                    semi_major_axis: 1.8017408357435853,
                    eccentricity: 0.2835625947849,
                    inclination: 0.20268281899795498,
                    ascending_node_longitude: 0.007951791820480579,
                    periapsis_argument: 1.2450647642589465,
                    mean_anomaly: 0.44080487866260754,
                },
                uncertainty: None,
                covariance: None,
            },
            rms: 66.97479294218962,
        },
        ExpectedResult {
            orbit: OrbitalElements::Keplerian {
                elements: KeplerianElements {
                    reference_epoch: 60672.2443617134,
                    semi_major_axis: 3.2199380906766213,
                    eccentricity: 0.062419209987893194,
                    inclination: 0.1829771029880662,
                    ascending_node_longitude: 0.030775930195269436,
                    periapsis_argument: 1.9053705720361855,
                    mean_anomaly: 4.980622835162594,
                },
                uncertainty: None,
                covariance: None,
            },
            rms: 0.5739558188544629,
        },
        ExpectedResult {
            orbit: OrbitalElements::Keplerian {
                elements: KeplerianElements {
                    reference_epoch: 60465.26777915681,
                    semi_major_axis: 2.1874983806597665,
                    eccentricity: 0.2025641449546787,
                    inclination: 0.11906245182869549,
                    ascending_node_longitude: 3.091806395999626,
                    periapsis_argument: 2.479324830682491,
                    mean_anomaly: 4.934465465979427,
                },
                uncertainty: None,
                covariance: None,
            },
            rms: 18.96375448637929,
        },
    ]
}

/// Load the DE440 in-house Horizon reader explicitly.
///
/// Unlike [`common::load_ephem`], this always resolves to the legacy DE
/// binary reader, regardless of whether `ephem-anise` is also enabled —
/// `JPLEphem::new` (which `common::load_ephem` follows) prefers ANISE
/// whenever it is compiled in, even alongside `ephem-builtin`. This file's
/// oracle values reproduce the Horizon reader's arithmetic bit-for-bit and
/// only agree with ANISE to `~1e-8`, so the tests need this specific reader
/// pinned explicitly rather than left to compile-time backend precedence.
fn load_horizon_ephem() -> JPLEphem {
    let source: EphemFileSource = "horizon:DE440"
        .try_into()
        .expect("failed to parse the Horizon ephemeris source");
    JPLEphem::from_builtin(source)
        .expect("failed to load the DE440 ephemeris through the in-house Horizon reader")
}

fn build_test_fixtures() -> (JPLEphem, Ut1Provider, ObsDataset, IODParams) {
    let ut1_provider = common::load_ut1();

    let jpl_ephem: JPLEphem = load_horizon_ephem();

    let (obs_dataset, errors) = ObsDataset::from_mpc_80_col_files(&[
        "tests/data/2015AB.obs",
        "tests/data/8467.obs",
        "tests/data/33803.obs",
    ]);
    if !errors.is_empty() {
        panic!("Failed to load observation datasets: {:?}", errors);
    }

    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()
        .unwrap();

    (jpl_ephem, ut1_provider, obs_dataset, iod_params)
}

fn assert_iod_results(mut full_orbit: FullOrbitResult, test_epsilon: f64, test_max_relative: f64) {
    // K09R05F
    let expected = &expected_results()[0];
    let best_orbit = full_orbit.remove(&"K09R05F".into()).unwrap().unwrap();
    assert!(approx_equal(
        &expected.orbit,
        best_orbit.orbital_elements(),
        test_epsilon
    ));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        expected.rms,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // 8467
    let expected = &expected_results()[1];
    let best_orbit = full_orbit.remove(&8467_u32.into()).unwrap().unwrap();
    assert!(approx_equal(
        &expected.orbit,
        best_orbit.orbital_elements(),
        test_epsilon
    ));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        expected.rms,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // 33803
    let expected = &expected_results()[2];
    let best_orbit = full_orbit.remove(&33803_u32.into()).unwrap().unwrap();
    assert!(approx_equal(
        &expected.orbit,
        best_orbit.orbital_elements(),
        test_epsilon
    ));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        expected.rms,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}

#[test]
fn test_gauss_iod() {
    let test_epsilon = 1e-11;
    let test_max_relative = 1e-11;

    let (jpl_ephem, ut1_provider, obs_dataset, iod_params) = build_test_fixtures();

    let full_orbit = obs_dataset
        .fit_full_iod(
            &jpl_ephem,
            &ut1_provider,
            &iod_params,
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    assert_iod_results(full_orbit, test_epsilon, test_max_relative);
}

#[test]
#[cfg(feature = "parallel")]
fn test_gauss_iod_parallel() {
    let test_epsilon = 1e-11;
    let test_max_relative = 1e-11;

    let (jpl_ephem, ut1_provider, obs_dataset, iod_params) = build_test_fixtures();

    let full_orbit = obs_dataset
        .fit_full_iod_parallel(
            &jpl_ephem,
            &ut1_provider,
            &iod_params,
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    assert_iod_results(full_orbit, test_epsilon, test_max_relative);
}
