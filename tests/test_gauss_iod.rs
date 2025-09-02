#![cfg(feature = "jpl-download")]

mod common;

use approx::assert_relative_eq;
use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::error_models::ErrorModel;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::observations::observations_ext::ObservationIOD;
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::orbit_type::keplerian_element::KeplerianElements;
use outfit::orbit_type::OrbitalElements;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::common::approx_equal;

fn run_iod(
    env_state: &mut Outfit,
    traj_set: &mut TrajectorySet,
    traj_number: &ObjectNumber,
) -> Result<(Option<GaussResult>, f64), OutfitError> {
    let obs = traj_set.get_mut(traj_number).unwrap();
    let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

    let default = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(obs.len())
        .max_triplets(30)
        .build()?;

    obs.estimate_best_orbit(env_state, &ErrorModel::FCCT14, &mut rng, &default)
}

#[test]

fn test_gauss_iod() {
    let test_max_relative = 1e-13;
    let test_epsilon = 5. * f64::EPSILON;

    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let path_file = Utf8Path::new("tests/data/2015AB.obs");
    let mut traj_set = TrajectorySet::new_from_80col(&mut env_state, path_file);

    let path_file = Utf8Path::new("tests/data/8467.obs");
    traj_set.add_from_80col(&mut env_state, path_file);

    let path_file = Utf8Path::new("tests/data/33803.obs");
    traj_set.add_from_80col(&mut env_state, path_file);

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("K09R05F".into()),
    )
    .unwrap();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 57049.22904452732,
        semi_major_axis: 1.8017634341924542,
        eccentricity: 0.28360400982137396,
        inclination: 0.20267485730439427,
        ascending_node_longitude: 0.00810182022710516,
        periapsis_argument: 1.2445523487100616,
        mean_anomaly: 0.44069989140091426,
    });

    let best_orbit_unwrapped = best_orbit.unwrap();
    let orbit = best_orbit_unwrapped.get_orbit();

    dbg!("expected orbit: {:?}", &expected_orbit);
    dbg!("actual orbit: {:?}", orbit);

    dbg!("expected rms: {:?}", 47.67954270293223);
    dbg!("best rms: {:?}", best_rms);

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        47.67954270293223,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("8467".into()),
    )
    .unwrap();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60672.24113100201,
        semi_major_axis: 3.189546977249391,
        eccentricity: 0.05434034666134485,
        inclination: 0.18343383575588465,
        ascending_node_longitude: 0.03253594968161228,
        periapsis_argument: 2.0197545218038355,
        mean_anomaly: 4.85070383704545,
    });

    let best_orbit_unwrapped = best_orbit.unwrap();
    let orbit = best_orbit_unwrapped.get_orbit();

    dbg!("expected orbit: {:?}", &expected_orbit);
    dbg!("actual orbit: {:?}", orbit);

    dbg!("expected rms: {:?}", 0.5509275597328892);
    dbg!("best rms: {:?}", best_rms);

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        0.5509275597328892,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("33803".into()),
    )
    .unwrap();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60465.26778016307,
        semi_major_axis: 2.192136202201971,
        eccentricity: 0.2042936374305811,
        inclination: 0.1189651192106584,
        ascending_node_longitude: 3.091130251223283,
        periapsis_argument: 2.4714439663661487,
        mean_anomaly: 4.9466622638827324,
    });

    let best_orbit_unwrapped = best_orbit.unwrap();
    let orbit = best_orbit_unwrapped.get_orbit();

    dbg!("expected orbit: {:?}", &expected_orbit);
    dbg!("actual orbit: {:?}", orbit);

    dbg!("expected rms: {:?}", 6.319395087742966);
    dbg!("best rms: {:?}", best_rms);

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        6.319395085728921,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}
