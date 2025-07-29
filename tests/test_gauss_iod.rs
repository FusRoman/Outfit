#![cfg(feature = "jpl-download")]

mod common;

use crate::common::assert_orbit_close;
use approx::assert_relative_eq;
use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::error_models::ErrorModel;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::keplerian_element::KeplerianElements;
use outfit::observations::observations_ext::ObservationIOD;
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use rand::rngs::StdRng;
use rand::SeedableRng;

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
    let test_epsilon = 1e-16;

    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let path_file = Utf8Path::new("tests/data/2015AB.obs");
    let mut traj_set = TrajectorySet::new_from_80col(&mut env_state, &path_file);

    let path_file = Utf8Path::new("tests/data/8467.obs");
    traj_set.add_from_80col(&mut env_state, &path_file);

    let path_file = Utf8Path::new("tests/data/33803.obs");
    traj_set.add_from_80col(&mut env_state, &path_file);

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("K09R05F".into()),
    )
    .unwrap();

    let expected_orbit = GaussResult::CorrectedOrbit(KeplerianElements {
        reference_epoch: 57049.22904452732,
        semi_major_axis: 1.8017634341924607,
        eccentricity: 0.28360400982137435,
        inclination: 0.20267485730439433,
        ascending_node_longitude: 0.008101820227108794,
        periapsis_argument: 1.2445523487100687,
        mean_anomaly: 0.44069989140090704,
    });

    let best_orbit_unwrapped = best_orbit.unwrap();
    assert_orbit_close(&best_orbit_unwrapped, &expected_orbit, test_epsilon);

    assert_relative_eq!(best_rms, 47.679542690830374, epsilon = test_epsilon);

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("8467".into()),
    )
    .unwrap();

    let expected_orbit = GaussResult::CorrectedOrbit(KeplerianElements {
        reference_epoch: 60672.24113100201,
        semi_major_axis: 3.1895469772492726,
        eccentricity: 0.05434034666133621,
        inclination: 0.18343383575588396,
        ascending_node_longitude: 0.032535949681617486,
        periapsis_argument: 2.0197545218041637,
        mean_anomaly: 4.850703837045108,
    });

    let best_orbit_unwrapped = best_orbit.unwrap();
    assert_orbit_close(&best_orbit_unwrapped, &expected_orbit, test_epsilon);
    assert_relative_eq!(best_rms, 0.550927559734149, epsilon = test_epsilon);

    let (best_orbit, best_rms) = run_iod(
        &mut env_state,
        &mut traj_set,
        &ObjectNumber::String("33803".into()),
    )
    .unwrap();

    let best_orbit_unwrapped = best_orbit.unwrap();

    let expected_orbit = GaussResult::CorrectedOrbit(KeplerianElements {
        reference_epoch: 60465.26778016307,
        semi_major_axis: 2.1921362022018465,
        eccentricity: 0.2042936374305541,
        inclination: 0.1189651192106595,
        ascending_node_longitude: 3.091130251223301,
        periapsis_argument: 2.4714439663663255,
        mean_anomaly: 4.9466622638824855,
    });

    assert_orbit_close(&best_orbit_unwrapped, &expected_orbit, test_epsilon);
    assert_relative_eq!(best_rms, 6.3193950851085035, epsilon = test_epsilon);
}
