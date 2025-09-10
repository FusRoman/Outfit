use outfit::constants::ObjectNumber;
use outfit::error_models::ErrorModel;
use outfit::observations::trajectory_ext::{ObservationBatch, TrajectoryExt};
use outfit::outfit::Outfit;
use outfit::TrajectorySet;

use std::sync::Arc;

#[test]
fn test_trajectories_from_vec() {
    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    // ---------- Obj 33803 ----------
    let traj_id = vec![33803];
    let ra_deg = vec![359.7403333333333_f64];
    let dec_deg = vec![-0.5039444444444444_f64];
    let time = vec![43041.93878_f64];
    let observer = env_state.get_observer_from_mpc_code(&"049".to_string()); // Uppsala-Kvistaberg

    let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &time);

    let mut traj_set: TrajectorySet =
        TrajectorySet::new_from_vec(&mut env_state, &batch, observer).unwrap();

    assert_eq!(traj_set.len(), 1);
    let obs_33803 = traj_set.get(&ObjectNumber::Int(33803)).unwrap();
    assert_eq!(obs_33803.len(), 1);

    assert!((obs_33803[0].time - 43041.93878).abs() < 1e-12);

    let expected_ra_rad = 359.7403333333333_f64.to_radians();
    let expected_dec_rad = (-0.5039444444444444_f64).to_radians();
    assert!((obs_33803[0].ra - expected_ra_rad).abs() < 1e-12);
    assert!((obs_33803[0].dec - expected_dec_rad).abs() < 1e-12);

    assert_eq!(
        obs_33803[0].get_observer(&env_state).name,
        Some("Uppsala-Kvistaberg".to_string())
    );

    // ---------- Obj 8467 ----------
    let traj_id = vec![8467];
    let ra_deg = vec![14.62025_f64];
    let dec_deg = vec![9.987777777777778_f64];
    let time = vec![43785.35799_f64];
    let observer: Arc<_> = env_state.get_observer_from_mpc_code(&"675".to_string()); // Palomar Mountain

    let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &time);

    traj_set
        .add_from_vec(&mut env_state, &batch, observer)
        .unwrap();

    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get(&ObjectNumber::Int(8467)).unwrap();
    assert_eq!(obs_8467.len(), 1);
    assert!((obs_8467[0].time - 43785.35799).abs() < 1e-12);

    let expected_ra_rad = 14.62025_f64.to_radians();
    let expected_dec_rad = 9.987777777777778_f64.to_radians();
    assert!((obs_8467[0].ra - expected_ra_rad).abs() < 1e-12);
    assert!((obs_8467[0].dec - expected_dec_rad).abs() < 1e-12);

    assert_eq!(
        obs_8467[0].get_observer(&env_state).name,
        Some("Palomar Mountain".to_string())
    );
}
