use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;

#[test]
fn test_trajectories_from_vec() {
    let mut env_state = Outfit::new("horizon:DE440");
    let object_number = "33803";
    let ra = vec![359.7403333333333];
    let dec = vec![-0.5039444444444444];
    let time = vec![43041.93878];
    let observer = env_state.get_observer_from_mpc_code(&"049".to_string());

    let mut traj_set: TrajectorySet = TrajectorySet::new_from_vec(
        &mut env_state,
        object_number,
        &ra,
        0.5,
        &dec,
        0.5,
        &time,
        observer,
    );

    assert_eq!(traj_set.len(), 1);
    let obs_33803 = traj_set.get(&ObjectNumber::String("33803".into())).unwrap();
    assert_eq!(obs_33803.len(), 1);
    assert_eq!(obs_33803[0].time, 43041.93878);
    assert_eq!(obs_33803[0].ra, 359.7403333333333);
    assert_eq!(obs_33803[0].dec, -0.5039444444444444);
    assert_eq!(
        obs_33803[0].get_observer(&env_state).name,
        Some("Uppsala-Kvistaberg".to_string())
    );

    let object_number = "8467";
    let ra = vec![14.62025];
    let dec = vec![9.987777777777778];
    let time = vec![43785.35799];
    let observer = env_state
        .get_observer_from_mpc_code(&"675".to_string())
        .clone();

    traj_set.add_from_vec(
        &mut env_state,
        object_number,
        &ra,
        0.5,
        &dec,
        0.5,
        &time,
        observer,
    );

    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get(&ObjectNumber::String("8467".into())).unwrap();
    assert_eq!(obs_8467.len(), 1);
    assert_eq!(obs_8467[0].time, 43785.35799);
    assert_eq!(obs_8467[0].ra, 14.62025);
    assert_eq!(obs_8467[0].dec, 9.987777777777778);
    assert_eq!(
        obs_8467[0].get_observer(&env_state).name,
        Some("Palomar Mountain".to_string())
    );
}
