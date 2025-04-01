use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn test_load_traj_from_parquet() {
    let mut env_state = Outfit::new();
    let path_file = Utf8Path::new("tests/data/trajectories.parquet");

    let ztf_observer = env_state.get_observer_from_mpc_code(&"I41".into());

    let traj_set = TrajectorySet::new_from_parquet(&mut env_state, &path_file, ztf_observer, None);
    assert_eq!(traj_set.len(), 4);
    assert_eq!(traj_set.get(&ObjectNumber::Int(1)).unwrap().len(), 3);
}

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn test_large_parquet() {
    let mut env_state = Outfit::new();
    let path_file = Utf8Path::new("tests/data/test_from_fink.parquet");

    let ztf_observer = env_state.get_observer_from_mpc_code(&"I41".into());

    let mut traj_set =
        TrajectorySet::new_from_parquet(&mut env_state, &path_file, ztf_observer, None);

    assert_eq!(traj_set.len(), 2082);
    assert_eq!(traj_set.get(&ObjectNumber::Int(1)).unwrap().len(), 6);
    assert_eq!(traj_set.get(&ObjectNumber::Int(2)).unwrap().len(), 7);
    assert_eq!(traj_set.get(&ObjectNumber::Int(3)).unwrap().len(), 7);
    assert_eq!(traj_set.get(&ObjectNumber::Int(4)).unwrap().len(), 7);

    let path_file = Utf8Path::new("tests/data/trajectories.parquet");
    let rubin_observer = env_state.get_observer_from_mpc_code(&"X05".into());
    traj_set.add_from_parquet(&mut env_state, path_file, rubin_observer, None);

    assert_eq!(traj_set.len(), 2082);
    let traj = traj_set.get(&ObjectNumber::Int(1)).unwrap();
    assert_eq!(traj.len(), 9);
    assert_eq!(traj_set.get(&ObjectNumber::Int(2)).unwrap().len(), 10);
    assert_eq!(traj_set.get(&ObjectNumber::Int(3)).unwrap().len(), 10);
    assert_eq!(traj_set.get(&ObjectNumber::Int(4)).unwrap().len(), 10);

    let first_obs = traj.get(0).unwrap();
    assert_eq!(first_obs.ra, 33.4247141);
    assert_eq!(first_obs.dec, 23.5516817);
    assert_eq!(first_obs.time, 58789.138125000056);
    assert_eq!(
        first_obs.get_observer(&env_state).name,
        Some("Palomar Mountain--ZTF".into())
    );

    let second_obs = traj.get(6).unwrap();
    assert_eq!(second_obs.ra, 192.1234);
    assert_eq!(second_obs.dec, -5.6789);
    assert_eq!(second_obs.time, 59396.0);
    assert_eq!(
        second_obs.get_observer(&env_state).name,
        Some("Simonyi Survey Telescope, Rubin Observatory".into())
    );
}
