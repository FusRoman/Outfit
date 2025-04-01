use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;

#[tokio::test(flavor = "multi_thread", worker_threads = 2)]
async fn test_80col_reader() {
    let mut env_state = Outfit::new();

    let path_file = Utf8Path::new("tests/data/33803.obs");
    let mut traj_set = TrajectorySet::new_from_80col(&mut env_state, &path_file);

    let obs_33803 = traj_set.get(&ObjectNumber::String("33803".into())).unwrap();
    assert_eq!(traj_set.len(), 1);
    assert_eq!(obs_33803.len(), 1812);
    assert_eq!(obs_33803[0].time, 43041.93878);
    assert_eq!(obs_33803[0].ra, 359.7403333333333);
    assert_eq!(obs_33803[0].dec, -0.5039444444444444);
    assert_eq!(
        obs_33803[0].get_observer(&env_state).name,
        Some("Uppsala-Kvistaberg".to_string())
    );

    let path_file = Utf8Path::new("tests/data/8467.obs");
    traj_set.add_from_80col(&mut env_state, &path_file);
    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get(&ObjectNumber::String("8467".into())).unwrap();
    assert_eq!(obs_8467.len(), 3748);
    assert_eq!(obs_8467[0].time, 43785.35799);
    assert_eq!(obs_8467[0].ra, 14.62025);
    assert_eq!(obs_8467[0].dec, 9.987777777777778);
    assert_eq!(
        obs_8467[0].get_observer(&env_state).name,
        Some("Palomar Mountain".to_string())
    );

    let path_file = Utf8Path::new("tests/data/K25D50B.obs");
    traj_set.add_from_80col(&mut env_state, &path_file);
    assert_eq!(traj_set.len(), 3);
    let obs_k25 = traj_set
        .get(&ObjectNumber::String("K25D50B".into()))
        .unwrap();
    assert_eq!(obs_k25.len(), 20);
    assert_eq!(obs_k25[0].time, 60732.280490000005);
    assert_eq!(obs_k25[0].ra, 154.65650833333333);
    assert_eq!(obs_k25[0].dec, 29.973188888888888);
    assert_eq!(
        obs_k25[0].get_observer(&env_state).name,
        Some("Kitt Peak-Bok".to_string())
    );
    assert_eq!(
        obs_k25[19].get_observer(&env_state).name,
        Some("Steward Observatory, Kitt Peak-Spacewatch".to_string())
    );
}
