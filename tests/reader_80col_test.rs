use camino::Utf8Path;
use outfit::constants::ObjectNumber;
use outfit::error_models::ErrorModel;
use outfit::outfit::Outfit;
use outfit::trajectories::trajectory_file::TrajectoryFile;
use outfit::TrajectorySet;

#[test]
fn test_80col_reader() {
    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let path_file = Utf8Path::new("tests/data/33803.obs");
    let mut traj_set = TrajectorySet::new_from_80col(&mut env_state, path_file);

    let obs_33803 = traj_set.get(&ObjectNumber::String("33803".into())).unwrap();
    assert_eq!(traj_set.len(), 1);
    assert_eq!(obs_33803.len(), 129);
    assert_eq!(obs_33803[0].time, 60324.52016874074);
    assert_eq!(obs_33803[0].ra, 3.5491391785131814);
    assert_eq!(obs_33803[0].dec, -0.15949710761897423);
    assert_eq!(
        obs_33803[0].get_observer(&env_state).name,
        Some("Mt. Lemmon Survey".to_string())
    );

    let path_file = Utf8Path::new("tests/data/8467.obs");
    traj_set.add_from_80col(&mut env_state, path_file);
    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get(&ObjectNumber::String("8467".into())).unwrap();
    assert_eq!(obs_8467.len(), 61);
    assert_eq!(obs_8467[0].time, 60647.053230740734);
    assert_eq!(obs_8467[0].ra, 0.10365423161131723);
    assert_eq!(obs_8467[0].dec, 0.1400047372376524);
    assert_eq!(
        obs_8467[0].get_observer(&env_state).name,
        Some("ATLAS Chile, Rio Hurtado".to_string())
    );

    let path_file = Utf8Path::new("tests/data/K25D50B.obs");
    traj_set.add_from_80col(&mut env_state, path_file);
    assert_eq!(traj_set.len(), 3);
    let obs_k25 = traj_set
        .get(&ObjectNumber::String("K25D50B".into()))
        .unwrap();
    assert_eq!(obs_k25.len(), 20);
    assert_eq!(obs_k25[0].time, 60732.28129074074);
    assert_eq!(obs_k25[0].ra, 2.6992652800547146);
    assert_eq!(obs_k25[0].dec, 0.5231308334332919);
    assert_eq!(
        obs_k25[0].get_observer(&env_state).name,
        Some("Kitt Peak-Bok".to_string())
    );
    assert_eq!(
        obs_k25[19].get_observer(&env_state).name,
        Some("Steward Observatory, Kitt Peak-Spacewatch".to_string())
    );
}
