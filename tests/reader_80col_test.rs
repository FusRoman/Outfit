use camino::Utf8Path;
use outfit::constants::TrajectorySet;
use outfit::observations::observations::TrajectoryExt;

#[test]
fn test_80col_reader() {
    let path_file = Utf8Path::new("tests/data/33803.obs");
    let mut traj_set = TrajectorySet::from_80col(&path_file);

    let obs_33803 = traj_set.get("33803").unwrap();
    assert_eq!(traj_set.len(), 1);
    assert_eq!(obs_33803.len(), 1812);
    assert_eq!(obs_33803[0].time, 43041.93878);
    assert_eq!(obs_33803[0].ra, 359.7403333333333);
    assert_eq!(obs_33803[0].dec, -0.5039444444444444);
    assert_eq!(obs_33803[0].observer, "049".to_string());

    let path_file = Utf8Path::new("tests/data/8467.obs");
    traj_set.add_80col(&path_file);
    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get("8467").unwrap();
    assert_eq!(obs_8467.len(), 3748);
    assert_eq!(obs_8467[0].time, 43785.35799);
    assert_eq!(obs_8467[0].ra, 14.62025);
    assert_eq!(obs_8467[0].dec, 9.987777777777778);
    assert_eq!(obs_8467[0].observer, "675".to_string());

    let path_file = Utf8Path::new("tests/data/K25D50B.obs");
    traj_set.add_80col(&path_file);
    assert_eq!(traj_set.len(), 3);
    let obs_k25 = traj_set.get("K25D50B").unwrap();
    assert_eq!(obs_k25.len(), 20);
    assert_eq!(obs_k25[0].time, 60732.280490000005);
    assert_eq!(obs_k25[0].ra, 154.65650833333333);
    assert_eq!(obs_k25[0].dec, 29.973188888888888);
    assert_eq!(obs_k25[0].observer, "V00".to_string());
    assert_eq!(obs_k25[19].observer, "691".to_string());
}
