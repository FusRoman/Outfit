use camino::Utf8Path;
use outfit::constants::TrajectorySet;
use outfit::observations::trajectory_ext::TrajectoryExt;

#[test]
fn test_load_traj_from_parquet() {
    let path_file = Utf8Path::new("tests/data/trajectories.parquet");

    let ztf_observer = "I41".to_string();

    let traj_set = TrajectorySet::new_from_parquet(&path_file, ztf_observer);
    assert_eq!(traj_set.len(), 4);
    assert_eq!(traj_set.get("1").unwrap().len(), 3);
}

#[test]
fn test_large_parquet() {
    let path_file = Utf8Path::new("tests/data/test_from_fink.parquet");

    let ztf_observer = "I41".to_string();

    let mut traj_set = TrajectorySet::new_from_parquet(&path_file, ztf_observer);

    assert_eq!(traj_set.len(), 2082);
    assert_eq!(traj_set.get("1").unwrap().len(), 6);
    assert_eq!(traj_set.get("2").unwrap().len(), 7);
    assert_eq!(traj_set.get("3").unwrap().len(), 7);
    assert_eq!(traj_set.get("4").unwrap().len(), 7);

    let path_file = Utf8Path::new("tests/data/trajectories.parquet");
    let rubin_observer = "X05".to_string();
    traj_set.add_from_parquet(path_file, rubin_observer);

    assert_eq!(traj_set.len(), 2082);
    assert_eq!(traj_set.get("1").unwrap().len(), 9);
    assert_eq!(traj_set.get("2").unwrap().len(), 10);
    assert_eq!(traj_set.get("3").unwrap().len(), 10);
    assert_eq!(traj_set.get("4").unwrap().len(), 10);
}
