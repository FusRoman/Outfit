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

    let traj_set = TrajectorySet::new_from_parquet(&path_file, ztf_observer);

    println!("{:?}", traj_set);
}
