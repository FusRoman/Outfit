use camino::Utf8Path;
use outfit::constants::TrajectorySet;
use outfit::observations::trajectory_ext::TrajectoryExt;


#[test]
fn test_load_traj_from_parquet() {
    let path_file = Utf8Path::new("tests/data/trajectories.parquet");

    let traj_set = TrajectorySet::new_from_parquet(&path_file);

}