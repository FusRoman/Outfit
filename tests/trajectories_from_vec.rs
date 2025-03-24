use outfit::constants::TrajectorySet;
use outfit::observations::trajectory_ext::TrajectoryExt;

#[test]
fn test_trajectories_from_vec() {
    let object_number = "33803";
    let ra = vec![359.7403333333333];
    let dec = vec![-0.5039444444444444];
    let time = vec![43041.93878];
    let observer = "049";

    let mut traj_set: TrajectorySet =
        TrajectorySet::new_from_vec(object_number, &ra, &dec, &time, observer);

    assert_eq!(traj_set.len(), 1);
    let obs_33803 = traj_set.get("33803").unwrap();
    assert_eq!(obs_33803.len(), 1);
    assert_eq!(obs_33803[0].time, 43041.93878);
    assert_eq!(obs_33803[0].ra, 359.7403333333333);
    assert_eq!(obs_33803[0].dec, -0.5039444444444444);
    assert_eq!(obs_33803[0].observer, "049".to_string());

    let object_number = "8467";
    let ra = vec![14.62025];
    let dec = vec![9.987777777777778];
    let time = vec![43785.35799];
    let observer = "675";

    traj_set.add_from_vec(object_number, &ra, &dec, &time, observer);

    assert_eq!(traj_set.len(), 2);
    let obs_8467 = traj_set.get("8467").unwrap();
    assert_eq!(obs_8467.len(), 1);
    assert_eq!(obs_8467[0].time, 43785.35799);
    assert_eq!(obs_8467[0].ra, 14.62025);
    assert_eq!(obs_8467[0].dec, 9.987777777777778);
    assert_eq!(obs_8467[0].observer, "675".to_string());
}
