use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::error_models::ErrorModel;
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;

#[test]
fn test_read_ades() {
    let mut outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let mut traj_set = TrajectorySet::new_from_ades(
        &mut outfit,
        &Utf8Path::new("tests/data/example_ades.xml"),
        None,
        None,
    );
    assert_eq!(traj_set.len(), 4);
    assert_eq!(traj_set.get(&ObjectNumber::Int(1234457)).unwrap().len(), 1);

    traj_set.add_from_ades(
        &mut outfit,
        &Utf8Path::new("tests/data/example_ades2.xml"),
        None,
        None,
    );

    assert_eq!(traj_set.len(), 7);
    let traj = traj_set
        .get(&ObjectNumber::String("2016 RD34".into()))
        .unwrap();
    assert_eq!(traj.len(), 2);
    let obs = traj.get(0).unwrap().get_observer(&outfit);
    assert_eq!(
        *obs.name.as_ref().unwrap(),
        "University of Hawaii 88-inch telescope, Maunakea".to_string()
    );

    traj_set.add_from_ades(
        &mut outfit,
        &Utf8Path::new("tests/data/flat_ades.xml"),
        None,
        None,
    );
    assert_eq!(traj_set.len(), 41);

    let traj = traj_set
        .get(&ObjectNumber::String("D/1993 F2-W".into()))
        .unwrap();

    assert_eq!(traj.len(), 1);
    assert_eq!(
        traj.get(0).unwrap().get_observer(&outfit).name,
        Some("La Palma".into())
    );
}
