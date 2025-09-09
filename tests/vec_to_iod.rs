use std::f64::consts::PI;

use approx::assert_relative_eq;
use outfit::observations::trajectory_ext::gauss_result_for;
use outfit::{
    observations::trajectory_ext::ObservationBatch, time::jd_to_mjd, ErrorModel, IODParams, Outfit,
    TrajectoryExt, TrajectorySet,
};
use outfit::{KeplerianElements, ObjectNumber};
use rand::rngs::StdRng;
use rand::SeedableRng;

#[inline]
fn angle_abs_diff(a: f64, b: f64) -> f64 {
    let tau = 2.0 * PI;
    let mut d = (a - b) % tau;
    if d > PI {
        d -= tau;
    }
    if d < -PI {
        d += tau;
    }
    d.abs()
}

pub fn assert_keplerian_approx_eq(
    got: &KeplerianElements,
    exp: &KeplerianElements,
    abs_eps: f64,
    rel_eps: f64,
) {
    // Scalars (non-angular)
    assert_relative_eq!(
        got.reference_epoch,
        exp.reference_epoch,
        epsilon = abs_eps,
        max_relative = rel_eps
    );
    assert_relative_eq!(
        got.semi_major_axis,
        exp.semi_major_axis,
        epsilon = abs_eps,
        max_relative = rel_eps
    );
    assert_relative_eq!(
        got.eccentricity,
        exp.eccentricity,
        epsilon = abs_eps,
        max_relative = rel_eps
    );

    // Angles (radians), compare with wrap-around
    for (name, g, e) in [
        ("inclination", got.inclination, exp.inclination),
        (
            "ascending_node_longitude",
            got.ascending_node_longitude,
            exp.ascending_node_longitude,
        ),
        (
            "periapsis_argument",
            got.periapsis_argument,
            exp.periapsis_argument,
        ),
        ("mean_anomaly", got.mean_anomaly, exp.mean_anomaly),
    ] {
        let diff = angle_abs_diff(g, e);
        // Allow absolute OR relative tolerance (whichever is larger).
        let tol = abs_eps.max(rel_eps * e.abs());
        assert!(
            diff <= tol,
            "Angle {name:?} differs too much: |Δ| = {diff:.6e} > tol {tol:.6e} (got={g:.15}, exp={e:.15})"
        );
    }
}

#[test]
fn vec_to_iod() {
    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(12)
        .max_triplets(30)
        .build()
        .unwrap();

    // ---------- Obj 33803 ----------
    let object_number = "0";
    let ra_deg = vec![
        20.9191548, 20.6388309, 20.6187259, 20.6137886, 19.892738, 18.2218784,
    ];
    let dec_deg = vec![
        20.0550441, 20.1218532, 20.1264229, 20.1275173, 20.2977473, 20.7096409,
    ];
    let time = jd_to_mjd(&[
        2458789.6362963,
        2458790.6995023,
        2458790.7733333,
        2458790.791412,
        2458793.6896759,
        2458801.7863426,
    ]);
    let observer = env_state.get_observer_from_mpc_code(&"049".to_string()); // Uppsala-Kvistaberg

    let batch = ObservationBatch::from_degrees_owned(&ra_deg, &dec_deg, 0.5, 0.5, &time);

    let mut traj_set: TrajectorySet =
        TrajectorySet::new_from_vec(&mut env_state, object_number, &batch, observer);

    let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility
    let orbits = traj_set.estimate_all_orbits(&env_state, &mut rng, &params);

    let traj_orbit = gauss_result_for(&orbits, &ObjectNumber::String("0".into()))
        .unwrap()
        .unwrap();

    let exp_orbit = KeplerianElements {
        reference_epoch: 58793.18358331323,
        semi_major_axis: 2.6906927279469097,
        eccentricity: 0.2657002116169482,
        inclination: 0.26700339023349917,
        ascending_node_longitude: 0.22380629928125767,
        periapsis_argument: 6.242813222404994,
        mean_anomaly: 0.2435936126817615,
    };

    let got_orbit = traj_orbit.0.as_inner().as_keplerian().unwrap();

    assert_keplerian_approx_eq(got_orbit, &exp_orbit, 1e-6, 1e-4);
}
