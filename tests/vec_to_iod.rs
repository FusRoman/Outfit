#![cfg(feature = "jpl-download")]

use std::collections::HashMap;
use std::f64::consts::PI;
use std::hash::DefaultHasher;

use approx::assert_relative_eq;
use outfit::time::utc_jd_slice_to_tt_mjd;
use outfit::trajectories::trajectory_fit::{gauss_result_for, FullOrbitResult};
use outfit::ObservationIOD;
use outfit::{
    trajectories::batch_reader::ObservationBatch, ErrorModel, IODParams, Outfit, TrajectoryFile,
    TrajectorySet,
};
use outfit::{KeplerianElements, ObjectNumber};
use rand::rngs::StdRng;
use rand::SeedableRng;
use std::hash::{Hash, Hasher};

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

fn assert_orbit(
    orbits: &FullOrbitResult,
    object_number: &ObjectNumber,
    exp_orbit: &KeplerianElements,
    exp_rms: f64,
    abs_eps: f64,
    rel_eps: f64,
) {
    let traj_orbit = gauss_result_for(orbits, object_number).unwrap().unwrap();

    let got_orbit = traj_orbit.0.as_inner().as_keplerian().unwrap();
    let got_rms = traj_orbit.1;

    assert_keplerian_approx_eq(got_orbit, exp_orbit, abs_eps, rel_eps);

    assert_relative_eq!(got_rms, exp_rms, epsilon = abs_eps, max_relative = rel_eps);
}

/// Derive a deterministic per-object seed from a global base seed and an [`ObjectNumber`].
///
/// This is used only in tests to ensure reproducible random number sequences
/// per object, independent of the `HashMap` iteration order.
///
/// Arguments
/// -----------------
/// * `base`: The global base seed (constant for the whole test run).
/// * `obj`: The object identifier used to derive a unique seed.
///
/// Return
/// ----------
/// * A `u64` seed value that can be fed into [`StdRng::seed_from_u64`].
fn seed_for_object(base: u64, obj: &ObjectNumber) -> u64 {
    let mut h = DefaultHasher::new();
    obj.hash(&mut h);
    let oh = h.finish();
    base ^ oh.rotate_left(17).wrapping_mul(0x9E37_79B9_7F4A_7C15)
}

#[test]
fn vec_to_iod() {
    // Create the Outfit environment with JPL DE440 ephemerides and
    // the FCCT14 astrometric error model.
    let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

    // Build the parameters for the Initial Orbit Determination (IOD).
    // - n_noise_realizations: number of noisy clones per triplet
    // - noise_scale: scale factor for astrometric uncertainties
    // - max_obs_for_triplets: upper bound on how many observations are used
    // - max_triplets: maximum number of triplets explored
    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.0)
        .max_obs_for_triplets(12)
        .max_triplets(30)
        .build()
        .unwrap();

    // Trajectory IDs: link each RA/DEC/time entry to one of the 3 objects (0,1,2).
    let traj_id = vec![0, 1, 2, 1, 2, 1, 0, 0, 0, 1, 2, 1, 1, 0, 2, 2, 0, 2, 2];

    // Right ascension values in degrees.
    let ra_deg = vec![
        20.9191548, 33.4247141, 32.1435128, 33.4159091, 32.1347282, 33.3829299, 20.6388309,
        20.6187259, 20.6137886, 32.7525147, 31.4874917, 32.4518231, 32.4495403, 19.892738,
        30.6416348, 30.0938936, 18.2218784, 28.3859403, 28.3818327,
    ];

    // Declination values in degrees.
    let dec_deg = vec![
        20.0550441, 23.5516817, 26.5139615, 23.5525348, 26.5160622, 23.5555991, 20.1218532,
        20.1264229, 20.1275173, 23.6064063, 26.6622284, 23.6270392, 23.6272157, 20.2977473,
        26.830301, 26.9256271, 20.7096409, 27.1602652, 27.160642,
    ];

    // Convert Julian Dates (UTC) to Modified Julian Dates (TT).
    let time = utc_jd_slice_to_tt_mjd(&[
        2458789.6362963,
        2458789.638125,
        2458789.638125,
        2458789.6663773,
        2458789.6663773,
        2458789.7706481,
        2458790.6995023,
        2458790.7733333,
        2458790.791412,
        2458791.8445602,
        2458791.8445602,
        2458792.8514699,
        2458792.8590741,
        2458793.6896759,
        2458794.7996759,
        2458796.7965162,
        2458801.7863426,
        2458803.7699537,
        2458803.7875231,
    ]);

    // Resolve the observer from MPC code I41 (ZTF - Palomar).
    let observer = env_state.get_observer_from_mpc_code(&"I41".to_string());

    // Build a batch of observations from RA/DEC/time arrays, with fixed uncertainties.
    // MJD needs to be in TT time scale for the IOD.
    let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &time);

    // Build the trajectory set (group observations by object).
    let mut traj_set: TrajectorySet =
        TrajectorySet::new_from_vec(&mut env_state, &batch, observer).unwrap();

    // Collect object identifiers and sort them for deterministic iteration order.
    let mut keys: Vec<_> = traj_set.keys().cloned().collect();
    keys.sort();

    // Base seed for reproducibility across runs.
    let base_seed = 42u64;

    // Container for final orbit results.
    let mut orbits: FullOrbitResult = HashMap::default();

    // For each object (in stable sorted order), derive a sub-RNG
    // from the base seed and object identifier.
    // This ensures reproducibility independent of HashMap iteration order.
    for obj in keys {
        let obs = traj_set.get_mut(&obj).expect("exists");
        let mut sub_rng = StdRng::seed_from_u64(seed_for_object(base_seed, &obj));

        // Perform Initial Orbit Determination (IOD) for this trajectory.
        let res =
            obs.estimate_best_orbit(&env_state, &env_state.error_model, &mut sub_rng, &params);

        orbits.insert(obj.clone(), res);
    }

    dbg!(&orbits);

    // ----- Validation for object 0 -----
    let exp_orbit = KeplerianElements {
        reference_epoch: 58793.18441761835,
        semi_major_axis: 2.6800907148611213,
        eccentricity: 0.26569343266377593,
        inclination: 0.26583378821788056,
        ascending_node_longitude: 0.22346768893017788,
        periapsis_argument: 6.2312396358686435,
        mean_anomaly: 0.25070862689966067,
    };
    let exp_rms = 0.41763854460398925;

    // Check that the computed orbit for object 0 matches the expected orbit.
    assert_orbit(
        &orbits,
        &ObjectNumber::Int(0),
        &exp_orbit,
        exp_rms,
        1e-6,
        1e-6,
    );

    // ----- Validation for object 1 -----
    let exp_orb_1 = KeplerianElements {
        reference_epoch: 58791.33825502848,
        semi_major_axis: 2.703826713146142,
        eccentricity: 0.30117809906476944,
        inclination: 0.25266176525888534,
        ascending_node_longitude: 0.317864582295561,
        periapsis_argument: 5.476805665094764,
        mean_anomaly: 0.6842587006557854,
    };
    let exp_rms_1 = 0.11855689104894801;

    assert_orbit(
        &orbits,
        &ObjectNumber::Int(1),
        &exp_orb_1,
        exp_rms_1,
        1e-6,
        1e-6,
    );

    // ----- Validation for object 2 -----
    let exp_orb_2 = KeplerianElements {
        reference_epoch: 58796.291754201615,
        semi_major_axis: 2.644301062854281,
        eccentricity: 0.29340523662320855,
        inclination: 0.23464394041365277,
        ascending_node_longitude: 0.21103621973001643,
        periapsis_argument: 6.263448848536585,
        mean_anomaly: 0.2960051245443382,
    };
    let exp_rms_2 = 0.3095556572660671;

    assert_orbit(
        &orbits,
        &ObjectNumber::Int(2),
        &exp_orb_2,
        exp_rms_2,
        1e-6,
        1e-6,
    );
}
