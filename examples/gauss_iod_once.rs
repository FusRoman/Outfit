#![cfg(feature = "jpl-download")]
use std::env;

use camino::Utf8Path;
use outfit::constants::{ObjectNumber, TrajectorySet};
use outfit::error_models::ErrorModel;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::observations::observations_ext::ObservationIOD;
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use rand::rngs::StdRng;
use rand::SeedableRng;

use std::{thread, time::Duration};

/// Run Gauss IOD on a single trajectory and return the best orbit and RMS.
///
/// Arguments
/// -----------------
/// * `env_state`: The global environment (ephemeris, EOP, error model).
/// * `traj_set`: The trajectory container with parsed observations.
/// * `traj_number`: The object identifier present in `traj_set`.
///
/// Return
/// ----------
/// * `Ok((Option<GaussResult>, f64))` — the best orbit (if any) and its RMS in mas.
/// * `Err(OutfitError)` — if the IOD pipeline fails.
///
/// See also
/// ------------
/// * [`ObservationIOD::estimate_best_orbit`] – High-level Gauss IOD entry point.
/// * [`IODParams`] – Controls triplet selection and stochastic noise.
/// * [`ErrorModel::FCCT14`] – Default astrometric error model used here.
fn run_iod_once(
    env_state: &mut Outfit,
    traj_set: &mut TrajectorySet,
    traj_number: &ObjectNumber,
) -> Result<(Option<GaussResult>, f64), OutfitError> {
    let obs = traj_set
        .get_mut(traj_number)
        .expect("trajectory not found in set");
    let mut rng = StdRng::seed_from_u64(42);

    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(obs.len())
        .max_triplets(30)
        .build()?;

    obs.estimate_best_orbit(env_state, &ErrorModel::FCCT14, &mut rng, &params)
}

/// Minimal driver: load three test trajectories, run IOD for the requested object once.
/// Usage:
///   gauss_iod_once <OBJECT> [--verbose]
/// Example:
///   gauss_iod_once K09R05F --verbose
fn main() -> Result<(), OutfitError> {
    let mut args = env::args().skip(1).collect::<Vec<_>>();
    let verbose = if let Some(pos) = args.iter().position(|a| a == "--verbose") {
        args.remove(pos);
        true
    } else {
        false
    };

    let object = args
        .first()
        .cloned()
        .unwrap_or_else(|| "K09R05F".to_string());
    let obj = ObjectNumber::String(object);

    // Warm environment (will read DE440 from cache if already downloaded).
    let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;

    // Parse observations (adjust paths if needed).
    let mut set = TrajectorySet::new_from_80col(&mut env, Utf8Path::new("tests/data/2015AB.obs"));
    set.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
    set.add_from_80col(&mut env, Utf8Path::new("tests/data/33803.obs"));

    // Run IOD once for the requested object.
    let (best, rms) = run_iod_once(&mut env, &mut set, &obj)?;

    thread::sleep(Duration::from_millis(500)); // pause 0,5 s

    // Run IOD once for the requested object.
    let (best2, rms2) = run_iod_once(&mut env, &mut set, &obj)?;

    if verbose {
        eprintln!("[gauss_iod_once] object = {obj:?}, rms(mas) = {rms}");
        if let Some(gr) = best {
            eprintln!("[gauss_iod_once] orbit = {:?}", gr.get_orbit());
        } else {
            eprintln!("[gauss_iod_once] no orbit found");
        }

        eprintln!("[gauss_iod_once] object = {obj:?}, rms(mas) = {rms2}");
        if let Some(gr) = best2 {
            eprintln!("[gauss_iod_once] orbit = {:?}", gr.get_orbit());
        } else {
            eprintln!("[gauss_iod_once] no orbit found");
        }
    }

    Ok(())
}
