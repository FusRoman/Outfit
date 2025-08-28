#![cfg_attr(not(feature = "jpl-download"), allow(dead_code))]
#[cfg(not(feature = "jpl-download"))]
compile_error!("This benchmark requires the `jpl-download` feature: run with `cargo bench --features jpl-download`.");

use std::cell::RefCell;

use camino::Utf8Path;
use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
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

/// Run Gauss IOD on a single trajectory and return the best orbit + RMS.
/// This mirrors your test helper `run_iod`, with a fixed RNG seed for reproducibility.
fn run_iod(
    env_state: &mut Outfit,
    traj_set: &mut TrajectorySet,
    traj_number: &ObjectNumber,
) -> Result<(Option<GaussResult>, f64), OutfitError> {
    // Safety: unwrap is fine in bench context — data is known to be present.
    let obs = traj_set.get_mut(traj_number).expect("trajectory not found");
    let mut rng = StdRng::seed_from_u64(42);

    // Keep the same parameters as in your test to benchmark the exact pipeline you care about.
    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(obs.len())
        .max_triplets(30)
        .build()?;

    obs.estimate_best_orbit(env_state, &ErrorModel::FCCT14, &mut rng, &params)
}

/// Prepare the environment: Outfit with DE440 + three trajectories loaded from 80-col files.
/// This function performs file IO and may hit the on-disk ephemeris cache (first run only).
fn prepare_env() -> (Outfit, TrajectorySet, [ObjectNumber; 3]) {
    // Note: the first call may download DE440 if not cached; subsequent runs use the cache.
    let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).expect("Outfit init");

    let mut set = TrajectorySet::new_from_80col(&mut env, Utf8Path::new("tests/data/2015AB.obs"));
    set.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
    set.add_from_80col(&mut env, Utf8Path::new("tests/data/33803.obs"));

    let ids = [
        ObjectNumber::String("K09R05F".into()),
        ObjectNumber::String("8467".into()),
        ObjectNumber::String("33803".into()),
    ];

    (env, set, ids)
}

/// End-to-end benchmark: parses observation files and runs IOD on the three objects.
/// This captures the user-visible latency of a full IOD session (after ephemeris cache is warm).
fn bench_gauss_iod_e2e(c: &mut Criterion) {
    c.bench_function("gauss_iod_e2e_all", |b| {
        b.iter_batched(
            // Per-iteration setup (includes file IO).
            prepare_env,
            // Measured code path.
            |(mut env, mut set, ids)| {
                for id in &ids {
                    let res = run_iod(&mut env, &mut set, id).expect("IOD");
                    black_box(res);
                }
            },
            BatchSize::SmallInput,
        )
    });
}

/// Core benchmark: focuses on the Gauss IOD computation for a single object,
/// reusing a persistent `Outfit` and rebuilding the trajectory set per iteration.
/// If `TrajectorySet` implements `Clone`, you can replace the rebuild with a `clone()`.
fn bench_gauss_iod_core(c: &mut Criterion) {
    let (env0, _set_once, _ids) = prepare_env();
    let env = RefCell::new(env0);
    let target = ObjectNumber::String("K09R05F".into());
    let obs_path = Utf8Path::new("tests/data/2015AB.obs");

    c.bench_function("gauss_iod_core_single_refcell", |b| {
        b.iter_batched(
            // Setup: build a fresh TrajectorySet using a *temporary* mutable borrow of env.
            || {
                let mut env_borrow = env.borrow_mut();
                TrajectorySet::new_from_80col(&mut env_borrow, obs_path)
            },
            // Measure: borrow env again (the previous borrow is dropped at end of setup).
            |mut set| {
                let mut env_borrow = env.borrow_mut();
                let res = run_iod(&mut env_borrow, &mut set, &target).expect("IOD");
                black_box(res);
            },
            BatchSize::SmallInput,
        )
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default()
        // Tweak as needed; end-to-end runs are heavy, so keep sample size modest.
        .sample_size(30)
        .warm_up_time(std::time::Duration::from_secs(5))
        .measurement_time(std::time::Duration::from_secs(25))
        .with_plots();
    targets = bench_gauss_iod_e2e, bench_gauss_iod_core
);
criterion_main!(benches);
