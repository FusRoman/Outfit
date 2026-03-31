//! Benchmark Gauss Initial Orbit Determination (IOD) on survey-like
//! trajectories loaded from a Parquet file.
//!
//! # Overview
//!
//! This benchmark loads a Parquet dataset, builds a `TrajectorySet`, and
//! evaluates the performance of Gauss IOD across
//! increasing batch sizes of *distinct* trajectories (survey-like scenario).
//!
//! The goal is to estimate the average cost per orbit-fit when processing
//! large collections of independent tracks, similar to real survey pipelines.
//!
//! # Usage
//! ```bash
//! cargo bench --bench iod_batch_diff --features jpl-download
//! ```
//!
//! The `jpl-download` feature is required so that the DE440 ephemeris can be
//! fetched if not already cached locally.

#![cfg_attr(not(feature = "jpl-download"), allow(dead_code))]

use std::cell::RefCell;

use camino::Utf8Path;
use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use outfit::constants::ObjectNumber;
use outfit::error_models::ErrorModel;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::observations::observations_ext::ObservationIOD;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use outfit::trajectories::trajectory_file::TrajectoryFile;
use outfit::TrajectorySet;
use rand::rngs::StdRng;
use rand::SeedableRng;

/// Run Gauss Initial Orbit Determination on a single trajectory.
///
/// This function retrieves the observation list for a given trajectory,
/// constructs reasonable IOD parameters (noise realizations, triplet limits,
/// etc.) and calls [`ObservationIOD::estimate_best_orbit`].
///
/// # Arguments
/// * `env_state` – Mutable Outfit runtime environment (ephemerides, settings…)
/// * `traj_set` – The global trajectory set loaded from the Parquet file
/// * `traj_number` – The ID of the trajectory to solve
///
/// # Returns
/// A tuple `(GaussResult, rms)` if a valid orbit is found, or a wrapped
/// [`OutfitError`] otherwise.
///
/// Note that failures are **expected** for some trajectories in survey-like
/// datasets (insufficient geometry, poor conditioning, etc.).
fn run_iod(
    env_state: &mut Outfit,
    traj_set: &mut TrajectorySet,
    traj_number: &ObjectNumber,
) -> Result<(GaussResult, f64), OutfitError> {
    // Retrieve the mutable list of observations for this trajectory.
    let obs = traj_set.get_mut(traj_number).expect("trajectory not found");

    // Use a deterministic RNG so the benchmark remains reproducible.
    let mut rng = StdRng::seed_from_u64(42);

    // Build the parameters controlling the Gauss IOD solver.
    let params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(obs.len())
        .max_triplets(30)
        .build()?;

    // Perform the actual orbit estimation.
    obs.estimate_best_orbit(env_state, &ErrorModel::FCCT14, &mut rng, &params)
}

/// Prepare the Outfit environment and the full `TrajectorySet` from the
/// Fink-like Parquet test file.
///
/// This function:
/// 1. Loads DE440 via Outfit,
/// 2. Initializes the observer (ZTF Palomar I41),
/// 3. Reads all trajectories from the Parquet file,
/// 4. Returns the environment, the trajectory set, and all available IDs.
///
/// This preparation is done **once** and reused across benchmark iterations.
fn prepare_env_parquet() -> (Outfit, TrajectorySet, Vec<ObjectNumber>) {
    // Create the Outfit environment using DE440.
    let mut env_state =
        Outfit::new("horizon:DE440", ErrorModel::FCCT14).expect("Outfit init (DE440)");

    let path_file = std::env::var("OUTFIT_IOD_BATCH_DIFF_PARQUET")
        .unwrap_or_else(|_| "tests/data/test_from_fink.parquet".to_string());
    let path_file = Utf8Path::new(&path_file);

    // Using the Palomar/ZTF observer (MPC code I41).
    let ztf_observer = env_state.get_observer_from_mpc_code(&"I41".into());

    // Load all trajectories from the Parquet file.
    //
    // The last two parameters are the RA/Dec clipping radii used when building
    // tangent-plane residuals. They are intentionally loose in this benchmark.
    let traj_set =
        TrajectorySet::new_from_parquet(&mut env_state, path_file, ztf_observer, 0.5, 0.5, None)
            .expect("load TrajectorySet from Parquet");

    // Collect all trajectory IDs. Assumes the underlying representation exposes `.keys()`.
    let ids: Vec<ObjectNumber> = traj_set.keys().cloned().collect();

    (env_state, traj_set, ids)
}

/// Benchmark the scaling of Gauss IOD across increasingly large batches of
/// **distinct trajectories**.
///
/// Unlike micro-benchmarks which refit the same orbit repeatedly, this test
/// mimics a survey scenario: different objects, heterogeneous geometries,
/// and varying observation counts.
///
/// This provides a more realistic estimation of the throughput one can expect
/// in daily survey operations.
fn bench_gauss_iod_parquet_batch_scaling(c: &mut Criterion) {
    // Prepare environment and dataset once; reuse across all benchmark runs.
    let (env0, traj0, all_ids) = prepare_env_parquet();
    let env = RefCell::new(env0);
    let traj = RefCell::new(traj0);

    // Batch sizes (numbers of distinct trajectories per iteration).
    // These remain moderate so the benchmark does not take too long.
    let batch_sizes = [
        1_usize, 4, 16, 64, 256, 512, 1024, 2048, 4096, 8192, 16384, //, 32768, 65536, 131072,
    ];

    let mut group = c.benchmark_group("gauss_iod_parquet_batch_scaling");
    group.sample_size(10);

    for &batch in &batch_sizes {
        if batch > all_ids.len() {
            continue;
        }

        let name = format!("parquet_batch_{batch}");
        group.bench_function(&name, |b| {
            b.iter_batched(
                || (), // No heavy per-iteration setup needed
                |_| {
                    let mut env_borrow = env.borrow_mut();
                    let mut traj_borrow = traj.borrow_mut();

                    // Counters to ensure useful work is not optimized away.
                    let mut ok_count: usize = 0;
                    let mut err_count: usize = 0;

                    // Fit orbits for the first `batch` trajectory IDs.
                    for id in all_ids.iter().take(batch) {
                        match run_iod(&mut env_borrow, &mut traj_borrow, id) {
                            Ok(result) => {
                                black_box(result);
                                ok_count += 1;
                            }
                            Err(err) => {
                                // Failures may occur depending on observation geometry.
                                black_box(&err);
                                err_count += 1;
                            }
                        }
                    }

                    // Black-box the counters to prevent loop elimination.
                    black_box((ok_count, err_count));
                },
                BatchSize::SmallInput,
            )
        });
    }

    group.finish();
}

#[cfg(feature = "parallel")]
mod iod_bench_parallel {

    use super::*;
    use outfit::TrajectoryFit;

    /// Benchmark Gauss IOD using the parallel batch API over subsets of trajectories.
    ///
    /// For each `batch_size`, we build a `TrajectorySet` restricted to the first
    /// `batch_size` objects (according to the ID list returned by `prepare_env_parquet`),
    /// and run `estimate_all_orbits_in_batches_parallel` on that subset.
    ///
    /// This way, the x-axis ("batch size") represents the *number of trajectories*
    /// processed, just like in the sequential benchmark.
    pub fn bench_gauss_iod_parquet_parallel_batches(c: &mut Criterion) {
        // Prepare environment, full dataset, and the full list of IDs once.
        let (env0, traj0, all_ids) = prepare_env_parquet();
        let env = RefCell::new(env0);

        // Same batch sizes as the sequential benchmark so that comparisons are easy.
        let batch_sizes = [
            1_usize, 4, 16, 64, 256, 512, 1024, 2048, 4096, 8192,
            16384, // 32768, 65536, 131072,
        ];

        let mut group = c.benchmark_group("gauss_iod_parquet_parallel_batches");
        group.sample_size(10);

        for &batch_size in &batch_sizes {
            // Skip if we don't have enough trajectories.
            if batch_size > all_ids.len() {
                continue;
            }

            let name = format!("parquet_parallel_batch_{batch_size}");

            // ---- Build a subset TrajectorySet with exactly `batch_size` trajectories ----
            //
            // This is done once, *outside* Criterion's inner loop so that we are
            // measuring the IOD cost only (not the cost of building the subset).
            let selected_ids: Vec<ObjectNumber> =
                all_ids.iter().take(batch_size).cloned().collect();

            let selected_set: std::collections::HashSet<ObjectNumber> =
                selected_ids.iter().cloned().collect();

            // Clone the full set and retain only the selected trajectories.
            let mut traj_subset = traj0.clone();
            traj_subset.retain(|obj_id, _obs| selected_set.contains(obj_id));

            let traj = RefCell::new(traj_subset);

            // ---- Benchmark on the restricted subset ----
            group.bench_function(&name, |b| {
                b.iter_batched(
                    || {
                        // Per-iteration setup: deterministic RNG + IOD parameters.
                        let rng = StdRng::seed_from_u64(42);

                        let chunk_size = batch_size.min(64);
                        let params = IODParams::builder()
                            .n_noise_realizations(10)
                            .noise_scale(1.1)
                            .max_triplets(30)
                            .batch_size(chunk_size)
                            .build()
                            .expect("build IODParams");

                        (rng, params)
                    },
                    |(mut rng, params)| {
                        let env_borrow = env.borrow();
                        let mut traj_borrow = traj.borrow_mut();

                        let full_result = traj_borrow.estimate_all_orbits_in_batches_parallel(
                            &env_borrow,
                            &mut rng,
                            &params,
                        );

                        // Black-box the result to avoid dead-code elimination.
                        black_box(full_result);
                    },
                    BatchSize::SmallInput,
                )
            });
        }

        group.finish();
    }
}

#[cfg(feature = "parallel")]
criterion_group!(
    gauss_parquet_benches,
    bench_gauss_iod_parquet_batch_scaling,
    iod_bench_parallel::bench_gauss_iod_parquet_parallel_batches,
);

#[cfg(not(feature = "parallel"))]
criterion_group!(gauss_parquet_benches, bench_gauss_iod_parquet_batch_scaling,);

criterion_main!(gauss_parquet_benches);
