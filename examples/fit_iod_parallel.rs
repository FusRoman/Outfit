//! Minimal walkthrough of **parallel** batch Initial Orbit Determination.
//!
//! `fit_full_iod` processes every trajectory in a dataset one after another.
//! When a batch holds many objects — the common case for survey-scale
//! processing — `fit_full_iod_parallel` (feature `parallel`, built on
//! [Rayon](https://docs.rs/rayon)) distributes trajectories across CPU cores
//! instead. Both produce **exactly the same result**: each trajectory's
//! random seed is derived deterministically from the caller's seed and the
//! trajectory's own identifier, independent of execution order or thread
//! count, so switching to the parallel version never changes what you get,
//! only how long it takes to get it.
//!
//! Run it with:
//! ```text
//! cargo run --release --example fit_iod_parallel --features parallel
//! ```
//! (`--release` matters here — this example times the fit; a debug build
//! would make both runs slow enough that the comparison is meaningless.)

use std::time::Instant;

use hifitime::ut1::Ut1Provider;
use outfit::{FitIOD, IODParams, JPLEphem, OutfitError};
use photom::io::polars::{ContiguousChoice, FromPolarsArgs};
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use polars::lazy::dsl::{col, lit};
use polars::lazy::frame::LazyFrame;
use rand::{rngs::StdRng, SeedableRng};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Load a batch of many trajectories.
    //
    // Parallelism only pays off when there is more than one independent
    // trajectory to spread across cores, so this example loads every
    // trajectory in the sample dataset that has at least 3 observations
    // (the minimum IOD needs), rather than a single object.
    let path_data = "tests/data/test_data_traj_str.parquet";
    let lf = LazyFrame::scan_parquet(path_data.into(), Default::default())
        .expect("scan_parquet must succeed")
        .filter(col("traj_id").is_not_null())
        .filter(
            col("traj_id")
                .count()
                .over([col("traj_id")])
                .gt_eq(lit(3u32)),
        );
    let polars_args = FromPolarsArgs {
        error_model: Some(ObsErrorModel::FCCT14),
        do_rechunk: Some(false),
        contiguous_choice: Some(ContiguousChoice::ContiguousTraj),
    };
    let obs_dataset = ObsDataset::from_lazy(lf, polars_args).unwrap();
    let n_traj = obs_dataset.iter_traj_id().unwrap().count();
    println!("Loaded {n_traj} trajectories\n");

    // Step 2 — Load the ephemeris, Earth-orientation data, and IOD tuning.
    let jpl_ephem: JPLEphem = "naif:DE440".try_into()?;
    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("failed to download UT1 Earth-orientation data");
    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()?;

    // Step 3 — Fit sequentially, timed.
    //
    // Both calls below seed their internal per-trajectory RNGs from the same
    // `StdRng::seed_from_u64(42)`, so their outputs are directly comparable.
    let start = Instant::now();
    let sequential = obs_dataset.clone().fit_full_iod(
        &jpl_ephem,
        &ut1_provider,
        &iod_params,
        ObsErrorModel::FCCT14,
        &mut StdRng::seed_from_u64(42),
    )?;
    let sequential_elapsed = start.elapsed();

    // Step 4 — Fit the same batch in parallel, timed.
    //
    // Only the function name changes; every argument is identical.
    let start = Instant::now();
    let parallel = obs_dataset.fit_full_iod_parallel(
        &jpl_ephem,
        &ut1_provider,
        &iod_params,
        ObsErrorModel::FCCT14,
        &mut StdRng::seed_from_u64(42),
    )?;
    let parallel_elapsed = start.elapsed();

    println!("Sequential: {sequential_elapsed:.2?} for {n_traj} trajectories");
    println!("Parallel:   {parallel_elapsed:.2?} for {n_traj} trajectories");
    println!(
        "Speedup:    {:.1}x\n",
        sequential_elapsed.as_secs_f64() / parallel_elapsed.as_secs_f64().max(1e-9)
    );

    // Step 5 — Confirm both runs agree exactly, trajectory by trajectory.
    //
    // This is the determinism guarantee in practice: same seed in, same
    // orbit out, regardless of how the work was scheduled across threads.
    let mut mismatches = 0usize;
    for (traj_id, seq_result) in sequential.iter() {
        let par_result = parallel.get(traj_id).expect("same trajectory set");
        match (seq_result, par_result) {
            (Ok(seq_orbit), Ok(par_orbit)) => {
                if seq_orbit.orbit_quality() != par_orbit.orbit_quality() {
                    mismatches += 1;
                }
            }
            (Err(_), Err(_)) => {}
            _ => mismatches += 1,
        }
    }
    println!("Trajectories where sequential and parallel disagree: {mismatches} / {n_traj}");

    Ok(())
}
