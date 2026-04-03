//! Benchmarks for `GaussObs::prelim_orbit` (single-threaded).
//!
//! # Overview
//!
//! This benchmark suite exercises the Gauss preliminary orbit determination
//! routine on a fixed synthetic fixture (`GaussObs`), in different usage
//! patterns:
//!
//! * Single call (micro-benchmark of the core algorithm).
//! * Fixed-size batches (e.g. 100 calls) to isolate algorithmic cost.
//! * Noisy realizations (Monte Carlo–like scenario).
//! * Scaling with batch size to study throughput and amortized overhead.
//!
//! All benchmarks are **single-threaded** by construction to avoid any
//! interference from Rayon or downstream parallelism.
//!
//! # Example commands
//!
//! ```bash
//! cargo bench --bench gauss_prelim_orbit --features jpl-download
//! cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/single_call
//! cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/batch_100
//! cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/noisy_single_call
//! cargo bench gauss_prelim_orbit -- gauss_prelim_orbit_batch_scaling/batch_size_256
//! ```
//!
//! Tip: you can also enforce single-threading from the command line with:
//!
//! ```bash
//! RAYON_NUM_THREADS=1 cargo bench --bench gauss_prelim_orbit
//! ```

#![cfg_attr(not(feature = "jpl-download"), allow(dead_code))]

use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use nalgebra::{Matrix3, Vector3};

use outfit::initial_orbit_determination::gauss::GaussObs;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use rand::SeedableRng;

/// Build the global Outfit state (ephemerides, reference frames, etc.).
///
/// This helper constructs an [`Outfit`] environment using the DE440 ephemeris
/// and the FCCT14 astrometric error model, matching the production setup.
///
/// Arguments
/// -----------------
/// * None.
///
/// Return
/// ----------
/// * A `Result` containing the initialized [`Outfit`] environment, or an
///   [`OutfitError`] if initialization fails.
///
/// See also
/// ------------
/// * [`make_fixture_gaussobs`] – Provides a deterministic Gauss fixture.
/// * [`bench_prelim_orbit`] – Uses this state for all prelim-orbit benchmarks.
/// * [`bench_prelim_orbit_vs_batch_size`] – Reuses the same state for
///   batch-scaling experiments.
fn build_state() -> Result<Outfit, OutfitError> {
    // Using FCCT14 as in production; adjust if you want to compare error models.
    use outfit::error_models::ErrorModel;
    Outfit::new("horizon:DE440", ErrorModel::FCCT14)
}

/// Create a deterministic `GaussObs` test fixture (angles in radians, times in MJD TT).
///
/// The fixture encodes:
/// * Three observation indices (0, 1, 2),
/// * Right ascensions and declinations (radians),
/// * Observation times (MJD TT),
/// * Heliocentric positions of the observer at each epoch (AU).
///
/// Arguments
/// -----------------
/// * None.
///
/// Return
/// ----------
/// * A [`GaussObs`] instance with fixed values suitable for reproducible
///   benchmarking.
///
/// See also
/// ------------
/// * [`build_state`] – Returns the global Outfit environment.
/// * [`bench_prelim_orbit`] – Benchmarks `prelim_orbit` on this fixture.
/// * [`bench_prelim_orbit_vs_batch_size`] – Reuses the same fixture to explore
///   performance scaling with batch size.
fn make_fixture_gaussobs() -> GaussObs {
    let idx_obs = Vector3::new(0, 1, 2);
    let ra = Vector3::new(
        1.689_468_098_510_894_5,
        1.689_861_452_091_062_9,
        1.752_645_090_442_272_3,
    );
    let dec = Vector3::new(
        1.082_598_452_265_743_7,
        0.943_679_018_934_623_1,
        0.827_517_321_571_201_4,
    );
    let time = Vector3::new(
        57_028.454_047_592_59,
        57_049.231_857_592_59,
        57_063.959_487_592_59,
    );

    let observer_helio_position = Matrix3::new(
        -0.264_135_633_607_079,
        -0.588_973_552_650_573_5,
        -0.774_192_148_350_372,
        0.869_046_620_910_086,
        0.724_011_718_791_646,
        0.561_510_219_548_918_2,
        0.376_746_685_666_572_5,
        0.313_873_420_677_094,
        0.243_444_791_401_658_5,
    );

    GaussObs::with_observer_position(idx_obs, ra, dec, time, observer_helio_position)
}

/// Force Rayon (if used by downstream code) to a single worker thread.
///
/// This function sets the `RAYON_NUM_THREADS` environment variable to `"1"`
/// **before** any Rayon thread pool is created, ensuring that all parallel
/// code in dependencies is effectively single-threaded.
///
/// Arguments
/// -----------------
/// * None.
///
/// Return
/// ----------
/// * Nothing. The function is used for its side-effect on the process
///   environment.
///
/// See also
/// ------------
/// * [`bench_prelim_orbit`] – Calls this before any benchmark runs.
/// * [`bench_prelim_orbit_vs_batch_size`] – Also enforces single-threaded
///   execution for scaling benchmarks.
fn force_single_thread_pool() {
    // Pure env-var approach: works even if we do not depend on Rayon directly.
    // Must be set very early, before any Rayon usage.
    std::env::set_var("RAYON_NUM_THREADS", "1");

    // If you prefer hard enforcement and have Rayon as a dev-dependency,
    // you could uncomment this block (but it is optional):
    //
    // #[cfg(any())]
    // {
    //     let _ = rayon::ThreadPoolBuilder::new()
    //         .num_threads(1)
    //         .build_global();
    //     // Ignore the error if the global pool was already built.
    // }
}

/// Benchmark various usage patterns of `GaussObs::prelim_orbit` on a fixed fixture.
///
/// This benchmark group (`gauss_prelim_orbit`) contains:
///
/// 1. **single_call** – A single call to `prelim_orbit`, micro-benchmarking
///    the solver itself.
/// 2. **batch_100** – Reuse the same input and call `prelim_orbit` 100 times
///    per iteration, focusing on algorithmic cost.
/// 3. **noisy_single_call** – Generate noisy realizations of the same geometry
///    (Monte Carlo–style) and call `prelim_orbit` on each noisy instance.
///
/// Single-threaded execution is enforced to avoid interference from any
/// Rayon-based parallelism in downstream code.
///
/// Arguments
/// -----------------
/// * `c` - Criterion handle used to register the benchmarks.
///
/// Return
/// ----------
/// * Nothing. The function registers benchmark cases in the given
///   [`Criterion`] context.
///
/// See also
/// ------------
/// * [`bench_prelim_orbit_vs_batch_size`] – Explores scaling vs. batch size.
/// * [`make_fixture_gaussobs`] – Provides the deterministic Gauss fixture.
/// * [`build_state`] – Constructs the shared Outfit environment.
fn bench_prelim_orbit(c: &mut Criterion) {
    // Ensure single-threaded execution before anything else touches Rayon.
    force_single_thread_pool();

    let mut group = c.benchmark_group("gauss_prelim_orbit");

    // Build global state once (outside hot loops).
    let state = build_state().expect("Outfit state");

    // Deterministic fixture.
    let gauss = make_fixture_gaussobs();

    // 1) Single call
    group.bench_function("single_call", |b| {
        b.iter(|| {
            let res = gauss.prelim_orbit(black_box(&state), &IODParams::default());
            match res {
                Ok(GaussResult::PrelimOrbit(_)) | Ok(GaussResult::CorrectedOrbit(_)) => {}
                Err(e) => panic!("prelim_orbit failed: {e:?}"),
            }
        })
    });

    // 2) Batch 100: reuse same input, measure algorithmic cost only.
    group.bench_function("batch_100", |b| {
        b.iter_batched(
            || gauss.clone(),
            |g| {
                for _ in 0..100 {
                    let res = g.prelim_orbit(&state, &IODParams::default());
                    black_box(&res);
                }
            },
            BatchSize::SmallInput,
        )
    });

    // 3) Noisy: generate noisy realizations then run prelim_orbit on each.
    group.bench_function("noisy_single_call", |b| {
        b.iter_batched(
            || gauss.clone(),
            |g| {
                // 0.3" ≈ 1.454e-6 rad; we use ~1.5e-6 rad as a representative astrometric noise.
                let sigma_rad = 1.5e-6_f64;
                let mut rng = rand::rngs::StdRng::seed_from_u64(42);

                // Per-angle standard deviations (same for RA and Dec here).
                let sigma_vec = Vector3::zeros().add_scalar(sigma_rad);

                let noisy =
                    g.generate_noisy_realizations(&sigma_vec, &sigma_vec, 100, 1.0, &mut rng);

                for gg in noisy {
                    let res = gg.prelim_orbit(&state, &IODParams::default());
                    black_box(&res);
                }
            },
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

/// Benchmark the scaling of `prelim_orbit` with the number of calls per batch.
///
/// This group (`gauss_prelim_orbit_batch_scaling`) explores how the cost per
/// orbit-fit behaves as you increase the batch size:
///
/// * For small batches, overheads (Criterion, function calls, etc.) dominate.
/// * For large batches, the raw cost of the preliminary orbit algorithm
///   becomes more apparent, yielding a stable "time per fit".
///
/// This is particularly useful to derive a "time per orbit-fit" curve for
/// performance notes or publications.
///
/// Arguments
/// -----------------
/// * `c` - Criterion handle used to register the benchmarks.
///
/// Return
/// ----------
/// * Nothing. The function registers the batch-scaling benchmarks.
///
/// See also
/// ------------
/// * [`bench_prelim_orbit`] – Baseline single-call and fixed-batch benchmarks.
/// * [`make_fixture_gaussobs`] – Fixed geometry used for all fits.
/// * [`build_state`] – Shared Outfit environment reused across batches.
fn bench_prelim_orbit_vs_batch_size(c: &mut Criterion) {
    // Ensure single-threaded execution before anything else touches Rayon.
    force_single_thread_pool();

    let mut group = c.benchmark_group("gauss_prelim_orbit_batch_scaling");

    // Build global state once (outside hot loops).
    let state = build_state().expect("Outfit state");

    // Deterministic fixture.
    let gauss = make_fixture_gaussobs();

    // Different batch sizes to explore scaling behaviour.
    let batch_sizes = [1_usize, 4, 16, 64, 128, 256, 512];

    for &batch in &batch_sizes {
        let name = format!("batch_size_{batch}");
        group.bench_function(&name, |b| {
            b.iter_batched(
                || gauss.clone(),
                |g| {
                    for _ in 0..batch {
                        let res = g.prelim_orbit(&state, &IODParams::default());
                        black_box(&res);
                    }
                },
                BatchSize::SmallInput,
            )
        });
    }

    group.finish();
}

criterion_group!(
    gauss_benches,
    bench_prelim_orbit,
    bench_prelim_orbit_vs_batch_size
);
criterion_main!(gauss_benches);
