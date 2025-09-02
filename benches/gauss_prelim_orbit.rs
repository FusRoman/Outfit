//! Benchmarks for GaussObs::prelim_orbit (single-threaded)
//!
//! Exemples d'exécution :
//!   cargo bench --bench gauss_prelim_orbit --features jpl-download
//!   cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/single_call
//!   cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/batch_100
//!   cargo bench gauss_prelim_orbit -- gauss_prelim_orbit/noisy_batch_100
//!
//! Astuce : vous pouvez aussi lancer en ligne de commande avec
//!   RAYON_NUM_THREADS=1 cargo bench --bench gauss_prelim_orbit

#![cfg_attr(not(feature = "jpl-download"), allow(dead_code))]

use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use nalgebra::{Matrix3, Vector3};

use outfit::initial_orbit_determination::gauss::GaussObs;
use outfit::initial_orbit_determination::gauss_result::GaussResult;
use outfit::initial_orbit_determination::IODParams;
use outfit::outfit::Outfit;
use outfit::outfit_errors::OutfitError;
use rand::SeedableRng;

/// Build global Outfit state (ephemerides, frames, etc.)
/// Note: keep this outside the hot loops.
fn build_state() -> Result<Outfit, OutfitError> {
    // English in-code comments per user preference:
    // Using FCCT14 as in production; adjust if you want to compare error models.
    use outfit::error_models::ErrorModel;
    Outfit::new("horizon:DE440", ErrorModel::FCCT14)
}

/// Deterministic GaussObs fixture (angles in radians, times in MJD TT)
fn make_fixture_gaussobs() -> GaussObs {
    let idx_obs = Vector3::new(0, 1, 2);
    let ra = Vector3::new(1.6894680985108945, 1.6898614520910629, 1.7526450904422723);
    let dec = Vector3::new(
        1.0825984522657437,
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
        .expect("GaussObs fixture")
}

/// Force Rayon (if used by downstream code) to a single worker thread.
/// Must be called before any Rayon pool is created by dependencies.
fn force_single_thread_pool() {
    // Pure env-var approach: works even if we don't depend on rayon directly.
    // Must be set very early, before any rayon usage.
    std::env::set_var("RAYON_NUM_THREADS", "1");

    // If you prefer hard enforcement and have rayon as a dev-dependency,
    // you could uncomment this block (but it's optional):
    //
    // #[cfg(any())]
    // {
    //     let _ = rayon::ThreadPoolBuilder::new()
    //         .num_threads(1)
    //         .build_global();
    //     // Ignore the error if the global pool was already built.
    // }
}

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

    // 2) Batch 100: reuse same input, measure algorithmic cost only
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

    // 3) Noisy batch 100: generate 100 noisy realizations then run prelim_orbit
    group.bench_function("noisy_batch_100", |b| {
        b.iter_batched(
            || gauss.clone(),
            |g| {
                // 0.3" ≈ 1.454e-6 rad; use your production noise scale if different.
                let sigma_rad = 1.5e-6_f64;
                let mut rng = rand::rngs::StdRng::seed_from_u64(42);
                let noisy = g.generate_noisy_realizations(
                    &(Vector3::zeros().add_scalar(1.0) * sigma_rad),
                    &(Vector3::zeros().add_scalar(1.0) * sigma_rad),
                    sigma_rad as usize,
                    1.0,
                    &mut rng,
                );

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

criterion_group!(gauss_benches, bench_prelim_orbit);
criterion_main!(gauss_benches);
