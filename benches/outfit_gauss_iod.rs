// benches/outfit_gauss_iod.rs

#[cfg(feature = "jpl-download")]
mod benches_impl {
    use std::cell::RefCell;

    use camino::Utf8Path;
    use criterion::{black_box, BatchSize, Criterion};
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
    /// This function:
    /// - retrieves the mutable observation list for the given trajectory,
    /// - builds a set of IOD parameters (noise realizations, triplet limits, …),
    /// - calls [`ObservationIOD::estimate_best_orbit`] and returns the result.
    ///
    /// Arguments
    /// -----------------
    /// * `env_state` - Mutable Outfit environment (DE ephemerides, settings, etc.).
    /// * `traj_set` - Global trajectory set containing the observations.
    /// * `traj_number` - ID of the trajectory to solve.
    ///
    /// Return
    /// ----------
    /// * A `Result` containing `(GaussResult, rms)` if a valid orbit is found,
    ///   or an [`OutfitError`] if the IOD fails.
    ///
    /// See also
    /// ------------
    /// * [`prepare_env`] - Helper to initialize the environment and trajectories.
    /// * [`bench_gauss_iod_e2e`] - End-to-end benchmark on three trajectories.
    /// * [`bench_gauss_iod_batch_scaling`] - Batch-scaling benchmark on the same set.
    fn run_iod(
        env_state: &mut Outfit,
        traj_set: &mut TrajectorySet,
        traj_number: &ObjectNumber,
    ) -> Result<(GaussResult, f64), OutfitError> {
        let obs = traj_set.get_mut(traj_number).expect("trajectory not found");
        let mut rng = StdRng::seed_from_u64(42);

        let params = IODParams::builder()
            .n_noise_realizations(10)
            .noise_scale(1.1)
            .max_obs_for_triplets(obs.len())
            .max_triplets(30)
            .build()?;

        obs.estimate_best_orbit(env_state, &ErrorModel::FCCT14, &mut rng, &params)
    }

    /// Prepare the Outfit environment (DE440) and a small set of test trajectories.
    ///
    /// This helper:
    /// - instantiates an [`Outfit`] environment using the DE440 ephemeris,
    /// - loads three 80-column observation files into a [`TrajectorySet`],
    /// - returns the environment, the trajectory set, and the three IDs.
    ///
    /// The three objects are:
    /// - `K09R05F` (corresponding to `2015AB.obs`),
    /// - `8467`,
    /// - `33803`.
    ///
    /// Arguments
    /// -----------------
    /// * None.
    ///
    /// Return
    /// ----------
    /// * A tuple `(env, set, ids)` where:
    ///   - `env` is the initialized [`Outfit`] environment,
    ///   - `set` is the [`TrajectorySet`] populated with the three objects,
    ///   - `ids` is the fixed array of [`ObjectNumber`] identifiers.
    ///
    /// See also
    /// ------------
    /// * [`run_iod`] - Core IOD routine using this prepared environment.
    /// * [`bench_gauss_iod_core`] - Reuses `env` but reloads trajectories each iteration.
    /// * [`bench_gauss_iod_batch_scaling`] - Uses these same objects in round-robin.
    fn prepare_env() -> (Outfit, TrajectorySet, [ObjectNumber; 3]) {
        let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).expect("Outfit init");

        let mut set =
            TrajectorySet::new_from_80col(&mut env, Utf8Path::new("tests/data/2015AB.obs"));
        set.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
        set.add_from_80col(&mut env, Utf8Path::new("tests/data/33803.obs"));

        let ids = [
            ObjectNumber::String("K09R05F".into()),
            ObjectNumber::String("8467".into()),
            ObjectNumber::String("33803".into()),
        ];

        (env, set, ids)
    }

    /// End-to-end benchmark on three real trajectories (full pipeline).
    ///
    /// This benchmark measures the cost of running Gauss IOD on all three
    /// trajectories, including environment and trajectory set initialization.
    ///
    /// Arguments
    /// -----------------
    /// * `c` - Criterion handle used to register the benchmark.
    ///
    /// See also
    /// ------------
    /// * [`bench_gauss_iod_core`] - Reuses the environment and reloads data per iteration.
    /// * [`bench_gauss_iod_batch_scaling`] - Scales the number of orbit fits per batch.
    pub fn bench_gauss_iod_e2e(c: &mut Criterion) {
        c.bench_function("gauss_iod_e2e_all", |b| {
            b.iter_batched(
                prepare_env,
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

    /// Core Gauss IOD benchmark: reuse environment, reload one trajectory per iteration.
    ///
    /// This benchmark keeps the [`Outfit`] environment in a shared [`RefCell`],
    /// and repeatedly:
    /// - rebuilds a `TrajectorySet` from `2015AB.obs`,
    /// - runs Gauss IOD only on `K09R05F`.
    ///
    /// This isolates the cost of:
    /// - loading a single trajectory from 80-column format,
    /// - performing one orbit determination for this trajectory.
    ///
    /// Arguments
    /// -----------------
    /// * `c` - Criterion handle used to register the benchmark.
    ///
    /// See also
    /// ------------
    /// * [`bench_gauss_iod_e2e`] - End-to-end benchmark over three trajectories.
    /// * [`bench_gauss_iod_batch_scaling`] - Varies the number of fits per batch.
    pub fn bench_gauss_iod_core(c: &mut Criterion) {
        let (env0, _set_once, _ids) = prepare_env();
        let env = RefCell::new(env0);
        let target = ObjectNumber::String("K09R05F".into());
        let obs_path = Utf8Path::new("tests/data/2015AB.obs");

        c.bench_function("gauss_iod_core_single_refcell", |b| {
            b.iter_batched(
                || {
                    let mut env_borrow = env.borrow_mut();
                    TrajectorySet::new_from_80col(&mut env_borrow, obs_path)
                },
                |mut set| {
                    let mut env_borrow = env.borrow_mut();
                    let res = run_iod(&mut env_borrow, &mut set, &target).expect("IOD");
                    black_box(res);
                },
                BatchSize::SmallInput,
            )
        });
    }

    /// Batch-scaling benchmark: N orbit fits per batch using the same 3 trajectories.
    ///
    /// This benchmark is useful to build "time per orbit-fit vs batch size"
    /// plots. It reuses the three real objects (K09R05F, 8467, 33803) and
    /// cycles through them in round-robin fashion for each batch.
    ///
    /// Arguments
    /// -----------------
    /// * `c` - Criterion handle used to register the benchmark.
    ///
    /// See also
    /// ------------
    /// * [`bench_gauss_iod_e2e`] - Baseline for full end-to-end performance.
    /// * [`bench_gauss_iod_core`] - Focuses on a single trajectory reload + fit.
    pub fn bench_gauss_iod_batch_scaling(c: &mut Criterion) {
        // Test several batch sizes: reuse the same three objects in round-robin.
        let batch_sizes = [1_usize, 2, 3, 6, 12, 24, 48, 96];

        for &batch in &batch_sizes {
            let name = format!("gauss_iod_batch_{batch}");
            c.bench_function(&name, |b| {
                b.iter_batched(
                    prepare_env,
                    |(mut env, mut set, ids)| {
                        for i in 0..batch {
                            let id = &ids[i % ids.len()];
                            let res = run_iod(&mut env, &mut set, id).expect("IOD");
                            black_box(res);
                        }
                    },
                    BatchSize::SmallInput,
                )
            });
        }
    }
}

#[cfg(feature = "jpl-download")]
use criterion::{criterion_group, criterion_main, Criterion};

#[cfg(feature = "jpl-download")]
criterion_group!(
    name = benches;
    config = Criterion::default()
        .sample_size(30)
        .warm_up_time(std::time::Duration::from_secs(5))
        .measurement_time(std::time::Duration::from_secs(25))
        .with_plots();
    targets =
        benches_impl::bench_gauss_iod_e2e,
        benches_impl::bench_gauss_iod_core,
        benches_impl::bench_gauss_iod_batch_scaling
);

#[cfg(feature = "jpl-download")]
criterion_main!(benches);

#[cfg(not(feature = "jpl-download"))]
fn main() {
    eprintln!(
        "This benchmark requires the `jpl-download` feature. \
         Run with: `cargo bench --features jpl-download`"
    );
}
