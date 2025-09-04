// benches/outfit_gauss_iod.rs

#[cfg(feature = "jpl-download")]
mod benches_impl {
    use std::cell::RefCell;

    use camino::Utf8Path;
    use criterion::{black_box, BatchSize, Criterion};
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
    fn run_iod(
        env_state: &mut Outfit,
        traj_set: &mut TrajectorySet,
        traj_number: &ObjectNumber,
    ) -> Result<(Option<GaussResult>, f64), OutfitError> {
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

    /// Prepare environment with DE440 and 3 trajectories.
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

    /// End-to-end benchmark.
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

    /// Core benchmark.
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
    targets = benches_impl::bench_gauss_iod_e2e, benches_impl::bench_gauss_iod_core
);

#[cfg(feature = "jpl-download")]
criterion_main!(benches);

// Fallback quand la feature est absente : fournit une main au crate.
#[cfg(not(feature = "jpl-download"))]
fn main() {
    eprintln!(
        "This benchmark requires the `jpl-download` feature. \
         Run with: `cargo bench --features jpl-download`"
    );
}
