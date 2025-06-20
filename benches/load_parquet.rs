use camino::Utf8Path;
use criterion::Throughput;
use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion};
use outfit::constants::TrajectorySet;
use outfit::error_models::ErrorModel;
use outfit::observations::trajectory_ext::TrajectoryExt;
use outfit::outfit::Outfit;

fn bench_load_parquet(c: &mut Criterion) {
    let mut outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
    let path = Utf8Path::new("tests/data/test_from_fink.parquet");
    let ztf_observer = outfit.get_observer_from_mpc_code(&"I41".into());

    let mut batch_group = c.benchmark_group("batch_sizes");

    for batch_size in [
        Some(2),
        Some(4),
        Some(8),
        Some(16),
        Some(32),
        Some(64),
        Some(128),
        Some(256),
        Some(512),
        Some(1024),
        Some(2048),
        Some(4096),
    ]
    .iter()
    {
        batch_group.throughput(Throughput::Elements(1));
        batch_group.bench_with_input(
            BenchmarkId::from_parameter(format!("{:?}", batch_size)),
            batch_size,
            |b, batch_size| {
                b.iter(|| {
                    let _ = TrajectorySet::new_from_parquet(
                        &mut outfit,
                        &path,
                        ztf_observer.clone(),
                        0.5,
                        0.5,
                        *batch_size,
                    );
                })
            },
        );
    }
}

criterion_group!(benches, bench_load_parquet);
criterion_main!(benches);
