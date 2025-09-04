use criterion::{black_box, criterion_group, criterion_main, BatchSize, Criterion};
use rand::rngs::StdRng;
use rand::{Rng, SeedableRng};

use outfit::orbit_type::equinoctial_element::EquinoctialElements;

/// Uniform random in [0, 2π)
#[inline]
fn rand_angle(rng: &mut StdRng) -> f64 {
    let two_pi = std::f64::consts::TAU;
    rng.random::<f64>() * two_pi
}

/// Build equinoctial elements with only h,k,λ set (others neutral).
#[inline]
fn make_equinoctial(h: f64, k: f64, lambda: f64) -> EquinoctialElements {
    EquinoctialElements {
        reference_epoch: 59000.0,
        semi_major_axis: 1.0, // not used here
        eccentricity_sin_lon: h,
        eccentricity_cos_lon: k,
        tan_half_incl_sin_node: 0.0,
        tan_half_incl_cos_node: 0.0,
        mean_longitude: lambda,
    }
}

/// Adjust λ so that λ >= ϖ, mimicking the production path in solve_two_body_problem.
#[inline]
fn align_lambda_with_periapsis(lambda: f64, w: f64) -> f64 {
    let mut lam = lambda.rem_euclid(std::f64::consts::TAU);
    let w_mod = w.rem_euclid(std::f64::consts::TAU);
    if lam < w_mod {
        lam += std::f64::consts::TAU;
    }
    lam
}

/// Typical regime: e ∈ [0.0, 0.7]
fn bench_typical(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0xDEADBEEF);
    let samples = 10_000usize;

    c.bench_function("solve_kepler_equation/typical_e<=0.7", |b| {
        b.iter_batched(
            || {
                // Pre-generate inputs to avoid RNG cost in the timed section
                (0..samples)
                    .map(|_| {
                        let e = rng.random_range(0.0..=0.7);
                        let w = rand_angle(&mut rng);
                        let lambda = align_lambda_with_periapsis(rand_angle(&mut rng), w);
                        let h = e * w.sin();
                        let k = e * w.cos();
                        (h, k, lambda, w)
                    })
                    .collect::<Vec<_>>()
            },
            |cases| {
                // Benchmark only the solver calls
                for (h, k, lambda, w) in cases {
                    let equ = make_equinoctial(h, k, lambda);
                    let f = equ
                        .solve_kepler_equation(black_box(lambda), black_box(w))
                        .unwrap();
                    black_box(f);
                }
            },
            BatchSize::LargeInput,
        )
    });
}

/// High-eccentricity (still elliptic): e ∈ [0.7, 0.9]
fn bench_high_e(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0xBADF00D);
    let samples = 10_000usize;

    c.bench_function("solve_kepler_equation/high_e_0.7..0.9", |b| {
        b.iter_batched(
            || {
                (0..samples)
                    .map(|_| {
                        let e = rng.random_range(0.7..0.9);
                        let w = rand_angle(&mut rng);
                        let lambda = align_lambda_with_periapsis(rand_angle(&mut rng), w);
                        let h = e * w.sin();
                        let k = e * w.cos();
                        (h, k, lambda, w)
                    })
                    .collect::<Vec<_>>()
            },
            |cases| {
                for (h, k, lambda, w) in cases {
                    let equ = make_equinoctial(h, k, lambda);
                    let _ = equ.solve_kepler_equation(black_box(lambda), black_box(w));
                }
            },
            BatchSize::LargeInput,
        )
    });
}

/// Near-circular regime: e ≈ 1e-12
fn bench_near_circular(c: &mut Criterion) {
    let mut rng = StdRng::seed_from_u64(0xFEEDFACE);
    let samples = 10_000usize;
    let e = 1e-12;

    c.bench_function("solve_kepler_equation/near_circular_e=1e-12", |b| {
        b.iter_batched(
            || {
                (0..samples)
                    .map(|_| {
                        let w = rand_angle(&mut rng);
                        let lambda = align_lambda_with_periapsis(rand_angle(&mut rng), w);
                        let h = e * w.sin();
                        let k = e * w.cos();
                        (h, k, lambda, w)
                    })
                    .collect::<Vec<_>>()
            },
            |cases| {
                for (h, k, lambda, w) in cases {
                    let equ = make_equinoctial(h, k, lambda);
                    let f = equ
                        .solve_kepler_equation(black_box(lambda), black_box(w))
                        .unwrap();
                    black_box(f);
                }
            },
            BatchSize::LargeInput,
        )
    });
}

/// Fixed “stress” case (from a previous failing input), useful for stability profiling.
fn bench_fixed_stress(c: &mut Criterion) {
    // Example numbers; feel free to replace with your own worst-case input.
    let e = 0.75_f64;
    let w = -2.823_013_355_485_587_6_f64;
    let lambda = 5.930_860_541_086_263_f64;
    let h = e * w.sin();
    let k = e * w.cos();
    let lambda = align_lambda_with_periapsis(lambda, w);
    let equ = make_equinoctial(h, k, lambda);

    c.bench_function("solve_kepler_equation/fixed_stress_case", |b| {
        b.iter(|| {
            let f = equ.solve_kepler_equation(black_box(lambda), black_box(w));
            black_box(f.ok());
        })
    });
}

criterion_group!(
    name = benches;
    config = Criterion::default();
    targets = bench_typical, bench_high_e, bench_near_circular, bench_fixed_stress
);
criterion_main!(benches);
