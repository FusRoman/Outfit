//! Criterion benchmarks for [`propagate_universal`], the two-body
//! universal-variable propagator used by fink-fat's Kalman filter for every
//! per-object state update (tens of thousands of calls in production).
//!
//! Fixtures reuse the exact `(position, velocity, dt)` triples already
//! validated in `src/kepler/propagation.rs::tests_propagation_edge_cases`
//! and `tests_propagation` (real fink-fat state) — each fixture below names
//! its corresponding test for traceability.
//!
//! Groups
//! ------
//! * `propagate_universal/scenarios` — end-to-end cost per orbital regime /
//!   gap length actually seen in production.
//! * `propagate_universal/solver_kind` — `NewtonRaphson` vs `BrentDecker` vs
//!   `Auto`, to quantify the cost of the automatic fallback.
//! * `propagate_universal/warm_start_chain` — a 20-step daily-cadence chain
//!   (as a Kalman filter would run it), cold-started every step vs
//!   warm-started from the previous `psy`.
//! * `kepler/components` — isolated cost of `prelim_kepuni` and `s_funct`,
//!   the two functions called on every Newton/Brent-Dekker iteration.
//!
//! Run with `cargo bench --bench propagate_universal`. To profile with
//! `perf` (matching the `tracking_analysis` workflow), build once with
//! `cargo bench --bench propagate_universal --no-run`, then `perf record`
//! directly against the resulting `target/release/deps/propagate_universal-*`
//! executable, optionally filtering to one Criterion group by passing its
//! name as an extra argument (Criterion's own CLI filter).

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use nalgebra::Vector3;
use outfit::kepler::{
    prelim_elliptic, prelim_hyperbolic, s_funct, ParabolicPrelimMethod, SolverKind, SolverParams,
    SolverType, UniversalKeplerParams,
};
use outfit::{kepler::propagate_universal, GAUSS_GRAV};

/// Reference MJD epoch (arbitrary; only `dt` matters for two-body dynamics).
const T0: f64 = 60000.0;

fn default_solver_type() -> SolverType {
    SolverType {
        kind: SolverKind::Auto,
        params: SolverParams {
            convergency: 2.220446049250313e-14,
            psi_guess: None,
            max_iter_prelim_kepuni: 20,
            parabolic_solving_method: ParabolicPrelimMethod::Cardano,
        },
    }
}

/// One `propagate_universal` input scenario, named after the corresponding
/// test in `src/kepler/propagation.rs`.
struct Scenario {
    name: &'static str,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    dt: f64,
}

fn scenarios() -> Vec<Scenario> {
    vec![
        // tests_propagation::test_propag — real fink-fat Kalman-filter state.
        Scenario {
            name: "real_fink_fat_state",
            position: Vector3::new(
                -8.264_959_160_036_185e-1,
                3.919_660_608_486_096_3e-1,
                2.229_919_607_182_842_5e-2,
            ),
            velocity: Vector3::new(
                -5.447_367_111_934_2e-3,
                -2.107_596_146_728_544e-2,
                1.560_811_152_125_889_6e-3,
            ),
            dt: 1.992476899875328,
        },
        // tests_propagation_edge_cases::test_quasi_circular_orbit — a=2.3 AU, e=1e-4.
        Scenario {
            name: "quasi_circular",
            position: Vector3::new(0.21053202031103074, 2.2808559935794808, 0.20574336652469302),
            velocity: Vector3::new(
                -0.011240896923812553,
                0.0009286862636099264,
                0.0012093812594121695,
            ),
            dt: 1.0,
        },
        // tests_propagation_edge_cases::test_high_eccentricity_near_perihelion — a=2.5 AU, e=0.95.
        Scenario {
            name: "high_eccentricity_near_perihelion",
            position: Vector3::new(
                -0.07334498751332834,
                -0.33895907499571465,
                -0.028392141622553987,
            ),
            velocity: Vector3::new(
                0.016602035346734694,
                -0.03512248571568977,
                -0.00855803655134319,
            ),
            dt: 2.0,
        },
        // tests_propagation_edge_cases::test_near_parabolic_elliptic — alpha = -1e-6.
        Scenario {
            name: "near_parabolic_elliptic",
            position: Vector3::new(-39.52202041820489, -31.796429108199263, -8.899233495929744),
            velocity: Vector3::new(
                -0.0021940078334573457,
                -0.0024735491628538825,
                -0.0007479533467381286,
            ),
            dt: 5.0,
        },
        // tests_propagation_edge_cases::test_near_parabolic_hyperbolic — alpha = +1e-6.
        Scenario {
            name: "near_parabolic_hyperbolic",
            position: Vector3::new(-39.52295104602465, -31.79684024886987, -8.899322047228939),
            velocity: Vector3::new(
                0.002194072590815172,
                0.0024735663723493223,
                0.0007479554224776113,
            ),
            dt: 5.0,
        },
        // tests_propagation_edge_cases::test_hyperbolic_orbit — a=-1.5 AU, e=2.0.
        Scenario {
            name: "hyperbolic",
            position: Vector3::new(1.249244222099196, 0.79545767943562, 0.31874549259903184),
            velocity: Vector3::new(
                0.011307076273210306,
                -0.019273992715882825,
                -0.00941294513552841,
            ),
            dt: 10.0,
        },
        // tests_propagation_edge_cases::test_gap_35_days_ztf_lsst_cadence.
        Scenario {
            name: "gap_35_days",
            position: Vector3::new(-0.5081279057987266, 1.9247110895634725, 0.35100377520249654),
            velocity: Vector3::new(
                -0.01246619702904591,
                -0.0016710189605835082,
                0.00028431118257845765,
            ),
            dt: 35.0,
        },
        // tests_propagation_edge_cases::test_gap_400_days_multi_revolution.
        Scenario {
            name: "gap_400_days_multi_revolution",
            position: Vector3::new(1.2047460759903845, 0.6820351942449735, -0.09343868944527914),
            velocity: Vector3::new(
                -0.006617740910480954,
                0.009290304236214792,
                0.0007113219381064047,
            ),
            dt: 400.0,
        },
    ]
}

fn bench_scenarios(c: &mut Criterion) {
    let mut group = c.benchmark_group("propagate_universal/scenarios");
    for scenario in scenarios() {
        group.bench_function(scenario.name, |b| {
            b.iter(|| {
                black_box(propagate_universal(
                    black_box(&scenario.position),
                    black_box(&scenario.velocity),
                    black_box(T0),
                    black_box(T0 + scenario.dt),
                    black_box(default_solver_type()),
                ))
            })
        });
    }
    group.finish();
}

fn bench_solver_kind(c: &mut Criterion) {
    let mut group = c.benchmark_group("propagate_universal/solver_kind");
    let all_scenarios = scenarios();
    let representative = [
        &all_scenarios[0], // real_fink_fat_state: typical short arc
        &all_scenarios[7], // gap_400_days_multi_revolution: hardest convergence case
    ];

    for scenario in representative {
        for kind in [
            SolverKind::NewtonRaphson,
            SolverKind::BrentDecker,
            SolverKind::Auto,
        ] {
            let mut solver_type = default_solver_type();
            solver_type.kind = kind;
            let id = format!("{}/{kind:?}", scenario.name);
            group.bench_function(id, |b| {
                b.iter(|| {
                    black_box(propagate_universal(
                        black_box(&scenario.position),
                        black_box(&scenario.velocity),
                        black_box(T0),
                        black_box(T0 + scenario.dt),
                        black_box(solver_type),
                    ))
                })
            });
        }
    }
    group.finish();
}

/// 20-step daily-cadence chain, as a Kalman filter would run it, starting
/// from the real fink-fat state (`tests_propagation::test_propag`).
fn bench_warm_start_chain(c: &mut Criterion) {
    let mut group = c.benchmark_group("propagate_universal/warm_start_chain");

    let start_position = Vector3::new(
        -8.264_959_160_036_185e-1,
        3.919_660_608_486_096_3e-1,
        2.229_919_607_182_842_5e-2,
    );
    let start_velocity = Vector3::new(
        -5.447_367_111_934_2e-3,
        -2.107_596_146_728_544e-2,
        1.560_811_152_125_889_6e-3,
    );
    const STEPS: usize = 20;
    const STEP_DT: f64 = 1.0;

    group.bench_function("cold", |b| {
        b.iter(|| {
            let mut position = start_position;
            let mut velocity = start_velocity;
            let mut t = T0;
            for _ in 0..STEPS {
                let solver_type = default_solver_type(); // psi_guess always None
                let result = propagate_universal(
                    black_box(&position),
                    black_box(&velocity),
                    t,
                    t + STEP_DT,
                    solver_type,
                )
                .expect("chain step should converge");
                position = result.r1;
                velocity = result.v1;
                t += STEP_DT;
            }
            black_box((position, velocity))
        })
    });

    group.bench_function("warm", |b| {
        b.iter(|| {
            let mut position = start_position;
            let mut velocity = start_velocity;
            let mut t = T0;
            let mut psi_guess = None;
            for _ in 0..STEPS {
                let mut solver_type = default_solver_type();
                solver_type.params.psi_guess = psi_guess;
                let result = propagate_universal(
                    black_box(&position),
                    black_box(&velocity),
                    t,
                    t + STEP_DT,
                    solver_type,
                )
                .expect("chain step should converge");
                position = result.r1;
                velocity = result.v1;
                psi_guess = Some(result.psy);
                t += STEP_DT;
            }
            black_box((position, velocity))
        })
    });

    group.finish();
}

/// Isolated cost of the two functions called on every Newton/Brent-Dekker
/// iteration inside `propagate_universal`: the initial-guess heuristic
/// (`prelim_kepuni`, exercised here via its public `prelim_elliptic` /
/// `prelim_hyperbolic` branches) and the Stumpff-like series `s_funct`.
fn bench_components(c: &mut Criterion) {
    let mu = GAUSS_GRAV * GAUSS_GRAV;
    let mut group = c.benchmark_group("kepler/components");

    let typical_elliptic_params = UniversalKeplerParams {
        dt: 1.992476899875328,
        r0: 0.9150028121122286,
        sig0: -0.21648695899581272,
        mu,
        alpha: -0.554947829724638, // -1/a, a~1.8 AU
        e0: 0.5005485966968782,
        solver_type: default_solver_type(),
    };
    let near_parabolic_elliptic_params = UniversalKeplerParams {
        alpha: -1e-6,
        ..typical_elliptic_params
    };
    let typical_hyperbolic_params = UniversalKeplerParams {
        alpha: 0.6666666666666666, // matches the hyperbolic scenario above (a=-1.5 AU)
        r0: 1.5148744695454394,
        sig0: 0.05,
        e0: 2.0,
        ..typical_elliptic_params
    };

    group.bench_function("prelim_elliptic/typical", |b| {
        b.iter(|| black_box(prelim_elliptic(black_box(&typical_elliptic_params))))
    });
    group.bench_function("prelim_elliptic/near_parabolic", |b| {
        b.iter(|| black_box(prelim_elliptic(black_box(&near_parabolic_elliptic_params))))
    });
    group.bench_function("prelim_hyperbolic/typical", |b| {
        b.iter(|| black_box(prelim_hyperbolic(black_box(&typical_hyperbolic_params))))
    });

    // s_funct: small-|beta| power-series branch vs large-|beta|
    // halving-and-duplication branch (see src/kepler/stumpff.rs).
    group.bench_function("s_funct/small_beta", |b| {
        b.iter(|| black_box(s_funct(black_box(2.1870299805300064), black_box(-0.5549))))
    });
    group.bench_function("s_funct/large_beta", |b| {
        b.iter(|| black_box(s_funct(black_box(1258.2351352918577), black_box(-0.5549))))
    });

    group.finish();
}

criterion_group!(
    benches,
    bench_scenarios,
    bench_solver_kind,
    bench_warm_start_chain,
    bench_components
);
criterion_main!(benches);
