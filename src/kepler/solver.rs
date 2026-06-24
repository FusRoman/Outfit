//! Newton–Raphson solver for the universal Kepler equation.

use super::orbit_type::OrbitType;
use super::params::UniversalKeplerParams;
use super::preliminary_guess::prelim_kepuni;
use super::stumpff::s_funct;

/// Solve the **universal Kepler equation** with a safeguarded Newton iteration.
///
/// Goal
/// -----
/// Find the universal anomaly `ψ` such that the residual
///
/// ```text
/// f(ψ) = r0·s1(ψ, α) + sig0·s2(ψ, α) + mu·s3(ψ, α) − dt = 0
/// ```
///
/// vanishes, where `(s0, s1, s2, s3)` are the Stumpff-like auxiliary functions
/// (see [`s_funct`](crate::kepler::s_funct)) for the energy parameter `α = 2E`.
///
/// Scientific background
/// ---------------------
/// The **universal-variable formulation** unifies elliptic and hyperbolic motion
/// (and the parabolic limit) in a single equation in the universal anomaly `ψ`.
/// The Stumpff-like functions `(s0..s3)` appear in the **Lagrange f–g series**,
/// enabling position/velocity propagation across conic regimes by combining
/// geometry (`r0`, `sig0`) with dynamics (`mu`, `α`) and the time of flight `dt`.
///
/// Supported regimes
/// -----------------
/// * **Elliptic** (`α < 0`) and **Hyperbolic** (`α > 0`) motions.
/// * **Parabolic** (`α = 0`) is **not** handled here and returns `None`.
///
/// Numerical strategy
/// ------------------
/// * Initial guess:
///   * Use [`prelim_kepuni`](crate::kepler::prelim_kepuni) unless an explicit `psi_guess` is provided (warm start).
/// * Newton–Raphson with safeguards:
///   * **Derivative guard**: if `f′(ψ)` is tiny/non-finite, shrink `ψ` (back-off).
///   * **Step limiter**: cap `|Δψ|` by a multiple of `1 + |ψ|` (stability on hyperbolas).
///   * **Sign-change damping**: if the update would flip `sign(ψ)`, halve `ψ` instead
///     to avoid zero-crossing oscillations.
/// * Convergence checks:
///   * **Residual criterion**: `|f(ψ)| ≤ tol_f` with a scale-aware tolerance
///     `tol_f = atol + rtol·|dt|` (robust for both small and large times of flight).
///   * **Step criterion**: `|Δψ| ≤ tol_step` (absolute) or `≤ tol_step·(1 + |ψ|)` (relative).
///
/// Arguments
/// -----------------
/// * `params` – Packed parameters for the universal formulation:
///   - `r0` (initial radius), `sig0` (radial-velocity proxy), `mu` (GM),
///     `alpha = 2E`, and `dt` (time of flight).
/// * `convergency` – Optional absolute tolerance on `ψ` (default: `100·ε`).
/// * `psi_guess` – Optional initial guess for `ψ` (skips [`prelim_kepuni`](crate::kepler::prelim_kepuni) if provided).
///
/// Return
/// ----------
/// * `Some((psi, s0, s1, s2, s3))` on convergence with `s⋅` consistent with the final `ψ`.
/// * `None` on failure to converge within the iteration budget or for `α = 0` (parabolic).
///
/// Notes
/// ------
/// * The step limiter `|Δψ| ≤ k·(1 + |ψ|)` with `k = 2` curbs runaway updates,
///   especially useful in **hyperbolic** cases where `ψ` can grow rapidly.
/// * The residual tolerance mixes absolute and relative components to keep
///   sensitivity proportional to the scale of `dt`.
///
/// See also
/// ------------
/// * [`s_funct`](crate::kepler::s_funct) – Stumpff-like auxiliary functions used to evaluate `f(ψ)` and `f′(ψ)`.
/// * [`prelim_kepuni`](crate::kepler::prelim_kepuni) – Heuristic initial guess for the universal anomaly.
/// * [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams) – Input container for the universal-variable formulation.
#[inline(always)]
pub fn solve_kepuni_with_guess(
    params: &UniversalKeplerParams,
    convergency: Option<f64>,
    psi_guess: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    println!("\n solve_kepuni_with_guess CALL ");
    println!("Input parameters: \n");
    println!("params : \n{params:?}\n");
    println!("convergency: {convergency:?}");
    println!("psi_guess: {psi_guess:?}");
    println!("----\n");

    const MAX_NEWTON_ITERATIONS: usize = 50;

    // Absolute step tolerance (controls convergence on delta_psi).
    let step_tolerance = convergency.unwrap_or(100.0 * f64::EPSILON);

    println!("Newton step tolerance: {}", step_tolerance);

    // Scale-aware residual tolerance for |f(psi)|: tol_f = atol + rtol * |dt|.
    // This keeps the criterion meaningful for both small and large flight times.
    let residual_tolerance = {
        let ten_times_epsilon = 10.0 * f64::EPSILON;
        ten_times_epsilon + ten_times_epsilon * params.dt.abs()
    };

    // Maximum allowed step length relative to the current |psi|.
    // Prevents overshooting, especially on hyperbolic trajectories.
    let max_relative_step = 2.0;

    // Orbit-type gate: parabolic motion is not supported by this solver.

    println!("detected orbit type: {:?}", params.orbit_type());

    match params.orbit_type() {
        OrbitType::Parabolic => return None,
        OrbitType::Elliptic | OrbitType::Hyperbolic => {}
    }

    // Shorthands to avoid repeated field loads in the inner loop.
    let initial_radius = params.r0;
    let radial_velocity_proxy = params.sig0;
    let gravitational_parameter = params.mu;
    let energy_parameter = params.alpha;
    let time_of_flight = params.dt;

    // Initial guess for psi: warm start, or the heuristic preliminary solver.
    let mut universal_anomaly = if let Some(guess) = psi_guess {
        guess
    } else {
        // May return None if the guess cannot be built consistently.
        prelim_kepuni(params, step_tolerance)?
    };

    println!("initial computed universal anomaly: {}", universal_anomaly);
    println!("starting newton loop ...");

    for _ in 0..MAX_NEWTON_ITERATIONS {
        println!("\nNewton universal anomaly: {}", universal_anomaly);

        // Defensive: if psi became non-finite (rare in practice), back off and continue.
        if !universal_anomaly.is_finite() {
            universal_anomaly = 0.5; // arbitrary finite reset that lets the iteration recover
        }

        // Stumpff functions at the current psi.
        // Derivative chain: d(s1)/d(psi) = s0, d(s2)/d(psi) = s1, d(s3)/d(psi) = s2.
        let (s0, s1, s2, s3) = s_funct(universal_anomaly, energy_parameter);

        // Residual and derivative:
        // f(psi)  = r0*s1 + sig0*s2 + mu*s3 - dt
        // f'(psi) = r0*s0 + sig0*s1 + mu*s2
        let residual = (((initial_radius * s1) + (radial_velocity_proxy * s2))
            + (gravitational_parameter * s3))
            - time_of_flight;
        let residual_derivative =
            ((initial_radius * s0) + (radial_velocity_proxy * s1)) + (gravitational_parameter * s2);

        println!("residual: {}", residual);
        println!("residual derivative: {}", residual_derivative);
        println!(
            "first exit loop condition: \n{} <= {}",
            residual.abs(),
            residual_tolerance
        );

        // --- Convergence on residual
        if residual.abs() <= residual_tolerance {
            return Some((universal_anomaly, s0, s1, s2, s3));
        }

        // --- Derivative safeguard
        // If f' is too small or non-finite, a Newton step would be unreliable -> shrink psi.
        if !residual_derivative.is_finite() || residual_derivative.abs() < 10.0 * f64::EPSILON {
            universal_anomaly *= 0.5;
            continue;
        }

        // --- Newton step (with step limiting)
        let mut newton_step = -residual / residual_derivative;

        // Limit the step magnitude: |delta_psi| <= max_relative_step * (1 + |psi|).
        let max_step_length = max_relative_step * (1.0 + universal_anomaly.abs());

        println!(
            "second exit loop condition: \n{} > {}",
            newton_step.abs(),
            max_step_length
        );

        if newton_step.abs() > max_step_length {
            newton_step = newton_step.signum() * max_step_length;
        }

        // --- Sign-change damping
        // If the update would cross zero (flip the sign of psi), halve psi instead
        // to avoid oscillations around 0 (common instability pattern).
        let candidate_anomaly = universal_anomaly + newton_step;
        universal_anomaly = if candidate_anomaly * universal_anomaly < 0.0 {
            0.5 * universal_anomaly
        } else {
            candidate_anomaly
        };

        // --- Convergence on step size
        // Absolute criterion: when |delta_psi| is already tiny, reuse current s.. to avoid one extra call.
        let step_is_small_in_absolute_terms = newton_step.abs() <= step_tolerance;

        println!(
            "third exit loop condition: \n{} <= {}",
            newton_step.abs(),
            step_tolerance
        );

        if step_is_small_in_absolute_terms {
            return Some((universal_anomaly, s0, s1, s2, s3));
        }

        // Relative-only success path: ensure strict consistency by recomputing s.. at the final psi.
        let step_is_small_in_relative_terms =
            newton_step.abs() <= step_tolerance * (1.0 + universal_anomaly.abs());

        println!(
            "fourth exit loop condition: \n{} <= {}",
            newton_step.abs(),
            step_tolerance * (1.0 + universal_anomaly.abs())
        );
        if step_is_small_in_relative_terms {
            let (s0_final, s1_final, s2_final, s3_final) =
                s_funct(universal_anomaly, energy_parameter);
            return Some((universal_anomaly, s0_final, s1_final, s2_final, s3_final));
        }
    }

    println!("No convergency");

    // No convergence within MAX_NEWTON_ITERATIONS.
    None
}

/// Solve the universal Kepler equation with the legacy two-argument signature.
///
/// This is a thin wrapper around [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) that omits the
/// optional warm-start argument (`psi_guess`). It exists for backward
/// compatibility with earlier code paths that only provide the orbital
/// parameters and a convergence tolerance.
///
/// Behavior
/// --------
/// * Delegates directly to [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) with `psi_guess = None`.
/// * Otherwise identical in functionality and return type.
/// * Use [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) if you want to supply a warm-start
///   value for the universal anomaly `ψ` to speed up iterative refinement.
///
/// Arguments
/// ---------
/// * `params` – Packed orbital parameters for the universal-variable formulation.
/// * `convergency` – Optional absolute tolerance on `ψ` steps
///   (default inside [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess): `100 × ε`).
///
/// Return
/// ------
/// * `Some((psi, s0, s1, s2, s3))` if the solver converged.
/// * `None` if convergence fails or if the orbit is parabolic (`alpha = 0`).
///
/// See also
/// --------
/// * [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) – Extended variant with warm-start support.
/// * [`s_funct`](crate::kepler::s_funct) – Computes the Stumpff functions `(s0..s3)`.
pub fn solve_kepuni(
    params: &UniversalKeplerParams,
    convergency: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    // Delegate directly, disabling warm-start by passing None for psi_guess.
    solve_kepuni_with_guess(params, convergency, None)
}

#[cfg(test)]
mod tests_solve_kepuni {
    use approx::assert_relative_eq;
    use proptest::prelude::*;

    use super::{solve_kepuni, UniversalKeplerParams};

    const MU: f64 = 1.0;
    const CONTR: f64 = 1e-12;

    fn make_params(
        dt: f64,
        r0: f64,
        sig0: f64,
        mu: f64,
        alpha: f64,
        e0: f64,
    ) -> UniversalKeplerParams {
        UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu,
            alpha,
            e0,
        }
    }

    #[test]
    fn test_solve_kepuni_returns_none_for_alpha_zero() {
        let params = make_params(1.0, 1.0, 0.1, MU, 0.0, 0.1);
        let res = solve_kepuni(&params, Some(CONTR));
        assert!(res.is_none());
    }

    #[test]
    fn test_solve_kepuni_elliptical_nominal() {
        let params = make_params(0.1, 1.0, 0.1, MU, -1.0, 0.5);
        let res = solve_kepuni(&params, Some(CONTR));
        assert!(
            res.is_some(),
            "solve_kepuni should converge for elliptical orbit"
        );
        let (psi, s0, s1, s2, s3) = res.unwrap();

        assert!(psi.is_finite());
        assert!(s0.is_finite());
        assert!(s1.is_finite());
        assert!(s2.is_finite());
        assert!(s3.is_finite());
    }

    #[test]
    fn test_solve_kepuni_hyperbolic_nominal() {
        let params = make_params(0.1, 2.0, -0.1, MU, 1.0, 1.5);
        let res = solve_kepuni(&params, Some(CONTR));
        assert!(
            res.is_some(),
            "solve_kepuni should converge for hyperbolic orbit"
        );
        let (psi, s0, s1, s2, s3) = res.unwrap();

        assert!(psi.is_finite());
        assert!(s0.is_finite());
        assert!(s1.is_finite());
        assert!(s2.is_finite());
        assert!(s3.is_finite());
    }

    #[test]
    fn test_solve_kepuni_large_dt_still_converges() {
        let params = make_params(10.0, 1.0, 0.1, MU, -1.0, 0.5);
        let res = solve_kepuni(&params, Some(CONTR));
        assert!(res.is_some(), "solve_kepuni should converge for long dt");
    }

    #[test]
    fn test_solve_kepuni_no_convergency_param_uses_default() {
        let params = make_params(0.5, 1.0, 0.2, MU, -1.0, 0.3);
        let res = solve_kepuni(&params, None);
        assert!(
            res.is_some(),
            "solve_kepuni should converge even with default tolerance"
        );
    }

    #[test]
    fn test_solve_kepuni_real_value() {
        let params = make_params(
            -20.765849999996135,
            1.3803870211345761,
            3.701_354_484_003_874_8E-3,
            2.959_122_082_855_911_5E-4,
            -1.642_158_377_771_140_7E-4,
            0.283_599_599_137_344_5,
        );

        let (psi, s0, s1, s2, s3) = solve_kepuni(&params, None).unwrap();

        assert_eq!(psi, -15.327414893041848);
        assert_eq!(s0, 0.9807723505583343);
        assert_eq!(s1, -15.229051668919967);
        assert_eq!(s2, 117.0876676813769);
        assert_eq!(s3, -598.9874390519309);

        let params2 = UniversalKeplerParams {
            alpha: 1.642_158_377_771_140_7E-4,
            ..params
        };
        let (psi, s0, s1, s2, s3) = solve_kepuni(&params2, None).unwrap();

        assert_eq!(psi, -15.1324122746124);
        assert_eq!(s0, 1.0188608766146905);
        assert_eq!(s1, -15.227430038021337);
        assert_eq!(s2, 114.854187452308);
        assert_eq!(s3, -578.615100072754);
    }

    #[test]
    fn test_solve_kepuni_invariant_residual_elliptical() {
        let params = make_params(0.1, 1.0, 0.1, MU, -1.0, 0.5);
        let res = solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");
        let (_, _, s1, s2, s3) = res;

        let residual = params.r0 * s1 + params.sig0 * s2 + MU * s3 - params.dt;
        assert_relative_eq!(residual, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn test_solve_kepuni_invariant_residual_hyperbolic() {
        let params = make_params(0.1, 2.0, -0.1, MU, 1.0, 1.5);
        let res = solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");
        let (_, _, s1, s2, s3) = res;

        let residual = params.r0 * s1 + params.sig0 * s2 + MU * s3 - params.dt;
        assert_relative_eq!(residual, 0.0, max_relative = 1e-10);
    }

    #[test]
    fn test_solve_kepuni_known_values_elliptical() {
        let params = make_params(0.5, 1.0, 0.2, MU, -1.0, 0.3);
        let (psi, _, _, _, _) = solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");

        let expected_psi = 0.47761843287737277;
        assert_relative_eq!(psi, expected_psi, epsilon = 1e-15);
    }

    #[test]
    fn test_solve_kepuni_known_values_hyperbolic() {
        let params = make_params(0.3, 2.0, -0.1, MU, 1.0, 1.5);
        let (psi, _, _, _, _) = solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");

        let expected_psi = 0.14972146123530983;
        assert_relative_eq!(psi, expected_psi, epsilon = 1e-15);
    }

    fn arb_common_params() -> impl Strategy<Value = (f64, f64, f64, f64, f64)> {
        (
            prop_oneof![-5.0..-0.1, 0.1..5.0],
            0.1..3.0f64,
            -0.5..0.5f64,
            0.5..2.0f64,
            0.01..2.0f64,
        )
    }

    proptest! {
        #[test]
        fn prop_solve_kepuni_elliptical(
            (dt, r0, sig0, mu, e0) in arb_common_params(),
            alpha in -5.0..-0.1f64
        ) {
            let params = make_params(dt, r0, sig0, mu, alpha, e0);
            let res = solve_kepuni(&params, Some(1e-12));

            if let Some((psi, s0, s1, s2, s3)) = res {
                prop_assert!(psi.is_finite());
                prop_assert!(s0.is_finite());
                prop_assert!(s1.is_finite());
                prop_assert!(s2.is_finite());
                prop_assert!(s3.is_finite());

                let residual = params.r0 * s1 + params.sig0 * s2 + mu * s3 - params.dt;
                prop_assert!(residual.abs() < 1e-8,
                    "Residual too large for elliptical case: {}", residual);
            }
        }
    }

    proptest! {
        #[test]
        fn prop_solve_kepuni_hyperbolic(
            (dt, r0, sig0, mu, e0) in arb_common_params(),
            alpha in 0.1..5.0f64
        ) {
            let params = make_params(dt, r0, sig0, mu, alpha, e0);
            let res = solve_kepuni(&params, Some(1e-12));

            if let Some((psi, s0, s1, s2, s3)) = res {
                prop_assert!(psi.is_finite());
                prop_assert!(s0.is_finite());
                prop_assert!(s1.is_finite());
                prop_assert!(s2.is_finite());
                prop_assert!(s3.is_finite());

                let residual = params.r0 * s1 + params.sig0 * s2 + mu * s3 - params.dt;
                prop_assert!(residual.abs() < 1e-8,
                    "Residual too large for hyperbolic case: {}", residual);
            }
        }
    }
}
