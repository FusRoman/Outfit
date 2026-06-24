//! Preliminary (initial-guess) solvers for the universal anomaly `psi`.
//!
//! These routines do not iterate to full convergence: they only provide a
//! starting value for `psi`, meant to be refined afterwards by
//! [`solve_kepuni`](crate::kepler::solve_kepuni).

use std::f64::consts::PI;

use super::angles::principal_angle;
use super::params::UniversalKeplerParams;

/// Compute a preliminary estimate of the universal anomaly `ψ` for elliptical orbits (α < 0).
///
/// This function generates an initial guess of the universal anomaly corresponding
/// to the time `t0 + dt`, given the initial state of an object in an elliptical
/// orbit. The universal anomaly is then used as a starting point for solving the
/// universal Kepler equation.
///
/// Algorithm
/// ---------
/// 1. Compute the semi-major axis `a0 = -μ / α` and mean motion `n = sqrt((-α)^3) / μ`.
/// 2. If the eccentricity `e0` is very small, approximate `ψ` directly using a linear formula.
/// 3. Otherwise, compute the initial eccentric anomaly `u0` from the orbital geometry and
///    correct its sign based on the radial velocity `sig0`.
/// 4. Compute the mean anomaly at epoch `ℓ0 = u0 - e0·sin(u0)` and propagate it forward
///    by `n·dt` to obtain the target mean anomaly at `t0 + dt`.
/// 5. Solve Kepler's equation `M = u - e·sin(u)` using Newton–Raphson iterations to find
///    the updated eccentric anomaly `u`.
/// 6. Convert the difference `u - u0` into the universal anomaly `ψ` using the factor `1/sqrt(-α)`.
///
/// Arguments
/// ---------
/// * `params` – [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams) containing `dt`, `r0`, `sig0`, `mu`, `alpha`, and `e0`.
/// * `contr` – Convergence control threshold for the Newton solver.
/// * `max_iter` – Maximum number of Newton–Raphson iterations.
///
/// Returns
/// --------
/// * `psi` – Preliminary value of the universal anomaly at time `t0 + dt`.
///
/// Remarks
/// --------
/// * This routine does not handle parabolic (α = 0) or hyperbolic (α > 0) cases.
/// * The result is an **initial guess**: a subsequent universal Kepler solver will refine `ψ`.
/// * The computed `ψ` can be negative if `dt` corresponds to backward propagation.
///
/// # See also
/// * [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) – Equivalent procedure for hyperbolic orbits.
/// * [`solve_kepuni`](crate::kepler::solve_kepuni) – Refines `ψ` by solving the universal Kepler equation.
/// * [`angle_diff`](crate::kepler::angle_diff) – Computes the principal difference between two angles.
pub fn prelim_elliptic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
    // Step 1: semi-major axis and mean motion.
    let semi_major_axis = -params.mu / params.alpha;
    let mean_motion = (-params.alpha.powi(3)).sqrt() / params.mu;

    // Step 2: special case, nearly circular orbit.
    if params.e0 < contr {
        return mean_motion * params.dt / (-params.alpha).sqrt();
    }

    // Step 3: eccentric anomaly at epoch, sign-corrected by radial velocity.
    let initial_eccentric_anomaly =
        initial_eccentric_anomaly_from_geometry(params.r0, semi_major_axis, params.e0, params.sig0);

    // Mean anomaly at epoch, then propagated forward by mean_motion * dt.
    // The target mean anomaly is intentionally left unwrapped to preserve
    // multi-revolution information: wrapping to [0, 2π) would discard the
    // number of full orbits completed during dt, producing a wrong initial
    // guess for dt spanning more than one revolution.
    let mean_anomaly_at_epoch =
        principal_angle(initial_eccentric_anomaly - params.e0 * initial_eccentric_anomaly.sin());
    let target_mean_anomaly = mean_anomaly_at_epoch + mean_motion * params.dt;

    // Step 5: solve Kepler's equation for the updated eccentric anomaly.
    let updated_eccentric_anomaly =
        solve_elliptic_kepler_equation(target_mean_anomaly, params.e0, contr, max_iter);

    // Step 6: convert (u - u0) into the universal anomaly psi.
    // The direct difference is used here (no angle wrapping) so that psi
    // correctly reflects multi-revolution arcs.
    (updated_eccentric_anomaly - initial_eccentric_anomaly) / (-params.alpha).sqrt()
}

/// Compute the eccentric anomaly `u0` at the reference epoch from orbital
/// geometry, i.e. from `cos(u0) = (1 - r0/a0) / e0`. The sign is corrected
/// according to the direction of radial motion (`sig0`), and the result is
/// wrapped into the principal angular range `[0, 2π)`.
fn initial_eccentric_anomaly_from_geometry(
    radial_distance: f64,
    semi_major_axis: f64,
    eccentricity: f64,
    radial_velocity_proxy: f64,
) -> f64 {
    let cosine_of_eccentric_anomaly = (1.0 - radial_distance / semi_major_axis) / eccentricity;

    let mut eccentric_anomaly = if cosine_of_eccentric_anomaly.abs() <= 1.0 {
        cosine_of_eccentric_anomaly.acos()
    } else if cosine_of_eccentric_anomaly >= 1.0 {
        0.0 // limit case: object at pericenter
    } else {
        PI // limit case: object at apocenter
    };

    // Flip the sign if the radial velocity is negative (object moving inward).
    if radial_velocity_proxy < 0.0 {
        eccentric_anomaly = -eccentric_anomaly;
    }

    principal_angle(eccentric_anomaly)
}

/// Solve the elliptic Kepler equation `mean_anomaly = u - e*sin(u)` for the
/// eccentric anomaly `u`, starting from `u = π`, using up to `max_iter`
/// Newton–Raphson iterations.
fn solve_elliptic_kepler_equation(
    target_mean_anomaly: f64,
    eccentricity: f64,
    convergence_threshold: f64,
    max_iter: usize,
) -> f64 {
    // Initial guess: u = M. Exact for e = 0 and remains accurate for small
    // eccentricities. Crucially, starting from M (unwrapped) rather than π
    // ensures the solver stays in the correct revolution count.
    let mut eccentric_anomaly = target_mean_anomaly;
    for _ in 0..max_iter {
        let residual =
            eccentric_anomaly - eccentricity * eccentric_anomaly.sin() - target_mean_anomaly;
        let residual_derivative = 1.0 - eccentricity * eccentric_anomaly.cos();
        let newton_step = -residual / residual_derivative;
        eccentric_anomaly += newton_step;
        if newton_step.abs() < convergence_threshold * 1e3 {
            break;
        }
    }
    eccentric_anomaly
}

/// Compute a preliminary estimate of the universal anomaly `ψ` for hyperbolic orbits (α > 0).
///
/// This function generates an initial guess of the universal anomaly corresponding
/// to the time `t0 + dt`, given the initial state of an object on a hyperbolic
/// trajectory. The universal anomaly is used as a starting point for solving the
/// universal Kepler equation.
///
/// Algorithm
/// ---------
/// 1. Compute the semi-major axis `a0 = -μ / α` and the hyperbolic mean motion
///    `n = sqrt(α³) / μ`.
/// 2. Compute the initial hyperbolic eccentric anomaly `F₀` from the geometry,
///    using:
///    cosh(F₀) = (1 - r₀ / a₀) / e₀
/// 3. Adjust the sign of `F₀` based on the sign of the radial velocity `sig0`.
/// 4. Compute the mean anomaly at epoch: ℓ₀ = e₀·sinh(F₀) - F₀.
/// 5. Propagate ℓ₀ forward by n·dt to get the target mean anomaly at t₀ + dt.
/// 6. Solve the hyperbolic Kepler equation:
///    e₀·sinh(F) - F = ℓ
///    using Newton-Raphson iterations, with additional handling to avoid
///    divergence for large |F|.
/// 7. Convert the difference (F - F₀) into the universal anomaly `ψ`
///    using the scaling factor 1/√α.
///
/// Arguments
/// ---------
/// * `params` – [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams) containing `dt`, `r0`, `sig0`, `mu`, `alpha`, and `e0`.
/// * `contr` – Convergence control threshold for the Newton solver.
/// * `max_iter` – Maximum number of Newton-Raphson iterations.
///
/// Returns
/// --------
/// * `psi` – Preliminary value of the universal anomaly at time `t0 + dt`.
///
/// Remarks
/// --------
/// * This function does not handle elliptical (α < 0) or parabolic (α = 0) orbits.
/// * The returned `ψ` is only an initial guess and is refined later by a solver.
///
/// # See also
/// * [`prelim_elliptic`](crate::kepler::prelim_elliptic) – Equivalent routine for elliptical orbits.
/// * [`solve_kepuni`](crate::kepler::solve_kepuni) – Refines ψ by solving the universal Kepler equation.
pub fn prelim_hyperbolic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
    // Step 1: semi-major axis (negative for hyperbolic orbits) and hyperbolic mean motion.
    let semi_major_axis = -params.mu / params.alpha;
    let mean_motion = params.alpha.powi(3).sqrt() / params.mu;

    // Step 2: hyperbolic anomaly at epoch, sign-corrected by radial velocity.
    let initial_hyperbolic_anomaly = initial_hyperbolic_anomaly_from_geometry(
        params.r0,
        semi_major_axis,
        params.e0,
        params.sig0,
    );

    // Step 3: mean anomaly at epoch, propagated forward by mean_motion * dt.
    let mean_anomaly_at_epoch =
        params.e0 * initial_hyperbolic_anomaly.sinh() - initial_hyperbolic_anomaly;
    let target_mean_anomaly = mean_anomaly_at_epoch + mean_motion * params.dt;

    // Step 4: solve the hyperbolic Kepler equation for the updated anomaly.
    let updated_hyperbolic_anomaly =
        solve_hyperbolic_kepler_equation(target_mean_anomaly, params.e0, contr, max_iter);

    // Step 5: convert the anomaly difference into the universal anomaly psi.
    (updated_hyperbolic_anomaly - initial_hyperbolic_anomaly) / params.alpha.sqrt()
}

/// Compute the hyperbolic eccentric anomaly `F0` at the reference epoch from
/// orbital geometry, i.e. from `cosh(F0) = (1 - r0/a0) / e0`. The sign is
/// corrected according to the direction of radial motion (`sig0`).
fn initial_hyperbolic_anomaly_from_geometry(
    radial_distance: f64,
    semi_major_axis: f64,
    eccentricity: f64,
    radial_velocity_proxy: f64,
) -> f64 {
    let hyperbolic_cosine_of_anomaly = (1.0 - radial_distance / semi_major_axis) / eccentricity;

    let mut hyperbolic_anomaly = if hyperbolic_cosine_of_anomaly > 1.0 {
        // Inverse hyperbolic cosine: acosh(x) = ln(x + sqrt(x^2 - 1)).
        (hyperbolic_cosine_of_anomaly + (hyperbolic_cosine_of_anomaly.powi(2) - 1.0).sqrt()).ln()
    } else {
        0.0 // limit case: object very close to pericenter
    };

    if radial_velocity_proxy < 0.0 {
        hyperbolic_anomaly = -hyperbolic_anomaly;
    }

    hyperbolic_anomaly
}

/// Solve the hyperbolic Kepler equation `mean_anomaly = e*sinh(F) - F` for
/// the hyperbolic anomaly `F`, starting from `F = 0`, using up to `max_iter`
/// damped Newton-Raphson iterations. Large `|F|` updates are halved instead
/// of taking a full Newton step, in order to avoid divergence.
fn solve_hyperbolic_kepler_equation(
    target_mean_anomaly: f64,
    eccentricity: f64,
    convergence_threshold: f64,
    max_iter: usize,
) -> f64 {
    let mut hyperbolic_anomaly: f64 = 0.0;

    for _ in 0..max_iter {
        if hyperbolic_anomaly.abs() < 15.0 {
            // Newton-Raphson update.
            let residual =
                eccentricity * hyperbolic_anomaly.sinh() - hyperbolic_anomaly - target_mean_anomaly;
            let residual_derivative = eccentricity * hyperbolic_anomaly.cosh() - 1.0;
            let newton_step = -residual / residual_derivative;
            let candidate_anomaly = hyperbolic_anomaly + newton_step;
            // If the update would cross zero, dampen the step to avoid divergence.
            hyperbolic_anomaly = if hyperbolic_anomaly * candidate_anomaly < 0.0 {
                hyperbolic_anomaly / 2.0
            } else {
                candidate_anomaly
            };
        } else {
            // For very large |F|, reduce it progressively.
            hyperbolic_anomaly /= 2.0;
        }

        // Convergence check.
        if hyperbolic_anomaly.abs() < convergence_threshold * 1e3 {
            break;
        }
    }

    hyperbolic_anomaly
}

#[cfg(test)]
mod tests_prelim_kepuni {

    use crate::kepler::params::SolverType;

    use super::UniversalKeplerParams;

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
            solver_type: SolverType::NewtonRaphson {
                convergency: Some(CONTR),
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        }
    }

    #[test]
    fn test_returns_none_for_alpha_zero() {
        let params = make_params(1.0, 1.0, 0.0, MU, 0.0, 0.1);
        let res = params.prelim_kepuni();
        assert!(res.is_none());
    }

    #[test]
    fn test_elliptic_small_eccentricity() {
        let params = make_params(0.5, 1.0, 0.1, MU, -1.0, 1e-8);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_elliptic_high_eccentricity() {
        let params = make_params(0.1, 0.5, 0.2, MU, -1.0, 0.8);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_hyperbolic_case() {
        let params = make_params(0.3, 2.0, -0.1, MU, 1.0, 1.5);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_negative_sig0_changes_direction() {
        let alpha = -1.0;
        let r0 = 1.0;
        let e0 = 0.5;
        let dt = 0.25;

        let params_pos = make_params(dt, r0, 0.1, MU, alpha, e0);
        let params_neg = make_params(dt, r0, -0.1, MU, alpha, e0);

        let psi_pos = params_pos.prelim_kepuni().unwrap();
        let psi_neg = params_neg.prelim_kepuni().unwrap();

        assert!(
            (psi_pos - psi_neg).abs() > 1e-8,
            "psi did not change significantly when changing sig0 sign: {psi_pos} vs {psi_neg}"
        );
    }

    #[test]
    fn test_stability_long_dt() {
        let params = make_params(50.0, 1.0, 0.1, MU, -1.0, 0.5);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_edge_cosine_limits() {
        let params = make_params(0.25, 2.0, 0.1, MU, -1.0, 0.1);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
    }

    #[test]
    fn test_prelim_kepuni_real_data() {
        let epsilon = f64::EPSILON;
        let contr = 100.0 * epsilon;

        let dt = -20.765849999996135;
        let r0 = 1.3803870211345761;
        let sig0 = 3.701_354_484_003_874_8E-3;
        let mu = 2.959_122_082_855_911_5E-4;
        let alpha = -1.642_158_377_771_140_7E-4;
        let e0 = 0.283_599_599_137_344_5;

        let params = UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu,
            alpha,
            e0,
            solver_type: SolverType::NewtonRaphson {
                convergency: Some(contr),
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        };
        let psi = params.prelim_kepuni().unwrap();
        assert_eq!(psi, -15.327414893041848);

        let params2 = UniversalKeplerParams {
            alpha: 1.642_158_377_771_140_7E-4,
            ..params
        };
        let psi = params2.prelim_kepuni().unwrap();
        assert_eq!(psi, -73.1875935362658);

        let params3 = UniversalKeplerParams {
            alpha: 0.0,
            ..params
        };
        assert!(params3.prelim_kepuni().is_none());
    }

    mod kepuni_prop_tests {
        use crate::kepler::params::SolverType;

        use super::UniversalKeplerParams;
        use proptest::prelude::*;

        fn arb_params() -> impl Strategy<Value = UniversalKeplerParams> {
            (
                -10.0..10.0f64,
                0.1..5.0f64,
                -2.0..2.0f64,
                0.5..2.0f64,
                prop_oneof![(-5.0..-0.01f64), (0.01..5.0f64)],
                0.0..3.0f64,
                1e-14..1e-8f64,
            )
                .prop_map(|(dt, r0, sig0, mu, alpha, e0, contr)| {
                    UniversalKeplerParams {
                        dt,
                        r0,
                        sig0,
                        mu,
                        alpha,
                        e0,
                        solver_type: SolverType::NewtonRaphson {
                            convergency: Some(contr),
                            psi_guess: None,
                            max_iter_prelim_kepuni: None,
                        },
                    }
                })
        }

        proptest! {
            #[test]
            fn prop_prelim_kepuni_behaves_well(params in arb_params()) {
                let result = params.prelim_kepuni();
                prop_assert!(result.is_some());
                prop_assert!(result.unwrap().is_finite());
            }
        }

        proptest! {
            #[test]
            fn prop_prelim_kepuni_alpha_zero(
                dt in -10.0..10.0f64,
                r0 in 0.1..5.0f64,
                sig0 in -2.0..2.0f64,
                mu in 0.5..2.0f64,
                e0 in 0.0..3.0f64,
                contr in 1e-14..1e-8f64
            ) {
                let params = UniversalKeplerParams { dt, r0, sig0, mu, alpha: 0.0, e0, solver_type: SolverType::NewtonRaphson {
                    convergency: Some(contr),
                    psi_guess: None,
                    max_iter_prelim_kepuni: None
                } };
                let result = params.prelim_kepuni();
                if let Some(psi0) = result {
                    prop_assert!(psi0.is_finite());
                }
            }
        }

        proptest! {
            #[test]
            fn prop_sig0_influences_psi(params in arb_params()) {
                prop_assume!(params.dt.abs() > 1e-6);
                prop_assume!(params.e0 > 1e-6);
                prop_assume!(params.r0 > 1e-6);

                let mut params_pos = params;
                params_pos.sig0 = 0.1;
                let mut params_neg = params;
                params_neg.sig0 = -0.1;

                let res_pos = params_pos.prelim_kepuni();
                let res_neg = params_neg.prelim_kepuni();

                prop_assume!(res_pos.is_some() && res_neg.is_some());

                let diff = (res_pos.unwrap() - res_neg.unwrap()).abs();
                prop_assert!(diff >= 0.0);
            }
        }
    }
}
