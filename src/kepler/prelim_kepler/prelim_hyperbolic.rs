use crate::kepler::UniversalKeplerParams;

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
pub fn prelim_hyperbolic(params: &UniversalKeplerParams) -> f64 {
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
    let updated_hyperbolic_anomaly = solve_hyperbolic_kepler_equation(
        target_mean_anomaly,
        params.e0,
        params.solver_type.params.convergency,
        params.solver_type.params.max_iter_prelim_kepuni,
    );

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
