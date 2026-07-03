use std::f64::consts::PI;

use crate::kepler::{principal_angle, UniversalKeplerParams};

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
pub fn prelim_elliptic(params: &UniversalKeplerParams) -> f64 {
    let contr = params.solver_type.params.convergency;
    let max_iter = params.solver_type.params.max_iter_prelim_kepuni;

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
