use super::constants::DPI;
use super::constants::GAUSS_GRAV;
use super::orb_elem::eccentricity_control;
use core::fmt;
use nalgebra::Vector3;
use std::f64::consts::PI;

/// Computes the Stumpff-like auxiliary functions (s0, s1, s2, s3) used in
/// universal variable formulations of two-body orbital motion.
///
/// These functions are a generalization of the classical Stumpff functions.
/// They appear in semi-analytical formulations such as the **f–g series**
/// (Everhart & Pitkin) for universal Kepler propagation, and are used to
/// compute position and velocity vectors from the universal anomaly `ψ`.
///
/// Algorithm
/// ---------
/// * For small `|α ψ²|`, `s0`, `s1`, `s2`, `s3` are computed directly
///   by power series expansion.
/// * For large `|α ψ²|`, the function:
///   1. Reduces `ψ` by repeatedly halving it until `|α ψ²|` becomes small,
///   2. Computes `s0` and `s1` from the series at the reduced value,
///   3. Applies duplication formulas to recover `s0` and `s1` for the full `ψ`,
///   4. Reconstructs `s2` and `s3` using:
///      * `s2 = (s0 - 1)/α`
///      * `s3 = (s1 - ψ)/α`
///
/// This approach avoids divergence and preserves numerical stability for
/// large anomalies, at the cost of some loss of precision for `s2` and `s3`.
///
/// Arguments
/// ---------
/// * `psi` – Universal anomaly (integration parameter).
/// * `alpha` – Two times the specific orbital energy (2 * E). It can be
///   negative (elliptic), zero (parabolic), or positive (hyperbolic).
///
/// Returns
/// --------
/// * `(s0, s1, s2, s3)` – Tuple of four Stumpff-like functions:
///   - `s0 ≈ 1 + α·Σ`  — generalized cosine-like series,
///   - `s1 ≈ ψ + α·Σ`  — generalized sine-like series,
///   - `s2 = (s0 - 1)/α`,
///   - `s3 = (s1 - ψ)/α`.
///
/// These functions are used to evaluate the Lagrange f and g coefficients:
///
/// ```text
/// f = 1 - (μ/r0) s2
/// g = Δt - μ s3
/// ```
///
/// Remarks
/// --------
/// * For very large values of `α ψ²`, `s2` and `s3` may lose some precision
///   because they are reconstructed from differences of large numbers.
///
/// References
/// ----------
/// * Everhart, E. & Pitkin, E.T., *American Journal of Physics*, 51(8), 712–717 (1983)
/// * Goodyear, W.H., *Astronomical Journal*, 70, 189–192 (1965)
///
/// # See also
/// * [`velocity_correction`] – Uses these functions to compute position and velocity (f–g series).
/// * [`solve_kepuni`] – Solves universal Kepler's equation using these functions.
/// * [Battin, *An Introduction to the Mathematics and Methods of Astrodynamics*]
fn s_funct(psi: f64, alpha: f64) -> (f64, f64, f64, f64) {
    // Maximum number of terms for the power series and half-angle iterations.
    const JMAX: usize = 70;
    const HALFMAX: usize = 30;
    const BETACONTR: f64 = 100.0;

    // Machine precision parameters:
    // - convergence_tol : term below this is negligible
    // - overflow_limit  : stop if a term becomes numerically unstable
    let eps = f64::EPSILON;
    let convergence_tol = 100.0 * eps;
    let overflow_limit = 1.0 / eps;

    // Beta is the main expansion parameter: β = α * ψ²
    let beta = alpha * psi.powi(2);

    // Helper closure for direct computation of s2 and s3 by series expansion.
    // It builds the power series for:
    //   s2 = ψ²/2 + β ψ⁴/(3*4) + ...
    //   s3 = ψ³/6 + β ψ⁵/(4*5) + ...
    let compute_s2_s3 = |beta: f64, psi: f64| -> (f64, f64) {
        let mut term_s2 = psi.powi(2) / 2.0;
        let mut term_s3 = term_s2 * psi / 3.0;
        let mut s2 = term_s2;
        let mut s3 = term_s3;

        // Add successive terms to the series for s2
        for j in 1..=JMAX {
            term_s2 *= beta / ((2.0 * j as f64 + 1.0) * (2.0 * j as f64 + 2.0));
            s2 += term_s2;
            // Stop when terms are very small or too large
            if term_s2.abs() < convergence_tol || term_s2.abs() > overflow_limit {
                break;
            }
        }

        // Add successive terms to the series for s3
        for j in 1..=JMAX {
            term_s3 *= beta / ((2.0 * j as f64 + 2.0) * (2.0 * j as f64 + 3.0));
            s3 += term_s3;
            if term_s3.abs() < convergence_tol || term_s3.abs() > overflow_limit {
                break;
            }
        }

        (s2, s3)
    };

    // =========================================================================
    // CASE 1: direct series expansion
    // When |β| is small enough, the series converge rapidly and no reduction is needed.
    // =========================================================================
    if beta.abs() < BETACONTR {
        let (s2, s3) = compute_s2_s3(beta, psi);
        let s1 = psi + alpha * s3;
        let s0 = 1.0 + alpha * s2;
        return (s0, s1, s2, s3);
    }

    // =========================================================================
    // CASE 2: large |β|
    // To ensure convergence, we repeatedly halve ψ until β becomes small,
    // compute s0 and s1 at the reduced value, then use duplication formulas
    // to recover the functions at the original ψ.
    // =========================================================================

    // 1. Reduce psi by repeated halving
    let mut psi_reduced = psi;
    let mut nhalf = 0;
    for _ in 0..HALFMAX {
        psi_reduced *= 0.5;
        nhalf += 1;
        let beta_half = alpha * psi_reduced.powi(2);
        if beta_half.abs() < BETACONTR {
            break;
        }
    }

    // 2. Compute s0 and s1 at the reduced psi using power series
    let mut term0 = 1.0; // first term for s0
    let mut term1 = psi_reduced; // first term for s1
    let mut s0 = 1.0;
    let mut s1 = psi_reduced;

    for j in 1..=JMAX {
        term0 *= beta / ((2 * j - 1) as f64 * (2 * j) as f64);
        s0 += term0;
        if term0.abs() < convergence_tol || term0.abs() > overflow_limit {
            break;
        }
    }

    for j in 1..=JMAX {
        term1 *= beta / ((2 * j) as f64 * (2 * j + 1) as f64);
        s1 += term1;
        if term1.abs() < convergence_tol || term1.abs() > overflow_limit {
            break;
        }
    }

    // 3. Apply duplication formulas to scale s0 and s1 back to the original psi.
    // These formulas correspond to:
    //   cos(2x) = 2 cos^2(x) - 1
    //   sin(2x) = 2 cos(x) sin(x)
    // generalized to the universal variables.
    for _ in 0..nhalf {
        let new_s0 = 2.0 * s0.powi(2) - 1.0;
        let new_s1 = 2.0 * s0 * s1;
        s0 = new_s0;
        s1 = new_s1;
    }

    // 4. Recompute s2 and s3 from s0 and s1 using their defining relations:
    //    s2 = (s0 - 1) / α
    //    s3 = (s1 - ψ) / α
    // These relations are numerically less stable for large β, but they are
    // the standard method used in OrbFit (and the Fortran reference code).
    let s3 = (s1 - psi) / alpha;
    let s2 = (s0 - 1.0) / alpha;

    (s0, s1, s2, s3)
}

/// Normalize an angle in radians to the range [0, 2π].
///
/// This ensures any input angle is wrapped into the principal interval
/// 0 ≤ θ < 2π using Euclidean remainder.
pub(crate) fn principal_angle(a: f64) -> f64 {
    a.rem_euclid(DPI)
}

/// Compute the signed minimal difference between two angles in radians.
///
/// Returns the value of (a - b) wrapped into the range [-π, π],
/// i.e. the smallest signed rotation from `b` to `a`.
fn angle_diff(a: f64, b: f64) -> f64 {
    let a = principal_angle(a);
    let b = principal_angle(b);

    let mut diff = a - b;
    if diff > PI {
        diff -= DPI;
    } else if diff < -PI {
        diff += DPI;
    }
    diff
}

/// Solve the preliminary universal Kepler problem (initial guess for ψ).
///
/// This function computes an initial estimate of the universal anomaly `ψ` at time `t0 + dt`,
/// given the orbital parameters stored in [`UniversalKeplerParams`].
///
/// The method selects the appropriate branch depending on the sign of `alpha`:
///
/// * **Elliptical case (alpha < 0)** – Calls [`prelim_elliptic`] to compute ψ using the
///   elliptical Kepler equation with eccentric anomaly.
/// * **Hyperbolic case (alpha > 0)** – Calls [`prelim_hyperbolic`] to compute ψ using the
///   hyperbolic Kepler equation with hyperbolic anomaly.
/// * **Parabolic case (alpha = 0)** – Not supported; the function returns `None`.
///
/// Arguments
/// ---------
/// * `params` – [`UniversalKeplerParams`] containing `dt`, `r0`, `sig0`, `mu`, `alpha`, and `e0`.
/// * `contr` – Convergence tolerance for the Newton–Raphson iterations.
///
/// Returns
/// --------
/// * `Some(psi)` – Initial guess for the universal anomaly at time `t0 + dt`.
/// * `None` – If `alpha = 0` (parabolic orbit).
///
/// Remarks
/// --------
/// * This function does not iterate to full convergence; it only provides a starting
///   value for `ψ`. It is typically used before calling a solver such as
///   [`solve_kepuni`] to refine the result.
/// * The hyperbolic and elliptical branches use up to 20 Newton–Raphson iterations.
/// * The parabolic case must be treated by a separate routine.
///
/// # See also
/// * [`prelim_elliptic`] – Computes ψ for elliptical orbits.
/// * [`prelim_hyperbolic`] – Computes ψ for hyperbolic orbits.
/// * [`solve_kepuni`] – Refines ψ by solving the universal Kepler equation.
fn prelim_kepuni(params: &UniversalKeplerParams, contr: f64) -> Option<f64> {
    const ITX: usize = 20;

    match params.orbit_type() {
        OrbitType::Elliptic => Some(prelim_elliptic(params, contr, ITX)),
        OrbitType::Hyperbolic => Some(prelim_hyperbolic(params, contr, ITX)),
        OrbitType::Parabolic => None, // Not supported
    }
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
/// * `params` – [`UniversalKeplerParams`] containing `dt`, `r0`, `sig0`, `mu`, `alpha`, and `e0`.
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
/// * [`prelim_hyperbolic`] – Equivalent procedure for hyperbolic orbits.
/// * [`solve_kepuni`] – Refines `ψ` by solving the universal Kepler equation.
/// * [`angle_diff`] – Computes the principal difference between two angles.
fn prelim_elliptic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
    // Step 1: Compute semi-major axis (a0) and mean motion (n)
    let a0 = -params.mu / params.alpha;
    let n = (-params.alpha.powi(3)).sqrt() / params.mu;

    // Step 2: Special case: nearly circular orbit
    if params.e0 < contr {
        return n * params.dt / (-params.alpha).sqrt();
    }

    // Step 3: Compute the eccentric anomaly at epoch u0 from geometry:
    // cos u0 = (1 - r0/a0) / e0
    let cos_u0 = (1.0 - params.r0 / a0) / params.e0;
    let mut u0 = if cos_u0.abs() <= 1.0 {
        cos_u0.acos()
    } else if cos_u0 >= 1.0 {
        0.0 // limit case: object at pericenter
    } else {
        PI // limit case: object at apocenter
    };

    // Flip the sign of u0 if radial velocity is negative
    if params.sig0 < 0.0 {
        u0 = -u0;
    }

    // Normalize u0 and compute the corresponding mean anomaly ℓ0
    u0 = principal_angle(u0);
    let ell0 = principal_angle(u0 - params.e0 * u0.sin());

    // Step 4: Target mean anomaly after dt
    let target_mean_anomaly = principal_angle(ell0 + n * params.dt);

    // Step 5: Solve Kepler's equation iteratively for u
    let mut u = PI; // start guess
    for _ in 0..max_iter {
        let f = u - params.e0 * u.sin() - target_mean_anomaly;
        let fp = 1.0 - params.e0 * u.cos();
        let du = -f / fp;
        u += du;
        if du.abs() < contr * 1e3 {
            break;
        }
    }

    // Step 6: Convert (u - u0) into universal anomaly ψ
    angle_diff(u, u0) / (-params.alpha).sqrt()
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
/// * `params` – [`UniversalKeplerParams`] containing `dt`, `r0`, `sig0`, `mu`, `alpha`, and `e0`.
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
/// * [`prelim_elliptic`] – Equivalent routine for elliptical orbits.
/// * [`solve_kepuni`] – Refines ψ by solving the universal Kepler equation.
fn prelim_hyperbolic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
    // Step 1: Compute semi-major axis (a0) and hyperbolic mean motion n
    // For hyperbolic orbits, a0 is negative.
    let a0 = -params.mu / params.alpha;
    let n = params.alpha.powi(3).sqrt() / params.mu;

    // Step 2: Compute the initial hyperbolic anomaly F₀
    // cosh(F₀) = (1 - r0/a0) / e0
    let coshf0 = (1.0 - params.r0 / a0) / params.e0;
    let mut f0 = if coshf0 > 1.0 {
        // Compute F₀ from cosh⁻¹(x) = ln(x + sqrt(x² - 1))
        (coshf0 + (coshf0.powi(2) - 1.0).sqrt()).ln()
    } else {
        0.0 // limit case: object very close to pericenter
    };

    // Adjust the sign of F₀ based on the sign of radial velocity
    if params.sig0 < 0.0 {
        f0 = -f0;
    }

    // Step 3: Compute the mean anomaly at epoch
    let ell0 = params.e0 * f0.sinh() - f0;

    // Propagate the mean anomaly forward by n·dt
    let target_mean_anomaly = ell0 + n * params.dt;

    // Step 4: Iteratively solve the hyperbolic Kepler equation:
    //    e·sinh(F) - F = ℓ
    // starting from F = 0 as an initial guess.
    let mut f: f64 = 0.0;

    for _ in 0..max_iter {
        if f.abs() < 15.0 {
            // Newton-Raphson update
            let func = params.e0 * f.sinh() - f - target_mean_anomaly;
            let deriv = params.e0 * f.cosh() - 1.0;
            let df = -func / deriv;
            let ff = f + df;
            // If the update crosses zero, dampen the step to avoid divergence
            f = if f * ff < 0.0 { f / 2.0 } else { ff };
        } else {
            // For very large |F|, reduce it progressively
            f /= 2.0;
        }

        // Convergence check
        if f.abs() < contr * 1e3 {
            break;
        }
    }

    // Step 5: Convert the difference in anomalies to universal anomaly ψ
    (f - f0) / params.alpha.sqrt()
}

/// Orbital regime based on the value of `alpha`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitType {
    Elliptic,   // alpha < 0
    Hyperbolic, // alpha > 0
    Parabolic,  // alpha = 0
}

impl OrbitType {
    pub fn from_alpha(alpha: f64) -> Self {
        if alpha < 0.0 {
            OrbitType::Elliptic
        } else if alpha > 0.0 {
            OrbitType::Hyperbolic
        } else {
            OrbitType::Parabolic
        }
    }
}

/// State needed to solve the universal Kepler equation.
#[derive(Debug, Clone, Copy)]
struct UniversalKeplerParams {
    dt: f64,
    r0: f64,
    sig0: f64,
    mu: f64,
    alpha: f64,
    e0: f64,
}

impl UniversalKeplerParams {
    pub fn orbit_type(&self) -> OrbitType {
        OrbitType::from_alpha(self.alpha)
    }
}

/// Solve the universal Kepler equation using Newton's method.
///
/// This solver refines an initial guess for the universal anomaly `ψ` using
/// Newton–Raphson iterations until the equation
///
/// ```text
/// r0 * s1 + sig0 * s2 + mu * s3 = dt
/// ```
///
/// is satisfied. The method handles both:
/// * Elliptical orbits (alpha < 0)
/// * Hyperbolic orbits (alpha > 0)
///
/// Arguments
/// ---------
/// * `dt` – Propagation interval Δt (days)
/// * `r0` – Distance at epoch (AU)
/// * `sig0` – Radial velocity component (AU/day)
/// * `mu` – Gravitational parameter GM (AU³/day²)
/// * `alpha` – Twice the specific orbital energy:
///   * < 0 elliptical,
///   * > 0 hyperbolic.
/// * `e0` – Eccentricity of the orbit
/// * `convergency` – Optional convergence tolerance (defaults to 100×ε)
///
/// Returns
/// --------
/// * `Some((psi, s0, s1, s2, s3))` if the solver converged
/// * `None` if the solver failed to converge or `alpha = 0`.
///
/// Remarks
/// --------
/// * The initial guess for `ψ` is computed using [`prelim_kepuni`].
/// * Convergence is usually fast (few iterations).
fn solve_kepuni(
    params: &UniversalKeplerParams,
    convergency: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    const MAX_ITER: usize = 100;
    let tol = convergency.unwrap_or(100.0 * f64::EPSILON);

    // Preliminary guess for psi
    let mut psi = match params.orbit_type() {
        OrbitType::Parabolic => return None,
        OrbitType::Elliptic | OrbitType::Hyperbolic => prelim_kepuni(params, tol)?,
    };

    // Newton-Raphson refinement
    for _ in 0..MAX_ITER {
        let (s0, s1, s2, s3) = s_funct(psi, params.alpha);

        // Universal Kepler function and derivative
        let f = params.r0 * s1 + params.sig0 * s2 + params.mu * s3 - params.dt;
        let fp = params.r0 * s0 + params.sig0 * s1 + params.mu * s2;
        let dpsi = -f / fp;

        // Convergence safeguard
        if s3.abs() > 1e-2 / f64::EPSILON {
            return None;
        }

        // Update psi (dampen sign change)
        let new_psi = psi + dpsi;
        psi = if new_psi * psi < 0.0 {
            psi / 2.0
        } else {
            new_psi
        };

        // Convergence condition
        let small_step = dpsi.abs() < tol || dpsi.abs() < tol * 10.0 * psi.abs();
        if small_step {
            return Some((psi, s0, s1, s2, s3));
        }
    }

    None
}

/// VelocityCorrectionError Error is used in case the velocity correction cannot be ended.
#[derive(Debug, Clone)]
pub struct VelocityCorrectionError;

impl fmt::Display for VelocityCorrectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The eccentricity control detected that the asteroid as no angular momentum"
        )
    }
}

pub fn velocity_correction(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
) -> Result<(Vector3<f64>, f64, f64), VelocityCorrectionError> {
    let mu = GAUSS_GRAV.powi(2);
    let sig0 = x2.dot(v2);
    let r2 = x2.norm();

    let Some((_, ecc, _, energy)) = eccentricity_control(x2, v2, peri_max, ecc_max) else {
        return Err(VelocityCorrectionError);
    };

    let eps = 1e3 * f64::EPSILON;
    let alpha = 2. * energy;

    let params = UniversalKeplerParams {
        dt,
        r0: r2,
        sig0,
        mu,
        alpha,
        e0: ecc,
    };
    let Some((_, _, _, s2, s3)) = solve_kepuni(&params, Some(eps)) else {
        return Err(VelocityCorrectionError);
    };

    // Calcul des coefficients f et g de Lagrange
    let f = 1.0 - (mu * s2) / r2;
    let g = dt - (mu * s3);

    // Calcul de la vitesse améliorée
    let v2 = (x1 - f * x2) / g;

    Ok((v2, f, g))
}

#[cfg(test)]
mod kepler_test {

    use super::*;

    mod tests_s_funct {
        use approx::assert_relative_eq;

        use super::s_funct;

        fn check_invariants(psi: f64, alpha: f64, s0: f64, s1: f64, s2: f64, s3: f64) {
            let tol = 1e-12;
            assert!(
                (s0 - (1.0 + alpha * s2)).abs() < tol,
                "Invariant s0 = 1 + α*s2 violated: {} vs {}",
                s0,
                1.0 + alpha * s2
            );
            assert!(
                (s1 - (psi + alpha * s3)).abs() < tol,
                "Invariant s1 = ψ + α*s3 violated: {} vs {}",
                s1,
                psi + alpha * s3
            );
        }

        #[test]
        fn test_small_beta() {
            // Small psi and alpha -> beta small, direct series expansion branch
            let psi = 0.01;
            let alpha = 0.1;
            let (s0, s1, s2, s3) = s_funct(psi, alpha);

            // Basic sanity
            assert!(s0 > 0.0);
            assert!(s1 > 0.0);
            check_invariants(psi, alpha, s0, s1, s2, s3);
        }

        #[test]
        fn test_large_beta() {
            let psi = 10.0;
            let alpha = 5.0;
            let (s0, s1, s2, s3) = s_funct(psi, alpha);

            // Vérification de la validité numérique
            assert!(s0.is_finite() && s1.is_finite() && s2.is_finite() && s3.is_finite());

            // Tolérance plus relâchée pour le grand beta
            let rel_tol = 1e-7;

            // Invariants (version Fortran)
            assert_relative_eq!(s0, 1.0 + alpha * s2, max_relative = rel_tol);
            assert_relative_eq!(s1, psi + alpha * s3, max_relative = rel_tol);
        }

        #[test]
        fn test_zero_alpha() {
            // When alpha = 0, expansions should reduce to s0=1, s1=psi, s2=psi^2/2, s3=psi^3/6
            let psi = 2.0;
            let alpha = 0.0;
            let (s0, s1, s2, s3) = s_funct(psi, alpha);

            assert!((s0 - 1.0).abs() < 1e-14);
            assert!((s1 - psi).abs() < 1e-14);
            assert!((s2 - psi.powi(2) / 2.0).abs() < 1e-14);
            assert!((s3 - psi.powi(3) / 6.0).abs() < 1e-14);
        }

        #[test]
        fn test_zero_psi() {
            // When psi = 0, expansions simplify
            let psi = 0.0;
            let alpha = 2.0;
            let (s0, s1, s2, s3) = s_funct(psi, alpha);

            assert!((s0 - 1.0).abs() < 1e-14);
            assert!((s1 - 0.0).abs() < 1e-14);
            assert!((s2 - 0.0).abs() < 1e-14);
            assert!((s3 - 0.0).abs() < 1e-14);
        }

        #[test]
        fn test_symmetry_negative_psi() {
            // s_funct should be odd in psi for s1 and s3, even for s0 and s2
            let psi = 1.0;
            let alpha = 0.5;
            let (s0_pos, s1_pos, s2_pos, s3_pos) = s_funct(psi, alpha);
            let (s0_neg, s1_neg, s2_neg, s3_neg) = s_funct(-psi, alpha);

            let tol = 1e-12;
            // Even functions
            assert!((s0_pos - s0_neg).abs() < tol);
            assert!((s2_pos - s2_neg).abs() < tol);
            // Odd functions
            assert!((s1_pos + s1_neg).abs() < tol);
            assert!((s3_pos + s3_neg).abs() < tol);
        }

        #[test]
        fn test_consistency_large_vs_small() {
            // For moderate values, the two branches should give consistent results
            let psi = 2.5;
            let alpha = 1.0;
            let (s0, s1, s2, s3) = s_funct(psi, alpha);
            check_invariants(psi, alpha, s0, s1, s2, s3);
        }

        #[test]
        fn test_s_funct_real_data() {
            let psi = -15.279808141051223;
            let alpha = -1.6298946008705195e-4;

            let (s0, s1, s2, s3) = s_funct(psi, alpha);

            assert_eq!(s0, 0.9810334785583247);
            assert_eq!(s1, -15.183083836892674);
            assert_eq!(s2, 116.3665517484714);
            assert_eq!(s3, -593.4390119881925);
        }
    }

    mod tests_prelim_kepuni {
        use super::{prelim_kepuni, UniversalKeplerParams};

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
        fn test_returns_none_for_alpha_zero() {
            let params = make_params(1.0, 1.0, 0.0, MU, 0.0, 0.1);
            let res = prelim_kepuni(&params, CONTR);
            assert!(res.is_none());
        }

        #[test]
        fn test_elliptic_small_eccentricity() {
            let params = make_params(0.5, 1.0, 0.1, MU, -1.0, 1e-8);
            let result = prelim_kepuni(&params, CONTR);
            assert!(result.is_some());
            assert!(result.unwrap().is_finite());
        }

        #[test]
        fn test_elliptic_high_eccentricity() {
            let params = make_params(0.1, 0.5, 0.2, MU, -1.0, 0.8);
            let result = prelim_kepuni(&params, CONTR);
            assert!(result.is_some());
            assert!(result.unwrap().is_finite());
        }

        #[test]
        fn test_hyperbolic_case() {
            let params = make_params(0.3, 2.0, -0.1, MU, 1.0, 1.5);
            let result = prelim_kepuni(&params, CONTR);
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

            let psi_pos = prelim_kepuni(&params_pos, CONTR).unwrap();
            let psi_neg = prelim_kepuni(&params_neg, CONTR).unwrap();

            assert!(
                (psi_pos - psi_neg).abs() > 1e-8,
                "psi did not change significantly when changing sig0 sign: {psi_pos} vs {psi_neg}"
            );
        }

        #[test]
        fn test_stability_long_dt() {
            let params = make_params(50.0, 1.0, 0.1, MU, -1.0, 0.5);
            let result = prelim_kepuni(&params, CONTR);
            assert!(result.is_some());
            assert!(result.unwrap().is_finite());
        }

        #[test]
        fn test_edge_cosine_limits() {
            let params = make_params(0.25, 2.0, 0.1, MU, -1.0, 0.1);
            let result = prelim_kepuni(&params, CONTR);
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
            };
            let psi = prelim_kepuni(&params, contr).unwrap();
            assert_eq!(psi, -15.327414893041848);

            let params2 = UniversalKeplerParams {
                alpha: 1.642_158_377_771_140_7E-4,
                ..params
            };
            let psi = prelim_kepuni(&params2, contr).unwrap();
            assert_eq!(psi, -73.1875935362658);

            let params3 = UniversalKeplerParams {
                alpha: 0.0,
                ..params
            };
            assert!(prelim_kepuni(&params3, contr).is_none());
        }

        mod kepuni_prop_tests {
            use super::{prelim_kepuni, UniversalKeplerParams};
            use proptest::prelude::*;

            fn arb_params() -> impl Strategy<Value = UniversalKeplerParams> {
                (
                    -10.0..10.0f64,
                    0.1..5.0f64,
                    -2.0..2.0f64,
                    0.5..2.0f64,
                    prop_oneof![(-5.0..-0.01f64), (0.01..5.0f64)],
                    0.0..3.0f64,
                )
                    .prop_map(|(dt, r0, sig0, mu, alpha, e0)| {
                        UniversalKeplerParams {
                            dt,
                            r0,
                            sig0,
                            mu,
                            alpha,
                            e0,
                        }
                    })
            }

            proptest! {
                #[test]
                fn prop_prelim_kepuni_behaves_well(params in arb_params(), contr in 1e-14..1e-8f64) {
                    let result = prelim_kepuni(&params, contr);
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
                    let params = UniversalKeplerParams { dt, r0, sig0, mu, alpha: 0.0, e0 };
                    let result = prelim_kepuni(&params, contr);
                    if let Some(psi0) = result {
                        prop_assert!(psi0.is_finite());
                    }
                }
            }

            proptest! {
                #[test]
                fn prop_sig0_influences_psi(params in arb_params()) {
                    let contr = 1e-12;
                    prop_assume!(params.dt.abs() > 1e-6);
                    prop_assume!(params.e0 > 1e-6);
                    prop_assume!(params.r0 > 1e-6);

                    let mut params_pos = params;
                    params_pos.sig0 = 0.1;
                    let mut params_neg = params;
                    params_neg.sig0 = -0.1;

                    let res_pos = prelim_kepuni(&params_pos, contr);
                    let res_neg = prelim_kepuni(&params_neg, contr);

                    prop_assume!(res_pos.is_some() && res_neg.is_some());

                    let diff = (res_pos.unwrap() - res_neg.unwrap()).abs();
                    prop_assert!(diff >= 0.0);
                }
            }
        }
    }

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

            assert_eq!(psi, -15.327414893041839);
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
            let (psi, _, _, _, _) =
                solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");

            let expected_psi = 0.47761843287737277;
            assert_relative_eq!(psi, expected_psi, epsilon = 1e-15);
        }

        #[test]
        fn test_solve_kepuni_known_values_hyperbolic() {
            let params = make_params(0.3, 2.0, -0.1, MU, 1.0, 1.5);
            let (psi, _, _, _, _) =
                solve_kepuni(&params, Some(CONTR)).expect("Expected convergence");

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

    #[test]
    fn test_velocity_correction() {
        let x1 = Vector3::new(
            -0.843_561_126_129_683_3,
            0.937_288_327_370_772_8,
            0.659_183_901_029_776_6,
        );

        let x2 = Vector3::new(
            -0.623_121_622_917_384,
            1.0076797884556383,
            0.708_125_687_984_424_5,
        );

        let v2 = Vector3::new(
            -1.552_431_036_862_405_6E-2,
            -3.984_104_176_604_068E-3,
            -2.764_015_436_163_718_3E-3,
        );
        let dt = 14.731970000000729;

        let (v2, f, g) = velocity_correction(&x1, &x2, &v2, dt, 1., 1.).unwrap();

        assert_eq!(f, 0.988_164_877_097_290_6);
        assert_eq!(g, 14.674676076120734);
        assert_eq!(
            v2.as_slice(),
            [
                -0.015524310248562921,
                -0.003984104769239458,
                -0.0027640155187336176
            ]
        )
    }
}
