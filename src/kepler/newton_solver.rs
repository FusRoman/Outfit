#![allow(clippy::doc_overindented_list_items)]
//! Newton–Raphson solver for the universal Kepler equation.
//!
//! This module provides a safeguarded Newton–Raphson solver as the primary
//! solver for the universal Kepler equation. For a more robust but slower
//! alternative, see [`brent_dekker`](crate::kepler::brent_dekker).
//!
//! # Algorithm overview
//!
//! 1. **Initial guess** ([`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni)) — heuristic estimate of  $ \psi_0 $ ,
//!    or a caller-supplied warm start.
//! 2. **Newton–Raphson iteration** — repeated application of
//!     $ \psi_{n+1} = \psi_n - \frac{f(\psi_n)}{f'(\psi_n)} $
//!    where
//!     $ f(\psi) = r_0 \cdot s_1(\psi,\alpha) + \sigma_0 \cdot s_2(\psi,\alpha)
//!               + \mu \cdot s_3(\psi,\alpha) - \Delta t $
//!     $ f'(\psi) = r_0 \cdot s_0(\psi,\alpha) + \sigma_0 \cdot s_1(\psi,\alpha)
//!               + \mu \cdot s_2(\psi,\alpha) $
//! 3. **Safeguards** — derivative guard, step limiter, sign-change damping.
//!
//! # Supported regimes
//!
//! * **Elliptic** ( $ \alpha < 0 $ ) and **Hyperbolic** ( $ \alpha > 0 $ ) motions.
//! * **Parabolic** ( $ \alpha = 0 $ ) is not handled and returns `None`.

use crate::kepler::UniversalKeplerSolution;

use super::params::UniversalKeplerParams;
use super::stumpff::s_funct;

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Maximum number of Newton–Raphson iterations before declaring non-convergence.
const MAX_NEWTON_ITERATIONS: usize = 50;

/// Maximum allowed step magnitude relative to the current  $ |\psi| $ :
///  $ |\Delta\psi| \leq k \cdot (1 + |\psi|) $ .
///
/// Prevents runaway updates on hyperbolic trajectories where  $ \psi $  can grow
/// rapidly.
const MAX_RELATIVE_STEP_FACTOR: f64 = 2.0;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Solve the **universal Kepler equation** using a safeguarded Newton–Raphson
/// iteration, with an optional warm-start guess.
///
/// # Goal
///
/// Find the universal anomaly  $ \psi^* $  such that the residual
///
///  $ f(\psi) = r_0 \cdot s_1(\psi, \alpha) + \sigma_0 \cdot s_2(\psi, \alpha)
///            + \mu \cdot s_3(\psi, \alpha) - \Delta t = 0 $
///
/// vanishes. The universal anomaly  $ \psi $  parametrises the along-track
/// progress of a body regardless of orbit type, unifying the classical
/// anomalies (eccentric, hyperbolic, parabolic) into a single variable.
///
/// # Scientific background
///
/// The **universal-variable formulation** unifies elliptic and hyperbolic
/// motion in a single equation. The Stumpff-like functions  $ (s_0, s_1, s_2,
/// s_3) $  appear in the **Lagrange f–g series**, enabling position and velocity
/// propagation across conic regimes by combining geometry ( $ r_0 $ ,  $ \sigma_0 $ )
/// with dynamics ( $ \mu $ ,  $ \alpha $ ) and the time of flight  $ \Delta t $ .
///
/// The energy parameter  $ \alpha = -\mu / (2a) $  encodes the orbit type:
///
/// | Regime | Condition | Shape |
/// |---|---|---|
/// | Elliptic |  $ \alpha < 0 $  | Closed orbit |
/// | Parabolic |  $ \alpha = 0 $  | Escape at exact escape velocity |
/// | Hyperbolic |  $ \alpha > 0 $  | Open, unbound trajectory |
///
/// # Numerical strategy
///
/// The solver pipeline proceeds in four steps:
///
/// 1. **Tolerance setup** — two tolerances are derived:
///    - Step tolerance:  $ \varepsilon_\psi $  (default  $ 100\varepsilon $ ), controls
///      convergence on  $ |\Delta\psi| $ .
///    - Residual tolerance:  $ \varepsilon_f = 10\varepsilon \cdot (1 + |\Delta t|) $ ,
///      scale-aware to remain meaningful for both small and large flight times.
/// 2. **Orbit-type guard** — parabolic orbits ( $ \alpha = 0 $ ) are rejected
///    immediately.
/// 3. **Initial guess** — uses the caller-supplied `psi_guess` if provided
///    (warm start), otherwise falls back to [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni).
/// 4. **Newton–Raphson loop** with three safeguards applied at each step:
///    - **Derivative guard** — if  $ |f'(\psi)| $  is too small or non-finite,
///       $ \psi $  is halved to escape a flat region.
///    - **Step limiter** —  $ |\Delta\psi| $  is capped at
///       $ k \cdot (1 + |\psi|) $  with  $ k = 2 $ , preventing overshooting on
///      hyperbolic trajectories.
///    - **Sign-change damping** — if the Newton step would cross zero,  $ \psi $
///      is halved instead to avoid oscillations around  $ \psi = 0 $ .
///
/// # Convergence criteria
///
/// The iteration terminates as soon as **any** of the following holds:
///
/// - **Residual**:  $ |f(\psi)| \leq \varepsilon_f $
/// - **Absolute step**:  $ |\Delta\psi| \leq \varepsilon_\psi $
/// - **Relative step**:  $ |\Delta\psi| \leq \varepsilon_\psi \cdot (1 + |\psi|) $
///
/// # Comparison with Brent–Dekker
///
/// | Property | Newton–Raphson | Brent–Dekker |
/// |---|---|---|
/// | Convergence guarantee | No (can diverge) | Yes (within bracket) |
/// | Convergence rate | Quadratic (near root) | Superlinear |
/// | Requires derivative | Yes | No |
/// | Requires bracket | No | Yes |
///
/// Newton–Raphson is preferred when a good initial guess is available and
/// quadratic convergence is desirable. Use the Brent-Dekker method as a fallback when this solver fails to converge.
///
/// # Arguments
///
/// * `params` – Packed orbital parameters for the universal-variable
///   formulation:
///   - `r0` — initial radius  $ r_0 $ ,
///   - `sig0` — radial velocity proxy  $ \sigma_0 = \mathbf{r}_0 \cdot \dot{\mathbf{r}}_0 / \sqrt{\mu} $ ,
///   - `mu` — gravitational parameter  $ \mu $ ,
///   - `alpha` — energy parameter  $ \alpha $ ,
///   - `dt` — time of flight  $ \Delta t $ .
/// * `convergency` – Optional absolute tolerance on  $ |\Delta\psi| $
///   (default:  $ 100\varepsilon \approx 2.2 \times 10^{-14} $ ).
/// * `psi_guess` – Optional warm-start value for  $ \psi_0 $ . If `None`,
///   [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) is used to construct the initial guess.
///
/// # Return
///
/// * `Some(UniversalKeplerSolution)` — the converged universal anomaly
///    $ \psi^* $  together with the evaluated Stumpff functions
///    $ (s_0, s_1, s_2, s_3) $  at  $ \psi^* $ .
/// * `None` if any of the following conditions occur:
///   - the orbit is parabolic ( $ \alpha = 0 $ ),
///   - [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) fails to produce a finite initial guess,
///   - the Newton loop exhausts maximum iteration without meeting any
///     convergence criterion.
///
/// # See also
///
/// * [`solve_kepuni`] – Thin wrapper without warm-start argument.
/// * [`s_funct`](crate::kepler::s_funct) – Stumpff auxiliary functions.
/// * [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) – Heuristic initial guess.
pub fn solve_kepuni_with_guess(params: &UniversalKeplerParams) -> Option<UniversalKeplerSolution> {
    let residual_tolerance = compute_residual_tolerance(params.mu.sqrt() * params.dt);

    let psi_initial = params
        .solver_type
        .params
        .psi_guess
        .map_or_else(|| params.prelim_kepuni(), Some)?;

    run_newton(psi_initial, params, residual_tolerance)
}

/// Solve the universal Kepler equation without a warm-start guess.
///
/// This is a thin wrapper around [`solve_kepuni_with_guess`] that omits the
/// optional warm-start argument. It exists for call sites that only have the
/// orbital parameters and a convergence tolerance available.
///
/// # Behavior
///
/// Delegates directly to [`solve_kepuni_with_guess`] with `psi_guess = None`.
/// The initial guess is then provided by [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni).
///
/// # Arguments
///
/// * `params` – Packed orbital parameters for the universal-variable
///   formulation.
/// * `convergency` – Optional absolute tolerance on  $ |\Delta\psi| $
///   (default:  $ 100\varepsilon $ ).
///
/// # Return
///
/// * `Some(UniversalKeplerSolution)` on convergence.
/// * `None` if convergence fails or the orbit is parabolic ( $ \alpha = 0 $ ).
///
/// # See also
///
/// * [`solve_kepuni_with_guess`] – Extended variant with warm-start support.
pub fn solve_kepuni(params: &UniversalKeplerParams) -> Option<UniversalKeplerSolution> {
    solve_kepuni_with_guess(params)
}

// ---------------------------------------------------------------------------
// Internal pipeline steps
// ---------------------------------------------------------------------------

/// Compute the scale-aware residual tolerance.
///
///  $ \varepsilon_f = 10\varepsilon \cdot (1 + |\sqrt{\mu}\Delta t|) $
///
/// Mixing an absolute floor  $ 10\varepsilon $  with a relative component
///  $ 10\varepsilon \cdot |\sqrt{\mu}\Delta t| $  keeps the criterion
/// meaningful across the full range of flight times. The residual is
/// compared against  $ \sqrt{\mu}\Delta t $  (not raw  $ \Delta t $ ), so the
/// tolerance is scaled accordingly.
#[inline]
fn compute_residual_tolerance(sqrt_mu_dt: f64) -> f64 {
    10.0 * f64::EPSILON * (1.0 + sqrt_mu_dt.abs())
}

/// Evaluate the universal Kepler residual  $ f(\psi) $  and its derivative
///  $ f'(\psi) $  at the given universal anomaly.
///
///  $ f(\psi) = r_0 \cdot s_1 + \sigma_0 \cdot s_2 + s_3 - \sqrt{\mu} \cdot \Delta t $
///  $ f'(\psi) = r_0 \cdot s_0 + \sigma_0 \cdot s_1 + s_2 $
///
/// The derivative chain of the Stumpff functions gives
///  $ \mathrm{d}s_k / \mathrm{d}\psi = s_{k-1} $ . `alpha` here is the
/// reciprocal semi-major-axis convention ( $ \alpha = -1/a $ ), so no extra
/// `mu` factor is needed on the `s2`/`s3` terms — see
/// [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams).
#[inline]
fn kepler_residual_and_derivative(
    psi: f64,
    params: &UniversalKeplerParams,
) -> (f64, f64, (f64, f64, f64, f64)) {
    let stumpff @ (s0, s1, s2, s3) = s_funct(psi, params.alpha);

    let residual = params.r0 * s1 + params.sig0 * s2 + s3 - params.mu.sqrt() * params.dt;
    let derivative = params.r0 * s0 + params.sig0 * s1 + s2;

    (residual, derivative, stumpff)
}

/// Run the core Newton–Raphson iteration loop.
///
/// Returns `Some(UniversalKeplerSolution)` on convergence, or `None` if
/// maximum iteration is exhausted without meeting any convergence
/// criterion.
fn run_newton(
    psi_initial: f64,
    params: &UniversalKeplerParams,
    residual_tolerance: f64,
) -> Option<UniversalKeplerSolution> {
    let mut psi = psi_initial;

    (0..MAX_NEWTON_ITERATIONS).find_map(|_| {
        if !psi.is_finite() {
            psi = 0.5;
            return None;
        }

        let (residual, derivative, stumpff) = kepler_residual_and_derivative(psi, params);

        if residual.abs() <= residual_tolerance {
            return Some(UniversalKeplerSolution::from_raw(psi, stumpff));
        }

        if is_derivative_degenerate(derivative) {
            psi *= 0.5;
            return None;
        }

        let newton_step = compute_safeguarded_step(psi, residual, derivative);
        let psi_candidate = apply_sign_change_damping(psi, newton_step);

        psi = psi_candidate;

        check_step_convergence(
            psi,
            newton_step,
            stumpff,
            params.solver_type.params.convergency,
            params.alpha,
        )
    })
}

// ---------------------------------------------------------------------------
// Newton step helpers
// ---------------------------------------------------------------------------

/// Return `true` if the derivative is too small or non-finite to produce a
/// reliable Newton step.
#[inline]
fn is_derivative_degenerate(derivative: f64) -> bool {
    !derivative.is_finite() || derivative.abs() < 10.0 * f64::EPSILON
}

/// Compute the Newton step  $ \Delta\psi = -f(\psi) / f'(\psi) $ , capped to
/// prevent overshooting.
///
/// The step magnitude is limited to
///  $ |\Delta\psi| \leq k \cdot (1 + |\psi|) $  with  $ k = $  [`MAX_RELATIVE_STEP_FACTOR`].
#[inline]
fn compute_safeguarded_step(psi: f64, residual: f64, derivative: f64) -> f64 {
    let raw_step = -residual / derivative;
    let max_step = MAX_RELATIVE_STEP_FACTOR * (1.0 + psi.abs());
    raw_step.clamp(-max_step, max_step)
}

/// Apply sign-change damping to the Newton candidate.
///
/// If applying `step` to `psi` would flip the sign of $\psi$ (i.e., cross
/// zero), halve $\psi$ instead. This prevents oscillations around $\psi = 0$
/// that are a common instability pattern in the Newton iteration.
#[inline]
fn apply_sign_change_damping(psi: f64, step: f64) -> f64 {
    let psi_candidate = psi + step;
    if psi_candidate * psi < 0.0 {
        0.5 * psi
    } else {
        psi_candidate
    }
}

/// Check convergence on the step size after updating $\psi$.
///
/// Two criteria are tested in order:
///
/// - **Absolute**: $|\Delta\psi| \leq \varepsilon_\psi$
///   — reuses the Stumpff values already computed at the previous $\psi$
///   to avoid one extra call.
/// - **Relative**: $|\Delta\psi| \leq \varepsilon_\psi \cdot (1 + |\psi|)$
///   — recomputes Stumpff functions at the updated $\psi$ to ensure strict
///   consistency of the returned solution.
///
/// Returns `Some(UniversalKeplerSolution)` if either criterion is met,
/// `None` otherwise.
#[inline]
fn check_step_convergence(
    psi: f64,
    step: f64,
    stumpff_prev: (f64, f64, f64, f64),
    step_tolerance: f64,
    alpha: f64,
) -> Option<UniversalKeplerSolution> {
    let step_abs = step.abs();

    // Absolute criterion: reuse already-computed Stumpff values.
    if step_abs <= step_tolerance {
        return Some(UniversalKeplerSolution::from_raw(psi, stumpff_prev));
    }

    // Relative criterion: recompute Stumpff at the updated psi for consistency.
    if step_abs <= step_tolerance * (1.0 + psi.abs()) {
        let stumpff_final = s_funct(psi, alpha);
        return Some(UniversalKeplerSolution::from_raw(psi, stumpff_final));
    }

    None
}
