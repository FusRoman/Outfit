#![allow(clippy::doc_overindented_list_items)]
//! Brent–Dekker solver for the universal Kepler equation.
//!
//! This module provides a bracket-based root-finding approach as an alternative
//! to the Newton–Raphson solver in [`newton_solver`](crate::kepler::newton_solver).
//! The Brent–Dekker method guarantees convergence when a valid bracket is
//! supplied, making it more robust than Newton's method in pathological cases
//! (e.g., near-parabolic orbits, large eccentricities).
//!
//! # Algorithm overview
//!
//! 1. **Bracket construction** ([`bracket_kepler_root`]) — finds an interval
//!     $ [\psi_a, \psi_b] $  such that  $ f(\psi_a) \cdot f(\psi_b) < 0 $ , where
//!
//!     $ f(\psi) = r_0 \cdot s_1(\psi, \alpha) + \sigma_0 \cdot s_2(\psi, \alpha)
//!               + \mu \cdot s_3(\psi, \alpha) - \Delta t $
//!
//! 2. **Root isolation** ([`solve_kepuni_brent_dekker`]) — applies the
//!    Brent–Dekker iteration within that bracket to converge on  $ \psi^* $ .
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

/// Maximum number of iterations for the bracket expansion phase.
const MAX_BRACKET_ITERATIONS: usize = 60;

/// Maximum number of iterations for the Brent–Dekker root-finding phase.
const MAX_BRENT_ITERATIONS: usize = 100;

/// Golden ratio  $ \phi \approx 1.618 $ , used as the geometric expansion factor
/// during bracket construction.
const PHI: f64 = 1.618_033_988_749_895;

// ---------------------------------------------------------------------------
// Residual evaluation
// ---------------------------------------------------------------------------

/// Evaluate the universal Kepler residual  $ f(\psi) $ .
///
///  $ f(\psi) = r_0 \cdot s_1(\psi,\alpha) + \sigma_0 \cdot s_2(\psi,\alpha)
///            + \mu \cdot s_3(\psi,\alpha) - \Delta t $
///
/// A root of this function corresponds to the universal anomaly  $ \psi^* $  at
/// time  $ \Delta t $ .
#[inline(always)]
fn kepler_residual(psi: f64, params: &UniversalKeplerParams) -> f64 {
    let (_, s1, s2, s3) = s_funct(psi, params.alpha);
    params.r0 * s1 + params.sig0 * s2 + params.mu * s3 - params.dt
}

// ---------------------------------------------------------------------------
// Bracket construction
// ---------------------------------------------------------------------------

/// Compute the initial half-width of the bracket search interval.
///
/// When  $ |\psi_0| $  is large enough, the half-width is set proportional to
///  $ |\psi_0| $  to scale the search with the magnitude of the guess. Otherwise
/// an absolute fallback of  $ 1.0 $  is used to avoid a degenerate zero-width
/// interval.
#[inline]
fn initial_bracket_half_width(psi0: f64) -> f64 {
    if psi0.abs() > 1.0 {
        psi0.abs()
    } else {
        1.0
    }
}

/// Expand the bracket endpoint whose residual magnitude is smaller.
///
/// The side with the smaller  $ |f| $  is closer to a root, so expanding it
/// biases the search toward the most promising direction. The expansion
/// grows the interval by a factor  $ \phi $  (golden ratio) at each call.
///
/// Returns the updated `(psi_lo, f_lo, psi_hi, f_hi)` tuple.
#[inline]
fn expand_bracket_toward_root(
    psi_lo: f64,
    f_lo: f64,
    psi_hi: f64,
    f_hi: f64,
    params: &UniversalKeplerParams,
) -> (f64, f64, f64, f64) {
    let interval_width = psi_hi - psi_lo;

    if f_lo.abs() < f_hi.abs() {
        let new_psi_lo = psi_lo - PHI * interval_width;
        (
            new_psi_lo,
            kepler_residual(new_psi_lo, params),
            psi_hi,
            f_hi,
        )
    } else {
        let new_psi_hi = psi_hi + PHI * interval_width;
        (
            psi_lo,
            f_lo,
            new_psi_hi,
            kepler_residual(new_psi_hi, params),
        )
    }
}

/// Return `true` if `(psi_lo, psi_hi)` is a valid bracket, i.e. the residuals
/// at both endpoints have opposite signs (or one is exactly zero).
#[inline]
fn is_valid_bracket(f_lo: f64, f_hi: f64) -> bool {
    f_lo * f_hi <= 0.0
}

/// Find an interval  $ [\psi_a, \psi_b] $  that brackets a root of the universal
/// Kepler residual.
///
/// # Strategy
///
/// Starting from an initial guess  $ \psi_0 $  (from [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) or the
/// caller), the algorithm expands the search interval geometrically until a
/// sign change is detected:
///
/// - The interval is grown by a factor  $ \phi $  (golden ratio,  $ \approx 1.618 $ )
///   at each step, alternating the expansion direction to remain centered on
///   the initial guess.
/// - If the initial guess already brackets the root (one side negative, one
///   positive), the interval is returned immediately.
///
/// # Arguments
///
/// * `psi0` – Initial guess for  $ \psi $  (center of the search).
/// * `params` – Orbital parameters (used to evaluate the residual).
///
/// # Return
///
/// * `Some((psi_lo, psi_hi))` such that  $ f(\psi_{lo}) \cdot f(\psi_{hi}) \leq 0 $ .
/// * `None` if no bracket is found within [`MAX_BRACKET_ITERATIONS`].
fn bracket_kepler_root(psi0: f64, params: &UniversalKeplerParams) -> Option<(f64, f64)> {
    let half_width = initial_bracket_half_width(psi0);

    let init_psi_lo = psi0 - half_width;
    let init_psi_hi = psi0 + half_width;

    let init_state = (
        init_psi_lo,
        kepler_residual(init_psi_lo, params),
        init_psi_hi,
        kepler_residual(init_psi_hi, params),
    );

    (0..MAX_BRACKET_ITERATIONS)
        .try_fold(init_state, |state, _| {
            let (psi_lo, f_lo, psi_hi, f_hi) = state;

            if is_valid_bracket(f_lo, f_hi) {
                // Signal early exit by returning Err (bracket found).
                return Err((psi_lo, psi_hi));
            }

            Ok(expand_bracket_toward_root(
                psi_lo, f_lo, psi_hi, f_hi, params,
            ))
        })
        .err()
}

// ---------------------------------------------------------------------------
// Brent–Dekker state and step logic
// ---------------------------------------------------------------------------

/// Internal state of the Brent–Dekker iteration.
///
/// By convention, `b` is always the *better* approximation, i.e.
///  $ |f(b)| \leq |f(a)| $  is maintained as an invariant throughout the
/// iteration.
struct BrentState {
    /// Current left endpoint of the bracket.
    a: f64,
    /// Residual at `a`.
    f_a: f64,
    /// Current best approximation (right endpoint of the bracket).
    b: f64,
    /// Residual at `b`.
    f_b: f64,
    /// Previous value of `b`, used for IQI and progress guards.
    c: f64,
    /// Residual at `c`.
    f_c: f64,
    /// Step length used in the previous iteration (for the progress guard).
    prev_step_length: f64,
    /// `true` if the previous step was a bisection.
    prev_was_bisection: bool,
}

impl BrentState {
    /// Initialize the Brent–Dekker state from a valid bracket
    /// `(psi_lo, psi_hi)`.
    ///
    /// Ensures from the start that `b` holds the better approximation by
    /// swapping `a` and `b` when  $ |f(a)| < |f(b)| $ .
    fn from_bracket(psi_lo: f64, psi_hi: f64, params: &UniversalKeplerParams) -> Self {
        let (mut a, mut f_a) = (psi_lo, kepler_residual(psi_lo, params));
        let (mut b, mut f_b) = (psi_hi, kepler_residual(psi_hi, params));

        if f_a.abs() < f_b.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut f_a, &mut f_b);
        }

        Self {
            a,
            f_a,
            b,
            f_b,
            c: a,
            f_c: f_a,
            prev_step_length: (psi_hi - psi_lo).abs(),
            prev_was_bisection: true,
        }
    }

    /// Return `true` if the bracket is sufficiently tight or the residual at
    /// `b` is within `tolerance` of zero.
    #[inline]
    fn has_converged(&self, tolerance: f64) -> bool {
        self.f_b.abs() <= tolerance || 0.5 * (self.b - self.a).abs() <= tolerance
    }
}

/// Attempt an inverse quadratic interpolation (IQI) step using the three
/// points  $ (a, f_a) $ ,  $ (b, f_b) $ ,  $ (c, f_c) $ .
///
/// IQI fits a parabola through the three points in the  $ f $ -direction and
/// evaluates it at  $ f = 0 $ . It provides superlinear convergence when the
/// function is smooth near the root.
///
/// Returns `None` if any two of the three  $ f $ -values are too close to
/// distinguish (the denominator would be numerically zero).
#[inline]
fn try_inverse_quadratic_interpolation(
    a: f64,
    f_a: f64,
    b: f64,
    f_b: f64,
    c: f64,
    f_c: f64,
) -> Option<f64> {
    let denominators_are_distinct =
        (f_a - f_c).abs() > f64::EPSILON && (f_b - f_c).abs() > f64::EPSILON;

    denominators_are_distinct.then(|| {
        let term_a = a * f_b * f_c / ((f_a - f_b) * (f_a - f_c));
        let term_b = b * f_a * f_c / ((f_b - f_a) * (f_b - f_c));
        let term_c = c * f_a * f_b / ((f_c - f_a) * (f_c - f_b));
        term_a + term_b + term_c
    })
}

/// Compute a secant step using the two points  $ (a, f_a) $  and  $ (b, f_b) $ .
///
/// The secant method interpolates linearly between the two points and
/// evaluates the interpolant at  $ f = 0 $ .
#[inline]
fn secant_step(a: f64, f_a: f64, b: f64, f_b: f64) -> f64 {
    b - f_b * (b - a) / (f_b - f_a)
}

/// Compute a bisection step between `a` and `b`.
#[inline]
fn bisection_step(a: f64, b: f64) -> f64 {
    0.5 * (a + b)
}

/// Return `true` if the interpolation candidate `s` falls strictly inside the
/// interval  $ (3a/4 + b/4,\; b) $ , i.e. within the inner three-quarters of the
/// bracket  $ [a, b] $ .
///
/// This condition prevents the interpolation step from landing too close to
/// `a`, which would give poor progress.
#[inline]
fn is_candidate_inside_bracket(s: f64, a: f64, b: f64) -> bool {
    let three_quarter_point = (3.0 * a + b) / 4.0;

    if three_quarter_point < b {
        s > three_quarter_point && s < b
    } else {
        s > b && s < three_quarter_point
    }
}

/// Return `true` if the interpolation candidate `s` is making sufficient
/// progress compared to the reference step length.
///
/// The interpolation step is accepted only when  $ |s - b| < \text{reference} / 2 $ .
/// Otherwise the method falls back to bisection to guarantee monotone bracket
/// shrinkage.
#[inline]
fn is_candidate_making_progress(s: f64, b: f64, reference_step_length: f64) -> bool {
    (s - b).abs() < 0.5 * reference_step_length
}

/// Select the next  $ \psi $  candidate and update the bisection flag.
///
/// The interpolation step (IQI or secant) is used only when it:
///
/// 1. Falls strictly inside the bracket (`is_candidate_inside_bracket`), and
/// 2. Makes sufficient progress over the previous step (`is_candidate_making_progress`).
///
/// Otherwise the method falls back to bisection.
///
/// Returns `(next_psi, was_bisection)`.
fn select_next_candidate(state: &BrentState) -> (f64, bool) {
    let interpolated = try_inverse_quadratic_interpolation(
        state.a, state.f_a, state.b, state.f_b, state.c, state.f_c,
    )
    .unwrap_or_else(|| secant_step(state.a, state.f_a, state.b, state.f_b));

    let reference_step_length = if state.prev_was_bisection {
        (state.b - state.c).abs()
    } else {
        state.prev_step_length
    };

    let use_interpolation = is_candidate_inside_bracket(interpolated, state.a, state.b)
        && is_candidate_making_progress(interpolated, state.b, reference_step_length);

    if use_interpolation {
        (interpolated, false)
    } else {
        (bisection_step(state.a, state.b), true)
    }
}

/// Update the bracket endpoints after evaluating the residual at `next_psi`.
///
/// The endpoint whose residual has the same sign as `f_next` is replaced by
/// `next_psi`. Then the invariant  $ |f(b)| \leq |f(a)| $  is restored by
/// swapping `a` and `b` if necessary.
fn update_bracket(state: &mut BrentState, next_psi: f64, f_next: f64) {
    state.c = state.b;
    state.f_c = state.f_b;

    if state.f_a * f_next < 0.0 {
        state.b = next_psi;
        state.f_b = f_next;
    } else {
        state.a = next_psi;
        state.f_a = f_next;
    }

    // Restore invariant: b is always the better approximation.
    if state.f_a.abs() < state.f_b.abs() {
        std::mem::swap(&mut state.a, &mut state.b);
        std::mem::swap(&mut state.f_a, &mut state.f_b);
    }
}

// ---------------------------------------------------------------------------
// Public solver
// ---------------------------------------------------------------------------

/// Solve the **universal Kepler equation** using the Brent–Dekker method.
///
/// # Goal
///
/// Find the universal anomaly $\psi^*$ such that the residual
///
/// $$f(\psi) = r_0 \cdot s_1(\psi, \alpha) + \sigma_0 \cdot s_2(\psi, \alpha)
///            + \mu \cdot s_3(\psi, \alpha) - \Delta t = 0$$
///
/// vanishes. The universal anomaly $\psi$ parametrises the along-track
/// progress of a body regardless of orbit type, unifying the classical
/// anomalies (eccentric, hyperbolic, parabolic) into a single variable.
///
/// # Scientific background
///
/// See [`newton_solver`](crate::kepler::newton_solver) for the physical
/// interpretation of $\psi$, $\alpha$, and the Stumpff functions
/// $(s_0, s_1, s_2, s_3)$.
///
/// The energy parameter $\alpha = -\mu / (2a)$ encodes the orbit type:
///
/// | Regime | Condition | Shape |
/// |---|---|---|
/// | Elliptic | $\alpha < 0$ | Closed orbit |
/// | Parabolic | $\alpha = 0$ | Escape at exact escape velocity |
/// | Hyperbolic | $\alpha > 0$ | Open, unbound trajectory |
///
/// # Numerical strategy
///
/// The solver pipeline proceeds in four steps:
///
/// 1. **Tolerance** — defaults to $100\varepsilon$ if not provided.
/// 2. **Orbit-type guard** — parabolic orbits ($\alpha = 0$) are rejected
///    immediately.
/// 3. **Bracket construction** ([`bracket_kepler_root`]) — starting from the
///    heuristic guess $\psi_0$ returned by [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni), the search
///    interval is expanded geometrically (by factor $\phi$, the golden ratio)
///    until a sign change is detected:
///    $$f(\psi_\text{lo}) \cdot f(\psi_\text{hi}) \leq 0$$
/// 4. **Brent–Dekker iteration** ([`run_brent_dekker`]) — within the bracket,
///    each step selects between:
///    - **Inverse quadratic interpolation (IQI)** — superlinear convergence
///      when three distinct points are available and the step falls inside the
///      bracket.
///    - **Secant step** — used when IQI is not applicable.
///    - **Bisection** — fallback when interpolation fails to make sufficient
///      progress, guaranteeing that the bracket width decreases monotonically.
///
/// # Comparison with Newton–Raphson
///
/// | Property | Newton–Raphson | Brent–Dekker |
/// |---|---|---|
/// | Convergence guarantee | No (can diverge) | Yes (within bracket) |
/// | Convergence rate | Quadratic (near root) | Superlinear |
/// | Requires derivative | Yes | No |
/// | Requires bracket | No | Yes |
///
/// Brent–Dekker is preferred when robustness matters more than raw speed, or
/// as a fallback when Newton's method fails to converge (e.g., near-parabolic
/// orbits, large eccentricities, or ill-conditioned initial conditions).
///
/// # Arguments
///
/// * `params` – Packed orbital parameters for the universal-variable
///   formulation:
///   - `r0` — initial radius $r_0$,
///   - `sig0` — radial velocity proxy $\sigma_0 = \mathbf{r}_0 \cdot \dot{\mathbf{r}}_0 / \sqrt{\mu}$,
///   - `mu` — gravitational parameter $\mu$,
///   - `alpha` — energy parameter $\alpha$,
///   - `dt` — time of flight $\Delta t$.
/// * `convergency` – Optional absolute tolerance on $|f(\psi^*)|$ and
///   $|\psi_b - \psi_a|$ (default: $100\varepsilon \approx 2.2 \times 10^{-14}$).
///
/// # Return
///
/// * `Some(UniversalKeplerSolution)` — the converged universal anomaly
///   $\psi^*$ together with the evaluated Stumpff functions
///   $(s_0, s_1, s_2, s_3)$ at $\psi^*$.
/// * `None` if any of the following conditions occur:
///   - the orbit is parabolic ($\alpha = 0$),
///   - [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) fails to produce a finite initial guess,
///   - no sign change is found within [`MAX_BRACKET_ITERATIONS`] expansions,
///   - the Brent–Dekker loop exhausts [`MAX_BRENT_ITERATIONS`] without
///     meeting `tolerance`.
///
/// # See also
///
/// * [`bracket_kepler_root`] – Bracket construction helper.
/// * [`run_brent_dekker`] – Core iteration loop.
/// * [`s_funct`](crate::kepler::s_funct) – Stumpff auxiliary functions.
/// * [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) – Initial guess used to
///   seed the bracket search.
/// * [`solve_kepuni_newton`](crate::kepler::newton_solver::solve_kepuni_newton) –
///   Alternative Newton–Raphson solver.
pub fn solve_kepuni_brent_dekker(
    params: &UniversalKeplerParams,
) -> Option<UniversalKeplerSolution> {
    let psi_initial_guess = params
        .solver_type
        .params
        .psi_guess
        .map_or_else(|| params.prelim_kepuni(), Some)?;

    let (psi_lo, psi_hi) = bracket_kepler_root(psi_initial_guess, params)?;

    run_brent_dekker(psi_lo, psi_hi, params)
}

// ---------------------------------------------------------------------------
// Internal pipeline steps
// ---------------------------------------------------------------------------

/// Run the core Brent–Dekker iteration loop starting from the bracket
/// `(psi_lo, psi_hi)`.
///
/// Returns `Some(UniversalKeplerSolution)` on convergence, or `None` if
/// [`MAX_BRENT_ITERATIONS`] is exhausted without meeting `tolerance`.
fn run_brent_dekker(
    psi_lo: f64,
    psi_hi: f64,
    params: &UniversalKeplerParams,
) -> Option<UniversalKeplerSolution> {
    let mut state = BrentState::from_bracket(psi_lo, psi_hi, params);

    (0..MAX_BRENT_ITERATIONS).find_map(|_| {
        if state.has_converged(params.solver_type.params.convergency) {
            return Some(build_solution(state.b, params));
        }

        let (next_psi, was_bisection) = select_next_candidate(&state);
        let f_next = kepler_residual(next_psi, params);

        state.prev_step_length = (state.b - state.c).abs();
        state.prev_was_bisection = was_bisection;
        update_bracket(&mut state, next_psi, f_next);

        None
    })
}

/// Build a [`UniversalKeplerSolution`] from the converged universal anomaly
/// `psi_star`.
///
/// Recomputes the Stumpff functions at $\psi^*$ to populate the solution
/// fields.
#[inline]
fn build_solution(psi_star: f64, params: &UniversalKeplerParams) -> UniversalKeplerSolution {
    let stumpff = s_funct(psi_star, params.alpha);
    UniversalKeplerSolution::from_raw(psi_star, stumpff)
}
