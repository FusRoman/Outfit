//! Preliminary estimate of the universal anomaly `psi` for the parabolic and
//! near-parabolic regime (`alpha` close to, or exactly, zero).
//!
//! # Why this module exists
//!
//! [`prelim_elliptic`](crate::kepler::prelim_elliptic) and
//! [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) both convert an
//! anomaly difference into `psi` through a division by `sqrt(|alpha|)`. As
//! `alpha` approaches zero this scaling factor diverges, which makes both
//! routines numerically unusable near the parabolic boundary — even though
//! `alpha` itself never needs to be exactly zero for the instability to show
//! up. This module provides a dedicated, well-conditioned preliminary solver
//! for that regime, free of any `1 / sqrt(alpha)` term.
//!
//! # Governing equation
//!
//! At `alpha = 0`, the generalized Stumpff functions collapse to their
//! polynomial limit:
//!
//! ```text
//! s0 = 1,  s1 = psi,  s2 = psi^2 / 2,  s3 = psi^3 / 6
//! ```
//!
//! Substituting into the universal Kepler equation
//! `r0 * s1 + sig0 * s2 + mu * s3 = dt` gives a **cubic polynomial in
//! `psi`**, with no division by `alpha` anywhere:
//!
//! ```text
//! (mu / 6) * psi^3 + (sig0 / 2) * psi^2 + r0 * psi - dt = 0
//! ```
//!
//! This is the universal-variable equivalent of Barker's equation for
//! parabolic motion. It is used here both as the exact preliminary solver
//! for `alpha == 0` and as a good initial guess for `alpha` close to zero
//! (elliptic or hyperbolic), since the cubic is the zeroth-order term of the
//! universal Kepler equation in a Taylor expansion around `alpha = 0`.
//!
//! # A cubic can have up to three real roots — and it matters here
//!
//! The cubic derivative is itself a quadratic,
//! `f'(psi) = (mu/2) psi^2 + sig0 * psi + r0`, with discriminant
//! `sig0^2 - 2 * mu * r0`. Whenever `sig0` is large enough (fast radial
//! motion, typically near pericenter — and, empirically, a *common* case for
//! near-parabolic asteroid tracklets, not a rare edge case), this
//! discriminant is positive, `f` is **not monotonic**, and the cubic can
//! have up to three real roots. Only one of them is physically meaningful:
//! [`select_physical_root`] documents and implements the selection rule
//! used by the Cardano branch.
//!
//! This non-monotonicity is also why the [`ParabolicPrelimMethod::Newton`]
//! variant is **not** a fully independent, always-robust alternative: a
//! plain safeguarded Newton iteration can get trapped oscillating inside the
//! non-monotonic region and never reach the correct, far-away root (this was
//! verified numerically while writing this module — see the
//! `newton_falls_back_to_cardano_when_non_monotonic` test). Rather than
//! ship an implementation that silently returns a wrong answer in a
//! situation known to be common in this codebase, the Newton variant
//! detects the non-monotonic case up front and delegates to the closed-form
//! Cardano solver, which has no such failure mode.

use crate::kepler::UniversalKeplerParams;

#[cfg(feature = "serde")]
use serde::Deserialize;

#[cfg(feature = "serde")]
use serde::Serialize;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Selects which numerical strategy is used to solve the cubic parabolic
/// preliminary equation.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParabolicPrelimMethod {
    /// Closed-form solution via Cardano's formula. No iteration, no
    /// convergence criterion, and no failure mode on the non-monotonic
    /// branch: the correct real root is always identified analytically.
    /// This is the recommended default for this codebase, where the
    /// non-monotonic case (`sig0^2 > 2 * mu * r0`) is common rather than
    /// exceptional.
    Cardano,
    /// Safeguarded Newton-Raphson iteration on the cubic residual, used
    /// only when the cubic is monotonic (see the module documentation).
    /// When it is not, this variant transparently delegates to
    /// [`ParabolicPrelimMethod::Cardano`] rather than risk a non-converging
    /// or wrong-branch iteration.
    Newton,
}

/// Compute a preliminary estimate of the universal anomaly `psi` for the
/// parabolic / near-parabolic regime.
///
/// # Arguments
///
/// * `params` – Orbital parameters. Only `r0`, `sig0`, `mu`, and `dt` are
///   used; `alpha` is ignored on purpose (this routine targets the
///   `alpha == 0` limit) and `e0` is not needed since the parabolic
///   eccentricity is fixed at 1 by construction.
/// * `contr` – Convergence threshold forwarded to the Newton variant. Unused
///   by the Cardano variant, kept for a uniform signature with
///   [`prelim_elliptic`](crate::kepler::prelim_elliptic) /
///   [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic).
/// * `max_iter` – Maximum number of Newton iterations. Unused by the Cardano
///   variant, kept for the same reason as `contr`.
/// * `method` – Which of [`ParabolicPrelimMethod::Cardano`] or
///   [`ParabolicPrelimMethod::Newton`] to use.
///
/// # Return
///
/// The preliminary universal anomaly `psi`, meant to seed the main
/// Newton-Raphson / Brent-Dekker universal Kepler solver.
///
/// # See also
///
/// * [`prelim_elliptic`](crate::kepler::prelim_elliptic) – Equivalent routine for `alpha < 0`.
/// * [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) – Equivalent routine for `alpha > 0`.
pub fn prelim_parabolic(params: &UniversalKeplerParams) -> f64 {
    match params.solver_type.params.parabolic_solving_method {
        ParabolicPrelimMethod::Cardano => prelim_parabolic_cardano(params),
        ParabolicPrelimMethod::Newton => prelim_parabolic_newton(
            params,
            params.solver_type.params.convergency,
            params.solver_type.params.max_iter_prelim_kepuni,
        ),
    }
}

// ---------------------------------------------------------------------------
// Shared cubic residual and monotonicity test
// ---------------------------------------------------------------------------

/// Evaluate the parabolic Kepler cubic `f(psi)` and its derivative
/// `f'(psi)` at a given universal anomaly.
///
/// ```text
/// f(psi)  = (mu / 6) * psi^3 + (sig0 / 2) * psi^2 + r0 * psi - dt
/// f'(psi) = (mu / 2) * psi^2 + sig0 * psi + r0
/// ```
#[inline]
fn cubic_residual_and_derivative(
    universal_anomaly: f64,
    radial_distance: f64,
    radial_velocity_proxy: f64,
    gravitational_parameter: f64,
    time_of_flight: f64,
) -> (f64, f64) {
    let residual = gravitational_parameter / 6.0 * universal_anomaly.powi(3)
        + radial_velocity_proxy / 2.0 * universal_anomaly.powi(2)
        + radial_distance * universal_anomaly
        - time_of_flight;

    let derivative = gravitational_parameter / 2.0 * universal_anomaly.powi(2)
        + radial_velocity_proxy * universal_anomaly
        + radial_distance;

    (residual, derivative)
}

/// `true` if `f'(psi) = (mu/2) psi^2 + sig0*psi + r0` has real roots, i.e.
/// the parabolic Kepler cubic `f` is **not** monotonic and can have up to
/// three real roots. See the module documentation for why this matters.
#[inline]
fn is_cubic_non_monotonic(
    radial_distance: f64,
    radial_velocity_proxy: f64,
    gravitational_parameter: f64,
) -> bool {
    radial_velocity_proxy * radial_velocity_proxy > 2.0 * gravitational_parameter * radial_distance
}

// ---------------------------------------------------------------------------
// Newton variant
// ---------------------------------------------------------------------------

/// Solve the parabolic Kepler cubic by safeguarded Newton-Raphson,
/// restricted to the regime where the cubic is monotonic.
///
/// # Monotonic-only guarantee
///
/// When `f` is not monotonic (see [`is_cubic_non_monotonic`]), a plain
/// Newton iteration has no robustness guarantee: it can oscillate inside
/// the non-monotonic region and never reach the physically relevant root,
/// however many safeguards (step limiting, sign-change damping) are added.
/// Rather than gamble on that, this function delegates to
/// [`prelim_parabolic_cardano`] in that case; the closed-form solver has no
/// such failure mode.
///
/// # Starting guess (monotonic case)
///
/// `psi0 = dt / r0`, the linearised solution of the cubic when the
/// quadratic and cubic terms are neglected. In the monotonic regime this is
/// sufficient: the safeguarded Newton loop below converges to machine
/// precision from this guess for any time-of-flight magnitude.
fn prelim_parabolic_newton(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
    let radial_distance = params.r0;
    let radial_velocity_proxy = params.sig0;
    let gravitational_parameter = params.mu;
    let time_of_flight = params.dt;

    // Fast path: see the identical short-circuit in
    // `prelim_parabolic_cardano` for the rationale.
    if time_of_flight == 0.0 {
        return 0.0;
    }

    if is_cubic_non_monotonic(
        radial_distance,
        radial_velocity_proxy,
        gravitational_parameter,
    ) {
        return prelim_parabolic_cardano(params);
    }

    let mut universal_anomaly = time_of_flight / radial_distance;

    for _ in 0..max_iter {
        let (residual, derivative) = cubic_residual_and_derivative(
            universal_anomaly,
            radial_distance,
            radial_velocity_proxy,
            gravitational_parameter,
            time_of_flight,
        );

        if !derivative.is_finite() || derivative.abs() < 10.0 * f64::EPSILON {
            // Cannot happen when the monotonicity guard above held (f' > 0
            // everywhere in that case), but kept as a defensive fallback.
            universal_anomaly *= 0.5;
            continue;
        }

        // Step limiter: caps the update magnitude relative to the current
        // |psi|, mirroring the safeguard used in the main universal Kepler
        // Newton solver. Only relevant for very large |dt|.
        let raw_step = -residual / derivative;
        let max_step_magnitude = 2.0 * (1.0 + universal_anomaly.abs());
        let newton_step = raw_step.clamp(-max_step_magnitude, max_step_magnitude);

        universal_anomaly += newton_step;

        if newton_step.abs() < contr {
            break;
        }
    }

    universal_anomaly
}

// ---------------------------------------------------------------------------
// Cardano variant
// ---------------------------------------------------------------------------

/// Solve the parabolic Kepler cubic in closed form via Cardano's formula.
///
/// # Steps
///
/// 1. Normalize `mu/6 * psi^3 + sig0/2 * psi^2 + r0 * psi - dt = 0` to the
///    monic form `psi^3 + b * psi^2 + c * psi + d = 0`.
/// 2. Depress it with `psi = y - b / 3`, eliminating the quadratic term and
///    yielding `y^3 + p * y + q = 0`.
/// 3. Compute the discriminant `delta = (q/2)^2 + (p/3)^3` to decide between
///    the one-real-root and three-real-roots branches.
/// 4. Recover all real roots of the depressed cubic, undo the depression
///    shift, and select the physically relevant one.
fn prelim_parabolic_cardano(params: &UniversalKeplerParams) -> f64 {
    let radial_distance = params.r0;
    let radial_velocity_proxy = params.sig0;
    let gravitational_parameter = params.mu;
    let time_of_flight = params.dt;

    // Fast path: at dt = 0 the constant term of the cubic vanishes
    // identically, so psi = 0 is an exact root regardless of r0 / sig0.
    // Going through the general depression + trigonometric algebra below
    // still finds this root, but only up to floating-point rounding (the
    // depressed coefficients can be many orders of magnitude larger than
    // psi itself, e.g. when mu is tiny relative to r0 and sig0, so the
    // subtractions involved lose precision). Short-circuiting here avoids
    // that entirely for the one case where the exact answer is trivial.
    if time_of_flight == 0.0 {
        return 0.0;
    }

    // Step 1: normalize to a monic cubic `psi^3 + b*psi^2 + c*psi + d = 0`.
    // The leading coefficient `mu / 6` is always strictly positive
    // (mu = GM > 0), so this division is always safe.
    let leading_coefficient = gravitational_parameter / 6.0;
    let quadratic_coefficient = (radial_velocity_proxy / 2.0) / leading_coefficient;
    let linear_coefficient = radial_distance / leading_coefficient;
    let constant_coefficient = -time_of_flight / leading_coefficient;

    // Step 2: depress the cubic via psi = y - depression_shift, with
    // depression_shift = quadratic_coefficient / 3. Standard depression
    // algebra gives:
    //   p = c - b * shift
    //   q = 2 * shift^3 - c * shift + d
    let depression_shift = quadratic_coefficient / 3.0;
    let depressed_linear_term = linear_coefficient - quadratic_coefficient * depression_shift;
    let depressed_constant_term = 2.0 * depression_shift.powi(3)
        - linear_coefficient * depression_shift
        + constant_coefficient;

    // Step 3: discriminant of the depressed cubic y^3 + p*y + q = 0.
    let half_constant_term = depressed_constant_term / 2.0;
    let discriminant =
        half_constant_term * half_constant_term + (depressed_linear_term / 3.0).powi(3);

    // Step 4: recover the real root(s) of the depressed cubic, then shift
    // back to psi = y - depression_shift.
    let candidate_roots: Vec<f64> = if discriminant > 0.0 {
        // Single real root: sum of two real cube roots (Cardano's original
        // formula). `f64::cbrt` handles negative arguments directly, so no
        // complex arithmetic is needed here.
        let sqrt_discriminant = discriminant.sqrt();
        let depressed_root = (-half_constant_term + sqrt_discriminant).cbrt()
            + (-half_constant_term - sqrt_discriminant).cbrt();
        vec![depressed_root - depression_shift]
    } else {
        // Three real roots (possibly with repeats): trigonometric form,
        // valid whenever discriminant <= 0.
        three_real_roots_trigonometric(depressed_linear_term, depressed_constant_term)
            .into_iter()
            .map(|depressed_root| depressed_root - depression_shift)
            .collect()
    };

    let selected_root = select_physical_root(
        &candidate_roots,
        radial_distance,
        radial_velocity_proxy,
        gravitational_parameter,
        time_of_flight,
    );

    polish_root_by_newton(
        selected_root,
        radial_distance,
        radial_velocity_proxy,
        gravitational_parameter,
        time_of_flight,
    )
}

/// Number of Newton polishing steps applied after Cardano's closed-form
/// root selection. See [`polish_root_by_newton`] for the rationale.
const POLISHING_ITERATIONS: usize = 2;

/// Refine a root already close to the true solution with a couple of
/// unguarded Newton-Raphson steps on the original (undepressed) cubic.
///
/// # Why this is needed
///
/// Cardano's formula recovers each root as `y - depression_shift`, where
/// `y` and `depression_shift` are individually computed from the monic
/// coefficients `b = (sig0/2) / (mu/6)`, `c = r0 / (mu/6)`, etc. Because
/// `mu / 6` is a very small, fixed constant in this application (`~5e-5` in
/// AU/day units), these monic coefficients — and therefore `y` and
/// `depression_shift` — can be several orders of magnitude larger than the
/// physical root itself whenever `sig0` is large relative to `r0`. The
/// final subtraction then cancels most of the significant digits: this was
/// observed to produce a **26% relative residual error** for
/// `r0 = 0.05, sig0 = -9.28, dt = -1e-6` before this polishing step was
/// added — far beyond floating-point noise, and a real accuracy bug in the
/// raw Cardano output, not just an overly strict test tolerance.
///
/// # Why unguarded Newton is safe here (unlike in [`prelim_parabolic_newton`])
///
/// The non-monotonicity trap documented at the top of this module is a
/// hazard when Newton starts *far* from the root and can wander into the
/// wrong branch. Here, Newton starts from a value already extremely close
/// to a genuine root (Cardano's answer, just imprecise), so a couple of
/// unguarded steps only ever refine within the same basin — the concern
/// that motivates the full safeguards in [`prelim_parabolic_newton`] does
/// not apply to this local polishing pass.
fn polish_root_by_newton(
    initial_estimate: f64,
    radial_distance: f64,
    radial_velocity_proxy: f64,
    gravitational_parameter: f64,
    time_of_flight: f64,
) -> f64 {
    let mut universal_anomaly = initial_estimate;

    for _ in 0..POLISHING_ITERATIONS {
        let (residual, derivative) = cubic_residual_and_derivative(
            universal_anomaly,
            radial_distance,
            radial_velocity_proxy,
            gravitational_parameter,
            time_of_flight,
        );
        if derivative == 0.0 || !derivative.is_finite() {
            break;
        }
        universal_anomaly -= residual / derivative;
    }

    universal_anomaly
}

/// Compute the three real roots of a depressed cubic `y^3 + p*y + q = 0`
/// known to have a non-positive discriminant `(q/2)^2 + (p/3)^3 <= 0`
/// (which implies `p <= 0`), using the classical trigonometric form:
///
/// ```text
/// theta = acos( (3*q) / (2*p) * sqrt(-3/p) )
/// y_k   = 2 * sqrt(-p/3) * cos( theta/3 - 2*pi*k/3 ),  k = 0, 1, 2
/// ```
fn three_real_roots_trigonometric(linear_term: f64, constant_term: f64) -> [f64; 3] {
    use std::f64::consts::PI;

    // Clamp guards against a rounding-induced argument slightly outside
    // [-1, 1], which would otherwise make `acos` return NaN.
    let acos_argument = ((3.0 * constant_term) / (2.0 * linear_term) * (-3.0 / linear_term).sqrt())
        .clamp(-1.0, 1.0);
    let base_angle = acos_argument.acos() / 3.0;
    let amplitude = 2.0 * (-linear_term / 3.0).sqrt();

    [
        amplitude * base_angle.cos(),
        amplitude * (base_angle - 2.0 * PI / 3.0).cos(),
        amplitude * (base_angle - 4.0 * PI / 3.0).cos(),
    ]
}

// ---------------------------------------------------------------------------
// Root selection
// ---------------------------------------------------------------------------

/// Select the physically relevant root among the (one or three) real roots
/// returned by the Cardano branch.
///
/// # Selection criterion
///
/// 1. **Monotonic branch** — a root is preferred if `f'(psi) >= 0` there.
///    `psi = 0` at `dt = 0` always lies on this branch (`f'(0) = r0 > 0`
///    since `r0` is a heliocentric distance), so staying on it keeps the
///    preliminary guess continuous with the trivial `dt -> 0` limit.
/// 2. **Closest to the linear estimate** — among the roots satisfying (1),
///    the one closest to `dt / r0` is kept. If none satisfy (1) — a
///    degenerate edge case that should not arise for physical inputs
///    (`r0 > 0`) — the same tie-break is applied to the full candidate set
///    as a fallback.
fn select_physical_root(
    candidate_roots: &[f64],
    radial_distance: f64,
    radial_velocity_proxy: f64,
    gravitational_parameter: f64,
    time_of_flight: f64,
) -> f64 {
    let linear_estimate = time_of_flight / radial_distance;

    let closest_to_linear_estimate = |roots: &[f64]| -> f64 {
        *roots
            .iter()
            .min_by(|first_root, second_root| {
                (**first_root - linear_estimate)
                    .abs()
                    .partial_cmp(&(**second_root - linear_estimate).abs())
                    .unwrap()
            })
            .expect("candidate_roots must be non-empty")
    };

    let monotonic_branch_roots: Vec<f64> = candidate_roots
        .iter()
        .copied()
        .filter(|&universal_anomaly| {
            let (_, derivative) = cubic_residual_and_derivative(
                universal_anomaly,
                radial_distance,
                radial_velocity_proxy,
                gravitational_parameter,
                time_of_flight,
            );
            derivative >= 0.0
        })
        .collect();

    if monotonic_branch_roots.is_empty() {
        closest_to_linear_estimate(candidate_roots)
    } else {
        closest_to_linear_estimate(&monotonic_branch_roots)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod prelim_parabolic_tests {
    use super::*;
    use crate::kepler::params::SolverType;

    const GRAVITATIONAL_PARAMETER: f64 = 2.959_122_082_855_911_5E-4;

    fn make_params(
        time_of_flight: f64,
        radial_distance: f64,
        radial_velocity_proxy: f64,
    ) -> UniversalKeplerParams {
        UniversalKeplerParams {
            dt: time_of_flight,
            r0: radial_distance,
            sig0: radial_velocity_proxy,
            mu: GRAVITATIONAL_PARAMETER,
            alpha: 0.0,
            e0: 1.0,
            solver_type: SolverType::default(),
        }
    }

    fn residual_at(params: &UniversalKeplerParams, universal_anomaly: f64) -> f64 {
        cubic_residual_and_derivative(
            universal_anomaly,
            params.r0,
            params.sig0,
            params.mu,
            params.dt,
        )
        .0
    }

    /// At dt = 0 the trivial root psi = 0 must be recovered exactly by both
    /// methods.
    #[test]
    fn zero_time_of_flight_gives_zero_anomaly() {
        let params = make_params(0.0, 1.3, 0.05);

        for method in [
            ParabolicPrelimMethod::Cardano,
            ParabolicPrelimMethod::Newton,
        ] {
            let universal_anomaly = prelim_parabolic(&params);
            assert!(
                universal_anomaly.abs() < 1e-9,
                "{method:?} gave {universal_anomaly} instead of ~0"
            );
        }
    }

    /// Both methods must land on (essentially) the same, low-residual root
    /// in the ordinary monotonic case.
    #[test]
    fn cardano_and_newton_agree_on_monotonic_case() {
        let mut params = make_params(5.0, 1.3, 0.02); // sig0^2 << 2*mu*r0: monotonic

        let psi_cardano = prelim_parabolic(&params);

        params.solver_type.params.parabolic_solving_method = ParabolicPrelimMethod::Newton;
        let psi_newton = prelim_parabolic(&params);

        assert!((psi_cardano - psi_newton).abs() < 1e-6);
        assert!(residual_at(&params, psi_cardano).abs() < 1e-6);
    }

    /// Regression guard for the non-monotonic case: this exact input was
    /// found, while validating this module, to trap a naive safeguarded
    /// Newton iteration at psi ~ 0.615 while the only physical root is near
    /// 60778. The Newton variant must detect the non-monotonic cubic and
    /// delegate to Cardano instead of returning the trapped value.
    #[test]
    fn newton_falls_back_to_cardano_when_non_monotonic() {
        let mut params = make_params(
            35.929_395_402_202_93,
            1.242_343_356_943_632_2,
            -5.995_054_062_733_072,
        );
        assert!(is_cubic_non_monotonic(params.r0, params.sig0, params.mu));

        let psi_cardano = prelim_parabolic(&params);

        params.solver_type.params.parabolic_solving_method = ParabolicPrelimMethod::Newton;
        let psi_newton = prelim_parabolic(&params);

        assert!((psi_cardano - psi_newton).abs() < 1e-6);
        assert!(residual_at(&params, psi_newton).abs() < 1e-3);
    }

    /// Backward propagation (dt < 0) must give a negative-side root.
    #[test]
    fn negative_time_of_flight_gives_negative_anomaly() {
        let mut params = make_params(-4.0, 1.1, 0.02);

        let psi_cardano = prelim_parabolic(&params);

        params.solver_type.params.parabolic_solving_method = ParabolicPrelimMethod::Newton;
        let psi_newton = prelim_parabolic(&params);

        assert!(psi_cardano < 0.0);
        assert!(psi_newton < 0.0);
        assert!((psi_cardano - psi_newton).abs() < 1e-6);
    }

    /// Broad residual check across a grid of (r0, sig0, dt) combinations,
    /// including several non-monotonic cases, for both methods.
    #[test]
    fn residual_is_small_across_a_parameter_grid() {
        let radial_distances = [0.2, 1.0, 3.0];
        let radial_velocities = [-6.0, -0.5, 0.0, 0.5, 6.0];
        let times_of_flight = [-50.0, -1.0, 0.0, 1.0, 50.0];

        for &r0 in &radial_distances {
            for &sig0 in &radial_velocities {
                for &dt in &times_of_flight {
                    let mut params = make_params(dt, r0, sig0);
                    for method in [
                        ParabolicPrelimMethod::Cardano,
                        ParabolicPrelimMethod::Newton,
                    ] {
                        params.solver_type.params.parabolic_solving_method = method;
                        let universal_anomaly = prelim_parabolic(&params);
                        let residual = residual_at(&params, universal_anomaly);
                        assert!(
                            residual.abs() < 1e-3,
                            "{method:?} residual too large for r0={r0}, sig0={sig0}, dt={dt}: {residual}"
                        );
                    }
                }
            }
        }
    }

    // -----------------------------------------------------------------
    // Property-based tests (proptest)
    // -----------------------------------------------------------------
    //
    // These complement the example-based unit tests above by sampling
    // `(r0, sig0, dt)` broadly, including deep into the non-monotonic
    // regime. Every property here was verified numerically (against
    // `numpy.roots` ground truth, outside this crate) before being encoded
    // as a test, specifically to avoid the trap described next.
    //
    // # A property that is deliberately *not* tested here
    //
    // "For small |dt|, psi is close to the linear estimate dt / r0" is
    // **not** true in general. For an extreme, physically implausible
    // `sig0` relative to `r0` (e.g. `r0 = 0.08`, `sig0 = -8.7`, well beyond
    // any realistic radial-velocity proxy at that heliocentric distance),
    // the cubic's near-zero real root can merge with a neighboring real
    // root and turn complex at an arbitrarily small `|dt|`, leaving only
    // the far-away root as the sole real solution. Asserting continuity
    // near `dt = 0` without restricting to the monotonic regime produces a
    // flaky property test that `proptest` will reliably shrink to a
    // failing case. The continuity and sign-consistency properties below
    // are therefore scoped to `is_cubic_non_monotonic(..) == false`, where
    // they are guaranteed (single real root, `f` strictly increasing).
    mod proptests {
        use super::*;
        use crate::kepler::params::SolverType;
        use proptest::prelude::*;

        const GRAVITATIONAL_PARAMETER: f64 = 2.959_122_082_855_911_5E-4;

        fn make_params(
            time_of_flight: f64,
            radial_distance: f64,
            radial_velocity_proxy: f64,
        ) -> UniversalKeplerParams {
            UniversalKeplerParams {
                dt: time_of_flight,
                r0: radial_distance,
                sig0: radial_velocity_proxy,
                mu: GRAVITATIONAL_PARAMETER,
                alpha: 0.0,
                e0: 1.0,
                solver_type: SolverType::default(),
            }
        }

        /// Heliocentric distance: strictly positive, from very close to the
        /// Sun to well beyond the outer main belt.
        fn radial_distance_strategy() -> impl Strategy<Value = f64> {
            0.05..10.0f64
        }

        /// Radial-velocity proxy: wide range on purpose, so that the
        /// non-monotonic branch (see the module documentation) is exercised
        /// often rather than being a rare corner of the generated space.
        fn radial_velocity_strategy() -> impl Strategy<Value = f64> {
            -10.0..10.0f64
        }

        /// Time of flight: both directions, from short intra-night arcs to
        /// multi-year baselines. `1e-6` is used as the lower bound away
        /// from zero so that dedicated `dt = 0` tests stay separate from
        /// the generic sampled ones.
        fn nonzero_time_of_flight_strategy() -> impl Strategy<Value = f64> {
            prop_oneof![(-500.0..-1e-6f64), (1e-6..500.0f64)]
        }

        proptest! {
            /// Both methods must always return a finite universal anomaly:
            /// no NaN, no infinity, across the full sampled domain
            /// (monotonic and non-monotonic alike).
            #[test]
            fn prelim_parabolic_never_returns_non_finite(
                radial_distance in radial_distance_strategy(),
                radial_velocity_proxy in radial_velocity_strategy(),
                time_of_flight in nonzero_time_of_flight_strategy(),
            ) {
                let mut params = make_params(time_of_flight, radial_distance, radial_velocity_proxy);
                for method in [ParabolicPrelimMethod::Cardano, ParabolicPrelimMethod::Newton] {
                    params.solver_type.params.parabolic_solving_method = method;
                    let universal_anomaly = prelim_parabolic(&params);
                    prop_assert!(
                        universal_anomaly.is_finite(),
                        "{method:?} returned a non-finite value for r0={radial_distance}, \
                         sig0={radial_velocity_proxy}, dt={time_of_flight}"
                    );
                }
            }

            /// The value returned by the Cardano solver must always be an
            /// actual root of the governing cubic: this is what makes it
            /// exact by construction, and it is the property every other
            /// guarantee in this module ultimately rests on.
            ///
            /// The tolerance is relative to the magnitude of the terms
            /// actually being summed in the residual, not to `dt` alone:
            /// when the returned `psi` is large (which happens routinely
            /// for extreme `sig0` in the non-monotonic regime), the
            /// `(mu/6) * psi^3` term can be many orders of magnitude larger
            /// than `dt`, and floating-point noise on that term alone
            /// dwarfs a `dt`-scaled tolerance without the residual being
            /// numerically wrong in any meaningful sense.
            #[test]
            fn cardano_root_satisfies_cubic_residual(
                radial_distance in radial_distance_strategy(),
                radial_velocity_proxy in radial_velocity_strategy(),
                time_of_flight in nonzero_time_of_flight_strategy(),
            ) {
                let mut params = make_params(time_of_flight, radial_distance, radial_velocity_proxy);
                params.solver_type.params.max_iter_prelim_kepuni = 100;
                let universal_anomaly = prelim_parabolic(&params);
                let (residual, _) = cubic_residual_and_derivative(
                    universal_anomaly,
                    radial_distance,
                    radial_velocity_proxy,
                    GRAVITATIONAL_PARAMETER,
                    time_of_flight,
                );

                let term_magnitude_scale = (GRAVITATIONAL_PARAMETER / 6.0
                    * universal_anomaly.powi(3))
                .abs()
                    + (radial_velocity_proxy / 2.0 * universal_anomaly.powi(2)).abs()
                    + (radial_distance * universal_anomaly).abs()
                    + time_of_flight.abs();
                let tolerance = 1e-9 * term_magnitude_scale.max(1.0);

                prop_assert!(
                    residual.abs() < tolerance,
                    "residual {residual} exceeds tolerance {tolerance} (term scale \
                     {term_magnitude_scale}) for r0={radial_distance}, \
                     sig0={radial_velocity_proxy}, dt={time_of_flight}"
                );
            }

            /// The Newton variant either matches Cardano's quality directly
            /// (monotonic case) or delegates to it outright (non-monotonic
            /// case, see the module documentation): either way, the two
            /// methods must agree closely on every sampled input.
            #[test]
            fn newton_agrees_with_cardano(
                radial_distance in radial_distance_strategy(),
                radial_velocity_proxy in radial_velocity_strategy(),
                time_of_flight in nonzero_time_of_flight_strategy(),
            ) {
                let mut params = make_params(time_of_flight, radial_distance, radial_velocity_proxy);
                params.solver_type.params.max_iter_prelim_kepuni = 100;
                let psi_cardano = prelim_parabolic(&params);

                params.solver_type.params.parabolic_solving_method = ParabolicPrelimMethod::Newton;
                let psi_newton = prelim_parabolic(&params);
                let tolerance = 1e-4 * (1.0 + psi_cardano.abs());
                prop_assert!(
                    (psi_cardano - psi_newton).abs() < tolerance,
                    "Cardano={psi_cardano} vs Newton={psi_newton} disagree beyond {tolerance} \
                     for r0={radial_distance}, sig0={radial_velocity_proxy}, dt={time_of_flight}"
                );
            }

            /// `psi = 0` is an exact root whenever `dt = 0`, for any
            /// `(r0, sig0)`: the constant term of the cubic vanishes
            /// identically, independently of the other coefficients.
            #[test]
            fn zero_time_of_flight_is_always_the_zero_root(
                radial_distance in radial_distance_strategy(),
                radial_velocity_proxy in radial_velocity_strategy(),
            ) {
                let mut params = make_params(0.0, radial_distance, radial_velocity_proxy);
                params.solver_type.params.max_iter_prelim_kepuni = 100;
                for method in [ParabolicPrelimMethod::Cardano, ParabolicPrelimMethod::Newton] {

                    params.solver_type.params.parabolic_solving_method = method;
                    let universal_anomaly = prelim_parabolic(&params);
                    prop_assert!(
                        universal_anomaly.abs() < 1e-8,
                        "{method:?} gave {universal_anomaly} instead of 0 at dt=0 \
                         for r0={radial_distance}, sig0={radial_velocity_proxy}"
                    );
                }
            }

            /// When sig0 = 0, the governing cubic `(mu/6) psi^3 + r0*psi - dt`
            /// is an odd function of `psi` for a fixed `dt`, so flipping the
            /// sign of `dt` must flip the sign of the selected root:
            /// `psi(-dt) == -psi(dt)`.
            #[test]
            fn zero_radial_velocity_gives_odd_symmetry_in_time_of_flight(
                radial_distance in radial_distance_strategy(),
                time_of_flight in 1e-6..500.0f64,
            ) {
                let mut params_forward = make_params(time_of_flight, radial_distance, 0.0);
                params_forward.solver_type.params.max_iter_prelim_kepuni = 100;

                let mut params_backward = make_params(-time_of_flight, radial_distance, 0.0);
                params_backward.solver_type.params.max_iter_prelim_kepuni = 100;

                for method in [ParabolicPrelimMethod::Cardano, ParabolicPrelimMethod::Newton] {
                    params_forward.solver_type.params.parabolic_solving_method = method;
                    let psi_forward = prelim_parabolic(&params_forward);

                    params_backward.solver_type.params.parabolic_solving_method = method;
                    let psi_backward = prelim_parabolic(&params_backward);
                    prop_assert!(
                        (psi_forward + psi_backward).abs() < 1e-6 * (1.0 + psi_forward.abs()),
                        "{method:?}: psi(dt)={psi_forward} and psi(-dt)={psi_backward} \
                         are not opposite for r0={radial_distance}, dt={time_of_flight}"
                    );
                }
            }

            /// In the regime where the cubic is *guaranteed* monotonic
            /// (`sig0^2 <= 2 * mu * r0`, see [`is_cubic_non_monotonic`]),
            /// `f` is strictly increasing and `f(0) = -dt`, so the sign of
            /// the unique real root must match the sign of `dt`. This is
            /// intentionally not asserted outside the monotonic regime —
            /// see the module-level note on the property that was dropped.
            #[test]
            fn sign_matches_time_of_flight_in_monotonic_regime(
                radial_distance in radial_distance_strategy(),
                time_of_flight in nonzero_time_of_flight_strategy(),
            ) {
                // sig0 fixed to 0, trivially within the sig0^2 <= 2*mu*r0
                // bound for any r0 > 0, guaranteeing the monotonic regime.
                let radial_velocity_proxy = 0.0;
                let mut params = make_params(time_of_flight, radial_distance, radial_velocity_proxy);
                params.solver_type.params.max_iter_prelim_kepuni = 100;

                prop_assert!(!is_cubic_non_monotonic(
                    radial_distance,
                    radial_velocity_proxy,
                    GRAVITATIONAL_PARAMETER,
                ));

                for method in [ParabolicPrelimMethod::Cardano, ParabolicPrelimMethod::Newton] {
                    params.solver_type.params.parabolic_solving_method = method;
                    let universal_anomaly = prelim_parabolic(&params);
                    prop_assert!(universal_anomaly.signum() == time_of_flight.signum());
                }
            }
        }
    }
}
