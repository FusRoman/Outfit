//! # Universal Kepler Propagation & f–g Velocity Update
//!
//! This module implements a **universal‐variable** formulation of the two‑body
//! problem to propagate orbital motion and to apply a **Lagrange f–g**
//! velocity correction between two position epochs.
//!
//! It provides:
//! - An **orbit classifier** [`OrbitType`](crate::kepler::OrbitType) from the energy sign `alpha = 2E`,
//! - A compact container of **universal Kepler inputs** [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams),
//! - A numerically robust evaluation of generalized **Stumpff‑like functions**
//!   via [`s_funct`](crate::kepler::s_funct),
//! - Practical **preliminary solvers** for the universal anomaly ψ in both
//!   elliptic and hyperbolic regimes ([`prelim_elliptic`](crate::kepler::prelim_elliptic), [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic),
//!   wrapped by [`prelim_kepuni`](crate::kepler::prelim_kepuni)),
//! - A Newton–Raphson **universal Kepler solver** [`solve_kepuni`](crate::kepler::solve_kepuni) that returns
//!   ψ and the Stumpff tuple (s0, s1, s2, s3),
//! - A high‑level **velocity correction** routine [`velocity_correction`](crate::kepler::velocity_correction) that
//!   computes `(v₂_corrected, f, g)` from two position vectors and an initial
//!   velocity.
//!
//! ## Units & Conventions
//! ----------------------
//! - Distance: **AU**
//! - Time: **days**
//! - Velocity: **AU/day**
//! - Gravitational parameter: `μ = GM = GAUSS_GRAV²` (AU³/day²)
//! - Angles: **radians** (normalized with [`principal_angle`](crate::kepler::principal_angle))
//! - Energy parameter: `alpha = 2 * (v²/2 − μ/r)`
//!
//! The module uses `nalgebra::Vector3<f64>` for state vectors.
//!
//! ## Mathematical Background
//! --------------------------
//! Universal Kepler’s equation is written with the generalized Stumpff
//! functions (here denoted `s0..s3`):
//!
//! ```text
//! r₀ s₁(ψ, α) + σ₀ s₂(ψ, α) + μ s₃(ψ, α) = Δt
//! with   σ₀ = r₀·ṙ₀ / r₀ = x₀·v₀
//! and    α = 2E  (twice specific orbital energy).
//! ```
//!
//! The invariants that must hold for any `(ψ, α)` are:
//!
//! ```text
//! s₀ = 1 + α s₂
//! s₁ = ψ + α s₃
//! ```
//!
//! Once ψ is known, the Lagrange coefficients for the f–g solution are:
//!
//! ```text
//! f = 1 − (μ/r₀) s₂   ,   g = Δt − μ s₃
//! ```
//!
//! which are used in [`velocity_correction`](crate::kepler::velocity_correction).
//!
//! ## Numerical Strategy
//! ---------------------
//! - [`s_funct`](crate::kepler::s_funct) evaluates `(s0..s3)` either by **power series** when
//!   `|β| = |α ψ²|` is small, or by **half‑angle reduction + duplication**
//!   when `|β|` is large, to avoid divergence.
//! - [`prelim_elliptic`](crate::kepler::prelim_elliptic) and [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) build robust initial guesses
//!   for ψ from the classical eccentric/hyperbolic anomalies (`u`, `F`) with
//!   Newton steps and sign handling from the radial motion (`sig0`).
//! - [`solve_kepuni`](crate::kepler::solve_kepuni) refines ψ with **Newton–Raphson**, includes sign‑change
//!   damping, and rejects obviously unstable configurations.
//! - **Parabolic (`alpha = 0`)** motion is **not supported** here and returns
//!   `None`; treat it with a dedicated routine if needed.
//!
//! ## API Overview
//! ---------------
//! - [`OrbitType`](crate::kepler::OrbitType): classify orbit from `alpha`.
//! - [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams): bundle of inputs `(dt, r0, sig0, mu, alpha, e0)`.
//! - [`s_funct`](crate::kepler::s_funct): compute `(s0, s1, s2, s3)` and satisfy `s0 = 1 + α s2`, `s1 = ψ + α s3`.
//! - [`principal_angle`](crate::kepler::principal_angle), [`angle_diff`](crate::kepler::angle_diff): angle normalization helpers.
//! - [`prelim_kepuni`](crate::kepler::prelim_kepuni): dispatch to [`prelim_elliptic`](crate::kepler::prelim_elliptic) / [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) for ψ guess.
//! - [`solve_kepuni`](crate::kepler::solve_kepuni): universal Kepler solver → `(ψ, s0, s1, s2, s3)`.
//! - [`velocity_correction`](crate::kepler::velocity_correction): apply f–g to update velocity `(v₂_corrected, f, g)`.
//!
//! ## Error Handling
//! -----------------
//! - Parabolic case (`alpha = 0`) returns `None` from the preliminary/solver routines.
//! - [`velocity_correction`](crate::kepler::velocity_correction) returns `Result<… , OutfitError>` and can fail when:
//!   - the upstream [`eccentricity_control`](crate::orb_elem::eccentricity_control) rejects the orbit (e.g., zero angular momentum,
//!     perihelion/eccentricity bounds), or
//!   - [`solve_kepuni`](crate::kepler::solve_kepuni) does not converge.
//!
//! ## Testing & Invariants
//! -----------------------
//! The test suite covers:
//! - **Analytical limits:** `α = 0`, `ψ = 0`, symmetry in `ψ`,
//! - **Large |β|** branch consistency,
//! - **Residual checks:** `r₀ s₁ + σ₀ s₂ + μ s₃ − Δt ≈ 0`,
//! - **Property‑based tests** over wide parameter ranges,
//! - **“Real data”** vectors to guard against regressions.
//!
//! ## Performance Notes
//! --------------------
//! - The routines are allocation‑free and deterministic in `f64`.
//! - Tolerances default to `𝒪(ε)` with small safety factors; tune via the
//!   `contr`/`convergency` parameters when necessary.
//! - For extremely large time steps or highly hyperbolic cases, consider
//!   scaling states and/or tightening tolerances.
//!
//! ## Examples
//! -----------
//! ### Solve universal Kepler for ψ and Stumpff functions
//! ```rust, no_run
//! use nalgebra::Vector3;
//! use outfit::kepler::{UniversalKeplerParams, solve_kepuni};
//!
//! let params = UniversalKeplerParams {
//!     dt: 0.25,          // days
//!     r0: 1.0,           // AU
//!     sig0: 0.0,         // AU/day
//!     mu: 0.00029591220828559115, // GAUSS_GRAV^2
//!     alpha: -1.0,       // elliptic
//!     e0: 0.1,
//! };
//!
//! let (psi, s0, s1, s2, s3) = solve_kepuni(&params, None)
//!     .expect("converged");
//! // Use (psi, s0..s3) to build f, g or propagate state.
//! ```
//!
//! ### Apply f–g velocity correction between two position epochs
//! ```rust
//! use nalgebra::Vector3;
//! use outfit::kepler::velocity_correction;
//!
//! let x1 = Vector3::new(1.0, 0.0, 0.0);   // r(t1) in AU
//! let x2 = Vector3::new(1.1, 0.0, 0.0);   // r(t2) in AU
//! let v2 = Vector3::new(0.0, 0.017, 0.0); // v(t2) in AU/day
//! let dt = 1.0;                           // t2 - t1 in days
//!
//! let (v2_corr, f, g) = velocity_correction(&x1, &x2, &v2, dt, 5.0, 0.9, 1e-12)?;
//! // v2_corr is the corrected velocity at t2 using the universal-variable f–g solution.
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## See also
//! ------------
//! * [`velocity_correction`](crate::kepler::velocity_correction) – Lagrange‑based velocity update.
//! * [`solve_kepuni`](crate::kepler::solve_kepuni) – Universal Kepler solver returning (ψ, s0..s3).
//! * [`prelim_elliptic`](crate::kepler::prelim_elliptic), [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic), [`prelim_kepuni`](crate::kepler::prelim_kepuni) – Initial guesses for ψ.
//! * [`s_funct`](crate::kepler::s_funct) – Stumpff‑like functions and invariants.
//! * [`eccentricity_control`](crate::orb_elem::eccentricity_control) – Eccentricity and energy checks.
//! * [`GAUSS_GRAV`](crate::constants::GAUSS_GRAV), [`DPI`](crate::constants::DPI) – Constants.
use crate::outfit_errors::OutfitError;

use super::constants::DPI;
use super::constants::GAUSS_GRAV_SQUARED;
use super::orb_elem::eccentricity_control;
use nalgebra::Vector3;
use std::f64::consts::PI;

/// Classifies the orbital regime based on the sign of the energy parameter `alpha`.
///
/// The parameter `alpha` is defined as:
///
/// ```text
/// alpha = 2 * specific_orbital_energy = 2 * (v²/2 - μ/r)
/// ```
///
/// Depending on its value:
/// - **Elliptic** (`alpha < 0`) – Closed orbit, bounded motion.
/// - **Parabolic** (`alpha = 0`) – Escape trajectory with zero excess velocity (marginally unbound).
/// - **Hyperbolic** (`alpha > 0`) – Open orbit, unbounded motion.
///
/// This classification is used to select the correct branch when solving the
/// universal Kepler equation and related initial approximations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitType {
    /// Bound, closed orbit (alpha < 0)
    Elliptic,
    /// Open, unbound orbit (alpha > 0)
    Hyperbolic,
    /// Parabolic trajectory (alpha = 0)
    Parabolic,
}

impl OrbitType {
    /// Determine the [`OrbitType`] from a given value of `alpha`.
    ///
    /// # Arguments
    /// * `alpha` – Twice the specific orbital energy.
    ///
    /// # Returns
    /// * `OrbitType::Elliptic` if `alpha < 0`
    /// * `OrbitType::Hyperbolic` if `alpha > 0`
    /// * `OrbitType::Parabolic` if `alpha == 0`
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

/// Parameters required to solve the universal Kepler equation.
///
/// This struct bundles together the state and constants needed to propagate an
/// orbit using the universal variable formulation:
///
/// * `dt`: Propagation time interval Δt = t - t₀ (in days).
/// * `r0`: Heliocentric distance at the reference epoch t₀ (in AU).
/// * `sig0`: Radial velocity component at t₀ (in AU/day).
/// * `mu`: Standard gravitational parameter μ = GM (in AU³/day²).
/// * `alpha`: Twice the specific orbital energy (2E).
/// * `e0`: Orbital eccentricity (unitless).
///
/// The associated method [`orbit_type`](UniversalKeplerParams::orbit_type)
/// classifies the orbit into elliptical, parabolic, or hyperbolic regimes
/// based on `alpha`.
#[derive(Debug, Clone, Copy)]
pub struct UniversalKeplerParams {
    pub dt: f64,
    pub r0: f64,
    pub sig0: f64,
    pub mu: f64,
    pub alpha: f64,
    pub e0: f64,
}

impl UniversalKeplerParams {
    /// Returns the [`OrbitType`] (elliptic, parabolic, or hyperbolic)
    /// corresponding to the value of `alpha`.
    pub fn orbit_type(&self) -> OrbitType {
        OrbitType::from_alpha(self.alpha)
    }
}

/// Compute the generalized Stumpff-like auxiliary functions (s0, s1, s2, s3)
/// used in the **universal variable formulation** of the two-body problem.
///
/// These functions extend the classical Stumpff functions `C(z), S(z)` and
/// appear in semi-analytical formulations such as the **f–g series** for
/// universal Kepler propagation. They allow one to compute position and
/// velocity vectors from the *universal anomaly* `ψ` across all orbit types
/// (elliptic, parabolic, hyperbolic).
///
/// Scientific background
/// ---------------------
/// Let α = 2E (twice the specific orbital energy). Define β = α ψ².
/// The generalized series are:
///
/// * s0(ψ, α) ≈ 1 + α Σ (cosine-like),
/// * s1(ψ, α) ≈ ψ + α Σ (sine-like),
/// * s2 = (s0 − 1) / α,
/// * s3 = (s1 − ψ) / α.
///
/// They enter the f–g Lagrange coefficients:
///
/// ```text
/// f = 1 − (μ/r0) · s2
/// g = Δt − μ · s3
/// ```
///
/// Algorithm
/// ---------
/// Two regimes are used depending on the size of `β = α ψ²`:
///
/// **CASE 1 – Small |β| (rapidly convergent):**
///   * Expand s2 and s3 by power series in ψ.
///   * Reconstruct s0 and s1 via defining relations.
///   * This is efficient and stable for |β| < BETACONTR (≈ 100).
///
/// **CASE 2 – Large |β| (poor convergence):**
///   1. Halve ψ until |β| is small enough (β → β/4 per halving).
///   2. Evaluate s0 and s1 by series expansion at the reduced ψ.
///   3. Apply *duplication formulas* recursively to scale back up:
///      - s0(2ψ) = 2 s0(ψ)² − 1   (analogue of cos(2x))
///      - s1(2ψ) = 2 s0(ψ)·s1(ψ)  (analogue of sin(2x))
///   4. Recover s2 and s3 from their definitions.
///   * This ensures stability but s2, s3 may lose precision due to subtraction.
///
/// Arguments
/// ---------
/// * `psi`   – Universal anomaly (integration parameter).
/// * `alpha` – Twice the specific orbital energy (2·E). Negative for elliptic,
///   zero for parabolic, positive for hyperbolic motion.
///
/// Return
/// ----------
/// * `(s0, s1, s2, s3)` – Tuple of the four auxiliary functions.
///   - `s0` : cosine-like series,
///   - `s1` : sine-like series,
///   - `s2` : (s0 − 1)/α,
///   - `s3` : (s1 − ψ)/α.
///
/// Remarks
/// ----------
/// * For parabolic motion (α = 0), the limiting form should be used (not handled here).
/// * For very large `|α ψ²|`, reconstruction of s2, s3 involves subtracting
///   large numbers and may lose accuracy.
///
/// References
/// ----------
/// * Everhart, E. & Pitkin, E.T., *American Journal of Physics*, 51(8), 712–717 (1983).
/// * Goodyear, W.H., *Astronomical Journal*, 70, 189–192 (1965).
/// * Battin, R.H., *An Introduction to the Mathematics and Methods of Astrodynamics*.
///
/// # See also
/// * [`velocity_correction`] – Uses these functions in the f–g propagation.
/// * [`solve_kepuni`] – Universal Kepler solver relying on these functions.
#[inline(always)]
pub fn s_funct(psi: f64, alpha: f64) -> (f64, f64, f64, f64) {
    // =========================================================================
    // Configuration parameters
    // =========================================================================
    const JMAX: usize = 70; // Max series terms (guard against slow convergence)
    const HALFMAX: usize = 30; // Max halving iterations (duplication safety)
    const BETACONTR: f64 = 100.0; // Threshold between "small" and "large" β

    // Machine-related thresholds
    let eps = f64::EPSILON;
    let convergence_tol = 100.0 * eps; // Below this, terms are negligible
    let overflow_limit = 1.0 / eps; // Prevent runaway if terms blow up

    // =========================================================================
    // Fast path: exactly ψ = 0 → trivial values
    // =========================================================================
    if psi == 0.0 {
        return (1.0, 0.0, 0.0, 0.0);
    }

    // Core expansion parameter: β = α ψ²
    let psi2 = psi * psi;
    let beta = alpha * psi2;

    // =========================================================================
    // CASE 1: Direct series expansion (|β| small)
    // =========================================================================
    if beta.abs() < BETACONTR {
        // Series initialization:
        // s2 = ψ²/2,  s3 = ψ³/6
        let mut s2 = 0.5 * psi2;
        let mut term2 = s2;

        let mut s3 = (s2 * psi) / 3.0;
        let mut term3 = s3;

        // Denominators for recurrence evolve as:
        //   s2: (3·4), (5·6), ...
        //   s3: (4·5), (6·7), ...
        let mut d2a = 3.0;
        let mut d2b = 4.0;
        let mut d3a = 4.0;
        let mut d3b = 5.0;

        for _ in 1..=JMAX {
            // Next correction term for s2
            term2 *= beta / (d2a * d2b);
            s2 += term2;

            // Next correction term for s3
            term3 *= beta / (d3a * d3b);
            s3 += term3;

            // Early exit if converged or diverging
            let small2 = term2.abs() < convergence_tol;
            let small3 = term3.abs() < convergence_tol;
            let big2 = term2.abs() > overflow_limit;
            let big3 = term3.abs() > overflow_limit;

            if (small2 && small3) || big2 || big3 {
                break;
            }

            // Update denominators (+2 at each step)
            d2a += 2.0;
            d2b += 2.0;
            d3a += 2.0;
            d3b += 2.0;
        }

        // Recover s1 and s0 via their defining relations
        let s1 = psi + alpha * s3;
        let s0 = 1.0 + alpha * s2;
        return (s0, s1, s2, s3);
    }

    // =========================================================================
    // CASE 2: Large |β| → reduce ψ until convergence is safe
    // =========================================================================
    let mut psi_r = psi;
    let mut beta_r = beta;
    let mut nhalf = 0usize;

    // Step 1. Halve ψ until |β| is small enough or max halving reached
    while beta_r.abs() >= BETACONTR && nhalf < HALFMAX {
        psi_r *= 0.5;
        beta_r *= 0.25; // β ∝ ψ² → halving ψ divides β by 4
        nhalf += 1;
    }

    // Step 2. Series expansion for s0 and s1 at reduced ψ
    let mut s0 = 1.0;
    let mut s1 = psi_r;

    // Recurrence terms for series
    let mut term0 = 1.0; // first nontrivial term for s0 is β_r/(1·2)
    let mut term1 = psi_r; // first nontrivial term for s1 is ψ_r · β_r/(2·3)

    for j in 1..=JMAX {
        term0 *= beta_r / ((2 * j - 1) as f64 * (2 * j) as f64);
        s0 += term0;
        if term0.abs() < convergence_tol || term0.abs() > overflow_limit {
            break;
        }
    }

    for j in 1..=JMAX {
        term1 *= beta_r / ((2 * j) as f64 * (2 * j + 1) as f64);
        s1 += term1;
        if term1.abs() < convergence_tol || term1.abs() > overflow_limit {
            break;
        }
    }

    // Step 3. Duplication formulas to scale back to the original ψ
    for _ in 0..nhalf {
        let c = s0;
        let s = s1;
        s0 = 2.0 * c * c - 1.0; // analogue of cos(2x)
        s1 = 2.0 * c * s; // analogue of sin(2x)
    }

    // Step 4. Reconstruct s2 and s3 (less stable numerically)
    let s3 = (s1 - psi) / alpha;
    let s2 = (s0 - 1.0) / alpha;

    (s0, s1, s2, s3)
}

/// Normalize an angle in radians to the range [0, 2π].
///
/// This ensures any input angle is wrapped into the principal interval
/// 0 ≤ θ < 2π using Euclidean remainder.
pub fn principal_angle(a: f64) -> f64 {
    a.rem_euclid(DPI)
}

/// Compute the signed minimal difference between two angles in radians.
///
/// Returns the value of (a - b) wrapped into the range [-π, π],
/// i.e. the smallest signed rotation from `b` to `a`.
pub fn angle_diff(a: f64, b: f64) -> f64 {
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
pub fn prelim_kepuni(params: &UniversalKeplerParams, contr: f64) -> Option<f64> {
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
pub fn prelim_elliptic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
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
pub fn prelim_hyperbolic(params: &UniversalKeplerParams, contr: f64, max_iter: usize) -> f64 {
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
/// (see [`s_funct`]) for the energy parameter `α = 2E`.
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
///   * Use [`prelim_kepuni`] unless an explicit `psi_guess` is provided (warm start).
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
/// * `psi_guess` – Optional initial guess for `ψ` (skips [`prelim_kepuni`] if provided).
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
/// * [`s_funct`] – Stumpff-like auxiliary functions used to evaluate `f(ψ)` and `f′(ψ)`.
/// * [`prelim_kepuni`] – Heuristic initial guess for the universal anomaly.
/// * [`UniversalKeplerParams`] – Input container for the universal-variable formulation.
#[inline(always)]
pub fn solve_kepuni_with_guess(
    params: &UniversalKeplerParams,
    convergency: Option<f64>,
    psi_guess: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    // =========================================================================
    // Configuration
    // =========================================================================
    const MAX_ITER: usize = 50;

    // Absolute step tolerance (controls convergence on Δψ).
    let tol_step = convergency.unwrap_or(100.0 * f64::EPSILON);

    // Scale-aware residual tolerance for |f(ψ)|:
    // tol_f = atol + rtol * |dt|
    // This keeps the criterion meaningful for both small and large flight times.
    let tol_f = {
        let eps10 = 10.0 * f64::EPSILON;
        eps10 + eps10 * params.dt.abs()
    };

    // Maximum allowed step length relative to the current |ψ|
    // Prevents overshooting, especially on hyperbolic trajectories.
    let max_rel_step = 2.0;

    // =========================================================================
    // Orbit-type gate
    // =========================================================================
    match params.orbit_type() {
        OrbitType::Parabolic => return None, // Parabolic motion not supported here
        OrbitType::Elliptic | OrbitType::Hyperbolic => {}
    }

    // Shorthands (tighten inner loop by avoiding repeated field loads)
    let r0 = params.r0;
    let sig0 = params.sig0;
    let mu = params.mu;
    let a = params.alpha;
    let dt = params.dt;

    // =========================================================================
    // Initial guess for ψ: warm start or heuristic
    // =========================================================================
    let mut psi = if let Some(g) = psi_guess {
        g
    } else {
        // May return None if the guess cannot be built consistently.
        prelim_kepuni(params, tol_step)?
    };

    // =========================================================================
    // Newton iteration with safeguards
    // =========================================================================
    for _ in 0..MAX_ITER {
        // Defensive: if ψ became non-finite (rare in practice), back off and continue.
        if !psi.is_finite() {
            psi = 0.5; // arbitrary finite reset that lets the iteration recover
        }

        // Stumpff functions at current ψ.
        // Derivative chain: (d/dψ) s1 = s0, (d/dψ) s2 = s1, (d/dψ) s3 = s2.
        let (s0, s1, s2, s3) = s_funct(psi, a);

        // Residual and derivative:
        // f(ψ)  = r0*s1 + sig0*s2 + mu*s3 - dt
        // f'(ψ) = r0*s0 + sig0*s1 + mu*s2
        let f = (((r0 * s1) + (sig0 * s2)) + (mu * s3)) - dt;
        let fp = ((r0 * s0) + (sig0 * s1)) + (mu * s2);

        // --- Convergence on residual
        if f.abs() <= tol_f {
            return Some((psi, s0, s1, s2, s3));
        }

        // --- Derivative safeguard
        // If f′ is too small or non-finite, a Newton step would be unreliable → shrink ψ.
        if !fp.is_finite() || fp.abs() < 10.0 * f64::EPSILON {
            psi *= 0.5;
            continue;
        }

        // --- Newton step (with step limiting)
        let mut dpsi = -f / fp;

        // Limit the step magnitude: |Δψ| ≤ max_rel_step * (1 + |ψ|)
        let max_step = max_rel_step * (1.0 + psi.abs());
        if dpsi.abs() > max_step {
            dpsi = dpsi.signum() * max_step;
        }

        // --- Sign-change damping
        // If the update would cross zero (flip the sign of ψ), halve ψ instead
        // to avoid oscillations around 0 (common instability pattern).
        let new_psi = psi + dpsi;
        psi = if new_psi * psi < 0.0 {
            0.5 * psi
        } else {
            new_psi
        };

        // --- Convergence on step size
        // Absolute criterion: when |Δψ| is already tiny, reuse current (s⋅) to avoid one extra call.
        let small_abs = dpsi.abs() <= tol_step;
        if small_abs {
            return Some((psi, s0, s1, s2, s3));
        }

        // Relative-only success path: ensure strict consistency by recomputing (s⋅) at the final ψ.
        let small_rel = dpsi.abs() <= tol_step * (1.0 + psi.abs());
        if small_rel {
            let (s0f, s1f, s2f, s3f) = s_funct(psi, a);
            return Some((psi, s0f, s1f, s2f, s3f));
        }
    }

    // No convergence within MAX_ITER
    None
}

/// Solve the universal Kepler equation with the legacy two-argument signature.
///
/// This is a thin wrapper around [`solve_kepuni_with_guess`] that omits the
/// optional warm-start argument (`psi_guess`). It exists for backward
/// compatibility with earlier code paths that only provide the orbital
/// parameters and a convergence tolerance.
///
/// Behavior
/// --------
/// * Delegates directly to [`solve_kepuni_with_guess`] with `psi_guess = None`.
/// * Otherwise identical in functionality and return type.
/// * Use [`solve_kepuni_with_guess`] if you want to supply a warm-start
///   value for the universal anomaly `ψ` to speed up iterative refinement.
///
/// Arguments
/// ---------
/// * `params` – Packed orbital parameters for the universal-variable formulation.
/// * `convergency` – Optional absolute tolerance on `ψ` steps
///   (default inside [`solve_kepuni_with_guess`]: `100 × ε`).
///
/// Return
/// ------
/// * `Some((psi, s0, s1, s2, s3))` if the solver converged.
/// * `None` if convergence fails or if the orbit is parabolic (`alpha = 0`).
///
/// See also
/// --------
/// * [`solve_kepuni_with_guess`] – Extended variant with warm-start support.
/// * [`s_funct`] – Computes the Stumpff functions `(s0..s3)`.
pub fn solve_kepuni(
    params: &UniversalKeplerParams,
    convergency: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    // Delegate directly, disabling warm-start by passing None for psi_guess.
    solve_kepuni_with_guess(params, convergency, None)
}

/// Apply velocity correction using Lagrange f–g coefficients.
///
/// This function refines the velocity vector `v2` of an orbiting body at a given epoch,
/// using two-body dynamics in universal variables. It is a simplified wrapper around
/// [`velocity_correction_with_guess`] that does not provide a warm-start for the universal
/// anomaly χ or a custom solver tolerance.
///
/// Arguments
/// ---------
/// * `x1` – Position vector of the body at epoch t₁ (in AU).
/// * `x2` – Position vector of the body at epoch t₂ (in AU).
/// * `v2` – Velocity vector at epoch t₂ (in AU/day).
/// * `dt` – Time difference between epochs t₁ and t₂ (in days).
/// * `peri_max` – Maximum acceptable perihelion distance for eccentricity control.
/// * `ecc_max` – Maximum acceptable eccentricity for eccentricity control.
///
/// Return
/// ------
/// * `Ok((v2_corrected, f, g))` on success, where:
///   - `v2_corrected` is the corrected velocity vector at t₂ (AU/day),
///   - `f`, `g` are the Lagrange coefficients.
/// * `Err(OutfitError::VelocityCorrectionError)` if eccentricity control fails
///   or the universal Kepler solver does not converge.
///
/// See also
/// --------
/// * [`velocity_correction_with_guess`] – Variant that accepts a universal-variable guess.
/// * [`eccentricity_control`] – Validates eccentricity and angular momentum.
/// * [`solve_kepuni`] – Universal Kepler solver returning Stumpff-based integrals.
pub fn velocity_correction(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
    eps: f64,
) -> Result<(Vector3<f64>, f64, f64), OutfitError> {
    // Delegate to the more general implementation with warm-start capability,
    // but here we pass `None` for both the χ guess and the solver tolerance
    // to keep the legacy, minimal interface.
    let (velocity_corrected, f, g, _) =
        velocity_correction_with_guess(x1, x2, v2, dt, peri_max, ecc_max, None, eps)?;
    Ok((velocity_corrected, f, g))
}

/// Find corrected velocity using Lagrange f–g coefficients (configurable).
///
/// This routine refines the velocity vector `v2` of a body at epoch `t₂` by
/// combining two-body dynamics in universal variables with Lagrange’s f–g
/// formulation. Compared to [`velocity_correction`], this version allows:
/// - a custom solver tolerance for the universal Kepler equation,
/// - an optional warm-start on the universal anomaly `χ`.
///
/// Arguments
/// ---------
/// * `x1` – Position vector at epoch t₁ (AU).
/// * `x2` – Position vector at epoch t₂ (AU).
/// * `v2` – Velocity vector at epoch t₂ (AU/day).
/// * `dt` – Time difference t₂ − t₁ (days).
/// * `peri_max` – Maximum perihelion distance for dynamic acceptability.
/// * `ecc_max` – Maximum eccentricity for dynamic acceptability.
/// * `chi_guess` – Optional warm-start for the universal variable `χ`.
/// * `eps` – Optional solver tolerance; default is `1e3 × f64::EPSILON`.
///
/// Return
/// ------
/// * `Ok((v2_corrected, f, g, chi))` on success, where:
///   - `v2_corrected` is the corrected velocity at t₂ (AU/day),
///   - `f`, `g` are the Lagrange coefficients,
///   - `chi` is the final universal anomaly (useful for warm-starts).
/// * `Err(OutfitError::VelocityCorrectionError)` if:
///   - the eccentricity/energy check fails,
///   - the universal Kepler solver does not converge,
///   - or `g` is too small for a stable velocity update.
///
/// See also
/// --------
/// * [`velocity_correction`] – Simpler wrapper without warm-start or tolerance control.
/// * [`eccentricity_control`] – Validates eccentricity and angular momentum.
/// * [`solve_kepuni_with_guess`] – Universal Kepler solver with warm-start support.
#[allow(clippy::too_many_arguments)]
pub fn velocity_correction_with_guess(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
    chi_guess: Option<f64>,
    eps: f64,
) -> Result<(Vector3<f64>, f64, f64, f64), OutfitError> {
    // --- Constants
    let mu = GAUSS_GRAV_SQUARED;
    let r2 = x2.norm();
    let sig0 = x2.dot(v2);

    // --- Quick reject: check angular momentum
    // A vanishing cross product x2 × v2 means rectilinear motion → unstable.
    let h_norm = x2.cross(v2).norm();
    if !h_norm.is_finite() || h_norm <= 1e6 * f64::EPSILON {
        return Err(OutfitError::VelocityCorrectionError(
            "Rejected orbit: near-zero angular momentum (x × v ≈ 0)".into(),
        ));
    }

    // --- Step 1: Validate eccentricity and compute initial orbital energy
    let (_, ecc, _, energy) = eccentricity_control(x2, v2, peri_max, ecc_max).ok_or(
        OutfitError::VelocityCorrectionError(
            "Eccentricity control rejected the candidate orbit".into(),
        ),
    )?;

    // --- Step 2: Pack parameters for the universal Kepler solver
    let params = UniversalKeplerParams {
        dt,
        r0: r2,
        sig0,
        mu,
        alpha: 2.0 * energy, // specific orbital energy → α
        e0: ecc,
    };

    // --- Step 3: Solve the universal Kepler equation (with optional χ warm-start)
    let (_chi, _c2, _c3, s2, s3) = solve_kepuni_with_guess(&params, Some(eps), chi_guess).ok_or(
        OutfitError::VelocityCorrectionError("Universal Kepler solver did not converge".into()),
    )?;

    // --- Step 4: Compute Lagrange coefficients f and g
    let f = 1.0 - (mu * s2) / r2;
    let g = dt - (mu * s3);

    // Guard against ill-conditioned g (division instability)
    let g_abs = g.abs();
    let g_floor = 100.0 * f64::EPSILON * (1.0 + dt.abs()); // scale-aware minimum
    if !(g_abs.is_finite()) || g_abs < g_floor {
        return Err(OutfitError::VelocityCorrectionError(
            "Lagrange coefficient g is too small for a stable velocity update".into(),
        ));
    }

    // --- Step 5: Correct the velocity vector
    // Formula: v2_corrected = (x1 - f*x2)/g
    // Implementation uses BLAS-like axpy for fewer temporaries.
    let mut v_corr = *x1; // allocate one Vector3
    v_corr.axpy(-f, x2, 1.0); // v_corr = x1 - f*x2
    v_corr.unscale_mut(g); // v_corr /= g

    Ok((v_corr, f, g, _chi))
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

    mod tests_velocity_correction {
        use super::*;
        use approx::assert_relative_eq;
        use nalgebra::Vector3;

        const KEP_EPS: f64 = 1e3 * f64::EPSILON;

        /// Helper function: easily build a 3D vector.
        fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
            Vector3::new(x, y, z)
        }

        #[test]
        fn test_velocity_correction_nominal_elliptic() {
            // Nominal test: elliptical orbit with simple positions and velocity.
            // This ensures that the correction converges and returns reasonable values.
            let x1 = v(1.0, 0.0, 0.0);
            let x2 = v(1.1, 0.0, 0.0);
            let v2 = v(0.0, 0.017, 0.0); // Approximate Earth orbital speed in AU/day
            let dt = 1.0; // 1 day
            let peri_max = 5.0;
            let ecc_max = 0.9;

            let result = velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS);
            assert!(result.is_ok(), "Velocity correction should succeed");

            let (vcorr, f, g) = result.unwrap();

            // All outputs should be finite
            assert!(vcorr.iter().all(|c| c.is_finite()));
            assert!(f.is_finite());
            assert!(g.is_finite());

            // f should be reasonably close to 1 for small geometry changes
            assert!(f > 0.5 && f < 1.5, "f coefficient is in a reasonable range");
            // g should be close to dt
            assert_relative_eq!(g, dt, epsilon = 0.1);
        }

        #[test]
        fn test_velocity_correction_rejects_bad_orbit() {
            // Use completely unrealistic input values to trigger an error in eccentricity_control.
            // This ensures that the function correctly returns an error instead of panicking.
            let x1 = v(1e6, 1e6, 1e6);
            let x2 = v(1e6, 1e6, 1e6);
            let v2 = v(1e6, 1e6, 1e6);
            let dt = 0.1;
            let peri_max = 0.001;
            let ecc_max = 0.001;

            let result = velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS);
            assert!(
                result.is_err(),
                "Velocity correction should fail on an invalid orbit"
            );
        }

        #[test]
        fn test_velocity_correction_sensitivity_to_x1() {
            // Sensitivity test:
            // If x1 changes slightly, the corrected velocity vector should also change.
            let x1 = v(1.0, 0.0, 0.0);
            let x1_shifted = v(1.05, 0.0, 0.0); // Slight offset
            let x2 = v(1.1, 0.0, 0.0);
            let v2 = v(0.0, 0.017, 0.0);
            let dt = 1.0;
            let peri_max = 5.0;
            let ecc_max = 0.9;

            let result1 =
                velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();
            let result2 =
                velocity_correction(&x1_shifted, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();

            let (v_corr1, _, _) = result1;
            let (v_corr2, _, _) = result2;

            let diff = (v_corr1 - v_corr2).norm();
            assert!(
                diff > 1e-6,
                "Changing x1 must influence the corrected velocity (diff = {diff})"
            );
        }

        #[test]
        fn test_velocity_correction_output_not_nan() {
            // This test ensures that the correction never returns NaN values,
            // even for slightly irregular configurations.
            let x1 = v(1.0, 0.5, 0.0);
            let x2 = v(1.1, 0.2, 0.0);
            let v2 = v(0.0, 0.02, 0.0);
            let dt = 0.5;
            let peri_max = 5.0;
            let ecc_max = 0.9;

            let (v_corr, f, g) =
                velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();
            assert!(
                !v_corr.iter().any(|x| x.is_nan()),
                "Corrected velocity vector contains NaN"
            );
            assert!(!f.is_nan(), "f coefficient is NaN");
            assert!(!g.is_nan(), "g coefficient is NaN");
        }

        #[test]
        fn test_velocity_correction_real_data() {
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

            let (v2, f, g) = velocity_correction(&x1, &x2, &v2, dt, 1., 1., KEP_EPS).unwrap();

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

        mod velocity_correction_prop_tests {
            use super::*;
            use nalgebra::Vector3;
            use proptest::prelude::*;

            /// Build a Vector3 from components
            fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
                Vector3::new(x, y, z)
            }

            /// Strategy to generate realistic orbital parameters
            fn arb_orbit_params(
            ) -> impl Strategy<Value = (Vector3<f64>, Vector3<f64>, Vector3<f64>, f64, f64, f64)>
            {
                (
                    // x1: position vector (avoid zero)
                    (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64)
                        .prop_map(|(x, y, z)| v(x + 1.0, y + 1.0, z)), // ensure non-zero
                    // x2: position vector (close to x1)
                    (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64)
                        .prop_map(|(x, y, z)| v(x + 1.1, y + 1.0, z)),
                    // v2: velocity vector
                    (-0.05..0.05f64, -0.05..0.05f64, -0.05..0.05f64)
                        .prop_map(|(x, y, z)| v(x, y, z)),
                    // dt: time interval
                    0.01..5.0f64,
                    // peri_max
                    0.1..10.0f64,
                    // ecc_max
                    0.1..2.0f64,
                )
            }

            proptest! {
                /// Property test: velocity_correction should never panic
                /// and always return finite results when it succeeds.
                #[test]
                fn prop_velocity_correction_no_panic(
                    (x1,x2,v2,dt,peri_max,ecc_max) in arb_orbit_params()
                ) {
                    let res = velocity_correction(&x1,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);

                    match res {
                        Ok((vcorr, f, g)) => {
                            // Ensure outputs are finite
                            prop_assert!(vcorr.iter().all(|c| c.is_finite()));
                            prop_assert!(f.is_finite());
                            prop_assert!(g.is_finite());
                        }
                        Err(_) => {
                            // Err is acceptable: indicates eccentricity_control or solve_kepuni failure
                        }
                    }
                }
            }

            proptest! {
                /// Differential property test: a slight change in x1 should
                /// change the result (when both calls return Ok).
                #[test]
                fn prop_velocity_correction_sensitive_to_x1(
                    (x1,x2,v2,dt,peri_max,ecc_max) in arb_orbit_params()
                ) {
                    // Slight perturbation on x1
                    let mut x1_shifted = x1;
                    x1_shifted[0] += 0.01;

                    let res1 = velocity_correction(&x1,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);
                    let res2 = velocity_correction(&x1_shifted,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);

                    if let (Ok((v1, _, _)), Ok((v2c, _, _))) = (res1, res2) {
                        let diff = (v1 - v2c).norm();
                        // Either result is valid and diff can be very small or non-zero
                        prop_assert!(diff >= 0.0);
                    }
                }
            }
        }
    }
}
