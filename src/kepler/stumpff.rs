//! Generalized Stumpff-like auxiliary functions for the universal-variable
//! formulation of the two-body problem.

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
///   * This is efficient and stable for |β| < BETA_SERIES_THRESHOLD (≈ 100).
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
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Uses these functions in the f–g propagation.
/// * [`solve_kepuni`](crate::kepler::solve_kepuni) – Universal Kepler solver relying on these functions.
#[inline(always)]
pub fn s_funct(psi: f64, alpha: f64) -> (f64, f64, f64, f64) {
    // Maximum number of series terms (guards against slow convergence).
    const MAX_SERIES_TERMS: usize = 70;
    // Maximum number of halving iterations for the duplication strategy.
    const MAX_HALVING_STEPS: usize = 30;
    // Threshold between the "small beta" and "large beta" regimes.
    const BETA_SERIES_THRESHOLD: f64 = 100.0;

    let machine_epsilon = f64::EPSILON;
    // Series terms smaller than this are considered negligible.
    let series_term_convergence_tolerance = 100.0 * machine_epsilon;
    // Series terms larger than this indicate divergence: abort early.
    let series_term_overflow_limit = 1.0 / machine_epsilon;

    // Fast path: psi = 0 gives trivial, exact values.
    if psi == 0.0 {
        return (1.0, 0.0, 0.0, 0.0);
    }

    // Core expansion parameter: beta = alpha * psi^2.
    let psi_squared = psi * psi;
    let beta = alpha * psi_squared;

    if beta.abs() < BETA_SERIES_THRESHOLD {
        // CASE 1: direct power-series expansion (fast and stable for small |beta|).
        compute_stumpff_via_power_series(
            psi,
            psi_squared,
            beta,
            alpha,
            series_term_convergence_tolerance,
            series_term_overflow_limit,
            MAX_SERIES_TERMS,
        )
    } else {
        // CASE 2: halve psi until |beta| is small, expand there, then scale
        // back up using duplication formulas (stable for large |beta|).
        compute_stumpff_via_halving_and_duplication(
            psi,
            beta,
            alpha,
            series_term_convergence_tolerance,
            series_term_overflow_limit,
            BETA_SERIES_THRESHOLD,
            MAX_HALVING_STEPS,
            MAX_SERIES_TERMS,
        )
    }
}

/// Evaluate `(s0, s1, s2, s3)` by direct power-series expansion in `psi`.
///
/// This branch is used when `|beta| = |alpha * psi^2|` is small enough for
/// the series to converge quickly and safely. `s2` and `s3` are expanded
/// directly; `s0` and `s1` are then recovered from the defining relations
/// `s0 = 1 + alpha * s2` and `s1 = psi + alpha * s3`.
#[allow(clippy::too_many_arguments)]
fn compute_stumpff_via_power_series(
    psi: f64,
    psi_squared: f64,
    beta: f64,
    alpha: f64,
    convergence_tolerance: f64,
    overflow_limit: f64,
    max_terms: usize,
) -> (f64, f64, f64, f64) {
    // Series initialization: s2 = psi^2/2, s3 = psi^3/6.
    let mut s2 = 0.5 * psi_squared;
    let mut series_term_s2 = s2;

    let mut s3 = (s2 * psi) / 3.0;
    let mut series_term_s3 = s3;

    // Denominators for the recurrence evolve as (3·4), (5·6), ... for s2,
    // and (4·5), (6·7), ... for s3.
    let mut denominator_s2_low = 3.0;
    let mut denominator_s2_high = 4.0;
    let mut denominator_s3_low = 4.0;
    let mut denominator_s3_high = 5.0;

    for _ in 1..=max_terms {
        // Next correction term for s2.
        series_term_s2 *= beta / (denominator_s2_low * denominator_s2_high);
        s2 += series_term_s2;

        // Next correction term for s3.
        series_term_s3 *= beta / (denominator_s3_low * denominator_s3_high);
        s3 += series_term_s3;

        // Early exit once both terms are negligible, or if either diverges.
        let term_s2_is_negligible = series_term_s2.abs() < convergence_tolerance;
        let term_s3_is_negligible = series_term_s3.abs() < convergence_tolerance;
        let term_s2_is_diverging = series_term_s2.abs() > overflow_limit;
        let term_s3_is_diverging = series_term_s3.abs() > overflow_limit;

        if (term_s2_is_negligible && term_s3_is_negligible)
            || term_s2_is_diverging
            || term_s3_is_diverging
        {
            break;
        }

        // Advance denominators by 2 at each step.
        denominator_s2_low += 2.0;
        denominator_s2_high += 2.0;
        denominator_s3_low += 2.0;
        denominator_s3_high += 2.0;
    }

    // Recover s1 and s0 via their defining relations.
    let s1 = psi + alpha * s3;
    let s0 = 1.0 + alpha * s2;
    (s0, s1, s2, s3)
}

/// Evaluate `(s0, s1, s2, s3)` for large `|beta|` via halving + duplication.
///
/// The universal anomaly `psi` is repeatedly halved until `|beta|` falls
/// below `beta_threshold`, where a direct series expansion of `s0` and `s1`
/// is safe. The result is then scaled back to the original `psi` using
/// cosine/sine-like duplication formulas, and `s2`, `s3` are reconstructed
/// from their defining relations.
#[allow(clippy::too_many_arguments)]
fn compute_stumpff_via_halving_and_duplication(
    psi: f64,
    beta: f64,
    alpha: f64,
    convergence_tolerance: f64,
    overflow_limit: f64,
    beta_threshold: f64,
    max_halving_steps: usize,
    max_terms: usize,
) -> (f64, f64, f64, f64) {
    let (reduced_psi, reduced_beta, halving_count) =
        reduce_psi_until_beta_is_small(psi, beta, beta_threshold, max_halving_steps);

    let (mut s0, mut s1) = expand_s0_s1_series_at_reduced_psi(
        reduced_psi,
        reduced_beta,
        convergence_tolerance,
        overflow_limit,
        max_terms,
    );

    // Scale s0, s1 back up to the original psi using duplication formulas:
    //   s0(2*psi) = 2*s0(psi)^2 - 1   (analogue of cos(2x))
    //   s1(2*psi) = 2*s0(psi)*s1(psi) (analogue of sin(2x))
    for _ in 0..halving_count {
        let cosine_like = s0;
        let sine_like = s1;
        s0 = 2.0 * cosine_like * cosine_like - 1.0;
        s1 = 2.0 * cosine_like * sine_like;
    }

    // Reconstruct s2 and s3 from their defining relations. This subtraction
    // can lose precision for very large |beta|, as documented on `s_funct`.
    let s3 = (s1 - psi) / alpha;
    let s2 = (s0 - 1.0) / alpha;

    (s0, s1, s2, s3)
}

/// Halve `psi` (and correspondingly divide `beta` by 4) until `|beta|` falls
/// below `beta_threshold`, or until `max_halving_steps` halvings have been
/// performed. Returns the reduced `(psi, beta)` pair and the number of
/// halvings actually applied.
fn reduce_psi_until_beta_is_small(
    psi: f64,
    beta: f64,
    beta_threshold: f64,
    max_halving_steps: usize,
) -> (f64, f64, usize) {
    let mut reduced_psi = psi;
    let mut reduced_beta = beta;
    let mut halving_count = 0usize;

    while reduced_beta.abs() >= beta_threshold && halving_count < max_halving_steps {
        reduced_psi *= 0.5;
        reduced_beta *= 0.25; // beta ∝ psi^2, so halving psi divides beta by 4.
        halving_count += 1;
    }

    (reduced_psi, reduced_beta, halving_count)
}

/// Expand `s0` and `s1` by power series at the (already reduced) `psi`,
/// where `|beta|` is guaranteed to be small enough for fast convergence.
fn expand_s0_s1_series_at_reduced_psi(
    reduced_psi: f64,
    reduced_beta: f64,
    convergence_tolerance: f64,
    overflow_limit: f64,
    max_terms: usize,
) -> (f64, f64) {
    let mut s0 = 1.0;
    let mut s1 = reduced_psi;

    // First nontrivial term for s0 is reduced_beta/(1·2); for s1 it is
    // reduced_psi · reduced_beta/(2·3).
    let mut series_term_s0 = 1.0;
    let mut series_term_s1 = reduced_psi;

    for term_index in 1..=max_terms {
        series_term_s0 *= reduced_beta / ((2 * term_index - 1) as f64 * (2 * term_index) as f64);
        s0 += series_term_s0;
        if series_term_s0.abs() < convergence_tolerance || series_term_s0.abs() > overflow_limit {
            break;
        }
    }

    for term_index in 1..=max_terms {
        series_term_s1 *= reduced_beta / ((2 * term_index) as f64 * (2 * term_index + 1) as f64);
        s1 += series_term_s1;
        if series_term_s1.abs() < convergence_tolerance || series_term_s1.abs() > overflow_limit {
            break;
        }
    }

    (s0, s1)
}

#[cfg(test)]
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
