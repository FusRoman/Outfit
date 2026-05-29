//! Assembly and resolution of the weighted least-squares normal equations for
//! differential orbit correction.
//!
//! This module implements step 5 of the two-body differential-correction
//! pipeline.  Given a set of astrometric observations with their partial
//! derivatives, it builds the normal-equation system, inverts it, and returns
//! the orbital-element correction together with uncertainty matrices and the
//! normalised RMS residual.
//!
//! ## Pipeline
//!
//! ```text
//! [ObsAndElementPartials]  →  ObservationEquation  →  solve_weighted_least_squares()
//!   (∂α/∂elem, ∂δ/∂elem)      + residuals                     ↓
//!   per observation            + weights          element_correction = Γ · GᵀWξ
//!                                                 OrbitalUncertainty { Γ, GᵀWG }
//!                                                 normalised_rms = √(ξᵀWξ / num_measurements)
//! ```
//!
//! ## Normal equations
//!
//! The design matrix **G** has shape `(2m × 6)`: for observation `i`, row `2i`
//! carries `∂α/∂elem` and row `2i+1` carries `∂δ/∂elem`.
//!
//! The per-observation weight sub-matrix is
//!
//! ```text
//! W_i = [ weight_ra     weight_cross ]
//!        [ weight_cross  weight_dec   ]
//! ```
//!
//! The normal matrix and right-hand side accumulate as
//!
//! ```text
//! GᵀWG[j,k] = Σ_i  g_α[i,j]·weight_ra ·g_α[i,k]  +  g_δ[i,j]·weight_dec·g_δ[i,k]
//!                +  weight_cross·(g_δ[i,j]·g_α[i,k]  +  g_α[i,j]·g_δ[i,k])
//!
//! GᵀWξ[j]   = Σ_i  (g_α[i,j]·weight_ra  + g_δ[i,j]·weight_cross)·ξ_α
//!                +  (g_α[i,j]·weight_cross + g_δ[i,j]·weight_dec )·ξ_δ
//! ```
//!
//! The element correction is then \\( \delta x = \Gamma \cdot G^\top W \xi \\)
//! where \\( \Gamma = (G^\top W G)^{-1} \\), and the normalised RMS is
//!
//! \\[
//!   \text{normalised\_rms} = \sqrt{\frac{\xi^\top W \xi}{n_{\text{measurements}}}}
//! \\]

use nalgebra::{Cholesky, Matrix6, Vector6, QR};

use crate::outfit_errors::OutfitError;

// ─────────────────────────────────────────────────────────────────────────────
// Public types
// ─────────────────────────────────────────────────────────────────────────────

/// Covariance and normal matrices produced by the differential-correction step.
#[derive(Debug, Clone)]
pub struct OrbitalUncertainty {
    /// Normal matrix \\( G^\top W G \\) (6×6).
    pub normal_matrix: Matrix6<f64>,
    /// Covariance matrix \\( \Gamma = (G^\top W G)^{-1} \\) (6×6).
    pub covariance: Matrix6<f64>,
    /// `true` if matrix inversion succeeded, `false` otherwise.
    pub inversion_succeeded: bool,
}

/// Per-observation input for assembling the normal equations.
///
/// Groups the partial derivatives (from
/// `ObsAndElementPartials`), the astrometric
/// residuals, and the statistical weights for a single optical observation.
#[derive(Debug, Clone)]
pub struct ObservationEquation {
    /// \\( \partial\alpha / \partial\text{elem} \\) — vector of 6 partial
    /// derivatives of the right ascension with respect to the equinoctial
    /// elements `(a, h, k, p, q, λ)`.
    pub d_ra_d_elem: Vector6<f64>,
    /// \\( \partial\delta / \partial\text{elem} \\) — vector of 6 partial
    /// derivatives of the declination.
    pub d_dec_d_elem: Vector6<f64>,
    /// Right-ascension residual \\( \xi_\alpha = \alpha_{\text{obs}} - \alpha_{\text{calc}} \\),
    /// wrapped to \\( (-\pi, \pi] \\).
    ///
    /// Use [`angular_diff`] to compute this residual.
    pub residual_ra: f64,
    /// Declination residual \\( \xi_\delta = \delta_{\text{obs}} - \delta_{\text{calc}} \\).
    pub residual_dec: f64,
    /// Right-ascension weight \\( w_\alpha = 1/\sigma_\alpha^2 \\).
    pub weight_ra: f64,
    /// Declination weight \\( w_\delta = 1/\sigma_\delta^2 \\).
    pub weight_dec: f64,
    /// Cross-correlation weight term \\( w_{\alpha\delta} \\).
    ///
    /// Set to `0.0` for uncorrelated optical observations (the common case).
    pub weight_cross: f64,
    /// Whether this observation contributes to the fit.
    ///
    /// When `false` the entry is skipped; its contribution to the normal
    /// matrix, right-hand side, and normalised RMS is exactly zero.
    pub active: bool,
}

impl ObservationEquation {
    /// Creates an uncorrelated entry (`weight_cross = 0`), the standard case
    /// for optical astrometry with independent RA/Dec uncertainties.
    ///
    /// # Arguments
    ///
    /// - `d_ra_d_elem` — partial derivatives \\( \partial\alpha/\partial\text{elem} \\)
    ///   from `ObsAndElementPartials`.
    /// - `d_dec_d_elem` — partial derivatives \\( \partial\delta/\partial\text{elem} \\).
    /// - `residual_ra` — angular difference \\( \alpha_{\text{obs}} - \alpha_{\text{calc}} \\);
    ///   compute with [`angular_diff`].
    /// - `residual_dec` — difference \\( \delta_{\text{obs}} - \delta_{\text{calc}} \\).
    /// - `sigma_ra` — right-ascension uncertainty in radians.
    /// - `sigma_dec` — declination uncertainty in radians.
    /// - `active` — pass `false` to exclude this observation from the fit.
    #[allow(clippy::too_many_arguments)]
    pub fn uncorrelated(
        d_ra_d_elem: Vector6<f64>,
        d_dec_d_elem: Vector6<f64>,
        residual_ra: f64,
        residual_dec: f64,
        sigma_ra: f64,
        sigma_dec: f64,
        active: bool,
    ) -> Self {
        Self {
            d_ra_d_elem,
            d_dec_d_elem,
            residual_ra,
            residual_dec,
            weight_ra: 1.0 / sigma_ra.powi(2),
            weight_dec: 1.0 / sigma_dec.powi(2),
            weight_cross: 0.0,
            active,
        }
    }
}

/// Output of [`solve_weighted_least_squares`]: the orbital-element correction,
/// uncertainty matrices, normalised RMS, and observation count.
#[derive(Debug, Clone)]
pub struct DifferentialCorrectionResult {
    /// Orbital-element correction \\( \delta x = \Gamma \cdot G^\top W \xi \\).
    ///
    /// Components for which `free_elements[j] == false` are set to `0`.
    pub element_correction: Vector6<f64>,
    /// Covariance and normal matrices from the inversion step.
    pub uncertainty: OrbitalUncertainty,
    /// Dimensionless normalised RMS
    /// \\( \sqrt{\xi^\top W \xi \,/\, \text{num\_measurements}} \\).
    pub normalised_rms: f64,
    /// Number of scalar measurements used (2 per active optical observation).
    pub num_measurements: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// Utilities
// ─────────────────────────────────────────────────────────────────────────────

/// Computes the angular difference `a − b` wrapped to \\( (-\pi, \pi] \\).
///
/// # Arguments
///
/// - `a` — minuend angle, in radians (any value).
/// - `b` — subtrahend angle, in radians (any value).
///
/// # Returns
///
/// The difference `a − b` reduced to the half-open interval \\( (-\pi, \pi] \\).
///
/// # Examples
///
/// ```
/// use outfit::differential_orbit_correction::least_square::angular_diff;
/// use std::f64::consts::PI;
///
/// // Ordinary difference
/// let d = angular_diff(0.1, 0.3);
/// assert!((d - (-0.2)).abs() < 1e-15);
///
/// // Wrapping near 2π
/// let d = angular_diff(0.1, 2.0 * PI - 0.1);
/// assert!((d - 0.2).abs() < 1e-14);
/// ```
pub fn angular_diff(a: f64, b: f64) -> f64 {
    use std::f64::consts::{PI, TAU};
    let mut d = a - b;
    // Reduce to (−π, π] with at most a few iterations (typical for nearby orbits)
    while d > PI {
        d -= TAU;
    }
    while d < -PI {
        d += TAU;
    }
    d
}

// ─────────────────────────────────────────────────────────────────────────────
// Normal-equation assembly and resolution
// ─────────────────────────────────────────────────────────────────────────────

/// Assembles the normal matrix \\( G^\top W G \\), solves the linear system,
/// and returns the differential correction \\( \delta x \\).
///
/// Active entries in `equations` contribute to the normal matrix, the
/// right-hand side, and the weighted sum of squared residuals.  Fixed elements
/// (those with `free_elements[j] == false`) have their row and column zeroed
/// in the normal matrix and are forced to `element_correction[j] = 0`.
///
/// # Arguments
///
/// - `equations` — one [`ObservationEquation`] per observation (order is irrelevant).
/// - `free_elements` — six-element boolean mask: `free_elements[j] = true`
///   means element `j` is free (solved for), `false` means it is held fixed.
///
/// # Returns
///
/// A [`DifferentialCorrectionResult`] containing the correction vector,
/// the [`OrbitalUncertainty`] matrices, the normalised RMS, and the
/// measurement count.
///
pub fn solve_weighted_least_squares(
    equations: &[ObservationEquation],
    free_elements: &[bool; 6],
) -> Result<DifferentialCorrectionResult, OutfitError> {
    // ─── 1. Count scalar measurements ───────────────────────────────────────
    let num_measurements: usize = equations.iter().filter(|e| e.active).count() * 2;

    // ─── 2. Accumulate normal matrix (GᵀWG) and right-hand side (GᵀWξ) ─────
    let mut normal_mat = Matrix6::<f64>::zeros();
    let mut right_hand_side = Vector6::<f64>::zeros();
    let mut weighted_sq_sum = 0.0_f64;

    for eq in equations.iter().filter(|e| e.active) {
        let partials_ra = &eq.d_ra_d_elem;
        let partials_dec = &eq.d_dec_d_elem;
        let weight_ra = eq.weight_ra;
        let weight_dec = eq.weight_dec;
        let weight_cross = eq.weight_cross;
        let residual_ra = eq.residual_ra;
        let residual_dec = eq.residual_dec;

        // normal_mat[j,k] += partials_ra[j]·weight_ra ·partials_ra[k]
        //                   + partials_dec[j]·weight_dec·partials_dec[k]
        //                   + weight_cross·(partials_dec[j]·partials_ra[k]
        //                                 + partials_ra[j]·partials_dec[k])
        for j in 0..6 {
            for k in 0..6 {
                normal_mat[(j, k)] += partials_ra[j] * weight_ra * partials_ra[k]
                    + partials_dec[j] * weight_dec * partials_dec[k]
                    + weight_cross
                        * (partials_dec[j] * partials_ra[k] + partials_ra[j] * partials_dec[k]);
            }

            // right_hand_side[j] += (partials_ra[j]·weight_ra  + partials_dec[j]·weight_cross)·ξ_α
            //                     + (partials_ra[j]·weight_cross + partials_dec[j]·weight_dec )·ξ_δ
            right_hand_side[j] += (partials_ra[j] * weight_ra + partials_dec[j] * weight_cross)
                * residual_ra
                + (partials_ra[j] * weight_cross + partials_dec[j] * weight_dec) * residual_dec;
        }

        // Q += weight_ra·ξ_α² + weight_dec·ξ_δ² + 2·weight_cross·ξ_α·ξ_δ
        weighted_sq_sum += weight_ra * residual_ra * residual_ra
            + weight_dec * residual_dec * residual_dec
            + 2.0 * weight_cross * residual_ra * residual_dec;
    }

    // ─── 3. Apply the free_elements mask ────────────────────────────────────
    //
    // Fixed elements (free_elements[j] == false) have their row and column
    // zeroed in the normal matrix, with the diagonal entry set to 1 to keep
    // the matrix invertible.  Their right-hand-side entry is also zeroed so
    // they contribute nothing to the correction.
    for j in 0..6 {
        if !free_elements[j] {
            for k in 0..6 {
                normal_mat[(j, k)] = 0.0;
                normal_mat[(k, j)] = 0.0;
            }
            normal_mat[(j, j)] = 1.0; // prevents singularity
            right_hand_side[j] = 0.0;
        }
    }

    // ─── 4. Inversion: Cholesky first, QR as fallback ───────────────────────
    let (covariance_mat, inversion_succeeded) = invert_normal_matrix(normal_mat);

    // ─── 5. Correction δx = Γ · GᵀWξ ────────────────────────────────────────
    //    Fixed components (free_elements == false) are forced to zero.
    let mut element_correction = if inversion_succeeded {
        covariance_mat * right_hand_side
    } else {
        Vector6::zeros()
    };

    for j in 0..6 {
        if !free_elements[j] {
            element_correction[j] = 0.0;
        }
    }

    // ─── 6. Normalised RMS ───────────────────────────────────────────────────
    let normalised_rms = if num_measurements > 0 {
        (weighted_sq_sum / num_measurements as f64).sqrt()
    } else {
        0.0
    };

    Ok(DifferentialCorrectionResult {
        element_correction,
        uncertainty: OrbitalUncertainty {
            normal_matrix: normal_mat,
            covariance: covariance_mat,
            inversion_succeeded,
        },
        normalised_rms,
        num_measurements,
    })
}

/// Attempts to invert the 6×6 normal matrix using Cholesky decomposition,
/// falling back to QR decomposition if the matrix is not positive-definite.
///
/// Returns `(covariance, inversion_succeeded)` where `covariance` is the
/// inverse (or a zero matrix on failure).
fn invert_normal_matrix(m: Matrix6<f64>) -> (Matrix6<f64>, bool) {
    // Attempt 1: Cholesky — preferred as it exploits symmetry
    if let Some(chol) = Cholesky::new(m) {
        return (chol.inverse(), true);
    }

    // Attempt 2: QR — more robust, handles ill-conditioned cases
    let qr = QR::new(m);
    match qr.try_inverse() {
        Some(inv) => (inv, true),
        None => (Matrix6::zeros(), false),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Covariance rescaling
// ─────────────────────────────────────────────────────────────────────────────

/// Applies a posterior rescaling factor \\( \mu \\) to the covariance and
/// normal matrices in `uncertainty`.
///
/// The rescaling factor is computed as follows:
///
/// - If \\( n_{\text{free}} \ge n_{\text{measurements}} \\): \\( \mu = 1 \\)
///   (under-determined system — no rescaling applied).
/// - If \\( \text{normalised\_rms} > 1 \\) (poor fit):
///   \\( \mu = \text{normalised\_rms} \cdot \sqrt{n_{\text{measurements}} / (n_{\text{measurements}} - n_{\text{free}})} \\).
/// - Otherwise (acceptable fit):
///   \\( \mu = \sqrt{n_{\text{measurements}} / (n_{\text{measurements}} - n_{\text{free}})} \\).
///
/// After rescaling, `covariance` is multiplied by \\( \mu^2 \\) and
/// `normal_matrix` is divided by \\( \mu^2 \\).
///
/// # Arguments
///
/// - `uncertainty` — the uncertainty matrices to rescale, modified in place.
/// - `num_free_params` — number of free parameters (number of `true` entries
///   in `free_elements`, typically 6).
/// - `num_measurements` — number of scalar measurements used (2 × number of
///   active optical observations).
/// - `normalised_rms` — dimensionless normalised RMS from
///   [`solve_weighted_least_squares`].
pub fn rescale_covariance(
    uncertainty: &mut OrbitalUncertainty,
    num_free_params: usize,
    num_measurements: usize,
    normalised_rms: f64,
) {
    let mu = if num_free_params < num_measurements {
        let factor =
            ((num_measurements as f64) / (num_measurements - num_free_params) as f64).sqrt();
        if normalised_rms > 1.0 {
            normalised_rms * factor
        } else {
            factor
        }
    } else {
        1.0
    };
    let mu2 = mu * mu;
    uncertainty.covariance *= mu2;
    uncertainty.normal_matrix /= mu2;
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod least_square_tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    const ALL_FREE: [bool; 6] = [true; 6];

    // ── Helpers ───────────────────────────────────────────────────────────────

    /// Generates `n` observations using an identity-block partial derivative
    /// pattern.
    ///
    /// The system requires at least 3 observations (`num_measurements = 6`) to
    /// be determined.  Observation `k % 6` sets `partials_ra[k] = 1` and
    /// `partials_dec[(k+1) % 6] = 1`, which produces a diagonal normal matrix.
    fn make_identity_equations(n: usize, sigma: f64) -> Vec<ObservationEquation> {
        (0..n)
            .map(|i| {
                let k = i % 6;
                let mut partials_ra = Vector6::zeros();
                let mut partials_dec = Vector6::zeros();
                partials_ra[k] = 1.0;
                partials_dec[(k + 1) % 6] = 1.0;
                ObservationEquation::uncorrelated(
                    partials_ra,
                    partials_dec,
                    0.0,
                    0.0,
                    sigma,
                    sigma,
                    true,
                )
            })
            .collect()
    }

    // ── Unit tests ─────────────────────────────────────────────────────────────

    /// Zero residuals must yield a zero correction and `normalised_rms = 0`.
    #[test]
    fn test_zero_residuals_give_zero_correction() {
        let equations = make_identity_equations(6, 1e-5);
        let result = solve_weighted_least_squares(&equations, &ALL_FREE).unwrap();
        assert_abs_diff_eq!(result.normalised_rms, 0.0, epsilon = 1e-15);
        for j in 0..6 {
            assert_abs_diff_eq!(result.element_correction[j], 0.0, epsilon = 1e-15);
        }
        assert!(result.uncertainty.inversion_succeeded);
        assert_eq!(result.num_measurements, 12);
    }

    /// Verify that \\( \Gamma \cdot G^\top W G \approx I \\).
    #[test]
    fn test_covariance_times_normal_is_identity() {
        let equations = make_identity_equations(6, 1e-5);
        let result = solve_weighted_least_squares(&equations, &ALL_FREE).unwrap();
        let product = result.uncertainty.covariance * result.uncertainty.normal_matrix;
        let diff = product - Matrix6::identity();
        assert!(
            diff.norm() < 1e-10,
            "Γ·GᵀWG must be close to I, error={}",
            diff.norm()
        );
    }

    /// Inactive observations (`active = false`) must not affect the result.
    #[test]
    fn test_rejected_observations_have_no_contribution() {
        let sigma = 1e-5;
        // One inactive entry with a non-zero residual
        let rejected = ObservationEquation::uncorrelated(
            Vector6::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            Vector6::new(1.0, 1.0, 1.0, 1.0, 1.0, 1.0),
            1e-4,
            1e-4,
            sigma,
            sigma,
            false, // inactive
        );
        // Fill with enough active zero-residual observations
        let mut equations = make_identity_equations(6, sigma);
        equations.push(rejected);
        let result = solve_weighted_least_squares(&equations, &ALL_FREE).unwrap();
        // The zero residuals on the active observations dominate → correction ≈ 0
        assert_abs_diff_eq!(result.normalised_rms, 0.0, epsilon = 1e-15);
        assert_eq!(result.num_measurements, 12); // 6 active obs × 2
    }

    /// Setting `free_elements[0] = false` must force `element_correction[0] = 0`.
    #[test]
    fn test_fixed_element_has_zero_correction() {
        let sigma = 1e-5;
        // Non-zero residuals to produce a non-trivial correction
        let equations: Vec<ObservationEquation> = (0..6)
            .map(|i| {
                let k = i % 6;
                let mut partials_ra = Vector6::zeros();
                let mut partials_dec = Vector6::zeros();
                partials_ra[k] = 1.0;
                partials_dec[(k + 1) % 6] = 1.0;
                ObservationEquation::uncorrelated(
                    partials_ra,
                    partials_dec,
                    1e-4,
                    1e-4,
                    sigma,
                    sigma,
                    true,
                )
            })
            .collect();

        let mut free_elements = [true; 6];
        free_elements[0] = false; // fix the first element

        let result = solve_weighted_least_squares(&equations, &free_elements).unwrap();
        assert_abs_diff_eq!(result.element_correction[0], 0.0, epsilon = 1e-15);
        // The other free elements may have non-zero corrections
    }

    /// The correction must equal the residual magnitude for a well-conditioned
    /// diagonal system.
    ///
    /// For 6 observations with `partials_ra = e_j` and
    /// `partials_dec = e_{(j+1) mod 6}`, uniform residuals `ξ_α = ξ_δ = r`,
    /// and equal weights `weight_ra = weight_dec = 1/σ²`, the normal matrix is
    /// diagonal with `2/σ²` on the diagonal, so
    /// `element_correction[j] = r` for all `j`.
    #[test]
    fn test_correction_magnitude_matches_residual() {
        let sigma = 1e-5;
        let r = 1e-5_f64;
        let equations: Vec<ObservationEquation> = (0..6)
            .map(|i| {
                let k = i % 6;
                let mut partials_ra = Vector6::zeros();
                let mut partials_dec = Vector6::zeros();
                partials_ra[k] = 1.0;
                partials_dec[(k + 1) % 6] = 1.0;
                ObservationEquation::uncorrelated(
                    partials_ra,
                    partials_dec,
                    r,
                    r,
                    sigma,
                    sigma,
                    true,
                )
            })
            .collect();

        let result = solve_weighted_least_squares(&equations, &ALL_FREE).unwrap();

        // GᵀWξ[j] = r/σ² × (count of j in partials_ra or partials_dec) = 2·r/σ²
        // GᵀWG[j,j] = 2/σ²  →  element_correction[j] = (2r/σ²) / (2/σ²) = r
        for j in 0..6 {
            assert_abs_diff_eq!(result.element_correction[j], r, epsilon = 1e-12);
        }
    }

    /// `num_measurements` must count only active observations (2 scalars each).
    #[test]
    fn test_num_measurements_counts_active_only() {
        let sigma = 1e-5;
        let mut equations = make_identity_equations(6, sigma);
        // Add 3 inactive observations
        for _ in 0..3 {
            let mut e = equations[0].clone();
            e.active = false;
            equations.push(e);
        }
        let result = solve_weighted_least_squares(&equations, &ALL_FREE).unwrap();
        assert_eq!(result.num_measurements, 12); // 6 active × 2
    }

    // ── Tests for angular_diff ────────────────────────────────────────────────

    #[test]
    fn test_angular_diff_basic() {
        use std::f64::consts::{PI, TAU};
        // Simple case
        assert_abs_diff_eq!(angular_diff(0.5, 0.3), 0.2, epsilon = 1e-15);
        // Positive wrap: 0.1 − (2π − 0.1) = 0.2 − 2π → +0.2 after correction
        assert_abs_diff_eq!(angular_diff(0.1, TAU - 0.1), 0.2, epsilon = 1e-14);
        // Negative wrap
        assert_abs_diff_eq!(angular_diff(TAU - 0.1, 0.1), -0.2, epsilon = 1e-14);
        // Exactly π must remain in (−π, π]
        let d = angular_diff(PI, 0.0);
        assert!(d > -PI && d <= PI);
    }

    // ── Tests for rescale_covariance ──────────────────────────────────────────

    #[test]
    fn test_rescale_identity_when_underdetermined() {
        let mut uncertainty = OrbitalUncertainty {
            normal_matrix: Matrix6::identity(),
            covariance: Matrix6::identity(),
            inversion_succeeded: true,
        };
        // num_free_params == num_measurements → mu = 1, no change
        rescale_covariance(&mut uncertainty, 6, 6, 1.5);
        assert_abs_diff_eq!(
            (uncertainty.covariance - Matrix6::identity()).norm(),
            0.0,
            epsilon = 1e-15
        );
    }

    #[test]
    fn test_rescale_good_fit() {
        // normalised_rms = 0.5 ≤ 1 → mu = sqrt(12/6) = sqrt(2)
        let mut uncertainty = OrbitalUncertainty {
            normal_matrix: Matrix6::identity(),
            covariance: Matrix6::identity(),
            inversion_succeeded: true,
        };
        rescale_covariance(&mut uncertainty, 6, 12, 0.5);
        let mu2 = 2.0_f64; // sqrt(12/6)² = 2
        assert_abs_diff_eq!(
            (uncertainty.covariance - Matrix6::identity() * mu2).norm(),
            0.0,
            epsilon = 1e-14
        );
        assert_abs_diff_eq!(
            (uncertainty.normal_matrix - Matrix6::identity() / mu2).norm(),
            0.0,
            epsilon = 1e-14
        );
    }

    #[test]
    fn test_rescale_bad_fit() {
        // normalised_rms = 2.0 > 1 → mu = 2.0 * sqrt(12/6) = 2.0 * sqrt(2)
        let mut uncertainty = OrbitalUncertainty {
            normal_matrix: Matrix6::identity(),
            covariance: Matrix6::identity(),
            inversion_succeeded: true,
        };
        rescale_covariance(&mut uncertainty, 6, 12, 2.0);
        let mu = 2.0 * 2.0_f64.sqrt();
        let mu2 = mu * mu;
        assert_abs_diff_eq!(
            (uncertainty.covariance - Matrix6::identity() * mu2).norm(),
            0.0,
            epsilon = 1e-14
        );
    }

    /// Oracle test:
    ///
    /// The dataset is 4 synthetic optical observations with diagonal weights (correl = 0).
    ///
    /// Key oracle values (tolerance 1e-10):
    ///   csinor = 2.0758197845135249   (= weighted_sq_sum in Rust)
    ///   dx0    = [1.6756756756756757e-5, 2.3513513513513514e-5,
    ///            -2.2972972972972989e-5, 3.5405405405405403e-5,
    ///             3.2702702702702714e-5,-4.5945945945945951e-5]
    #[test]
    fn test_oracle_min_sol() {
        use nalgebra::Vector6;

        // ── observations ────────────────────────────────────────────────────
        // Matching the Fortran driver exactly: g rows are standard basis
        // vectors (or mixed for obs 4), weights = 1/sigma^2.

        let obs1 = ObservationEquation::uncorrelated(
            Vector6::new(1.0, 0.0, 0.0, 0.0, 0.0, 0.0), // g_ra  = e_1
            Vector6::new(0.0, 1.0, 0.0, 0.0, 0.0, 0.0), // g_dec = e_2
            2.0e-5,                                     // residual_ra
            3.0e-5,                                     // residual_dec
            1.0e-5,                                     // sigma_ra
            1.0e-5,                                     // sigma_dec
            true,
        );
        let obs2 = ObservationEquation::uncorrelated(
            Vector6::new(0.0, 0.0, 1.0, 0.0, 0.0, 0.0), // g_ra  = e_3
            Vector6::new(0.0, 0.0, 0.0, 1.0, 0.0, 0.0), // g_dec = e_4
            -1.0e-5,                                    // residual_ra
            5.0e-5,                                     // residual_dec
            2.0e-5,                                     // sigma_ra
            1.5e-5,                                     // sigma_dec
            true,
        );
        let obs3 = ObservationEquation::uncorrelated(
            Vector6::new(0.0, 0.0, 0.0, 0.0, 1.0, 0.0), // g_ra  = e_5
            Vector6::new(0.0, 0.0, 0.0, 0.0, 0.0, 1.0), // g_dec = e_6
            4.0e-5,                                     // residual_ra
            -2.0e-5,                                    // residual_dec
            1.5e-5,                                     // sigma_ra
            2.0e-5,                                     // sigma_dec
            true,
        );
        let obs4 = ObservationEquation::uncorrelated(
            Vector6::new(0.5, 0.5, 0.5, 0.5, 0.5, 0.5), // g_ra  (mixed)
            Vector6::new(0.5, -0.5, 0.5, -0.5, 0.5, -0.5), // g_dec (mixed)
            1.0e-5,                                     // residual_ra
            1.0e-5,                                     // residual_dec
            1.0e-5,                                     // sigma_ra
            1.0e-5,                                     // sigma_dec
            true,
        );

        let equations = [obs1, obs2, obs3, obs4];
        let free = [true; 6];

        let result = solve_weighted_least_squares(&equations, &free)
            .expect("solve_weighted_least_squares failed");

        // ── Oracle constants ───────────────────────────
        let oracle_dx = [
            1.6756756756756757e-5_f64,
            2.3513513513513514e-5,
            -2.297_297_297_297_299e-5,
            3.5405405405405403e-5,
            3.2702702702702714e-5,
            -4.594_594_594_594_595e-5,
        ];
        // csinor = sqrt(Q/nused) = Rust normalised_rms
        let oracle_csinor: f64 = 2.075_819_784_513_525;

        let eps = 1e-10;
        for (i, &expected) in oracle_dx.iter().enumerate() {
            assert_abs_diff_eq!(result.element_correction[i], expected, epsilon = eps);
        }
        assert_abs_diff_eq!(result.normalised_rms, oracle_csinor, epsilon = eps);
    }
}
