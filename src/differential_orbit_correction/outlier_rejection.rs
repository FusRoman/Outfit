//! Projection-based outlier rejection for the differential-correction loop.
//!
//! After each Newton–Raphson iteration, each observation's chi-squared
//! contribution is recomputed using the *projected* residual variance — the
//! variance of the residual that is **not** explained by the parameter
//! correction.  This accounts for the covariance of the fitted orbit and
//! avoids systematically over-rejecting observations whose apparent residual
//! is partly absorbed by the correction.
//!
//! ## Projected residual variance
//!
//! For a single observation with partial-derivative rows
//! \\( g_\alpha \\) and \\( g_\delta \\) (each a \\( 1 \times 6 \\) vector),
//! the projected 2×2 residual covariance is
//!
//! \\[
//!   V = W^{-1} - g \cdot \Gamma \cdot g^\top
//! \\]
//!
//! where \\( W^{-1} \\) is the observation covariance and
//! \\( \Gamma = (G^\top W G)^{-1} \\) is the orbit covariance.
//!
//! The per-observation chi-squared is
//!
//! \\[
//!   \chi^2_i = \xi_i^\top \cdot V^{-1} \cdot \xi_i
//! \\]
//!
//! An observation is **rejected** when \\( \chi^2_i > \chi^2_{\text{reject}} \\)
//! and **re-admitted** when it was previously rejected (but not
//! [`ObsSelection::ForcedOut`]) and \\( \chi^2_i \le \chi^2_{\text{recover}} \\).
//!
//! ## Configuration
//!
//! See [`OutlierRejectionConfig`] for the two thresholds and
//! [`update_observation_selection`] for the function that applies them.

use nalgebra::{Matrix2, Vector2};

use crate::differential_orbit_correction::{
    least_square::{ObservationEquation, OrbitalUncertainty},
    obs_fit_data::{ObsFitData, ObsSelection},
};

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Tuning parameters for the projection-based outlier rejection step.
///
/// Both thresholds are expressed as dimensionless chi-squared values.
#[derive(Debug, Clone)]
pub struct OutlierRejectionConfig {
    /// Chi-squared threshold above which an observation is rejected.
    ///
    /// An active observation whose projected \\( \chi^2_i \\) exceeds this
    /// value is moved to [`ObsSelection::Rejected`].
    ///
    /// Default: `25.0` (≈ 5σ for uncorrelated RA/Dec, 2 degrees of freedom).
    pub chi_squared_rejection_threshold: f64,

    /// Chi-squared threshold at or below which a previously rejected
    /// observation is re-admitted.
    ///
    /// An observation in [`ObsSelection::Rejected`] whose projected
    /// \\( \chi^2_i \\) drops back to this value is moved to
    /// [`ObsSelection::Active`].  [`ObsSelection::ForcedOut`] observations
    /// are never re-admitted regardless of their chi-squared value.
    ///
    /// Default: `9.0` (≈ 3σ for uncorrelated RA/Dec, 2 degrees of freedom).
    pub chi_squared_recovery_threshold: f64,
}

impl Default for OutlierRejectionConfig {
    fn default() -> Self {
        Self {
            chi_squared_rejection_threshold: 25.0,
            chi_squared_recovery_threshold: 9.0,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Updates the selection flag of each observation using a projection-based
/// chi-squared criterion and returns the updated collection.
///
/// For each observation the function:
///
/// 1. Computes the projected 2×2 residual covariance
///    \\( V = W^{-1} - g \cdot \Gamma \cdot g^\top \\).
/// 2. Inverts \\( V \\) using its analytic 2×2 formula.
/// 3. Evaluates \\( \chi^2_i = \xi_i^\top V^{-1} \xi_i \\).
/// 4. Applies the rejection / recovery decision.
///
/// # Arguments
///
/// - `obs_fit_data` — per-observation fit data (immutable).
/// - `equations` — per-observation linearised equations from the last
///   [`single_iteration`](super::single_iteration::single_iteration) call.
///   Must have the same length as `obs_fit_data`.
/// - `uncertainty` — orbit covariance \\( \Gamma \\) from the last iteration.
/// - `config` — rejection and recovery thresholds.
///
/// # Returns
///
/// A tuple `(updated, num_changes)` where `updated` is a new `Vec<ObsFitData>`
/// with the selection flags updated and `num_changes` is the number of
/// selection changes (rejections + recoveries) made in this call.  The caller
/// uses `num_changes` to detect when the selection has stabilised
/// (i.e. when it equals `0`).
///
/// # Panics
///
/// Panics if `obs_fit_data.len() != equations.len()`.
pub fn update_observation_selection(
    obs_fit_data: &[ObsFitData],
    equations: &[ObservationEquation],
    uncertainty: &OrbitalUncertainty,
    config: &OutlierRejectionConfig,
) -> (Vec<ObsFitData>, usize) {
    assert_eq!(
        obs_fit_data.len(),
        equations.len(),
        "obs_fit_data and equations must have the same length"
    );

    let covariance = &uncertainty.covariance;

    obs_fit_data.iter().zip(equations.iter()).fold(
        (Vec::with_capacity(obs_fit_data.len()), 0_usize),
        |(mut acc, changes), (fit_data, eq)| {
            // ForcedOut observations are permanently excluded — never touch them.
            if fit_data.selection == ObsSelection::ForcedOut {
                acc.push(fit_data.clone());
                return (acc, changes);
            }

            // ── 1. Observation measurement covariance W⁻¹ (diagonal, 2×2) ─
            //
            // W_i = diag(w_α, w_δ)  →  W_i⁻¹ = diag(σ_α², σ_δ²)
            //
            // Cross-weight is zero for the standard uncorrelated case.
            let var_ra = fit_data.sigma_ra * fit_data.sigma_ra;
            let var_dec = fit_data.sigma_dec * fit_data.sigma_dec;
            let cov_cross = -fit_data.sigma_ra * fit_data.sigma_dec * eq.weight_cross
                / (eq.weight_ra * eq.weight_dec); // = 0 when weight_cross = 0

            // ── 2. Projection term g · Γ · gᵀ (2×2) ──────────────────────
            //
            //   g = [ g_α ]  (2×6)
            //       [ g_δ ]
            //
            //   proj[0,0] = g_αᵀ · Γ · g_α
            //   proj[1,1] = g_δᵀ · Γ · g_δ
            //   proj[0,1] = g_αᵀ · Γ · g_δ  (= proj[1,0] by symmetry)
            let g_alpha = &eq.d_ra_d_elem;
            let g_delta = &eq.d_dec_d_elem;

            let gamma_g_alpha = covariance * g_alpha; // 6-vector
            let gamma_g_delta = covariance * g_delta; // 6-vector

            let proj_aa = g_alpha.dot(&gamma_g_alpha);
            let proj_dd = g_delta.dot(&gamma_g_delta);
            let proj_ad = g_alpha.dot(&gamma_g_delta);

            // ── 3. Projected residual variance V = W⁻¹ − g·Γ·gᵀ ──────────
            let projected_var = Matrix2::new(
                var_ra - proj_aa,
                cov_cross - proj_ad,
                cov_cross - proj_ad,
                var_dec - proj_dd,
            );

            // ── 4. Invert V analytically (2×2) ────────────────────────────
            //
            // det(V) = V[0,0]·V[1,1] − V[0,1]²
            let det = projected_var[(0, 0)] * projected_var[(1, 1)]
                - projected_var[(0, 1)] * projected_var[(0, 1)];

            // Use a relative singularity check: |det| < ε × max diagonal²
            let scale = projected_var[(0, 0)].abs().max(projected_var[(1, 1)].abs());
            let singular_threshold = f64::EPSILON * scale * scale;
            if det.abs() < singular_threshold || scale == 0.0 {
                // Singular or numerically degenerate projected covariance — skip.
                acc.push(fit_data.clone());
                return (acc, changes);
            }

            let inv_projected_var = Matrix2::new(
                projected_var[(1, 1)] / det,
                -projected_var[(0, 1)] / det,
                -projected_var[(1, 0)] / det,
                projected_var[(0, 0)] / det,
            );

            // ── 5. Chi-squared ─────────────────────────────────────────────
            let residual = Vector2::new(fit_data.residual_ra, fit_data.residual_dec);
            let chi_squared = residual.dot(&(inv_projected_var * residual));

            // ── 6. Rejection / recovery decision ──────────────────────────
            let new_selection = match fit_data.selection {
                ObsSelection::Active if chi_squared > config.chi_squared_rejection_threshold => {
                    Some(ObsSelection::Rejected)
                }
                ObsSelection::Rejected if chi_squared <= config.chi_squared_recovery_threshold => {
                    Some(ObsSelection::Active)
                }
                _ => None,
            };

            let (updated_fit_data, delta) = match new_selection {
                Some(new_sel) => {
                    let mut updated = fit_data.clone();
                    updated.selection = new_sel;
                    (updated, 1)
                }
                None => (fit_data.clone(), 0),
            };

            acc.push(updated_fit_data);
            (acc, changes + delta)
        },
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod outlier_rejection_tests {
    use super::*;
    use crate::differential_orbit_correction::least_square::ObservationEquation;
    use nalgebra::{Matrix6, Vector6};

    // ── Helpers ───────────────────────────────────────────────────────────────

    fn identity_uncertainty() -> OrbitalUncertainty {
        OrbitalUncertainty {
            normal_matrix: Matrix6::identity(),
            covariance: Matrix6::identity(),
            inversion_succeeded: true,
        }
    }

    fn zero_uncertainty() -> OrbitalUncertainty {
        OrbitalUncertainty {
            normal_matrix: Matrix6::zeros(),
            covariance: Matrix6::zeros(),
            inversion_succeeded: false,
        }
    }

    /// Build a diagonal observation equation with given partial derivatives and
    /// residuals.
    fn make_eq(
        g_ra: Vector6<f64>,
        g_dec: Vector6<f64>,
        res_ra: f64,
        res_dec: f64,
        sigma: f64,
        active: bool,
    ) -> ObservationEquation {
        ObservationEquation::uncorrelated(g_ra, g_dec, res_ra, res_dec, sigma, sigma, active)
    }

    // ── Tests ─────────────────────────────────────────────────────────────────

    /// An observation with zero residuals must never be rejected.
    #[test]
    fn test_zero_residual_never_rejected() {
        let sigma = 1e-5;
        let eq = make_eq(Vector6::zeros(), Vector6::zeros(), 0.0, 0.0, sigma, true);
        let uncertainty = zero_uncertainty(); // Γ = 0 → proj_var = W⁻¹

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = 0.0;
        fit_data.residual_dec = 0.0;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );
        assert_eq!(changes, 0);
        assert_eq!(updated[0].selection, ObsSelection::Active);
    }

    /// A large residual (>> 5σ) must be rejected.
    #[test]
    fn test_large_residual_is_rejected() {
        let sigma = 1e-5;
        // Γ = 0 → no projection term → V = W⁻¹ = diag(σ², σ²)
        // residual = 100σ → χ² = (100σ/σ)² + (100σ/σ)² = 20000 >> 25
        let eq = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            100.0 * sigma,
            100.0 * sigma,
            sigma,
            true,
        );
        let uncertainty = zero_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = 100.0 * sigma;
        fit_data.residual_dec = 100.0 * sigma;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        assert_eq!(changes, 1);
        assert_eq!(updated[0].selection, ObsSelection::Rejected);
    }

    /// A previously rejected observation with small residuals must be
    /// re-admitted.
    #[test]
    fn test_small_residual_recovers_rejected_observation() {
        let sigma = 1e-5;
        let eq = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            0.5 * sigma, // small residual
            0.5 * sigma,
            sigma,
            false, // currently inactive
        );
        let uncertainty = zero_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = 0.5 * sigma;
        fit_data.residual_dec = 0.5 * sigma;
        fit_data.selection = ObsSelection::Rejected;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        // χ² = 0.5² + 0.5² = 0.5 ≤ 9 → recovered
        assert_eq!(changes, 1);
        assert_eq!(updated[0].selection, ObsSelection::Active);
    }

    /// A [`ObsSelection::ForcedOut`] observation must never change state.
    #[test]
    fn test_forced_out_is_never_changed() {
        let sigma = 1e-5;
        let eq = make_eq(Vector6::zeros(), Vector6::zeros(), 0.0, 0.0, sigma, false);
        let uncertainty = zero_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.selection = ObsSelection::ForcedOut;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        assert_eq!(changes, 0);
        assert_eq!(updated[0].selection, ObsSelection::ForcedOut);
    }

    /// Verify the projected variance computation against a hand-calculated
    /// oracle.
    ///
    /// Setup:
    /// - σ_α = σ_δ = 1e-5  →  W⁻¹ = diag(1e-10, 1e-10)
    /// - g_α = e₁, g_δ = e₂  (standard basis vectors)
    /// - Γ = identity
    /// - proj[0,0] = g_αᵀ · I · g_α = 1,  proj[1,1] = 1,  proj[0,1] = 0
    /// - V = diag(1e-10 − 1, 1e-10 − 1)  → det = (1e-10 − 1)²
    ///
    /// With residual_ra = residual_dec = 1e-5 = σ:
    /// - χ² = ξ² / (σ² − 1) + ξ² / (σ² − 1)   (negative denom → singular)
    ///
    /// Because Γ = I and the partials are unit vectors, the projected variance
    /// becomes negative (the orbit fully explains the residual).  The function
    /// must detect the singular / negative-determinant case and skip the
    /// observation without changing its state.
    #[test]
    fn test_singular_projected_variance_skipped() {
        let sigma = 1e-5;
        // g_α = e₀, g_δ = e₁ so proj = diag(1, 1)
        let mut g_ra = Vector6::zeros();
        g_ra[0] = 1.0;
        let mut g_dec = Vector6::zeros();
        g_dec[1] = 1.0;

        let eq = make_eq(g_ra, g_dec, sigma, sigma, sigma, true);
        // Γ = identity → proj_aa = 1 >> σ² → negative V diagonal → singular
        let uncertainty = identity_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = sigma;
        fit_data.residual_dec = sigma;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        // Cannot invert → observation state unchanged
        assert_eq!(changes, 0);
        assert_eq!(updated[0].selection, ObsSelection::Active);
    }

    /// Returns `0` changes when no threshold is crossed by any observation.
    #[test]
    fn test_no_change_when_within_thresholds() {
        let sigma = 1e-5;
        // residual ≈ 2σ → χ² ≈ 4+4 = 8 < 25 → no rejection
        let eq = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            2.0 * sigma,
            2.0 * sigma,
            sigma,
            true,
        );
        let uncertainty = zero_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = 2.0 * sigma;
        fit_data.residual_dec = 2.0 * sigma;

        let (updated, changes) = update_observation_selection(
            &[fit_data],
            &[eq],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        assert_eq!(changes, 0);
        assert_eq!(updated[0].selection, ObsSelection::Active);
    }

    /// Custom thresholds are respected.
    #[test]
    fn test_custom_thresholds() {
        let sigma = 1e-5;
        // residual = 3σ → χ² = 9+9 = 18
        // Default threshold 25 → no rejection
        // Custom threshold 16 → rejection
        let eq = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            3.0 * sigma,
            3.0 * sigma,
            sigma,
            true,
        );
        let uncertainty = zero_uncertainty();

        let mut fit_data = ObsFitData::new(sigma, sigma);
        fit_data.residual_ra = 3.0 * sigma;
        fit_data.residual_dec = 3.0 * sigma;

        let config = OutlierRejectionConfig {
            chi_squared_rejection_threshold: 16.0,
            chi_squared_recovery_threshold: 4.0,
        };

        let (updated, changes) =
            update_observation_selection(&[fit_data], &[eq], &uncertainty, &config);

        assert_eq!(changes, 1);
        assert_eq!(updated[0].selection, ObsSelection::Rejected);
    }

    /// `num_changes` counts correctly when multiple observations change state.
    #[test]
    fn test_multiple_changes_counted_correctly() {
        let sigma = 1e-5;
        let uncertainty = zero_uncertainty();

        // obs 0: active, small residual → no change
        let eq0 = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            sigma,
            sigma,
            sigma,
            true,
        );
        let mut fd0 = ObsFitData::new(sigma, sigma);
        fd0.residual_ra = sigma;
        fd0.residual_dec = sigma;

        // obs 1: active, large residual → reject
        let eq1 = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            10.0 * sigma,
            10.0 * sigma,
            sigma,
            true,
        );
        let mut fd1 = ObsFitData::new(sigma, sigma);
        fd1.residual_ra = 10.0 * sigma;
        fd1.residual_dec = 10.0 * sigma;

        // obs 2: rejected, tiny residual → recover
        let eq2 = make_eq(
            Vector6::zeros(),
            Vector6::zeros(),
            0.1 * sigma,
            0.1 * sigma,
            sigma,
            false,
        );
        let mut fd2 = ObsFitData::new(sigma, sigma);
        fd2.residual_ra = 0.1 * sigma;
        fd2.residual_dec = 0.1 * sigma;
        fd2.selection = ObsSelection::Rejected;

        let (updated, changes) = update_observation_selection(
            &[fd0, fd1, fd2],
            &[eq0, eq1, eq2],
            &uncertainty,
            &OutlierRejectionConfig::default(),
        );

        assert_eq!(changes, 2);
        assert_eq!(updated[0].selection, ObsSelection::Active);
        assert_eq!(updated[1].selection, ObsSelection::Rejected);
        assert_eq!(updated[2].selection, ObsSelection::Active);
    }
}
