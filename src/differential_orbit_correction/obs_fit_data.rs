//! Per-observation statistical fit data: uncertainties, biases, residuals and
//! selection flag.
//!
//! [`ObsFitData`] carries the statistical metadata that accompanies an
//! astrometric observation through the differential-correction loop without
//! modifying the raw observation itself.
//!
//! ## Typical lifecycle
//!
//! 1. Build one [`ObsFitData`] per observation from the error model
//!    (σ_α, σ_δ and optional biases).
//! 2. Pass a `&[ObsFitData]` (together with the matching `&[Observation]`) to
//!    [`single_iteration`](super::single_iteration::single_iteration).
//! 3. Inspect the [`super::single_iteration::SingleIterationResult::updated_obs_fit_data`] returned to
//!    read the updated residuals and `chi` values.

// ─────────────────────────────────────────────────────────────────────────────
// Selection flag
// ─────────────────────────────────────────────────────────────────────────────

/// Participation flag for a single observation in the differential-correction
/// fit.
///
/// | Value | Variant |
/// |---|---|
/// | inactive | [`ObsSelection::Rejected`] |
/// | active | [`ObsSelection::Active`] |
/// | excluded | [`ObsSelection::ForcedOut`] |
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObsSelection {
    /// Observation contributes to the fit.
    Active,
    /// Automatically rejected by the outlier-rejection step.
    Rejected,
    /// Manually excluded; never re-activated.
    ForcedOut,
}

impl ObsSelection {
    /// Returns `true` if this observation should contribute to the normal
    /// equations (i.e. it is [`Active`](ObsSelection::Active)).
    #[inline]
    pub fn is_active(self) -> bool {
        matches!(self, ObsSelection::Active)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main struct
// ─────────────────────────────────────────────────────────────────────────────

/// Statistical fit data for a single astrometric observation.
///
/// Holds the per-observation uncertainties, systematic biases, residuals and
/// selection flag used during the differential-correction loop.  It is *separate*
/// from the raw [`photom::observation_dataset::observation::Observation`] so
/// that the loop can update residuals and selection flags without mutating
/// immutable observation data.
///
/// ## Units
///
/// All angular quantities (sigmas, biases, residuals) are in **radians**.
///
/// ## Default biases
///
/// Set `bias_ra = 0.0` and `bias_dec = 0.0` unless a catalogue or night-block
/// debiasing step has produced non-zero values.
#[derive(Debug, Clone)]
pub struct ObsFitData {
    /// Right-ascension uncertainty σ_α \[rad\].
    pub sigma_ra: f64,
    /// Declination uncertainty σ_δ \[rad\].
    pub sigma_dec: f64,
    /// Right-ascension systematic bias \[rad\] (subtracted from the observed
    /// RA before computing residuals).
    pub bias_ra: f64,
    /// Declination systematic bias \[rad\].
    pub bias_dec: f64,
    /// Right-ascension residual
    /// \\( \xi_\alpha = \alpha_{\text{obs}} - \text{bias}_\alpha - \alpha_{\text{calc}} \\)
    /// \[rad\], set by [`single_iteration`](super::single_iteration::single_iteration).
    pub residual_ra: f64,
    /// Declination residual
    /// \\( \xi_\delta = \delta_{\text{obs}} - \text{bias}_\delta - \delta_{\text{calc}} \\)
    /// \[rad\], set by [`single_iteration`](super::single_iteration::single_iteration).
    pub residual_dec: f64,
    /// Whether this observation participates in the fit.
    pub selection: ObsSelection,
    /// \\( \sqrt{\chi^2} \\) contribution of this observation, filled after
    /// each call to [`single_iteration`](super::single_iteration::single_iteration).
    pub chi: f64,
}

impl ObsFitData {
    /// Constructs an active, unbiased entry from the observation uncertainties.
    ///
    /// This is the standard constructor for the first iteration: residuals and
    /// `chi` are initialised to `0.0` and will be populated by
    /// [`single_iteration`](super::single_iteration::single_iteration).
    ///
    /// # Arguments
    ///
    /// - `sigma_ra` — right-ascension uncertainty \[rad\].
    /// - `sigma_dec` — declination uncertainty \[rad\].
    pub fn new(sigma_ra: f64, sigma_dec: f64) -> Self {
        Self {
            sigma_ra,
            sigma_dec,
            bias_ra: 0.0,
            bias_dec: 0.0,
            residual_ra: 0.0,
            residual_dec: 0.0,
            selection: ObsSelection::Active,
            chi: 0.0,
        }
    }

    /// Returns `true` if this observation is [`Active`](ObsSelection::Active).
    #[inline]
    pub fn is_active(&self) -> bool {
        self.selection.is_active()
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod obs_fit_data_tests {
    use super::*;

    #[test]
    fn test_new_is_active_with_zero_residuals() {
        let fit = ObsFitData::new(1e-5, 1.5e-5);
        assert!(fit.is_active());
        assert_eq!(fit.residual_ra, 0.0);
        assert_eq!(fit.residual_dec, 0.0);
        assert_eq!(fit.bias_ra, 0.0);
        assert_eq!(fit.bias_dec, 0.0);
        assert_eq!(fit.chi, 0.0);
        assert_eq!(fit.sigma_ra, 1e-5);
        assert_eq!(fit.sigma_dec, 1.5e-5);
    }

    #[test]
    fn test_rejected_is_not_active() {
        let mut fit = ObsFitData::new(1e-5, 1e-5);
        fit.selection = ObsSelection::Rejected;
        assert!(!fit.is_active());
    }

    #[test]
    fn test_forced_out_is_not_active() {
        let mut fit = ObsFitData::new(1e-5, 1e-5);
        fit.selection = ObsSelection::ForcedOut;
        assert!(!fit.is_active());
    }

    #[test]
    fn test_obs_selection_is_active() {
        assert!(ObsSelection::Active.is_active());
        assert!(!ObsSelection::Rejected.is_active());
        assert!(!ObsSelection::ForcedOut.is_active());
    }
}
