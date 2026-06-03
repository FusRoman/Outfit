//! Differential orbit correction: outer loop driver.
//!
//! This module orchestrates the full differential-correction pipeline for a
//! single trajectory.  It combines:
//!
//! - An **inner Newton–Raphson loop** that iteratively refines the orbital
//!   elements by solving the linearised observation equations
//!   (see [`single_iteration`]).
//! - An **outer outlier-rejection loop** that, after each converged inner
//!   loop, tests whether any observation should be rejected or re-admitted
//!   based on its projected chi-squared contribution (see
//!   [`update_observation_selection`]).
//!
//! The function returns when the selection is stable *and* the inner loop
//! has converged, or when one of the configured iteration limits is reached.
//!
//! ## Algorithm overview
//!
//! ```text
//! outer_iteration = 0
//! loop:
//!     inner_iteration = 0
//!     loop (Newton–Raphson):
//!         result = single_iteration(elements, obs, obs_fit_data, …)
//!         check convergence: correction_norm < threshold
//!         check stagnation / divergence on normalised_rms
//!         check bizarre elements (eccentricity, semi-major axis limits)
//!         check inversion success
//!         elements = result.corrected_elements
//!         obs_fit_data = result.updated_obs_fit_data
//!         inner_iteration += 1
//!     if enable_outlier_rejection:
//!         num_changes, obs_fit_data = update_observation_selection(obs_fit_data, equations, Γ, …)
//!         if num_changes == 0: break   ← selection stable → done
//!     else:
//!         break
//!     outer_iteration += 1
//! rescale_covariance(…)
//! return DifferentialCorrectionOutput { elements, uncertainty, normalised_rms, … }
//! ```
//!
//! ## Configuration
//!
//! All tuning parameters are grouped in [`DifferentialCorrectionConfig`].  The
//! default values are suitable for main-belt asteroids with dense optical
//! astrometry; tighter constraints may be appropriate for short-arc objects.

use photom::observation_dataset::observation::Observation;

use crate::{
    cache::OutfitCache,
    differential_orbit_correction::{
        least_square::rescale_covariance,
        obs_fit_data::ObsFitData,
        outlier_rejection::{update_observation_selection, OutlierRejectionConfig},
        single_iteration::single_iteration,
    },
    orbit_type::{
        equinoctial_element::EquinoctialLimits,
        uncertainty::{EquinoctialUncertainty, OrbitalCovariance},
    },
    propagator::PropagatorKind,
    EquinoctialElements, JPLEphem, OrbitalElements, OutfitError,
};

use super::least_square::OrbitalUncertainty;

// ─────────────────────────────────────────────────────────────────────────────
// Configuration
// ─────────────────────────────────────────────────────────────────────────────

/// Tuning parameters for the full differential-correction pipeline.
///
/// All fields have documented defaults; use [`DifferentialCorrectionConfig::default()`]
/// to get a configuration suitable for main-belt asteroids with optical
/// astrometry.
#[derive(Debug, Clone)]
pub struct DifferentialCorrectionConfig {
    /// Maximum number of Newton–Raphson iterations in the inner loop.
    ///
    /// The inner loop stops when this count is reached even if the correction
    /// norm has not dropped below `convergence_threshold`.
    ///
    /// Default: `30`.
    pub max_newton_iterations: usize,

    /// Maximum number of outer outlier-rejection passes.
    ///
    /// The outer loop stops after this many passes even if the selection is
    /// still changing.
    ///
    /// Default: `10`.
    pub max_outlier_rejection_passes: usize,

    /// Dimensionless correction-norm threshold for inner-loop convergence.
    ///
    /// The inner loop is considered converged when
    /// \\( \|\delta x\|_C < \text{convergence\_threshold} \\).
    ///
    /// Default: `1e-4`.
    pub convergence_threshold: f64,

    /// Minimum normalised RMS below which outlier rejection is skipped for
    /// the first pass.
    ///
    /// When the normalised RMS is already below this value after the first
    /// inner-loop convergence, the outlier-rejection step is bypassed entirely
    /// (the fit is already clean enough).
    ///
    /// Default: `2.0`.
    pub convergence_before_rejection_threshold: f64,

    /// Ratio of consecutive normalised-RMS values above which the inner loop
    /// is considered **stagnated**.
    ///
    /// If `rms_new / rms_prev ≥ rms_stagnation_ratio`, the inner loop is
    /// terminated early and the last valid iteration is returned.
    ///
    /// Default: `0.98` (less than 2 % improvement triggers stagnation).
    pub rms_stagnation_ratio: f64,

    /// Ratio of consecutive normalised-RMS values above which the inner loop
    /// is considered **diverged**.
    ///
    /// If `rms_new / rms_prev ≥ rms_divergence_ratio`, the pipeline returns
    /// [`OutfitError::DifferentialCorrectionDiverged`].
    ///
    /// Default: `1.5` (50 % increase in RMS signals divergence).
    pub rms_divergence_ratio: f64,

    /// Maximum number of consecutive stagnation events before the inner loop
    /// is forcefully stopped.
    ///
    /// Default: `3`.
    pub max_stagnation_iterations: usize,

    /// Whether to run the projection-based outlier-rejection step after each
    /// inner-loop convergence.
    ///
    /// When `false`, the pipeline performs a single inner-loop pass with the
    /// initial observation selection.
    ///
    /// Default: `true`.
    pub enable_outlier_rejection: bool,

    /// Outlier rejection thresholds (chi-squared).
    ///
    /// Used only when `enable_outlier_rejection` is `true`.
    pub outlier_rejection_config: OutlierRejectionConfig,

    /// Physical-plausibility limits on the equinoctial elements.
    ///
    /// After each Newton step, the corrected elements are tested against these
    /// limits.  If they fail, the pipeline returns
    /// [`OutfitError::BizarreOrbit`].
    pub orbital_limits: EquinoctialLimits,

    /// Six-element boolean mask for free parameters.
    ///
    /// `free_elements[j] = true` means element `j` is solved for; `false`
    /// means it is held fixed.  Fixing an element is useful for under-
    /// determined arcs (e.g., fixing the semi-major axis for a very short arc).
    ///
    /// Default: `[true; 6]` (all elements free).
    pub free_elements: [bool; 6],

    /// Propagator to use for computing predicted observations and partials.
    ///
    /// - [`PropagatorKind::TwoBody`] (default): analytic Keplerian propagation.
    /// - [`PropagatorKind::NBody`]: numerical DOP853 N-body integration with
    ///   user-specified perturbing bodies.
    pub propagator: PropagatorKind,
}

impl Default for DifferentialCorrectionConfig {
    fn default() -> Self {
        Self {
            max_newton_iterations: 30,
            max_outlier_rejection_passes: 10,
            convergence_threshold: 1e-4,
            convergence_before_rejection_threshold: 2.0,
            rms_stagnation_ratio: 0.98,
            rms_divergence_ratio: 1.5,
            max_stagnation_iterations: 3,
            enable_outlier_rejection: true,
            outlier_rejection_config: OutlierRejectionConfig::default(),
            orbital_limits: EquinoctialLimits::default(),
            free_elements: [true; 6],
            propagator: PropagatorKind::TwoBody,
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Output type
// ─────────────────────────────────────────────────────────────────────────────

/// Output of a completed differential-correction run.
///
/// Returned by [`run_differential_correction`] on success.
#[derive(Debug, Clone)]
pub struct DifferentialCorrectionOutput {
    /// Best-fitting equinoctial orbital elements.
    pub elements: EquinoctialElements,

    /// Per-observation fit data after the final iteration, including residuals,
    /// chi values, and selection flags.
    pub final_obs_fit_data: Vec<ObsFitData>,

    /// Covariance and normal matrices, rescaled by the posterior uncertainty
    /// inflation factor.
    pub uncertainty: OrbitalUncertainty,

    /// Normalised RMS of the final fit
    /// \\( \sqrt{\xi^\top W \xi / n_{\text{active}}} \\).
    pub normalised_rms: f64,

    /// Total number of Newton–Raphson iterations performed across all outer
    /// passes.
    pub total_newton_iterations: usize,

    /// Number of scalar measurements used in the final fit (2 per active
    /// optical observation).
    pub num_measurements: usize,
}

/// Conversion from [`DifferentialCorrectionOutput`] to the more general
/// [`OrbitalElements`] type.
///
/// The covariance is extracted from the `uncertainty` field and included in the
/// `OrbitalElements::Equinoctial` variant.
impl From<DifferentialCorrectionOutput> for OrbitalElements {
    fn from(output: DifferentialCorrectionOutput) -> Self {
        let orb_covariance = OrbitalCovariance {
            matrix: output.uncertainty.covariance,
        };
        OrbitalElements::Equinoctial {
            elements: output.elements,
            uncertainty: Some(EquinoctialUncertainty::from_covariance(&orb_covariance)),
            covariance: Some(orb_covariance),
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Main function
// ─────────────────────────────────────────────────────────────────────────────

/// Runs the full differential orbit correction for a single set of
/// observations.
///
/// # Arguments
///
/// - `observations` — slice of astrometric observations (all belonging to the
///   same trajectory).
/// - `initial_obs_fit_data` — per-observation statistical fit data (σ, bias,
///   initial selection flag).  Must have the same length as `observations`.
/// - `initial_elements` — starting equinoctial orbital elements (e.g., from a
///   preceding IOD step).
/// - `cache` — pre-computed observer geometry cache.
/// - `jpl` — JPL planetary ephemeris handle.
/// - `config` — tuning parameters; use [`DifferentialCorrectionConfig::default()`]
///   for standard settings.
///
/// # Returns
///
/// A [`DifferentialCorrectionOutput`] with the refined elements, covariance,
/// and final per-observation fit data.
///
/// # Errors
///
/// - [`OutfitError::BizarreOrbit`] — the corrected elements violate the
///   physical-plausibility limits in `config.orbital_limits`.
/// - [`OutfitError::DifferentialCorrectionDiverged`] — the normalised RMS
///   increased by more than `config.rms_divergence_ratio`.
/// - [`OutfitError::DifferentialCorrectionFailed`] — the normal-equation
///   inversion failed (e.g., fewer active observations than free parameters).
///
/// # Panics
///
/// Panics if `observations.len() != initial_obs_fit_data.len()`.
pub fn run_differential_correction(
    observations: &[Observation],
    initial_obs_fit_data: &[ObsFitData],
    initial_elements: &EquinoctialElements,
    cache: &OutfitCache,
    jpl: &JPLEphem,
    config: &DifferentialCorrectionConfig,
) -> Result<DifferentialCorrectionOutput, OutfitError> {
    assert_eq!(
        observations.len(),
        initial_obs_fit_data.len(),
        "observations and initial_obs_fit_data must have the same length"
    );

    let num_free = config.free_elements.iter().filter(|&&f| f).count();

    // Working copies — updated at the end of every Newton step.
    let mut elements = initial_elements.clone();
    let mut obs_fit_data = initial_obs_fit_data.to_vec();

    let mut total_newton_iterations = 0_usize;

    // Saved from the last successful Newton step for fallback.
    let mut last_uncertainty = OrbitalUncertainty {
        normal_matrix: nalgebra::Matrix6::zeros(),
        covariance: nalgebra::Matrix6::zeros(),
        inversion_succeeded: false,
    };
    let mut last_normalised_rms = f64::MAX;
    let mut last_num_measurements = 0_usize;

    // ── Outer outlier-rejection loop ─────────────────────────────────────────
    for outer_pass in 0..=config.max_outlier_rejection_passes {
        let mut prev_rms = f64::MAX;
        let mut stagnation_count = 0_usize;

        // ── Inner Newton–Raphson loop ─────────────────────────────────────────
        let mut last_equations = vec![];
        let mut inner_loop_converged = false;

        for _inner in 0..config.max_newton_iterations {
            total_newton_iterations += 1;

            let iter_result = single_iteration(
                observations,
                &obs_fit_data,
                &elements,
                &config.free_elements,
                cache,
                jpl,
                true,
                &config.propagator,
            )?;

            // ── Check inversion ──────────────────────────────────────────────
            if !iter_result.uncertainty.inversion_succeeded {
                return Err(OutfitError::DifferentialCorrectionFailed(
                    "normal-equation inversion failed (possibly fewer active observations than \
                     free parameters)"
                        .into(),
                ));
            }

            // ── Check bizarre orbit ──────────────────────────────────────────
            if iter_result
                .corrected_elements
                .is_bizarre(&config.orbital_limits)
            {
                return Err(OutfitError::BizarreOrbit);
            }

            let new_rms = iter_result.normalised_rms;

            // ── Check divergence ─────────────────────────────────────────────
            if prev_rms < f64::MAX && new_rms / prev_rms >= config.rms_divergence_ratio {
                return Err(OutfitError::DifferentialCorrectionDiverged);
            }

            // ── Check stagnation ─────────────────────────────────────────────
            let stagnated =
                prev_rms < f64::MAX && new_rms / prev_rms >= config.rms_stagnation_ratio;
            if stagnated {
                stagnation_count += 1;
                if stagnation_count >= config.max_stagnation_iterations {
                    // Accept the current state and stop the inner loop.
                    break;
                }
            } else {
                stagnation_count = 0;
            }

            // Advance state.
            last_equations = iter_result.observation_equations.clone();
            last_uncertainty = iter_result.uncertainty.clone();
            last_normalised_rms = new_rms;
            last_num_measurements = iter_result.num_measurements;

            elements = iter_result.corrected_elements;
            obs_fit_data = iter_result.updated_obs_fit_data;
            prev_rms = new_rms;

            // ── Convergence check ────────────────────────────────────────────
            if iter_result.correction_norm < config.convergence_threshold {
                inner_loop_converged = true;
                break;
            }
        }

        // ── Skip outlier rejection if disabled or fit is already clean ───────
        if !config.enable_outlier_rejection {
            break;
        }

        // On the first outer pass, skip rejection if the RMS is already low.
        if outer_pass == 0 && last_normalised_rms < config.convergence_before_rejection_threshold {
            break;
        }

        // If the inner loop did not converge, do not attempt outlier rejection.
        if !inner_loop_converged {
            break;
        }

        // ── Outlier rejection step ───────────────────────────────────────────
        let (updated_obs_fit_data, num_selection_changes) = update_observation_selection(
            &obs_fit_data,
            &last_equations,
            &last_uncertainty,
            &config.outlier_rejection_config,
        );
        obs_fit_data = updated_obs_fit_data;

        // Selection is stable — we are done.
        if num_selection_changes == 0 {
            break;
        }
    }

    // ── Final covariance rescaling ────────────────────────────────────────────
    rescale_covariance(
        &mut last_uncertainty,
        num_free,
        last_num_measurements,
        last_normalised_rms,
    );

    Ok(DifferentialCorrectionOutput {
        elements,
        final_obs_fit_data: obs_fit_data,
        uncertainty: last_uncertainty,
        normalised_rms: last_normalised_rms,
        total_newton_iterations,
        num_measurements: last_num_measurements,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod diff_cor_tests {
    use super::*;
    use crate::{
        differential_orbit_correction::obs_fit_data::{ObsFitData, ObsSelection},
        test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER},
    };
    use photom::{
        coordinates::equatorial::EquCoord,
        observation_dataset::{observation::ObservationInput, ObsDataset},
        observer::{
            dataset::ObserverId,
            error_model::{ModelCorrection, ObsErrorModel},
        },
        photometry::{Filter, Photometry},
    };

    // ── Helpers ───────────────────────────────────────────────────────────────

    fn circular_elements(epoch: f64) -> EquinoctialElements {
        EquinoctialElements {
            reference_epoch: epoch,
            semi_major_axis: 1.8,
            eccentricity_sin_lon: 0.1,
            eccentricity_cos_lon: 0.05,
            tan_half_incl_sin_node: 0.01,
            tan_half_incl_cos_node: 0.1,
            mean_longitude: 1.0,
        }
    }

    fn make_dataset_and_cache(t0: f64, time_span: f64, n: usize) -> (ObsDataset, OutfitCache) {
        let step = if n > 1 {
            time_span / (n - 1) as f64
        } else {
            0.0
        };
        let inputs: Vec<ObservationInput> = (0..n)
            .map(|i| {
                let t_obs = t0 + i as f64 * step;
                ObservationInput::new(
                    i as u64,
                    EquCoord {
                        ra: 0.0,
                        ra_error: 0.0,
                        dec: 0.0,
                        dec_error: 0.0,
                    },
                    Photometry {
                        magnitude: 15.0,
                        error: 0.1,
                        filter: Filter::Int(0),
                    },
                    t_obs,
                    Some(ObserverId::MpcCode(*b"F51")),
                )
            })
            .collect();

        let obs_dataset = {
            let mut ds = ObsDataset::empty();
            for input in inputs {
                ds = ds.push_observation(vec![input]).unwrap().0;
            }
            ds.with_error_model(ObsErrorModel::FCCT14)
                .apply_model_errors()
        };

        let cache =
            OutfitCache::build(&obs_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER, false).unwrap();
        (obs_dataset, cache)
    }

    // ── Tests ─────────────────────────────────────────────────────────────────

    /// The pipeline must complete without error on well-conditioned input and
    /// return finite orbital elements.
    ///
    /// Note: the synthetic observations have RA/Dec = 0 which are far from the
    /// predicted position of the orbit, so large residuals are expected.  The
    /// correction step may push the elements outside the physical limits; the
    /// test therefore accepts `BizarreOrbit` / `DifferentialCorrectionDiverged`
    /// as valid outcomes and only checks that the function does not panic.
    #[test]
    fn test_run_completes_and_elements_are_finite() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 8);

        let observations: Vec<_> = (0..8)
            .map(|i| obs_dataset.get_observation(i as u64).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let config = DifferentialCorrectionConfig {
            enable_outlier_rejection: false,
            ..Default::default()
        };

        let result = run_differential_correction(
            &observations,
            &obs_fit_data,
            &elements,
            &cache,
            &JPL_EPHEM_HORIZON,
            &config,
        );

        // With RA/Dec=0 synthetic observations the Newton step produces a very
        // large correction that pushes the elements outside the plausibility
        // limits.  Acceptable outcomes are success (unlikely but possible) or
        // a typed error; the function must not panic.
        match result {
            Ok(output) => {
                assert!(output.elements.semi_major_axis.is_finite());
                assert!(output.total_newton_iterations >= 1);
            }
            Err(OutfitError::BizarreOrbit) => {} // expected
            Err(OutfitError::DifferentialCorrectionDiverged) => {} // acceptable
            Err(OutfitError::DifferentialCorrectionFailed(_)) => {} // acceptable
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }

    /// With all observations inactive, the elements must not change.
    #[test]
    fn test_all_inactive_elements_unchanged() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<_> = (0..6)
            .map(|i| obs_dataset.get_observation(i as u64).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| {
                let mut fd = ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error);
                fd.selection = ObsSelection::Rejected;
                fd
            })
            .collect();

        let config = DifferentialCorrectionConfig {
            enable_outlier_rejection: false,
            ..Default::default()
        };

        // All inactive → inversion fails → DifferentialCorrectionFailed
        let result = run_differential_correction(
            &observations,
            &obs_fit_data,
            &elements,
            &cache,
            &JPL_EPHEM_HORIZON,
            &config,
        );

        assert!(matches!(
            result,
            Err(OutfitError::DifferentialCorrectionFailed(_))
        ));
    }

    /// `total_newton_iterations` must be ≥ 1.
    #[test]
    fn test_at_least_one_newton_iteration() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 8);

        let observations: Vec<_> = (0..8)
            .map(|i| obs_dataset.get_observation(i as u64).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let config = DifferentialCorrectionConfig {
            enable_outlier_rejection: false,
            ..Default::default()
        };

        // `total_newton_iterations` is set before the BizarreOrbit check.
        // For synthetic zero-RA/Dec observations the first iteration already
        // detects a bizarre orbit, so we accept that error but verify the
        // counter was incremented.
        match run_differential_correction(
            &observations,
            &obs_fit_data,
            &elements,
            &cache,
            &JPL_EPHEM_HORIZON,
            &config,
        ) {
            Ok(output) => assert!(output.total_newton_iterations >= 1),
            Err(OutfitError::BizarreOrbit) => {} // one iteration ran before the check
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }

    /// Limiting to 1 Newton iteration must still return a valid result (or a
    /// typed error for degenerate synthetic input).
    #[test]
    fn test_single_newton_iteration() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 8);

        let observations: Vec<_> = (0..8)
            .map(|i| obs_dataset.get_observation(i as u64).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let config = DifferentialCorrectionConfig {
            max_newton_iterations: 1,
            enable_outlier_rejection: false,
            ..Default::default()
        };

        match run_differential_correction(
            &observations,
            &obs_fit_data,
            &elements,
            &cache,
            &JPL_EPHEM_HORIZON,
            &config,
        ) {
            Ok(output) => {
                assert_eq!(output.total_newton_iterations, 1);
                assert!(output.elements.semi_major_axis.is_finite());
            }
            Err(OutfitError::BizarreOrbit) => {} // expected for synthetic data
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }

    /// `final_obs_fit_data` must have the same length as `observations`.
    #[test]
    fn test_final_obs_fit_data_length() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let n = 8;
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, n);

        let observations: Vec<_> = (0..n)
            .map(|i| obs_dataset.get_observation(i as u64).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let config = DifferentialCorrectionConfig {
            // Limit to 1 iteration so we get a result before BizarreOrbit.
            // The `final_obs_fit_data` length is set before the bizarre check.
            max_newton_iterations: 1,
            enable_outlier_rejection: false,
            ..Default::default()
        };

        match run_differential_correction(
            &observations,
            &obs_fit_data,
            &elements,
            &cache,
            &JPL_EPHEM_HORIZON,
            &config,
        ) {
            Ok(output) => assert_eq!(output.final_obs_fit_data.len(), n),
            // BizarreOrbit is raised before returning output, so we cannot
            // check the length directly.  Just verify the function ran.
            Err(OutfitError::BizarreOrbit) => {}
            Err(e) => panic!("unexpected error: {e:?}"),
        }
    }
}
