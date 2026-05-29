//! One Newton–Raphson iteration of the differential orbit correction.
//!
//! A single call to [`single_iteration`] performs one full pass:
//!
//! 1. For each active observation, compute the predicted (RA, DEC) and the
//!    partial derivatives ∂(α,δ)/∂(elements) via the two-body propagator.
//! 2. Compute the astrometric residuals and store them in the returned
//!    [`ObsFitData`] entries.
//! 3. Assemble the normal equations and solve for the element correction
//!    δx = Γ · GᵀWξ.
//! 4. Optionally apply δx to the input elements, returning the corrected
//!    [`EquinoctialElements`].
//! 5. Compute the correction norm ‖δx‖_C = √(δxᵀ · C · δx) where C is the
//!    normal matrix GᵀWG.
//!
//! ## Functional style
//!
//! All outputs are returned in [`SingleIterationResult`]; inputs are borrowed
//! immutably.  The caller obtains new values for the orbital elements and the
//! per-observation fit data by destructuring the result.
//!
//! ```text
//! let result = single_iteration(&obs, &obs_fit_data, &elements, ...)?;
//! let elements    = result.corrected_elements;
//! let obs_fit_data = result.updated_obs_fit_data;
//! ```
//!
//! ## Observation failures
//!
//! If the orbit propagator fails for an individual observation (e.g. light-time
//! non-convergence), that observation is silently skipped and treated as inactive
//! for this iteration.  All other observations are processed normally.

use photom::observation_dataset::observation::Observation;

use crate::{
    cache::OutfitCache,
    differential_orbit_correction::{
        least_square::{
            angular_diff, solve_weighted_least_squares, ObservationEquation, OrbitalUncertainty,
        },
        obs_fit_data::ObsFitData,
    },
    observation_ephemeris::ObservationEphemeris,
    propagator::PropagatorKind,
    EquinoctialElements, JPLEphem, OutfitError,
};

// ─────────────────────────────────────────────────────────────────────────────
// Result type
// ─────────────────────────────────────────────────────────────────────────────

/// Output of a single differential-correction iteration.
///
/// This struct bundles every quantity produced by [`single_iteration`]:
/// the corrected orbital elements, updated per-observation fit data,
/// the per-observation linearised equations, convergence diagnostics,
/// and uncertainty matrices.
///
/// ## Correction norm
///
/// `correction_norm` is the dimensionless scalar
///
/// \\[
///   \|\delta x\|_C = \sqrt{\delta x^\top \cdot C \cdot \delta x}
/// \\]
///
/// where \\(C = G^\top W G\\) is the normal matrix.  It provides a unit-free
/// measure of the step size in the parameter space and drives the convergence
/// check.
///
/// ## Observation equations
///
/// `observation_equations` mirrors the per-observation linearised system built
/// during this iteration.  The `i`-th entry contains the partial derivatives
/// and weights for the `i`-th input observation.  Downstream steps such as
/// outlier rejection use these to compute the projected residual variance
/// \\( g \cdot \Gamma \cdot g^\top \\) without rerunning the propagator.
#[derive(Debug, Clone)]
pub struct SingleIterationResult {
    /// Orbital elements after applying the correction δx (or unchanged if
    /// `apply_correction` was `false`).
    pub corrected_elements: EquinoctialElements,
    /// Per-observation fit data with updated residuals and chi values.
    ///
    /// The `i`-th entry corresponds to the `i`-th input observation.
    /// Observations that failed to propagate keep their residuals from the
    /// previous iteration (or `0.0` on the first call).
    pub updated_obs_fit_data: Vec<ObsFitData>,
    /// Per-observation linearised equations built during this iteration.
    ///
    /// The `i`-th entry holds the partial derivatives
    /// \\( \partial\alpha/\partial\text{elem} \\),
    /// \\( \partial\delta/\partial\text{elem} \\), residuals, and weights for
    /// the `i`-th input observation.  Inactive observations carry a zero
    /// placeholder equation.
    pub observation_equations: Vec<ObservationEquation>,
    /// Dimensionless correction norm \\( \|\delta x\|_C \\).
    pub correction_norm: f64,
    /// Normalised RMS residual \\( \sqrt{\xi^\top W \xi / n} \\).
    pub normalised_rms: f64,
    /// Covariance and normal matrices from the inversion step.
    pub uncertainty: OrbitalUncertainty,
    /// Number of scalar measurements used (2 per active optical observation).
    pub num_measurements: usize,
}

// ─────────────────────────────────────────────────────────────────────────────
// Public function
// ─────────────────────────────────────────────────────────────────────────────

/// Performs one Newton–Raphson iteration of the differential orbit correction.
///
/// # Arguments
///
/// - `observations` — slice of astrometric observations (order is preserved in
///   the returned [`SingleIterationResult::updated_obs_fit_data`]).
/// - `obs_fit_data` — per-observation statistical fit data (σ, bias, selection
///   flag).  Must have the same length as `observations`.
/// - `elements` — current equinoctial orbital elements.
/// - `free_elements` — six-element boolean mask; `true` means the corresponding
///   element is solved for, `false` means it is held fixed.
/// - `cache` — pre-computed observer geometry cache.
/// - `jpl` — JPL planetary ephemeris handle.
/// - `apply_correction` — if `true`, the element correction δx is applied to
///   produce [`SingleIterationResult::corrected_elements`]; if `false` only
///   the covariance is computed (matrix-only mode).
///
/// # Errors
///
/// Returns [`OutfitError`] only if the normal-equation solver itself fails
/// (e.g. the normal matrix is identically zero because all observations are
/// inactive).  Per-observation propagation failures are handled gracefully by
/// skipping the failing observation.
///
/// # Panics
///
/// Panics if `observations.len() != obs_fit_data.len()`.
#[allow(clippy::too_many_arguments)]
pub fn single_iteration(
    observations: &[Observation],
    obs_fit_data: &[ObsFitData],
    elements: &EquinoctialElements,
    free_elements: &[bool; 6],
    cache: &OutfitCache,
    jpl: &JPLEphem,
    apply_correction: bool,
    propagator: &PropagatorKind,
) -> Result<SingleIterationResult, OutfitError> {
    assert_eq!(
        observations.len(),
        obs_fit_data.len(),
        "observations and obs_fit_data must have the same length"
    );

    // ── 1. Build observation equations ───────────────────────────────────────
    //
    // For each observation, attempt to compute predicted (RA, DEC) and element
    // partials.  Failures are logged and the observation is marked inactive.

    // Collect (ObservationEquation, updated ObsFitData) pairs.
    let (equations, updated_obs_fit_data): (Vec<ObservationEquation>, Vec<ObsFitData>) =
        observations
            .iter()
            .zip(obs_fit_data.iter())
            .map(|(obs, fit_data)| {
                // Inactive observations: build a zero-weight placeholder and keep
                // the existing residuals unchanged.
                if !fit_data.is_active() {
                    return (
                        ObservationEquation::uncorrelated(
                            nalgebra::Vector6::zeros(),
                            nalgebra::Vector6::zeros(),
                            0.0,
                            0.0,
                            1.0, // dummy sigma — weight ignored because active=false
                            1.0,
                            false,
                        ),
                        fit_data.clone(),
                    );
                }

                // Propagate orbit and compute predicted position + element partials.
                let partials_result = match propagator {
                    PropagatorKind::TwoBody => {
                        obs.compute_obs_and_partials_2body(cache, jpl, elements)
                    }
                    PropagatorKind::NBody(nbody_config) => {
                        obs.compute_obs_and_partials_nbody(cache, jpl, elements, nbody_config)
                    }
                };
                match partials_result {
                    Ok(partials) => {
                        // Residual RA: angular difference in (−π, π] accounting for wrapping.
                        // The observed RA is corrected for the catalogue bias before differencing.
                        let obs_ra_debiased = obs.equ_coord().ra - fit_data.bias_ra;
                        let residual_ra = angular_diff(obs_ra_debiased, partials.ra);

                        // Residual Dec: simple difference (no wrapping needed).
                        let obs_dec_debiased = obs.equ_coord().dec - fit_data.bias_dec;
                        let residual_dec = obs_dec_debiased - partials.dec;

                        // chi = normalised residual magnitude for this observation.
                        let chi = ((residual_ra / fit_data.sigma_ra).powi(2)
                            + (residual_dec / fit_data.sigma_dec).powi(2))
                        .sqrt();

                        let eq = ObservationEquation::uncorrelated(
                            partials.d_ra_d_elem,
                            partials.d_dec_d_elem,
                            residual_ra,
                            residual_dec,
                            fit_data.sigma_ra,
                            fit_data.sigma_dec,
                            true,
                        );

                        let updated = ObsFitData {
                            residual_ra,
                            residual_dec,
                            chi,
                            ..fit_data.clone()
                        };

                        (eq, updated)
                    }
                    Err(err) => {
                        eprintln!(
                        "single_iteration: propagation failed for observation at MJD {:.4}: {}. \
                         Observation skipped for this iteration.",
                        obs.mjd_tt(),
                        err
                    );
                        // Keep previous residuals; mark as inactive for this
                        // iteration only (do not change the selection flag).
                        (
                            ObservationEquation::uncorrelated(
                                nalgebra::Vector6::zeros(),
                                nalgebra::Vector6::zeros(),
                                0.0,
                                0.0,
                                1.0, // dummy sigma
                                1.0,
                                false,
                            ),
                            fit_data.clone(),
                        )
                    }
                }
            })
            .unzip();

    // ── 2. Solve normal equations ─────────────────────────────────────────────
    let ls_result = solve_weighted_least_squares(&equations, free_elements)?;

    // ── 3. Correction norm ‖δx‖_C = √(δxᵀ · C · δx) ─────────────────────────
    let dx = &ls_result.element_correction;
    let c = &ls_result.uncertainty.normal_matrix;
    let correction_norm = (dx.dot(&(c * dx))).sqrt();

    // ── 4. Apply correction to elements ──────────────────────────────────────
    let corrected_elements = if apply_correction {
        let mut new_coord = [
            elements.semi_major_axis,
            elements.eccentricity_sin_lon,
            elements.eccentricity_cos_lon,
            elements.tan_half_incl_sin_node,
            elements.tan_half_incl_cos_node,
            elements.mean_longitude,
        ];
        for j in 0..6 {
            if free_elements[j] {
                new_coord[j] += dx[j];
            }
        }
        EquinoctialElements {
            reference_epoch: elements.reference_epoch,
            semi_major_axis: new_coord[0],
            eccentricity_sin_lon: new_coord[1],
            eccentricity_cos_lon: new_coord[2],
            tan_half_incl_sin_node: new_coord[3],
            tan_half_incl_cos_node: new_coord[4],
            mean_longitude: new_coord[5],
        }
    } else {
        elements.clone()
    };

    Ok(SingleIterationResult {
        corrected_elements,
        updated_obs_fit_data,
        observation_equations: equations,
        correction_norm,
        normalised_rms: ls_result.normalised_rms,
        uncertainty: ls_result.uncertainty,
        num_measurements: ls_result.num_measurements,
    })
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod single_iteration_tests {
    use super::*;
    use crate::{
        differential_orbit_correction::obs_fit_data::ObsSelection,
        test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER},
    };
    use approx::assert_abs_diff_eq;
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

    /// Build a dataset with `n` observations spread over `time_span` days
    /// starting at `t0` (MJD TT), all from observatory F51 (Pan-STARRS).
    fn make_dataset_and_cache(t0: f64, time_span: f64, n: usize) -> (ObsDataset, OutfitCache) {
        let step = if n > 1 {
            time_span / (n - 1) as f64
        } else {
            0.0
        };

        // Use the propagator to get realistic (RA, DEC) for each epoch
        // so that residuals start near zero.
        let inputs: Vec<ObservationInput> = (0..n)
            .map(|i| {
                let t_obs = t0 + i as f64 * step;
                // Dummy RA/Dec — will be replaced after cache is built.
                // We use 0.0 here and fix them in the test bodies as needed.
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

    /// When all observations are inactive, the correction must be zero and the
    /// corrected elements must be unchanged.
    #[test]
    fn test_all_inactive_gives_zero_correction() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<Observation> = (0..6)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData {
                sigma_ra: obs.equ_coord().ra_error,
                sigma_dec: obs.equ_coord().dec_error,
                bias_ra: 0.0,
                bias_dec: 0.0,
                residual_ra: 0.0,
                residual_dec: 0.0,
                selection: ObsSelection::Rejected, // all inactive
                chi: 0.0,
            })
            .collect();

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &[true; 6],
            &cache,
            &JPL_EPHEM_HORIZON,
            true,
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        // No active observations → correction is zero
        assert_abs_diff_eq!(result.correction_norm, 0.0, epsilon = 1e-15);
        assert_abs_diff_eq!(result.normalised_rms, 0.0, epsilon = 1e-15);
        assert_eq!(result.num_measurements, 0);

        // Elements must be unchanged (delta applied was zero)
        assert_abs_diff_eq!(
            result.corrected_elements.semi_major_axis,
            elements.semi_major_axis,
            epsilon = 1e-15
        );
    }

    /// When `apply_correction = false`, the returned elements must be
    /// identical to the input elements regardless of the residuals.
    #[test]
    fn test_matonly_does_not_change_elements() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<Observation> = (0..6)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &[true; 6],
            &cache,
            &JPL_EPHEM_HORIZON,
            false, // matonly
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        assert_abs_diff_eq!(
            result.corrected_elements.semi_major_axis,
            elements.semi_major_axis,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            result.corrected_elements.mean_longitude,
            elements.mean_longitude,
            epsilon = 1e-15
        );
    }

    /// Updated weights slice must have the same length as the input.
    #[test]
    fn test_updated_weights_length_matches_input() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 4);

        let observations: Vec<Observation> = (0..4)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &[true; 6],
            &cache,
            &JPL_EPHEM_HORIZON,
            true,
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        assert_eq!(result.updated_obs_fit_data.len(), observations.len());
    }

    /// A fixed element must not change even when the correction is non-zero for
    /// free elements.
    #[test]
    fn test_fixed_element_not_corrected() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<Observation> = (0..6)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        // Fix element 0 (semi_major_axis)
        let mut free = [true; 6];
        free[0] = false;

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &free,
            &cache,
            &JPL_EPHEM_HORIZON,
            true,
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        assert_abs_diff_eq!(
            result.corrected_elements.semi_major_axis,
            elements.semi_major_axis,
            epsilon = 1e-15,
        );
    }

    /// The selection flag of inactive observations must be preserved unchanged
    /// in `updated_weights`.
    #[test]
    fn test_inactive_selection_flag_preserved() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<Observation> = (0..6)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let mut obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        // Mark observations 1 and 3 as rejected
        obs_fit_data[1].selection = ObsSelection::Rejected;
        obs_fit_data[3].selection = ObsSelection::ForcedOut;

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &[true; 6],
            &cache,
            &JPL_EPHEM_HORIZON,
            true,
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        assert_eq!(
            result.updated_obs_fit_data[1].selection,
            ObsSelection::Rejected
        );
        assert_eq!(
            result.updated_obs_fit_data[3].selection,
            ObsSelection::ForcedOut
        );
    }

    /// Active observations must have their residuals updated (non-NaN).
    #[test]
    fn test_active_observations_residuals_are_finite() {
        let t0 = 59000.0;
        let elements = circular_elements(t0);
        let (obs_dataset, cache) = make_dataset_and_cache(t0, 365.0, 6);

        let observations: Vec<Observation> = (0..6)
            .map(|i| obs_dataset.get_observation(i).unwrap().clone())
            .collect();

        let obs_fit_data: Vec<ObsFitData> = observations
            .iter()
            .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
            .collect();

        let result = single_iteration(
            &observations,
            &obs_fit_data,
            &elements,
            &[true; 6],
            &cache,
            &JPL_EPHEM_HORIZON,
            true,
            &PropagatorKind::TwoBody,
        )
        .unwrap();

        for fit in &result.updated_obs_fit_data {
            assert!(fit.residual_ra.is_finite());
            assert!(fit.residual_dec.is_finite());
            assert!(fit.chi.is_finite());
        }
    }
}
