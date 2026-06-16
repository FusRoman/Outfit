use nalgebra::Matrix2;
use photom::observation_dataset::observation::Observation;

pub mod diff_cor;
pub mod least_square;
pub mod obs_dataset_api;
pub mod obs_fit_data;
pub mod outlier_rejection;
pub mod single_iteration;

pub use diff_cor::{
    run_differential_correction, DifferentialCorrectionConfig, DifferentialCorrectionOutput,
};
pub use least_square::{
    angular_diff, rescale_covariance, solve_weighted_least_squares, DifferentialCorrectionResult,
    ObservationEquation, OrbitalUncertainty,
};
pub use obs_dataset_api::FitLSQ;
pub use obs_fit_data::{ObsFitData, ObsSelection};
pub use outlier_rejection::{update_observation_selection, OutlierRejectionConfig};
pub use single_iteration::{single_iteration, SingleIterationResult};

use crate::{
    cache::OutfitCache, constants::FitOrbitResult,
    initial_orbit_determination::obs_dataset_api::run_iod_on_observations, IODParams, JPLEphem,
    OutfitError,
};

/// Returns the 2×2 diagonal weight matrix for a single astrometric observation.
///
/// The weight matrix encodes the per-axis measurement precision:
///
/// ```text
/// W = diag(1/σ_α², 1/σ_δ²)
/// ```
///
/// When `rejected` is `true` a zero matrix is returned, which effectively
/// removes the observation from all subsequent linear-algebra steps (normal
/// matrix accumulation, residual summation, etc.).
///
/// # Arguments
///
/// - `obs` — the observation whose RA/Dec errors define the weights.
/// - `rejected` — pass `true` to exclude the observation from the fit.
///
/// # Returns
///
/// A symmetric 2×2 matrix.  The off-diagonal entries are always `0` (no
/// cross-correlation between RA and Dec is assumed).
pub fn observation_weight(obs: &Observation, rejected: bool) -> Matrix2<f64> {
    if rejected {
        return Matrix2::zeros();
    }

    let sa2 = obs.equ_coord().ra_error.powf(2.0);
    let sd2 = obs.equ_coord().dec_error.powf(2.0);
    Matrix2::new(1.0 / sa2, 0.0, 0.0, 1.0 / sd2)
}

pub fn differential_correction(
    observations: &[Observation],
    cache: &OutfitCache,
    jpl: &JPLEphem,
    iod_params: &IODParams,
    diff_cor_config: &DifferentialCorrectionConfig,
    initial_orbits: Option<&FitOrbitResult>,
    rng: &mut impl rand::Rng,
) -> Result<FitOrbitResult, OutfitError> {
    // ── Obtain the starting equinoctial elements ──────────────────────────────
    let initial_orbit: FitOrbitResult = match initial_orbits {
        Some(FitOrbitResult::DifferentialCorrection(_)) => {
            return Err(OutfitError::DifferentialCorrectionFailed(
                "differential_correction expects an IOD result as initial orbit, \
             got a DifferentialCorrection result"
                    .to_string(),
            ));
        }
        Some(orbit) => orbit.clone(),
        None => run_iod_on_observations(observations, cache, jpl, iod_params, rng)?,
    };

    let equinoctial_elements = initial_orbit
        .orbital_elements()
        .to_equinoctial()?
        .as_equinoctial()
        .ok_or(OutfitError::InvalidConversion(
            "Conversion to equinoctial elements failed".to_string(),
        ))?;

    // ── Build per-observation fit data from the error-model uncertainties ─────
    let obs_fit_data: Vec<ObsFitData> = observations
        .iter()
        .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
        .collect();

    // ── Run the differential-correction loop ──────────────────────────────────
    match run_differential_correction(
        observations,
        &obs_fit_data,
        &equinoctial_elements,
        cache,
        jpl,
        diff_cor_config,
    ) {
        Ok(dc_output) => {
            // ── Package the result ────────────────────────────────────────────────────
            let normalised_rms = dc_output.normalised_rms;
            Ok(FitOrbitResult::DifferentialCorrection((
                dc_output.into(),
                normalised_rms,
            )))
        }
        Err(_) => Ok(initial_orbit),
    }
}
