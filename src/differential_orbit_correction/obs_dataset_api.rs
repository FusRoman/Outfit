//! Least-squares orbit fitting entry points for [`ObsDataset`].
//!
//! This module exposes the [`FitLSQ`] trait, which adds differential orbit
//! correction methods directly to
//! [`photom::observation_dataset::ObsDataset`].
//!
//! The main entry point is [`FitLSQ::fit_lsq`], which runs the full
//! differential-correction pipeline for every trajectory in the dataset.
//! When no initial orbit is provided for a trajectory, a preliminary orbit is
//! first derived via the Gauss IOD method and then fed into the
//! differential-correction loop.
//!
//! # Type aliases
//!
//! - [`crate::constants::Chi2`] — scalar quality metric (normalised RMS after
//!   the final iteration).
//! - [`crate::FullOrbitResult`] — batch result map keyed by trajectory ID.

use std::collections::HashMap;

use crate::{
    cache::OutfitCache,
    constants::FitOrbitResult,
    differential_orbit_correction::{
        diff_cor::{run_differential_correction, DifferentialCorrectionConfig},
        obs_fit_data::ObsFitData,
    },
    initial_orbit_determination::{obs_dataset_api::run_iod_on_observations, IODParams},
    FullOrbitResult, JPLEphem, OrbitalElements, OutfitError,
};
use hifitime::ut1::Ut1Provider;
use photom::{
    observation_dataset::{observation::Observation, ObsDataset},
    observer::error_model::{ModelCorrection, ObsErrorModel},
    TrajId,
};
use rand::{rngs::SmallRng, SeedableRng};

// ─────────────────────────────────────────────────────────────────────────────
// Trait definition
// ─────────────────────────────────────────────────────────────────────────────

/// Extension trait that adds differential orbit correction to [`ObsDataset`].
///
/// Import this trait to call [`fit_lsq`](FitLSQ::fit_lsq) on any
/// [`ObsDataset`] value.
pub trait FitLSQ {
    /// Run the differential-correction pipeline for **every** trajectory in
    /// the dataset.
    ///
    /// For each trajectory:
    ///
    /// 1. If `initial_orbits` contains an entry for that trajectory, it is
    ///    converted to equinoctial elements and used as the starting point.
    /// 2. Otherwise a preliminary orbit is first derived via the Gauss IOD
    ///    method (using `iod_params` and `iod_error_model`) and the result is
    ///    used as the starting point.
    /// 3. The differential-correction loop is run using `diff_cor_config`.
    ///
    /// # Arguments
    ///
    /// - `jpl` — JPL planetary ephemeris.
    /// - `ut1_provider` — UT1 time-scale data for Earth orientation.
    /// - `error_model` — astrometric error model applied to every observation
    ///   before the fit.
    /// - `iod_params` — IOD tuning parameters used when no initial orbit is
    ///   provided for a trajectory.
    /// - `diff_cor_config` — differential-correction tuning parameters.
    /// - `initial_orbits` — optional map from trajectory ID to a known initial
    ///   orbit.  Trajectories absent from this map (or when the map is `None`)
    ///   will be initialised via IOD.
    /// - `rng` — source of randomness for IOD noise sampling.
    ///
    /// # Quality metric — normalised RMS
    ///
    /// Each successful entry in the returned map carries a scalar quality
    /// metric (`Chi2`, accessible as the second element of
    /// [`FitOrbitResult::DifferentialCorrection`]).  This is the **normalised
    /// RMS** of the final fit, defined as:
    ///
    /// \\[ \text{normalised\_rms} = \sqrt{\frac{\xi^\top W \xi}{n_{\text{active}}}} \\]
    ///
    /// where **ξ** is the vector of residuals (observed minus computed
    /// positions in RA and Dec), **W** the diagonal weight matrix (inverse
    /// of the per-observation variances), and *n*_active the total number of
    /// active scalar measurements (2 per optical observation).
    ///
    /// This is the **square root of the reduced chi²** (chi² per degree of
    /// freedom, with the number of degrees of freedom approximated by
    /// *n*_active):
    ///
    /// - **≈ 1.0** — residuals are consistent with the reported observation
    ///   uncertainties; the fit is statistically good.
    /// - **> 1.0** — residuals exceed the expected noise; the orbit does not
    ///   fit well, or the uncertainties are under-estimated.
    /// - **< 1.0** — residuals are smaller than the noise; the uncertainties
    ///   may be over-estimated, or the fit is over-constrained.
    ///
    /// Unlike the raw chi² (which grows with the number of observations), the
    /// normalised RMS is directly comparable across trajectories with
    /// different numbers of observations.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if the shared observer cache cannot be built.
    /// Individual trajectory failures are stored as `Err(…)` entries inside
    /// the returned [`FullOrbitResult`].
    #[allow(clippy::too_many_arguments)]
    fn fit_lsq(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        error_model: ObsErrorModel,
        iod_params: &IODParams,
        diff_cor_config: &DifferentialCorrectionConfig,
        initial_orbits: Option<&FullOrbitResult>,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError>;
}

// ─────────────────────────────────────────────────────────────────────────────
// Implementation
// ─────────────────────────────────────────────────────────────────────────────

impl FitLSQ for ObsDataset {
    fn fit_lsq(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        error_model: ObsErrorModel,
        iod_params: &IODParams,
        diff_cor_config: &DifferentialCorrectionConfig,
        initial_orbits: Option<&FullOrbitResult>,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError> {
        // ── 1. Apply the error model and build the shared position cache ──────
        let corrected_dataset = self
            .with_error_model(error_model)
            .apply_model_errors()
            .apply_batch_rms_correction(iod_params.gap_max);

        let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider, true)?;

        // Draw a base seed for deterministic per-trajectory RNGs.
        let base_seed: u64 = rng.random();

        // ── 2. Iterate over trajectories ──────────────────────────────────────
        let mut result_map = HashMap::with_hasher(ahash::RandomState::new());

        let traj_ids: Vec<TrajId> = corrected_dataset
            .iter_traj_id()
            .ok_or(OutfitError::NoTrajectoryIndex)?
            .cloned()
            .collect();

        for traj_id in &traj_ids {
            let traj_seed = base_seed ^ traj_id.stable_hash();
            let mut local_rng = SmallRng::seed_from_u64(traj_seed);

            let outcome = run_differential_correction_for_trajectory(
                traj_id,
                &corrected_dataset,
                &cache,
                jpl,
                ut1_provider,
                iod_params,
                diff_cor_config,
                initial_orbits,
                &mut local_rng,
            );

            result_map.insert(traj_id.clone(), outcome);
        }

        Ok(result_map)
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Internal helpers
// ─────────────────────────────────────────────────────────────────────────────

/// Runs the full pipeline (IOD if needed + differential correction) for a
/// single trajectory.
///
/// Returns `Ok(FitOrbitResult::DifferentialCorrection(…))` on success or an
/// [`OutfitError`] on any failure (IOD failure, bizarre orbit, divergence, …).
#[allow(clippy::too_many_arguments)]
fn run_differential_correction_for_trajectory(
    traj_id: &TrajId,
    dataset: &ObsDataset,
    cache: &OutfitCache,
    jpl: &JPLEphem,
    _ut1_provider: &Ut1Provider,
    iod_params: &IODParams,
    diff_cor_config: &DifferentialCorrectionConfig,
    initial_orbits: Option<&FullOrbitResult>,
    rng: &mut impl rand::Rng,
) -> Result<FitOrbitResult, OutfitError> {
    // ── Collect and sort observations for this trajectory ─────────────────────
    let materialized = dataset
        .materialize_trajectory(traj_id)
        .ok_or_else(|| OutfitError::TrajectoryIdNotFound(traj_id.clone()))?;

    let mut obs_refs: Vec<&Observation> = materialized.collect_into_vec();
    obs_refs.sort_by(|a, b| a.mjd_tt().total_cmp(&b.mjd_tt()));
    let observations: Vec<Observation> = obs_refs.into_iter().cloned().collect();

    // ── Obtain the starting equinoctial elements ──────────────────────────────
    let initial_equinoctial = match initial_orbits.and_then(|map| map.get(traj_id)) {
        Some(Ok(orbital_elements)) => {
            // Caller provided an initial orbit — convert to equinoctial.
            orbital_elements.orbital_elements().to_equinoctial()?
        }
        Some(Err(_)) | None => {
            // No initial orbit — run IOD directly on the already-corrected
            // observations, reusing the cache that was built for the full
            // dataset.  This avoids reconstructing an ObsDataset and
            // rebuilding the cache.
            let iod_result = run_iod_on_observations(&observations, cache, jpl, iod_params, rng)?;
            iod_result.orbital_elements().to_equinoctial()?
        }
    };

    println!("\n\nTrajectory {traj_id}: initial equinoctial elements: {initial_equinoctial:?}\n\n");

    // ── Build per-observation fit data from the error-model uncertainties ─────
    let obs_fit_data: Vec<ObsFitData> = observations
        .iter()
        .map(|obs| ObsFitData::new(obs.equ_coord().ra_error, obs.equ_coord().dec_error))
        .collect();

    // ── Run the differential-correction loop ──────────────────────────────────
    let dc_output = run_differential_correction(
        &observations,
        &obs_fit_data,
        &initial_equinoctial,
        cache,
        jpl,
        diff_cor_config,
    )?;

    println!("\ndc_output for trajectory {traj_id}: {dc_output:?}\n");

    // ── Package the result ────────────────────────────────────────────────────
    let orbital_elements = OrbitalElements::Equinoctial(dc_output.elements);
    Ok(FitOrbitResult::DifferentialCorrection((
        orbital_elements,
        dc_output.normalised_rms,
    )))
}
