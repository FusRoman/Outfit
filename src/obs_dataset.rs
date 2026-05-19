use ahash::AHashMap;
use hifitime::ut1::Ut1Provider;
use photom::{
    observation_dataset::{observation::Observation, ObsDataset},
    observer::error_model::{ModelCorrection, ObsErrorModel},
    TrajId,
};
use rand::{rngs::SmallRng, SeedableRng};

use crate::{
    cache::OutfitCache, trajectory::TrajectoryFit, GaussResult, IODParams, JPLEphem, OutfitError,
};

/// Type alias for the RMS of normalized residuals from an IOD fit.
/// This is a single scalar value representing the overall fit quality of the IOD solution.
pub type IODRMS = f64;

/// Full batch orbit determination results.
///
/// Each entry maps an [`TrajId`] to the outcome of a full
/// Initial Orbit Determination (IOD) attempt on its set of observations.
///
/// Internally, this is implemented as:
///
/// ```ignore
/// HashMap<TrajId, Result<(GaussResult, IODRMS), OutfitError>, RandomState>
/// ```
///
/// Return semantics
/// -----------------
/// * `Ok((GaussResult, IODRMS))` – a successful IOD with its RMS of normalized residuals.
/// * `Err(OutfitError)` – a failure isolated to that object.
pub type FullOrbitResult = AHashMap<TrajId, Result<(GaussResult, IODRMS), OutfitError>>;

pub trait FitIOD {
    fn fit_full_iod(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError>;

    fn fit_iod(
        self,
        traj: impl Into<TrajId>,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<(GaussResult, IODRMS), OutfitError>;
}

impl FitIOD for ObsDataset {
    fn fit_iod(
        self,
        traj: impl Into<TrajId>,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<(GaussResult, IODRMS), OutfitError> {
        let corrected_dataset = self
            .with_error_model(error_model)
            .apply_batch_rms_correction(params.gap_max);

        let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider)?;

        fit_single_traj(&traj.into(), &corrected_dataset, &cache, jpl, params, rng)
    }

    fn fit_full_iod(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError> {
        let corrected_dataset = self
            .with_error_model(error_model)
            .apply_model_errors()
            .apply_batch_rms_correction(params.gap_max);

        let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider)?;

        // Draw a single base seed from the caller's RNG.
        // All per-trajectory RNGs are derived from this seed → deterministic
        // regardless of trajectory ordering or parallelism.
        let base_seed: u64 = rng.random();

        let results: FullOrbitResult = corrected_dataset
            .iter_traj_id()
            .into_iter()
            .flatten()
            .map(|traj_id| {
                // Derive a per-trajectory seed from the base seed and the trajectory ID.
                // Two trajectories with different IDs always get independent sequences.
                let traj_seed = base_seed ^ traj_id.stable_hash();
                let mut local_rng = SmallRng::seed_from_u64(traj_seed);

                let result = fit_single_traj(
                    &traj_id,
                    &corrected_dataset,
                    &cache,
                    jpl,
                    params,
                    &mut local_rng,
                );
                (traj_id.clone(), result)
            })
            .collect();

        Ok(results)
    }
}

fn fit_single_traj(
    traj: &TrajId,
    corrected_dataset: &ObsDataset,
    cache: &OutfitCache,
    jpl: &JPLEphem,
    params: &IODParams,
    rng: &mut impl rand::Rng,
) -> Result<(GaussResult, IODRMS), OutfitError> {
    let materialized_traj = corrected_dataset
        .materialize_trajectory(traj)
        .ok_or_else(|| OutfitError::TrajectoryIdNotFound(traj.clone()))?;

    let mut obs_vec_refs: Vec<&Observation> = materialized_traj.collect_into_vec();
    obs_vec_refs.sort_by(|a, b| a.mjd_tt().total_cmp(&b.mjd_tt()));

    obs_vec_refs.estimate_best_orbit(cache, jpl, params, rng)
}
