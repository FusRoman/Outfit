use ahash::AHashMap;
use hifitime::ut1::Ut1Provider;
use photom::{
    observation_dataset::{observation::Observation, ObsDataset},
    observer::error_model::ModelCorrection,
    TrajId,
};

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
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError>;

    fn fit_iod(
        self,
        traj: impl Into<TrajId>,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
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
        rng: &mut impl rand::Rng,
    ) -> Result<(GaussResult, IODRMS), OutfitError> {
        let corrected_dataset = self
            .apply_model_errors()
            .apply_batch_rms_correction(params.gap_max);

        let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider)?;

        fit_single_traj(&traj.into(), &corrected_dataset, &cache, jpl, params, rng)
    }

    fn fit_full_iod(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError> {
        let corrected_dataset = self
            .apply_model_errors()
            .apply_batch_rms_correction(params.gap_max);

        let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider)?;

        let results: FullOrbitResult = corrected_dataset
            .iter_traj_id()
            .into_iter()
            .flatten()
            .map(|traj| {
                println!("Fitting IOD for Trajectory ID: {}", traj);
                let result = fit_single_traj(traj, &corrected_dataset, &cache, jpl, params, rng);
                (traj.clone(), result)
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

    let mut obs_vec_refs: Vec<&Observation> = materialized_traj.iter().collect();
    obs_vec_refs.sort_by(|a, b| a.mjd_tt().total_cmp(&b.mjd_tt()));

    obs_vec_refs.estimate_best_orbit(cache, jpl, params, rng)
}
