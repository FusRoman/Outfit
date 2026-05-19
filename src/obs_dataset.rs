//! IOD entry points that operate on [`ObsDataset`](photom::observation_dataset::ObsDataset).
//!
//! This module exposes the [`FitIOD`] trait, which adds Initial Orbit
//! Determination methods directly to
//! [`photom::observation_dataset::ObsDataset`]. Two entry points are provided:
//!
//! - [`FitIOD::fit_iod`] — run IOD for a **single** named trajectory.
//! - [`FitIOD::fit_full_iod`] — run IOD for **every** trajectory in the
//!   dataset and collect the results in a [`FullOrbitResult`] map.
//!
//! Both methods apply an astrometric error model, batch-RMS corrections, and
//! build a per-observation position cache before invoking the Gauss method.
//!
//! # Type aliases
//!
//! - [`IODRMS`] — scalar quality metric (RMS of normalised residuals).
//! - [`FullOrbitResult`] — batch result map keyed by trajectory ID.

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

/// Extension trait that adds Initial Orbit Determination methods to
/// [`ObsDataset`].
///
/// Import this trait to call [`fit_iod`](FitIOD::fit_iod) or
/// [`fit_full_iod`](FitIOD::fit_full_iod) on any [`ObsDataset`] value.
pub trait FitIOD {
    /// Run the Gauss IOD pipeline for **every** trajectory in the dataset.
    ///
    /// Applies the given astrometric `error_model`, builds a shared
    /// per-observation position cache, and then attempts a preliminary orbit
    /// determination for each trajectory independently.  Each trajectory gets
    /// its own deterministic random seed derived from `rng`, so results are
    /// reproducible regardless of the order in which trajectories are
    /// processed.
    ///
    /// # Arguments
    ///
    /// - `jpl` — JPL planetary ephemeris used for heliocentric Earth positions
    ///   and light-time corrections.
    /// - `ut1_provider` — UT1 time-scale data for Earth orientation corrections.
    /// - `params` — IOD tuning parameters (triplet selection, noise
    ///   realizations, RMS window, …).
    /// - `error_model` — astrometric error model assigned to every observation
    ///   before the fit.
    /// - `rng` — source of randomness for Monte-Carlo noise sampling; a single
    ///   seed is drawn here and then per-trajectory seeds are derived from it.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if the shared cache cannot be built (e.g., an
    /// observation has no associated observer ID or the JPL ephemeris is out of
    /// range).  Individual trajectory failures are **not** propagated as errors;
    /// they are stored as `Err(…)` entries in the returned [`FullOrbitResult`].
    fn fit_full_iod(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError>;

    /// Run the Gauss IOD pipeline for a **single** named trajectory.
    ///
    /// Filters the dataset to the observations belonging to `traj`, applies
    /// `error_model`, builds a position cache, and then searches for the
    /// best-fitting preliminary orbit via the Gauss method with Monte-Carlo
    /// noise sampling.
    ///
    /// # Arguments
    ///
    /// - `traj` — trajectory identifier; anything that converts into
    ///   [`photom::TrajId`] (e.g., a `&str` or `String`).
    /// - `jpl` — JPL planetary ephemeris.
    /// - `ut1_provider` — UT1 time-scale data.
    /// - `params` — IOD tuning parameters.
    /// - `error_model` — astrometric error model.
    /// - `rng` — source of randomness for noise sampling.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if:
    /// - `traj` is not found in the dataset ([`OutfitError::TrajectoryIdNotFound`]),
    /// - the position cache cannot be built, or
    /// - no viable orbit is found after exhausting all triplet candidates
    ///   ([`OutfitError::NoViableOrbit`]).
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
                    traj_id,
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
