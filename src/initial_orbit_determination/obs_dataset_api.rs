//! IOD entry points that operate on [`ObsDataset`](photom::observation_dataset::ObsDataset).
//!
//! This module exposes the [`crate::FitIOD`] trait, which adds Initial Orbit
//! Determination methods directly to
//! [`photom::observation_dataset::ObsDataset`]. Two entry points are provided:
//!
//! - [`crate::FitIOD::fit_iod`] — run IOD for a **single** named trajectory.
//! - [`crate::FitIOD::fit_full_iod`] — run IOD for **every** trajectory in the
//!   dataset and collect the results in a [`crate::FullOrbitResult`] map.
//!
//! Both methods apply an astrometric error model, batch-RMS corrections, and
//! build a per-observation position cache before invoking the Gauss method.
//!
//! # Type aliases
//!
//! - [`crate::IODRMS`] — scalar quality metric (RMS of normalised residuals).
//! - [`crate::FullOrbitResult`] — batch result map keyed by trajectory ID.

use hifitime::ut1::Ut1Provider;
use photom::{
    observation_dataset::{observation::Observation, ObsDataset},
    observer::error_model::{ModelCorrection, ObsErrorModel},
    TrajId,
};
use rand::{rngs::SmallRng, SeedableRng};
use std::collections::HashMap;

use crate::{
    cache::OutfitCache, constants::FitOrbitResult, trajectory::TrajectoryFit, FullOrbitResult,
    IODParams, JPLEphem, OutfitError,
};

#[cfg(feature = "parallel")]
use rayon::iter::ParallelIterator;

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

    /// Parallel version of [`FitIOD::fit_full_iod`] that processes trajectories concurrently.
    ///
    /// Enabled with the `parallel` feature flag, which also enables the `rayon`
    /// dependency.  Trajectories are processed in parallel using Rayon, but each
    /// trajectory's RNG is still deterministically derived from the caller's `rng`
    /// to ensure reproducibility regardless of processing order or parallelism.
    #[cfg(feature = "parallel")]
    fn fit_full_iod_parallel(
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
    ) -> Result<FitOrbitResult, OutfitError>;
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
    ) -> Result<FitOrbitResult, OutfitError> {
        let (corrected_dataset, cache, _) =
            prepare_iod(self, jpl, ut1_provider, params, error_model, rng)?;

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
        let (corrected_dataset, cache, base_seed) =
            prepare_iod(self, jpl, ut1_provider, params, error_model, rng)?;

        // Default is not implemented for RandomState so a collect cannot be used directly.
        // Instead, we create a new map which initialize correctly the RandomState
        // then for each trajectory, insert the results into it.
        corrected_dataset
            .iter_traj_id()
            .ok_or(OutfitError::NoTrajectoryIndex)
            .map(|iter| {
                let mut map = HashMap::with_hasher(ahash::RandomState::new());
                iter.map(|traj_id| {
                    process_traj(traj_id, &corrected_dataset, &cache, jpl, params, base_seed)
                })
                .for_each(|(k, v)| {
                    map.insert(k, v);
                });
                map
            })
    }

    #[cfg(feature = "parallel")]
    fn fit_full_iod_parallel(
        self,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        params: &IODParams,
        error_model: ObsErrorModel,
        rng: &mut impl rand::Rng,
    ) -> Result<FullOrbitResult, OutfitError> {
        let (corrected_dataset, cache, base_seed) =
            prepare_iod(self, jpl, ut1_provider, params, error_model, rng)?;

        let new_map = || HashMap::with_hasher(ahash::RandomState::new());

        // parallel iteration over the trajectory Ids
        // For each trajectory, we process it and get a map of results for that trajectory.
        // Fold produce a single map for each thread, and then we reduce all the maps into a single one.
        corrected_dataset
            .par_iter_traj_id()
            .ok_or(OutfitError::NoTrajectoryIndex)
            .map(|iter| {
                iter.map(|traj_id| {
                    process_traj(traj_id, &corrected_dataset, &cache, jpl, params, base_seed)
                })
                .fold(new_map, |mut map, (k, v)| {
                    map.insert(k, v);
                    map
                })
                .reduce(new_map, |mut a, b| {
                    a.extend(b);
                    a
                })
            })
    }
}

fn fit_single_traj(
    traj: &TrajId,
    corrected_dataset: &ObsDataset,
    cache: &OutfitCache,
    jpl: &JPLEphem,
    params: &IODParams,
    rng: &mut impl rand::Rng,
) -> Result<FitOrbitResult, OutfitError> {
    let materialized_traj = corrected_dataset
        .materialize_trajectory(traj)
        .ok_or_else(|| OutfitError::TrajectoryIdNotFound(traj.clone()))?;

    let mut obs_vec_refs: Vec<&Observation> = materialized_traj.collect_into_vec();
    obs_vec_refs.sort_by(|a, b| a.mjd_tt().total_cmp(&b.mjd_tt()));

    obs_vec_refs.estimate_best_orbit(cache, jpl, params, rng)
}

/// Run IOD directly on a pre-sorted, pre-corrected slice of observations using
/// an already-built position cache.
///
/// This is the low-level entry point used when the caller has already applied
/// the error model and built the [`OutfitCache`].  It avoids the cost of
/// reconstructing an [`ObsDataset`] and rebuilding the cache.
///
/// # Arguments
///
/// - `observations` — slice of observations, **sorted by MJD** (ascending).
/// - `cache` — position cache already built for these observations.
/// - `jpl` — JPL ephemeris.
/// - `params` — IOD tuning parameters.
/// - `rng` — random-number generator for noise realisations.
pub(crate) fn run_iod_on_observations(
    observations: &[Observation],
    cache: &OutfitCache,
    jpl: &JPLEphem,
    params: &IODParams,
    rng: &mut impl rand::Rng,
) -> Result<FitOrbitResult, OutfitError> {
    let mut refs: Vec<&Observation> = observations.iter().collect();
    refs.sort_by(|a, b| a.mjd_tt().total_cmp(&b.mjd_tt()));
    refs.estimate_best_orbit(cache, jpl, params, rng)
}

fn prepare_iod(
    dataset: ObsDataset,
    jpl: &JPLEphem,
    ut1_provider: &Ut1Provider,
    params: &IODParams,
    error_model: ObsErrorModel,
    rng: &mut impl rand::Rng,
) -> Result<(ObsDataset, OutfitCache, u64), OutfitError> {
    let corrected_dataset = dataset
        .with_error_model(error_model)
        .apply_model_errors()
        .apply_batch_rms_correction(params.gap_max);

    let cache = OutfitCache::build(&corrected_dataset, jpl, ut1_provider, true)?;

    // Draw a single base seed from the caller's RNG.
    // All per-trajectory RNGs are derived from this seed → deterministic
    // regardless of trajectory ordering or parallelism.
    let base_seed: u64 = rng.random();

    Ok((corrected_dataset, cache, base_seed))
}

fn process_traj(
    traj_id: &TrajId,
    corrected_dataset: &ObsDataset,
    cache: &OutfitCache,
    jpl: &JPLEphem,
    params: &IODParams,
    base_seed: u64,
) -> (TrajId, Result<FitOrbitResult, OutfitError>) {
    let traj_seed = base_seed ^ traj_id.stable_hash();
    let mut local_rng = SmallRng::seed_from_u64(traj_seed);
    let result = fit_single_traj(
        traj_id,
        corrected_dataset,
        cache,
        jpl,
        params,
        &mut local_rng,
    );
    (traj_id.clone(), result)
}
