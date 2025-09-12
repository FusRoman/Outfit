//! # Batch Gauss IOD over Trajectory Sets
//!
//! Run a full **Gauss-based Initial Orbit Determination (IOD)** over a
//! [`TrajectorySet`], collect **per-object outcomes**, and expose convenience
//! helpers to query or extract solutions and summarize observation counts.
//!
//! ## Overview
//! -----------------
//! A [`TrajectorySet`] maps each [`ObjectNumber`] to its time-ordered
//! [`Observations`]. This module implements the [`TrajectoryFit`] trait on
//! `TrajectorySet`, providing:
//!
//! * `estimate_all_orbits` – run the Gauss IOD pipeline on **every object**,
//! * `estimate_all_orbits_with_cancel` – same, with **cooperative cancellation**,
//! * `total_observations` / `number_of_trajectories` – quick set-level metrics,
//! * `obs_count_stats` – summary statistics on observation counts,
//! * `gauss_result_for` / `take_gauss_result` – ergonomic access to results.
//!
//! All objects are processed with the same [`Outfit`] state (ephemerides, error
//! model, frames), a caller-provided RNG, and a single [`IODParams`] configuration.
//!
//! ## Result Model
//! -----------------
//! Batch outcomes are returned as a [`FullOrbitResult`]:
//!
//! ```text
//! ObjectNumber → Result<(GaussResult, rms: f64), OutfitError>
//! ```
//!
//! * `Ok((GaussResult, rms))` – the best preliminary/corrected orbit and its RMS
//!   of normalized astrometric residuals,
//! * `Err(OutfitError)` – a failure **isolated** to that object (other objects
//!   continue to be processed).
//!
//! Use [`gauss_result_for`] to **borrow** a solution and its RMS, or
//! [`take_gauss_result`] to **move** them out of the map.
//!
//! ## Execution Modes
//! -----------------
//! ### Progress UI (feature: `progress`)
//! When compiled with the `progress` feature, `estimate_all_orbits` renders a
//! live progress bar (via `indicatif`) and reports per-iteration timing via
//! a lightweight moving average to help diagnose throughput bottlenecks.
//!
//! ### Cooperative cancellation
//! `estimate_all_orbits_with_cancel` periodically calls a user-provided
//! closure `should_cancel()` based on **wall-clock intervals** (not iteration
//! counts) to keep cancellation latency stable even if some objects are slow.
//!
//! ## Performance Notes
//! -----------------
//! * The loop walks the underlying map once; overall time scales with the number
//!   of objects × the cost of `ObservationIOD::estimate_best_orbit` (triplet
//!   enumeration, scoring, optional correction).
//! * Results are accumulated in a `HashMap` that uses `ahash::RandomState`,
//!   matching the default hasher used elsewhere in the crate.
//! * No mutation of the observations themselves; only per-object IOD is performed.
//!
//! ## Error Semantics
//! -----------------
//! * Failures are **per-object**: an error for one object does **not** abort
//!   the batch.
//! * The returned map always contains **one entry per processed object**,
//!   each entry being either `Ok((GaussResult, rms))` or `Err(OutfitError)`.
//!
//! ## Examples
//! -----------------
//! Minimal end-to-end run (no progress UI):
//!
//! ```rust,no_run
//! use rand::SeedableRng;
//! use outfit::{Outfit, TrajectorySet};
//! use outfit::initial_orbit_determination::IODParams;
//! use outfit::trajectories::trajectory_fit::TrajectoryFit;
//!
//! # fn demo(mut trajs: TrajectorySet) -> Result<(), outfit::outfit_errors::OutfitError> {
//! let state = Outfit::new("horizon:DE440", outfit::error_models::ErrorModel::FCCT14)?;
//! let mut rng = rand::rngs::StdRng::from_os_rng();
//! let params = IODParams::builder().max_triplets(32).build()?;
//!
//! let results = trajs.estimate_all_orbits(&state, &mut rng, &params);
//! for (obj, res) in &results {
//!     match res {
//!         Ok((g, rms)) => eprintln!("{obj:?}: orbit={} rms={:.4}", g, rms),
//!         Err(e)       => eprintln!("{obj:?}: error={e}"),
//!     }
//! }
//! # Ok(()) }
//! ```
//!
//! Cooperative cancellation (poll every ~20 ms):
//!
//! ```rust,no_run
//! use std::sync::atomic::{AtomicBool, Ordering};
//! use outfit::trajectories::trajectory_fit::TrajectoryFit;
//!
//! # fn cancelable(mut trajs: outfit::TrajectorySet,
//! #                state: &outfit::Outfit,
//! #                rng: &mut impl rand::Rng,
//! #                params: &outfit::IODParams,
//! # ) -> outfit::trajectories::trajectory_fit::FullOrbitResult {
//! let stop = AtomicBool::new(false);
//! // … flip `stop` from another thread / signal handler …
////!
//! trajs.estimate_all_orbits_with_cancel(state, rng, params, || stop.load(Ordering::Relaxed))
//! # }
//! ```
//!
//! Quick stats for logging/reporting:
//!
//! ```rust
//! use outfit::trajectories::trajectory_fit::TrajectoryFit;
//!
//! fn summarize(set: &outfit::TrajectorySet) {
//!     let n_obj = set.number_of_trajectories();
//!     let n_obs = set.total_observations();
//!     if let Some(stats) = set.obs_count_stats() {
//!         eprintln!("Trajectories: {n_obj}, Observations: {n_obs}");
//!         eprintln!("{:#}", stats);
//!     }
//! }
//! ```
//!
//! ## See also
//! ------------
//! * [`TrajectoryFit`] – Trait implemented by `TrajectorySet` for batch IOD and stats.
//! * [`ObservationIOD::estimate_best_orbit`] – Per-object Gauss IOD (called internally).
//! * [`GaussResult`] – Preliminary/corrected orbit container.
//! * [`IODParams`] – Tuning for triplet generation, scoring, correction.
//! * [`gauss_result_for`] / [`take_gauss_result`] – Accessors for the `FullOrbitResult` map.
use std::{collections::HashMap, fmt};

use ahash::RandomState;
use rand::Rng;

use crate::{
    constants::Observations, GaussResult, IODParams, ObjectNumber, ObservationIOD, Outfit,
    OutfitError, TrajectorySet,
};

use std::time::{Duration, Instant};

#[cfg(feature = "progress")]
use super::progress_bar::IterTimer;
#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};

/// Full batch orbit determination results.
///
/// Each entry maps an [`ObjectNumber`] to the outcome of a full
/// Initial Orbit Determination (IOD) attempt on its set of observations.
/// The result type is:
///
/// * `Ok((Option<GaussResult>, f64))` – successful execution of the IOD pipeline:
///   * `Option<GaussResult>` – the best preliminary or corrected orbit found
///     for this object (or `None` if no viable orbit was selected),
///   * `f64` – the RMS of normalized astrometric residuals for that solution.
/// * `Err(OutfitError)` – a failure during orbit estimation for this object.
///   Errors are per-object and do not abort the rest of the batch.
///
/// Internally, this is implemented as:
///
/// ```ignore
/// HashMap<ObjectNumber, Result<(Option<GaussResult>, f64), OutfitError>, RandomState>
/// ```
pub type FullOrbitResult =
    HashMap<ObjectNumber, Result<(GaussResult, f64), OutfitError>, RandomState>;

/// Borrow a Gauss solution (if any) and its RMS for a given key.
///
/// Arguments
/// -----------------
/// * `all`: The map of all IOD outcomes.
/// * `key`: The object identifier.
///
/// Return
/// ----------
/// * `Ok(Some((&GaussResult, f64)))` – a solution is present for the key.
/// * `Ok(None)` – key absent OR present but no acceptable solution (`None`).
/// * `Err(&OutfitError)` – the IOD attempt failed for that key.
///
/// See also
/// ------------
/// * [`GaussResult`] – Gauss IOD output structure.
pub fn gauss_result_for<'a>(
    all: &'a FullOrbitResult,
    key: &ObjectNumber,
) -> Result<Option<(&'a GaussResult, f64)>, &'a OutfitError> {
    match all.get(key) {
        None => Ok(None),
        Some(Err(e)) => Err(e),
        Some(Ok((g, rms))) => Ok(Some((g, *rms))),
    }
}

/// Take ownership of the solution for `key`, removing it from the map.
pub fn take_gauss_result(
    all: &mut FullOrbitResult,
    key: &ObjectNumber,
) -> Result<Option<(GaussResult, f64)>, OutfitError> {
    match all.remove(key) {
        None => Ok(None),
        Some(Err(e)) => Err(e),
        Some(Ok((g, rms))) => Ok(Some((g, rms))),
    }
}

/// Summary statistics for per-trajectory observation counts.
///
/// Each [`TrajectorySet`] entry (one object) has an associated
/// [`Observations`] container. This structure stores basic distribution
/// statistics on the **number of observations per trajectory**, as
/// returned by [`obs_count_stats`](crate::trajectories::trajectory_fit::TrajectoryFit::obs_count_stats).
///
/// Fields
/// -----------------
/// * `min` – smallest number of observations in any trajectory.
/// * `p25` – 25th percentile (first quartile) of observation counts.
/// * `median` – 50th percentile (second quartile).
/// * `p95` – 95th percentile, indicating the upper tail of the distribution.
/// * `max` – largest number of observations in any trajectory.
///
/// Percentiles are computed using the *nearest-rank* method:
/// the index is `round(q × (N-1))` for quantile `q ∈ [0,1]`, clamped to valid range.
/// This convention makes results stable even for small sample sizes.
///
/// Display
/// -----------------
/// * `format!("{}", stats)` – compact single-line summary, e.g.:
///   ```text
///   min=2, p25=4, median=8, p95=15, max=20
///   ```
///
/// * `format!("{:#}", stats)` – pretty multi-line table, e.g.:
///   ```text
///   Observation count per trajectory — summary
///   -----------------------------------------
///   min    : 2
///   p25    : 4
///   median : 8
///   p95    : 15
///   max    : 20
///   ```
///
/// See also
/// ------------
/// * [`obs_count_stats`](crate::trajectories::trajectory_fit::TrajectoryFit::obs_count_stats) – Computes these statistics from a [`TrajectorySet`].
#[derive(Debug, Clone, Copy)]
pub struct ObsCountStats {
    pub min: usize,
    pub p25: usize,
    pub median: usize,
    pub p95: usize,
    pub max: usize,
}

impl fmt::Display for ObsCountStats {
    /// Compact by default; pretty multi-line when using the alternate flag (`{:#}`).
    ///
    /// See also
    /// ------------
    /// * [`obs_count_stats`](crate::trajectories::trajectory_fit::TrajectoryFit::obs_count_stats) – Builder of these summary statistics.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            // Pretty, multi-line, aligned output (ASCII-only for portability).
            writeln!(f, "Observation count per trajectory — summary")?;
            writeln!(f, "-----------------------------------------")?;
            writeln!(f, "min    : {}", self.min)?;
            writeln!(f, "p25    : {}", self.p25)?;
            writeln!(f, "median : {}", self.median)?;
            writeln!(f, "p95    : {}", self.p95)?;
            write!(f, "max    : {}", self.max)
        } else {
            // Compact single-line for logs and quick prints.
            write!(
                f,
                "min={}, p25={}, median={}, p95={}, max={}",
                self.min, self.p25, self.median, self.p95, self.max
            )
        }
    }
}

pub trait TrajectoryFit {
    /// Estimate an initial orbit for **every trajectory** in the set and collect the results.
    ///
    /// This method runs the Gauss-based IOD pipeline on each `(ObjectNumber → Observations)`
    /// pair contained in the trajectory set. All objects are processed with the same
    /// configuration (`state`, `error_model`, `params`) and random number generator.
    /// Results are aggregated into a [`FullOrbitResult`] map.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: The global environment providing ephemerides, constants, and reference frames.
    /// * `rng`: Random number generator used for noisy triplet realizations (e.g., \[`StdRng`\]).
    /// * `params`: Parameters controlling triplet generation, scoring, and correction loops.
    ///
    /// Return
    /// ----------
    /// * A [`FullOrbitResult`] mapping each object to either:
    ///   * `Ok((Option<GaussResult>, f64))` – a valid solution with its RMS,
    ///   * `Err(OutfitError)` – an error diagnostic for that object.
    ///
    /// Notes
    /// ----------
    /// * The method iterates in place on the current set.  
    /// * Observations themselves are not modified, but computation time scales
    ///   with the number of trajectories and candidate triplets.  
    /// * Errors are isolated: one object failing does not prevent others from being processed.
    ///
    /// See also
    /// ------------
    /// * [`ObservationIOD::estimate_best_orbit`] – Per-trajectory IOD with best-orbit selection.
    /// * [`GaussResult`] – Variants for preliminary or corrected orbit solutions.
    /// * [`IODParams`] – Tuning parameters for IOD batch execution.
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult;

    /// Count the total number of [`Observation`](crate::observations::Observation) entries across all trajectories.
    ///
    /// This method iterates once over all values in the [`TrajectorySet`],
    /// summing the length of each [`Observations`] container.
    ///
    /// Return
    /// ----------
    /// * The total number of observations across all objects.
    fn total_observations(&self) -> usize;

    /// Compute distribution statistics for the number of observations per trajectory.
    ///
    /// Each trajectory (one object in the [`TrajectorySet`]) has an associated
    /// [`Observations`] container. This function collects their sizes and computes:
    ///
    /// * `min` – smallest number of observations in any trajectory,
    /// * `p25` – 25th percentile (first quartile),
    /// * `median` – 50th percentile (second quartile),
    /// * `p95` – 95th percentile (upper tail indicator),
    /// * `max` – largest number of observations in any trajectory.
    ///
    /// Percentiles are computed using the *nearest-rank* method:
    /// the index is `round(q × (N-1))` for quantile `q ∈ [0,1]`, clamped to valid range.
    /// This makes results robust even for small datasets.
    ///
    /// Return
    /// ----------
    /// * `None` if the set is empty.
    /// * `Some(ObsCountStats)` containing the summary statistics otherwise.
    ///
    /// See also
    /// ------------
    /// * [`total_observations`](crate::trajectories::trajectory_fit::TrajectoryFit::total_observations) – Sum of all observations across trajectories.
    fn obs_count_stats(&self) -> Option<ObsCountStats>;

    /// Return the number of distinct trajectories (objects) in the set.
    ////
    /// This is simply the number of keys in the underlying map.
    ///
    /// Return
    /// ------
    /// * The number of distinct trajectories (objects) in the set.
    ///
    /// See also
    /// ------------
    /// * [`total_observations`](crate::trajectories::trajectory_fit::TrajectoryFit::total_observations) – Sum of all observations across trajectories.
    /// * [`obs_count_stats`](crate::trajectories::trajectory_fit::TrajectoryFit::obs_count_stats) – Statistics on the number of observations per trajectory.
    fn number_of_trajectories(&self) -> usize;

    fn estimate_all_orbits_with_cancel<F>(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
        should_cancel: F,
    ) -> FullOrbitResult
    where
        F: FnMut() -> bool;
}

impl TrajectoryFit for TrajectorySet {
    #[cfg(feature = "progress")]
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult {
        let total = self.len() as u64;
        let pb = ProgressBar::new(total.max(1));
        pb.set_style(
            ProgressStyle::with_template(
                "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) \
             | {per_sec} | ETA {eta_precise} | {msg}",
            )
            .expect("indicatif template"),
        );
        pb.enable_steady_tick(Duration::from_millis(200));

        let mut results: FullOrbitResult = HashMap::default();
        let mut it_timer = IterTimer::new(0.2);

        for (obj, observations) in self.iter_mut() {
            // Time of the previous loop (zero for the first, acceptable).

            use super::progress_bar::fmt_dur;
            let last = it_timer.tick();
            let avg = it_timer.avg();
            pb.set_message(format!("last: {}, avg: {}", fmt_dur(last), fmt_dur(avg)));

            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);

            pb.inc(1);
        }

        pb.finish_and_clear();
        results
    }

    /// Cooperative cancellation version: the loop periodically calls `should_cancel()`
    /// based on a wall-clock timer (not on iteration count).
    #[cfg(feature = "progress")]
    fn estimate_all_orbits_with_cancel<F>(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
        mut should_cancel: F,
    ) -> FullOrbitResult
    where
        F: FnMut() -> bool,
    {
        let total = self.len() as u64;
        let pb = ProgressBar::new(total.max(1));
        pb.set_style(
        ProgressStyle::with_template(
            "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) | {per_sec} | ETA {eta_precise} | {msg}",
        )
        .expect("indicatif template"),
    );
        pb.enable_steady_tick(Duration::from_millis(200));

        let mut results: FullOrbitResult = HashMap::default();
        let mut it_timer = IterTimer::new(0.2);

        // --- Timer-based polling configuration
        // Keep the cancellation latency roughly constant regardless of iteration cost.
        const POLL_INTERVAL: Duration = Duration::from_millis(20);
        let mut last_poll = Instant::now();

        for (obj, observations) in self.iter_mut() {
            // --- Cancellation poll based on elapsed time
            // Only acquire the GIL / call should_cancel() if enough time has elapsed.
            if last_poll.elapsed() >= POLL_INTERVAL {
                if should_cancel() {
                    pb.set_message("Interrupted");
                    pb.disable_steady_tick();
                    pb.finish_and_clear();
                    break; // or early-return
                }
                last_poll = Instant::now();
            }

            // Progress message (kept as-is)
            use super::progress_bar::fmt_dur;
            let last = it_timer.tick();
            let avg = it_timer.avg();
            pb.set_message(format!("last: {}, avg: {}", fmt_dur(last), fmt_dur(avg)));

            // Core work
            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);

            pb.inc(1);
        }

        pb.disable_steady_tick();
        pb.finish_and_clear();
        results
    }

    #[cfg(not(feature = "progress"))]
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult {
        // Output map using the same fast hasher as TrajectorySet.
        let mut results: FullOrbitResult = HashMap::default();

        for (obj, observations) in self.iter_mut() {
            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);
        }

        results
    }

    #[cfg(not(feature = "progress"))]
    fn estimate_all_orbits_with_cancel<F>(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
        mut should_cancel: F,
    ) -> FullOrbitResult
    where
        F: FnMut() -> bool,
    {
        let mut results: FullOrbitResult = HashMap::default();

        // --- Timer configuration
        let poll_interval = Duration::from_millis(20); // intervalle minimal
        let mut last_poll = Instant::now();

        for (obj, observations) in self.iter_mut() {
            // --- Vérifie si l’intervalle est écoulé
            if last_poll.elapsed() >= poll_interval {
                if should_cancel() {
                    break; // Stoppe la boucle
                }
                last_poll = Instant::now();
            }

            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);
        }

        results
    }

    #[inline]
    fn total_observations(&self) -> usize {
        self.values().map(|obs: &Observations| obs.len()).sum()
    }

    #[inline]
    fn number_of_trajectories(&self) -> usize {
        self.len()
    }

    fn obs_count_stats(&self) -> Option<ObsCountStats> {
        // Collect sizes (one pass, O(N))
        let mut counts: Vec<usize> = self.values().map(|obs| obs.len()).collect();
        if counts.is_empty() {
            return None;
        }

        // Sort once, O(N log N). `unstable` is fine since we only need order.
        counts.sort_unstable();

        #[inline]
        fn q_index(n: usize, q: f64) -> usize {
            // Nearest-rank on [0, n-1] using linear index; robust for small n.
            let pos = q * (n as f64 - 1.0);
            let idx = pos.round() as isize;
            idx.clamp(0, (n as isize) - 1) as usize
        }

        let n = counts.len();
        let min = counts[0];
        let max = counts[n - 1];
        let p25 = counts[q_index(n, 0.25)];
        let median = counts[q_index(n, 0.50)];
        let p95 = counts[q_index(n, 0.95)];

        Some(ObsCountStats {
            min,
            p25,
            median,
            p95,
            max,
        })
    }
}
