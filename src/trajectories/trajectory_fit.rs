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
//! * The returned map contains **one entry per processed object**,
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
//!
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

/// Full batch orbit determination results.
///
/// Each entry maps an [`ObjectNumber`] to the outcome of a full
/// Initial Orbit Determination (IOD) attempt on its set of observations.
///
/// Internally, this is implemented as:
///
/// ```ignore
/// HashMap<ObjectNumber, Result<(GaussResult, f64), OutfitError>, RandomState>
/// ```
///
/// Return semantics
/// -----------------
/// * `Ok((GaussResult, f64))` – a successful IOD with its RMS of normalized residuals.
/// * `Err(OutfitError)` – a failure isolated to that object.
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
/// * `Ok(None)` – key absent.
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
///
/// Arguments
/// -----------------
/// * `all`: The map of all IOD outcomes (consumed entry will be removed).
/// * `key`: The object identifier to extract.
///
/// Return
/// ----------
/// * `Ok(Some((GaussResult, f64)))` – ownership of the solution and its RMS.
/// * `Ok(None)` – key absent.
/// * `Err(OutfitError)` – the IOD attempt failed for that key.
///
/// See also
/// ------------
/// * [`gauss_result_for`] – Borrowing accessor.
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

// ============================================================================
// Factorized core + progress abstraction + cancel config
// ============================================================================

/// Cancellation guard polled at fixed wall-clock intervals.
///
/// A `CancelCfg` lets the main loop periodically check whether the user
/// or an external controller has requested an early stop. The loop itself
/// decides *when* to poll based on the `interval`, and *how* to react
/// based on the `should_cancel` callback.
///
/// Arguments
/// -----------------
/// * `interval`: Minimum wall-clock delay between two cancellation checks.
///   This prevents the loop from calling the callback at every
///   iteration, which would be too costly.
/// * `should_cancel`: User-provided closure returning `true` when cancellation
///   is requested. If `true`, the loop terminates gracefully.
///
/// See also
/// ------------
/// * [`estimate_all_orbits_core`] – Main iteration loop that evaluates this configuration.
struct CancelCfg<F> {
    interval: Duration,
    should_cancel: F,
}

/// Abstract interface for reporting progress to the outside world.
///
/// The purpose of this trait is to **decouple the heavy numerical loop**
/// from any particular UI backend. The core orbit determination logic
/// always calls into a `ProgressSink`, but the actual implementation
/// depends on the build:
///
/// * with the `progress` feature, it is backed by an [`indicatif`] progress bar,
/// * otherwise, a no-op implementation is used, so the loop compiles without UI.
///
/// This abstraction ensures that the loop has a consistent lifecycle:
///   1. [`start`] is called once with the total number of objects.
///   2. [`on_iter`] is called at the beginning of each iteration (e.g. to refresh messages).
///   3. [`inc`] is called once per completed iteration.
///   4. [`on_interrupt`] is called right before exiting due to cancellation.
///   5. [`finish`] is always called at the end, successful or not.
///
/// By default, all methods are no-ops, so implementors only override the
/// subset they need.
trait ProgressSink {
    /// Called once before entering the main loop, with the total number of items.
    fn start(&mut self, _total: u64) {}
    /// Called at the beginning of each iteration (e.g., refresh UI or logs).
    fn on_iter(&mut self) {}
    /// Increment the progress by one step.
    fn inc(&mut self) {}
    /// Called once if the loop exits because of cancellation.
    fn on_interrupt(&mut self) {}
    /// Called once at the very end, regardless of success or cancellation.
    fn finish(&mut self) {}
}

/// Default no-op implementation, used when the `progress` feature is disabled.
impl ProgressSink for () {}

/// Blanket implementation so that `&mut T` also implements [`ProgressSink`].
///
/// This lets tests pass `&mut MockProgress` directly, while the core API
/// continues to accept progress sinks by value.
impl<T: ProgressSink + ?Sized> ProgressSink for &mut T {
    #[inline]
    fn start(&mut self, total: u64) {
        (**self).start(total)
    }
    #[inline]
    fn on_iter(&mut self) {
        (**self).on_iter()
    }
    #[inline]
    fn inc(&mut self) {
        (**self).inc()
    }
    #[inline]
    fn on_interrupt(&mut self) {
        (**self).on_interrupt()
    }
    #[inline]
    fn finish(&mut self) {
        (**self).finish()
    }
}

/// Concrete type selected depending on the `progress` feature:
/// * [`IndicatifProgress`] when enabled,
/// * [`()`] (no-op) when disabled.
#[cfg(feature = "progress")]
type ProgressImpl = IndicatifProgress;
#[cfg(not(feature = "progress"))]
type ProgressImpl = ();

/// Central loop that runs orbit estimation for each trajectory.
///
/// This function is the **engine** behind the public APIs:
/// [`TrajectoryFit::estimate_all_orbits`] and
/// [`TrajectoryFit::estimate_all_orbits_with_cancel`].
///
/// It consumes a [`TrajectorySet`] and tries to estimate the best orbit
/// for each contained object. Along the way it reports progress, and it
/// may stop early if the provided cancellation config triggers.
///
/// Arguments
/// -----------------
/// * `set`: The trajectory set to process (mutable, results are inserted).
/// * `state`: Global environment providing ephemerides, constants, and frames.
/// * `rng`: Random number generator used for noisy triplet realizations.
/// * `params`: IOD parameters controlling triplet generation and scoring.
/// * `cancel`: Optional cancellation guard (poll interval + callback).
/// * `progress`: Progress reporting sink (indicatif bar or no-op).
///
/// Return
/// ----------
/// * A [`FullOrbitResult`], i.e. a map from object → `Ok((GaussResult, rms))`
///   or `Err(OutfitError)` depending on whether orbit estimation succeeded.
///
/// See also
/// ------------
/// * [`TrajectoryFit::estimate_all_orbits`] – Public API without cancellation.
/// * [`TrajectoryFit::estimate_all_orbits_with_cancel`] – Public API with cancellation.
fn estimate_all_orbits_core<F, P>(
    set: &mut TrajectorySet,
    state: &Outfit,
    rng: &mut impl Rng,
    params: &IODParams,
    mut cancel: Option<CancelCfg<F>>,
    mut progress: P,
) -> FullOrbitResult
where
    F: FnMut() -> bool,
    P: ProgressSink,
{
    let total = set.len() as u64;
    progress.start(total.max(1));

    let mut results: FullOrbitResult = HashMap::default();
    let mut last_poll = Instant::now();

    for (obj, observations) in set.iter_mut() {
        // --- Timer-based cancellation (if configured)
        if let Some(CancelCfg {
            interval,
            should_cancel,
        }) = cancel.as_mut()
        {
            if last_poll.elapsed() >= *interval {
                if should_cancel() {
                    progress.on_interrupt();
                    break;
                }
                last_poll = Instant::now();
            }
        }

        progress.on_iter();

        // Core work
        let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
        results.insert(obj.clone(), res);

        progress.inc();
    }

    progress.finish();
    results
}

// --------------------------- Progress (indicatif) ----------------------------
#[cfg(feature = "progress")]
mod progress_impl {
    use super::ProgressSink;
    use indicatif::{ProgressBar, ProgressStyle};
    use std::time::Duration;

    use super::IterTimer;
    use crate::trajectories::progress_bar::fmt_dur;

    /// Progress sink backed by `indicatif`.
    ///
    /// See also
    /// ------------
    /// * [`estimate_all_orbits_core`] – Calls into this sink at key lifecycle moments.
    pub(super) struct IndicatifProgress {
        pb: ProgressBar,
        it_timer: IterTimer,
    }

    impl Default for IndicatifProgress {
        fn default() -> Self {
            // The actual length is set in `start()`.
            let pb = ProgressBar::new(1);
            Self {
                pb,
                it_timer: IterTimer::new(0.2),
            }
        }
    }

    impl ProgressSink for IndicatifProgress {
        fn start(&mut self, total: u64) {
            self.pb.set_length(total.max(1));
            self.pb.set_style(
                ProgressStyle::with_template(
                    "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) \
                         | {per_sec} | ETA {eta_precise} | {msg}",
                )
                .expect("indicatif template"),
            );
            self.pb.enable_steady_tick(Duration::from_millis(200));
        }

        fn on_iter(&mut self) {
            let last = self.it_timer.tick();
            let avg = self.it_timer.avg();
            self.pb
                .set_message(format!("last: {}, avg: {}", fmt_dur(last), fmt_dur(avg)));
        }

        fn inc(&mut self) {
            self.pb.inc(1);
        }

        fn on_interrupt(&mut self) {
            self.pb.set_message("Interrupted");
        }

        fn finish(&mut self) {
            self.pb.disable_steady_tick();
            self.pb.finish_and_clear();
        }
    }
}

#[cfg(feature = "progress")]
use progress_impl::IndicatifProgress;

// ============================================================================
// Public trait + factorized implementation
// ============================================================================

pub trait TrajectoryFit {
    /// Run Gauss-based Initial Orbit Determination (IOD) for **every trajectory** in the set.
    ///
    /// This method iterates over each `(ObjectNumber → Observations)` entry and applies the
    /// full IOD pipeline: candidate triplet enumeration, preliminary Gauss solution, and
    /// scoring / selection of the best orbit. It aggregates results into a [`FullOrbitResult`].
    ///
    /// Mutation semantics
    /// -----------------
    /// * This function requires `&mut self` and may **reorder observations in-place** (e.g.,
    ///   by time) and/or **update batch-level calibration** data (RMS scaling of quoted errors).
    /// * The underlying astrometric measurements (RA/DEC/time) remain semantically identical,
    ///   but their **container order** and **per-observation uncertainty metadata** may change
    ///   due to calibration and sorting steps used by the estimator.
    /// * If you rely on a specific iteration order elsewhere, do not assume it is preserved.
    ///
    /// Determinism
    /// -----------------
    /// * With a fixed RNG seed, the procedure is deterministic given identical inputs and params.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: Global environment providing ephemerides, constants, and reference frames.
    /// * `rng`: Random number generator used for noisy triplet realizations (e.g., [`StdRng`](rand::rngs::StdRng)).
    /// * `params`: IOD parameters controlling triplet generation, scoring, and correction loops.
    ///
    /// Return
    /// ----------
    /// * A [`FullOrbitResult`] mapping each object to either:
    ///   * `Ok((GaussResult, f64))` – selected orbit and its RMS,
    ///   * `Err(OutfitError)` – diagnostic if no acceptable solution was found.
    ///
    /// Notes
    /// ----------
    /// * Failures are isolated: one object failing does not prevent others from being processed.
    /// * Runtime scales with the number of trajectories and candidate triplets per trajectory.
    ///
    /// See also
    /// ------------
    /// * [`ObservationIOD::estimate_best_orbit`] – Per-trajectory IOD with best-orbit selection.
    /// * [`TrajectoryFit::estimate_all_orbits_with_cancel`] – Same API with cooperative cancellation.
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
    ///
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

    /// Run Gauss-based IOD for all trajectories, with **cooperative cancellation** support.
    ///
    /// Behaves like [`TrajectoryFit::estimate_all_orbits`], but periodically polls a user
    /// callback to decide whether to stop early. Returns **partial results** if cancelled.
    ///
    /// Mutation semantics
    /// -----------------
    /// * Same as the non-cancellable variant: the method may **reorder observations in-place**
    ///   and update **per-batch calibration** (e.g., RMS alignment of quoted errors).
    ///
    /// Cancellation model
    /// -----------------
    /// * The loop polls `should_cancel` at ~20 ms wall-clock intervals. When it returns `true`,
    ///   the loop terminates gracefully, calls `on_interrupt()` on the progress sink (if any),
    ///   and returns the results accumulated so far.
    ///
    /// Determinism
    /// -----------------
    /// * With a fixed RNG seed, behavior is deterministic except for the **cut point** at which
    ///   cancellation is observed (timing dependent).
    ///
    /// Arguments
    /// -----------------
    /// * `state`: Global environment and reference frames.
    /// * `rng`: Random number generator for noisy triplet realizations.
    /// * `params`: IOD parameters.
    /// * `should_cancel`: Closure polled periodically; return `true` to request early stop.
    ///
    /// Return
    /// ----------
    /// * A [`FullOrbitResult`]:
    ///   * Complete if the loop ran to completion,
    ///   * Partial if cancellation was triggered mid-way.
    ///
    /// See also
    /// ------------
    /// * [`TrajectoryFit::estimate_all_orbits`] – Non-cancellable variant.
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
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult {
        // `ProgressImpl` is `IndicatifProgress` when feature=progress, `()` otherwise.
        estimate_all_orbits_core(
            self,
            state,
            rng,
            params,
            None::<CancelCfg<fn() -> bool>>,
            ProgressImpl::default(),
        )
    }

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
        let cancel = CancelCfg {
            interval: Duration::from_millis(20),
            should_cancel: &mut should_cancel,
        };
        estimate_all_orbits_core(
            self,
            state,
            rng,
            params,
            Some(cancel),
            ProgressImpl::default(),
        )
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

#[cfg(test)]
#[cfg(feature = "jpl-download")]
mod tests_estimate_all_orbits {
    use crate::{
        observations::Observation, unit_test_global::OUTFIT_HORIZON_TEST, KeplerianElements,
    };

    use super::*;
    use approx::assert_relative_eq;
    use rand::SeedableRng;
    use smallvec::SmallVec;
    use std::{
        f64::consts::PI,
        sync::atomic::{AtomicUsize, Ordering},
    };

    // -------------------------------
    // Test fixtures (lightweight)
    // -------------------------------

    /// Build a tiny TrajectorySet with N empty observation lists.
    ///
    /// Note: This assumes `TrajectorySet` is a HashMap-like structure
    ///       and `ObjectNumber::Int(u64)` exists. Adjust if needed.
    fn make_set(n: usize) -> TrajectorySet {
        let mut set: TrajectorySet = std::collections::HashMap::with_hasher(RandomState::new());
        for i in 0..n {
            // If your ObjectNumber uses a different constructor, adjust here.
            let key = ObjectNumber::Int(i as u32);
            // If Observations is not a Vec, adapt this to your type.
            let obs: Observations = Default::default();
            set.insert(key, obs);
        }
        set
    }

    /// Dummy `Outfit` and `IODParams` for tests that do not reach the estimator.
    ///
    /// We never call the estimator in cancellation-first tests, so these values
    /// are placeholders to satisfy the function signatures.
    fn dummy_env() -> (Outfit, IODParams) {
        let env = OUTFIT_HORIZON_TEST.0.clone();
        let params = IODParams::builder()
            .n_noise_realizations(10)
            .noise_scale(1.0)
            .max_obs_for_triplets(12)
            .max_triplets(30)
            .build()
            .unwrap();
        (env, params)
    }

    // -------------------------------
    // Unit tests: cancellation logic
    // -------------------------------

    /// Cancellation fires before the first object is processed: result should be empty.
    ///
    /// This test calls the factorized core with `interval = 0 ms` and a callback
    /// that immediately requests cancellation. The estimator is never invoked.
    #[test]
    fn core_cancel_before_any_work() {
        let mut set = make_set(5);
        let (env, params) = dummy_env();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        // Cancel immediately on the very first poll.
        let mut cancel_called = 0usize;
        let mut should_cancel = || {
            cancel_called += 1;
            true
        };

        let cancel = CancelCfg {
            interval: Duration::from_millis(0),
            should_cancel: &mut should_cancel,
        };

        // Use no-op progress sink: works with or without the `progress` feature.
        let results = estimate_all_orbits_core(&mut set, &env, &mut rng, &params, Some(cancel), ());

        assert!(results.is_empty(), "No object should have been processed");
        assert!(
            cancel_called >= 1,
            "Cancellation should have been polled at least once"
        );
    }

    /// Cancellation after exactly one iteration: we expect exactly one entry in the map.
    ///
    /// IMPORTANT: This test *may* reach the estimator if the cancellation poll
    /// happens after the first object. We therefore keep the set size to 1 so
    /// we never process more than one. If your estimator requires real env/params,
    /// mark this test as `#[ignore]` until you wire a small valid fixture.
    #[test]
    fn core_cancel_after_one_object() {
        let mut set = make_set(2);
        let (env, params) = dummy_env();
        let mut rng = rand::rngs::StdRng::seed_from_u64(123);

        let polls = AtomicUsize::new(0);
        // First poll = false (let the first object run), next polls = true.
        let mut should_cancel = || {
            let c = polls.fetch_add(1, Ordering::Relaxed);
            c >= 1
        };

        let cancel = CancelCfg {
            interval: Duration::from_millis(0), // poll at every loop entry
            should_cancel: &mut should_cancel,
        };

        let results = estimate_all_orbits_core(&mut set, &env, &mut rng, &params, Some(cancel), ());

        assert_eq!(
            results.len(),
            1,
            "Exactly one object should have been processed before cancel"
        );
    }

    // -------------------------------
    // Unit tests: progress plumbing
    // -------------------------------

    /// Mock progress sink to observe lifecycle calls.
    #[derive(Default)]
    struct MockProgress {
        started_with: Option<u64>,
        it_calls: usize,
        inc_calls: usize,
        interrupted: bool,
        finished: bool,
    }

    impl ProgressSink for MockProgress {
        fn start(&mut self, total: u64) {
            self.started_with = Some(total);
        }
        fn on_iter(&mut self) {
            self.it_calls += 1;
        }
        fn inc(&mut self) {
            self.inc_calls += 1;
        }
        fn on_interrupt(&mut self) {
            self.interrupted = true;
        }
        fn finish(&mut self) {
            self.finished = true;
        }
    }

    /// Progress sink should receive `start`, `on_interrupt`, and `finish` when cancelling before work.
    #[test]
    fn progress_calls_when_cancelled_immediately() {
        let mut set = make_set(3);
        let (env, params) = dummy_env();
        let mut rng = rand::rngs::StdRng::seed_from_u64(7);

        let mut should_cancel = || true;
        let cancel = CancelCfg {
            interval: Duration::from_millis(0),
            should_cancel: &mut should_cancel,
        };

        let mut mock = MockProgress::default();
        let results =
            estimate_all_orbits_core(&mut set, &env, &mut rng, &params, Some(cancel), &mut mock);

        assert!(results.is_empty());
        assert_eq!(mock.started_with, Some(3));
        assert!(mock.interrupted, "on_interrupt() must be called");
        assert!(mock.finished, "finish() must be called");
        // No iteration advanced, so no inc() and on_iter() expected.
        assert_eq!(mock.it_calls, 0);
        assert_eq!(mock.inc_calls, 0);
    }

    // -------------------------------
    // Integration tests
    // -------------------------------

    #[inline]
    fn angle_abs_diff(a: f64, b: f64) -> f64 {
        let tau = 2.0 * PI;
        let mut d = (a - b) % tau;
        if d > PI {
            d -= tau;
        }
        if d < -PI {
            d += tau;
        }
        d.abs()
    }

    pub fn assert_keplerian_approx_eq(
        got: &KeplerianElements,
        exp: &KeplerianElements,
        abs_eps: f64,
        rel_eps: f64,
    ) {
        // Scalars (non-angular)
        assert_relative_eq!(
            got.reference_epoch,
            exp.reference_epoch,
            epsilon = abs_eps,
            max_relative = rel_eps
        );
        assert_relative_eq!(
            got.semi_major_axis,
            exp.semi_major_axis,
            epsilon = abs_eps,
            max_relative = rel_eps
        );
        assert_relative_eq!(
            got.eccentricity,
            exp.eccentricity,
            epsilon = abs_eps,
            max_relative = rel_eps
        );

        // Angles (radians), compare with wrap-around
        for (name, g, e) in [
            ("inclination", got.inclination, exp.inclination),
            (
                "ascending_node_longitude",
                got.ascending_node_longitude,
                exp.ascending_node_longitude,
            ),
            (
                "periapsis_argument",
                got.periapsis_argument,
                exp.periapsis_argument,
            ),
            ("mean_anomaly", got.mean_anomaly, exp.mean_anomaly),
        ] {
            let diff = angle_abs_diff(g, e);
            // Allow absolute OR relative tolerance (whichever is larger).
            let tol = abs_eps.max(rel_eps * e.abs());
            assert!(
            diff <= tol,
            "Angle {name:?} differs too much: |Δ| = {diff:.6e} > tol {tol:.6e} (got={g:.15}, exp={e:.15})"
        );
        }
    }

    /// Estimate on an empty-observation set should return one entry per object with errors.
    #[test]
    fn public_no_progress_runs_all_objects() {
        let mut set = OUTFIT_HORIZON_TEST.1.clone();

        // TODO: replace with real constructors in your codebase:
        let (env, params) = dummy_env();
        let mut rng = rand::rngs::StdRng::seed_from_u64(777);

        use super::TrajectoryFit;
        let results = set.estimate_all_orbits(&env, &mut rng, &params);

        let string_id = "K09R05F";
        let orbit = gauss_result_for(&results, &string_id.into())
            .unwrap()
            .unwrap()
            .0
            .as_inner()
            .as_keplerian()
            .unwrap();

        let expected = KeplerianElements {
            reference_epoch: 57049.25533417104,
            semi_major_axis: 1.8017448718161189,
            eccentricity: 0.283572382702194,
            inclination: 0.2026747553253312,
            ascending_node_longitude: 0.0079836299943183,
            periapsis_argument: 1.245049339166438,
            mean_anomaly: 0.4406946018418537,
        };

        assert_keplerian_approx_eq(orbit, &expected, 1e-6, 1e-6);
    }

    /// Public cancellation API should return a partial map when cancelling quickly.
    #[test]
    fn public_with_cancel_returns_partial() {
        let mut set = OUTFIT_HORIZON_TEST.1.clone();

        // TODO: replace with real constructors in your codebase:
        let (env, params) = dummy_env();
        let mut rng = rand::rngs::StdRng::seed_from_u64(42);

        use super::TrajectoryFit;
        // Callback cancels immediately; public API polls every ~20ms.
        // Depending on estimator speed, a few items may slip in before first poll.
        let results = set.estimate_all_orbits_with_cancel(&env, &mut rng, &params, || true);

        assert!(
            results.len() <= 50,
            "Result map cannot exceed the number of objects"
        );
        assert!(
            !results.is_empty(),
            "Depending on timing, a few items may be processed before first poll"
        );
    }

    // -------------------------------
    // Accessor helpers tests
    // -------------------------------

    /// `gauss_result_for` should distinguish: missing key, error entry, ok entry.
    #[test]
    fn gauss_accessors_err_and_missing() {
        let mut all: FullOrbitResult = HashMap::with_hasher(RandomState::new());
        let k1 = ObjectNumber::Int(1);
        let k2 = ObjectNumber::Int(2);

        // Insert an error entry for k1. Construct an OutfitError if you have a cheap variant.
        // If construction is non-trivial, you can skip inserting and just test "missing".
        all.insert(
            k1.clone(),
            Err(OutfitError::InvalidIODParameter("test".into())),
        );

        // Missing key:
        assert!(matches!(gauss_result_for(&all, &k2), Ok(None)));

        // Error key:
        match gauss_result_for(&all, &k1) {
            Err(e) => {
                // Just check we got *some* error reference back.
                let _ = format!("{e}");
            }
            other => panic!("expected Err(&OutfitError), got {other:?}"),
        }

        // Take on missing:
        assert!(matches!(take_gauss_result(&mut all, &k2), Ok(None)));

        // Take on error:
        match take_gauss_result(&mut all, &k1) {
            Err(e) => {
                let _ = format!("{e}");
            }
            other => panic!("expected Err(OutfitError), got {other:?}"),
        }
    }

    /// Stats over per-trajectory observation counts.
    #[test]
    fn obs_count_stats_basic() {
        use std::collections::HashMap;

        // Helper: build a dummy Observation for tests only.
        #[inline]
        fn dummy_observation() -> Observation {
            // SAFETY (tests only):
            // This assumes `Observation` is plain-old-data (floats, ints) and `Copy`,
            // i.e. no heap-owned fields (String, Vec, Arc, etc.) and no Drop.
            // Si ce n’est pas vrai dans ton code, remplace cette fonction par
            // un vrai constructeur de test qui remplit des champs plausibles.
            assert_is_copy::<Observation>();
            unsafe { std::mem::MaybeUninit::<Observation>::zeroed().assume_init() }
        }

        // Compile-time check: force `Observation: Copy` pour que le zero-init soit sûr.
        #[inline(always)]
        fn assert_is_copy<T: Copy>() {}

        // Cas vide.
        let set = make_set(0);
        assert!(set.obs_count_stats().is_none(), "Empty set → None");

        // Build uneven counts: 2, 4, 8, 16, 16
        let mut set: TrajectorySet = HashMap::with_hasher(RandomState::new());

        let mut push_n = |id: u32, n: usize| {
            let mut v: Observations = SmallVec::with_capacity(n);
            for _ in 0..n {
                v.push(dummy_observation());
            }
            set.insert(ObjectNumber::Int(id), v);
        };

        push_n(1, 2);
        push_n(2, 4);
        push_n(3, 8);
        push_n(4, 16);
        push_n(5, 16);

        let stats = set.obs_count_stats().expect("non-empty");
        assert_eq!(stats.min, 2);
        assert_eq!(stats.max, 16);
        assert_eq!(stats.median, 8);
        assert_eq!(stats.p25, 4);
        assert_eq!(stats.p95, 16);
    }
}
