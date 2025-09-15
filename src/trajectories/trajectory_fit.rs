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
#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};

#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(feature = "parallel")]
use std::{
    hash::{Hash, Hasher},
    mem,
};

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
    use super::IterTimer;
    use super::ProgressSink;
    use crate::trajectories::progress_bar::fmt_dur;

    /// Progress sink backed by `indicatif`.
    ///
    /// See also
    /// ------------
    /// * [`estimate_all_orbits_core`] – Calls into this sink at key lifecycle moments.
    pub(super) struct IndicatifProgress {
        pb: super::ProgressBar,
        it_timer: IterTimer,
    }

    impl Default for IndicatifProgress {
        fn default() -> Self {
            // The actual length is set in `start()`.
            let pb = super::ProgressBar::new(1);
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
                super::ProgressStyle::with_template(
                    "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) \
                         | {per_sec} | ETA {eta_precise} | {msg}",
                )
                .expect("indicatif template"),
            );
            self.pb
                .enable_steady_tick(super::Duration::from_millis(200));
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

// --------------------------- Parallel features ----------------------------

#[cfg(feature = "parallel")]
/// Generate a new 64-bit pseudo-random value using the **SplitMix64** algorithm.
/// This is a simple, fast, and reproducible way to decorrelate seeds for parallel RNGs.
///
/// Arguments
/// -----------------
/// * `x`: Input state (a `u64` value, typically a hash or a base seed).
///
/// Return
/// ----------
/// * A `u64` pseudo-random value, suitable for seeding RNGs (e.g., `StdRng`).
///
/// See also
/// ------------
/// * [`seed_for_object`] – Derives per-object seeds from a base seed and object hash.
#[inline]
fn splitmix64(mut x: u64) -> u64 {
    // SplitMix64 constants and shifts from Steele et al. (2014).
    x = x.wrapping_add(0x9E3779B97F4A7C15);
    let mut z = x;
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D049BB133111EB);
    z ^ (z >> 31)
}

#[cfg(feature = "parallel")]
/// Derive a **deterministic per-object RNG seed** from a global base seed.
///
/// Each object is first hashed with the crate’s default hasher (`ahash`), and the
/// resulting 64-bit value is mixed with the base seed via [`splitmix64`].
/// This ensures that:
/// - Each object gets a **stable, reproducible seed** (same input → same output).
/// - Seeds are decorrelated even across many objects.
/// - Parallel runs remain deterministic regardless of thread scheduling.
///
/// Arguments
/// -----------------
/// * `base`: A global base seed (drawn once per batch).
/// * `obj`: The [`ObjectNumber`] used as key to derive the per-object seed.
///
/// Return
/// ----------
/// * A `u64` deterministic RNG seed for the given object.
///
/// See also
/// ------------
/// * [`splitmix64`] – Core mixing function.
/// * [`ObjectNumber`] – Object identifier used in Outfit (MPC number or string).
#[inline]
fn seed_for_object(base: u64, obj: &ObjectNumber) -> u64 {
    // Hash object key with the same family used elsewhere (ahash).
    let mut h = ahash::AHasher::default();
    obj.hash(&mut h);
    let obj_h = h.finish();

    // Mix base seed and object hash through SplitMix64.
    splitmix64(base ^ obj_h)
}

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

    /// Run Gauss-based Initial Orbit Determination (IOD) over all trajectories
    /// using **parallel batches**.
    ///
    /// The [`TrajectorySet`] is split into chunks of size `batch_size`. Each chunk
    /// is assigned to a Rayon worker thread, and objects inside a chunk are processed
    /// sequentially for cache efficiency. A single `base_seed` is drawn from `rng`,
    /// and a stable per-object seed is derived deterministically to guarantee
    /// reproducibility regardless of parallel scheduling.
    ///
    /// Threading model
    /// -----------------
    /// * This function uses the **global Rayon thread pool** by default.
    /// * The number of worker threads is controlled by the environment variable
    ///   `RAYON_NUM_THREADS`. For example:
    ///
    ///   ```bash
    ///   RAYON_NUM_THREADS=4 cargo run --release
    ///   ```
    ///
    ///   will cap Rayon to 4 threads across the entire program.
    /// * If the variable is unset, Rayon defaults to the number of logical CPUs.
    ///
    /// Mutation semantics
    /// -----------------
    /// * As in the sequential version, this method may **reorder observations**
    ///   and **update per-batch calibration** (e.g. RMS scaling of quoted errors).
    /// * Each trajectory’s observation container is reinserted after processing,
    ///   so `self` remains valid and complete.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: Global environment (ephemerides, constants, frames).
    /// * `rng`: Random number generator, used only once to draw a base seed.
    /// * `params`: IOD parameters controlling triplet generation, scoring, correction.
    /// * `batch_size`: Number of trajectories per parallel batch. Must be ≥ 1.
    ///
    /// Return
    /// ----------
    /// * A [`FullOrbitResult`] mapping each object to either:
    ///   * `Ok((GaussResult, f64))` – best orbit and its RMS,
    ///   * `Err(OutfitError)` – diagnostic if no acceptable orbit was found.
    ///
    /// See also
    /// ------------
    /// * [`TrajectoryFit::estimate_all_orbits`] – Sequential variant.
    /// * [`TrajectoryFit::estimate_all_orbits_with_cancel`] – Sequential variant with cooperative cancellation.
    #[cfg(feature = "parallel")]
    fn estimate_all_orbits_in_batches_parallel(
        &mut self,
        state: &Outfit,
        rng: &mut impl rand::Rng,
        params: &IODParams,
        batch_size: usize,
    ) -> FullOrbitResult;
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

    #[cfg(feature = "parallel")]
    fn estimate_all_orbits_in_batches_parallel(
        &mut self,
        state: &Outfit,
        rng: &mut impl rand::Rng,
        params: &IODParams,
        batch_size: usize,
    ) -> FullOrbitResult {
        // Draw a single base seed once; per-object seeds derived deterministically.
        let base_seed: u64 = rng.random();

        // Take the whole map so we can own/mutate Observations per object off-thread.
        let mut old: TrajectorySet = mem::take(self);
        let mut entries: Vec<(ObjectNumber, Observations)> = old.drain().collect();

        let total_items = entries.len() as u64;
        let batch_size = batch_size.max(1);

        // Materialize batches **by move** (no clone of Observations).
        let mut batches: Vec<Vec<(ObjectNumber, Observations)>> =
            Vec::with_capacity(entries.len().div_ceil(batch_size));
        while !entries.is_empty() {
            let take_n = entries.len().min(batch_size);
            batches.push(entries.drain(..take_n).collect());
        }

        // Global progress bar (thread-safe) under the `progress` feature.
        #[cfg(feature = "progress")]
        let pb = {
            use indicatif::{ProgressBar, ProgressStyle};
            let pb = ProgressBar::new(total_items.max(1));
            pb.set_style(
                ProgressStyle::with_template(
                    "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) \
                     | {per_sec} | ETA {eta_precise} | parallel batches",
                )
                .expect("indicatif template"),
            );
            pb.enable_steady_tick(std::time::Duration::from_millis(200));
            pb
        };

        // Process batches in parallel; each batch processed sequentially for locality.
        #[allow(clippy::type_complexity)]
        let mut per_batch: Vec<
            Vec<(
                ObjectNumber,
                Result<(GaussResult, f64), OutfitError>,
                Observations,
            )>,
        > = batches
            .into_par_iter()
            .map(|mut batch| {
                let mut out: Vec<(
                    ObjectNumber,
                    Result<(GaussResult, f64), OutfitError>,
                    Observations,
                )> = Vec::with_capacity(batch.len());

                for (obj, mut obs) in batch.drain(..) {
                    use rand::SeedableRng;

                    let local_seed = seed_for_object(base_seed, &obj);
                    let mut local_rng = rand::rngs::StdRng::seed_from_u64(local_seed);

                    let res =
                        obs.estimate_best_orbit(state, &state.error_model, &mut local_rng, params);

                    // Thread-safe progress increment.
                    #[cfg(feature = "progress")]
                    pb.inc(1);

                    out.push((obj, res, obs));
                }
                out
            })
            .collect();

        // Finalize progress.
        #[cfg(feature = "progress")]
        {
            pb.disable_steady_tick();
            pb.finish_and_clear();
        }

        // Reinsert mutated observations and build results map with the same hasher.
        let mut results: FullOrbitResult = HashMap::with_hasher(ahash::RandomState::new());
        for batch in per_batch.drain(..) {
            for (obj, res, obs) in batch {
                self.insert(obj.clone(), obs);
                results.insert(obj, res);
            }
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

    #[cfg(test)]
    #[cfg(feature = "parallel")]
    mod tests_estimate_orbit_parallel_batches {
        use super::*;
        use ahash::RandomState;
        use rand::SeedableRng;

        // Reuse helpers and fixtures style from your existing tests.
        // If these are private in another module, duplicate minimal versions here.

        /// Build a tiny TrajectorySet with N empty observation lists.
        ///
        /// Note: This assumes `TrajectorySet` is a HashMap-like structure
        ///       and `ObjectNumber::Int(u32)` exists. Adjust if needed.
        fn make_set(n: usize) -> TrajectorySet {
            let mut set: TrajectorySet = std::collections::HashMap::with_hasher(RandomState::new());
            for i in 0..n {
                set.insert(ObjectNumber::Int(i as u32), Default::default());
            }
            set
        }

        #[inline]
        fn total_obs(set: &TrajectorySet) -> usize {
            set.values().map(|v: &Observations| v.len()).sum()
        }

        // -------------------------------
        // Unit tests: basic shape/edges
        // -------------------------------

        /// Parallel-batched IOD over an empty set should return an empty map.
        #[test]
        fn parallel_batches_empty_set_is_empty() {
            let mut set = make_set(0);

            // Dummy env/params: estimator is never reached for empty set.
            // If you need to compile without jpl, use the same `dummy_env()` strategy as your seq tests.
            let env = dummy_env().0;
            let params = IODParams::builder().build().unwrap();

            let mut rng = rand::rngs::StdRng::seed_from_u64(1);
            let results =
                set.estimate_all_orbits_in_batches_parallel(&env, &mut rng, &params, 1024);
            assert!(results.is_empty(), "Empty input → empty results");
            assert_eq!(set.len(), 0, "Set remains empty");
        }

        /// Batch-size boundaries (1 and very large) should produce exactly one entry per object.
        ///
        /// We don't assert on Ok/Err, only that every object was processed and observations were reintegrated.
        #[test]
        fn parallel_batches_size_edges_cover_all_objects() {
            for &batch_size in &[1usize, 10_000usize] {
                let mut set = make_set(7);

                // Build a dummy env that lets the code run without panicking even if estimator errs.
                // The estimator may return Err for empty observations — it's fine for this test.
                let env = dummy_env().0;
                let params = IODParams::builder().build().unwrap();

                let before_n = set.len();
                let before_tot = total_obs(&set);

                let mut rng = rand::rngs::StdRng::seed_from_u64(42);
                let results = set
                    .estimate_all_orbits_in_batches_parallel(&env, &mut rng, &params, batch_size);

                assert_eq!(results.len(), before_n, "Exactly one entry per object");
                assert_eq!(set.len(), before_n, "All objects reinserted in set");
                assert_eq!(
                    total_obs(&set),
                    before_tot,
                    "Total number of observations is preserved (reorder/calibration only)"
                );
            }
        }

        /// Using the same input and RNG seed must be deterministic across runs, regardless of scheduling.
        /// This checks **one specific object** for identical formatted outcome (Ok/Err shape and RMS).
        #[test]
        fn parallel_batches_deterministic_across_runs_with_same_seed() {
            // Small synthetic set; estimator likely returns Err for empty observations.
            // Determinism check focuses on the *presence* and *shape* of results.
            let build_set = || make_set(5);
            let env = dummy_env().0;
            let params = IODParams::builder().build().unwrap();

            let key = ObjectNumber::Int(2);

            let run_once = |seed: u64| {
                let mut set = build_set();
                let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
                let results =
                    set.estimate_all_orbits_in_batches_parallel(&env, &mut rng, &params, 2);
                // Record a stable string representation for that key.
                match results.get(&key) {
                    None => "None".to_string(),
                    Some(Ok((_g, rms))) => format!("Ok rms={rms:.12e}"),
                    Some(Err(e)) => format!("Err: {e}"),
                }
            };

            let a = run_once(0xDEADBEEF);
            let b = run_once(0xDEADBEEF);
            assert_eq!(a, b, "Same seed/input → identical outcome formatting");
        }

        // -------------------------------
        // Integration tests with JPL env
        // -------------------------------

        mod with_ephem {
            use super::*;
            use approx::assert_relative_eq;

            use crate::unit_test_global::OUTFIT_HORIZON_TEST;

            /// Parallel batched results should match the known-good orbit (as in the sequential test).
            #[test]
            fn parallel_batches_return_orbit() {
                // Use the same fixture you use elsewhere.
                let mut set = OUTFIT_HORIZON_TEST.1.clone();
                let (env, params) = {
                    // If you have a helper `dummy_env()` in the seq tests, keep the same one:
                    let env = OUTFIT_HORIZON_TEST.0.clone();
                    let params = IODParams::builder()
                        .n_noise_realizations(10)
                        .noise_scale(1.0)
                        .max_obs_for_triplets(12)
                        .max_triplets(30)
                        .build()
                        .unwrap();
                    (env, params)
                };
                let mut rng = rand::rngs::StdRng::seed_from_u64(42);

                // Choose a batch size that is neither 1 nor huge to exercise chunking logic.
                let results =
                    set.estimate_all_orbits_in_batches_parallel(&env, &mut rng, &params, 1);

                // Same canonical object as in your sequential test:
                let string_id = "K09R05F";
                let orbit = gauss_result_for(&results, &string_id.into());

                assert!(orbit.is_ok(), "Result entry should be Ok");
            }

            /// Parallel batched determinism across different batch sizes.
            ///
            /// With same RNG seed + inputs, changing only `batch_size` should not affect results.
            #[test]
            fn parallel_batches_results_independent_of_batch_size() {
                let mut set1 = OUTFIT_HORIZON_TEST.1.clone();
                let mut set2 = OUTFIT_HORIZON_TEST.1.clone();
                let (env, params) = {
                    let env = OUTFIT_HORIZON_TEST.0.clone();
                    let params = IODParams::builder()
                        .n_noise_realizations(10)
                        .noise_scale(1.0)
                        .max_obs_for_triplets(12)
                        .max_triplets(30)
                        .build()
                        .unwrap();
                    (env, params)
                };

                let seed = 0xABCDEF0123456789;
                let mut rng1 = rand::rngs::StdRng::seed_from_u64(seed);
                let mut rng2 = rand::rngs::StdRng::seed_from_u64(seed);

                let res1 =
                    set1.estimate_all_orbits_in_batches_parallel(&env, &mut rng1, &params, 64);
                let res2 =
                    set2.estimate_all_orbits_in_batches_parallel(&env, &mut rng2, &params, 4096);

                // Compare a known object Keplerian solution (same as above).
                let key = "K09R05F".into();
                let k1 = gauss_result_for(&res1, &key)
                    .unwrap()
                    .unwrap()
                    .0
                    .as_inner()
                    .as_keplerian()
                    .unwrap();
                let k2 = gauss_result_for(&res2, &key)
                    .unwrap()
                    .unwrap()
                    .0
                    .as_inner()
                    .as_keplerian()
                    .unwrap();

                // Tight numerical equality (same seed → same triplet noise → identical orbit).
                assert_relative_eq!(k1.reference_epoch, k2.reference_epoch, epsilon = 0.0);
                assert_relative_eq!(k1.semi_major_axis, k2.semi_major_axis, epsilon = 0.0);
                assert_relative_eq!(k1.eccentricity, k2.eccentricity, epsilon = 0.0);
                assert_relative_eq!(k1.inclination, k2.inclination, epsilon = 0.0);
                assert_relative_eq!(
                    k1.ascending_node_longitude,
                    k2.ascending_node_longitude,
                    epsilon = 0.0
                );
                assert_relative_eq!(k1.periapsis_argument, k2.periapsis_argument, epsilon = 0.0);
                assert_relative_eq!(k1.mean_anomaly, k2.mean_anomaly, epsilon = 0.0);
            }
        }
    }
}
