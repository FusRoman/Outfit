//! # IOD Triplet Index Generator (lazy, windowed by time span)
//!
//! Streams **reduced indices** `(first, middle, last)` of time-feasible IOD triplets,
//! after sorting and downsampling the input observations. This is intended to be the
//! *index-level* counterpart of a Gauss triplet generator, allowing clients to compute
//! their own heuristics (e.g., a spacing weight) before materializing actual `GaussObs`.
//!
//! ## What “reduced indices” means
//! The generator first downsamples the full observation set to a smaller, **time-sorted**
//! subset (the “reduced” view). Indices yielded by the iterator refer to this reduced
//! view; use `selected_original_indices()` to map **reduced → original** indices.
//!
//! ## Feasibility constraints
//! For each yielded triplet `(first, middle, last)` with `first < middle < last`, the
//! following time-span constraint holds:
//!
//! ```text
//! dt_min ≤ t[last] - t[first] ≤ dt_max
//! ```
//!
//! where `t[*]` are the **reduced** epochs (same units as your observation times, e.g. TT/MJD).
//!
//! ## Why this generator?
//! - **Lazy**: no intermediate `Vec` of triplets; you can start consuming immediately.
//! - **Better complexity**: a per-anchor two-pointer window on `last` reduces the
//!   effective cost towards ~`O(n²)` in typical time distributions (vs. `O(n³)` brute force).
//! - **No overlapping borrows**: the generator **owns** its internal buffers (times,
//!   mapping), so you can iterate without holding long-lived borrows of the input.
//!
//! ## Typical flow
//! 1. Build with `TripletIndexGenerator::from_observations(...)` (sorts + downsamples).
//! 2. Iterate lazily over `(i, j, k)` **reduced** indices.
//! 3. If needed, map to originals via `gen.selected_original_indices()[i]`.
//! 4. (Optional) Compute a heuristic (e.g., spacing weight) and keep the best-K with a heap.
//! 5. Only then materialize full `GaussObs` from the original observations.
//!
//! ## Invariants per anchor
//! - `first < middle < last` always holds.
//! - For the current `first`, the **feasible window** for `last` is
//!   `[last_lower_bound_reduced_idx, last_upper_bound_reduced_idx]` such that
//!   `dt_min ≤ t[last] - t[first] ≤ dt_max`.
//! - The initial `middle` is `first + 1` and the initial `last` is
//!   `max(last_lower_bound_reduced_idx, middle + 1)`.
//!
//! ## Complexity
//! - Time: typically ~`O(n²)` (one window sweep per anchor).
//! - Space: `O(1)` per yielded triplet (the generator holds only the reduced times and mapping).
//!
//! ## Notes & pitfalls
//! - `dt_min`/`dt_max` must be in the **same time units** as your observation epochs.
//! - If `dt_min > dt_max` or there are fewer than 3 reduced observations, the iterator is empty.
//! - The generator **does not** impose any ordering by a heuristic (e.g., “optimal interval”);
//!   it only ensures time-feasibility. Best-first strategies should be layered on top.
//!
//! ## See also
//! - A Gauss triplet generator that yields `GaussObs` and can be combined with Monte-Carlo
//!   perturbations (`realizations_iter`).
//! - An IOD search routine (e.g., `estimate_best_orbit`) that consumes indices for early-stop.

use crate::observations::{triplets_iod::downsample_uniform_with_edges_indices, Observations};

/// Stream-only generator of **reduced indices** `(first, middle, last)`
/// for time-feasible IOD triplets.
///
/// Arguments
/// -----------------
/// * `dt_min` – Minimum allowed time span between **first** and **last**.
/// * `dt_max` – Maximum allowed time span between **first** and **last**.
/// * `max_triplets_to_yield` – Optional cap on the number of yielded triplets (use `usize::MAX` for no cap).
///
/// Return
/// ----------
/// * Implements `Iterator<Item = (usize, usize, usize)>` where the tuple contains
///   **reduced** indices `(first, middle, last)`.
///
/// See also
/// ------------
/// * [`selected_original_indices`](TripletIndexGenerator::selected_original_indices) – Map reduced indices back to original ones.
/// * A `GaussObs`-level generator if you need full triplets instead of indices.
pub struct TripletIndexGenerator {
    /// Map reduced index → original index (owned).
    reduced_to_original_index: Vec<usize>,
    /// Epochs of the reduced observations (same units as input times; owned).
    reduced_epochs_tt_mjd: Vec<f64>,

    /// Current **first** (anchor) index in reduced space.
    first_reduced_idx: usize,
    /// Current **middle** index in reduced space.
    middle_reduced_idx: usize,
    /// Current **last** index in reduced space.
    last_reduced_idx: usize,

    /// Lower bound (inclusive) for `last` given the current `first`.
    last_lower_bound_reduced_idx: usize,
    /// Upper bound (inclusive) for `last` given the current `first`.
    last_upper_bound_reduced_idx: usize,

    /// Number of reduced observations.
    reduced_len: usize,

    /// Time-window constraints on `(first, last)`.
    dt_min: f64,
    dt_max: f64,

    /// Count of triplets yielded so far (monotonic).
    yielded_triplets_count: usize,
    /// Hard cap on the number of triplets to yield.
    max_triplets_to_yield: usize,
}

impl TripletIndexGenerator {
    /// Build a generator from a full observation set:
    /// - Sorts by time (in-place),
    /// - Downsamples to at most `max_obs_for_triplets`,
    /// - Caches reduced epochs and the reduced→original mapping (both **owned**),
    /// - Positions on the first feasible window if any.
    ///
    /// Arguments
    /// -----------------
    /// * `observations` – The full set; will be **sorted by time in place**.
    /// * `dt_min`, `dt_max` – Time-span constraints on `(first, last)`.
    /// * `max_obs_for_triplets` – Downsampling cap (uniform with edges).
    /// * `max_triplets_to_yield` – Upper bound on yielded triplets (`usize::MAX` for no cap).
    ///
    /// Return
    /// ----------
    /// * A `TripletIndexGenerator` positioned at the first feasible window,
    ///   or “empty” if fewer than 3 reduced observations remain.
    ///
    /// See also
    /// ------------
    /// * [`TripletIndexGenerator::selected_original_indices`]
    /// * [`TripletIndexGenerator::reduced_times`]
    pub fn from_observations(
        observations: &mut Observations,
        dt_min: f64,
        dt_max: f64,
        max_obs_for_triplets: usize,
        max_triplets_to_yield: usize,
    ) -> Self {
        // 1) Sort by epoch (ascending)
        observations.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

        // 2) Downsample → keep indices only (no long borrows)
        let reduced_to_original_index =
            downsample_uniform_with_edges_indices(observations.len(), max_obs_for_triplets);

        // 3) Cache reduced epochs (owned) aligned with reduced indices
        let reduced_epochs_tt_mjd: Vec<f64> = reduced_to_original_index
            .iter()
            .map(|&orig| observations[orig].time)
            .collect();

        let reduced_len = reduced_epochs_tt_mjd.len();

        // Initialize; the precise window is set by `init_last_window_for_first`
        let mut gen = Self {
            reduced_to_original_index,
            reduced_epochs_tt_mjd,
            first_reduced_idx: 0,
            middle_reduced_idx: 1,
            last_reduced_idx: 2,
            last_lower_bound_reduced_idx: 2,
            last_upper_bound_reduced_idx: 1, // sentinel; will be recomputed
            reduced_len,
            dt_min,
            dt_max,
            yielded_triplets_count: 0,
            max_triplets_to_yield,
        };

        if gen.reduced_len >= 3 {
            gen.init_last_window_for_first();
        }
        gen
    }

    /// Access the reduced→original index mapping.
    ///
    /// Return
    /// ----------
    /// * A slice of original indices; `selected_original_indices()[r]` maps reduced `r` → original.
    ///
    /// See also
    /// ------------
    /// * [`TripletIndexGenerator::reduced_times`]
    pub fn selected_original_indices(&self) -> &[usize] {
        &self.reduced_to_original_index
    }

    /// Access the reduced epochs (TT/MJD).
    ///
    /// Return
    /// ----------
    /// * A slice of epochs aligned with reduced indices.
    pub fn reduced_times(&self) -> &[f64] {
        &self.reduced_epochs_tt_mjd
    }

    /// Recompute the feasible `last` window `[lower, upper]` for the current `first`.
    ///
    /// Invariants
    /// -----------------
    /// * `lower` is the earliest `last` with `t[last] - t[first] ≥ dt_min`.
    /// * `upper` is the latest  `last` with `t[last] - t[first] ≤ dt_max`.
    /// * `middle = first + 1`.
    /// * `last   = max(lower, middle + 1)`.
    fn init_last_window_for_first(&mut self) {
        let first = self.first_reduced_idx;

        // Lower bound (earliest last satisfying the min span; need one middle in (first, last))
        let mut lower = first + 2;
        while lower < self.reduced_len
            && (self.reduced_epochs_tt_mjd[lower] - self.reduced_epochs_tt_mjd[first]) < self.dt_min
        {
            lower += 1;
        }

        // Upper bound (latest last satisfying the max span)
        let mut upper = lower.saturating_sub(1).max(first + 1);
        while (upper + 1) < self.reduced_len
            && (self.reduced_epochs_tt_mjd[upper + 1] - self.reduced_epochs_tt_mjd[first])
                <= self.dt_max
        {
            upper += 1;
        }

        self.last_lower_bound_reduced_idx = lower;
        self.last_upper_bound_reduced_idx = upper;

        self.middle_reduced_idx = first + 1;
        self.last_reduced_idx = self
            .last_lower_bound_reduced_idx
            .max(self.middle_reduced_idx + 1);
    }

    /// Move to the next `first` anchor and refresh its feasible `last` window.
    ///
    /// Return
    /// ----------
    /// * `true` if there are still enough reduced observations to form a triplet; `false` otherwise.
    fn advance_first_anchor(&mut self) -> bool {
        self.first_reduced_idx += 1;
        if self.first_reduced_idx + 2 >= self.reduced_len {
            return false;
        }
        self.init_last_window_for_first();
        true
    }

    /// Whether the current `(first, middle, last)` window is empty.
    ///
    /// Return
    /// ----------
    /// * `true` if the time-feasible window is empty or invalid; `false` otherwise.
    fn last_window_is_empty(&self) -> bool {
        self.last_lower_bound_reduced_idx >= self.reduced_len
            || self.last_lower_bound_reduced_idx <= self.first_reduced_idx + 1
            || self.last_upper_bound_reduced_idx <= self.first_reduced_idx + 1
            || self.last_lower_bound_reduced_idx > self.last_upper_bound_reduced_idx
    }
}

impl Iterator for TripletIndexGenerator {
    type Item = (usize, usize, usize); // reduced indices (first, middle, last)

    fn next(&mut self) -> Option<Self::Item> {
        // Respect the optional global cap.
        if self.yielded_triplets_count >= self.max_triplets_to_yield {
            return None;
        }

        while self.first_reduced_idx + 2 < self.reduced_len {
            // If the current last-window is empty, move to the next anchor.
            if self.last_window_is_empty() {
                if !self.advance_first_anchor() {
                    return None;
                }
                continue;
            }

            // If `middle` reached the upper bound, switch to the next `first`.
            if self.middle_reduced_idx >= self.last_upper_bound_reduced_idx {
                if !self.advance_first_anchor() {
                    return None;
                }
                continue;
            }

            // Ensure `last` is within the feasible window and respects `middle < last`.
            if self.last_reduced_idx < self.last_lower_bound_reduced_idx
                || self.last_reduced_idx <= self.middle_reduced_idx
            {
                self.last_reduced_idx = self
                    .last_lower_bound_reduced_idx
                    .max(self.middle_reduced_idx + 1);
            }

            // If `last` exceeded the window, advance `middle` and reset `last`.
            if self.last_reduced_idx > self.last_upper_bound_reduced_idx {
                self.middle_reduced_idx += 1;
                self.last_reduced_idx = self
                    .last_lower_bound_reduced_idx
                    .max(self.middle_reduced_idx + 1);
                continue;
            }

            // We have a feasible triplet (first, middle, last).
            let i = self.first_reduced_idx;
            let j = self.middle_reduced_idx;
            let k = self.last_reduced_idx;

            // Prepare the next candidate by advancing `last`.
            self.last_reduced_idx += 1;

            self.yielded_triplets_count += 1;
            return Some((i, j, k));
        }

        // No more anchors → enumeration complete.
        None
    }
}
