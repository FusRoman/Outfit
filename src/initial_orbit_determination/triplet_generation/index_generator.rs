//! # IOD Triplet Index Generator (lazy, windowed by time span)
//!
//! Streams index triplets `(i, j, k)` — called **anchor**, **middle**, **last** — that
//! satisfy a time-span constraint on the outer pair:
//!
//! $$dt\_min \leq t_k - t_i \leq dt\_max \quad \text{with } i < j < k$$
//!
//! Indices refer to a *downsampled* view of the input observations (the
//! **reduced set**). They map directly to positions in the input slice.
//!
//! ## Input contract
//!
//! The observation slice passed to [`TripletIndexGenerator::from_observations`]
//! **must already be sorted in ascending time order**. No internal sorting is
//! performed. Violating this contract produces silently incorrect triplets.
//!
//! ## Algorithm
//!
//! For each anchor `i`, a two-pointer sweep finds the valid window
//! $[lo, hi]$ for `k` such that the time-span constraint holds.
//! `j` then ranges over `(i, hi)` and `k` over `(j, hi]`.
//! This gives ~$O(n^2)$ complexity versus $O(n^3)$ brute-force.
//!
//! ## Typical usage
//!
//! ```text
//! // observations must be sorted by ascending epoch before this call.
//! let gen = TripletIndexGenerator::from_observations(&obs, dt_min, dt_max, 200, usize::MAX);
//! for (i, j, k) in gen {
//!     // i, j, k index directly into the (downsampled view of the) input slice.
//! }
//! ```
//!
//! ## Notes
//! - `dt_min`/`dt_max` must share the same time unit as the observation epochs (TT/MJD).
//! - Fewer than 3 reduced observations or `dt_min > dt_max` → empty iterator.
//! - No ordering by heuristic is imposed; layer a best-K heap on top if needed.

use photom::observation_dataset::observation::Observation;

// ---------------------------------------------------------------------------
// Downsampling
// ---------------------------------------------------------------------------

/// Select `max_keep` indices from `0..n` uniformly, always including `0` and `n-1`.
///
/// This reduces the $O(n^3)$ triplet explosion while preserving the full time span.
///
/// Behavior
/// -----------------
/// | Condition          | Result                              |
/// |--------------------|-------------------------------------|
/// | `n == 0`           | `[]`                                |
/// | `max_keep >= n`    | `[0, 1, …, n-1]` (identity)        |
/// | `max_keep <= 3`    | `[0, n/2, n-1]`                    |
/// | otherwise          | `max_keep` uniformly spaced indices |
///
/// Arguments
/// -----------------
/// * `n`        – Total number of points.
/// * `max_keep` – Maximum number of indices to return.
///
/// Return
/// ----------
/// Indices in **strictly ascending** order.
pub(crate) fn downsample_uniform_with_edges(n: usize, max_keep: usize) -> Vec<usize> {
    match n {
        0 => vec![],
        _ if max_keep >= n => (0..n).collect(),
        _ if max_keep <= 3 => vec![0, n / 2, n - 1],
        _ => (0..max_keep)
            .map(|i| i * (n - 1) / (max_keep - 1))
            .collect(),
    }
}

// ---------------------------------------------------------------------------
// Feasible window for a fixed anchor
// ---------------------------------------------------------------------------

/// Valid range of `last` indices for a fixed anchor `i`.
///
/// `lo` is the smallest index `k > i` with $t_k - t_i \geq dt\_min$.
/// `hi` is the largest  index `k > i` with $t_k - t_i \leq dt\_max$.
///
/// The window is **empty** when `lo > hi` or no such `k` exists.
#[derive(Debug, Clone, Copy)]
struct LastWindow {
    lo: usize,
    hi: usize,
}

impl LastWindow {
    fn compute(anchor: usize, epochs: &[f64], dt_min: f64, dt_max: f64) -> Self {
        let n = epochs.len();
        let t0 = epochs[anchor];

        let mut lo = anchor + 2;
        while lo < n && epochs[lo] - t0 < dt_min {
            lo += 1;
        }

        let mut hi = lo.saturating_sub(1).max(anchor + 1);
        while hi + 1 < n && epochs[hi + 1] - t0 <= dt_max {
            hi += 1;
        }

        Self { lo, hi }
    }

    fn is_empty(&self, anchor: usize, n: usize) -> bool {
        self.lo >= n || self.lo > self.hi || self.hi <= anchor + 1
    }
}

// ---------------------------------------------------------------------------
// TripletIndexGenerator
// ---------------------------------------------------------------------------

/// Lazy iterator over time-feasible IOD triplet indices `(anchor, middle, last)`.
///
/// # Input contract
///
/// The observation slice passed to [`from_observations`](Self::from_observations)
/// **must be sorted in ascending time order** before construction.
///
/// # Index space
///
/// Indices yielded by the iterator refer directly to positions in the
/// (downsampled) input slice — no remapping is needed.
///
/// See the [module-level documentation](self) for the algorithm and usage.
pub struct TripletIndexGenerator {
    /// Epochs of the reduced set (TT/MJD), extracted from the downsampled positions.
    epochs: Vec<f64>,

    // --- iteration state ---
    anchor: usize,
    middle: usize,
    last: usize,
    window: LastWindow,

    n: usize,
    dt_min: f64,
    dt_max: f64,

    /// Remaining triplets allowed before the iterator stops (`usize::MAX` = no cap).
    remaining: usize,
}

impl TripletIndexGenerator {
    /// Construct directly from a reduced epoch vector.
    ///
    /// `epochs` must be in **ascending** order.
    ///
    /// Arguments
    /// -----------------
    /// * `epochs`           – Reduced epochs in ascending order.
    /// * `dt_min`, `dt_max` – Time-span bounds on `(anchor, last)`.
    /// * `cap`              – Maximum triplets to yield (`usize::MAX` for no limit).
    pub fn new(epochs: Vec<f64>, dt_min: f64, dt_max: f64, cap: usize) -> Self {
        let n = epochs.len();
        let window = if n >= 3 {
            LastWindow::compute(0, &epochs, dt_min, dt_max)
        } else {
            LastWindow { lo: n, hi: 0 }
        };

        Self {
            n,
            epochs,
            dt_min,
            dt_max,
            anchor: 0,
            middle: 1,
            last: window.lo.max(2),
            window,
            remaining: cap,
        }
    }

    /// Build from a time-sorted observation slice, with optional downsampling.
    ///
    /// The input slice **must already be sorted in ascending time order**.
    /// Downsampling via `downsample_uniform_with_edges` is applied, always
    /// preserving the first and last observations.
    ///
    /// Arguments
    /// -----------------
    /// * `observations`      – Time-sorted observation slice (ascending epoch).
    /// * `dt_min`, `dt_max`  – Time-span bounds on `(anchor, last)`.
    /// * `max_reduced`       – Downsampling cap (uniform with endpoints).
    /// * `cap`               – Maximum triplets to yield (`usize::MAX` for no limit).
    pub fn from_observations(
        observations: &[&Observation],
        dt_min: f64,
        dt_max: f64,
        max_reduced: usize,
        cap: usize,
    ) -> Self {
        let keep = downsample_uniform_with_edges(observations.len(), max_reduced);
        let epochs: Vec<f64> = keep.iter().map(|&i| observations[i].mjd_tt()).collect();
        Self::new(epochs, dt_min, dt_max, cap)
    }

    /// Reduced epochs (TT/MJD), aligned with the indices yielded by the iterator.
    pub fn reduced_times(&self) -> &[f64] {
        &self.epochs
    }

    fn advance_anchor(&mut self) -> bool {
        self.anchor += 1;
        if self.anchor + 2 >= self.n {
            return false;
        }
        self.window = LastWindow::compute(self.anchor, &self.epochs, self.dt_min, self.dt_max);
        self.middle = self.anchor + 1;
        self.last = self.window.lo.max(self.middle + 1);
        true
    }

    #[inline]
    fn reset_last_for_middle(&mut self) {
        self.last = self.window.lo.max(self.middle + 1);
    }
}

impl Iterator for TripletIndexGenerator {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        if self.remaining == 0 {
            return None;
        }

        loop {
            if self.anchor + 2 >= self.n {
                return None;
            }

            if self.window.is_empty(self.anchor, self.n) {
                if !self.advance_anchor() {
                    return None;
                }
                continue;
            }

            if self.middle >= self.window.hi {
                if !self.advance_anchor() {
                    return None;
                }
                continue;
            }

            if self.last <= self.middle {
                self.reset_last_for_middle();
            }

            if self.last > self.window.hi {
                self.middle += 1;
                self.reset_last_for_middle();
                continue;
            }

            let triplet = (self.anchor, self.middle, self.last);
            self.last += 1;
            self.remaining -= 1;
            return Some(triplet);
        }
    }
}

#[cfg(test)]
mod triplet_generator_tests {
    use super::*;
    use proptest::prelude::*;

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------

    /// Build a generator directly from a sorted epoch slice (no observations needed).
    fn gen_from_epochs(epochs: Vec<f64>, dt_min: f64, dt_max: f64) -> TripletIndexGenerator {
        TripletIndexGenerator::new(epochs, dt_min, dt_max, usize::MAX)
    }

    /// Collect all triplets and verify every invariant inline.
    fn collect_and_validate(
        epochs: &[f64],
        dt_min: f64,
        dt_max: f64,
    ) -> Vec<(usize, usize, usize)> {
        let gen = gen_from_epochs(epochs.to_vec(), dt_min, dt_max);
        let mut out = Vec::new();
        for (i, j, k) in gen {
            assert!(i < j, "first < middle violated: ({i},{j},{k})");
            assert!(j < k, "middle < last violated: ({i},{j},{k})");
            let span = epochs[k] - epochs[i];
            assert!(
                span >= dt_min - 1e-12,
                "span {span} < dt_min {dt_min}: ({i},{j},{k})"
            );
            assert!(
                span <= dt_max + 1e-12,
                "span {span} > dt_max {dt_max}: ({i},{j},{k})"
            );
            out.push((i, j, k));
        }
        out
    }

    /// Brute-force reference: all (i,j,k) with i<j<k and dt_min ≤ t[k]-t[i] ≤ dt_max.
    fn brute_force(epochs: &[f64], dt_min: f64, dt_max: f64) -> Vec<(usize, usize, usize)> {
        let n = epochs.len();
        let mut out = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    let span = epochs[k] - epochs[i];
                    if span >= dt_min - 1e-12 && span <= dt_max + 1e-12 {
                        out.push((i, j, k));
                    }
                }
            }
        }
        out
    }

    // -----------------------------------------------------------------------
    // downsample_uniform_with_edges_indices
    // -----------------------------------------------------------------------

    #[test]
    fn downsample_empty() {
        assert!(downsample_uniform_with_edges(0, 10).is_empty());
    }

    #[test]
    fn downsample_no_op_when_max_ge_n() {
        let result = downsample_uniform_with_edges(5, 10);
        assert_eq!(result, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn downsample_exact_n() {
        let result = downsample_uniform_with_edges(5, 5);
        assert_eq!(result, vec![0, 1, 2, 3, 4]);
    }

    #[test]
    fn downsample_max_keep_3() {
        // Always returns [0, mid, n-1].
        let result = downsample_uniform_with_edges(9, 3);
        assert_eq!(result, vec![0, 4, 8]);
    }

    #[test]
    fn downsample_max_keep_le_3_small_n() {
        let result = downsample_uniform_with_edges(3, 2);
        // Edge case: max_keep ≤ 3 branch → [0, mid, n-1] = [0, 1, 2]
        assert_eq!(result[0], 0);
        assert_eq!(*result.last().unwrap(), 2);
    }

    #[test]
    fn downsample_endpoints_always_present() {
        for n in 4..=20 {
            for max_keep in 3..=n {
                let result = downsample_uniform_with_edges(n, max_keep);
                assert_eq!(
                    result[0], 0,
                    "first endpoint missing for n={n} max={max_keep}"
                );
                assert_eq!(
                    *result.last().unwrap(),
                    n - 1,
                    "last endpoint missing for n={n} max={max_keep}"
                );
            }
        }
    }

    #[test]
    fn downsample_length_respects_max_keep() {
        for n in 4..=30 {
            for max_keep in 3..n {
                let result = downsample_uniform_with_edges(n, max_keep);
                assert!(
                    result.len() <= max_keep,
                    "len={} > max_keep={max_keep} for n={n}",
                    result.len()
                );
            }
        }
    }

    #[test]
    fn downsample_strictly_increasing() {
        let result = downsample_uniform_with_edges(100, 10);
        for w in result.windows(2) {
            assert!(w[0] < w[1], "not strictly increasing: {:?}", result);
        }
    }

    // -----------------------------------------------------------------------
    // TripletIndexGenerator — edge cases
    // -----------------------------------------------------------------------

    #[test]
    fn generator_empty_on_fewer_than_3_obs() {
        for n in 0..=2 {
            let epochs: Vec<f64> = (0..n).map(|i| i as f64).collect();
            let triplets = collect_and_validate(&epochs, 0.0, 10.0);
            assert!(triplets.is_empty(), "expected empty for n={n}");
        }
    }

    #[test]
    fn generator_empty_when_dt_min_gt_dt_max() {
        let epochs = vec![0.0, 1.0, 2.0, 3.0];
        let triplets = collect_and_validate(&epochs, 5.0, 2.0);
        assert!(triplets.is_empty());
    }

    #[test]
    fn generator_empty_when_all_spans_below_dt_min() {
        // All spans ≤ 2, dt_min = 10.
        let epochs = vec![0.0, 1.0, 2.0, 3.0];
        let triplets = collect_and_validate(&epochs, 10.0, 100.0);
        assert!(triplets.is_empty());
    }

    #[test]
    fn generator_empty_when_all_spans_above_dt_max() {
        // All spans ≥ 3, dt_max = 1.
        let epochs = vec![0.0, 10.0, 20.0, 30.0];
        let triplets = collect_and_validate(&epochs, 0.0, 1.0);
        assert!(triplets.is_empty());
    }

    #[test]
    fn generator_single_feasible_triplet() {
        // Only (0,1,2) has span 2.0 ∈ [2, 2].
        let epochs = vec![0.0, 1.0, 2.0];
        let triplets = collect_and_validate(&epochs, 2.0, 2.0);
        assert_eq!(triplets, vec![(0, 1, 2)]);
    }

    #[test]
    fn generator_matches_brute_force_small() {
        let epochs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let dt_min = 1.5;
        let dt_max = 3.5;

        let mut got = collect_and_validate(&epochs, dt_min, dt_max);
        let mut expected = brute_force(&epochs, dt_min, dt_max);

        got.sort();
        expected.sort();
        assert_eq!(got, expected);
    }

    #[test]
    fn generator_matches_brute_force_no_constraint() {
        // dt_min = 0, dt_max = ∞ → all (i<j<k) must appear.
        let epochs: Vec<f64> = (0..6).map(|i| i as f64).collect();

        let mut got = collect_and_validate(&epochs, 0.0, f64::MAX);
        let mut expected = brute_force(&epochs, 0.0, f64::MAX);

        got.sort();
        expected.sort();
        assert_eq!(got, expected);
    }

    #[test]
    fn generator_matches_brute_force_equal_spacing() {
        let epochs: Vec<f64> = (0..7).map(|i| i as f64 * 2.0).collect();
        let dt_min = 3.0;
        let dt_max = 9.0;

        let mut got = collect_and_validate(&epochs, dt_min, dt_max);
        let mut expected = brute_force(&epochs, dt_min, dt_max);

        got.sort();
        expected.sort();
        assert_eq!(got, expected);
    }

    #[test]
    fn generator_no_duplicates() {
        let epochs: Vec<f64> = (0..8).map(|i| i as f64).collect();
        let mut triplets = collect_and_validate(&epochs, 1.0, 6.0);
        triplets.sort();
        triplets.dedup();
        let all = collect_and_validate(&epochs, 1.0, 6.0);
        assert_eq!(triplets.len(), all.len(), "duplicates detected");
    }

    // -----------------------------------------------------------------------
    // max_triplets_to_yield cap
    // -----------------------------------------------------------------------

    #[test]
    fn generator_respects_max_triplets_cap() {
        let epochs: Vec<f64> = (0..10).map(|i| i as f64).collect();
        let cap = 5;
        let gen = TripletIndexGenerator::new(epochs, 1.0, 20.0, cap);
        let count = gen.count();
        assert_eq!(count, cap);
    }

    #[test]
    fn generator_cap_zero_yields_nothing() {
        let epochs = vec![0.0, 1.0, 2.0, 3.0];
        let gen = TripletIndexGenerator::new(epochs, 0.0, 10.0, 0);
        assert_eq!(gen.count(), 0);
    }

    // -----------------------------------------------------------------------
    // reduced_to_original mapping
    // -----------------------------------------------------------------------

    #[test]
    fn reduced_times_match_input_epochs() {
        let epochs: Vec<f64> = (0..5).map(|i| i as f64).collect();
        let gen = TripletIndexGenerator::new(epochs.clone(), 0.0, 10.0, usize::MAX);
        assert_eq!(gen.reduced_times(), epochs.as_slice());
    }

    #[test]
    fn reduced_times_aligned_with_mapping() {
        let epochs = vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let gen = TripletIndexGenerator::new(epochs.clone(), 0.0, 20.0, usize::MAX);

        // reduced_times length must match the number of epochs provided.
        assert_eq!(gen.reduced_times().len(), epochs.len());
    }

    // -----------------------------------------------------------------------
    // EpochSortedIndices
    // -----------------------------------------------------------------------

    // (tested indirectly through from_observations; direct access is private)

    // -----------------------------------------------------------------------
    // Proptest: invariants hold for random inputs
    // -----------------------------------------------------------------------

    proptest! {
        /// For random sorted epoch vectors and random dt windows, every yielded
        /// triplet must satisfy i < j < k and dt_min ≤ t[k]-t[i] ≤ dt_max.
        #[test]
        fn prop_all_invariants_hold(
            // Generate 3..=12 strictly increasing epochs starting near 0.
            n in 3usize..=12,
            steps in prop::collection::vec(0.1f64..5.0, 11),
            dt_min in 0.0f64..10.0,
            dt_range in 0.0f64..20.0,
        ) {
            let dt_max = dt_min + dt_range;
            // Build epochs as cumulative sum of steps (strictly increasing).
            let mut epochs = vec![0.0f64];
            for &s in steps.iter().take(n - 1) {
                epochs.push(epochs.last().unwrap() + s);
            }

            let _triplets = collect_and_validate(&epochs, dt_min, dt_max);
            // collect_and_validate asserts invariants inline.
        }

        /// Generator output matches brute-force for small random inputs.
        #[test]
        fn prop_matches_brute_force(
            n in 3usize..=8,
            steps in prop::collection::vec(0.5f64..3.0, 7),
            dt_min in 0.0f64..5.0,
            dt_range in 0.0f64..10.0,
        ) {
            let dt_max = dt_min + dt_range;
            let mut epochs = vec![0.0f64];
            for &s in steps.iter().take(n - 1) {
                epochs.push(epochs.last().unwrap() + s);
            }

            let mut got = collect_and_validate(&epochs, dt_min, dt_max);
            let mut expected = brute_force(&epochs, dt_min, dt_max);
            got.sort();
            expected.sort();
            prop_assert_eq!(got, expected);
        }

        /// No duplicate triplets regardless of input.
        #[test]
        fn prop_no_duplicates(
            n in 3usize..=10,
            steps in prop::collection::vec(0.1f64..4.0, 9),
            dt_min in 0.0f64..5.0,
            dt_range in 0.0f64..15.0,
        ) {
            let dt_max = dt_min + dt_range;
            let mut epochs = vec![0.0f64];
            for &s in steps.iter().take(n - 1) {
                epochs.push(epochs.last().unwrap() + s);
            }

            let mut triplets = collect_and_validate(&epochs, dt_min, dt_max);
            let total = triplets.len();
            triplets.sort();
            triplets.dedup();
            prop_assert_eq!(triplets.len(), total, "duplicate triplets found");
        }

        /// Cap is always respected.
        #[test]
        fn prop_cap_respected(
            n in 3usize..=12,
            steps in prop::collection::vec(0.5f64..3.0, 11),
            cap in 0usize..=20,
        ) {
            let mut epochs = vec![0.0f64];
            for &s in steps.iter().take(n - 1) {
                epochs.push(epochs.last().unwrap() + s);
            }
            let gen = TripletIndexGenerator::new(
                epochs, 0.0, f64::MAX, cap,
            );
            prop_assert!(gen.count() <= cap);
        }

        /// downsample always returns endpoints and respects the length cap.
        #[test]
        fn prop_downsample_endpoints_and_length(
            n in 1usize..=50,
            max_keep in 3usize..=50,
        ) {
            let result = downsample_uniform_with_edges(n, max_keep);
            prop_assert!(result.len() <= max_keep.min(n).max(3));
            if n >= 1 {
                prop_assert_eq!(result[0], 0);
                prop_assert_eq!(*result.last().unwrap(), n - 1);
            }
        }

        /// downsample output is strictly increasing.
        #[test]
        fn prop_downsample_strictly_increasing(
            n in 2usize..=50,
            max_keep in 3usize..=50,
        ) {
            let result = downsample_uniform_with_edges(n, max_keep);
            for w in result.windows(2) {
                prop_assert!(w[0] < w[1]);
            }
        }
    }
}
