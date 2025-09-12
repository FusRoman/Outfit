//! # Triplet selection for Gauss Initial Orbit Determination (IOD)
//!
//! Utilities to **thin**, **enumerate**, and **rank** triplets of astrometric
//! observations prior to running the Gauss IOD solver.
//!
//! ## What this module does
//!
//! - **Downsample** dense time series to a representative subset while keeping endpoints,
//! - **Generate** all valid triplets under configurable time-span constraints,
//! - **Score** triplets with a spacing-based weight that favors near-uniform timing,
//! - **Select** only the best candidates using a heap-based pruning strategy.
//!
//! The output is a small set of [`GaussObs`] triplets, ready to be passed to
//! the Gauss preliminary-orbit routine.
//!
//! ## Workflow overview
//!
//! 1. **Sort in time** (in-place).
//! 2. **Downsample** to at most *N* points with uniform coverage (always keep first/last).
//! 3. **Enumerate** triplets `(i, j, k)` such that `dt_min ≤ t_k − t_i ≤ dt_max`.
//! 4. **Score** each triplet with [`triplet_weight`], using `optimal_interval_time` as target spacing.
//! 5. **Keep top-K** (lowest weight) using a `BinaryHeap` configured as a **max-heap**.
//! 6. **Materialize** the winners as [`GaussObs`] with precomputed observer positions.
//!
//! ## Complexity & performance
//!
//! - Naïve enumeration scales as **O(n³)**.  
//!   Downsampling plus early pruning reduces candidates dramatically for large `n`.
//! - The heap keeps memory bounded to `O(max_triplet)`.
//! - Weighting is **O(1)** per triplet.
//!
//! ## Configuration knobs (typical values)
//!
//! - `dt_min` *(days)*: **0.03** — minimum triplet span,
//! - `dt_max` *(days)*: **150.0** — maximum triplet span,
//! - `optimal_interval_time` *(days)*: **20.0** — ideal spacing for `(t₂−t₁)` and `(t₃−t₂)`,
//! - `max_obs_for_triplets`: **100** — post-thinning observation cap,
//! - `max_triplet`: **10** — number of best triplets to return.
//!
//! ## Example
//!
//! ```rust,no_run
//! use outfit::constants::Observations;
//! use outfit::observations::triplets_iod::generate_triplets;
//!
//! // Your observations loaded elsewhere:
//! let mut obs: Observations = unimplemented!();
//!
//! // Build best triplets for Gauss IOD
//! let triplets = generate_triplets(
//!     &mut obs,
//!     0.03,   // dt_min [days]
//!     150.0,  // dt_max [days]
//!     20.0,   // optimal_interval_time [days]
//!     100,    // max_obs_for_triplets
//!     10      // max_triplet
//! );
//!
//! // Now pass `triplets` to your Gauss solver…
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## Guarantees & edge cases
//!
//! - Endpoints are **always preserved** by the downsampler.
//! - If the input has fewer than three observations, triplet generation yields **no result**.
//! - Triplet weights are strictly **lower-is-better**; ties are allowed.
//!
//! ## See also
//!
//! - [`generate_triplets`] – Core function assembling and ranking triplets.
//! - [`triplet_weight`] – Spacing-based scoring rule.
//! - [`GaussObs`] – Triplet container consumed by the Gauss IOD solver.
//! - [`crate::observations::observations_ext::ObservationsExt::compute_triplets`] – Higher-level wrapper.
use nalgebra::{Matrix3, Vector3};
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::constants::Observations;
use crate::initial_orbit_determination::gauss::GaussObs;
use crate::observations::triplets_generator::TripletIndexGenerator;
use crate::observations::Observation;

/// Internal structure holding a weighted observation triplet during selection.
///
/// This container stores the indices `(i, j, k)` of a candidate triplet together
/// with its scalar `weight`. The type implements ordering so it can be placed in a
/// `BinaryHeap` acting as a **max-heap** on `weight`: the *worst* (largest weight)
/// candidate sits at the top and is pruned first.
///
/// Fields
/// -----------------
/// * `weight` – Triplet score (lower is better).
/// * `i`, `j`, `k` – Indices of the first, second, and third observation.
///
/// Remarks
/// -----------------
/// * Ordering treats larger `weight` as “greater”, enabling a max-heap.
/// * `NaN` comparisons are mapped to `Ordering::Equal` (see `Ord` impl),
///   which effectively groups `NaN` weights and prevents panics. Avoid `NaN`
///   weights in normal operation.
///
/// See also
/// ------------
/// * [`generate_triplets`] – Produces and ranks triplets using this structure.
/// * [`triplet_weight`] – Scoring function used to produce `weight`.
#[derive(Clone, Copy, Debug)]
struct WeightedTriplet {
    weight: f64,
    first_idx: usize,
    middle_idx: usize,
    last_idx: usize,
}

impl PartialEq for WeightedTriplet {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
            && self.first_idx == other.first_idx
            && self.middle_idx == other.middle_idx
            && self.last_idx == other.last_idx
    }
}
impl Eq for WeightedTriplet {}

impl PartialOrd for WeightedTriplet {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        // IMPORTANT: **do not** reverse. We want a normal max-heap by weight.
        Some(self.cmp(other))
    }
}
impl Ord for WeightedTriplet {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.weight.partial_cmp(&other.weight) {
            Some(ord) => ord,
            None => (self.first_idx, self.middle_idx, self.last_idx).cmp(&(
                other.first_idx,
                other.middle_idx,
                other.last_idx,
            )),
        }
    }
}

/// Compute a spacing-based weight for a time-ordered observation triplet.
///
/// The score penalizes **unbalanced** or **far-from-optimal** spacings between
/// consecutive observations. A lower score indicates a more desirable (nearly uniform)
/// spacing close to the target interval `dtw`.
///
/// Formula
/// -----------------
/// For consecutive gaps `dt₁ = t₂ − t₁`, `dt₂ = t₃ − t₂`, define:
/// ```text
/// s(dt, dtw) =  { dtw/dt,   if dt ≤ dtw
///               { 1 + dt/dtw, otherwise
/// ```
/// The total weight is:
/// ```text
/// weight = s(dt₁, dtw) + s(dt₂, dtw)
/// ```
/// This function computes the same value while precomputing `1/dtw` once
/// to reduce divisions.
///
/// Arguments
/// -----------------
/// * `time1` – Epoch of the first observation (must be finite).
/// * `time2` – Epoch of the middle observation (must be finite).
/// * `time3` – Epoch of the last observation (must be finite).
/// * `dtw`   – Target per-gap spacing (same unit as the epochs, strictly positive).
///
/// Return
/// ----------
/// * A scalar weight (`f64`) where **smaller is better** (ideal gaps ≈ `dtw`).
///
/// Remarks
/// -------------
/// * Assumes strictly increasing times (`time1 < time2 < time3`). If not guaranteed,
///   sort beforehand or validate upstream.
/// * Complexity: **O(1)**; uses one reciprocal of `dtw` and two cheap branches.
/// * Very small gaps (`dt → 0⁺`) yield large penalties (`∞` limit). If you prefer
///   to cap this behavior, clamp `dt` to an epsilon upstream (changes semantics).
///
/// See also
/// ------------
/// * [`triplet_weight_with_inv`] – Same logic but takes `inv_dtw` precomputed (hot loop friendly).
/// * [`TripletIndexGenerator`] – Streams feasible `(i, j, k)` triplets by time window.
/// * [`generate_triplets`] – Uses this weight to keep the best-K candidates.
#[inline(always)]
pub fn triplet_weight(time1: f64, time2: f64, time3: f64, dtw: f64) -> f64 {
    // Precompute reciprocal once (one division instead of two on average).
    let inv_dtw = dtw.recip();
    // Gaps are positive if times are sorted.
    let dt12 = time2 - time1;
    let dt23 = time3 - time2;
    s_gap(dt12, inv_dtw) + s_gap(dt23, inv_dtw)
}

/// Compute the same spacing-based weight as [`triplet_weight`], but taking
/// the reciprocal of the target spacing (`inv_dtw = 1.0 / dtw`) as input.
///
/// This variant is intended for **hot loops** where `dtw` is constant and can be
/// inverted once upstream, eliminating one division per call.
///
/// Arguments
/// -----------------
/// * `time1`   – Epoch of the first observation (finite).
/// * `time2`   – Epoch of the middle observation (finite).
/// * `time3`   – Epoch of the last observation (finite).
/// * `inv_dtw` – Precomputed reciprocal of the target spacing (`1.0 / dtw`, finite and positive).
///
/// Return
/// ----------
/// * The same scalar weight as [`triplet_weight`] (`f64`), **lower is better**.
///
/// Remarks
/// -------------
/// * Assumes `time1 < time2 < time3`.
/// * Prefer this function inside inner loops when `dtw` is constant across many triplets.
/// * Complexity: **O(1)**; avoids the extra division by reusing `inv_dtw`.
///
/// See also
/// ------------
/// * [`triplet_weight`] – Convenience version that accepts `dtw` directly.
#[inline(always)]
pub fn triplet_weight_with_inv(time1: f64, time2: f64, time3: f64, inv_dtw: f64) -> f64 {
    // Gaps are positive if times are sorted.
    let dt12 = time2 - time1;
    let dt23 = time3 - time2;
    s_gap(dt12, inv_dtw) + s_gap(dt23, inv_dtw)
}

/// Piecewise per-gap penalty used by the triplet spacing weight.
///
/// For the ratio `r = dt / dtw`, the penalty is:
/// ```text
/// s(dt, dtw) =  { 1 / r,     if r ≤ 1      (too small ⇢ penalize by inverse)
///               { 1 + r,     otherwise     (too large ⇢ linear growth)
/// ```
/// Passing `inv_dtw = 1 / dtw` avoids a division in the common path.
///
/// Arguments
/// -----------------
/// * `dt`      – The consecutive time gap (positive if epochs are strictly increasing).
/// * `inv_dtw` – Reciprocal of the target spacing (`1.0 / dtw`, positive).
///
/// Return
/// ----------
/// * The per-gap penalty `s(dt, dtw)` (`f64`).
///
/// Remarks
/// -------------
/// * Complexity: **O(1)**, one multiply and a single branch.
/// * As `dt → 0⁺`, `s(dt, dtw) → +∞`; clamp upstream if you need bounded penalties.
/// * Numeric domains must be finite; `NaN` or non-finite inputs propagate unspecified results.
///
/// See also
/// ------------
/// * [`triplet_weight`] – Sums two `s_gap` values for `(t₂−t₁)` and `(t₃−t₂)`.
/// * [`triplet_weight_with_inv`] – Same, using a precomputed `inv_dtw`.
#[inline(always)]
fn s_gap(dt: f64, inv_dtw: f64) -> f64 {
    // r = dt / dtw
    let r = dt * inv_dtw;
    if r <= 1.0 {
        r.recip()
    } else {
        // `1.0 + r` avoids an extra division.
        1.0 + r
    }
}

/// Downsample observation indices while preserving endpoints and temporal coverage.
///
/// If the input has more than `max_keep` points, this routine selects a subset of indices
/// **uniformly in time** while always including the **first** and **last** observation.
/// This reduces the `O(n³)` triplet explosion without losing global time-span information.
///
/// Behavior
/// -----------------
/// * If `n == 0`: returns `[]`.
/// * If `max_keep >= n`: returns all indices `0..n`.
/// * If `max_keep <= 3`: returns `[0, mid, n-1]` (with `mid = n/2`).
///   For `n < 3`, indices may repeat (callers should handle deduplication if needed).
/// * Otherwise: returns `max_keep` indices distributed uniformly between `1` and `n-2`,
///   plus the endpoints `0` and `n-1`.
///
/// Arguments
/// -----------------
/// * `n` – Total number of observations.
/// * `max_keep` – Maximum number of indices to return.
///
/// Return
/// ----------
/// * A `Vec<usize>` with the selected indices in **ascending** order.
///
/// Remarks
/// -------------
/// * Complexity: **O(max_keep)** after the trivial cases.
/// * The selection is **index-uniform** over the span `[1, n-2]`; if strictly
///   time-uniform selection is required, pre-sort observations by time first
///   (as done in [`generate_triplets`]).
///
/// See also
/// ------------
/// * [`generate_triplets`] – Uses this function before triplet enumeration.
pub(crate) fn downsample_uniform_with_edges_indices(n: usize, max_keep: usize) -> Vec<usize> {
    match n {
        0 => Vec::new(),
        _ if max_keep <= 3 => {
            let mid = n / 2;
            vec![0, mid, n - 1]
        }
        _ if max_keep >= n => (0..n).collect(),
        _ => {
            let slots = max_keep - 2;
            std::iter::once(0)
                .chain((0..slots).map(move |i| {
                    let fraction = (i + 1) as f64 / (slots + 1) as f64;
                    // distribute indices uniformly between 1 and n-2
                    1 + (fraction * (n - 2) as f64).floor() as usize
                }))
                .chain(std::iter::once(n - 1))
                .collect()
        }
    }
}

/// Generate and select **best-K** triplets of astrometric observations
/// for Gauss Initial Orbit Determination (IOD), using a **lazy index stream**
/// and a bounded **max-heap** on a spacing weight.
///
/// Overview
/// -----------------
/// This routine constructs good candidate triplets `(first, middle, last)` as inputs
/// to the **Gauss method** while avoiding the `O(n³)` blow-up:
///
/// 1. **Sort & downsample (in `TripletIndexGenerator`)** – Observations are sorted by epoch
///    and uniformly thinned (keeping endpoints) to at most `max_obs_for_triplets`.
/// 2. **Time-feasible enumeration (lazy)** – Indices `(i, j, k)` are streamed by
///    `TripletIndexGenerator`, constrained by:
///    `dt_min ≤ t[k] − t[i] ≤ dt_max` with `i < j < k`.
/// 3. **Weight scoring** – Each feasible triplet receives a weight via [`triplet_weight`],
///    favoring near-uniform spacing around `optimal_interval_time`.
/// 4. **Best-K selection** – A bounded **max-heap** retains only the `max_triplet`
///    lowest-weight candidates (the heap’s `peek()` is the current **worst**).
/// 5. **Materialization** – Only for the selected indices, we re-borrow `observations`
///    immutably and build [`GaussObs`] with precomputed observer heliocentric columns
///    (via [`Observation::get_observer_helio_position`]).
///
/// Design notes
/// -----------------
/// * Enumeration and scoring happen on **reduced indices** (owned by the generator),
///   so there are **no overlapping borrows** of `observations`.
/// * The function avoids cloning the generator’s internal buffers (times, mapping);
///   we read them through short-lived immutable borrows between `next()` calls.
/// * The final `Vec<GaussObs>` is **sorted by increasing weight** (best first).
///
/// Arguments
/// -----------------
/// * `observations` – Mutable set of astrometric observations; epochs are sorted **in-place**.
/// * `dt_min` – Minimum allowed time span (same units as `Observation::time`) between first and last.
/// * `dt_max` – Maximum allowed time span between first and last.
/// * `optimal_interval_time` – Target per-gap spacing used by [`triplet_weight`].
/// * `max_obs_for_triplets` – Downsampling cap (uniform with edges).
/// * `max_triplet` – Number `K` of best triplets to return (heap capacity).
///
/// Return
/// ----------
/// * A `Vec<GaussObs>` of length `≤ max_triplet`, sorted by **ascending** heuristic weight.
///
/// Complexity
/// -----------------
/// * Enumeration: typically ~`O(n²)` thanks to the per-anchor time window in
///   [`TripletIndexGenerator`].
/// * Selection: `O(n log K)` due to the bounded heap.
/// * Space: `O(1)` per yielded triplet during enumeration; only the final `K` are materialized.
///
/// See also
/// ------------
/// * [`TripletIndexGenerator`] – Streams time-feasible reduced indices lazily.
/// * [`triplet_weight`] – Heuristic favoring evenly spaced triplets around a target gap.
/// * [`GaussObs::realizations_iter`] – Lazy Monte-Carlo perturbations per triplet.
/// * [`ObservationsExt::compute_triplets`](crate::observations::observations_ext::ObservationsExt::compute_triplets) – Typical high-level wrapper.
pub fn generate_triplets(
    observations: &mut Observations,
    dt_min: f64,
    dt_max: f64,
    optimal_interval_time: f64,
    max_obs_for_triplets: usize,
    max_triplet: u32,
) -> Vec<GaussObs> {
    if max_triplet == 0 {
        return Vec::new();
    }

    // --- Phase 1: enumerate feasible reduced indices & keep best-K by weight (no &Observation borrows).
    let mut index_gen = TripletIndexGenerator::from_observations(
        observations,
        dt_min,
        dt_max,
        max_obs_for_triplets,
        usize::MAX, // scan all feasible triplets; the heap does best-K filtering
    );

    let k_cap = max_triplet as usize;
    let mut heap: BinaryHeap<WeightedTriplet> = BinaryHeap::with_capacity(k_cap.saturating_add(1));

    // Bounded push: maintain the K smallest weights in a BinaryHeap (max-heap).
    let mut push_best_k = |cand: WeightedTriplet| {
        if !cand.weight.is_finite() {
            return; // guard against NaN/Inf
        }
        if heap.len() < k_cap {
            heap.push(cand);
        } else if let Some(worst) = heap.peek() {
            if cand.weight < worst.weight {
                heap.pop();
                heap.push(cand);
            }
        }
    };

    // Consume the reduced-index stream. After each `next()`, take a short immutable
    // borrow of the times to compute the weight (no overlap with the next `next()`).

    let inv_dtw = optimal_interval_time.recip(); // precompute once
    while let Some((i, j, k)) = index_gen.next() {
        let times = index_gen.reduced_times();
        let w = triplet_weight_with_inv(times[i], times[j], times[k], inv_dtw);
        push_best_k(WeightedTriplet {
            weight: w,
            first_idx: i,
            middle_idx: j,
            last_idx: k,
        });
    }

    // Best-K by ascending weight.
    let mut best_reduced = heap.into_sorted_vec();
    best_reduced.sort_by(|a, b| a.weight.partial_cmp(&b.weight).unwrap()); // defensive

    // --- Phase 2: materialize GaussObs for the selected indices (immutable borrows now safe).
    let mapping = index_gen.selected_original_indices();

    best_reduced
        .into_iter()
        .map(|wt| {
            let (i, j, k) = (wt.first_idx, wt.middle_idx, wt.last_idx);

            // reduced → original indices
            let oi = mapping[i];
            let oj = mapping[j];
            let ok = mapping[k];

            // Immutable borrows of the original observations occur only here.
            let o1: &Observation = &observations[oi];
            let o2: &Observation = &observations[oj];
            let o3: &Observation = &observations[ok];

            // Observer 3×3 matrix (columns = heliocentric observer positions at each epoch).
            let observer_matrix: Matrix3<f64> = Matrix3::from_columns(&[
                o1.get_observer_helio_position(),
                o2.get_observer_helio_position(),
                o3.get_observer_helio_position(),
            ]);

            GaussObs::with_observer_position(
                Vector3::new(oi, oj, ok),
                Vector3::new(o1.ra, o2.ra, o3.ra),
                Vector3::new(o1.dec, o2.dec, o3.dec),
                Vector3::new(o1.time, o2.time, o3.time),
                observer_matrix,
            )
        })
        .collect()
}

#[cfg(test)]
mod triplets_iod_tests {

    #[cfg(feature = "jpl-download")]
    use approx::assert_relative_eq;

    use super::*;

    #[cfg(feature = "jpl-download")]
    pub(crate) fn assert_gauss_obs_approx_eq(a: &GaussObs, b: &GaussObs, tol: f64) {
        assert_eq!(a.idx_obs, b.idx_obs);
        assert_relative_eq!(a.ra, b.ra, max_relative = tol);
        assert_relative_eq!(a.dec, b.dec, max_relative = tol);
        assert_relative_eq!(a.time, b.time, max_relative = tol);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_compute_triplets() {
        use camino::Utf8Path;

        use crate::{
            trajectories::trajectory_file::TrajectoryFile, unit_test_global::OUTFIT_HORIZON_TEST, TrajectorySet
        };

        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let mut traj_set =
            TrajectorySet::new_from_80col(&mut env_state, Utf8Path::new("tests/data/2015AB.obs"));

        let traj_number = crate::constants::ObjectNumber::String("K09R05F".into());
        let traj_len = traj_set
            .get(&traj_number)
            .expect("Failed to get trajectory")
            .len();

        let traj_mut = traj_set
            .get_mut(&traj_number)
            .expect("Failed to get trajectory");

        let triplets = generate_triplets(traj_mut, 0.03, 150.0, 20.0, traj_len, 10);

        assert_eq!(
            triplets.len(),
            10,
            "Expected 10 triplets, got {}",
            triplets.len()
        );

        let expected_triplets = GaussObs {
            idx_obs: [[23, 24, 33]].into(),
            ra: [[1.6893715963476699, 1.689861452091063, 1.7527345385664372]].into(),
            dec: [[1.082468037385525, 0.9436790189346231, 0.8273762407899986]].into(),
            time: [[57028.479297592596, 57049.2318575926, 57063.97711759259]].into(),
            observer_helio_position: [
                [-0.2645666171486676, 0.8689351643673471, 0.3766996211112465],
                [-0.5889735526502539, 0.7240117187952059, 0.3138734206791042],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        assert_gauss_obs_approx_eq(&triplets[0], &expected_triplets, 1e-12);

        let expected_triplet = GaussObs {
            idx_obs: [[21, 25, 33]].into(),
            ra: [[1.6894680985108947, 1.6898894500811472, 1.7527345385664372]].into(),
            dec: [[1.0825984522657437, 0.9435805047946215, 0.8273762407899986]].into(),
            time: [[57028.45404759259, 57049.245147592585, 57063.97711759259]].into(),
            observer_helio_position: [
                [-0.26413563361674103, 0.8690466209095019, 0.3767466856686271],
                [-0.5891631852172257, 0.7238872516832191, 0.3138186516545291],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        assert_gauss_obs_approx_eq(&triplets[9], &expected_triplet, 1e-12);
    }

    mod downsampling_observations_tests {
        use nalgebra::Vector3;

        use super::*;

        fn make_obs(n: usize) -> Observations {
            (0..n)
                .map(|i| Observation {
                    observer: 0,
                    ra: 0.0,
                    dec: 0.0,
                    error_ra: 0.0,
                    error_dec: 0.0,
                    time: i as f64,
                    observer_earth_position: Vector3::zeros(),
                    observer_helio_position: Vector3::zeros(),
                })
                .collect()
        }

        #[test]
        fn returns_all_when_max_keep_ge_n() {
            let n = 5;
            let indices = downsample_uniform_with_edges_indices(n, 5);
            assert_eq!(indices, vec![0, 1, 2, 3, 4]);

            let indices = downsample_uniform_with_edges_indices(n, 10);
            assert_eq!(indices, vec![0, 1, 2, 3, 4]);
        }

        #[test]
        fn empty_input_returns_empty() {
            assert!(downsample_uniform_with_edges_indices(0, 0).is_empty());
            assert!(downsample_uniform_with_edges_indices(0, 10).is_empty());
        }

        #[test]
        fn max_keep_less_than_three_returns_first_middle_last() {
            let n = 10;
            let mid = n / 2;
            for max_keep in [0, 1, 2] {
                let indices = downsample_uniform_with_edges_indices(n, max_keep);
                assert_eq!(indices, vec![0, mid, n - 1]);
            }
        }

        #[test]
        fn max_keep_three_exactly_returns_first_middle_last() {
            let n = 10;
            let indices = downsample_uniform_with_edges_indices(n, 3);
            assert_eq!(indices, vec![0, n / 2, n - 1]);

            let n = 3;
            let indices = downsample_uniform_with_edges_indices(n, 3);
            assert_eq!(indices, vec![0, 1, 2]);
        }

        #[test]
        fn downsampling_uniformity_for_general_case() {
            let n = 10;
            let max_keep = 5;
            let indices = downsample_uniform_with_edges_indices(n, max_keep);

            assert_eq!(indices.len(), max_keep);
            assert_eq!(indices.first().unwrap(), &0);
            assert_eq!(indices.last().unwrap(), &(n - 1));

            // Indices doivent être strictement croissants
            assert!(indices.windows(2).all(|w| w[1] > w[0]));
        }

        #[test]
        fn works_with_large_data() {
            let n = 1000;
            let max_keep = 100;
            let indices = downsample_uniform_with_edges_indices(n, max_keep);

            assert_eq!(indices.len(), max_keep);
            assert_eq!(indices.first().unwrap(), &0);
            assert_eq!(indices.last().unwrap(), &(n - 1));
        }

        #[test]
        fn indices_match_observations() {
            let obs = make_obs(10);
            let max_keep = 5;
            let indices = downsample_uniform_with_edges_indices(obs.len(), max_keep);

            // Vérifie que les indices correspondent bien aux temps dans obs
            let times: Vec<_> = indices.iter().map(|&i| obs[i].time).collect();
            assert!(times.windows(2).all(|w| w[1] > w[0]));
        }
    }
}
