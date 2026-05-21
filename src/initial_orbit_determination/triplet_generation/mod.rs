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
//! ```rust,ignore
//! use outfit::initial_orbit_determination::triplet_generation::generate_triplets;
//!
//! // Your observations loaded elsewhere (sorted by ascending epoch):
//! let obs: Vec<&Observation> = unimplemented!();
//! let cache: &OutfitCache = unimplemented!();
//!
//! // Build best triplets for Gauss IOD
//! let triplets = generate_triplets(
//!     &obs,
//!     cache,
//!     0.03,   // dt_min [days]
//!     150.0,  // dt_max [days]
//!     20.0,   // optimal_interval_time [days]
//!     100,    // max_obs_for_triplets
//!     10      // max_triplet
//! );
//!
//! // Now pass `triplets` to your Gauss solver…
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
//! - `ObservationsExt::compute_triplets` – Higher-level wrapper (internal API).

pub mod index_generator;

use nalgebra::{Matrix3, Vector3};
use photom::observation_dataset::observation::Observation;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::cache::OutfitCache;
use crate::initial_orbit_determination::gauss::GaussObs;
use crate::initial_orbit_determination::triplet_generation::index_generator::TripletIndexGenerator;
use crate::IODParams;

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

/// Generate and select **best-K** triplets of astrometric observations
/// for Gauss Initial Orbit Determination (IOD), using a **lazy index stream**
/// and a bounded **max-heap** on a spacing weight.
///
/// Overview
/// -----------------
/// This routine constructs good candidate triplets `(first, middle, last)` as inputs
/// to the **Gauss method** while avoiding the `O(n³)` blow-up:
///
/// 1. **Downsample (in `TripletIndexGenerator`)** – Observations are uniformly thinned
///    (keeping endpoints) to at most `max_obs_for_triplets`. The input slice is assumed
///    to be **already sorted by ascending epoch**.
/// 2. **Time-feasible enumeration (lazy)** – Indices `(i, j, k)` are streamed by
///    `TripletIndexGenerator`, constrained by:
///    `dt_min ≤ t[k] − t[i] ≤ dt_max` with `i < j < k`.
/// 3. **Weight scoring** – Each feasible triplet receives a weight via [`triplet_weight`],
///    favoring near-uniform spacing around `optimal_interval_time`.
/// 4. **Best-K selection** – A bounded **max-heap** retains only the `max_triplet`
///    lowest-weight candidates (the heap's `peek()` is the current **worst**).
/// 5. **Materialization** – Only for the selected indices, we re-borrow `observations`
///    immutably and build [`GaussObs`] with precomputed observer heliocentric columns.
///
/// # Design notes
///
/// - The input slice must be **sorted in ascending time order** before this call.
/// - Enumeration and scoring happen on reduced epoch indices owned by the generator,
///   so there are **no overlapping borrows** of `observations`.
/// - The final `Vec<GaussObs>` is **sorted by increasing weight** (best first).
///
/// # Arguments
///
/// - `observations`         – Time-sorted observation slice (ascending epoch).
/// - `cache`                – Precomputed heliocentric observer positions.
/// - `params`               – IOD parameters controlling time bounds, downsampling cap,
///   optimal spacing, and the best-K limit.
///
/// # Return
///
/// - A `Vec<GaussObs>` of length `≤ max_triplet`, sorted by **ascending** heuristic weight.
///
/// # Complexity
///
/// * Enumeration: typically ~$O(n^2)$ thanks to the per-anchor time window in
///   [`TripletIndexGenerator`].
/// * Selection: $O(n \log K)$ due to the bounded heap.
/// * Space: $O(1)$ per yielded triplet during enumeration; only the final `K` are materialized.
///
/// # See also
///
/// * [`TripletIndexGenerator`] – Streams time-feasible reduced indices lazily.
/// * [`triplet_weight`] – Heuristic favoring evenly spaced triplets around a target gap.
/// * [`GaussObs::realizations_iter`] – Lazy Monte-Carlo perturbations per triplet.
pub fn generate_triplets(
    observations: &[&Observation],
    cache: &OutfitCache,
    params: &IODParams,
) -> Vec<GaussObs> {
    if params.max_triplets == 0 || observations.len() < 3 {
        return Vec::new();
    }

    // --- Phase 1: build the reduced-index stream and score all feasible triplets.
    let mut index_gen = TripletIndexGenerator::from_observations(
        observations,
        params.dt_min,
        params.dt_max_triplet,
        params.max_obs_for_triplets,
        usize::MAX,
    );

    let k_cap = params.max_triplets as usize;
    let inv_dtw = params.optimal_interval_time.recip();
    let best_k = collect_best_k_triplets(&mut index_gen, k_cap, inv_dtw);

    // --- Phase 2: materialize GaussObs for the selected indices.
    // Indices in WeightedTriplet refer directly into the downsampled view of
    // `observations` — no remapping is needed.
    best_k
        .into_iter()
        .map(|wt| build_gauss_obs(cache, observations, wt))
        .collect()
}

/// Consume the triplet stream and retain the `max_triplets` best candidates by ascending weight.
///
/// Uses a bounded max-heap: when the heap is full, a new candidate replaces the
/// current worst only if its weight is strictly smaller.
///
/// Returns candidates sorted by ascending weight.
fn collect_best_k_triplets(
    gen: &mut TripletIndexGenerator,
    max_triplets: usize,
    inv_optimal_interval: f64,
) -> Vec<WeightedTriplet> {
    let mut heap: BinaryHeap<WeightedTriplet> =
        BinaryHeap::with_capacity(max_triplets.saturating_add(1));

    while let Some((first, middle, last)) = gen.next() {
        let times = gen.reduced_times();
        let weight = triplet_weight_with_inv(
            times[first],
            times[middle],
            times[last],
            inv_optimal_interval,
        );

        if !weight.is_finite() {
            continue;
        }

        if heap.len() < max_triplets {
            heap.push(WeightedTriplet {
                weight,
                first_idx: first,
                middle_idx: middle,
                last_idx: last,
            });
        } else if heap.peek().is_some_and(|worst| weight < worst.weight) {
            heap.pop();
            heap.push(WeightedTriplet {
                weight,
                first_idx: first,
                middle_idx: middle,
                last_idx: last,
            });
        }
    }

    // Max-heap → ascending weight order.
    let mut result = heap.into_vec();
    result.sort_unstable_by(|a, b| a.weight.partial_cmp(&b.weight).unwrap_or(Ordering::Equal));
    result
}

/// Materialize a single [`GaussObs`] from a [`WeightedTriplet`].
///
/// Indices in `wt` map directly into `observations` — no remapping is needed
/// since [`TripletIndexGenerator`] works on the input slice as-is.
fn build_gauss_obs(
    cache: &OutfitCache,
    observations: &[&Observation],
    wt: WeightedTriplet,
) -> GaussObs {
    let o1 = &observations[wt.first_idx];
    let o2 = &observations[wt.middle_idx];
    let o3 = &observations[wt.last_idx];

    let observer_matrix = Matrix3::from_columns(&[
        *cache.get_helio_position(o1.index()),
        *cache.get_helio_position(o2.index()),
        *cache.get_helio_position(o3.index()),
    ]);

    let (o1_ra, o1_dec) = (o1.equ_coord().ra, o1.equ_coord().dec);
    let (o2_ra, o2_dec) = (o2.equ_coord().ra, o2.equ_coord().dec);
    let (o3_ra, o3_dec) = (o3.equ_coord().ra, o3.equ_coord().dec);

    GaussObs::with_observer_position(
        Vector3::new(wt.first_idx, wt.middle_idx, wt.last_idx),
        Vector3::new(o1_ra, o2_ra, o3_ra),
        Vector3::new(o1_dec, o2_dec, o3_dec),
        Vector3::new(o1.mjd_tt(), o2.mjd_tt(), o3.mjd_tt()),
        observer_matrix.map(|x| x.into_inner()),
    )
}

#[cfg(test)]
mod triplets_iod_tests {
    use super::*;

    #[test]
    fn test_compute_triplets() {
        use crate::cache::OutfitCache;
        use crate::test_fixture::{DATASET_2015AB, JPL_EPHEM_HORIZON, UT1_PROVIDER};
        use crate::IODParams;
        use photom::observer::error_model::{ModelCorrection, ObsErrorModel};

        // The error model must be set before building the cache.
        let dataset = DATASET_2015AB
            .clone()
            .with_error_model(ObsErrorModel::FCCT14)
            .apply_batch_rms_correction(30.0);

        // Build the cache from the real 2015AB dataset.
        let cache = OutfitCache::build(&dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

        // Pick a trajectory with enough observations.
        let traj = dataset
            .materialize_trajectory("K09R05F")
            .unwrap()
            .collect_into_vec();

        assert!(
            traj.len() >= 3,
            "trajectory must have at least 3 observations"
        );

        let params = IODParams {
            dt_min: 0.03,
            dt_max_triplet: 150.0,
            optimal_interval_time: 20.0,
            max_obs_for_triplets: traj.len(),
            max_triplets: 10,
            ..Default::default()
        };

        let triplets = generate_triplets(&traj, &cache, &params);

        // We should get at least one triplet back.
        assert!(!triplets.is_empty(), "expected at least one triplet");

        // No more than max_triplets.
        assert!(
            triplets.len() <= params.max_triplets as usize,
            "got {} triplets, expected ≤ {}",
            triplets.len(),
            params.max_triplets
        );

        // Triplets must be sorted by ascending weight (best first).
        // We verify this by re-computing weights and checking order.
        for window in triplets.windows(2) {
            let t1 = &window[0].time;
            let t2 = &window[1].time;
            let w1 = triplet_weight(t1[0], t1[1], t1[2], params.optimal_interval_time);
            let w2 = triplet_weight(t2[0], t2[1], t2[2], params.optimal_interval_time);
            assert!(
                w1 <= w2 + 1e-12,
                "triplets not sorted by ascending weight: w1={w1} > w2={w2}"
            );
        }

        // Each triplet's indices must be strictly increasing.
        for t in &triplets {
            assert!(
                t.idx_obs[0] < t.idx_obs[1] && t.idx_obs[1] < t.idx_obs[2],
                "triplet indices not strictly increasing: {:?}",
                t.idx_obs
            );
        }
    }

    mod downsampling_observations_tests {

        use crate::initial_orbit_determination::triplet_generation::index_generator::downsample_uniform_with_edges;

        #[test]
        fn returns_all_when_max_keep_ge_n() {
            let n = 5;
            let indices = downsample_uniform_with_edges(n, 5);
            assert_eq!(indices, vec![0, 1, 2, 3, 4]);

            let indices = downsample_uniform_with_edges(n, 10);
            assert_eq!(indices, vec![0, 1, 2, 3, 4]);
        }

        #[test]
        fn empty_input_returns_empty() {
            assert!(downsample_uniform_with_edges(0, 0).is_empty());
            assert!(downsample_uniform_with_edges(0, 10).is_empty());
        }

        #[test]
        fn max_keep_less_than_three_returns_first_middle_last() {
            let n = 10;
            let mid = n / 2;
            for max_keep in [1, 2, 3] {
                let indices = downsample_uniform_with_edges(n, max_keep);
                assert_eq!(indices, vec![0, mid, n - 1]);
            }
        }

        #[test]
        fn max_keep_three_exactly_returns_first_middle_last() {
            let n = 10;
            let indices = downsample_uniform_with_edges(n, 3);
            assert_eq!(indices, vec![0, n / 2, n - 1]);

            let n = 3;
            let indices = downsample_uniform_with_edges(n, 3);
            assert_eq!(indices, vec![0, 1, 2]);
        }

        #[test]
        fn downsampling_uniformity_for_general_case() {
            let n = 10;
            let max_keep = 5;
            let indices = downsample_uniform_with_edges(n, max_keep);

            assert_eq!(indices.len(), max_keep);
            assert_eq!(indices.first().unwrap(), &0);
            assert_eq!(indices.last().unwrap(), &(n - 1));

            // Indices must be strictly increasing
            assert!(indices.windows(2).all(|w| w[1] > w[0]));
        }

        #[test]
        fn works_with_large_data() {
            let n = 1000;
            let max_keep = 100;
            let indices = downsample_uniform_with_edges(n, max_keep);

            assert_eq!(indices.len(), max_keep);
            assert_eq!(indices.first().unwrap(), &0);
            assert_eq!(indices.last().unwrap(), &(n - 1));
        }
    }
}
