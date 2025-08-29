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
//! )?;
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
use crate::observations::Observation;
use crate::outfit_errors::OutfitError;

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
#[derive(Debug)]
struct WeightedTriplet {
    weight: f64,
    i: usize,
    j: usize,
    k: usize,
}

// Implements ordering so that the BinaryHeap acts as a max-heap based on weight.
// The "worst" triplet (largest weight) will be at the top, making pruning efficient.
impl Eq for WeightedTriplet {}
impl PartialEq for WeightedTriplet {
    fn eq(&self, other: &Self) -> bool {
        self.weight == other.weight
    }
}
impl Ord for WeightedTriplet {
    fn cmp(&self, other: &Self) -> Ordering {
        // Compare by weight (descending). Largest weight is considered "greater".
        self.weight
            .partial_cmp(&other.weight)
            .unwrap_or(Ordering::Equal)
    }
}
impl PartialOrd for WeightedTriplet {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Compute a spacing-based weight for a time-ordered observation triplet.
///
/// The score penalizes **unbalanced** or **far-from-optimal** spacings between
/// consecutive observations. A lower score means a more desirable (nearly uniform)
/// spacing close to the target interval `dtw`.
///
/// Formula
/// -----------------
/// For consecutive time gaps `dt₁ = t₂ − t₁`, `dt₂ = t₃ − t₂`,
/// the per-gap penalty is:
/// ```text
/// s(dt, dtw) = { dtw/dt,  if dt ≤ dtw
///              { 1 + dt/dtw, otherwise
/// ```
/// and the total weight is:
/// ```text
/// weight = s(dt₁, dtw) + s(dt₂, dtw)
/// ```
///
/// Arguments
/// -----------------
/// * `time1`, `time2`, `time3` – Observation times (MJD), with `time1 < time2 < time3`.
/// * `dtw` – Desired interval `[days]` between consecutive observations.
///
/// Return
/// ----------
/// * A scalar score (`f64`), **lower is better**.
///
/// Remarks
/// -------------
/// * Complexity: **O(1)**.
/// * Behavior is monotonic in each gap: penalties grow when a gap shrinks well below
///   `dtw` or grows far above it.
/// * Assumes strictly increasing times; if not guaranteed by callers, sort beforehand.
///
/// See also
/// ------------
/// * [`generate_triplets`] – Uses this function to rank candidate triplets.
pub fn triplet_weight(time1: f64, time2: f64, time3: f64, dtw: f64) -> f64 {
    fn s3dtw(dt: f64, dtw: f64) -> f64 {
        if dt <= dtw {
            dtw / dt
        } else {
            1.0 + dt / dtw
        }
    }

    let dt1 = time2 - time1;
    let dt2 = time3 - time2;

    s3dtw(dt1, dtw) + s3dtw(dt2, dtw)
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
fn downsample_uniform_with_edges_indices(n: usize, max_keep: usize) -> Vec<usize> {
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

/// Generate and select optimized triplets of astrometric observations
/// for Gauss initial orbit determination (IOD).
///
/// Overview
/// -----------------
/// This function prepares the candidate triplets of observations that serve
/// as input to the **Gauss method**. It reduces the combinatorial explosion
/// of all possible `(i, j, k)` choices by combining temporal constraints and
/// a spacing-based weight function:
///
/// 1. **Temporal sorting** – Observations are sorted chronologically (in-place).
/// 2. **Downsampling (optional)** – If more than `max_obs_for_triplets` points are present,
///    a uniformly spaced subset is selected (always retaining the first and last).
/// 3. **Triplet enumeration** – All index triplets `(i, j, k)` are generated from
///    the reduced set, subject to the constraint:
///    `dt_min <= (t_k − t_i) <= dt_max`.
/// 4. **Weight scoring** – Each valid triplet is assigned a weight via [`triplet_weight`],
///    favoring evenly spaced observation times close to `optimal_interval_time`.
/// 5. **Heap selection** – A max-heap retains only the `max_triplet` lowest-weight
///    triplets (best candidates).
/// 6. **Conversion** – The selected triplets are converted into [`GaussObs`] structures,
///    with precomputed observer heliocentric positions (using [`Observation::get_observer_helio_position`]).
///
/// Arguments
/// -----------------
/// * `observations` – Mutable vector of astrometric observations.  
///   Sorted **in-place** by observation epoch before processing.
/// * `dt_min` – Minimum allowed timespan `[days]` between first and last observation in a triplet.
/// * `dt_max` – Maximum allowed timespan `[days]` between first and last observation in a triplet.
/// * `optimal_interval_time` – Ideal spacing `[days]` between consecutive observations (used in weighting).
/// * `max_obs_for_triplets` – Maximum number of observations to keep after downsampling.  
///   Larger datasets are uniformly thinned in time (keeping endpoints).
/// * `max_triplet` – Maximum number of best-scoring triplets to return.
///
/// Return
/// ----------
/// * `Result<Vec<GaussObs>, OutfitError>` – At most `max_triplet` triplets, sorted by ascending weight
///   (best geometric distribution first).
///
/// Remarks
/// -------------
/// * Complexity without pruning grows as **O(n³)**; downsampling and weight filtering
///   are critical for large datasets (hundreds or thousands of observations).
/// * The returned [`GaussObs`] can be passed directly to [`GaussObs::prelim_orbit`]
///   to attempt orbit estimation.
/// * This is usually the **first step** in Gauss IOD workflows,
///   often wrapped by [`ObservationsExt::compute_triplets`](crate::observations::observations_ext::ObservationsExt::compute_triplets).
///
/// See also
/// -------------
/// * [`GaussObs`] – Triplet structure used for Gauss orbit determination.
/// * [`triplet_weight`] – Weighting function for scoring observation spacing.
/// * [`ObservationsExt::compute_triplets`](crate::observations::observations_ext::ObservationsExt::compute_triplets) – Higher-level wrapper around this function.
pub fn generate_triplets(
    observations: &mut Observations,
    dt_min: f64,
    dt_max: f64,
    optimal_interval_time: f64,
    max_obs_for_triplets: usize,
    max_triplet: u32,
) -> Result<Vec<GaussObs>, OutfitError> {
    // 1. Sort observations by time (in-place)
    observations.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

    // 2. Downsample the observations (returns a list of indices)
    let selected_indices =
        downsample_uniform_with_edges_indices(observations.len(), max_obs_for_triplets);

    // 3. Build a reduced view using only selected indices
    let reduced: Vec<&Observation> = selected_indices.iter().map(|&i| &observations[i]).collect();
    let n = reduced.len();

    // Use a max-heap to keep only the best triplets (lowest weights)
    let mut heap = BinaryHeap::with_capacity((max_triplet + 1) as usize);

    // Capture reduced by reference to avoid moving it into the closure
    let reduced_ref = &reduced;

    // 4. Generate triplets using iterators
    (0..n)
        .flat_map(|i| {
            (i + 1..n).flat_map(move |j| {
                let reduced_ref = reduced_ref;
                (j + 1..n).filter_map(move |k| {
                    let reduced_ref = reduced_ref;

                    // Check time span constraints
                    let dt13 = reduced_ref[k].time - reduced_ref[i].time;
                    if dt13 < dt_min || dt13 > dt_max {
                        return None;
                    }

                    // Compute triplet weight
                    let weight = triplet_weight(
                        reduced_ref[i].time,
                        reduced_ref[j].time,
                        reduced_ref[k].time,
                        optimal_interval_time,
                    );

                    Some(WeightedTriplet { weight, i, j, k })
                })
            })
        })
        // 5. Maintain a bounded heap of best triplets
        .for_each(|wt| {
            if heap.len() < max_triplet as usize {
                heap.push(wt);
            } else if let Some(top) = heap.peek() {
                // Replace the worst triplet if a better one is found
                if wt.weight < top.weight {
                    heap.pop();
                    heap.push(wt);
                }
            }
        });

    // Extract and sort the best triplets from the heap
    let mut best: Vec<_> = heap.into_sorted_vec();
    best.sort_by(|a, b| a.weight.partial_cmp(&b.weight).unwrap());

    // 6. Convert selected triplets into GaussObs
    best.into_iter()
        .map(|wt| {
            let WeightedTriplet { i, j, k, .. } = wt;
            let idx1 = selected_indices[i];
            let idx2 = selected_indices[j];
            let idx3 = selected_indices[k];

            GaussObs::with_observer_position(
                Vector3::new(idx1, idx2, idx3),
                Vector3::new(reduced[i].ra, reduced[j].ra, reduced[k].ra),
                Vector3::new(reduced[i].dec, reduced[j].dec, reduced[k].dec),
                Vector3::new(reduced[i].time, reduced[j].time, reduced[k].time),
                Matrix3::from_columns(&[
                    reduced[i].get_observer_helio_position(),
                    reduced[j].get_observer_helio_position(),
                    reduced[k].get_observer_helio_position(),
                ]),
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
            constants::TrajectorySet, error_models::ErrorModel,
            observations::trajectory_ext::TrajectoryExt, outfit::Outfit,
        };

        let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
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

        let triplets = generate_triplets(traj_mut, 0.03, 150.0, 20.0, traj_len, 10).unwrap();

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
