//! Triplet selection for initial orbit determination (Gauss method).
//!
//! This module provides tools to:
//! - Downsample a large set of observations to a representative subset,
//! - Generate triplets of observations suitable for the Gauss IOD method,
//! - Score triplets using a weight function that favors well-spaced times,
//! - Select the best triplets using a heap-based pruning algorithm.

use std::cmp::Ordering;
use std::collections::BinaryHeap;

use crate::constants::Observations;
use crate::initial_orbit_determination::gauss::GaussObs;
use crate::observations::Observation;
use crate::observers::Observer;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;

/// Internal structure used to store a weighted observation triplet during
/// the selection process.
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

/// Compute the "weight" of a triplet based on time spacing.
///
/// The weight is low when the intervals between observations are close
/// to the `dtw` (optimal interval), and larger when the intervals are
/// unbalanced or far from `dtw`.
///
/// # Arguments
/// * `time1`, `time2`, `time3` - times of the three observations (MJD).
/// * `dtw` - Optimal desired time spacing between consecutive observations.
///
/// # Returns
/// A floating-point score (lower is better).
fn triplet_weight(time1: f64, time2: f64, time3: f64, dtw: f64) -> f64 {
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

/// Downsample a set of observations by returning a list of indices.
///
/// If the dataset is large, only `max_keep` points are kept, distributed
/// uniformly in time. The first and last observations are always included.
///
/// * If `max_keep < 3`, the first, middle, and last indices are returned.
/// * If `max_keep >= n`, all indices are returned.
///
/// This ensures a reasonable coverage of the time span without processing
/// every single observation.
///
/// # Arguments
/// * `n` - Total number of observations.
/// * `max_keep` - Maximum number of indices to return.
///
/// # Returns
/// A vector of indices of the selected observations.
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

/// Generate and select high-quality triplets of astrometric observations
/// for initial orbit determination using the Gauss method.
///
/// This function implements the preprocessing stage of the IOD pipeline:
///
/// 1. **Temporal sorting** – Observations are sorted chronologically.
/// 2. **Downsampling (optional)** – If the dataset is very large, only
///    `max_obs_for_triplets` observations are retained using a uniform
///    time distribution (always keeping the first and last).
/// 3. **Triplet enumeration** – All index triplets `(i, j, k)` are generated
///    from the reduced set, restricted to those where the time span between
///    the first and last observation satisfies:
///    `dt_min <= (t_k - t_i) <= dt_max`.
/// 4. **Weight scoring** – Each valid triplet is assigned a weight using
///    [`triplet_weight`], which favors evenly spaced observation times.
/// 5. **Selection** – A max-heap is used to retain only the `max_triplet`
///    triplets with the lowest weights (i.e., most promising for Gauss IOD).
/// 6. **Output** – The selected triplets are converted into [`GaussObs`]
///    structures, including their observer positions computed from [`Outfit`].
///
/// This procedure drastically reduces the combinatorial explosion that
/// arises when working with dense observation sets and ensures that the
/// best geometrically spaced triplets are passed to the Gauss solver.
///
/// # Arguments
///
/// * `observations` - Mutable vector of astrometric observations. The vector
///   is sorted in-place by observation time before processing.
/// * `state` - Global [`Outfit`] context used to compute the observer's
///   position at the time of each observation.
/// * `dt_min` - Minimum time span (in days) allowed between the first and
///   last observation of a triplet. Default: `0.03`.
/// * `dt_max` - Maximum time span (in days) allowed between the first and
///   last observation of a triplet. Default: `150.0`.
/// * `optimal_interval_time` - Ideal time spacing (in days) between
///   consecutive observations in a triplet, used when computing the
///   triplet weight. Default: `20.0`.
/// * `max_obs_for_triplets` - Maximum number of observations to keep after
///   downsampling. If the input dataset is larger, a subset is selected
///   uniformly in time. Default: `100`.
/// * `max_triplet` - Maximum number of best triplets to return. Default: `10`.
///
/// # Returns
///
/// A `Vec<GaussObs>` containing at most `max_triplet` triplets, sorted in
/// ascending order of weight (best first).
///
/// # Notes
///
/// * This function is typically the first step of the Gauss initial orbit
///   determination process.
/// * Downsampling is critical when the number of observations is very large
///   (hundreds or thousands), as the total number of combinations grows as `O(n³)`.
///
/// # See also
/// * [`GaussObs`] – Representation of a triplet ready for Gauss IOD.
/// * [`triplet_weight`] – The scoring function used to rank triplets.
pub(crate) fn generate_triplets(
    observations: &mut Observations,
    state: &Outfit,
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
    let precomputed_observers: Vec<&Observer> =
        reduced.iter().map(|obs| obs.get_observer(state)).collect();
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
                state,
                nalgebra::Vector3::new(idx1, idx2, idx3),
                nalgebra::Vector3::new(reduced[i].ra, reduced[j].ra, reduced[k].ra),
                nalgebra::Vector3::new(reduced[i].dec, reduced[j].dec, reduced[k].dec),
                nalgebra::Vector3::new(reduced[i].time, reduced[j].time, reduced[k].time),
                [
                    precomputed_observers[i],
                    precomputed_observers[j],
                    precomputed_observers[k],
                ],
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

        let triplets =
            generate_triplets(traj_mut, &env_state, 0.03, 150.0, 20.0, traj_len, 10).unwrap();

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
            observer_position: [
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
            observer_position: [
                [-0.26413563361674103, 0.8690466209095019, 0.3767466856686271],
                [-0.5891631852172257, 0.7238872516832191, 0.3138186516545291],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        assert_gauss_obs_approx_eq(&triplets[9], &expected_triplet, 1e-12);
    }

    mod downsampling_observations_tests {
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
