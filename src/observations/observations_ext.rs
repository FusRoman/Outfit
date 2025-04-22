use itertools::Itertools;
use nalgebra::Vector3;

use crate::{constants::Observations, initial_orbit_determination::gauss::GaussObs};

/// Calculate the weight of a triplet of observations based on the time difference.
///
/// Arguments
/// ---------
/// * `time`: A reference to a `Vector3<f64>` containing the time of the three observations.
/// * `dtw`: A `f64` representing an optimal interval time requested between the observations of the triplet.
///
/// Return
/// ------
/// * A `f64` representing the weight of the triplet.
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

pub(crate) trait ObservationsExt {
    fn compute_triplets(
        &mut self,
        dt_min: Option<f64>,
        dt_max: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplet: Option<u32>,
    ) -> Vec<GaussObs>;
}

impl ObservationsExt for Observations {

    /// Compute triplets of observations.
    /// 
    /// This function computes triplets of observations from the given set of observations.
    /// It filters the observations based on the time difference between the first and last observation in the triplet,
    /// and sorts the triplets based on their weights.
    /// 
    /// Arguments
    /// ---------
    /// * `dt_min`: An optional minimum time difference between the first and last observation in the triplet. (default: 0.03)
    /// * `dt_max`: An optional maximum time difference between the first and last observation in the triplet. (default: 150.0)
    /// * `optimal_interval_time`: An optional optimal interval time requested between the observations of the triplet. (default: 20.0)
    /// * `max_triplet`: An optional maximum number of triplets to return. (default: 10)
    /// 
    /// Return
    /// ------
    /// * A vector of `GaussObs` representing the computed triplets of observations.
    fn compute_triplets(
        &mut self,
        dt_min: Option<f64>,
        dt_max: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplet: Option<u32>,
    ) -> Vec<GaussObs> {
        self.sort_by(|a, b| {
            a.time.partial_cmp(&b.time).expect(
                format!(
                    "Error sorting observations for object number: {:?} vs {:?}",
                    a.time, b.time
                )
                .as_str(),
            )
        });

        self.iter()
            .enumerate()
            .tuple_combinations::<(_, _, _)>()
            .filter_map(|(obs1, obs2, obs3)| {
                let dt13 = obs3.1.time - obs1.1.time;
                if dt13 < dt_min.unwrap_or(0.03) || dt13 > dt_max.unwrap_or(150.) {
                    return None;
                }

                Some((
                    (obs1.0, obs2.0, obs3.0),
                    triplet_weight(
                        obs1.1.time,
                        obs2.1.time,
                        obs3.1.time,
                        optimal_interval_time.unwrap_or(20.0),
                    ),
                ))
            })
            .sorted_by(|(_, w1), (_, w2)| {
                w1.partial_cmp(w2).expect("Error sorting triplet weights")
            })
            .take(max_triplet.unwrap_or(10) as usize)
            .map(|((idx1, idx2, idx3), _)| {
                let obs1 = &self[idx1];
                let obs2 = &self[idx2];
                let obs3 = &self[idx3];

                GaussObs::new(
                    Vector3::new(obs1.ra, obs2.ra, obs3.ra),
                    Vector3::new(obs1.dec, obs2.dec, obs3.dec),
                    Vector3::new(obs1.time, obs2.time, obs3.time),
                )
            })
            .collect::<Vec<GaussObs>>()
    }
}

#[cfg(test)]
mod test_obs_ext {
    use camino::Utf8Path;

    use crate::{
        constants::TrajectorySet, observations::trajectory_ext::TrajectoryExt, outfit::Outfit,
    };

    use super::*;

    #[test]
    fn test_compute_triplets() {
        let mut env_state = Outfit::new("horizon:DE440");
        let mut traj_set =
            TrajectorySet::new_from_80col(&mut env_state, &Utf8Path::new("tests/data/2015AB.obs"));

        let triplets = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory")
            .compute_triplets(Some(0.03), Some(150.), None, Some(10));

        assert_eq!(
            triplets.len(),
            10,
            "Expected 10 triplets, got {}",
            triplets.len()
        );

        assert_eq!(
            triplets[0],
            GaussObs::new(
                [[96.7938625, 96.82192916666668, 100.42429166666668]].into(),
                [[62.020849999999996, 54.068825000000004, 47.405166666666666]].into(),
                [[57028.479297592596, 57049.2318575926, 57063.97711759259]].into(),
            )
        );

        assert_eq!(
            triplets[9],
            GaussObs::new(
                [[96.79939166666666, 96.82353333333334, 100.42429166666668]].into(),
                [[62.02832222222222, 54.063180555555554, 47.405166666666666]].into(),
                [[57028.45404759259, 57049.245147592585, 57063.97711759259]].into(),
            )
        );
    }
}
