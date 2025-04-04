use std::collections::HashMap;

use itertools::Itertools;
use nalgebra::Vector3;

use crate::{
    constants::{TrajectorySet, Triplets},
    outfit::Outfit,
};

use super::gauss::GaussObs;

fn triplets_selection<'a>(
    trajectory_set: &'a mut TrajectorySet,
    env_state: &'a mut Outfit,
    dt_min: Option<f64>,
    dt_max: Option<f64>,
    max_triplet: Option<u32>,
) -> Triplets<'a> {
    let mut result = HashMap::default();

    for (object_number, observations) in trajectory_set {
        observations.sort_by(|a, b| {
            a.time.partial_cmp(&b.time).expect(
                format!(
                    "Error sorting observations for object number {:?}: {:?} vs {:?}",
                    object_number, a.time, b.time
                )
                .as_str(),
            )
        });
        let mut triplets: Vec<GaussObs> = Vec::new();
        let mut nb_triplets = 0;

        for (o1, o2, o3) in observations.iter().tuple_combinations() {
            let dt13 = o3.time - o1.time;

            if dt13 >= dt_min.unwrap_or(0.03) && dt13 <= dt_max.unwrap_or(150.) {
                // Construction de la structure GaussObs
                let gauss = GaussObs::new(
                    Vector3::new(o1.ra, o2.ra, o3.ra),
                    Vector3::new(o1.dec, o2.dec, o3.dec),
                    Vector3::new(o1.time, o2.time, o3.time),
                    [
                        o1.get_observer(env_state),
                        o2.get_observer(env_state),
                        o3.get_observer(env_state),
                    ],
                );
                triplets.push(gauss);
                nb_triplets += 1;
                if nb_triplets >= max_triplet.unwrap_or(10) {
                    break;
                }
            }
        }
        if !triplets.is_empty() {
            result.insert(object_number.clone(), triplets);
        }
    }
    result
}

#[cfg(test)]
mod trajectory_ext_test {
    use camino::Utf8Path;

    use crate::{
        constants::ObjectNumber, observations::trajectory_ext::TrajectoryExt, outfit::Outfit,
    };

    use super::*;

    #[test]
    fn test_triplets_selection() {
        let mut env_state = Outfit::new();
        let mut traj_set =
            TrajectorySet::new_from_80col(&mut env_state, &Utf8Path::new("tests/data/8467.obs"));

        let triplets = triplets_selection(
            &mut traj_set,
            &mut env_state,
            Some(0.03),
            Some(150.),
            Some(10),
        );

        assert_eq!(
            triplets
                .get(&ObjectNumber::String("8467".into()))
                .unwrap()
                .len(),
            10
        );
        println!("{:?}", triplets);
    }
}
