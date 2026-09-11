mod common;

use approx::assert_abs_diff_eq;
use camino::Utf8Path;
use hifitime::ut1::Ut1Provider;
use nalgebra::Vector3;
use outfit::{cache::OutfitCache, IODParams, JPLEphem};
use photom::{
    observation_dataset::ObsDataset,
    observer::error_model::{ModelCorrection, ObsErrorModel},
};

const POSITION_EPSILON: f64 = 1e-12;

struct CacheFixture {
    ut1_provider: Ut1Provider,
    jpl_ephem: JPLEphem,
    default_params: IODParams,
}

impl CacheFixture {
    fn new() -> Self {
        let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
            .expect("Download of the JPL short time scale UT1 data failed");

        let jpl_ephem = common::load_ephem();

        Self {
            ut1_provider,
            jpl_ephem,
            default_params: IODParams::default(),
        }
    }

    fn build_cache<P: AsRef<Utf8Path>>(&self, paths: &[P]) -> (OutfitCache, ObsDataset) {
        let (obs_dataset, _) = ObsDataset::from_mpc_80_col_files(paths);

        let obs_dataset = obs_dataset
            .with_error_model(ObsErrorModel::FCCT14)
            .apply_batch_rms_correction(self.default_params.gap_max);

        let cache = OutfitCache::build(&obs_dataset, &self.jpl_ephem, &self.ut1_provider, false)
            .expect("Failed to build outfit cache");

        (cache, obs_dataset)
    }
}

/// Sample the cached heliocentric position of `trajectory_id` at its first,
/// middle and last observation.
///
/// # Returns
///
/// The three sampled positions, in AU, in observation order.
fn sample_helio_positions(
    cache: &OutfitCache,
    obs_dataset: &ObsDataset,
    trajectory_id: &str,
) -> [Vector3<f64>; 3] {
    let traj = obs_dataset
        .materialize_trajectory(trajectory_id)
        .unwrap()
        .collect_into_vec();

    let indices = [0, traj.len() / 2, traj.len() - 1];
    indices.map(|idx| {
        cache
            .get_centric(traj[idx].index())
            .helio_position
            .map(|x| x.into_inner())
    })
}

#[test]
fn test_cache_consistency() {
    let fixture = CacheFixture::new();

    let dataset_combinations: &[&[&str]] = &[
        &["tests/data/2015AB.obs"],
        &["tests/data/8467.obs", "tests/data/2015AB.obs"],
        &[
            "tests/data/2015AB.obs",
            "tests/data/8467.obs",
            "tests/data/33803.obs",
        ],
    ];

    // The cached heliocentric position of a given object must not depend on which
    // sibling datasets are loaded alongside it: every combination must reproduce
    // the positions obtained from the first combination, whatever the backend.
    let mut reference: Option<[Vector3<f64>; 3]> = None;
    for paths in dataset_combinations {
        let (cache, obs_dataset) = fixture.build_cache(paths);
        let sampled = sample_helio_positions(&cache, &obs_dataset, "K09R05F");
        match &reference {
            None => reference = Some(sampled),
            Some(reference) => {
                for (got, want) in sampled.iter().zip(reference.iter()) {
                    assert_abs_diff_eq!(got, want, epsilon = POSITION_EPSILON);
                }
            }
        }
    }

    // Sanity check on the reference itself: the object is between ~0.9 and ~1.6 AU
    // from the Sun over the arc.
    for position in reference.expect("at least one dataset combination") {
        assert!((0.8..2.0).contains(&position.norm()));
    }
}
