use approx::assert_abs_diff_eq;
use camino::Utf8Path;
use hifitime::ut1::Ut1Provider;
use nalgebra::Vector3;
use outfit::{
    cache::OutfitCache, jpl_ephem::download_jpl_file::EphemFileSource, IODParams, JPLEphem,
};
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

        let jpl_file: EphemFileSource = "horizon:DE440"
            .try_into()
            .expect("Failed to parse JPL ephemeris source");
        let jpl_ephem =
            JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Horizon");

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

fn assert_helio_position(cache: &OutfitCache, obs_dataset: &ObsDataset, trajectory_id: &str) {
    let traj = obs_dataset
        .materialize_trajectory(trajectory_id)
        .unwrap()
        .collect_into_vec();

    let checks: &[(usize, [f64; 3])] = &[
        (
            0,
            [0.996798968524259, -0.12232935537370689, -0.0530044254994447],
        ),
        (
            traj.len() / 2,
            [-0.2588304494454786, 0.8703635675926336, 0.3773300400916685],
        ),
        (
            traj.len() - 1,
            [-0.8383497659757538, 0.479843435538848, 0.20801659206288547],
        ),
    ];

    for (idx, expected) in checks {
        let obs = &traj[*idx];
        let helio_pos = cache
            .get_centric(obs.index())
            .helio_position
            .map(|x| x.into_inner());

        assert_abs_diff_eq!(
            helio_pos,
            Vector3::from(*expected),
            epsilon = POSITION_EPSILON
        );
    }
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

    for paths in dataset_combinations {
        let (cache, obs_dataset) = fixture.build_cache(paths);
        assert_helio_position(&cache, &obs_dataset, "K09R05F");
    }
}
