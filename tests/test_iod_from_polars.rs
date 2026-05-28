mod common;

use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::orbit_type::{keplerian_element::KeplerianElements, OrbitalElements};
use outfit::{FitIOD, IODParams, JPLEphem};
use photom::io::polars::ContiguousChoice;
use photom::{
    io::polars::FromPolarsArgs, observation_dataset::ObsDataset,
    observer::error_model::ObsErrorModel,
};
use polars::{
    lazy::{
        dsl::{col, lit},
        frame::LazyFrame,
    },
    prelude::NamedFrom,
    series::Series,
};
use rand::{rngs::StdRng, SeedableRng};

use crate::common::approx_equal;

#[test]
fn test_iod_from_polars() {
    let test_max_relative = 1e-11;
    let test_epsilon = 1e-11;

    let path_data = "tests/data/test_data_traj_str.parquet";

    let ids = Series::new("".into(), ["95777", "14226", "29757"]);
    let lf = LazyFrame::scan_parquet(path_data.into(), Default::default())
        .expect("scan_parquet must succeed")
        .filter(col("traj_id").is_in(lit(ids).implode(), true));

    let polars_args = FromPolarsArgs {
        error_model: Some(ObsErrorModel::FCCT14),
        do_rechunk: Some(false),
        contiguous_choice: Some(ContiguousChoice::ContiguousTraj),
    };

    let obs_dataset = ObsDataset::from_lazy(lf, polars_args).unwrap();

    let max_traj_size = obs_dataset
        .iter_traj_id()
        .unwrap()
        .map(|traj_id| obs_dataset.len_trajectory(traj_id).unwrap())
        .max()
        .unwrap();

    let default = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(max_traj_size)
        .max_triplets(30)
        .build()
        .unwrap();

    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("Download of the JPL short time scale UT1 data failed");

    let jpl_ephem: JPLEphem = "horizon:DE440"
        .try_into()
        .expect("Failed to load JPL ephemeris");

    let mut full_orbit = obs_dataset
        .fit_full_iod(
            &jpl_ephem,
            &ut1_provider,
            &default,
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    // --- traj 14226 ---
    let best_orbit = full_orbit.remove(&"14226".into()).unwrap().unwrap();
    let orbit = best_orbit.orbital_elements();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60894.372896385554,
        semi_major_axis: 0.5415009930884174,
        eccentricity: 0.9027228307040831,
        inclination: 0.31200939353818746,
        ascending_node_longitude: 5.550735343593096,
        periapsis_argument: 3.1638244350882596,
        mean_anomaly: 2.7888128618151495,
    });

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.02704195897369085,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // --- traj 29757 ---
    let best_orbit = full_orbit.remove(&"29757".into()).unwrap().unwrap();
    let orbit = best_orbit.orbital_elements();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60835.25573266984,
        semi_major_axis: 0.635955220245824,
        eccentricity: 0.5904849550180079,
        inclination: 0.2263529126300279,
        ascending_node_longitude: 4.366539949885583,
        periapsis_argument: 3.3107966723035602,
        mean_anomaly: 3.0157533331616966,
    });

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.025397381294328548,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // --- traj 95777 ---
    let best_orbit = full_orbit.remove(&"95777".into()).unwrap().unwrap();
    let orbit = best_orbit.orbital_elements();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60894.252965553926,
        semi_major_axis: 1.24701989952089,
        eccentricity: 0.2082069422196415,
        inclination: 0.08116316040114972,
        ascending_node_longitude: 2.49554922649176,
        periapsis_argument: 2.5470318525197477,
        mean_anomaly: 0.2983936748249412,
    });

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.010284390096535293,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}
