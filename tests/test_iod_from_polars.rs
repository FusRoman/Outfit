//! Exact-oracle non-regression test for the Gauss IOD pipeline on a
//! parquet-sourced dataset. The reference orbits reproduce the in-house
//! Horizon ephemeris reader's arithmetic bit-for-bit (`1e-11`), so the whole
//! file is gated on `ephem-builtin` and forces that reader explicitly rather
//! than following `JPLEphem::new`'s compile-time backend precedence — ANISE
//! only agrees with it to `~1e-8`.
#![cfg(feature = "ephem-builtin")]

mod common;

use approx::assert_relative_eq;
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
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

    let ut1_provider = common::load_ut1();

    let source: EphemFileSource = "horizon:DE440"
        .try_into()
        .expect("failed to parse the Horizon ephemeris source");
    let jpl_ephem: JPLEphem = JPLEphem::from_builtin(source)
        .expect("failed to load the DE440 ephemeris through the in-house Horizon reader");

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

    let expected_orbit = OrbitalElements::Keplerian {
        elements: KeplerianElements {
            reference_epoch: 60894.372896385554,
            semi_major_axis: 0.541500993374083,
            eccentricity: 0.9027228296226688,
            inclination: 0.3120093908106981,
            ascending_node_longitude: 5.550735343624791,
            periapsis_argument: 3.163824435180773,
            mean_anomaly: 2.788812862316506,
        },
        uncertainty: None,
        covariance: None,
    };

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.038243152464060015,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // --- traj 29757 ---
    let best_orbit = full_orbit.remove(&"29757".into()).unwrap().unwrap();
    let orbit = best_orbit.orbital_elements();

    let expected_orbit = OrbitalElements::Keplerian {
        elements: KeplerianElements {
            reference_epoch: 60835.25573266984,
            semi_major_axis: 0.6359552205081981,
            eccentricity: 0.5904849543618848,
            inclination: 0.22635291292788373,
            ascending_node_longitude: 4.366539950357365,
            periapsis_argument: 3.3107966717676875,
            mean_anomaly: 3.0157533336404514,
        },
        uncertainty: None,
        covariance: None,
    };

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.025397386915472865,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    // --- traj 95777 ---
    let best_orbit = full_orbit.remove(&"95777".into()).unwrap().unwrap();
    let orbit = best_orbit.orbital_elements();

    let expected_orbit = OrbitalElements::Keplerian {
        elements: KeplerianElements {
            reference_epoch: 60894.252965553926,
            semi_major_axis: 1.247019899458668,
            eccentricity: 0.20820694219067884,
            inclination: 0.08116316043057617,
            ascending_node_longitude: 2.4955492264292034,
            periapsis_argument: 2.547031852373195,
            mean_anomaly: 0.29839367498266756,
        },
        uncertainty: None,
        covariance: None,
    };

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_orbit.orbit_quality(),
        0.014544406342677106,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}
