mod common;

use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::initial_orbit_determination::IODParams;
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::obs_dataset::FitIOD;
use outfit::orbit_type::keplerian_element::KeplerianElements;
use outfit::orbit_type::OrbitalElements;
use outfit::JPLEphem;
use photom::observation_dataset::ObsDataset;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::common::approx_equal;

#[test]

fn test_gauss_iod() {
    let test_max_relative = 1e-11;
    let test_epsilon = 1e-11;

    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("Download of the JPL short time scale UT1 data failed");

    let jpl_file: EphemFileSource = "naif:DE440"
        .try_into()
        .expect("Failed to parse JPL ephemeris source");
    let jpl_ephem = JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Naif");

    let (obs_dataset, errors) = ObsDataset::from_mpc_80_col_files(&[
        "tests/data/2015AB.obs",
        "tests/data/8467.obs",
        "tests/data/33803.obs",
    ]);

    if !errors.is_empty() {
        panic!("Failed to load observation datasets: {:?}", errors);
    }

    let default = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130) // number of observation for the largest trajectory in the dataset
        .max_triplets(30)
        .build()
        .unwrap();

    let mut full_orbit = obs_dataset
        .fit_full_iod(
            &jpl_ephem,
            &ut1_provider,
            &default,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    let (best_orbit, best_rms) = full_orbit.remove(&"K09R05F".into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 57049.22904452732,
        semi_major_axis: 1.8017634341924542,
        eccentricity: 0.28360400982137396,
        inclination: 0.20267485730439427,
        ascending_node_longitude: 0.00810182022710516,
        periapsis_argument: 1.2445523487100616,
        mean_anomaly: 0.44069989140091426,
    });

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        47.67954270293223,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60672.24113100201,
        semi_major_axis: 3.189546977249391,
        eccentricity: 0.05434034666134485,
        inclination: 0.18343383575588465,
        ascending_node_longitude: 0.03253594968161228,
        periapsis_argument: 2.0197545218038355,
        mean_anomaly: 4.85070383704545,
    });

    let (best_orbit, best_rms) = full_orbit.remove(&"8467".into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        0.550927559734816,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60465.26778016307,
        semi_major_axis: 2.192136202201971,
        eccentricity: 0.2042936374305811,
        inclination: 0.1189651192106584,
        ascending_node_longitude: 3.091130251223283,
        periapsis_argument: 2.4714439663661487,
        mean_anomaly: 4.9466622638827324,
    });

    let (best_orbit, best_rms) = full_orbit.remove(&"33803".into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        6.319395085728921,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}
