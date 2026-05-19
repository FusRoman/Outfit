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
use photom::observer::error_model::ObsErrorModel;
use rand::rngs::StdRng;
use rand::SeedableRng;

use crate::common::approx_equal;

#[test]
fn test_gauss_iod() {
    let test_max_relative = 1e-11;
    let test_epsilon = 1e-11;

    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("Download of the JPL short time scale UT1 data failed");

    let jpl_file: EphemFileSource = "horizon:DE440"
        .try_into()
        .expect("Failed to parse JPL ephemeris source");
    let jpl_ephem = JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Horizon");

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
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    let (best_orbit, best_rms) = full_orbit.remove(&"K09R05F".into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 57049.2684537375,
        semi_major_axis: 1.801740835743616,
        eccentricity: 0.28356259478492557,
        inclination: 0.2026828189979528,
        ascending_node_longitude: 0.007951791820548622,
        periapsis_argument: 1.2450647642587158,
        mean_anomaly: 0.4408048786626789,
    });

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        66.97479288637471,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60672.2443617134,
        semi_major_axis: 3.2199380906809876,
        eccentricity: 0.0624192099888107,
        inclination: 0.1829771029880289,
        ascending_node_longitude: 0.030775930195064964,
        periapsis_argument: 1.9053705720223801,
        mean_anomaly: 4.980622835177979,
    });

    let (best_orbit, best_rms) = full_orbit.remove(&8467_u32.into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        0.5739558189489471,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );

    let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
        reference_epoch: 60465.26777915681,
        semi_major_axis: 2.1874983804796972,
        eccentricity: 0.20256414489486008,
        inclination: 0.11906245183260411,
        ascending_node_longitude: 3.0918063960305293,
        periapsis_argument: 2.4793248309745692,
        mean_anomaly: 4.934465465531324,
    });

    let (best_orbit, best_rms) = full_orbit.remove(&33803_u32.into()).unwrap().unwrap();
    let orbit = best_orbit.get_orbit();

    assert!(approx_equal(&expected_orbit, orbit, test_epsilon));
    assert_relative_eq!(
        best_rms,
        18.963755528781288,
        epsilon = test_epsilon,
        max_relative = test_max_relative
    );
}
