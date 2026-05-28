mod common;

use crate::common::approx_equal;
use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::{
    jpl_ephem::download_jpl_file::EphemFileSource,
    orbit_type::{equinoctial_element::EquinoctialElements, OrbitalElements},
    DifferentialCorrectionConfig, FitLSQ, IODParams, JPLEphem,
};
use photom::TrajId;
use photom::{observation_dataset::ObsDataset, observer::error_model::ObsErrorModel};
use rand::{rngs::StdRng, SeedableRng};

fn build_test_fixtures() -> (
    JPLEphem,
    Ut1Provider,
    ObsDataset,
    IODParams,
    DifferentialCorrectionConfig,
) {
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

    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()
        .unwrap();

    // NOTE: rms_divergence_ratio is raised above the default (1.5) to allow
    // the 2-body differential corrector to handle objects whose osculating
    // 2-body elements differ noticeably from the N-body solution (e.g. 8467).
    //
    // Root cause: Outfit uses pure 2-body (Keplerian) propagation during
    // differential correction.  For well-observed main-belt asteroids the
    // best-fit 2-body elements differ from the best-fit N-body elements. When
    // the IOD starting point happens to be close to the N-body solution, the
    // first Newton step moves the orbit toward the 2-body minimum — which is
    // further from the data in N-body residuals — causing the RMS to jump by
    // more than the default 1.5× threshold before the iteration has a chance
    // to settle.  Setting the ratio to 10 lets the iteration continue through
    // that transient increase.
    //
    // Long-term fix: implement N-body propagation + variational equations in
    // the differential corrector (tracked separately).
    let diff_cor_config = DifferentialCorrectionConfig {
        rms_divergence_ratio: 10.0,
        ..DifferentialCorrectionConfig::default()
    };

    (
        jpl_ephem,
        ut1_provider,
        obs_dataset,
        iod_params,
        diff_cor_config,
    )
}

/// Non-regression test for differential orbit correction.
///
/// Oracle values were captured from a known-good Outfit run with seed 42.
/// Tolerances:
///   - Non-regression (Outfit vs oracle): 1e-10 absolute
#[test]
fn test_diff_cor() {
    let nr_tol = 1e-10;

    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, diff_cor_config) = build_test_fixtures();

    let full_orbit = obs_dataset
        .fit_lsq(
            &jpl_ephem,
            &ut1_provider,
            ObsErrorModel::FCCT14,
            &iod_params,
            &diff_cor_config,
            None,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    // -------------------------------------------------------------------------
    // 2015 AB  (MPC packed designation: K09R05F)
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::from("K09R05F"))
            .expect("K09R05F (2015AB) not found in results")
            .as_ref()
            .expect("K09R05F (2015AB) should converge");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 57049.2684537375,
            semi_major_axis: 1.801837227645679,
            eccentricity_sin_lon: 0.26941036025991355,
            eccentricity_cos_lon: 0.08909600747061494,
            tan_half_incl_sin_node: 0.0008708024189761142,
            tan_half_incl_cos_node: 0.10166598640878513,
            mean_longitude: 1.6929834276945714,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), nr_tol),
            "K09R05F orbital elements differ from oracle beyond tolerance {nr_tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 1.272e0, max_relative = 1e-3);
    }

    // -------------------------------------------------------------------------
    // 33803
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(33803))
            .expect("33803 not found in results")
            .as_ref()
            .expect("33803 should converge");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 60465.26777915681,
            semi_major_axis: 2.190614169340076,
            eccentricity_sin_lon: -0.13393967896355405,
            eccentricity_cos_lon: 0.1533932583177835,
            tan_half_incl_sin_node: 0.002997272576917091,
            tan_half_incl_cos_node: -0.05948928702443621,
            mean_longitude: 4.224671691074116,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), nr_tol),
            "33803 orbital elements differ from oracle beyond tolerance {nr_tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 4.344e-1, max_relative = 1e-3);
    }

    // -------------------------------------------------------------------------
    // 8467
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(8467))
            .expect("8467 not found in results")
            .as_ref()
            .expect("8467 should converge with rms_divergence_ratio=10");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 60672.2443617134,
            semi_major_axis: 3.2073734821020743,
            eccentricity_sin_lon: 0.053597752212361474,
            eccentricity_cos_lon: -0.023229330026225303,
            tan_half_incl_sin_node: 0.0028890355813102732,
            tan_half_incl_cos_node: 0.09179492536540514,
            mean_longitude: 0.626741395885302,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), nr_tol),
            "8467 orbital elements differ from oracle beyond tolerance {nr_tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 3.450e-1, max_relative = 1e-3);
    }
}
