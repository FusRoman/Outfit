mod common;

use crate::common::approx_equal;
use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
};
use outfit::{
    orbit_type::{equinoctial_element::EquinoctialElements, OrbitalElements},
    propagator::{NBodyConfig, PropagatorKind},
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

    let jpl_ephem: JPLEphem = "horizon:DE440"
        .try_into()
        .expect("Failed to load JPL ephemeris");

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

/// N-body differential orbit correction test.
///
/// Runs the same dataset with Sun + Jupiter as perturbers and verifies:
///
/// 1. All three objects still converge.
/// 2. The fitted semi-major axes are physically consistent with the two-body
///    oracle (no sign flip, reasonable magnitude).
/// 3. For the shortest arc (8467, ~40 days) the N-body elements agree with
///    the two-body oracle within 5e-2 AU / 5e-2 (eccentricity components).
///    Jupiter's perturbation over 40 days at 3.2 AU is measurable but small.
/// 4. The orbit quality (normalised RMS) remains below a generous bound of 5.0
///    for all objects, confirming the fit converged to a good residual level.
///
/// The test intentionally does **not** demand tight element agreement for the
/// longer arcs (2015AB at ~5 years, 33803 at ~5 months), where the best-fit
/// N-body and two-body elements differ by physical amounts (Jovian
/// perturbations accumulate over the arc).  What it does verify is that the
/// propagator produces a self-consistent, converged orbit.
#[test]
fn test_diff_cor_nbody() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, _) = build_test_fixtures();

    // Perturbers: Sun (central body) + Jupiter (dominant perturber in the
    // main belt).  Using the same rms_divergence_ratio as the 2-body test
    // because the N-body corrector converges smoothly to its own minimum.
    let nbody_config = NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::PB(PlanetaryBary::Jupiter),
        ],
        abs_tol: 1e-12,
        rel_tol: 1e-12,
    };

    let diff_cor_config = DifferentialCorrectionConfig {
        rms_divergence_ratio: 10.0,
        propagator: PropagatorKind::NBody(nbody_config),
        ..DifferentialCorrectionConfig::default()
    };

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

    // 2-body oracle values (from test_diff_cor) used as reference below.
    let twobody_sma_k09r05f = 1.801837227645679_f64;
    let twobody_sma_33803 = 2.190614169340076_f64;
    let twobody_sma_8467 = 3.2073734821020743_f64;

    // -------------------------------------------------------------------------
    // 2015 AB  (K09R05F) — 5-year arc, larger N-body/2-body divergence expected
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::from("K09R05F"))
            .expect("K09R05F not found in N-body results")
            .as_ref()
            .expect("K09R05F should converge under N-body propagation");

        let elem = match orbit.orbital_elements() {
            OrbitalElements::Equinoctial(e) => e,
            _ => panic!("Expected equinoctial elements for K09R05F"),
        };

        // Semi-major axis must be positive and within 0.3 AU of the 2-body value.
        assert!(
            elem.semi_major_axis > 0.0,
            "K09R05F N-body a must be positive"
        );
        assert!(
            (elem.semi_major_axis - twobody_sma_k09r05f).abs() < 0.3,
            "K09R05F N-body a = {} differs from 2-body oracle {} by more than 0.3 AU",
            elem.semi_major_axis,
            twobody_sma_k09r05f
        );

        // Fit quality must remain physically reasonable.
        assert!(
            orbit.orbit_quality() < 5.0,
            "K09R05F N-body orbit quality {} exceeds bound 5.0",
            orbit.orbit_quality()
        );
    }

    // -------------------------------------------------------------------------
    // 33803 — 5-month arc
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(33803))
            .expect("33803 not found in N-body results")
            .as_ref()
            .expect("33803 should converge under N-body propagation");

        let elem = match orbit.orbital_elements() {
            OrbitalElements::Equinoctial(e) => e,
            _ => panic!("Expected equinoctial elements for 33803"),
        };

        assert!(
            elem.semi_major_axis > 0.0,
            "33803 N-body a must be positive"
        );
        assert!(
            (elem.semi_major_axis - twobody_sma_33803).abs() < 0.1,
            "33803 N-body a = {} differs from 2-body oracle {} by more than 0.1 AU",
            elem.semi_major_axis,
            twobody_sma_33803
        );
        assert!(
            orbit.orbit_quality() < 5.0,
            "33803 N-body orbit quality {} exceeds bound 5.0",
            orbit.orbit_quality()
        );
    }

    // -------------------------------------------------------------------------
    // 8467 — 40-day arc (tightest comparison with 2-body oracle)
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(8467))
            .expect("8467 not found in N-body results")
            .as_ref()
            .expect("8467 should converge under N-body propagation");

        let elem = match orbit.orbital_elements() {
            OrbitalElements::Equinoctial(e) => e,
            _ => panic!("Expected equinoctial elements for 8467"),
        };

        let twobody_8467 = EquinoctialElements {
            reference_epoch: 60672.2443617134,
            semi_major_axis: twobody_sma_8467,
            eccentricity_sin_lon: 0.053597752212361474,
            eccentricity_cos_lon: -0.023229330026225303,
            tan_half_incl_sin_node: 0.0028890355813102732,
            tan_half_incl_cos_node: 0.09179492536540514,
            mean_longitude: 0.626741395885302,
        };

        // Over a 40-day arc at 3.2 AU, Jovian perturbations produce element
        // changes well below 5e-2 in dimensionless units.
        let a_tol = 5e-2; // AU
        let e_tol = 5e-2; // eccentricity components (dimensionless)

        assert!(
            (elem.semi_major_axis - twobody_8467.semi_major_axis).abs() < a_tol,
            "8467 N-body a = {} differs from 2-body oracle {} by more than {} AU",
            elem.semi_major_axis,
            twobody_8467.semi_major_axis,
            a_tol
        );
        assert!(
            (elem.eccentricity_sin_lon - twobody_8467.eccentricity_sin_lon).abs() < e_tol,
            "8467 N-body h = {} differs from 2-body oracle {} by more than {}",
            elem.eccentricity_sin_lon,
            twobody_8467.eccentricity_sin_lon,
            e_tol
        );
        assert!(
            (elem.eccentricity_cos_lon - twobody_8467.eccentricity_cos_lon).abs() < e_tol,
            "8467 N-body k = {} differs from 2-body oracle {} by more than {}",
            elem.eccentricity_cos_lon,
            twobody_8467.eccentricity_cos_lon,
            e_tol
        );
        assert!(
            orbit.orbit_quality() < 5.0,
            "8467 N-body orbit quality {} exceeds bound 5.0",
            orbit.orbit_quality()
        );
    }
}

/// Strict non-regression test for the N-body differential corrector.
///
/// Reference values were captured from a deterministic run of `test_diff_cor_nbody`
/// with `--nocapture` and must remain reproducible to 1e-10.
#[test]
fn test_diff_cor_nbody_nonregression() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, _) = build_test_fixtures();

    let nbody_config = NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::PB(PlanetaryBary::Jupiter),
        ],
        abs_tol: 1e-12,
        rel_tol: 1e-12,
    };

    let diff_cor_config = DifferentialCorrectionConfig {
        rms_divergence_ratio: 10.0,
        propagator: PropagatorKind::NBody(nbody_config),
        ..DifferentialCorrectionConfig::default()
    };

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

    let tol = 1e-10_f64;

    // -------------------------------------------------------------------------
    // 8467 — reference N-body final equinoctial elements + orbit quality
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(8467))
            .expect("8467 not found")
            .as_ref()
            .expect("8467 should converge");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 60672.2443617134,
            semi_major_axis: 3.2064058028481552,
            eccentricity_sin_lon: 0.05300520970081429,
            eccentricity_cos_lon: -0.023197692700634636,
            tan_half_incl_sin_node: 0.002896813138792025,
            tan_half_incl_cos_node: 0.09181010554057693,
            mean_longitude: 0.6256995904459722,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "8467 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3486122845933199, epsilon = tol);
    }

    // -------------------------------------------------------------------------
    // 33803 — reference N-body final equinoctial elements + orbit quality
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::Int(33803))
            .expect("33803 not found")
            .as_ref()
            .expect("33803 should converge");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 60465.26777915681,
            semi_major_axis: 2.190348311458185,
            eccentricity_sin_lon: -0.13373910921857446,
            eccentricity_cos_lon: 0.15339157238172804,
            tan_half_incl_sin_node: 0.0029876412023201416,
            tan_half_incl_cos_node: -0.05950692044872062,
            mean_longitude: 4.224365041422834,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "33803 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.7034091187041202, epsilon = tol);
    }

    // -------------------------------------------------------------------------
    // K09R05F — reference N-body final equinoctial elements + orbit quality
    // -------------------------------------------------------------------------
    {
        let orbit = full_orbit
            .get(&TrajId::from("K09R05F"))
            .expect("K09R05F not found")
            .as_ref()
            .expect("K09R05F should converge");

        let expected = OrbitalElements::Equinoctial(EquinoctialElements {
            reference_epoch: 57049.2684537375,
            semi_major_axis: 1.8021517900042052,
            eccentricity_sin_lon: 0.2694922786015968,
            eccentricity_cos_lon: 0.08955282358108035,
            tan_half_incl_sin_node: 0.0008974287327937245,
            tan_half_incl_cos_node: 0.10167548786557225,
            mean_longitude: 1.6921653421358704,
        });

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "K09R05F N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3608868439717083, epsilon = tol);
    }
}
