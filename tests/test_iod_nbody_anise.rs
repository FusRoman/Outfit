//! Non-regression and cross-backend consistency tests for the Gauss IOD and
//! N-body differential correction pipelines on the **ANISE** backend.
//!
//! `tests/test_gauss_iod.rs` and `tests/test_diff_cor.rs` pin the same
//! reference trajectories (2015AB/K09R05F, 8467, 33803) to the in-house
//! Horizon reader's arithmetic bit-for-bit. This file exercises the same
//! pipelines through ANISE instead, and checks two properties for each
//! trajectory:
//!
//! 1. **Own-backend non-regression**: the fitted elements and orbit quality
//!    reproduce a captured, deterministic ANISE run bit-for-bit (`1e-10`).
//! 2. **Cross-backend consistency**: those same ANISE-fitted elements agree
//!    with the Horizon-backend oracle to a much looser tolerance. ANISE and
//!    the in-house reader parse the same DE440 data through independent code
//!    paths and agree on raw positions to only `~1e-8` AU (see the backend
//!    parity tests in `src/jpl_ephem`), so the fitted elements are expected
//!    to differ slightly — this check guards against a large, real
//!    disagreement rather than requiring identity.
#![cfg(feature = "ephem-anise")]

mod common;

use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
};
use outfit::orbit_type::{
    equinoctial_element::EquinoctialElements, keplerian_element::KeplerianElements, OrbitalElements,
};
use outfit::{
    propagator::{NBodyConfig, PropagatorKind},
    DifferentialCorrectionConfig, FitIOD, FitLSQ, IODParams, JPLEphem,
};
use photom::TrajId;
use photom::{observation_dataset::ObsDataset, observer::error_model::ObsErrorModel};
use rand::{rngs::StdRng, SeedableRng};

/// Tolerance for cross-backend element comparison (absolute, same units as
/// the element: AU for semi-major axis, radians for angles).
///
/// ANISE and the in-house Horizon reader agree on positions to `~1e-8` AU.
/// Empirically, the largest cross-backend element disagreement observed
/// across the fixtures in this file is `~3e-7` (Gauss IOD, 33803's mean
/// anomaly) — this tolerance keeps a comfortable margin over that without
/// being so loose it would miss a real regression.
const CROSS_BACKEND_ELEMENT_TOL: f64 = 1e-6;

/// Relative tolerance for cross-backend orbit-quality (RMS) comparison.
///
/// Observed cross-backend relative differences are `~2e-4` at most (Gauss
/// IOD, 33803); N-body orbit quality agrees far more closely (`~1e-6`)
/// since the differential corrector converges to a shared physical minimum
/// largely independent of ephemeris precision.
const CROSS_BACKEND_RMS_MAX_RELATIVE: f64 = 1e-3;

/// Tolerance for the own-ANISE non-regression oracle: this backend's
/// arithmetic is deterministic, so this should hold to numerical noise.
const ANISE_ORACLE_TOL: f64 = 1e-10;

fn load_anise_ephem() -> JPLEphem {
    let source: EphemFileSource = "naif:DE440"
        .try_into()
        .expect("failed to parse the NAIF ephemeris source");
    JPLEphem::from_anise(source).expect("failed to load the DE440 ephemeris through ANISE")
}

fn build_test_fixtures() -> (JPLEphem, Ut1Provider, ObsDataset, IODParams) {
    let ut1_provider = common::load_ut1();

    let jpl_ephem = load_anise_ephem();

    let (obs_dataset, errors) = ObsDataset::from_mpc_80_col_files(&[
        "tests/data/2015AB.obs",
        "tests/data/8467.obs",
        "tests/data/33803.obs",
    ]);
    if !errors.is_empty() {
        panic!("Failed to load observation datasets: {errors:?}");
    }

    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()
        .unwrap();

    (jpl_ephem, ut1_provider, obs_dataset, iod_params)
}

/// Compare two Keplerian element sets component-wise, each within `tol`.
fn keplerian_close(a: &KeplerianElements, b: &KeplerianElements, tol: f64) -> bool {
    (a.semi_major_axis - b.semi_major_axis).abs() < tol
        && (a.eccentricity - b.eccentricity).abs() < tol
        && (a.inclination - b.inclination).abs() < tol
        && (a.ascending_node_longitude - b.ascending_node_longitude).abs() < tol
        && (a.periapsis_argument - b.periapsis_argument).abs() < tol
        && (a.mean_anomaly - b.mean_anomaly).abs() < tol
}

/// Compare two equinoctial element sets component-wise, each within `tol`.
fn equinoctial_close(a: &EquinoctialElements, b: &EquinoctialElements, tol: f64) -> bool {
    (a.semi_major_axis - b.semi_major_axis).abs() < tol
        && (a.eccentricity_sin_lon - b.eccentricity_sin_lon).abs() < tol
        && (a.eccentricity_cos_lon - b.eccentricity_cos_lon).abs() < tol
        && (a.tan_half_incl_sin_node - b.tan_half_incl_sin_node).abs() < tol
        && (a.tan_half_incl_cos_node - b.tan_half_incl_cos_node).abs() < tol
        && (a.mean_longitude - b.mean_longitude).abs() < tol
}

#[test]
fn gauss_iod_anise_matches_horizon_within_tolerance() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params) = build_test_fixtures();

    let mut full_orbit = obs_dataset
        .fit_full_iod(
            &jpl_ephem,
            &ut1_provider,
            &iod_params,
            ObsErrorModel::FCCT14,
            &mut StdRng::seed_from_u64(42),
        )
        .unwrap();

    // ---- K09R05F ----
    {
        let best_orbit = full_orbit.remove(&"K09R05F".into()).unwrap().unwrap();
        let OrbitalElements::Keplerian {
            elements: anise, ..
        } = best_orbit.orbital_elements()
        else {
            panic!("expected Keplerian elements");
        };

        let anise_oracle = KeplerianElements {
            reference_epoch: 57049.268453737495,
            semi_major_axis: 1.8017408423118564,
            eccentricity: 0.283562596550713,
            inclination: 0.20268281960979226,
            ascending_node_longitude: 0.007951790269219904,
            periapsis_argument: 1.2450647682683995,
            mean_anomaly: 0.4408048750959458,
        };
        assert!(
            keplerian_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "K09R05F (ANISE, Gauss IOD): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            66.95944969907143,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = KeplerianElements {
            reference_epoch: 57049.2684537375,
            semi_major_axis: 1.8017408357435853,
            eccentricity: 0.2835625947849,
            inclination: 0.20268281899795498,
            ascending_node_longitude: 0.007951791820480579,
            periapsis_argument: 1.2450647642589465,
            mean_anomaly: 0.44080487866260754,
        };
        assert!(
            keplerian_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "K09R05F (ANISE, Gauss IOD): elements {anise:?} disagree with the Horizon-backend \
             oracle {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            66.97479294218962,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }

    // ---- 8467 ----
    {
        let best_orbit = full_orbit.remove(&8467_u32.into()).unwrap().unwrap();
        let OrbitalElements::Keplerian {
            elements: anise, ..
        } = best_orbit.orbital_elements()
        else {
            panic!("expected Keplerian elements");
        };

        let anise_oracle = KeplerianElements {
            reference_epoch: 60672.24436171338,
            semi_major_axis: 3.2199381014940345,
            eccentricity: 0.06241921076866068,
            inclination: 0.18297710293694305,
            ascending_node_longitude: 0.030775929741641336,
            periapsis_argument: 1.9053705337483742,
            mean_anomaly: 4.980622874255642,
        };
        assert!(
            keplerian_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "8467 (ANISE, Gauss IOD): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            0.5739558166838535,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = KeplerianElements {
            reference_epoch: 60672.2443617134,
            semi_major_axis: 3.2199380906766213,
            eccentricity: 0.062419209987893194,
            inclination: 0.1829771029880662,
            ascending_node_longitude: 0.030775930195269436,
            periapsis_argument: 1.9053705720361855,
            mean_anomaly: 4.980622835162594,
        };
        assert!(
            keplerian_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "8467 (ANISE, Gauss IOD): elements {anise:?} disagree with the Horizon-backend \
             oracle {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            0.5739558188544629,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }

    // ---- 33803 ----
    {
        let best_orbit = full_orbit.remove(&33803_u32.into()).unwrap().unwrap();
        let OrbitalElements::Keplerian {
            elements: anise, ..
        } = best_orbit.orbital_elements()
        else {
            panic!("expected Keplerian elements");
        };

        let anise_oracle = KeplerianElements {
            reference_epoch: 60465.26777915684,
            semi_major_axis: 2.1874982448104126,
            eccentricity: 0.2025641009051343,
            inclination: 0.11906245470092386,
            ascending_node_longitude: 3.091806419129847,
            periapsis_argument: 2.4793250506440367,
            mean_anomaly: 4.9344651309211205,
        };
        assert!(
            keplerian_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "33803 (ANISE, Gauss IOD): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            18.96454770398426,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = KeplerianElements {
            reference_epoch: 60465.26777915681,
            semi_major_axis: 2.1874983806597665,
            eccentricity: 0.2025641449546787,
            inclination: 0.11906245182869549,
            ascending_node_longitude: 3.091806395999626,
            periapsis_argument: 2.479324830682491,
            mean_anomaly: 4.934465465979427,
        };
        assert!(
            keplerian_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "33803 (ANISE, Gauss IOD): elements {anise:?} disagree with the Horizon-backend \
             oracle {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            best_orbit.orbit_quality(),
            18.96375448637929,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }
}

#[test]
fn nbody_diff_cor_anise_matches_horizon_within_tolerance() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params) = build_test_fixtures();

    let nbody_config = NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::PB(PlanetaryBary::Jupiter),
        ],
        abs_tol: 1e-12,
        rel_tol: 1e-12,
        ..NBodyConfig::default()
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

    // ---- 8467 ----
    {
        let orbit = full_orbit
            .get(&TrajId::Int(8467))
            .expect("8467 not found")
            .as_ref()
            .expect("8467 should converge");
        let OrbitalElements::Equinoctial {
            elements: anise, ..
        } = orbit.orbital_elements()
        else {
            panic!("expected equinoctial elements");
        };

        let anise_oracle = EquinoctialElements {
            reference_epoch: 60672.24436171338,
            semi_major_axis: 3.207570968343197,
            eccentricity_sin_lon: 0.053581766035026554,
            eccentricity_cos_lon: -0.023173717297720883,
            tan_half_incl_sin_node: 0.002888164539958646,
            tan_half_incl_cos_node: 0.09179508428588358,
            mean_longitude: 0.6266478295853958,
        };
        assert!(
            equinoctial_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "8467 (ANISE, N-body): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.34506775156504055,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = EquinoctialElements {
            reference_epoch: 60672.2443617134,
            semi_major_axis: 3.207570962866833,
            eccentricity_sin_lon: 0.05358176605779368,
            eccentricity_cos_lon: -0.02317371872117522,
            tan_half_incl_sin_node: 0.0028881645534142926,
            tan_half_incl_cos_node: 0.09179508427372225,
            mean_longitude: 0.6266478312948931,
        };
        assert!(
            equinoctial_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "8467 (ANISE, N-body): elements {anise:?} disagree with the Horizon-backend oracle \
             {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.3450677501261278,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }

    // ---- 33803 ----
    {
        let orbit = full_orbit
            .get(&TrajId::Int(33803))
            .expect("33803 not found")
            .as_ref()
            .expect("33803 should converge");
        let OrbitalElements::Equinoctial {
            elements: anise, ..
        } = orbit.orbital_elements()
        else {
            panic!("expected equinoctial elements");
        };

        let anise_oracle = EquinoctialElements {
            reference_epoch: 60465.26777915684,
            semi_major_axis: 2.1907001018465175,
            eccentricity_sin_lon: -0.1340027234565689,
            eccentricity_cos_lon: 0.1533800919526236,
            tan_half_incl_sin_node: 0.002998523851621838,
            tan_half_incl_cos_node: -0.05948593474189358,
            mean_longitude: 4.224751252995968,
        };
        assert!(
            equinoctial_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "33803 (ANISE, N-body): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.40877494651441326,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = EquinoctialElements {
            reference_epoch: 60465.26777915681,
            semi_major_axis: 2.1907001009003797,
            eccentricity_sin_lon: -0.13400272182064238,
            eccentricity_cos_lon: 0.1533800913537971,
            tan_half_incl_sin_node: 0.002998523845120412,
            tan_half_incl_cos_node: -0.05948593480260023,
            mean_longitude: 4.224751249408136,
        };
        assert!(
            equinoctial_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "33803 (ANISE, N-body): elements {anise:?} disagree with the Horizon-backend oracle \
             {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.4087745379548085,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }

    // ---- K09R05F ----
    {
        let orbit = full_orbit
            .get(&TrajId::from("K09R05F"))
            .expect("K09R05F not found")
            .as_ref()
            .expect("K09R05F should converge");
        let OrbitalElements::Equinoctial {
            elements: anise, ..
        } = orbit.orbital_elements()
        else {
            panic!("expected equinoctial elements");
        };

        let anise_oracle = EquinoctialElements {
            reference_epoch: 57049.268453737495,
            semi_major_axis: 1.8017539669836284,
            eccentricity_sin_lon: 0.269397024892828,
            eccentricity_cos_lon: 0.0886949233352466,
            tan_half_incl_sin_node: 0.0008307916904969703,
            tan_half_incl_cos_node: 0.10167059708040241,
            mean_longitude: 1.6936595604040174,
        };
        assert!(
            equinoctial_close(anise, &anise_oracle, ANISE_ORACLE_TOL),
            "K09R05F (ANISE, N-body): elements {anise:?} differ from the ANISE oracle beyond {ANISE_ORACLE_TOL}"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.6128051633992085,
            epsilon = ANISE_ORACLE_TOL
        );

        let horizon_oracle = EquinoctialElements {
            reference_epoch: 57049.2684537375,
            semi_major_axis: 1.8017539671066376,
            eccentricity_sin_lon: 0.26939702471014915,
            eccentricity_cos_lon: 0.0886949242110166,
            tan_half_incl_sin_node: 0.0008307916711973057,
            tan_half_incl_cos_node: 0.1016705971606363,
            mean_longitude: 1.6936595589339083,
        };
        assert!(
            equinoctial_close(anise, &horizon_oracle, CROSS_BACKEND_ELEMENT_TOL),
            "K09R05F (ANISE, N-body): elements {anise:?} disagree with the Horizon-backend oracle \
             {horizon_oracle:?} by more than {CROSS_BACKEND_ELEMENT_TOL} — expected close \
             agreement between backends, not this far apart"
        );
        assert_relative_eq!(
            orbit.orbit_quality(),
            0.6128100688234707,
            max_relative = CROSS_BACKEND_RMS_MAX_RELATIVE
        );
    }
}
