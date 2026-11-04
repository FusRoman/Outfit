mod common;

use hifitime::ut1::Ut1Provider;
use outfit::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
};
use outfit::{
    orbit_type::{equinoctial_element::EquinoctialElements, OrbitalElements},
    propagator::{NBodyConfig, PropagatorKind},
    DifferentialCorrectionConfig, FitLSQ, IODParams, JPLEphem,
};

// Used only by the exact-oracle tests gated on the in-house ephemeris reader.
#[cfg(feature = "ephem-builtin")]
use crate::common::approx_equal;
#[cfg(feature = "ephem-builtin")]
use approx::assert_relative_eq;
#[cfg(feature = "ephem-builtin")]
use outfit::orbit_type::uncertainty::{EquinoctialUncertainty, OrbitalCovariance};
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
    build_test_fixtures_with(common::load_ephem())
}

/// Load the DE440 in-house Horizon reader explicitly.
///
/// Unlike [`common::load_ephem`], this always resolves to the legacy DE
/// binary reader (`ephem-builtin`'s [`horizon`](outfit::jpl_ephem::horizon)
/// backend), regardless of whether `ephem-anise` is also enabled — `naif:`
/// sources would instead be resolved by ANISE whenever it is compiled in
/// (see [`outfit::JPLEphem::new`]'s backend precedence). The exact-oracle
/// tests in this file need this: their reference values were captured
/// against this specific reader, and only agree with it to `1e-10` — ANISE
/// and the built-in reader agree with each other to `~1e-8` only.
#[cfg(feature = "ephem-builtin")]
fn load_horizon_ephem() -> JPLEphem {
    let source: outfit::jpl_ephem::download_jpl_file::EphemFileSource = "horizon:DE440"
        .try_into()
        .expect("failed to parse the Horizon ephemeris source");
    JPLEphem::from_builtin(source)
        .expect("failed to load the DE440 ephemeris through the in-house Horizon reader")
}

/// Same as [`build_test_fixtures`], but with an explicitly provided
/// ephemeris backend instead of [`common::load_ephem`]'s compile-time
/// default.
fn build_test_fixtures_with(
    jpl_ephem: JPLEphem,
) -> (
    JPLEphem,
    Ut1Provider,
    ObsDataset,
    IODParams,
    DifferentialCorrectionConfig,
) {
    let ut1_provider = common::load_ut1();

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
// Exact-oracle reproduction of the in-house ephemeris reader; the ANISE pipeline
// is covered by `test_diff_cor_nbody` (physical bounds) and the library parity
// tests.
#[cfg(feature = "ephem-builtin")]
#[test]
fn test_diff_cor() {
    let nr_tol = 1e-10;

    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, diff_cor_config) =
        build_test_fixtures_with(load_horizon_ephem());

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 57049.2684537375,
                semi_major_axis: 1.8018372276456807,
                eccentricity_sin_lon: 0.2694103602599234,
                eccentricity_cos_lon: 0.08909600747061278,
                tan_half_incl_sin_node: 0.0008708024189767214,
                tan_half_incl_cos_node: 0.10166598640878237,
                mean_longitude: 1.6929834276945717,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 1.3935756195557644e-6,
                eccentricity_sin_lon: 2.399103572387907e-6,
                eccentricity_cos_lon: 9.38058462462098e-6,
                tan_half_incl_sin_node: 4.2486965578790797e-7,
                tan_half_incl_cos_node: 9.93805458900432e-7,
                mean_longitude: 1.569946253578505e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        1.9420530074202325e-12,
                        -3.736554279179659e-13,
                        1.250111986690657e-11,
                        -3.8069559981098713e-13,
                        5.495356214438088e-13,
                        -2.1061628709661727e-11,
                    ],
                    [
                        -3.7365542791796506e-13,
                        5.755697951044417e-12,
                        -8.919579569638933e-12,
                        6.829258005855371e-13,
                        -2.190283686529886e-12,
                        1.4156679660622977e-11,
                    ],
                    [
                        1.2501119866906569e-11,
                        -8.919579569638912e-12,
                        8.799536789967552e-11,
                        -3.1575631054090046e-12,
                        5.930188849726796e-12,
                        -1.4720731392959058e-10,
                    ],
                    [
                        -3.806955998109871e-13,
                        6.829258005855374e-13,
                        -3.1575631054090054e-12,
                        1.805142244093354e-13,
                        -3.575156211335803e-13,
                        5.229181990929756e-12,
                    ],
                    [
                        5.495356214438094e-13,
                        -2.190283686529882e-12,
                        5.930188849726797e-12,
                        -3.5751562113358023e-13,
                        9.876492901402984e-13,
                        -9.673289523061097e-12,
                    ],
                    [
                        -2.1061628709661724e-11,
                        1.4156679660622993e-11,
                        -1.4720731392959063e-10,
                        5.229181990929757e-12,
                        -9.673289523061097e-12,
                        2.4647312391251836e-10,
                    ],
                ]
                .into(),
            }),
        };

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60465.26777915681,
                semi_major_axis: 2.1906141693400434,
                eccentricity_sin_lon: -0.13393967896355333,
                eccentricity_cos_lon: 0.15339325831776995,
                tan_half_incl_sin_node: 0.002997272576917366,
                tan_half_incl_cos_node: -0.059489287024434304,
                mean_longitude: 4.224671691074105,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.1400421559880274e-5,
                eccentricity_sin_lon: 1.364670439651079e-5,
                eccentricity_cos_lon: 5.318530114156113e-6,
                tan_half_incl_sin_node: 3.4496877522531336e-7,
                tan_half_incl_cos_node: 8.503880052285361e-7,
                mean_longitude: 2.6643012050850683e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.579780429405884e-10,
                        -2.443785426075781e-10,
                        7.203221689136443e-11,
                        -1.8831696298318545e-12,
                        -6.327911237985101e-12,
                        4.3441160815079203e-10,
                    ],
                    [
                        -2.4437854260757805e-10,
                        1.8623254088574694e-10,
                        -6.032986816795786e-11,
                        7.99977386634886e-15,
                        -6.598752075467102e-13,
                        -3.582952843163607e-10,
                    ],
                    [
                        7.203221689136459e-11,
                        -6.032986816795789e-11,
                        2.8286762575185433e-11,
                        2.0398130597537942e-14,
                        1.4218640627193858e-13,
                        1.2758725519523657e-10,
                    ],
                    [
                        -1.8831696298318565e-12,
                        7.999773866348862e-15,
                        2.039813059753784e-14,
                        1.1900345588045277e-13,
                        2.643338264230069e-13,
                        3.756599803488474e-13,
                    ],
                    [
                        -6.327911237985105e-12,
                        -6.598752075467105e-13,
                        1.421864062719383e-13,
                        2.6433382642300696e-13,
                        7.231597594365688e-13,
                        2.6056879092312455e-12,
                    ],
                    [
                        4.344116081507928e-10,
                        -3.582952843163608e-10,
                        1.2758725519523657e-10,
                        3.7565998034884743e-13,
                        2.6056879092312455e-12,
                        7.098500911417747e-10,
                    ],
                ]
                .into(),
            }),
        };

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60672.2443617134,
                semi_major_axis: 3.2073734821108264,
                eccentricity_sin_lon: 0.05359775221539591,
                eccentricity_cos_lon: -0.023229330025440753,
                tan_half_incl_sin_node: 0.0028890355812592043,
                tan_half_incl_cos_node: 0.09179492536533305,
                mean_longitude: 0.6267413958895128,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.007583180075782994,
                eccentricity_sin_lon: 0.0024784066499032027,
                eccentricity_cos_lon: 0.0007443879843659475,
                tan_half_incl_sin_node: 4.277383427235611e-5,
                tan_half_incl_cos_node: 5.706392945200011e-5,
                mean_longitude: 0.0033339957721846693,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.750462006175217e-5,
                        1.872989807598382e-5,
                        5.604249242358418e-6,
                        -3.237007652220501e-7,
                        -4.2973184571546797e-7,
                        2.5046336679932686e-5,
                    ],
                    [
                        1.872989807598381e-5,
                        6.142499522284416e-6,
                        1.8071842874556456e-6,
                        -1.0560688804961835e-7,
                        -1.4092476242345648e-7,
                        8.250952978571766e-6,
                    ],
                    [
                        5.6042492423584035e-6,
                        1.8071842874556456e-6,
                        5.54113471268398e-7,
                        -3.147288674832421e-8,
                        -4.147174996587513e-8,
                        2.400571809616114e-6,
                    ],
                    [
                        -3.237007652220494e-7,
                        -1.0560688804961831e-7,
                        -3.14728867483242e-8,
                        1.8296008983589866e-9,
                        2.4353470981477625e-9,
                        -1.4137266553913023e-7,
                    ],
                    [
                        -4.2973184571546765e-7,
                        -1.4092476242345642e-7,
                        -4.147174996587511e-8,
                        2.4353470981477625e-9,
                        3.2562920445028455e-9,
                        -1.8928601559669537e-7,
                    ],
                    [
                        2.5046336679932676e-5,
                        8.250952978571766e-6,
                        2.400571809616116e-6,
                        -1.4137266553913026e-7,
                        -1.8928601559669537e-7,
                        1.1115527808945248e-5,
                    ],
                ]
                .into(),
            }),
        };

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
            OrbitalElements::Equinoctial { elements, .. } => elements,
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
            OrbitalElements::Equinoctial { elements, .. } => elements,
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
            OrbitalElements::Equinoctial { elements, .. } => elements,
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

/// Helper: dumps a converged N-body orbit as a ready-to-paste oracle block.
/// Run with `cargo test --test test_diff_cor dump_nbody_nonregression_oracle -- --ignored --nocapture`.
#[cfg(feature = "ephem-builtin")]
#[test]
#[ignore = "prints regenerated oracle literals for test_diff_cor_nbody_nonregression"]
fn dump_nbody_nonregression_oracle() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, _) =
        build_test_fixtures_with(load_horizon_ephem());

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

    for label in ["8467", "33803", "K09R05F"] {
        let key = match label {
            "K09R05F" => TrajId::from("K09R05F"),
            n => TrajId::Int(n.parse().unwrap()),
        };
        let orbit = full_orbit.get(&key).unwrap().as_ref().unwrap();
        let (elements, unc, cov) = match orbit.orbital_elements() {
            OrbitalElements::Equinoctial {
                elements,
                uncertainty,
                covariance,
            } => (
                elements,
                uncertainty.as_ref().unwrap(),
                covariance.as_ref().unwrap(),
            ),
            _ => unreachable!(),
        };
        println!("// ---- {label} ----");
        println!("let expected = OrbitalElements::Equinoctial {{");
        println!("    elements: EquinoctialElements {{");
        println!("        reference_epoch: {:?},", elements.reference_epoch);
        println!("        semi_major_axis: {:?},", elements.semi_major_axis);
        println!(
            "        eccentricity_sin_lon: {:?},",
            elements.eccentricity_sin_lon
        );
        println!(
            "        eccentricity_cos_lon: {:?},",
            elements.eccentricity_cos_lon
        );
        println!(
            "        tan_half_incl_sin_node: {:?},",
            elements.tan_half_incl_sin_node
        );
        println!(
            "        tan_half_incl_cos_node: {:?},",
            elements.tan_half_incl_cos_node
        );
        println!("        mean_longitude: {:?},", elements.mean_longitude);
        println!("    }},");
        println!("    uncertainty: Some(EquinoctialUncertainty {{");
        println!("        semi_major_axis: {:?},", unc.semi_major_axis);
        println!(
            "        eccentricity_sin_lon: {:?},",
            unc.eccentricity_sin_lon
        );
        println!(
            "        eccentricity_cos_lon: {:?},",
            unc.eccentricity_cos_lon
        );
        println!(
            "        tan_half_incl_sin_node: {:?},",
            unc.tan_half_incl_sin_node
        );
        println!(
            "        tan_half_incl_cos_node: {:?},",
            unc.tan_half_incl_cos_node
        );
        println!("        mean_longitude: {:?},", unc.mean_longitude);
        println!("    }}),");
        println!("    covariance: Some(OrbitalCovariance {{");
        println!("        matrix: [");
        for r in 0..6 {
            println!("            [");
            for c in 0..6 {
                println!("                {:?},", cov.matrix[(r, c)]);
            }
            println!("            ],");
        }
        println!("        ]");
        println!("        .into(),");
        println!("    }}),");
        println!("}};");
        println!("// orbit_quality = {:?}", orbit.orbit_quality());
        println!();
    }
}

/// Strict non-regression test for the N-body differential corrector.
///
/// Reference values were captured from a deterministic run of
/// `dump_nbody_nonregression_oracle` and must remain reproducible to 1e-10.
// Exact-oracle reproduction of the in-house ephemeris reader.
#[cfg(feature = "ephem-builtin")]
#[test]
fn test_diff_cor_nbody_nonregression() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, _) =
        build_test_fixtures_with(load_horizon_ephem());

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60672.2443617134,
                semi_major_axis: 3.207570962866833,
                eccentricity_sin_lon: 0.05358176605779368,
                eccentricity_cos_lon: -0.02317371872117522,
                tan_half_incl_sin_node: 0.0028881645534142926,
                tan_half_incl_cos_node: 0.09179508427372225,
                mean_longitude: 0.6266478312948931,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.007583438399464079,
                eccentricity_sin_lon: 0.002478450802258219,
                eccentricity_cos_lon: 0.0007444024242159036,
                tan_half_incl_sin_node: 4.2771800144422136e-5,
                tan_half_incl_cos_node: 5.706349678441748e-5,
                mean_longitude: 0.003334126710926486,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.7508537958466306e-5,
                        1.8730855691322925e-5,
                        5.604551395320573e-6,
                        -3.2369629351180094e-7,
                        -4.2974302166896236e-7,
                        2.5048160817696977e-5,
                    ],
                    [
                        1.8730855691322932e-5,
                        6.14271837921441e-6,
                        1.8072495578763175e-6,
                        -1.0560370749732765e-7,
                        -1.4092617085115345e-7,
                        8.251426418071773e-6,
                    ],
                    [
                        5.6045513953205804e-6,
                        1.8072495578763184e-6,
                        5.541349691785142e-7,
                        -3.147198305954118e-8,
                        -4.147223680785181e-8,
                        2.4007126818504506e-6,
                    ],
                    [
                        -3.236962935118007e-7,
                        -1.0560370749732763e-7,
                        -3.147198305954116e-8,
                        1.8294268875943892e-9,
                        2.4352133443717015e-9,
                        -1.413715091783988e-7,
                    ],
                    [
                        -4.297430216689621e-7,
                        -1.4092617085115345e-7,
                        -4.147223680785182e-8,
                        2.435213344371702e-9,
                        3.256242665265224e-9,
                        -1.8929196285398248e-7,
                    ],
                    [
                        2.5048160817696977e-5,
                        8.251426418071773e-6,
                        2.4007126818504498e-6,
                        -1.413715091783988e-7,
                        -1.8929196285398253e-7,
                        1.1116400924513468e-5,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "8467 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3450677501261278, epsilon = tol);
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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60465.26777915681,
                semi_major_axis: 2.1907001009003797,
                eccentricity_sin_lon: -0.13400272182064238,
                eccentricity_cos_lon: 0.1533800913537971,
                tan_half_incl_sin_node: 0.002998523845120412,
                tan_half_incl_cos_node: -0.05948593480260023,
                mean_longitude: 4.224751249408136,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.1404880736301646e-5,
                eccentricity_sin_lon: 1.3646865441089579e-5,
                eccentricity_cos_lon: 5.318472489569789e-6,
                tan_half_incl_sin_node: 3.4495331038910287e-7,
                tan_half_incl_cos_node: 8.503524434875798e-7,
                mean_longitude: 2.6641390855810307e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.581689193352973e-10,
                        -2.4444869402377964e-10,
                        7.205169375761403e-11,
                        -1.882944769574507e-12,
                        -6.327633385709223e-12,
                        4.3450961651635614e-10,
                    ],
                    [
                        -2.444486940237796e-10,
                        1.8623693636720508e-10,
                        -6.03298981972857e-11,
                        7.881985349776607e-15,
                        -6.599946297549495e-13,
                        -3.582770230881789e-10,
                    ],
                    [
                        7.205169375761404e-11,
                        -6.03298981972857e-11,
                        2.828614962231067e-11,
                        2.0470448604142446e-14,
                        1.4233196526440624e-13,
                        1.2757883512903262e-10,
                    ],
                    [
                        -1.882944769574508e-12,
                        7.881985349776709e-15,
                        2.047044860414234e-14,
                        1.1899278634840073e-13,
                        2.6430760135856e-13,
                        3.758769572772045e-13,
                    ],
                    [
                        -6.3276333857092236e-12,
                        -6.599946297549496e-13,
                        1.4233196526440609e-13,
                        2.6430760135856e-13,
                        7.230992781452975e-13,
                        2.605828924776723e-12,
                    ],
                    [
                        4.345096165163567e-10,
                        -3.58277023088179e-10,
                        1.2757883512903262e-10,
                        3.758769572772045e-13,
                        2.605828924776723e-12,
                        7.09763706732053e-10,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "33803 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.4087745379548085, epsilon = tol);
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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 57049.2684537375,
                semi_major_axis: 1.8017539671066376,
                eccentricity_sin_lon: 0.26939702471014915,
                eccentricity_cos_lon: 0.0886949242110166,
                tan_half_incl_sin_node: 0.0008307916711973057,
                tan_half_incl_cos_node: 0.1016705971606363,
                mean_longitude: 1.6936595589339083,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 8.918432978038933e-7,
                eccentricity_sin_lon: 1.4528052999506174e-6,
                eccentricity_cos_lon: 5.920319266943404e-6,
                tan_half_incl_sin_node: 2.792168578053201e-7,
                tan_half_incl_cos_node: 6.218091729545387e-7,
                mean_longitude: 9.896911634400027e-6,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        7.953844678377238e-13,
                        -1.0232494678801584e-13,
                        5.059210617758408e-12,
                        -1.506692655121542e-13,
                        2.4607703529115327e-13,
                        -8.522632811132329e-12,
                    ],
                    [
                        -1.0232494678801523e-13,
                        2.1106432395646035e-12,
                        -3.0905811444994442e-12,
                        2.4837164956545964e-13,
                        -7.926397586796095e-13,
                        4.808599014209342e-12,
                    ],
                    [
                        5.059210617758407e-12,
                        -3.0905811444994515e-12,
                        3.5050180222541284e-11,
                        -1.2295382103683321e-12,
                        2.4527485767349156e-12,
                        -5.855948314989829e-11,
                    ],
                    [
                        -1.5066926551215407e-13,
                        2.4837164956546045e-13,
                        -1.2295382103683325e-12,
                        7.796205368267633e-14,
                        -1.454002790606985e-13,
                        2.0297758780287726e-12,
                    ],
                    [
                        2.4607703529115327e-13,
                        -7.9263975867961e-13,
                        2.4527485767349156e-12,
                        -1.454002790606985e-13,
                        3.8664664757040743e-13,
                        -3.9888503465101905e-12,
                    ],
                    [
                        -8.522632811132327e-12,
                        4.808599014209346e-12,
                        -5.855948314989828e-11,
                        2.0297758780287726e-12,
                        -3.9888503465101905e-12,
                        9.794885989912263e-11,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "K09R05F N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.6128100688234707, epsilon = tol);
    }
}
