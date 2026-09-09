mod common;

use crate::common::approx_equal;
use approx::assert_relative_eq;
use hifitime::ut1::Ut1Provider;
use outfit::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
};
use outfit::orbit_type::uncertainty::{EquinoctialUncertainty, OrbitalCovariance};
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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 57049.2684537375,
                semi_major_axis: 1.801837227645679,
                eccentricity_sin_lon: 0.26941036025991355,
                eccentricity_cos_lon: 0.08909600747061494,
                tan_half_incl_sin_node: 0.0008708024189761142,
                tan_half_incl_cos_node: 0.10166598640878513,
                mean_longitude: 1.6929834276945714,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 1.3935756201273647e-6,
                eccentricity_sin_lon: 2.399103573371585e-6,
                eccentricity_cos_lon: 9.380584628466963e-6,
                tan_half_incl_sin_node: 4.2486965596206456e-7,
                tan_half_incl_cos_node: 9.938054593077774e-7,
                mean_longitude: 1.5699462542222023e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        1.942053009013369e-12,
                        -3.7365542822268565e-13,
                        1.250111987715944e-11,
                        -3.8069560012308287e-13,
                        5.495356218939393e-13,
                        -2.1061628726935973e-11,
                    ],
                    [
                        -3.736554282226888e-13,
                        5.7556979557643085e-12,
                        -8.919579576942644e-12,
                        6.829258011452513e-13,
                        -2.190283688325579e-12,
                        1.4156679672214094e-11,
                    ],
                    [
                        1.2501119877159442e-11,
                        -8.919579576942621e-12,
                        8.799536797183067e-11,
                        -3.157563107997367e-12,
                        5.930188854586023e-12,
                        -1.472073140503015e-10,
                    ],
                    [
                        -3.806956001230829e-13,
                        6.829258011452509e-13,
                        -3.157563107997368e-12,
                        1.8051422455732311e-13,
                        -3.5751562142662264e-13,
                        5.229181995216352e-12,
                    ],
                    [
                        5.495356218939391e-13,
                        -2.1902836883255787e-12,
                        5.930188854586025e-12,
                        -3.5751562142662264e-13,
                        9.876492909499423e-13,
                        -9.67328953098736e-12,
                    ],
                    [
                        -2.1061628726935976e-11,
                        1.4156679672214063e-11,
                        -1.472073140503015e-10,
                        5.229181995216351e-12,
                        -9.673289530987361e-12,
                        2.464731241146324e-10,
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
                semi_major_axis: 2.190614169340076,
                eccentricity_sin_lon: -0.13393967896355405,
                eccentricity_cos_lon: 0.1533932583177835,
                tan_half_incl_sin_node: 0.002997272576917091,
                tan_half_incl_cos_node: -0.05948928702443621,
                mean_longitude: 4.224671691074116,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.1400421559849134e-5,
                eccentricity_sin_lon: 1.364670439647764e-5,
                eccentricity_cos_lon: 5.318530114145479e-6,
                tan_half_incl_sin_node: 3.44968775225327e-7,
                tan_half_incl_cos_node: 8.503880052285401e-7,
                mean_longitude: 2.664301205078454e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.5797804293925557e-10,
                        -2.443785426064791e-10,
                        7.203221689097433e-11,
                        -1.883169629832777e-12,
                        -6.3279112379918766e-12,
                        4.3441160814862357e-10,
                    ],
                    [
                        -2.443785426064796e-10,
                        1.8623254088484216e-10,
                        -6.032986816763725e-11,
                        7.999773867024745e-15,
                        -6.598752075412107e-13,
                        -3.5829528431457476e-10,
                    ],
                    [
                        7.203221689097439e-11,
                        -6.032986816763721e-11,
                        2.8286762575072326e-11,
                        2.0398130597296797e-14,
                        1.4218640626998597e-13,
                        1.2758725519460455e-10,
                    ],
                    [
                        -1.883169629832779e-12,
                        7.99977386702494e-15,
                        2.0398130597296844e-14,
                        1.190034558804622e-13,
                        2.64333826423024e-13,
                        3.756599803475119e-13,
                    ],
                    [
                        -6.327911237991877e-12,
                        -6.598752075412104e-13,
                        1.4218640626998607e-13,
                        2.64333826423024e-13,
                        7.231597594365756e-13,
                        2.605687909220327e-12,
                    ],
                    [
                        4.3441160814862383e-10,
                        -3.582952843145747e-10,
                        1.2758725519460457e-10,
                        3.7565998034751195e-13,
                        2.6056879092203274e-12,
                        7.098500911382502e-10,
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
                semi_major_axis: 3.2073734821020743,
                eccentricity_sin_lon: 0.053597752212361474,
                eccentricity_cos_lon: -0.023229330026225303,
                tan_half_incl_sin_node: 0.0028890355813102732,
                tan_half_incl_cos_node: 0.09179492536540514,
                mean_longitude: 0.626741395885302,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.00758317975106881,
                eccentricity_sin_lon: 0.002478406542589576,
                eccentricity_cos_lon: 0.0007443879537814839,
                tan_half_incl_sin_node: 4.277383244080703e-5,
                tan_half_incl_cos_node: 5.706392699913953e-5,
                mean_longitude: 0.00333399562783862,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.750461513702002e-5,
                        1.8729896457450725e-5,
                        5.604248768814215e-6,
                        -3.2370073744381016e-7,
                        -4.297318085854602e-7,
                        2.504633450274609e-5,
                    ],
                    [
                        1.8729896457450735e-5,
                        6.1424989903508165e-6,
                        1.8071841318216132e-6,
                        -1.0560687892019813e-7,
                        -1.409247502206143e-7,
                        8.250952263039232e-6,
                    ],
                    [
                        5.604248768814217e-6,
                        1.807184131821612e-6,
                        5.541134257349846e-7,
                        -3.14728840772654e-8,
                        -4.14717463955493e-8,
                        2.4005716002617356e-6,
                    ],
                    [
                        -3.237007374438101e-7,
                        -1.0560687892019811e-7,
                        -3.147288407726542e-8,
                        1.8296007416742358e-9,
                        2.435346888714026e-9,
                        -1.4137265325860534e-7,
                    ],
                    [
                        -4.2973180858546056e-7,
                        -1.4092475022061433e-7,
                        -4.1471746395549346e-8,
                        2.4353468887140264e-9,
                        3.2562917645631254e-9,
                        -1.8928599918199224e-7,
                    ],
                    [
                        2.50463345027461e-5,
                        8.250952263039232e-6,
                        2.400571600261738e-6,
                        -1.4137265325860537e-7,
                        -1.8928599918199224e-7,
                        1.1115526846447033e-5,
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
#[test]
#[ignore = "prints regenerated oracle literals for test_diff_cor_nbody_nonregression"]
fn dump_nbody_nonregression_oracle() {
    let (jpl_ephem, ut1_provider, obs_dataset, iod_params, _) = build_test_fixtures();

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
                semi_major_axis: 3.2064595210068236,
                eccentricity_sin_lon: 0.05302053705961858,
                eccentricity_cos_lon: -0.023191345016486484,
                tan_half_incl_sin_node: 0.002896547023047195,
                tan_half_incl_cos_node: 0.09180979390814888,
                mean_longitude: 0.6257180509530422,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.007572768801577769,
                eccentricity_sin_lon: 0.0024777810677302815,
                eccentricity_cos_lon: 0.0007445454670854451,
                tan_half_incl_sin_node: 4.2789214029110725e-05,
                tan_half_incl_cos_node: 5.709011408307015e-05,
                mean_longitude: 0.0033348986376892224,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.7346827322149607e-05,
                        1.8699274495858503e-05,
                        5.5978341299500165e-06,
                        -3.2337218267235595e-07,
                        -4.293370052067704e-07,
                        2.5018395598006037e-05,
                    ],
                    [
                        1.8699274495858496e-05,
                        6.139399019602613e-06,
                        1.8071180883575584e-06,
                        -1.0561813508866851e-07,
                        -1.4095413666601356e-07,
                        8.251115621205382e-06,
                    ],
                    [
                        5.597834129950014e-06,
                        1.8071180883575576e-06,
                        5.543479525574835e-07,
                        -3.1491072612759017e-08,
                        -4.1500199526702724e-08,
                        2.401759196797288e-06,
                    ],
                    [
                        -3.2337218267235584e-07,
                        -1.0561813508866851e-07,
                        -3.149107261275904e-08,
                        1.830916837229046e-09,
                        2.437345023372592e-09,
                        -1.414619819340877e-07,
                    ],
                    [
                        -4.293370052067702e-07,
                        -1.4095413666601353e-07,
                        -4.1500199526702724e-08,
                        2.4373450233725913e-09,
                        3.2592811260179643e-09,
                        -1.8942436671705862e-07,
                    ],
                    [
                        2.5018395598006003e-05,
                        8.25111562120538e-06,
                        2.401759196797288e-06,
                        -1.4146198193408767e-07,
                        -1.894243667170586e-07,
                        1.1121548923661431e-05,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "8467 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3487536466316906, epsilon = tol);
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
                semi_major_axis: 2.1903762390743804,
                eccentricity_sin_lon: -0.13376571752518818,
                eccentricity_cos_lon: 0.1533923223347731,
                tan_half_incl_sin_node: 0.00298829610339966,
                tan_half_incl_cos_node: -0.05950520161506636,
                mean_longitude: 4.224408175658702,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.138754063212835e-05,
                eccentricity_sin_lon: 1.3645893495668419e-05,
                eccentricity_cos_lon: 5.318606266540278e-06,
                tan_half_incl_sin_node: 3.4497740513000285e-07,
                tan_half_incl_cos_node: 8.504305789982595e-07,
                mean_longitude: 2.6646418667811342e-05,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.574268942909411e-10,
                        -2.4416518119419675e-10,
                        7.19676662428883e-11,
                        -1.8832309847367074e-12,
                        -6.328435623539356e-12,
                        4.3408659975282814e-10,
                    ],
                    [
                        -2.441651811941968e-10,
                        1.8621040929512565e-10,
                        -6.032837370428415e-11,
                        8.144705735767642e-15,
                        -6.596163586244179e-13,
                        -3.583210794623511e-10,
                    ],
                    [
                        7.19676662428883e-11,
                        -6.032837370428418e-11,
                        2.8287572618481517e-11,
                        2.036214125139651e-14,
                        1.4220215216291498e-13,
                        1.276056880086939e-10,
                    ],
                    [
                        -1.8832309847367062e-12,
                        8.144705735767625e-15,
                        2.0362141251396485e-14,
                        1.190094100502301e-13,
                        2.643549446442591e-13,
                        3.754660829157307e-13,
                    ],
                    [
                        -6.328435623539354e-12,
                        -6.596163586244183e-13,
                        1.4220215216291485e-13,
                        2.6435494464425914e-13,
                        7.232321696953148e-13,
                        2.6058390871950087e-12,
                    ],
                    [
                        4.3408659975282783e-10,
                        -3.583210794623512e-10,
                        1.2760568800869388e-10,
                        3.754660829157308e-13,
                        2.6058390871950087e-12,
                        7.100316278202847e-10,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "33803 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.6707182200067897, epsilon = tol);
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
                semi_major_axis: 1.8034373217966757,
                eccentricity_sin_lon: 0.2698207464026321,
                eccentricity_cos_lon: 0.08956625181322,
                tan_half_incl_sin_node: 0.0008477163489886412,
                tan_half_incl_cos_node: 0.10175513615418677,
                mean_longitude: 1.6918910074340627,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 1.8231056185204266e-06,
                eccentricity_sin_lon: 2.6968871000583087e-06,
                eccentricity_cos_lon: 1.2086907109458045e-05,
                tan_half_incl_sin_node: 5.433069007447511e-07,
                tan_half_incl_cos_node: 1.1352970467179496e-06,
                mean_longitude: 2.0262209200273614e-05,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        3.323714096280747e-12,
                        -5.43891761052487e-13,
                        2.128926441258226e-11,
                        -5.1182720087428e-13,
                        8.589495682097479e-13,
                        -3.58656751229932e-11,
                    ],
                    [
                        -5.438917610524938e-13,
                        7.273200030460913e-12,
                        -1.1779921957670923e-11,
                        8.870676736169013e-13,
                        -2.79563820045302e-12,
                        1.8776749774202017e-11,
                    ],
                    [
                        2.1289264412582264e-11,
                        -1.1779921957670865e-11,
                        1.4609332347266744e-10,
                        -4.232757191041452e-12,
                        8.577236986870777e-12,
                        -2.448322501242045e-10,
                    ],
                    [
                        -5.118272008742799e-13,
                        8.870676736169012e-13,
                        -4.232757191041451e-12,
                        2.9518238839686677e-13,
                        -4.835600243574041e-13,
                        7.023402307507248e-12,
                    ],
                    [
                        8.589495682097476e-13,
                        -2.795638200453019e-12,
                        8.577236986870775e-12,
                        -4.83560024357404e-13,
                        1.2888993842864982e-12,
                        -1.405977300032291e-11,
                    ],
                    [
                        -3.5865675122993194e-11,
                        1.8776749774201972e-11,
                        -2.448322501242045e-10,
                        7.023402307507247e-12,
                        -1.4059773000322908e-11,
                        4.105571216756527e-10,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "K09R05F N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3495925675201415, epsilon = tol);
    }
}
