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
    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("Download of the JPL short time scale UT1 data failed");

    let jpl_ephem: JPLEphem = common::load_ephem();

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
#[cfg(feature = "ephem-builtin")]
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
// Exact-oracle reproduction of the in-house ephemeris reader.
#[cfg(feature = "ephem-builtin")]
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
                semi_major_axis: 3.2075709628598497,
                eccentricity_sin_lon: 0.05358176605531737,
                eccentricity_cos_lon: -0.023173718721778213,
                tan_half_incl_sin_node: 0.0028881645534556155,
                tan_half_incl_cos_node: 0.09179508427378168,
                mean_longitude: 0.6266478312914119,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.0075834386193465685,
                eccentricity_sin_lon: 0.0024784508747146784,
                eccentricity_cos_lon: 0.0007444024450377598,
                tan_half_incl_sin_node: 4.277180138435706e-5,
                tan_half_incl_cos_node: 5.7063498440961735e-5,
                mean_longitude: 0.00333412680819657,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.7508541293396986e-5,
                        1.873085678575111e-5,
                        5.604551716850951e-6,
                        -3.2369631231954327e-7,
                        -4.297430467785363e-7,
                        2.5048162288448977e-5,
                    ],
                    [
                        1.873085678575113e-5,
                        6.1427187383739545e-6,
                        1.8072496633933745e-6,
                        -1.0560371366948271e-7,
                        -1.4092617909138465e-7,
                        8.251426900729693e-6,
                    ],
                    [
                        5.60455171685095e-6,
                        1.8072496633933745e-6,
                        5.54135000178195e-7,
                        -3.147198487285182e-8,
                        -4.1472239228740826e-8,
                        2.4007128236500167e-6,
                    ],
                    [
                        -3.236963123195433e-7,
                        -1.0560371366948272e-7,
                        -3.147198487285183e-8,
                        1.829426993662889e-9,
                        2.4352134859800713e-9,
                        -1.4137151747287231e-7,
                    ],
                    [
                        -4.297430467785363e-7,
                        -1.4092617909138465e-7,
                        -4.147223922874081e-8,
                        2.435213485980071e-9,
                        3.2562428543216423e-9,
                        -1.892919739276469e-7,
                    ],
                    [
                        2.5048162288448987e-5,
                        8.251426900729693e-6,
                        2.400712823650017e-6,
                        -1.4137151747287231e-7,
                        -1.892919739276469e-7,
                        1.1116401573135049e-5,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "8467 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3450677500966914, epsilon = tol);
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
                semi_major_axis: 2.1907001009004063,
                eccentricity_sin_lon: -0.13400272182063963,
                eccentricity_cos_lon: 0.15338009135380934,
                tan_half_incl_sin_node: 0.0029985238451201465,
                tan_half_incl_cos_node: -0.05948593480260206,
                mean_longitude: 4.224751249408142,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.1404880736386688e-5,
                eccentricity_sin_lon: 1.3646865441167257e-5,
                eccentricity_cos_lon: 5.318472489593462e-6,
                tan_half_incl_sin_node: 3.449533103891244e-7,
                tan_half_incl_cos_node: 8.503524434876494e-7,
                mean_longitude: 2.664139085596166e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.581689193389379e-10,
                        -2.4444869402655374e-10,
                        7.20516937585694e-11,
                        -1.8829447695752316e-12,
                        -6.3276333857022975e-12,
                        4.3450961652176595e-10,
                    ],
                    [
                        -2.4444869402655374e-10,
                        1.862369363693252e-10,
                        -6.032989819801514e-11,
                        7.88198535010003e-15,
                        -6.599946297608992e-13,
                        -3.5827702309231364e-10,
                    ],
                    [
                        7.205169375856947e-11,
                        -6.032989819801514e-11,
                        2.8286149622562476e-11,
                        2.047044860400452e-14,
                        1.4233196526634721e-13,
                        1.2757883513045594e-10,
                    ],
                    [
                        -1.8829447695752336e-12,
                        7.881985350100153e-15,
                        2.047044860400451e-14,
                        1.1899278634841562e-13,
                        2.6430760135859734e-13,
                        3.758769572765907e-13,
                    ],
                    [
                        -6.327633385702298e-12,
                        -6.599946297608993e-13,
                        1.4233196526634721e-13,
                        2.6430760135859734e-13,
                        7.230992781454159e-13,
                        2.6058289247883292e-12,
                    ],
                    [
                        4.345096165217667e-10,
                        -3.582770230923137e-10,
                        1.2757883513045594e-10,
                        3.7587695727659063e-13,
                        2.605828924788329e-12,
                        7.097637067401174e-10,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "33803 N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.40877453793803264, epsilon = tol);
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
                semi_major_axis: 1.8017539671066383,
                eccentricity_sin_lon: 0.2693970247101432,
                eccentricity_cos_lon: 0.08869492421101689,
                tan_half_incl_sin_node: 0.0008307916711964613,
                tan_half_incl_cos_node: 0.10167059716063825,
                mean_longitude: 1.69365955893391,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 8.91843297803619e-7,
                eccentricity_sin_lon: 1.4528052999504664e-6,
                eccentricity_cos_lon: 5.920319266940421e-6,
                tan_half_incl_sin_node: 2.792168578051512e-7,
                tan_half_incl_cos_node: 6.218091729543099e-7,
                mean_longitude: 9.89691163439512e-6,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        7.953844678372347e-13,
                        -1.0232494678734991e-13,
                        5.059210617753894e-12,
                        -1.5066926551194571e-13,
                        2.4607703529074084e-13,
                        -8.52263281112491e-12,
                    ],
                    [
                        -1.0232494678735079e-13,
                        2.110643239564165e-12,
                        -3.0905811444958427e-12,
                        2.483716495651774e-13,
                        -7.926397586792164e-13,
                        4.808599014203309e-12,
                    ],
                    [
                        5.0592106177538936e-12,
                        -3.0905811444958334e-12,
                        3.505018022250596e-11,
                        -1.2295382103666956e-12,
                        2.4527485767320916e-12,
                        -5.855948314983968e-11,
                    ],
                    [
                        -1.5066926551194556e-13,
                        2.483716495651779e-13,
                        -1.2295382103666956e-12,
                        7.796205368258202e-14,
                        -1.4540027906052378e-13,
                        2.0297758780260734e-12,
                    ],
                    [
                        2.460770352907407e-13,
                        -7.926397586792161e-13,
                        2.452748576732091e-12,
                        -1.454002790605238e-13,
                        3.866466475701229e-13,
                        -3.988850346505511e-12,
                    ],
                    [
                        -8.522632811124905e-12,
                        4.808599014203313e-12,
                        -5.855948314983968e-11,
                        2.0297758780260734e-12,
                        -3.988850346505511e-12,
                        9.794885989902548e-11,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "K09R05F N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.6128100691261191, epsilon = tol);
    }
}
