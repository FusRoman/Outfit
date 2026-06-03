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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60672.2443617134,
                semi_major_axis: 3.2064058028481552,
                eccentricity_sin_lon: 0.05300520970081429,
                eccentricity_cos_lon: -0.023197692700634636,
                tan_half_incl_sin_node: 0.002896813138792025,
                tan_half_incl_cos_node: 0.09181010554057693,
                mean_longitude: 0.6256995904459722,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 0.007572375820104381,
                eccentricity_sin_lon: 0.0024777464468933156,
                eccentricity_cos_lon: 0.0007445419051153811,
                tan_half_incl_sin_node: 4.2789628256661375e-5,
                tan_half_incl_cos_node: 5.7090614265788426e-5,
                mean_longitude: 0.003334899745150928,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        5.73408755609015e-5,
                        1.869803895909995e-5,
                        5.597518961032454e-6,
                        -3.23358524526529e-7,
                        -4.293184367946207e-7,
                        2.5017097632757542e-5,
                    ],
                    [
                        1.8698038959099974e-5,
                        6.139227455092451e-6,
                        1.8070844259588129e-6,
                        -1.0561768228833672e-7,
                        -1.409534091855587e-7,
                        8.251003158870094e-6,
                    ],
                    [
                        5.597518961032457e-6,
                        1.8070844259588129e-6,
                        5.543426484728413e-7,
                        -3.149123145262512e-8,
                        -4.1500376105573355e-8,
                        2.4017490220016504e-6,
                    ],
                    [
                        -3.233585245265291e-7,
                        -1.0561768228833674e-7,
                        -3.149123145262512e-8,
                        1.8309522863432738e-9,
                        2.4373900210699776e-9,
                        -1.4146340017585484e-7,
                    ],
                    [
                        -4.2931843679462117e-7,
                        -1.4095340918555872e-7,
                        -4.1500376105573355e-8,
                        2.4373900210699772e-9,
                        3.259338237245045e-9,
                        -1.894260958072586e-7,
                    ],
                    [
                        2.501709763275752e-5,
                        8.251003158870094e-6,
                        2.401749022001649e-6,
                        -1.414634001758548e-7,
                        -1.894260958072586e-7,
                        1.1121556310207724e-5,
                    ],
                ]
                .into(),
            }),
        };

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 60465.26777915681,
                semi_major_axis: 2.190348311458185,
                eccentricity_sin_lon: -0.13373910921857446,
                eccentricity_cos_lon: 0.15339157238172804,
                tan_half_incl_sin_node: 0.0029876412023201416,
                tan_half_incl_cos_node: -0.05950692044872062,
                mean_longitude: 4.224365041422834,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 2.1385808329040844e-5,
                eccentricity_sin_lon: 1.3645976741407893e-5,
                eccentricity_cos_lon: 5.318680335330679e-6,
                tan_half_incl_sin_node: 3.4498285885877785e-7,
                tan_half_incl_cos_node: 8.504424577287931e-7,
                mean_longitude: 2.6647348205193914e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        4.573527978864727e-10,
                        -2.441395550149477e-10,
                        7.195967928874385e-11,
                        -1.8832832793515527e-12,
                        -6.328517179823655e-12,
                        4.340505064567979e-10,
                    ],
                    [
                        -2.4413955501494734e-10,
                        1.862126812270452e-10,
                        -6.032972260112714e-11,
                        8.172697046412589e-15,
                        -6.596088825542575e-13,
                        -3.58336068263501e-10,
                    ],
                    [
                        7.195967928874368e-11,
                        -6.032972260112711e-11,
                        2.8288360509433262e-11,
                        2.035203222614985e-14,
                        1.4220205527899883e-13,
                        1.2761184787385386e-10,
                    ],
                    [
                        -1.8832832793515527e-12,
                        8.17269704641279e-15,
                        2.0352032226149824e-14,
                        1.1901317290637546e-13,
                        2.6436375372263e-13,
                        3.754248073793805e-13,
                    ],
                    [
                        -6.328517179823653e-12,
                        -6.596088825542574e-13,
                        1.4220205527899883e-13,
                        2.6436375372262996e-13,
                        7.2325237390779e-13,
                        2.605903338872906e-12,
                    ],
                    [
                        4.340505064567979e-10,
                        -3.5833606826350087e-10,
                        1.2761184787385386e-10,
                        3.754248073793805e-13,
                        2.6059033388729057e-12,
                        7.100811663688513e-10,
                    ],
                ]
                .into(),
            }),
        };

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

        let expected = OrbitalElements::Equinoctial {
            elements: EquinoctialElements {
                reference_epoch: 57049.2684537375,
                semi_major_axis: 1.8021517900042052,
                eccentricity_sin_lon: 0.2694922786015968,
                eccentricity_cos_lon: 0.08955282358108035,
                tan_half_incl_sin_node: 0.0008974287327937245,
                tan_half_incl_cos_node: 0.10167548786557225,
                mean_longitude: 1.6921653421358704,
            },
            uncertainty: Some(EquinoctialUncertainty {
                semi_major_axis: 1.910876358918557e-6,
                eccentricity_sin_lon: 2.7271919973585478e-6,
                eccentricity_cos_lon: 1.2559941333300101e-5,
                tan_half_incl_sin_node: 6.143234310625764e-7,
                tan_half_incl_cos_node: 1.1476173256703189e-6,
                mean_longitude: 2.1064465635865037e-5,
            }),
            covariance: Some(OrbitalCovariance {
                matrix: [
                    [
                        3.651448459073842e-12,
                        -4.87907485491453e-13,
                        2.321298362132558e-11,
                        -3.7695250201166625e-13,
                        8.511532638002078e-13,
                        -3.91138523482157e-11,
                    ],
                    [
                        -4.879074854914533e-13,
                        7.437576190456506e-12,
                        -1.1647669978804286e-11,
                        9.359797430147383e-13,
                        -2.8577594338429333e-12,
                        1.853502993770551e-11,
                    ],
                    [
                        2.3212983621325566e-11,
                        -1.164766997880434e-11,
                        1.577521262959403e-10,
                        -3.47676746499932e-12,
                        8.610023673871895e-12,
                        -2.644913915663376e-10,
                    ],
                    [
                        -3.7695250201166625e-13,
                        9.359797430147385e-13,
                        -3.4767674649993202e-12,
                        3.7739327795249603e-13,
                        -5.048815271306508e-13,
                        5.7505636344116006e-12,
                    ],
                    [
                        8.511532638002078e-13,
                        -2.857759433842935e-12,
                        8.610023673871898e-12,
                        -5.048815271306507e-13,
                        1.3170255261786945e-12,
                        -1.4110008489365913e-11,
                    ],
                    [
                        -3.911385234821569e-11,
                        1.8535029937705585e-11,
                        -2.6449139156633765e-10,
                        5.750563634411601e-12,
                        -1.4110008489365913e-11,
                        4.437117125245391e-10,
                    ],
                ]
                .into(),
            }),
        };

        assert!(
            approx_equal(&expected, orbit.orbital_elements(), tol),
            "K09R05F N-body orbital elements differ from oracle beyond tolerance {tol}"
        );
        assert_relative_eq!(orbit.orbit_quality(), 0.3608868439717083, epsilon = tol);
    }
}
