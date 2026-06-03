mod common;

use outfit::{
    orbit_type::uncertainty::{EquinoctialUncertainty, KeplerianUncertainty, OrbitalCovariance},
    EquinoctialElements, KeplerianElements, OrbitalElements,
};

use crate::common::approx_equal;

#[test]
fn test_uncertainty_propagation() {
    let equinoctial_orbit = OrbitalElements::Equinoctial {
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

    let keplerian_orbit = equinoctial_orbit.to_keplerian().unwrap();

    let expected_keplerian_propag = OrbitalElements::Keplerian {
        elements: KeplerianElements {
            reference_epoch: 57049.2684537375,
            semi_major_axis: 1.8021517900042052,
            eccentricity: 0.2839820354128493,
            inclination: 0.20266238925780133,
            ascending_node_longitude: 0.008826172835575467,
            periapsis_argument: 1.2411480851756391,
            mean_anomaly: 0.4421910841246559,
        },
        uncertainty: Some(KeplerianUncertainty {
            semi_major_axis: 1.910876358918557e-6,
            eccentricity: 3.926080684435881e-6,
            inclination: 2.2639852329024065e-6,
            ascending_node_longitude: 6.113264876575711e-6,
            periapsis_argument: 4.049775340683106e-5,
            mean_anomaly: 2.2182426229638676e-5,
        }),
        covariance: Some(OrbitalCovariance {
            matrix: [
                [
                    3.651448459073842e-12,
                    6.857127156611333e-12,
                    1.6782354228854548e-12,
                    -3.781001511911568e-12,
                    -7.433110873463038e-11,
                    3.899825789832625e-11,
                ],
                [
                    6.857127156611329e-12,
                    1.5414109540700513e-11,
                    2.690953229794561e-15,
                    -2.0474618140821963e-12,
                    -1.2349406349235225e-10,
                    5.97243215927523e-11,
                ],
                [
                    1.6782354228854548e-12,
                    2.6909532297930087e-15,
                    5.1256291348001634e-12,
                    -9.989144038881854e-12,
                    -5.3024087432235095e-11,
                    3.518354634255312e-11,
                ],
                [
                    -3.781001511911568e-12,
                    -2.047461814082196e-12,
                    -9.989144038881855e-12,
                    3.7372007451174244e-11,
                    8.98813435388229e-11,
                    -6.947495524468516e-11,
                ],
                [
                    -7.433110873463033e-11,
                    -1.2349406349235207e-10,
                    -5.302408743223507e-11,
                    8.988134353882289e-11,
                    1.6400680310004965e-9,
                    -8.833005679743845e-10,
                ],
                [
                    3.8998257898326207e-11,
                    5.972432159275218e-11,
                    3.5183546342553095e-11,
                    -6.947495524468513e-11,
                    -8.833005679743845e-10,
                    4.920600334333619e-10,
                ],
            ]
            .into(),
        }),
    };

    approx_equal(&expected_keplerian_propag, &keplerian_orbit, 1e-10);
}
