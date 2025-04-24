use crate::{equinoctial_element::EquinoctialElements, kepler::principal_angle};

/// Keplerian orbital elements
/// Units:
/// * `reference_epoch`: MJD (Modified Julian Date)
/// * `semi_major_axis`: AU (Astronomical Units)
/// * `eccentricity`: unitless
/// * `inclination`: radians
/// * `ascending_node_longitude`: radians
/// * `periapsis_argument`: radians
/// * `mean_anomaly`: radians
#[derive(Debug, PartialEq)]
pub struct KeplerianElements {
    pub reference_epoch: f64,
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub ascending_node_longitude: f64,
    pub periapsis_argument: f64,
    pub mean_anomaly: f64,
}

impl KeplerianElements {
    pub (crate) fn from_equinoctial_internal(
        reference_epoch: f64,
        semi_major_axis: f64,
        eccentricity_sin_lon: f64,
        eccentricity_cos_lon: f64,
        tan_half_incl_sin_node: f64,
        tan_half_incl_cos_node: f64,
        mean_longitude: f64,
    ) -> Self {
        let eps = 1.0e-12; // small value for near-circular/near-equatorial tests
        let a = semi_major_axis;
        let ecc = (eccentricity_sin_lon.powi(2) + eccentricity_cos_lon.powi(2)).sqrt();

        // Compute dig = ω + Ω (or set to 0 if eccentricity ≈ 0)
        let dig = if ecc < eps {
            0.0
        } else {
            eccentricity_sin_lon.atan2(eccentricity_cos_lon)
        };

        let tgi2 = (tan_half_incl_sin_node.powi(2) + tan_half_incl_cos_node.powi(2)).sqrt();

        // Compute Ω (ascending node longitude)
        let omega_node = if tgi2 < eps {
            0.0
        } else {
            tan_half_incl_sin_node.atan2(tan_half_incl_cos_node)
        };

        let inclination = 2.0 * tgi2.atan(); // radians

        // Compute angular elements
        let periapsis_arg = principal_angle(dig - omega_node);
        let mean_anomaly = principal_angle(mean_longitude - dig);

        Self {
            reference_epoch: reference_epoch,
            semi_major_axis: a,
            eccentricity: ecc,
            inclination: inclination,
            ascending_node_longitude: omega_node,
            periapsis_argument: periapsis_arg,
            mean_anomaly: mean_anomaly,
        }
    }
}

impl From<KeplerianElements> for EquinoctialElements {
    fn from(k: KeplerianElements) -> Self {
        EquinoctialElements::from_kepler_internal(
            k.reference_epoch,
            k.semi_major_axis,
            k.eccentricity,
            k.inclination,
            k.ascending_node_longitude,
            k.periapsis_argument,
            k.mean_anomaly,
        )
    }
}

impl From<&KeplerianElements> for EquinoctialElements {
    fn from(k: &KeplerianElements) -> Self {
        EquinoctialElements::from_kepler_internal(
            k.reference_epoch,
            k.semi_major_axis,
            k.eccentricity,
            k.inclination,
            k.ascending_node_longitude,
            k.periapsis_argument,
            k.mean_anomaly,
        )
    }
}

#[cfg(test)]
mod test_keplerian_element {
    use super::*;
    use crate::equinoctial_element::EquinoctialElements;

    #[test]
    fn test_keplerian_conversion() {
        let kepler = KeplerianElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713,
            eccentricity: 0.2835591457,
            inclination: 0.2026738329,
            ascending_node_longitude: 0.0079559790,
            periapsis_argument: 1.2451951388,
            mean_anomaly: 0.4405458902,
        };

        let equinoctial: EquinoctialElements = kepler.into();

        assert_eq!(
            equinoctial,
            EquinoctialElements {
                reference_epoch: 0.0,
                semi_major_axis: 1.8017360713,
                eccentricity_sin_lon: 0.2693736809404963,
                eccentricity_cos_lon: 0.08856415260522467,
                tan_half_incl_sin_node: 0.0008089970142830734,
                tan_half_incl_cos_node: 0.10168201110394352,
                mean_longitude: 1.693697008,
            }
        );
    }
}
