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

impl From<KeplerianElements> for EquinoctialElements {
    fn from(kepler: KeplerianElements) -> Self {
        let KeplerianElements {
            reference_epoch,
            semi_major_axis,
            eccentricity,
            inclination,
            ascending_node_longitude,
            periapsis_argument,
            mean_anomaly,
        } = kepler;

        let dig = ascending_node_longitude + periapsis_argument;
        let h = eccentricity * dig.sin();
        let k = eccentricity * dig.cos();

        let tan_half_i = (inclination / 2.0).tan();
        let p = tan_half_i * ascending_node_longitude.sin();
        let q = tan_half_i * ascending_node_longitude.cos();

        let mut lambda = dig + mean_anomaly;
        lambda = principal_angle(lambda);

        EquinoctialElements {
            reference_epoch,
            semi_major_axis,
            eccentricity_sin_lon: h,
            eccentricity_cos_lon: k,
            tan_half_incl_sin_node: p,
            tan_half_incl_cos_node: q,
            mean_longitude: lambda,
        }
    }
}

#[cfg(test)]
mod test_equinoctial_element {
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
