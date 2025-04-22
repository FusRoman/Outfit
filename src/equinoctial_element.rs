use crate::{kepler::principal_angle, keplerian_element::KeplerianElements};

/// Equinoctial orbital elements.
/// Units:
/// - a: AU (astronomical units)
/// - h, k: dimensionless (related to eccentricity)
/// - p, q: dimensionless (related to inclination)
/// - lambda: radians (mean longitude)
#[derive(Debug)]
pub struct EquinoctialElements {
    pub reference_epoch: f64,        // Reference epoch (MJD)
    pub semi_major_axis: f64,        // Semi-major axis (AU)
    pub eccentricity_sin_lon: f64,   // h = e * sin(Ω + ω)
    pub eccentricity_cos_lon: f64,   // k = e * cos(Ω + ω)
    pub tan_half_incl_sin_node: f64, // p = tan(i/2) * sin(Ω)
    pub tan_half_incl_cos_node: f64, // q = tan(i/2) * cos(Ω)
    pub mean_longitude: f64,         // λ = Ω + ω + M
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
    use crate::keplerian_element::KeplerianElements;

    #[test]
    fn test_equinoctial_conversion() {
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
        dbg!(&equinoctial);
    }
}