use crate::{kepler::principal_angle, keplerian_element::KeplerianElements};

/// Equinoctial orbital elements.
/// Units:
/// - a: AU (astronomical units)
/// - h, k: dimensionless (related to eccentricity)
/// - p, q: dimensionless (related to inclination)
/// - lambda: radians (mean longitude)
#[derive(Debug, PartialEq)]
pub struct EquinoctialElements {
    pub reference_epoch: f64,        // Reference epoch (MJD)
    pub semi_major_axis: f64,        // Semi-major axis (AU)
    pub eccentricity_sin_lon: f64,   // h = e * sin(Ω + ω)
    pub eccentricity_cos_lon: f64,   // k = e * cos(Ω + ω)
    pub tan_half_incl_sin_node: f64, // p = tan(i/2) * sin(Ω)
    pub tan_half_incl_cos_node: f64, // q = tan(i/2) * cos(Ω)
    pub mean_longitude: f64,         // λ = Ω + ω + M
}

impl From<EquinoctialElements> for KeplerianElements {
    fn from(equinoctial: EquinoctialElements) -> Self {
        let EquinoctialElements {
            reference_epoch,
            semi_major_axis,
            eccentricity_sin_lon,
            eccentricity_cos_lon,
            tan_half_incl_sin_node,
            tan_half_incl_cos_node,
            mean_longitude,
        } = equinoctial;

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

        KeplerianElements {
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

#[cfg(test)]
mod test_equinoctial_element {
    use super::*;

    #[test]
    fn test_equinoctial_conversion() {
        let equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713,
            eccentricity_sin_lon: 0.2693736809404963,
            eccentricity_cos_lon: 0.08856415260522467,
            tan_half_incl_sin_node: 0.0008089970142830734,
            tan_half_incl_cos_node: 0.10168201110394352,
            mean_longitude: 1.693697008,
        };
        let kepler: KeplerianElements = equ.into();

        assert_eq!(
            kepler,
            KeplerianElements {
                reference_epoch: 0.0,
                semi_major_axis: 1.8017360713,
                eccentricity: 0.2835591457,
                inclination: 0.20267383289999996,
                ascending_node_longitude: 0.007955979,
                periapsis_argument: 1.2451951388,
                mean_anomaly: 0.4405458902000001,
            }
        );
    }
}
