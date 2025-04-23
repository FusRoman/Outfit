use core::f64;
use std::f64::consts::PI;

use nalgebra::Vector3;
use roots::{find_root_newton_raphson, SimpleConvergency};

use crate::{
    constants::{DPI, GAUSS_GRAV_SQUARED},
    kepler::principal_angle,
    keplerian_element::KeplerianElements,
    outfit_errors::OutfitError,
};

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

impl EquinoctialElements {
    fn solve_kepler_equation(
        &self,
        mean_longitude_t1: f64,
        longitude_of_periastre: f64,
    ) -> Result<f64, OutfitError> {
        // Fonction R(F) = F - k·sin(F) + h·cos(F) - lambda
        let f = |fval: f64| -> f64 {
            fval - self.eccentricity_cos_lon * fval.sin() + self.eccentricity_sin_lon * fval.cos()
                - mean_longitude_t1
        };

        // Dérivée R'(F)
        let df = |fval: f64| -> f64 {
            1.0 - self.eccentricity_cos_lon * fval.cos() - self.eccentricity_sin_lon * fval.sin()
        };

        // Point de départ initial (comme en Fortran : pig + pol)
        let x0 = PI + longitude_of_periastre;

        // Critère de convergence similaire à eps=1e-14
        let mut tol = SimpleConvergency {
            eps: f64::EPSILON * 1e2, // ~2e-14
            max_iter: 25,
        };

        // Appel à la méthode de Newton-Raphson
        Ok(find_root_newton_raphson(x0, &f, &df, &mut tol)?)
    }

    fn compute_cartesian_position_and_velocity(
        &self,
        mean_motion: f64,
        eccentric_anomaly: f64,
        eccentricity_pow2: f64,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let beta = 1. / (1. + (1. - eccentricity_pow2).sqrt());

        let beta_ecc_term = beta * self.eccentricity_sin_lon * self.eccentricity_cos_lon;

        let sin_ecc_anom = eccentric_anomaly.sin();
        let cos_ecc_anom = eccentric_anomaly.cos();

        let xe = self.semi_major_axis
            * ((1. - beta * self.eccentricity_sin_lon.powi(2)) * cos_ecc_anom
                + beta_ecc_term * sin_ecc_anom
                - self.eccentricity_cos_lon);

        let ye = self.semi_major_axis
            * ((1. - beta * self.eccentricity_cos_lon.powi(2)) * sin_ecc_anom
                + beta_ecc_term * cos_ecc_anom
                - self.eccentricity_sin_lon);

        let u = 1. + self.tan_half_incl_sin_node.powi(2) + self.tan_half_incl_cos_node.powi(2);
        let inv_u = 1.0 / u;

        let common_component =
            2. * self.tan_half_incl_sin_node * self.tan_half_incl_cos_node * inv_u;

        let f_vector = Vector3::new(
            (1. - self.tan_half_incl_sin_node.powi(2) + self.tan_half_incl_cos_node.powi(2))
                * inv_u,
            common_component,
            -2. * self.tan_half_incl_sin_node * inv_u,
        );

        let g_vector = Vector3::new(
            common_component,
            (1. + self.tan_half_incl_sin_node.powi(2) - self.tan_half_incl_cos_node.powi(2))
                * inv_u,
            2. * self.tan_half_incl_cos_node * inv_u,
        );

        let cartesian_position = xe * f_vector + ye * g_vector;

        let v_const = mean_motion * self.semi_major_axis.powi(2) / (xe.powi(2) + ye.powi(2)).sqrt();

        let v_xe = v_const
            * (beta_ecc_term * cos_ecc_anom
                - (1. - beta * self.eccentricity_sin_lon.powi(2)) * sin_ecc_anom);
        let v_ye = v_const
            * ((1. - beta * self.eccentricity_cos_lon.powi(2)) * cos_ecc_anom
                - beta_ecc_term * sin_ecc_anom);
        let cartesian_velocity = v_xe * f_vector + v_ye * g_vector;

        (cartesian_position, cartesian_velocity)
    }

    pub(crate) fn solve_two_body_problem(
        &self,
        t0: f64,
        t1: f64,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        let mean_motion = (GAUSS_GRAV_SQUARED / self.semi_major_axis.powi(3)).sqrt();
        let mut mean_longitude_t1 = self.mean_longitude + mean_motion * (t1 - t0);

        let eccentricity_pow2 =
            self.eccentricity_sin_lon.powi(2) + self.eccentricity_cos_lon.powi(2);
        let epsilon = f64::EPSILON * 1e2;

        let mut longitude_of_periastre = 0.0;
        if eccentricity_pow2 > epsilon {
            longitude_of_periastre =
                principal_angle(self.eccentricity_sin_lon.atan2(self.eccentricity_cos_lon));
        }

        mean_longitude_t1 = principal_angle(mean_longitude_t1);
        if mean_longitude_t1 < longitude_of_periastre {
            mean_longitude_t1 += DPI;
        }

        let eccentric_anomaly =
            self.solve_kepler_equation(mean_longitude_t1, longitude_of_periastre)?;

        Ok(self.compute_cartesian_position_and_velocity(
            mean_motion,
            eccentric_anomaly,
            eccentricity_pow2,
        ))
    }
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

    #[test]
    fn test_kepler_equation() {
        let equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.26937368090922720,
            eccentricity_cos_lon: 8.8564152600135601E-002,
            tan_half_incl_sin_node: 8.0899701663963020E-004,
            tan_half_incl_cos_node: 0.10168201109730375,
            mean_longitude: 1.6936970079414786,
        };

        let result = equ.solve_kepler_equation(1.8432075709935847, 1.2531511177826073);

        assert_eq!(
            result.unwrap(),
            2.0450042417470673,
            "Kepler equation solution does not match expected value"
        );
    }

    #[test]
    fn test_two_body_problem() {
        let equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.26937368090922720,
            eccentricity_cos_lon: 8.8564152600135601E-002,
            tan_half_incl_sin_node: 8.0899701663963020E-004,
            tan_half_incl_cos_node: 0.10168201109730375,
            mean_longitude: 1.6936970079414786,
        };

        let (pos, vel) = equ.solve_two_body_problem(0., 21.019733018845727).unwrap();
        assert_eq!(
            pos,
            Vector3::new(-0.9321264203108841, 1.0784562905421133, 0.22313456997634373)
        );
        assert_eq!(
            vel,
            Vector3::new(
                -0.013800441828595238,
                -0.007301622877053736,
                -0.001477839051396935
            )
        );
    }
}
