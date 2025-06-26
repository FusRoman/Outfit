use core::f64;
use std::f64::consts::PI;

use nalgebra::{Matrix3x6, Matrix6x3, Vector3};
use roots::{find_root_newton_raphson, SimpleConvergency};

use crate::{
    constants::{DPI, GAUSS_GRAV_SQUARED},
    kepler::principal_angle,
    keplerian_element::KeplerianElements,
    outfit_errors::OutfitError,
};

pub type TwoBodyResult = (
    Vector3<f64>,
    Vector3<f64>,
    Option<(Matrix6x3<f64>, Matrix6x3<f64>)>,
);

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
    /// Create a new instance of `EquinoctialElements` from Keplerian elements.
    ///
    /// Arguments
    /// ---------
    /// * `reference_epoch`: Reference epoch (MJD)
    /// * `semi_major_axis`: Semi-major axis (AU)
    /// * `eccentricity`: Eccentricity (dimensionless)
    /// * `inclination`: Inclination (radians)
    /// * `ascending_node_longitude`: Ascending node longitude (radians)
    /// * `periapsis_argument`: Argument of periapsis (radians)
    /// * `mean_anomaly`: Mean anomaly (radians)
    ///
    /// Return
    /// ------
    /// * A new instance of `EquinoctialElements`
    pub(crate) fn from_kepler_internal(
        reference_epoch: f64,
        semi_major_axis: f64,
        eccentricity: f64,
        inclination: f64,
        ascending_node_longitude: f64,
        periapsis_argument: f64,
        mean_anomaly: f64,
    ) -> Self {
        let dig = ascending_node_longitude + periapsis_argument;
        let h = eccentricity * dig.sin();
        let k = eccentricity * dig.cos();

        let tan_half_i = (inclination / 2.0).tan();
        let p = tan_half_i * ascending_node_longitude.sin();
        let q = tan_half_i * ascending_node_longitude.cos();

        let lambda = principal_angle(dig + mean_anomaly);

        Self {
            reference_epoch,
            semi_major_axis,
            eccentricity_sin_lon: h,
            eccentricity_cos_lon: k,
            tan_half_incl_sin_node: p,
            tan_half_incl_cos_node: q,
            mean_longitude: lambda,
        }
    }

    /// Find the eccentric anomaly from the mean longitude and the longitude of periapsis.
    /// Solve the general kepler equation using equinoctial elements and the Newton-Raphson method.
    ///
    /// Arguments
    /// ---------
    /// * `mean_longitude_t1`: The mean longitude at time t1.
    /// * `longitude_of_periastre`: The longitude of periapsis.
    ///
    /// Return
    /// ------
    /// * A `Result` containing the eccentric anomaly (radians) or an error if the solution fails.
    fn solve_kepler_equation(
        &self,
        mean_longitude_t1: f64,
        longitude_of_periastre: f64,
    ) -> Result<f64, OutfitError> {
        let f = |fval: f64| -> f64 {
            fval - self.eccentricity_cos_lon * fval.sin() + self.eccentricity_sin_lon * fval.cos()
                - mean_longitude_t1
        };

        let df = |fval: f64| -> f64 {
            1.0 - self.eccentricity_cos_lon * fval.cos() - self.eccentricity_sin_lon * fval.sin()
        };

        let x0 = PI + longitude_of_periastre;

        let mut tol = SimpleConvergency {
            eps: f64::EPSILON * 1e2, // ~2e-14
            max_iter: 25,
        };

        Ok(find_root_newton_raphson(x0, &f, &df, &mut tol)?)
    }

    fn compute_w_vector(&self, inv_u: f64) -> Vector3<f64> {
        Vector3::new(
            2. * self.tan_half_incl_sin_node * inv_u,
            -2. * self.tan_half_incl_cos_node * inv_u,
            (1. - self.tan_half_incl_sin_node.powi(2) - self.tan_half_incl_cos_node.powi(2))
                * inv_u,
        )
    }

    /// Computes partial derivatives of the position and velocity vectors
    /// with respect to equinoctial orbital elements.
    ///
    /// This function evaluates how the Cartesian state vector of an orbiting body
    /// (position and velocity in 3D space) changes when each of its orbital elements
    /// is perturbed. It returns the Jacobian matrices of partial derivatives:
    ///
    /// - The first matrix contains ∂(position)/∂(element)
    /// - The second matrix contains ∂(velocity)/∂(element)
    ///
    /// These derivatives are useful for orbit propagation, uncertainty estimation,
    /// sensitivity analysis, and orbital parameter fitting.
    ///
    /// Arguments
    /// ---------
    /// * `t0`, `t1`: Epoch and target times (in days)
    /// * `mean_motion`: √(GM / a³), in rad/day
    /// * `mean_longitude_t1`: λ(t1) = λ₀ + n·(t1 - t0)
    /// * `eccentric_anomaly`: generalized eccentric anomaly F
    /// * `inv_u`: 1 / (1 + p² + q²) — normalization of equinoctial frame
    /// * `beta`: helper factor = 1 / (1 + √(1 - h² - k²))
    /// * `sin_ecc_anom`, `cos_ecc_anom`: sin(F), cos(F)
    /// * `xe`, `ye`: position in orbital plane
    /// * `v_xe`, `v_ye`: velocity in orbital plane
    /// * `f_vector`, `g_vector`: equinoctial base vectors
    /// * `cart_position`, `cart_velocity`: 3D Cartesian state vectors
    ///
    /// Returns
    /// -------
    /// * `(dxde, dvde)`: Tuple of 6×3 matrices containing partial derivatives of
    ///   position and velocity (rows) with respect to each orbital element (columns)
    ///
    /// Units
    /// -----
    /// * All lengths: UA
    /// * Time: days
    /// * Velocity: UA/day
    fn compute_derivative(
        &self,
        t0: f64,
        t1: f64,
        mean_motion: f64,
        mean_longitude_t1: f64,
        eccentric_anomaly: f64,
        inv_u: f64,
        beta: f64,
        sin_ecc_anom: f64,
        cos_ecc_anom: f64,
        xe: f64,
        ye: f64,
        v_xe: f64,
        v_ye: f64,
        f_vector: Vector3<f64>,
        g_vector: Vector3<f64>,
        cart_position: Vector3<f64>,
        cart_velocity: Vector3<f64>,
    ) -> (Matrix6x3<f64>, Matrix6x3<f64>) {
        let w_vector = self.compute_w_vector(inv_u);

        let r = (xe.powi(2) + ye.powi(2)).sqrt();
        let inv_r = 1.0 / r;
        let inv_1_beta = 1.0 / (1. - beta);

        let tmp1 = mean_longitude_t1 - eccentric_anomaly;
        let tmp2 = beta + self.eccentricity_sin_lon.powi(2) * beta.powi(3) * inv_1_beta;
        let tmp3 =
            self.eccentricity_sin_lon * self.eccentricity_cos_lon * beta.powi(3) * inv_1_beta;
        let tmp4 = beta * self.eccentricity_sin_lon - sin_ecc_anom;
        let tmp5 = beta * self.eccentricity_cos_lon - cos_ecc_anom;
        let tmp6 = beta + self.eccentricity_cos_lon.powi(2) * beta.powi(3) * inv_1_beta;
        let tmp7 = 1. - r / self.semi_major_axis;
        let tmp8 = sin_ecc_anom - self.eccentricity_sin_lon;
        let tmp9 = cos_ecc_anom - self.eccentricity_cos_lon;
        let tmp10 = self.semi_major_axis * cos_ecc_anom * inv_r;
        let tmp11 = self.semi_major_axis * sin_ecc_anom * inv_r;
        let tmp12 = mean_motion * self.semi_major_axis.powi(2) * inv_r;

        // compute the derivatives of the position w.r.t to the equinoctial elements
        let dxde_pos1 =
            (cart_position - 3. * cart_velocity * (t1 - t0) / 2.) / self.semi_major_axis;

        let dx1de2_pos = -self.semi_major_axis
            * (tmp1 * tmp2 + self.semi_major_axis * cos_ecc_anom * tmp4 * inv_r);
        let dx2de2_pos = self.semi_major_axis
            * (tmp1 * tmp3 - 1.0 + self.semi_major_axis * cos_ecc_anom * tmp5 * inv_r);

        let dxde_pos2 = dx1de2_pos * f_vector + dx2de2_pos * g_vector;

        let dx1de3_pos = -self.semi_major_axis
            * (tmp1 * tmp3 + 1.0 - self.semi_major_axis * sin_ecc_anom * tmp4 * inv_r);
        let dx2de3_pos = self.semi_major_axis
            * (tmp1 * tmp6 - self.semi_major_axis * sin_ecc_anom * tmp5 * inv_r);
        let dxde_pos3 = dx1de3_pos * f_vector + dx2de3_pos * g_vector;
        let dxde_pos4 = 2.0
            * (self.tan_half_incl_cos_node * (ye * f_vector - xe * g_vector) - xe * w_vector)
            * inv_u;
        let dxde_pos5 = 2.0
            * (self.tan_half_incl_sin_node * (-ye * f_vector + xe * g_vector) + ye * w_vector)
            * inv_u;
        let dxde_pos6 = cart_velocity / mean_motion;

        // compute the derivatives of the velocity w.r.t to the equinoctial elements
        let dxde_vel1 = -(cart_velocity
            - 3. * GAUSS_GRAV_SQUARED * cart_position * (t1 - t0) / r.powi(3))
            / (2. * self.semi_major_axis);

        let dx4de2_vel = tmp12
            * (tmp7 * tmp2
                + self.semi_major_axis.powi(2) * tmp8 * tmp4 * inv_r.powi(2)
                + tmp10 * cos_ecc_anom);
        let dx5de2_vel = -tmp12
            * (tmp7 * tmp3 + self.semi_major_axis.powi(2) * tmp8 * tmp5 * inv_r.powi(2)
                - tmp10 * sin_ecc_anom);
        let dxde_vel2 = dx4de2_vel * f_vector + dx5de2_vel * g_vector;

        let dx4de3_vel = tmp12
            * (tmp7 * tmp3 + self.semi_major_axis.powi(2) * tmp9 * tmp4 * inv_r.powi(2)
                - tmp11 * cos_ecc_anom);
        let dx5de3_vel = -tmp12
            * (tmp7 * tmp6
                + self.semi_major_axis.powi(2) * tmp9 * tmp5 * inv_r.powi(2)
                + tmp11 * sin_ecc_anom);
        let dxde_vel3 = dx4de3_vel * f_vector + dx5de3_vel * g_vector;

        let dxde_vel4 = 2.0
            * (self.tan_half_incl_cos_node * (v_ye * f_vector - v_xe * g_vector) - v_xe * w_vector)
            * inv_u;
        let dxde_vel5 = 2.0
            * (self.tan_half_incl_sin_node * (-v_ye * f_vector + v_xe * g_vector)
                + v_ye * w_vector)
            * inv_u;
        let dxde_vel6 = -mean_motion * self.semi_major_axis.powi(3) * cart_position * inv_r.powi(3);

        let derivative_position_matrix = Matrix3x6::from_columns(&[
            dxde_pos1, dxde_pos2, dxde_pos3, dxde_pos4, dxde_pos5, dxde_pos6,
        ])
        .transpose();

        let derivative_velocity_matrix = Matrix3x6::from_columns(&[
            dxde_vel1, dxde_vel2, dxde_vel3, dxde_vel4, dxde_vel5, dxde_vel6,
        ])
        .transpose();

        (derivative_position_matrix, derivative_velocity_matrix)
    }

    /// Computes the Cartesian position and velocity vectors from equinoctial orbital elements.
    ///
    /// This function transforms the current equinoctial orbital state of the body into
    /// inertial Cartesian coordinates. It does so by solving for the position in the orbital plane
    /// using the generalized eccentric anomaly, and then projecting the result into the
    /// 3D inertial space using the equinoctial reference frame.
    ///
    /// Optionally, it can also compute the partial derivatives of the position and velocity vectors
    /// with respect to the six equinoctial orbital elements — useful for orbit fitting,
    /// sensitivity analysis, and uncertainty propagation.
    ///
    /// Arguments
    /// ---------
    ///
    /// * `mean_motion`:  
    ///   The mean motion (`n = √(GM / a³)`), expressed in **radians per day**.  
    ///   Represents the average angular rate of the body in its orbit.
    ///
    /// * `eccentric_anomaly`:  
    ///   Generalized eccentric anomaly \( F \), in **radians**.  
    ///   This is typically obtained by solving the time evolution of the mean longitude.
    ///
    /// * `eccentricity_pow2`:  
    ///   The squared scalar eccentricity (`e² = h² + k²`),  
    ///   used in computing intermediate factors such as the β parameter.
    ///
    /// * `compute_derivatives`:  
    ///   If `true`, the function also computes the Jacobian matrices of the Cartesian
    ///   state vector with respect to the equinoctial orbital elements.
    ///
    /// Returns
    /// -------
    ///
    /// A `TwoBodyResult` tuple with the following elements:
    ///
    /// - `position`: Cartesian position vector (in **astronomical units**, UA).  
    ///   Expressed in the inertial reference frame.
    ///
    /// - `velocity`: Cartesian velocity vector (in **UA/day**).  
    ///   Also expressed in the inertial frame.
    ///
    /// - `Option<(Matrix6x3<f64>, Matrix6x3<f64>)>`:  
    ///   If `compute_derivatives` is `true`, returns:
    ///     - `dxde`: 6×3 matrix of partial derivatives of position w.r.t. orbital elements  
    ///     - `dvde`: 6×3 matrix of partial derivatives of velocity w.r.t. orbital elements  
    ///   Otherwise, this value is `None`.
    ///
    /// Units
    /// -----
    ///
    /// - Position: **astronomical units** (UA)  
    /// - Velocity: **UA/day**  
    /// - Time: **days**  
    /// - Angles: **radians**
    fn compute_cartesian_position_and_velocity(
        &self,
        mean_motion: f64,
        eccentric_anomaly: f64,
        eccentricity_pow2: f64,
        compute_derivatives: bool,
    ) -> TwoBodyResult {
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

        if compute_derivatives {
            let (dxde_pos, dxde_vel) = self.compute_derivative(
                self.reference_epoch,
                self.reference_epoch,
                mean_motion,
                self.mean_longitude,
                eccentric_anomaly,
                beta,
                beta_ecc_term,
                sin_ecc_anom,
                cos_ecc_anom,
                xe,
                ye,
                v_xe,
                v_ye,
                f_vector,
                g_vector,
                cartesian_position,
                cartesian_velocity,
            );

            return (
                cartesian_position,
                cartesian_velocity,
                Some((dxde_pos, dxde_vel)),
            );
        }

        (cartesian_position, cartesian_velocity, None)
    }

    /// Solves the two-body orbital motion using equinoctial elements.
    ///
    /// This function computes the inertial Cartesian position and velocity of an orbiting body
    /// at a future time `t1`, starting from its equinoctial orbital elements at epoch `t0`.
    /// It uses a generalized form of Kepler's equation and a transformation to 3D Cartesian
    /// space through the equinoctial reference frame.
    ///
    /// Optionally, the function can also compute the Jacobian matrix of partial derivatives
    /// of the position and velocity vectors with respect to the six equinoctial orbital elements.
    ///
    /// Arguments
    /// ---------
    ///
    /// - `t0`:  
    ///   Epoch time, in days. This is the time at which the orbital elements are defined.
    ///
    /// - `t1`:  
    ///   Prediction time, in days. The function computes the state at this time.
    ///
    /// - `compute_derivatives`:  
    ///   If `true`, the function also returns the Jacobian matrices of the state vector
    ///   with respect to the orbital elements.
    ///
    /// Returns
    /// -------
    ///
    /// Returns a `Result<TwoBodyResult, OutfitError>` where:
    ///
    /// - `TwoBodyResult` is a tuple:
    ///   - `position`: 3D position vector `[x, y, z]` in astronomical units (UA)
    ///   - `velocity`: 3D velocity vector `[vx, vy, vz]` in UA/day
    ///   - `Option<(dxde, dvde)>`: If `compute_derivatives` is `true`, this contains:
    ///     - `dxde`: matrix of ∂(position)/∂(elements)
    ///     - `dvde`: matrix of ∂(velocity)/∂(elements)
    ///
    /// Units
    /// -----
    ///
    /// - Lengths: **astronomical units (UA)**  
    /// - Velocities: **UA/day**  
    /// - Time: **days**  
    /// - Angles: **radians**
    ///
    /// Notes
    /// -----
    ///
    /// - This method uses equinoctial orbital elements, which provide a non-singular
    ///   representation even for circular and equatorial orbits.
    /// - The generalized Kepler equation is solved numerically:
    ///
    ///   `F - k·sin(F) + h·cos(F) = λ(t1)`
    ///
    ///   where `F` is a generalized eccentric anomaly, and `λ(t1)` is the mean longitude
    ///   at time `t1`.
    ///
    /// - The solution is constrained to lie in the interval `[W, W + 2π]`
    ///   for numerical stability (`W` being the longitude of pericenter).
    /// - The position and velocity are projected from the orbital plane
    ///   using the equinoctial basis vectors `f`, `g`, and `w`.
    pub(crate) fn solve_two_body_problem(
        &self,
        t0: f64,
        t1: f64,
        compute_derivatives: bool,
    ) -> Result<TwoBodyResult, OutfitError> {
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
            compute_derivatives,
        ))
    }
}

impl From<EquinoctialElements> for KeplerianElements {
    fn from(equinoctial: EquinoctialElements) -> Self {
        KeplerianElements::from_equinoctial_internal(
            equinoctial.reference_epoch,
            equinoctial.semi_major_axis,
            equinoctial.eccentricity_sin_lon,
            equinoctial.eccentricity_cos_lon,
            equinoctial.tan_half_incl_sin_node,
            equinoctial.tan_half_incl_cos_node,
            equinoctial.mean_longitude,
        )
    }
}

impl From<&EquinoctialElements> for KeplerianElements {
    fn from(equinoctial: &EquinoctialElements) -> Self {
        KeplerianElements::from_equinoctial_internal(
            equinoctial.reference_epoch,
            equinoctial.semi_major_axis,
            equinoctial.eccentricity_sin_lon,
            equinoctial.eccentricity_cos_lon,
            equinoctial.tan_half_incl_sin_node,
            equinoctial.tan_half_incl_cos_node,
            equinoctial.mean_longitude,
        )
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

        let (pos, vel, _) = equ
            .solve_two_body_problem(0., 21.019733018845727, false)
            .unwrap();
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

    #[test]
    fn test_compute_derivative() {
        let equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713154256,
            eccentricity_sin_lon: 0.26937368090922720,
            eccentricity_cos_lon: 8.8564152600135601E-002,
            tan_half_incl_sin_node: 8.0899701663963020E-004,
            tan_half_incl_cos_node: 0.10168201109730375,
            mean_longitude: 1.6936970079414786,
        };

        let (dx_pos, dx_vel) = equ.compute_derivative(
            0.,
            21.019733018845727,
            0.00711286689122354,
            1.8432075709935847,
            2.0450042417470673,
            0.9897659332253373,
            0.5104763141856585,
            0.8896546935605525,
            -0.4566339083178991,
            -0.9323069355123041,
            1.10114506236264,
            -0.013799246261211583,
            -0.007451892523877908,
            Vector3::new(
                0.9999987044435599,
                0.00016283716950135768,
                -0.0016014353743016747,
            ),
            Vector3::new(
                0.00016283716950135768,
                0.979533162007115,
                0.2012827812119039,
            ),
            Vector3::new(-0.9321264203108841, 1.0784562905421133, 0.22313456997634373),
            Vector3::new(
                -0.013800441828595238,
                -0.007301622877053736,
                -0.001477839051396935,
            ),
        );

        assert_eq!(
            dx_pos,
            [
                [
                    -0.2758472919839214,
                    -0.5803614626760855,
                    -3.3051181917865815,
                    0.2246273101991508,
                    0.0017270780533123044,
                    -1.9402080820074667
                ],
                [
                    0.7263403095474552,
                    -2.2723053964839406,
                    -1.1670672177854213,
                    -0.18762099832127083,
                    -0.44020925155213336,
                    -1.0265372582837307
                ],
                [
                    0.1497057464344368,
                    -0.4659843688851336,
                    -0.23441565316351645,
                    1.8451739525659905,
                    2.1348385937023004,
                    -0.20776981686813492
                ]
            ]
            .into()
        );

        assert_eq!(
            dx_vel,
            [
                [
                    0.002222700614910293,
                    -0.005788282594204328,
                    0.018371322890135426,
                    -0.0014557385356716304,
                    -1.1693077165124217e-5,
                    0.012911021052381672
                ],
                [
                    0.0038856205602975087,
                    -0.015583165352767119,
                    -0.010403249849722409,
                    -0.0027777913132127417,
                    0.0029475300114507746,
                    -0.014937857749615903
                ],
                [
                    0.0007948174310456126,
                    -0.0031927019517180885,
                    -0.0021677860848341836,
                    0.027318414370803085,
                    -0.014453795161933127,
                    -0.003090669964614741
                ]
            ]
            .into()
        );
    }
}
