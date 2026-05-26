//! # Equinoctial orbital elements
//!
//! This module defines the [`crate::orbit_type::equinoctial_element::EquinoctialElements`] struct and its associated methods,
//! providing a **non-singular representation of orbital elements**.  
//! Equinoctial elements are particularly well suited for orbit determination and propagation
//! because they remain regular for small eccentricities and inclinations, unlike classical Keplerian elements.
//!
//! ## What are equinoctial elements?
//!
//! The six equinoctial elements are:
//!
//! 1. **a** – Semi-major axis (AU)
//! 2. **h** = e · sin(ω + Ω)
//! 3. **k** = e · cos(ω + Ω)
//! 4. **p** = tan(i / 2) · sin Ω
//! 5. **q** = tan(i / 2) · cos Ω
//! 6. **λ** = Mean longitude = Ω + ω + M
//!
//! where:
//! - `e` is the eccentricity,
//! - `i` is the inclination,
//! - `Ω` is the longitude of ascending node,
//! - `ω` is the argument of periapsis,
//! - `M` is the mean anomaly.
//!
//! These elements are **well-behaved even for nearly circular and equatorial orbits**.
//!
//! ## Provided functionality
//!
//! - **Conversion** between [`crate::orbit_type::equinoctial_element::EquinoctialElements`] and [`crate::orbit_type::keplerian_element::KeplerianElements`].
//! - **Propagation** of two-body orbits:
//!   * Solve the generalized Kepler equation
//!   * Compute inertial Cartesian position and velocity from equinoctial elements
//! - **Variational derivatives**:
//!   * Compute the Jacobian matrices of position and velocity with respect to the six elements
//!   * Useful for orbit fitting and uncertainty propagation
//!
//! ## Key methods
//!
//! - [`crate::orbit_type::equinoctial_element::EquinoctialElements::solve_kepler_equation`]  
//!   Solve the generalized Kepler equation in equinoctial form.
//!
//! - [`crate::orbit_type::equinoctial_element::EquinoctialElements::solve_two_body_problem`]  
//!   Propagate an orbit from `t0` to `t1` using the two-body approximation,
//!   returning position, velocity, and optionally Jacobians.
//!
//! - [`crate::orbit_type::equinoctial_element::EquinoctialElements::compute_cartesian_position_and_velocity`]  
//!   Compute position/velocity vectors from the equinoctial elements at a given epoch.
//!
//! - [`crate::orbit_type::equinoctial_element::EquinoctialElements::compute_w_vector`]  
//!   Compute the **w** basis vector of the equinoctial frame (normal to the orbital plane).
//!
//! ## Units
//!
//! - Lengths: **AU**
//! - Velocities: **AU/day**
//! - Angles: **radians**
//! - Time: **days** (usually Modified Julian Date, TDB/TT)
//!
//! ## Advantages of equinoctial elements
//!
//! * Avoid singularities when `e → 0` or `i → 0`
//! * Smooth derivatives, ideal for gradient-based fitting
//! * Directly compatible with least-squares adjustment and orbit uncertainty analysis
//!
//! ## Example
//!
//! ```rust, no_run
//! use outfit::orbit_type::equinoctial_element::EquinoctialElements;
//! use outfit::orbit_type::keplerian_element::KeplerianElements;
//!
//! // Define equinoctial elements
//! let equ = EquinoctialElements {
//!     reference_epoch: 59000.0,
//!     semi_major_axis: 1.2,
//!     eccentricity_sin_lon: 0.02,
//!     eccentricity_cos_lon: 0.03,
//!     tan_half_incl_sin_node: 0.001,
//!     tan_half_incl_cos_node: 0.05,
//!     mean_longitude: 1.0,
//! };
//!
//! // Propagate 10 days ahead using the two-body model
//! let (pos, vel, _) = equ.solve_two_body_problem(59000.0, 59010.0, false).unwrap();
//!
//! // Convert to classical Keplerian elements
//! let kep: KeplerianElements = equ.into();
//! ```
//!
//! ## See also
//!
//! - [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements)
//! - Milani & Gronchi, *Theory of Orbit Determination* (2010).
use core::f64;
use std::{f64::consts::PI, fmt};

use nalgebra::{Matrix3x6, Matrix6x3, Vector3};
use roots::{find_root_newton_raphson, SimpleConvergency};

use crate::{
    constants::{DPI, GAUSS_GRAV_SQUARED},
    kepler::principal_angle,
    orbit_type::keplerian_element::KeplerianElements,
    outfit_errors::OutfitError,
};

/// Type alias for the result of two-body propagation from equinoctial elements.
///
/// The tuple contains:
/// -   `Vector3<f64>`: Cartesian position vector (AU)
/// -   `Vector3<f64>`: Cartesian velocity vector (AU/day)
/// -   `Option<(Matrix6x3<f64>, Matrix6x3<f64>)>`: Optional Jacobian matrices of position and velocity with respect to the six equinoctial elements.
pub type TwoBodyResult = (
    Vector3<f64>,
    Vector3<f64>,
    Option<(Matrix6x3<f64>, Matrix6x3<f64>)>,
);

#[derive(Debug, Clone, PartialEq)]
pub struct EquinoctialLimits {
    pub eccentricity_limit: f64,
    min_semi_major_axis: f64,
    max_semi_major_axis: f64,
    min_periapsis_distance: f64,
    max_apoapsis_distance: f64,
}

impl Default for EquinoctialLimits {
    fn default() -> Self {
        Self {
            eccentricity_limit: 1.2, // Allow for some hyperbolic orbits (e > 1) but not extreme ones
            min_semi_major_axis: 1e-6, // 1e-6 AU ~ 150 km
            max_semi_major_axis: 1e4, // 10,000 AU
            min_periapsis_distance: 1e-6, // 1e-6 AU ~ 150 km
            max_apoapsis_distance: 1e4, // 10,000 AU
        }
    }
}

/// Equinoctial orbital elements.
/// Units:
/// - a: AU (astronomical units)
/// - h, k: dimensionless (related to eccentricity)
/// - p, q: dimensionless (related to inclination)
/// - lambda: radians (mean longitude)
#[derive(Debug, PartialEq, Clone)]
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
    /// Compute the squared orbital eccentricity from the equinoctial parameters `h` and `k`.
    ///
    /// Details
    /// -------
    ///
    /// In the equinoctial element formulation:
    ///
    /// * `h = e * sin(Ω + ω)`
    /// * `k = e * cos(Ω + ω)`
    ///
    /// where `e` is the eccentricity, `Ω` the longitude of ascending node,
    /// and `ω` the argument of periapsis.
    ///
    /// The squared eccentricity is simply:
    ///
    /// ```text
    /// e² = h² + k²
    /// ```
    ///
    /// Returns
    /// -------
    /// * The value of `e²` (dimensionless).
    pub fn squared_eccentricity(&self) -> f64 {
        self.eccentricity_sin_lon.powi(2) + self.eccentricity_cos_lon.powi(2)
    }

    /// Compute the orbital eccentricity from the equinoctial parameters `h` and `k`.
    ///
    /// Details
    /// -------
    ///
    /// The eccentricity `e` is obtained from:
    ///
    /// ```text
    /// e = sqrt(h² + k²)
    /// ```
    ///
    /// See also [`EquinoctialElements::squared_eccentricity`] for the squared value.
    ///
    /// Returns
    /// -------
    /// * The eccentricity `e` (dimensionless).
    ///
    /// See also
    /// --------
    /// * [`EquinoctialElements::squared_eccentricity`] – Returns `e²` directly.
    pub fn eccentricity(&self) -> f64 {
        self.squared_eccentricity().sqrt()
    }

    /// Check if the equinoctial elements represent a "bizarre" orbit based on specified limits.
    ///
    /// # Arguments
    ///
    /// - `limits` - An instance of `EquinoctialLimits` defining thresholds for eccentricity, semi-major axis, periapsis distance, and apoapsis distance.
    ///
    /// # Returns
    ///
    /// - `true` if the orbit is considered bizarre (e.g., too eccentric, too small/large semi-major axis, or extreme periapsis/apoapsis distances), otherwise `false`.
    pub fn is_bizarre(&self, limits: &EquinoctialLimits) -> bool {
        let e = self.eccentricity();
        let periapsis = self.semi_major_axis * (1. - e); // periapsis distance
        let apoapsis = self.semi_major_axis * (1. + e); // apoapsis distance

        e > limits.eccentricity_limit
            || self.semi_major_axis < limits.min_semi_major_axis
            || self.semi_major_axis > limits.max_semi_major_axis
            || periapsis < limits.min_periapsis_distance
            || apoapsis > limits.max_apoapsis_distance
    }

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
    pub fn solve_kepler_equation(
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

    /// Compute the **w** basis vector of the equinoctial reference frame.
    ///
    /// In the equinoctial orbital elements formalism, the three orthonormal
    /// basis vectors (**f**, **g**, **w**) define the orientation of the orbital
    /// plane in inertial space:
    ///
    /// * **f**: points roughly toward the pericenter,
    /// * **g**: 90° ahead of f in the plane,
    /// * **w**: normal to the orbital plane.
    ///
    /// This method computes the third vector **w**, which is orthogonal to both
    /// **f** and **g**, using the equinoctial parameters:
    /// - `tan_half_incl_sin_node` = tan(i/2) * sin(Ω)
    /// - `tan_half_incl_cos_node` = tan(i/2) * cos(Ω)
    ///
    /// The resulting vector ensures that {f, g, w} forms a right-handed
    /// orthonormal basis.
    ///
    /// Arguments
    /// ---------
    /// * `inv_u` – Reciprocal of `u`, where `u = 1 + p² + q²`, used to normalize
    ///   the equinoctial basis vectors.
    ///
    /// Returns
    /// --------
    /// * `Vector3<f64>` – The **w** vector of the equinoctial frame, normalized.
    pub fn compute_w_vector(&self, inv_u: f64) -> Vector3<f64> {
        Vector3::new(
            2. * self.tan_half_incl_sin_node * inv_u,
            -2. * self.tan_half_incl_cos_node * inv_u,
            (1. - self.tan_half_incl_sin_node.powi(2) - self.tan_half_incl_cos_node.powi(2))
                * inv_u,
        )
    }

    /// Computes the Jacobian matrices of position and velocity
    /// with respect to the equinoctial orbital elements.
    ///
    /// This function calculates how small variations in the six equinoctial
    /// elements affect the Cartesian state vector (3D position and 3D velocity)
    /// of an orbiting body. It returns two 6×3 matrices:
    ///
    /// * The first matrix contains the partial derivatives of position
    ///   with respect to each orbital element.
    /// * The second matrix contains the partial derivatives of velocity
    ///   with respect to each orbital element.
    ///
    /// The equinoctial elements are, in order:
    /// 1. a – semi-major axis (AU)
    /// 2. h – e * sin(ω + Ω)
    /// 3. k – e * cos(ω + Ω)
    /// 4. p – tan(i/2) * sin Ω
    /// 5. q – tan(i/2) * cos Ω
    /// 6. λ – mean longitude at epoch
    ///
    /// These derivatives are commonly used for:
    /// * Orbit propagation using variational equations,
    /// * Covariance propagation and orbit uncertainty estimation,
    /// * Sensitivity analysis,
    /// * Non-linear least squares fitting of orbital parameters.
    ///
    /// # Arguments
    ///
    /// * `t0`, `t1` – Epoch of the elements and target time (days)
    /// * `mean_motion` – Mean motion (sqrt(GM / a^3), in rad/day)
    /// * `mean_longitude_t1` – Mean longitude at `t1`
    /// * `eccentric_anomaly` – Generalized eccentric anomaly F (radians)
    /// * `inv_u` – 1 / (1 + p² + q²), normalizing factor for the equinoctial frame
    /// * `beta` – Auxiliary factor: 1 / (1 + sqrt(1 - h² - k²))
    /// * `sin_ecc_anom`, `cos_ecc_anom` – sin(F) and cos(F)
    /// * `xe`, `ye` – Position components in the orbital plane (AU)
    /// * `v_xe`, `v_ye` – Velocity components in the orbital plane (AU/day)
    /// * `f_vector`, `g_vector` – Base vectors of the equinoctial reference frame
    /// * `cart_position`, `cart_velocity` – Cartesian position and velocity vectors
    ///
    /// # Returns
    ///
    /// Returns a tuple `(dxde, dvde)` where:
    ///
    /// * `dxde` is a 6×3 matrix of partial derivatives of position
    /// * `dvde` is a 6×3 matrix of partial derivatives of velocity
    ///
    /// # Units
    ///
    /// * Lengths: astronomical units (AU)
    /// * Time: days
    /// * Velocities: AU/day
    ///
    /// # See also
    /// * [`EquinoctialElements::compute_w_vector`] – builds the third equinoctial base vector (W)
    /// * [`KeplerianElements`] – conversion between Keplerian and equinoctial elements
    #[allow(clippy::too_many_arguments)]
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
        // Compute the third base vector of the equinoctial reference frame.
        // f_vector and g_vector define the orbital plane, w_vector is perpendicular.
        let w_vector = self.compute_w_vector(inv_u);

        // Compute useful scalars for derivative formulas
        let r = (xe.powi(2) + ye.powi(2)).sqrt(); // radial distance in the orbital plane
        let inv_r = 1.0 / r; // 1/r
        let inv_1_beta = 1.0 / (1. - beta); // 1/(1 - beta) auxiliary term

        // Precompute repeated terms to simplify derivative expressions
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

        // =========================================================================
        // 1. Partial derivatives of position with respect to equinoctial elements
        // =========================================================================

        // Derivative of position with respect to semi-major axis (a)
        // The expression includes a correction for time difference (t1 - t0).
        let dxde_pos1 =
            (cart_position - 3. * cart_velocity * (t1 - t0) / 2.) / self.semi_major_axis;

        // Derivative with respect to element h (second column)
        let dx1de2_pos = -self.semi_major_axis
            * (tmp1 * tmp2 + self.semi_major_axis * cos_ecc_anom * tmp4 * inv_r);
        let dx2de2_pos = self.semi_major_axis
            * (tmp1 * tmp3 - 1.0 + self.semi_major_axis * cos_ecc_anom * tmp5 * inv_r);
        let dxde_pos2 = dx1de2_pos * f_vector + dx2de2_pos * g_vector;

        // Derivative with respect to element k (third column)
        let dx1de3_pos = -self.semi_major_axis
            * (tmp1 * tmp3 + 1.0 - self.semi_major_axis * sin_ecc_anom * tmp4 * inv_r);
        let dx2de3_pos = self.semi_major_axis
            * (tmp1 * tmp6 - self.semi_major_axis * sin_ecc_anom * tmp5 * inv_r);
        let dxde_pos3 = dx1de3_pos * f_vector + dx2de3_pos * g_vector;

        // Derivative with respect to element p (related to inclination and node)
        let dxde_pos4 = 2.0
            * (self.tan_half_incl_cos_node * (ye * f_vector - xe * g_vector) - xe * w_vector)
            * inv_u;

        // Derivative with respect to element q (related to inclination and node)
        let dxde_pos5 = 2.0
            * (self.tan_half_incl_sin_node * (-ye * f_vector + xe * g_vector) + ye * w_vector)
            * inv_u;

        // Derivative with respect to mean longitude (lambda)
        let dxde_pos6 = cart_velocity / mean_motion;

        // =========================================================================
        // 2. Partial derivatives of velocity with respect to equinoctial elements
        // =========================================================================

        // Derivative of velocity with respect to semi-major axis
        let dxde_vel1 = -(cart_velocity
            - 3. * GAUSS_GRAV_SQUARED * cart_position * (t1 - t0) / r.powi(3))
            / (2. * self.semi_major_axis);

        // Derivatives of velocity with respect to h
        let dx4de2_vel = tmp12
            * (tmp7 * tmp2
                + self.semi_major_axis.powi(2) * tmp8 * tmp4 * inv_r.powi(2)
                + tmp10 * cos_ecc_anom);
        let dx5de2_vel = -tmp12
            * (tmp7 * tmp3 + self.semi_major_axis.powi(2) * tmp8 * tmp5 * inv_r.powi(2)
                - tmp10 * sin_ecc_anom);
        let dxde_vel2 = dx4de2_vel * f_vector + dx5de2_vel * g_vector;

        // Derivatives of velocity with respect to k
        let dx4de3_vel = tmp12
            * (tmp7 * tmp3 + self.semi_major_axis.powi(2) * tmp9 * tmp4 * inv_r.powi(2)
                - tmp11 * cos_ecc_anom);
        let dx5de3_vel = -tmp12
            * (tmp7 * tmp6
                + self.semi_major_axis.powi(2) * tmp9 * tmp5 * inv_r.powi(2)
                + tmp11 * sin_ecc_anom);
        let dxde_vel3 = dx4de3_vel * f_vector + dx5de3_vel * g_vector;

        // Derivative of velocity with respect to p
        let dxde_vel4 = 2.0
            * (self.tan_half_incl_cos_node * (v_ye * f_vector - v_xe * g_vector) - v_xe * w_vector)
            * inv_u;

        // Derivative of velocity with respect to q
        let dxde_vel5 = 2.0
            * (self.tan_half_incl_sin_node * (-v_ye * f_vector + v_xe * g_vector)
                + v_ye * w_vector)
            * inv_u;

        // Derivative of velocity with respect to mean longitude
        let dxde_vel6 = -mean_motion * self.semi_major_axis.powi(3) * cart_position * inv_r.powi(3);

        // =========================================================================
        // 3. Assemble the Jacobian matrices
        // =========================================================================
        // Position derivatives: 6x3 matrix (transposed from 3x6)
        let derivative_position_matrix = Matrix3x6::from_columns(&[
            dxde_pos1, dxde_pos2, dxde_pos3, dxde_pos4, dxde_pos5, dxde_pos6,
        ])
        .transpose();

        // Velocity derivatives: 6x3 matrix
        let derivative_velocity_matrix = Matrix3x6::from_columns(&[
            dxde_vel1, dxde_vel2, dxde_vel3, dxde_vel4, dxde_vel5, dxde_vel6,
        ])
        .transpose();

        // Return the Jacobian matrices (d(r)/dE, d(v)/dE)
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
    ///       Otherwise, this value is `None`.
    ///
    /// Units
    /// -----
    ///
    /// - Position: **astronomical units** (UA)  
    /// - Velocity: **UA/day**  
    /// - Time: **days**  
    /// - Angles: **radians**
    #[allow(clippy::too_many_arguments)]
    pub fn compute_cartesian_position_and_velocity(
        &self,
        mean_motion: f64,
        eccentric_anomaly: f64,
        eccentricity_pow2: f64,
        compute_derivatives: bool,
        t0: f64,
        t1: f64,
        mean_longitude_t1: f64,
    ) -> TwoBodyResult {
        // -------------------------------------------------------------------------
        // 1. Compute auxiliary parameters
        // -------------------------------------------------------------------------
        // beta factor: depends on eccentricity, used to express coordinates in equinoctial form
        let beta = 1. / (1. + (1. - eccentricity_pow2).sqrt());

        // Combination of beta, sin(ω+Ω) and cos(ω+Ω), reused many times
        let beta_ecc_term = beta * self.eccentricity_sin_lon * self.eccentricity_cos_lon;

        // Precompute sine and cosine of the generalized eccentric anomaly
        let sin_ecc_anom = eccentric_anomaly.sin();
        let cos_ecc_anom = eccentric_anomaly.cos();

        // -------------------------------------------------------------------------
        // 2. Compute position in the orbital plane (x_e, y_e)
        // -------------------------------------------------------------------------
        // These expressions correspond to the equinoctial formulation
        // of the solution of Kepler's equation.
        let xe = self.semi_major_axis
            * ((1. - beta * self.eccentricity_sin_lon.powi(2)) * cos_ecc_anom
                + beta_ecc_term * sin_ecc_anom
                - self.eccentricity_cos_lon);

        let ye = self.semi_major_axis
            * ((1. - beta * self.eccentricity_cos_lon.powi(2)) * sin_ecc_anom
                + beta_ecc_term * cos_ecc_anom
                - self.eccentricity_sin_lon);

        // -------------------------------------------------------------------------
        // 3. Build equinoctial reference frame (f, g)
        // -------------------------------------------------------------------------
        // u = normalization factor = 1 + p^2 + q^2
        let u = 1. + self.tan_half_incl_sin_node.powi(2) + self.tan_half_incl_cos_node.powi(2);
        let inv_u = 1.0 / u;

        // Term common to several components of the f and g vectors
        let common_component =
            2. * self.tan_half_incl_sin_node * self.tan_half_incl_cos_node * inv_u;

        // f_vector points approximately toward periapsis
        let f_vector = Vector3::new(
            (1. - self.tan_half_incl_sin_node.powi(2) + self.tan_half_incl_cos_node.powi(2))
                * inv_u,
            common_component,
            -2. * self.tan_half_incl_sin_node * inv_u,
        );

        // g_vector is orthogonal to f_vector within the orbital plane
        let g_vector = Vector3::new(
            common_component,
            (1. + self.tan_half_incl_sin_node.powi(2) - self.tan_half_incl_cos_node.powi(2))
                * inv_u,
            2. * self.tan_half_incl_cos_node * inv_u,
        );

        // -------------------------------------------------------------------------
        // 4. Transform plane coordinates into 3D inertial position
        // -------------------------------------------------------------------------
        let cartesian_position = xe * f_vector + ye * g_vector;

        // -------------------------------------------------------------------------
        // 5. Compute velocity
        // -------------------------------------------------------------------------
        // v_const is a scaling factor for velocity components in the plane
        let v_const = mean_motion * self.semi_major_axis.powi(2) / (xe.powi(2) + ye.powi(2)).sqrt();

        // Velocity components in the orbital plane
        let v_xe = v_const
            * (beta_ecc_term * cos_ecc_anom
                - (1. - beta * self.eccentricity_sin_lon.powi(2)) * sin_ecc_anom);
        let v_ye = v_const
            * ((1. - beta * self.eccentricity_cos_lon.powi(2)) * cos_ecc_anom
                - beta_ecc_term * sin_ecc_anom);

        // Project velocity components onto inertial frame
        let cartesian_velocity = v_xe * f_vector + v_ye * g_vector;

        // -------------------------------------------------------------------------
        // 6. Optionally compute partial derivatives (Jacobian matrices)
        // -------------------------------------------------------------------------
        if compute_derivatives {
            let (dxde_pos, dxde_vel) = self.compute_derivative(
                t0,
                t1,
                mean_motion,
                mean_longitude_t1,
                eccentric_anomaly,
                inv_u,
                beta,
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

        // If derivatives are not requested, return position and velocity only
        (cartesian_position, cartesian_velocity, None)
    }

    /// Propagates an orbit (two-body problem) from equinoctial elements to a future state vector.
    ///
    /// This function computes the inertial Cartesian position and velocity of a body
    /// at a target epoch `t1`, starting from its equinoctial orbital elements defined at `t0`.
    ///
    /// Internally, the algorithm:
    /// 1. Computes the mean motion from the semi-major axis.
    /// 2. Advances the mean longitude to `t1`.
    /// 3. Solves the generalized Kepler equation to obtain the generalized eccentric anomaly `F`.
    /// 4. Converts the result from the orbital plane to 3D inertial coordinates using
    ///    the equinoctial basis vectors (`f`, `g`, `w`).
    ///
    /// If requested, it can also compute the Jacobian matrices of the state vector
    /// with respect to the six equinoctial elements.
    ///
    /// # Arguments
    ///
    /// * `t0` – Epoch of the orbital elements, in days.
    /// * `t1` – Target time for propagation, in days.
    /// * `compute_derivatives` – If `true`, also computes Jacobians
    ///   of position and velocity with respect to the elements.
    ///
    /// # Returns
    ///
    /// On success, returns a [`TwoBodyResult`] tuple:
    ///
    /// * `position`: 3D position in AU
    /// * `velocity`: 3D velocity in AU/day
    /// * `Option<(dxde, dvde)>`: Optional Jacobians if requested
    ///
    /// # Notes
    ///
    /// * Uses non-singular equinoctial elements, robust for nearly circular/equatorial orbits.
    /// * The generalized Kepler equation solved is:
    ///   `F - k·sin(F) + h·cos(F) = λ(t1)`
    /// * The solution domain is constrained to `[ω, ω + 2π]` for stability,
    ///   where ω is the longitude of periapsis.
    ///
    /// # Units
    /// * Lengths: AU
    /// * Velocities: AU/day
    /// * Angles: radians
    /// * Time: days
    ///
    /// # See also
    /// * [`EquinoctialElements::compute_cartesian_position_and_velocity`] – projects equinoctial elements to 3D
    /// * [`EquinoctialElements::solve_kepler_equation`] – solves the generalized Kepler equation
    /// * [`principal_angle`] – normalizes an angle to [0, 2π)
    pub fn solve_two_body_problem(
        &self,
        t0: f64,
        t1: f64,
        compute_derivatives: bool,
    ) -> Result<TwoBodyResult, OutfitError> {
        // ---------------------------------------------------------------------
        // 1. Compute mean motion n = sqrt(mu / a^3)
        // ---------------------------------------------------------------------
        let mean_motion = (GAUSS_GRAV_SQUARED / self.semi_major_axis.powi(3)).sqrt();

        // Advance the mean longitude to the target time:
        // λ(t1) = λ0 + n * (t1 - t0)
        let mut mean_longitude_t1 = self.mean_longitude + mean_motion * (t1 - t0);

        // ---------------------------------------------------------------------
        // 2. Compute squared eccentricity = h² + k²
        //    Used to detect circular orbits (avoid division by zero)
        // ---------------------------------------------------------------------
        let eccentricity_pow2 = self.squared_eccentricity();
        let epsilon = f64::EPSILON * 1e2; // threshold for near-circular orbit

        // ---------------------------------------------------------------------
        // 3. Compute longitude of periapsis (ω = atan2(h, k)) if e > 0
        // ---------------------------------------------------------------------
        let mut longitude_of_periastre = 0.0;
        if eccentricity_pow2 > epsilon {
            longitude_of_periastre =
                principal_angle(self.eccentricity_sin_lon.atan2(self.eccentricity_cos_lon));
        }

        // Normalize λ(t1) to [0, 2π)
        mean_longitude_t1 = principal_angle(mean_longitude_t1);

        // Ensure λ(t1) lies ahead of ω, so that F is in [ω, ω+2π]
        if mean_longitude_t1 < longitude_of_periastre {
            mean_longitude_t1 += DPI;
        }

        // ---------------------------------------------------------------------
        // 4. Solve the generalized Kepler equation for F
        //    F - k sin F + h cos F = λ(t1)
        // ---------------------------------------------------------------------
        let eccentric_anomaly =
            self.solve_kepler_equation(mean_longitude_t1, longitude_of_periastre)?;

        // ---------------------------------------------------------------------
        // 5. Compute final Cartesian position and velocity
        // ---------------------------------------------------------------------
        Ok(self.compute_cartesian_position_and_velocity(
            mean_motion,
            eccentric_anomaly,
            eccentricity_pow2,
            compute_derivatives,
            t0,
            t1,
            mean_longitude_t1,
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

impl fmt::Display for EquinoctialElements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        writeln!(
            f,
            "Equinoctial Elements @ epoch (MJD): {:.6}",
            self.reference_epoch
        )?;
        writeln!(f, "------------------------------------------------")?;
        writeln!(
            f,
            "  a   (semi-major axis)            = {:.6} AU",
            self.semi_major_axis
        )?;
        writeln!(
            f,
            "  h   (e·sin(Ω+ω))                 = {:.6}",
            self.eccentricity_sin_lon
        )?;
        writeln!(
            f,
            "  k   (e·cos(Ω+ω))                 = {:.6}",
            self.eccentricity_cos_lon
        )?;
        writeln!(
            f,
            "  p   (tan(i/2)·sinΩ)              = {:.6}",
            self.tan_half_incl_sin_node
        )?;
        writeln!(
            f,
            "  q   (tan(i/2)·cosΩ)              = {:.6}",
            self.tan_half_incl_cos_node
        )?;
        writeln!(
            f,
            "  λ   (mean longitude)             = {:.6} rad ({:.6}°)",
            self.mean_longitude,
            self.mean_longitude * rad_to_deg
        )
    }
}

#[cfg(test)]
mod test_equinoctial_element {
    use super::*;

    #[test]
    fn test_is_bizarre() {
        let equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713,
            eccentricity_sin_lon: 0.2693736809404963,
            eccentricity_cos_lon: 0.08856415260522467,
            tan_half_incl_sin_node: 0.0008089970142830734,
            tan_half_incl_cos_node: 0.10168201110394352,
            mean_longitude: 1.693697008,
        };
        let limits = EquinoctialLimits::default();
        assert!(!equ.is_bizarre(&limits));

        let bad_equ = EquinoctialElements {
            reference_epoch: 0.0,
            semi_major_axis: 0.000001, // too small
            eccentricity_sin_lon: 0.2693736809404963,
            eccentricity_cos_lon: 0.08856415260522467,
            tan_half_incl_sin_node: 0.0008089970142830734,
            tan_half_incl_cos_node: 0.10168201110394352,
            mean_longitude: 1.693697008 + 2. * std::f64::consts::PI, // out of range
        };
        assert!(bad_equ.is_bizarre(&limits));
    }

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
            eccentricity_sin_lon: 0.269_373_680_909_227_2,
            eccentricity_cos_lon: 8.856_415_260_013_56E-2,
            tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
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
            eccentricity_sin_lon: 0.269_373_680_909_227_2,
            eccentricity_cos_lon: 8.856_415_260_013_56E-2,
            tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
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
            eccentricity_sin_lon: 0.269_373_680_909_227_2,
            eccentricity_cos_lon: 8.856_415_260_013_56E-2,
            tan_half_incl_sin_node: 8.089_970_166_396_302E-4,
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

    #[cfg(test)]
    mod equinoctial_kepler_proptests {
        use super::*;
        use proptest::prelude::*;

        /// Small helper: wrap angle to [0, 2π).
        fn wrap_0_2pi(x: f64) -> f64 {
            let twopi = std::f64::consts::TAU;
            let y = x.rem_euclid(twopi);
            if y < 0.0 {
                y + twopi
            } else {
                y
            }
        }

        /// Reference Newton solver for the *classical* Kepler equation:
        /// E - e*sin(E) = M, with M ∈ (-π, π].
        ///
        /// Arguments
        /// -----------------
        /// * `e`: eccentricity in [0, 1).
        /// * `m`: mean anomaly wrapped to (-π, π].
        ///
        /// Return
        /// ----------
        /// * Eccentric anomaly E (radians).
        fn solve_kepler_classic_newton(e: f64, m: f64) -> f64 {
            // Initial guess: good for e <~ 0.9
            let (sm, cm) = m.sin_cos();
            let mut e_curr = m + e * sm * (1.0 + e * cm);

            // Two or three Newton iterations are enough in practice
            for _ in 0..5 {
                let (se, ce) = e_curr.sin_cos();
                let f = e_curr - e * se - m;
                let fp = 1.0 - e * ce;
                let delta = f / fp;
                e_curr -= delta;
                if delta.abs() < 1e-14 {
                    break;
                }
            }
            e_curr
        }

        /// Strategy to generate (e, ϖ, λ) in safe/realistic ranges.
        /// - e in [0, 0.9]
        /// - ϖ (longitude of periapsis) in (-π, π]
        /// - λ in [0, 2π), then optionally shifted so that λ >= ϖ,
        ///   mirroring the constraint used in `solve_two_body_problem`.
        fn strategy_e_w_lambda() -> impl Strategy<Value = (f64, f64, f64)> {
            let e_strat = prop_oneof![
                Just(0.0),
                Just(1e-12),
                0.0f64..=0.85f64,  // <-- borne à 0.85 au lieu de 0.9
                0.75f64..=0.85f64, // densifie près de la borne
            ];
            let w_strat = -std::f64::consts::PI..=std::f64::consts::PI;
            let lam_strat = 0.0f64..std::f64::consts::TAU;

            (e_strat, w_strat, lam_strat).prop_map(|(e, w, mut lam)| {
                lam = (lam).rem_euclid(std::f64::consts::TAU);
                if lam < (w).rem_euclid(std::f64::consts::TAU) {
                    lam += std::f64::consts::TAU;
                }
                (e, w, lam)
            })
        }

        /// Build an EquinoctialElements instance using only the relevant fields (h, k, λ).
        fn make_equinoctial_from_hk_lambda(h: f64, k: f64, lambda: f64) -> EquinoctialElements {
            EquinoctialElements {
                reference_epoch: 59000.0,
                semi_major_axis: 1.0, // not used by the solver
                eccentricity_sin_lon: h,
                eccentricity_cos_lon: k,
                tan_half_incl_sin_node: 0.0, // irrelevant here
                tan_half_incl_cos_node: 0.0, // irrelevant here
                mean_longitude: lambda,
            }
        }

        // Property 1: The returned solution F satisfies the equinoctial Kepler equation
        // within tight numerical tolerance: |F - k sin F + h cos F - λ| < 1e-12.
        proptest! {
            #[test]
            fn solution_satisfies_equation((e, w, lambda) in strategy_e_w_lambda()) {
                let h = e * w.sin();
                let k = e * w.cos();

                let equ = make_equinoctial_from_hk_lambda(h, k, lambda);
                let res = equ.solve_kepler_equation(lambda, w);

                prop_assume!(res.is_ok());
                let f_sol = res.unwrap();

                let residual = f_sol - k * f_sol.sin() + h * f_sol.cos() - lambda;
                prop_assert!(residual.abs() < 1e-12, "residual = {}", residual);
            }
        }

        // Property 2: In the quasi-circular limit (e ≈ 0), the solution is ≈ λ.
        proptest! {
            #[test]
            fn near_circular_gives_f_approx_lambda(
                w in -std::f64::consts::PI..=std::f64::consts::PI,
                lambda_base in 0.0f64..std::f64::consts::TAU
            ) {
                // e extremely small
                let e = 1e-14;
                let lambda = {
                    let mut l = wrap_0_2pi(lambda_base);
                    if l < wrap_0_2pi(w) { l += std::f64::consts::TAU; }
                    l
                };
                let h = e * w.sin();
                let k = e * w.cos();

                let equ = make_equinoctial_from_hk_lambda(h, k, lambda);
                let f_sol = equ.solve_kepler_equation(lambda, w).expect("NR should converge");

                // For e≈0, equation reduces to F ≈ λ
                let err = (f_sol - lambda).abs();
                prop_assert!(err < 1e-10, "F={}, λ={}, |F-λ|={}", f_sol, lambda, err);
            }
        }

        // Property 3: Consistency with classical Kepler:
        // Let M = λ - ϖ wrapped to (-π, π]; solving classical E - e sin E = M,
        // we must have F ≈ E + ϖ (mod 2π).
        proptest! {
            #[test]
            fn matches_classical_kepler_solution(
                (e, w, lambda) in strategy_e_w_lambda()
            ) {
                let h = e * w.sin();
                let k = e * w.cos();

                let equ = make_equinoctial_from_hk_lambda(h, k, lambda);
                let res = equ.solve_kepler_equation(lambda, w);

                // Si la prod ne converge pas, on rejette ce cas (pas d'échec de test)
                prop_assume!(res.is_ok());
                let f_sol = res.unwrap();

                // Classical mean anomaly
                let m = {
                    let pi = std::f64::consts::PI;
                    let y = (lambda - w + pi).rem_euclid(std::f64::consts::TAU) - pi;
                    if y <= -pi { y + std::f64::consts::TAU } else { y }
                };
                let e_classic = solve_kepler_classic_newton(e, m);

                // Compare F to (E + ϖ) modulo 2π avec min distance sur le cercle
                let got    = f_sol.rem_euclid(std::f64::consts::TAU);
                let target = (e_classic + w).rem_euclid(std::f64::consts::TAU);
                let diff = ((got - target + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU))
                         - std::f64::consts::PI;

                prop_assert!(diff.abs() < 1e-10, "F={}, target(E+ϖ)={}, diff={}", got, target, diff);
            }
        }
    }
}
