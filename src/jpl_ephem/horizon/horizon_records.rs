use nalgebra::Vector3;

use super::interpolation_result::InterpResult;

/// The HorizonRecord struct represents a record of Tchebycheff coefficients
/// to derive a celestial object's position over a specified time interval.
#[derive(Debug, PartialEq)]
pub struct HorizonRecord {
    pub start_jd: f64,
    pub end_jd: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl HorizonRecord {
    /// Create a new HorizonRecord instance.
    ///
    /// Arguments
    /// ---------
    /// * `start_jd`: The Julian date at the start of the time interval.
    /// * `end_jd`: The Julian date at the end of the time interval.
    /// * `coeffs`: A slice of Tchebycheff coefficients for the x, y, and z
    ///   coordinates.
    /// * `offset`: The offset in the coefficients array where the Tchebycheff
    ///   coefficients for this record start.
    /// * `n_subintervals`: The number of subintervals in this record.
    /// * `n_coeffs`: The number of Tchebycheff coefficients for each coordinate.
    ///
    /// Returns
    /// -------
    /// * A new HorizonRecord instance.
    pub fn new(
        start_jd: f64,
        end_jd: f64,
        coeffs: &[f64],
        offset: usize,
        n_subintervals: usize,
        n_coeffs: usize,
    ) -> Self {
        let mut x = Vec::with_capacity(n_coeffs);
        let mut y = Vec::with_capacity(n_coeffs);
        let mut z = Vec::with_capacity(n_coeffs);

        let base = offset - 3 + n_subintervals * n_coeffs * 3;

        for j in 0..n_coeffs {
            let idx = base + j;

            let x_val = coeffs.get(idx).copied().expect("x coeff out of bounds");
            let y_val = coeffs
                .get(idx + n_coeffs)
                .copied()
                .expect("y coeff out of bounds");
            let z_val = coeffs
                .get(idx + n_coeffs * 2)
                .copied()
                .expect("z coeff out of bounds");

            x.push(x_val);
            y.push(y_val);
            z.push(z_val);
        }

        HorizonRecord {
            start_jd,
            end_jd,
            x,
            y,
            z,
        }
    }

    /// Compute the position, velocity, and acceleration of a celestial object
    /// at a given time (tau) using Tchebycheff coefficients.
    ///
    /// Arguments
    /// ---------
    /// * `tau`: A normalized time value between 0 and 1, representing the
    ///   position within the time interval defined by `start_jd` and `end_jd`.
    /// * `compute_velocity`: A boolean flag indicating whether to compute
    ///   the velocity of the object.
    /// * `compute_acceleration`: A boolean flag indicating whether to compute
    ///   the acceleration of the object.
    /// * `n_subintervals`: The number of subintervals to use for the
    ///   interpolation. This affects the scaling of the velocity and
    ///   acceleration calculations.
    ///
    /// Returns
    /// -------
    /// * An `InterpResult` struct containing the computed position,
    ///   velocity, and acceleration of the celestial object.
    pub fn interpolate(
        &self,
        tau: f64,
        compute_velocity: bool,
        compute_acceleration: bool,
        n_subintervals: usize,
    ) -> InterpResult {
        let dt1 = tau as i64;
        let temp = n_subintervals as f64 * tau;
        let tc = 2.0 * (temp.rem_euclid(1.0) + dt1 as f64) - 1.0;
        let mut twot = 0.;
        let n_coeff = self.x.len();

        let mut tcheb = vec![0.0; n_coeff];
        tcheb[0] = 1.0;

        if tc != tcheb[1] {
            tcheb[1] = tc;
            twot = tc + tc;
        }

        if n_coeff > 1 {
            for i in 2..n_coeff {
                tcheb[i] = twot * tcheb[i - 1] - tcheb[i - 2];
            }
        }

        let vfac = if compute_velocity || compute_acceleration {
            (2.0 * n_subintervals as f64) / (self.end_jd - self.start_jd)
        } else {
            0.0
        };

        let afac = vfac * vfac;

        let mut tcheb_deriv = vec![0.0; n_coeff];
        tcheb_deriv[1] = 1.0;
        tcheb_deriv[2] = twot + twot;

        let mut tcheb_accel = vec![0.0; n_coeff];

        if compute_velocity {
            if n_coeff > 1 {
                for i in 3..n_coeff {
                    tcheb_deriv[i] = twot * tcheb_deriv[i - 1] + 2.0 * tcheb[i - 1]
                        - tcheb_deriv[i - 2];
                }
            }
        }

        if compute_acceleration {
            tcheb_accel[0] = 0.0;
            if n_coeff > 1 {
                tcheb_accel[1] = 0.0;
                for i in 2..n_coeff {
                    tcheb_accel[i] = 2.0 * tcheb[1] * tcheb_accel[i - 1] + 4.0 * tcheb_deriv[i - 1]
                        - tcheb_accel[i - 2];
                }
            }
        }

        let eval = |coeffs: &Vec<f64>, basis: &Vec<f64>| -> f64 {
            coeffs.iter().zip(basis.iter()).map(|(c, b)| c * b).sum()
        };

        let x = eval(&self.x, &tcheb);

        let y = eval(&self.y, &tcheb);

        let z = eval(&self.z, &tcheb);

        let velocity = if compute_velocity {
            Some(Vector3::new(
                vfac * eval(&self.x, &tcheb_deriv),
                vfac * eval(&self.y, &tcheb_deriv),
                vfac * eval(&self.z, &tcheb_deriv),
            ))
        } else {
            None
        };

        let acceleration = if compute_acceleration {
            Some(Vector3::new(
                afac * eval(&self.x, &tcheb_accel),
                afac * eval(&self.y, &tcheb_accel),
                afac * eval(&self.z, &tcheb_accel),
            ))
        } else {
            None
        };

        InterpResult {
            position: Vector3::new(x, y, z),
            velocity,
            acceleration,
        }
    }
}
