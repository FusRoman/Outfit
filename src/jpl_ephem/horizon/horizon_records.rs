use nalgebra::Vector3;

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

#[derive(Debug, PartialEq)]
pub struct InterpResult {
    pub position: Vector3<f64>,
    pub velocity: Option<Vector3<f64>>,
    pub acceleration: Option<Vector3<f64>>,
}

impl HorizonRecord {
    pub fn interpolate(
        &self,
        tau: f64,
        compute_velocity: bool,
        compute_acceleration: bool,
        n_subintervals: usize,
    ) -> InterpResult {
        let n_coeff = self.x.len();

        let mut tcheb = vec![0.0; n_coeff];
        tcheb[0] = 1.0;
        if n_coeff > 1 {
            tcheb[1] = 2.0 * tau - 1.0;
            for i in 2..n_coeff {
                tcheb[i] = 2.0 * tcheb[1] * tcheb[i - 1] - tcheb[i - 2];
            }
        }

        let vfac = if compute_velocity || compute_acceleration {
            (2.0 * n_subintervals as f64) / (self.end_jd - self.start_jd)
        } else {
            0.0
        };

        let afac = vfac * vfac;

        let mut tcheb_deriv = vec![0.0; n_coeff];
        let mut tcheb_accel = vec![0.0; n_coeff];

        if compute_velocity {
            tcheb_deriv[0] = 0.0;
            if n_coeff > 1 {
                tcheb_deriv[1] = 1.0;
                for i in 2..n_coeff {
                    tcheb_deriv[i] = 2.0 * tcheb[1] * tcheb_deriv[i - 1] + 2.0 * tcheb[i - 1]
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
