use nalgebra::Vector3;

use super::interpolation_result::InterpResult;

/// One interpolation segment from a Horizons DATA RECORD.
///
/// A `HorizonRecord` stores the Chebyshev coefficients needed to compute
/// the position (and optionally velocity and acceleration) of a body over
/// a **single sub-interval** of a DATA RECORD. The interval is defined
/// by `start_jd` and `end_jd` (Julian Dates, TDB).
///
/// Fields
/// --------
/// * `start_jd` — Julian Date (TDB) marking the beginning of the sub-interval.
/// * `end_jd` — Julian Date (TDB) marking the end of the sub-interval.
/// * `x`, `y`, `z` — Vectors of Chebyshev coefficients for each spatial
///   component of the body’s state.
///
/// Notes
/// --------
/// * Coefficients are parsed from the raw Horizons binary using the
///   IPT table (`offset`, `n_coeffs`, `n_subs`).
/// * A `HorizonRecord` is never constructed directly by users; it is
///   built internally by \[`extract_body_records`\] when parsing DATA RECORDS.
///
/// See also
/// --------
/// * [`InterpResult`] — container for evaluated position/velocity/acceleration.
#[derive(Debug, PartialEq, Clone)]
pub struct HorizonRecord {
    pub start_jd: f64,
    pub end_jd: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl HorizonRecord {
    /// Construct a new [`HorizonRecord`] from the raw coefficients of a DATA RECORD.
    ///
    /// Arguments
    /// -----------------
    /// * `start_jd` — Julian Date at the start of the sub-interval.
    /// * `end_jd` — Julian Date at the end of the sub-interval.
    /// * `coeffs` — Slice of Chebyshev coefficients extracted from the DATA RECORD.
    /// * `offset` — Offset in `coeffs` to the start of this body’s coefficients.
    /// * `n_subintervals` — Total number of sub-intervals in the parent DATA RECORD.
    /// * `n_coeffs` — Number of Chebyshev coefficients per coordinate axis.
    ///
    /// Return
    /// ----------
    /// * A new `HorizonRecord` instance containing the three coefficient arrays
    ///   (`x`, `y`, `z`) for this sub-interval.
    ///
    /// Panics
    /// ----------
    /// * If the computed index range exceeds the bounds of `coeffs`.
    ///
    /// Notes
    /// ----------
    /// * The indexing formula `offset - 3 + n_subintervals * n_coeffs * 3`
    ///   accounts for the fact that the first two values in each block are
    ///   `[JD_start, JD_end]` and the coefficients are stored contiguously
    ///   per body, per component.
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

    /// Interpolate the state of the body within this sub-interval.
    ///
    /// Evaluates the Chebyshev polynomial series at normalized time `tau`,
    /// returning position and optionally velocity and acceleration vectors.
    ///
    /// Arguments
    /// -----------------
    /// * `tau` — Normalized time ∈ \[0,1\], relative to [`start_jd`, `end_jd`].
    /// * `compute_velocity` — If `true`, compute the velocity vector.
    /// * `compute_acceleration` — If `true`, compute the acceleration vector.
    /// * `n_subintervals` — Number of sub-intervals in the parent DATA RECORD,
    ///   used for scaling derivatives.
    ///
    /// Return
    /// ----------
    /// * An [`InterpResult`] containing:
    ///   - `position` (`Vector3<f64>` in km),
    ///   - `velocity` (`Option<Vector3<f64>>`, km/day),
    ///   - `acceleration` (`Option<Vector3<f64>>`, km/day²).
    ///
    /// Notes
    /// ----------
    /// * The Chebyshev basis functions are generated up to `n_coeffs`.
    /// * Velocity and acceleration are scaled by factors derived from
    ///   the interval length and number of sub-intervals.
    /// * This method is low-level: most users will prefer
    ///   [`HorizonData::ephemeris`](super::horizon_data::HorizonData::ephemeris).
    ///
    /// See also
    /// --------
    /// * [`InterpResult`] — evaluated result container.
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

        if compute_velocity && n_coeff > 1 {
            for i in 3..n_coeff {
                tcheb_deriv[i] =
                    twot * tcheb_deriv[i - 1] + 2.0 * tcheb[i - 1] - tcheb_deriv[i - 2];
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
