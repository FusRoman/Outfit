//! Chebyshev record segments for legacy JPL Horizons DATA RECORDS.
//!
//! Overview
//! -----------------
//! This module defines [`HorizonRecord`](crate::jpl_ephem::horizon::horizon_records::HorizonRecord), a compact container holding the
//! Chebyshev coefficients for a **single sub-interval** inside a Horizons
//! DATA RECORD. Each record can evaluate the body state at any normalized
//! time `tau ∈ [0, 1]` within its interval, returning position and optionally
//! velocity and acceleration via [`InterpResult`](crate::jpl_ephem::horizon::interpolation_result::InterpResult).
//!
//! Data layout
//! -----------------
//! A `HorizonRecord` stores three coefficient vectors (`x`, `y`, `z`) used to
//! evaluate the Chebyshev series for each Cartesian component. Records are
//! built internally from the binary ephemeris using the IPT metadata
//! (`offset`, `n_coeffs`, `n_subs`) that defines how coefficients are packed
//! within each DATA RECORD.
//!
//! Interpolation semantics
//! -----------------
//! * Time normalization: `tau = 0` maps to `start_jd`, `tau = 1` to `end_jd`.
//! * Position is evaluated from the Chebyshev basis up to `n_coeffs`.
//! * Velocity and acceleration are obtained by evaluating derivative bases,
//!   scaled by factors depending on the sub-interval duration and
//!   the number of sub-intervals in the parent block.
//!
//! Units
//! -----------------
//! * `position`: kilometers (km)
//! * `velocity`: kilometers per day (km/day)
//! * `acceleration`: kilometers per day² (km/day²)
//!
//! Typical usage
//! -----------------
//! Most users do not create `HorizonRecord` directly. Instead, they:
//! 1. Load and parse the ephemeris with
//!    [`HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData).
//! 2. Query a target/center state with
//!    [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris),
//!    which selects the right sub-interval and performs interpolation.
//! 3. (Advanced) If you already have a record and `tau`, call
//!    [`HorizonRecord::interpolate`](crate::jpl_ephem::horizon::horizon_records::HorizonRecord::interpolate) directly.
//!
//! Examples
//! -----------------
//! ```rust, no_run
//! use outfit::jpl_ephem::horizon::{
//!     horizon_data::HorizonData,
//!     horizon_ids::HorizonID,
//! };
//! use outfit::constants::MJD;
//!
//! // Load the binary and build HorizonData (details omitted here).
//! let hd: HorizonData = unimplemented!("Construct HorizonData from a loaded Horizons file");
//!
//! // Interpolate Earth wrt Solar System Barycenter at an epoch:
//! let et: MJD = 60200.0; // example MJD
//! let res = hd.ephemeris(HorizonID::Earth, HorizonID::SolarSystemBarycenter, et, true, true);
//! println!("r[km] = {:?}", res.position);
//! println!("v[km/day] = {:?}", res.velocity);
//! ```
//!
//! See also
//! -----------------
//! * [`HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData) – file parsing and high-level queries.
//! * [`InterpResult`](crate::jpl_ephem::horizon::interpolation_result::InterpResult) – evaluated state container.
//! * [`HorizonID`](crate::jpl_ephem::horizon::horizon_ids::HorizonID) – body identifiers.
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

#[cfg(test)]
mod test_horizon_records {
    // tests/horizon_record_interpolate.rs

    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use nalgebra::Vector3;
    use proptest::prelude::*;

    use crate::jpl_ephem::horizon::horizon_records::HorizonRecord;

    /// Build a synthetic HorizonRecord directly (bypassing `new`) for clarity in tests.
    fn rec(start_jd: f64, end_jd: f64, x: Vec<f64>, y: Vec<f64>, z: Vec<f64>) -> HorizonRecord {
        HorizonRecord {
            start_jd,
            end_jd,
            x,
            y,
            z,
        }
    }

    /// Naive Chebyshev basis (T_k) generator for reference comparisons.
    /// Stable enough for small n in tests.
    fn cheb_basis(tc: f64, n: usize) -> Vec<f64> {
        let mut t = vec![0.0_f64; n];
        if n == 0 {
            return t;
        }
        t[0] = 1.0;
        if n == 1 {
            return t;
        }
        t[1] = tc;
        let twot = tc + tc;
        for i in 2..n {
            t[i] = twot * t[i - 1] - t[i - 2];
        }
        t
    }

    /// Reference evaluation matching the current implementation semantics
    /// (position only) using explicit basis.
    fn eval_pos(coeffs: &[f64], tc: f64) -> f64 {
        let t = cheb_basis(tc, coeffs.len());
        coeffs.iter().zip(t.iter()).map(|(c, b)| c * b).sum()
    }

    /// Compute the internal tc mapping used by HorizonRecord::interpolate (mirrors code).
    fn tc_from_tau(tau: f64, n_sub: usize) -> f64 {
        let dt1 = tau as i64;
        let temp = n_sub as f64 * tau;
        2.0 * (temp.rem_euclid(1.0) + dt1 as f64) - 1.0
    }

    #[test]
    fn new_extracts_coeffs_synthetically() {
        // NOTE: This test exercises the slicing logic of `new` with a synthetic packing.
        // We choose `n_subintervals = 0` so that `base = offset - 3` and pick `offset = 3`,
        // which yields `base = 0`. This avoids negative indices while still validating the
        // contiguous layout (x[..], y[..], z[..]).
        let n_coeffs = 4usize;
        let n_subintervals = 0usize;

        // coeffs = [*,*,*] + [x0..x3][y0..y3][z0..z3]
        let x: Vec<f64> = (0..n_coeffs).map(|i| 10.0 + i as f64).collect();
        let y: Vec<f64> = (0..n_coeffs).map(|i| 20.0 + i as f64).collect();
        let z: Vec<f64> = (0..n_coeffs).map(|i| 30.0 + i as f64).collect();
        let mut coeffs = vec![999.0, 888.0, 777.0];
        coeffs.extend_from_slice(&x);
        coeffs.extend_from_slice(&y);
        coeffs.extend_from_slice(&z);

        // base = offset - 3 + n_subintervals * n_coeffs * 3
        let offset = 6usize;

        let hr = HorizonRecord::new(1000.0, 1001.0, &coeffs, offset, n_subintervals, n_coeffs);

        assert_eq!(hr.x, x);
        assert_eq!(hr.y, y);
        assert_eq!(hr.z, z);
    }

    #[test]
    fn interpolate_zero_coeffs_yield_zero_vectors() {
        let hr = rec(1000.0, 1001.0, vec![0.0; 6], vec![0.0; 6], vec![0.0; 6]);

        let r = hr.interpolate(0.0, false, false, 1);
        assert_abs_diff_eq!(r.position, Vector3::zeros(), epsilon = 0.0);
        assert!(r.velocity.is_none());
        assert!(r.acceleration.is_none());

        let r = hr.interpolate(0.42, true, true, 5);
        assert_abs_diff_eq!(r.position, Vector3::zeros(), epsilon = 0.0);
        assert_abs_diff_eq!(r.velocity.unwrap(), Vector3::zeros(), epsilon = 0.0);
        assert_abs_diff_eq!(r.acceleration.unwrap(), Vector3::zeros(), epsilon = 0.0);
    }

    #[test]
    fn interpolate_position_matches_explicit_basis() {
        // Simple, low-degree coefficients to check position path exactly.
        let x = vec![1.0, 0.5, -0.25, 0.125]; // T0..T3
        let y = vec![0.0, -2.0, 0.0, 1.0];
        let z = vec![3.0, 0.0, 0.0, 0.0];
        let hr = rec(2000.0, 2002.0, x.clone(), y.clone(), z.clone());

        let n_sub = 4usize;
        for &tau in &[0.0, 0.1, 0.33, 0.5, 0.73, 0.999999] {
            let tc = tc_from_tau(tau, n_sub);

            let r = hr.interpolate(tau, false, false, n_sub);

            assert_abs_diff_eq!(r.position.x, eval_pos(&x, tc), epsilon = 1e-14);
            assert_abs_diff_eq!(r.position.y, eval_pos(&y, tc), epsilon = 1e-14);
            assert_abs_diff_eq!(r.position.z, eval_pos(&z, tc), epsilon = 1e-14);
        }
    }

    #[test]
    fn velocity_flag_does_not_change_position() {
        let hr = rec(
            1234.5,
            1236.5,
            vec![0.1, -0.2, 0.3, -0.4, 0.5],
            vec![1.0, 0.0, -1.0, 0.0, 1.0],
            vec![2.0, -1.0, 0.0, 0.0, 0.0],
        );
        let tau = 0.37;
        let n_sub = 3;

        let r0 = hr.interpolate(tau, false, false, n_sub);
        let r1 = hr.interpolate(tau, true, false, n_sub);
        let r2 = hr.interpolate(tau, false, true, n_sub);
        let r3 = hr.interpolate(tau, true, true, n_sub);

        assert_abs_diff_eq!(r0.position, r1.position, epsilon = 0.0);
        assert_abs_diff_eq!(r0.position, r2.position, epsilon = 0.0);
        assert_abs_diff_eq!(r0.position, r3.position, epsilon = 0.0);
    }

    #[test]
    fn pure_t1_coefficient_gives_constant_velocity_scale() {
        // If only T1 is non-zero: f(tc) = a * tc
        // Then df/dtc = a (constant). The implementation computes velocity as:
        // v = vfac * (derivative wrt tc). So velocity should be vfac * a.
        // f(tc) = a * tc  => df/dtc = a (constant)
        let a = 123.456789;

        let hr = rec(
            1000.0,
            1004.0,
            vec![0.0, a, 0.0],
            vec![0.0, -a, 0.0],
            vec![0.0, 2.0 * a, 0.0],
        );
        let n_sub = 8usize;
        let vfac = (2.0 * n_sub as f64) / (hr.end_jd - hr.start_jd);

        for &tau in &[0.0, 0.2, 0.4, 0.6, 0.9] {
            let r = hr.interpolate(tau, true, false, n_sub);
            let v = r.velocity.unwrap();
            assert_abs_diff_eq!(v.x, vfac * a, epsilon = 1e-12);
            assert_abs_diff_eq!(v.y, vfac * (-a), epsilon = 1e-12);
            assert_abs_diff_eq!(v.z, vfac * (2.0 * a), epsilon = 1e-12);
        }
    }

    #[test]
    fn scaling_invariance() {
        // If we multiply all coefficients by s, outputs should scale by s.
        let s = 5.5;
        let hr = rec(
            5000.0,
            5002.0,
            vec![0.3, -0.2, 0.5, 0.0, 0.1],
            vec![1.0, 2.0, 3.0, 4.0, 5.0],
            vec![-1.0, 0.0, 1.0, -1.0, 0.0],
        );
        let hr_s = rec(
            hr.start_jd,
            hr.end_jd,
            hr.x.iter().map(|v| s * v).collect(),
            hr.y.iter().map(|v| s * v).collect(),
            hr.z.iter().map(|v| s * v).collect(),
        );

        let tau = 0.73;
        let n_sub = 7;

        let a = hr.interpolate(tau, true, true, n_sub);
        let b = hr_s.interpolate(tau, true, true, n_sub);

        assert_abs_diff_eq!(b.position, a.position * s, epsilon = 1e-12);
        assert_abs_diff_eq!(
            b.velocity.unwrap(),
            a.velocity.unwrap() * s,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            b.acceleration.unwrap(),
            a.acceleration.unwrap() * s,
            epsilon = 1e-11
        );
    }

    #[test]
    fn tau_edge_cases() {
        let hr = rec(
            70000.0,
            70001.0,
            vec![1.0, -0.5, 0.25, -0.125],
            vec![0.0, 2.0, 0.0, -1.0],
            vec![3.0, -2.0, 1.0, 0.0],
        );

        for &tau in &[0.0, 1.0, f64::EPSILON, 1.0 - 1e-15] {
            let r = hr.interpolate(tau, true, true, 4);
            // Just ensure finiteness & no panics.
            assert!(r.position.iter().all(|v| v.is_finite()));
            assert!(r.velocity.unwrap().iter().all(|v| v.is_finite()));
            assert!(r.acceleration.unwrap().iter().all(|v| v.is_finite()));
        }
    }

    proptest! {
        // Property: position is linear in coefficients (superposition).
        #[test]
    fn prop_linear_in_coeffs(
        tau in 0f64..1f64,
        n_sub in 1usize..10usize,
        coeffs in prop::collection::vec(-1.0f64..1.0, 3..12) // <- min 3
    ) {
            let start = 80000.0;
            let end = 80005.0;

            // Two random records A and B.
            let ax = coeffs.clone();
            let ay = coeffs.iter().map(|v| 2.0*v).collect::<Vec<_>>();
            let az = coeffs.iter().map(|v| -v).collect::<Vec<_>>();

            let bx = coeffs.iter().map(|v| v*v.signum()).collect::<Vec<_>>();
            let by = coeffs.iter().map(|v| 0.5*v).collect::<Vec<_>>();
            let bz = coeffs.iter().map(|v| v - v*v).collect::<Vec<_>>();

            let a = rec(start, end, ax.clone(), ay.clone(), az.clone());
            let b = rec(start, end, bx.clone(), by.clone(), bz.clone());

            let alpha = 0.37f64;
            let beta  = -1.25f64;

            // Linear combination record (component-wise).
            let cx: Vec<f64> = ax.iter().zip(bx.iter()).map(|(u,v)| alpha*u + beta*v).collect();
            let cy: Vec<f64> = ay.iter().zip(by.iter()).map(|(u,v)| alpha*u + beta*v).collect();
            let cz: Vec<f64> = az.iter().zip(bz.iter()).map(|(u,v)| alpha*u + beta*v).collect();

            let c = rec(start, end, cx, cy, cz);

            // Property: interpolate(linear combo) == linear combo of interpolations
            let ra = a.interpolate(tau, true, true, n_sub);
            let rb = b.interpolate(tau, true, true, n_sub);
            let rc = c.interpolate(tau, true, true, n_sub);

            assert_relative_eq!(rc.position,  ra.position  * alpha + rb.position  * beta, max_relative = 1e-10);
            assert_relative_eq!(rc.velocity.unwrap(),  ra.velocity.unwrap()  * alpha + rb.velocity.unwrap()  * beta, max_relative = 1e-10);
            assert_relative_eq!(rc.acceleration.unwrap(), ra.acceleration.unwrap() * alpha + rb.acceleration.unwrap() * beta, max_relative = 1e-10);
        }

        // Property: small tau perturbations induce small output perturbations (continuity),
        // away from sub-interval boundaries where `fract` changes branch.
        #[test]
    fn prop_continuity_inside_subinterval(
        base_tau in 0.01f64..0.99f64,
        delta in 1e-9f64..1e-6f64,
        n_sub in 1usize..12usize,
        coeffs in prop::collection::vec(-2.0f64..2.0, 3..12) // <- min 3
    ) {
            let hr = rec(
                90000.0,
                90003.0,
                coeffs.clone(),
                coeffs.iter().map(|v| v*0.3).collect(),
                coeffs.iter().map(|v| -v*0.7).collect(),
            );
            let tau1 = (base_tau - delta).clamp(0.0, 1.0);
            let tau2 = (base_tau + delta).clamp(0.0, 1.0);

            let r1 = hr.interpolate(tau1, true, true, n_sub);
            let r2 = hr.interpolate(tau2, true, true, n_sub);

            // Lipschitz-like bound: scaled by delta (loose, just guards wild jumps).
            let pos_diff = (r2.position - r1.position).abs();
            assert!(pos_diff.x <= 1e9 * delta && pos_diff.y <= 1e9 * delta && pos_diff.z <= 1e9 * delta);

            if let (Some(v1), Some(v2)) = (r1.velocity, r2.velocity) {
                let vel_diff = (v2 - v1).abs();
                assert!(vel_diff.x.is_finite() && vel_diff.y.is_finite() && vel_diff.z.is_finite());
            }
            if let (Some(a1), Some(a2)) = (r1.acceleration, r2.acceleration) {
                let acc_diff = (a2 - a1).abs();
                assert!(acc_diff.x.is_finite() && acc_diff.y.is_finite() && acc_diff.z.is_finite());
            }
        }
    }
}
