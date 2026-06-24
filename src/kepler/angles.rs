//! Angle normalization utilities used by the preliminary Kepler solvers.
//!
//! These helpers keep angular quantities (eccentric anomaly, mean anomaly,
//! ...) inside well-defined ranges, so that differences between angles are
//! always taken along the shortest signed path.

use std::f64::consts::PI;

use crate::constants::DPI;

/// Normalize an angle in radians to the range `[0, 2π)`.
///
/// This ensures any input angle is wrapped into the principal interval
/// `0 ≤ θ < 2π` using Euclidean remainder.
pub fn principal_angle(angle_radians: f64) -> f64 {
    angle_radians.rem_euclid(DPI)
}

/// Compute the signed minimal difference between two angles in radians.
///
/// Returns the value of `(first_angle - second_angle)` wrapped into the range
/// `[-π, π]`, i.e. the smallest signed rotation from `second_angle` to
/// `first_angle`.
pub fn angle_diff(first_angle_radians: f64, second_angle_radians: f64) -> f64 {
    let normalized_first_angle = principal_angle(first_angle_radians);
    let normalized_second_angle = principal_angle(second_angle_radians);

    let mut signed_difference = normalized_first_angle - normalized_second_angle;
    if signed_difference > PI {
        signed_difference -= DPI;
    } else if signed_difference < -PI {
        signed_difference += DPI;
    }
    signed_difference
}
