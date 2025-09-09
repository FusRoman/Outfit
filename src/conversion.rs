//! # Parsing and conversion of astronomical coordinates
//!
//! This module provides utility functions for:
//!
//! - **Parsing textual Right Ascension (RA) and Declination (DEC)** from sexagesimal strings
//!   (formats like `HH MM SS.SS` or `±DD MM SS.SS`),
//! - **Estimating the positional accuracy** based on the number of decimals
//!   provided in the input strings,
//! - **Converting 3D Cartesian position vectors to equatorial coordinates (RA, DEC)**.
//!
//! ## Overview
//!
//! These functions are used to transform raw observational data (often stored
//! in MPC-style text files) into numerical coordinates that can be consumed
//! by orbit determination algorithms.
//!
//! ### Provided features
//!
//! * [`parse_ra_to_deg`](crate::conversion::parse_ra_to_deg) – Parse a sexagesimal RA string into degrees and accuracy in arcseconds.
//! * [`parse_dec_to_deg`](crate::conversion::parse_dec_to_deg) – Parse a sexagesimal DEC string into degrees and accuracy in arcseconds.
//! * [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) – Convert a 3D Cartesian position vector to right ascension and declination (in radians).
//!
//! ### Accuracy estimation
//!
//! Both parsing functions compute a **measurement accuracy** in arcseconds
//! directly from the number of digits present after the decimal point in
//! the seconds field. This is useful when weighting observations in the
//! orbit determination process.
//!
//! ### Units
//!
//! - RA is returned in **degrees** (0° ≤ RA < 360°).
//! - DEC is returned in **degrees** (−90° ≤ DEC ≤ +90°).
//! - Accuracy is returned in **arcseconds**.
//! - Cartesian inputs for [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) can be in any length unit,
//!   but the result angles are in **radians** and the norm is returned
//!   in the same unit as the input.
//!
//! ## Example
//!
//! ```rust,ignore
//! use outfit::conversion::{parse_ra_to_deg, parse_dec_to_deg, cartesian_to_radec};
//! use nalgebra::Vector3;
//!
//! // Parse RA and DEC strings
//! let (ra_deg, ra_acc) = parse_ra_to_deg("10 12 33.44").unwrap();
//! let (dec_deg, dec_acc) = parse_dec_to_deg("-20 33 10.5").unwrap();
//! println!("RA = {ra_deg} deg ± {ra_acc} arcsec");
//! println!("DEC = {dec_deg} deg ± {dec_acc} arcsec");
//!
//! // Convert Cartesian vector to RA/DEC (radians)
//! let pos = Vector3::new(1.0, 1.0, 0.5);
//! let (alpha, delta, rho) = cartesian_to_radec(pos);
//! println!("α = {alpha} rad, δ = {delta} rad, distance = {rho}");
//! ```
//!
//! These functions are typically called when reading observations from
//! MPC 80-column formatted files or other astrometric data sources.

use nalgebra::Vector3;

use crate::constants::{ArcSec, Degree, Radian, DPI};

/// Estimate the accuracy of a numeric string based on its decimal precision.
///
/// Arguments
/// ---------------
/// * `field`: a string slice containing the numeric value (e.g., `"56.78"`), typically the last component of an angle
/// * `factor`: a scale factor to apply to the accuracy (e.g., `1.0 / 60.0` for arcminutes, `1.0 / 3600.0` for arcseconds)
///
/// Return
/// ----------
/// * `Option<f64>`: the estimated accuracy scaled by `factor`, or `None` if the input is malformed
fn compute_accuracy(field: &str, factor: f64) -> Option<f64> {
    if let Some(dot_pos) = field.find('.') {
        let length = field.trim().len();
        let digits_after_dot = length - dot_pos - 1;
        let acc = 10f64.powi(-(digits_after_dot as i32)) * factor;
        Some(acc)
    } else {
        Some(1.0 * factor)
    }
}

/// Convert an angle from **arcseconds** to **radians**.
///
/// Arguments
/// -----------------
/// * `arcsec` — Angle in **arcseconds**.
///
/// Return
/// ----------
/// * Angle in **radians** (`Radian`).
#[inline]
pub fn arcsec_to_rad(arcsec: ArcSec) -> Radian {
    (arcsec / 3600.0).to_radians()
}

/// Parse a right ascension (RA) string and convert it to degrees, with an estimate of its accuracy.
///
/// This function parses a right ascension expressed in **sexagesimal hours**
/// (`HH MM SS.SS`) and converts it into degrees. The input must have exactly
/// three whitespace-separated fields: hours, minutes, and seconds (which may
/// include fractional seconds).
///
/// # Arguments
///
/// * `ra` – A string in the format:
///   * `"HH MM SS.SS"` (e.g., `"12 30 45.67"`).
///
/// # Returns
///
/// Returns `Some((ra_deg, accuracy_arcsec))` where:
/// * `ra_deg` – Right ascension in **degrees** (0° ≤ RA < 360°),
/// * `accuracy_arcsec` – Estimated accuracy of the RA, derived from the number
///   of decimals provided in the seconds field, in **arcseconds**.
///
/// Returns `None` if:
/// * The string does not have exactly 3 whitespace-separated components,
/// * Any component fails to parse as a floating-point number.
///
/// # Formula
///
/// ```text
/// RA(deg) = (hours + minutes / 60 + seconds / 3600) × 15
/// ```
///
/// # See also
/// * [`parse_dec_to_deg`] – Parses declination strings into degrees.
pub fn parse_ra_to_deg(ra: &str) -> Option<(Degree, ArcSec)> {
    let parts: Vec<&str> = ra.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }

    let h: f64 = parts[0].parse().ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s_raw = parts[2];
    let s: f64 = s_raw.parse().ok()?;

    let ra_deg = (h + m / 60.0 + s / 3600.0) * 15.0;
    let acc_arcsec = compute_accuracy(s_raw, 1.0 / 3600.0)?;
    Some((ra_deg, acc_arcsec))
}

/// Parse a declination (DEC) string and convert it to degrees, with an estimate of its accuracy.
///
/// This function parses a declination expressed in **sexagesimal degrees**
/// (`±DD MM SS.SS`) and converts it into degrees. The input must have exactly
/// three whitespace-separated fields: degrees (with sign), minutes, and seconds
/// (which may include fractional seconds).
///
/// # Arguments
///
/// * `dec` – A string in the format:
///   * `"±DD MM SS.SS"` (e.g., `"-23 26 45.1"` or `"+10 15 30"`).
///
/// # Returns
///
/// Returns `Some((dec_deg, accuracy_arcsec))` where:
/// * `dec_deg` – Declination in **degrees** (−90° ≤ DEC ≤ +90°),
/// * `accuracy_arcsec` – Estimated accuracy of the DEC, derived from the number
///   of decimals provided in the seconds field, in **arcseconds**.
///
/// Returns `None` if:
/// * The string does not have exactly 3 whitespace-separated components,
/// * Any component fails to parse as a floating-point number.
///
/// # Formula
///
/// ```text
/// DEC(deg) = sign × (degrees + minutes / 60 + seconds / 3600)
/// ```
///
/// # See also
/// * [`parse_ra_to_deg`] – Parses right ascension strings into degrees.
pub fn parse_dec_to_deg(dec: &str) -> Option<(Degree, ArcSec)> {
    let parts: Vec<&str> = dec.split_whitespace().collect();
    if parts.len() != 3 {
        return None;
    }

    let sign = if parts[0].starts_with('-') { -1.0 } else { 1.0 };
    let d: f64 = parts[0].trim_start_matches(&['-', '+'][..]).parse().ok()?;
    let m: f64 = parts[1].parse().ok()?;
    let s_raw = parts[2];
    let s: f64 = s_raw.parse().ok()?;

    let dec_deg = sign * (d + m / 60.0 + s / 3600.0);
    let acc_arcsec = compute_accuracy(s_raw, 1. / 3600.)?;
    Some((dec_deg, acc_arcsec))
}

/// Convert a 3D Cartesian position vector to right ascension and declination.
///
/// Given a position vector expressed in Cartesian coordinates (typically in an equatorial frame),
/// this function returns the corresponding right ascension (α), declination (δ), and norm (distance).
///
/// Arguments
/// ---------
/// * `cartesian_position`: 3D position vector in Cartesian coordinates [AU or any length unit].
///
/// Returns
/// --------
/// * Tuple `(α, δ, ρ)`:
///     - `α`: right ascension in radians, in the range [0, 2π).
///     - `δ`: declination in radians, in the range [−π/2, +π/2].
///     - `ρ`: Euclidean norm of the vector (distance to the origin).
///
/// Remarks
/// -------
/// * If the input vector has zero norm, the result is `(0.0, 0.0, 0.0)`.
/// * The RA computation uses `atan2` to preserve quadrant information.
/// * This function is used when converting inertial position vectors to observable angles.
///
/// # See also
/// * [`correct_aberration`](crate::observations::correct_aberration) – apply aberration correction before calling this if needed
pub fn cartesian_to_radec(cartesian_position: Vector3<f64>) -> (f64, f64, f64) {
    let pos_norm = cartesian_position.norm();
    if pos_norm == 0. {
        return (0.0, 0.0, pos_norm);
    }

    let delta = (cartesian_position.z / pos_norm).asin();

    let cos_delta = delta.cos();
    if cos_delta == 0.0 {
        return (0.0, delta, pos_norm);
    }

    let cos_alpha = cartesian_position.x / (pos_norm * cos_delta);
    let sin_alpha = cartesian_position.y / (pos_norm * cos_delta);
    let alpha = sin_alpha.atan2(cos_alpha);
    let alpha = if alpha < 0.0 { alpha + DPI } else { alpha };
    (alpha, delta, pos_norm)
}

#[cfg(test)]
mod observations_test {
    use super::*;

    #[test]
    fn test_ra_to_deg() {
        assert_eq!(
            parse_ra_to_deg("22 52 23.37"),
            Some((343.097375, 2.777777777777778e-6))
        );
        assert_eq!(
            parse_ra_to_deg("23 58 57.68"),
            Some((359.7403333333333, 2.777777777777778e-6))
        );
        assert_eq!(
            parse_ra_to_deg("04 41 04.77"),
            Some((70.269875, 2.777777777777778e-6))
        );
        assert_eq!(parse_ra_to_deg("1 2 3.4.5"), None);
        assert_eq!(parse_ra_to_deg("1 2"), None);
        assert_eq!(
            parse_ra_to_deg("06 50 13.370"),
            Some((102.55570833333333, 2.7777777777777776e-7))
        );
    }

    #[test]
    fn test_dec_to_deg() {
        assert_eq!(
            parse_dec_to_deg("-00 30 14.2"),
            Some((-0.5039444444444444, 2.777777777777778e-5))
        );
        assert_eq!(
            parse_dec_to_deg("+13 55 42.7"),
            Some((13.928527777777777, 2.777777777777778e-5))
        );
        assert_eq!(
            parse_dec_to_deg("89 15 50.2"),
            Some((89.26394444444445, 2.777777777777778e-5))
        );
        assert_eq!(parse_dec_to_deg("89 15 50.2.3"), None);
        assert_eq!(parse_dec_to_deg("89 15"), None);

        assert_eq!(
            parse_dec_to_deg("-14 47 05.4"),
            Some((-14.784833333333333, 2.777777777777778e-5))
        );
    }

    #[test]
    fn test_estimate_accuracy() {
        assert_eq!(
            compute_accuracy("23.3", 1. / 3600.),
            Some(2.777777777777778e-5)
        );
        assert_eq!(
            compute_accuracy("23", 1. / 3600.),
            Some(0.0002777777777777778)
        );
        assert_eq!(
            compute_accuracy("23.370", 1. / 3600.),
            Some(2.7777777777777776e-7)
        );
        assert_eq!(
            compute_accuracy("23.37", 1. / 3600.),
            Some(2.777777777777778e-6)
        );
    }
}
