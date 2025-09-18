//! # Parsing, conversion, and formatting of astronomical coordinates
//!
//! This module provides robust utilities to:
//!
//! - **Parse** textual Right Ascension (RA) and Declination (DEC) from common sexagesimal strings,
//! - **Estimate** per-measurement accuracy from the number of decimals in the input,
//! - **Convert** 3D Cartesian vectors to equatorial angles (RA, DEC),
//! - **Format** angles and vectors for human-friendly displays (sexagesimal H/M/S and D/M/S, AU vectors).
//!
//! ## Overview
//!
//! These helpers are typically used when ingesting MPC-style astrometric observations and
//! converting them into numerical values usable by orbit determination routines.
//!
//! ### Provided features
//!
//! **Parsing & accuracy**
//! -----------------
//! * [`parse_ra_to_deg`](crate::conversion::parse_ra_to_deg) — Parse a sexagesimal RA string (`"HH MM SS.SS"`) into **degrees** and **accuracy** (arcsec).
//! * [`parse_dec_to_deg`](crate::conversion::parse_dec_to_deg) — Parse a sexagesimal DEC string (`"±DD MM SS.SS"`) into **degrees** and **accuracy** (arcsec).
//! * [`arcsec_to_rad`](crate::conversion::arcsec_to_rad) — Convert **arcseconds** to **radians** (utility used by callers).
//!
//! **Angle & vector formatting**
//! -----------------
//! * [`ra_hms_prec`](crate::conversion::ra_hms_prec) — Convert RA (radians) to canonical `(HH, MM, SS.sss)` with rounding and carry, hours wrapped to `[0, 24)`.
//! * [`dec_sdms_prec`](crate::conversion::dec_sdms_prec) — Convert DEC (radians) to `(sign, DD, MM, SS.sss)` with rounding, carry, and clamping at `±90°`.
//! * [`fmt_vec3_au`](crate::conversion::fmt_vec3_au) — Format a `nalgebra::Vector3<f64>` as `[ x, y, z ] AU` with fixed decimal precision.
//!
//! **Cartesian → Equatorial**
//! -----------------
//! * [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) — Convert a 3D Cartesian position vector to `(α, δ, ρ)` where `α, δ` are in **radians**,
//!   and `ρ` is the input norm (same units as the input vector).
//!
//! ### Units
//!
//! - **RA** returned by [`parse_ra_to_deg`](crate::conversion::parse_ra_to_deg) is in **degrees** (`0° ≤ RA < 360°`).
//! - **DEC** returned by [`parse_dec_to_deg`](crate::conversion::parse_dec_to_deg) is in **degrees** (`−90° ≤ DEC ≤ +90°`).
//! - **Accuracy** estimates are in **arcseconds** (derived from the decimal precision of the input seconds).
//! - [`cartesian_to_radec`](crate::conversion::cartesian_to_radec) returns angles in **radians** and `ρ` in the **same unit** as the input vector.
//! - [`ra_hms_prec`](crate::conversion::ra_hms_prec) returns `(HH, MM, SS)` with `HH ∈ [0, 23]`, `MM ∈ [0, 59]`, `SS ∈ [0.0, 60.0)` after carry.
//! - [`dec_sdms_prec`](crate::conversion::dec_sdms_prec) returns `(sign, DD, MM, SS)` with `sign ∈ {'+','-'}`, `DD ∈ [0, 90]`, `MM ∈ [0, 59]`,
//!   `SS ∈ [0.0, 60.0)` after carry, and clamps to `90°00′00″` at the pole.
//! - [`fmt_vec3_au`](crate::conversion::fmt_vec3_au) prints components in **fixed-point** with exactly `prec` decimals, suffixed with `" AU"`.
//!
//! ### Accuracy estimation
//!
//! Both parsers compute an **accuracy hint** directly from the number of digits after the decimal point
//! in the *seconds* field. For instance, `"12 34 56.7"` implies `0.1″`, whereas `"12 34 56"` implies `1″`.
//! This is convenient when deriving per-observation weights for orbit determination.
//!
//! ## Examples
//!
//! ```rust,no_run
//! use outfit::conversion::{parse_ra_to_deg, parse_dec_to_deg, cartesian_to_radec,
//!                          ra_hms_prec, dec_sdms_prec, fmt_vec3_au};
//! use nalgebra::Vector3;
//!
//! // --- Parsing with accuracy ------------------------------------------------
//! let (ra_deg, ra_acc)   = parse_ra_to_deg("10 12 33.44").unwrap();   // deg, arcsec
//! let (dec_deg, dec_acc) = parse_dec_to_deg("-20 33 10.5").unwrap();  // deg, arcsec
//! println!("RA = {ra_deg:.6}° ± {ra_acc:.3}\"");
//! println!("DEC = {dec_deg:.6}° ± {dec_acc:.3}\"");
//!
//! // --- Cartesian → Equatorial (radians) ------------------------------------
//! let pos = Vector3::new(1.0, 1.0, 0.5);
//! let (alpha, delta, rho) = cartesian_to_radec(pos);
//! println!("α = {alpha} rad, δ = {delta} rad, ρ = {rho}");
//!
//! // --- Sexagesimal formatting helpers --------------------------------------
//! let (hh, mm, ss) = ra_hms_prec(alpha, 3);        // RA → (HH, MM, SS.sss)
//! let (sgn, d, m, s) = dec_sdms_prec(delta, 3);    // DEC → (sign, DD, MM, SS.sss)
//! println!("RA ≈ {hh:02}h{mm:02}m{ss:.3}s,  DEC ≈ {sgn}{d:02}°{m:02}'{s:.3}\"");
//!
//! // --- Vector formatting (AU) ----------------------------------------------
//! let r_geo = Vector3::new(0.1234567, -1.0, 2.0);
//! println!("{}", fmt_vec3_au(&r_geo, 6)); // → "[ 0.123457, -1.000000, 2.000000 ] AU"
//! ```
//!
//! ## See also
//!
//! - Sexagesimal string rendering of seconds: `fmt_ss` (in the observations display helpers).
//! - Reference frame conversions and aberration: `ref_system` and `observations` modules.
//! - MPC/ADES ingestion modules where these utilities are typically used.
use std::f64::consts::TAU;

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

/// Format a 3D vector (AU) with a configurable fixed decimal precision.
/// The output is rendered as:
/// `[ {x:.prec}, {y:.prec}, {z:.prec} ] AU`
///
/// Arguments
/// -----------------
/// * `v`: The position vector in **astronomical units (AU)**, expressed in the
///   **equatorial mean J2000** frame if you follow the crate’s convention.
/// * `prec`: Number of fractional digits for each component (fixed‐point).
///
/// Return
/// ----------
/// * A `String` like `"[ 0.123457, -1.000000, 0.000042 ] AU"` when `prec = 6`.
///
/// Notes
/// ----------
/// * Rounding uses Rust’s default `Display` formatting for `f64` (round half
///   away from zero).
/// * No thousands separator or scientific notation is used: components are
///   printed in **fixed‐point** with exactly `prec` decimals.
/// * The function does **not** sanitize non-finite values: `NaN`, `inf`, and
///   `-inf` will be forwarded as is.
/// * Units are **not converted**: call this function only if your vector is
///   already in AU. For other units (e.g. km), provide a separate formatter.
///
/// Examples
/// ----------
/// ```rust, ignore
/// use nalgebra::Vector3;
///
/// let v = Vector3::new(0.1234567, -1.0, 2.0);
/// assert_eq!(fmt_vec3_au(&v, 3), "[ 0.123, -1.000, 2.000 ] AU");
/// assert_eq!(fmt_vec3_au(&v, 6), "[ 0.123457, -1.000000, 2.000000 ] AU");
/// ```
///
/// See also
/// ------------
/// * [`fmt_ss`](crate::time::fmt_ss) – Zero-padded seconds string for sexagesimal outputs.
/// * [`ra_hms_prec`] – RA (radians) → `(HH, MM, SS.sss)` with carry.
/// * [`dec_sdms_prec`] – DEC (radians) → `(sign, DD, MM, SS.sss)` with carry.
pub fn fmt_vec3_au(v: &Vector3<f64>, prec: usize) -> String {
    let x = v.x;
    let y = v.y;
    let z = v.z;
    let p = prec; // capture for dynamic precision
    format!("[ {x:.p$}, {y:.p$}, {z:.p$} ] AU")
}

// --- Angle formatting --------------------------------------------------------

/// Convert a right ascension (radians) to sexagesimal **hours–minutes–seconds**
/// with rounding and carry handling.
///
/// The input angle is normalized to `[0, 2π)` (i.e., modulo one full turn),
/// then converted to **total seconds** in `[0h, 24h)`. Seconds are rounded to
/// `prec` fractional digits; potential overflows are carried into minutes, and
/// minutes into hours (hours wrap modulo 24).
///
/// Arguments
/// -----------------
/// * `rad`: Right ascension in **radians**. Any finite value is accepted;
///   negatives and values ≥ 2π are normalized to `[0, 2π)`.
/// * `prec`: Number of fractional digits for the **seconds** component.
///
/// Return
/// ----------
/// * A tuple `(HH, MM, SS)` where:
///   - `HH ∈ [0, 23]`,
///   - `MM ∈ [0, 59]`,
///   - `SS ∈ [0.0, 60.0)` after rounding and carry,
///     guaranteeing canonical sexagesimal components ready for display.
///
/// Notes
/// ----------
/// * Rounding uses `f64::round` (half away from zero), then carry is applied:
///   `59.9995s` at `prec = 3` becomes `60.000s → +1 min`.
/// * Final hours are wrapped modulo 24; e.g. a value very close to `2π` can
///   round to `24h00m00s` which is reported as `00h00m00s`.
/// * This function **does not** format strings. For zero-padded seconds like
///   `"SS.sss"`, combine with [`fmt_ss`](crate::time::fmt_ss).
///
/// See also
/// ------------
/// * [`fmt_ss`](crate::time::fmt_ss) – Enforce zero-padded second strings (`"SS.sss"`).
/// * [`dec_sdms_prec`] – Declination to `sign, DD, MM, SS.sss`.
pub fn ra_hms_prec(rad: f64, prec: usize) -> (u32, u32, f64) {
    // Normalize RA to [0, 2π)
    let mut a = rad % TAU;
    if a < 0.0 {
        a += TAU;
    }

    // Total seconds in [0h, 24h)
    let total_sec = a * 12.0 / std::f64::consts::PI * 3600.0;

    let mut h = (total_sec / 3600.0).floor() as u32;
    let rem = total_sec - (h as f64) * 3600.0;
    let mut m = (rem / 60.0).floor() as u32;
    let mut s = rem - (m as f64) * 60.0;

    // Round seconds to precision
    let pow = 10f64.powi(prec as i32);
    s = (s * pow).round() / pow;

    // Carry seconds -> minutes
    if s >= 60.0 {
        s -= 60.0;
        m += 1;
    }
    // Carry minutes -> hours
    if m >= 60 {
        m = 0;
        h = (h + 1) % 24;
    }
    (h, m, s)
}

/// Convert a declination (radians) to **signed** sexagesimal
/// **degrees–minutes–seconds** with rounding, carry, and polar clamping.
///
/// The sign is taken from the input (`'+'` for `rad ≥ 0`, `'-'` otherwise).
/// The absolute value is converted to **total arcseconds**, seconds are rounded
/// to `prec` fractional digits, and carry is applied seconds→minutes→degrees.
/// Final values are clamped to the physical pole at **±90°00′00″**.
///
/// Arguments
/// -----------------
/// * `rad`: Declination in **radians**. Typical physical range is
///   `[-π/2, +π/2]`, but any finite value is accepted; the absolute value is
///   used for the `DD/MM/SS` decomposition.
/// * `prec`: Number of fractional digits for the **seconds** component.
///
/// Return
/// ----------
/// * A tuple `(sign, DD, MM, SS)` where:
///   - `sign ∈ {'+', '-'}` reflects the input sign,
///   - `DD ∈ [0, 90]`, `MM ∈ [0, 59]`,
///   - `SS ∈ [0.0, 60.0)` after rounding and carry,
///     with a **final clamp** at the pole: if rounding would exceed `90°`,
///     the function returns `(sign, 90, 0, 0.0)`.
///
/// Notes
/// ----------
/// * Rounding uses `f64::round` (half away from zero), then carry is applied:
///   e.g. `59.9995″` at `prec = 3` becomes `60.000″ → +1′`.
/// * If carrying minutes produces `60′`, it becomes `+1°`.
/// * At the upper bound, values that round past `90°` are clamped to
///   `90°00′00.000″` to maintain a valid declination.
///
/// See also
/// ------------
/// * [`ra_hms_prec`] – Right ascension to `HH, MM, SS.sss`.
/// * [`fmt_ss`](crate::time::fmt_ss) – Zero-padded second string formatting for display.
pub fn dec_sdms_prec(rad: f64, prec: usize) -> (char, u32, u32, f64) {
    let sign = if rad < 0.0 { '-' } else { '+' };
    let deg_abs = rad.abs() * 180.0 / std::f64::consts::PI;

    let mut d = deg_abs.floor() as u32;
    let rem = deg_abs - d as f64;
    let mut m = (rem * 60.0).floor() as u32;
    let mut s = (rem * 3600.0) - (m as f64) * 60.0;

    let pow = 10f64.powi(prec as i32);
    s = (s * pow).round() / pow;

    // Carry seconds -> minutes
    if s >= 60.0 {
        s -= 60.0;
        m += 1;
    }
    // Carry minutes -> degrees
    if m >= 60 {
        m = 0;
        d += 1;
    }
    // Safety clamp near the pole
    if d > 90 {
        d = 90;
        m = 0;
        s = 0.0;
    }
    (sign, d, m, s)
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
