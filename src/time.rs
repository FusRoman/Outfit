//! # Time utilities for astronomical computations
//!
//! Utilities to convert between common astronomical time representations
//! (**calendar date, MJD, JD**), format epochs as **ISO-8601** on various
//! time scales, and compute the **Greenwich Mean Sidereal Time (GMST)**.
//!
//! Built on top of [`hifitime`](https://docs.rs/hifitime) for precise epoch handling.
//! Supported time scales: **UTC**, **TT**, and **UT1** (for GMST input).
//!
//! ## Provided conversions & helpers
//!
//! **Calendar / JD / MJD**
//! -----------------
//! * [`date_to_mjd`](crate::time::date_to_mjd) — Parse ISO-8601 calendar strings (`YYYY-MM-DDTHH:mm:ss`) to **MJD (UTC)**.
//! * [`mjd_to_jd`](crate::time::mjd_to_jd) / [`jd_to_mjd`](crate::time::jd_to_mjd) — Convert **MJD (UTC)** ↔ **JD (UTC)**.
//! * [`frac_date_to_mjd`](crate::time::frac_date_to_mjd) — Parse `YYYY MM DD.FFFFF` (UTC) to **MJD (TT)**.
//!
//! **Epoch & ISO formatting**
//! -----------------
//! * [`iso_tt_from_epoch`](crate::time::iso_tt_from_epoch) — Format an [`Epoch`](hifitime::Epoch) as `YYYY-MM-DDThh:mm:SS.sss TT`.
//! * [`iso_utc_from_epoch`](crate::time::iso_utc_from_epoch) — Format an [`Epoch`](hifitime::Epoch) as `YYYY-MM-DDThh:mm:SS.sssZ` (UTC).
//! * [`fmt_ss`](crate::time::fmt_ss) — Two-digit, fixed-precision seconds string (e.g., `"05.123"`).
//!
//! **Batch scale conversions**
//! -----------------
//! * [`utc_mjd_slice_to_tt`](crate::time::utc_mjd_slice_to_tt) — Map a slice of **MJD (UTC)** → **MJD (TT)**.
//! * [`utc_jd_slice_to_tt_mjd`](crate::time::utc_jd_slice_to_tt_mjd) — Map a slice of **JD (UTC)** → **MJD (TT)**.
//!
//! **Sidereal time**
//! -----------------
//! * [`gmst`](crate::time::gmst) — GMST angle (radians) for a given **MJD (UT1)**, normalized to `[0, 2π)`.
//!
//! ## Units & conventions
//!
//! - **Calendar strings** are interpreted as **UTC** unless otherwise stated.
//! - **JD / MJD** are expressed in **days** (floating-point). By convention,
//!   `MJD = JD − 2_400_000.5`, origin at 1858-11-17T00:00:00.
//! - **GMST** is returned in **radians**, normalized to `[0, 2π)`; input is **MJD (UT1)**.
//! - **Epoch constants**: [`T2000`](crate::constants::T2000) = J2000.0 (MJD 51544.5, TT).
//! - **Angle constants**: [`DPI`](crate::constants::DPI) = `2π`.
//!
//! ## Examples
//!
//! ```rust, no_run
//! use outfit::time::*;
//! use hifitime::{Epoch, TimeScale};
//!
//! // --- Calendar → MJD (UTC) ------------------------------------------------
//! let mjd_utc = date_to_mjd(&vec!["2021-01-01T00:00:00"]);
//! assert_eq!(mjd_utc, vec![59215.0]);
//!
//! // --- MJD (UTC) ↔ JD (UTC) ------------------------------------------------
//! let jd_utc = mjd_to_jd(&mjd_utc);
//! assert_eq!(jd_utc, vec![2459215.5]);
//! let mjd_utc_back = jd_to_mjd(&jd_utc);
//! assert_eq!(mjd_utc_back, vec![59215.0]);
//!
//! // --- Batch UTC → TT (days) -----------------------------------------------
//! let mjd_tt_vec = utc_mjd_slice_to_tt(&mjd_utc);          // MJD(UTC) → MJD(TT)
//! let mjd_tt_from_jd = utc_jd_slice_to_tt_mjd(&jd_utc);     // JD(UTC) → MJD(TT)
//!
//! // --- Fractional-date (UTC) to MJD (TT) -----------------------------------
//! let mjd_tt = frac_date_to_mjd("1976 09 20.93878").unwrap();
//!
//! // --- ISO formatting (TT / UTC) -------------------------------------------
//! let epoch_tt = Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT);
//! let iso_tt  = iso_tt_from_epoch(epoch_tt, 3);   // "YYYY-MM-DDThh:mm:SS.sss TT"
//! let iso_utc = iso_utc_from_epoch(epoch_tt, 3);  // "YYYY-MM-DDThh:mm:SS.sssZ"
//!
//! // --- GMST (UT1) -----------------------------------------------------------
//! let gmst_angle = gmst(57028.478514610404);  // radians in [0, 2π)
//! assert!((gmst_angle - 4.851925725092499).abs() < 1e-12);
//!
//! // --- Seconds formatter for sexagesimal displays ---------------------------
//! let ss = fmt_ss(5.1234, 3);  // → "05.123"
//! assert_eq!(ss, "05.123");
//! ```
//!
//! ## References
//!
//! - IAU 1982; IERS Conventions 1996/2000  
//! - *Explanatory Supplement to the Astronomical Almanac* (1992)  
//! - Vallado, *Fundamentals of Astrodynamics and Applications* (2007)
//!
//! ## See also
//! -------------
//! * [`hifitime::Epoch`] — core epoch type for parsing and conversions.
//! * [`gmst`](crate::time::gmst) — Greenwich Mean Sidereal Time in radians.
//! * [`date_to_mjd`](crate::time::date_to_mjd) / [`mjd_to_jd`](crate::time::mjd_to_jd) / [`jd_to_mjd`](crate::time::jd_to_mjd) — basic conversions.
//! * [`iso_tt_from_epoch`](crate::time::iso_tt_from_epoch) / [`iso_utc_from_epoch`](crate::time::iso_utc_from_epoch) — ISO renderers.
//! * [`utc_mjd_slice_to_tt`](crate::time::utc_mjd_slice_to_tt) / [`utc_jd_slice_to_tt_mjd`](crate::time::utc_jd_slice_to_tt_mjd) — batch scale changes.
use hifitime::{Epoch, TimeScale};
use std::str::FromStr;

use crate::constants::{DPI, T2000};

/// Format seconds with a forced two-digit integer part, e.g. `"SS.sss"`.
/// This helper guarantees at least two digits before the decimal separator,
/// which is convenient for sexagesimal angle rendering (RA/DEC seconds).
///
/// Arguments
/// -----------------
/// * `seconds`: The seconds value to format. **Must** be within `[0, 60)`
///   *after any carry has already been applied*. This function does not perform carry.
/// * `prec`: Decimal precision (number of fractional digits).
///
/// Return
/// ----------
/// * A `String` representing the seconds with zero-padded integer part
///   and `prec` fractional digits (rounded half away from zero).
///
/// Notes
/// ----------
/// * Rounding uses `f64::round`, i.e. half away from zero.
/// * This function assumes callers like [`ra_hms_prec`](crate::conversion::ra_hms_prec) / [`dec_sdms_prec`](crate::conversion::dec_sdms_prec)
///   have already normalized and applied any carry so that `0 ≤ seconds < 60`.
///
/// See also
/// ------------
/// * [`ra_hms_prec`](crate::conversion::ra_hms_prec) – Carry-safe RA (radians) → `HH, MM, SS.sss`.
/// * [`dec_sdms_prec`](crate::conversion::dec_sdms_prec) – Carry-safe DEC (radians) → `sign, DD, MM, SS.sss`.
pub fn fmt_ss(seconds: f64, prec: usize) -> String {
    let pow = 10u64.pow(prec as u32);
    let total = (seconds * pow as f64).round() as u64;
    let int = total / pow;
    let frac = total % pow;
    format!("{int:02}.{frac:0prec$}")
}

/// Format an [`Epoch`] as an **ISO-8601** string on the **TT** time scale.
/// This performs a Gregorian calendar breakdown in TT and renders
/// `YYYY-MM-DDThh:mm:SS.sss TT`.
///
/// Arguments
/// -----------------
/// * `epoch_tt`: The epoch to format (interpreted on the TT time scale).
/// * `sec_prec`: Number of fractional digits for the seconds component.
///
/// Return
/// ----------
/// * A `String` like `"2000-01-01T12:00:00.000 TT"`.
///
/// Notes
/// ----------
/// * The seconds field is formatted via [`fmt_ss`] to enforce a two-digit
///   integer part and `sec_prec` fractional digits.
/// * The output is **not** locale-dependent and always uses zero padding.
///
/// See also
/// ------------
/// * [`iso_utc_from_epoch`] – UTC rendering, accounting for leap seconds.
#[inline]
pub fn iso_tt_from_epoch(epoch_tt: Epoch, sec_prec: usize) -> String {
    // Many hifitime versions expose `to_gregorian(TimeScale) -> (i32,u8,u8,u8,u8,f64, …)`
    let (y, mo, d, h, mi, s, _) = epoch_tt.to_gregorian(TimeScale::TT);
    let s_str = fmt_ss(s.into(), sec_prec);
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s_str} TT")
}

/// Format an [`Epoch`] as an **ISO-8601** string on the **UTC** time scale (with `'Z'` suffix).
/// Leap seconds are handled by `hifitime`’s internal tables.
///
/// Arguments
/// -----------------
/// * `epoch_tt`: The epoch to format; it is converted to **UTC** internally.
/// * `sec_prec`: Number of fractional digits for the seconds component.
///
/// Return
/// ----------
/// * A `String` like `"1998-01-01T00:00:00.000Z"`.
///
/// Notes
/// ----------
/// * The conversion TT → UTC is handled by `hifitime`, including leap seconds.
/// * The seconds field is formatted via [`fmt_ss`] to enforce a two-digit
///   integer part and `sec_prec` fractional digits.
///
/// See also
/// ------------
/// * [`iso_tt_from_epoch`] – TT rendering for dynamical times.
#[inline]
pub fn iso_utc_from_epoch(epoch_tt: Epoch, sec_prec: usize) -> String {
    let (y, mo, d, h, mi, s, _) = epoch_tt.to_gregorian(TimeScale::UTC);
    let s_str = fmt_ss(s.into(), sec_prec);
    // Use trailing 'Z' for UTC
    format!("{y:04}-{mo:02}-{d:02}T{h:02}:{mi:02}:{s_str}Z")
}

/// Transformation from date in the format YYYY-MM-ddTHH:mm:ss to modified julian date (MJD)
///
/// Argument
/// --------
/// * `date`: a vector of date in the format YYYY-MM-ddTHH:mm:ss
///
/// Return
/// ------
/// * a vector of float representing the input date in modified julian date (MJD)
pub fn date_to_mjd(date: &Vec<&str>) -> Vec<f64> {
    date.iter()
        .map(|x| Epoch::from_str(x).unwrap().to_mjd_utc_days())
        .collect::<Vec<f64>>()
}

/// Transformation from modified julian date (MJD) in julian date (JD)
///
/// Argument
/// --------
/// * `mjd`: a vector of MJD
///
/// Return
/// ------
/// * a vector of jd
pub fn mjd_to_jd(mjd: &[f64]) -> Vec<f64> {
    mjd.iter()
        .map(|x| Epoch::from_mjd_utc(*x).to_jde_utc_days())
        .collect()
}

/// Transformation from julian date (JD) in modified julian date (MJD)
///
/// Argument
/// --------
/// * `jd`: a vector of JD
///
/// Return
/// ------
/// * a vector of MJD
pub fn jd_to_mjd(jd: &[f64]) -> Vec<f64> {
    jd.iter()
        .map(|x| Epoch::from_jde_utc(*x).to_mjd_utc_days())
        .collect()
}

/// Map a slice of MJD(UTC) into a newly allocated Vec of MJD(TT).
///
/// Arguments
/// -----------------
/// * `mjds_utc`: Slice of MJD expressed on the UTC time scale.
///
/// Return
/// ----------
/// * A `Vec<f64>` containing the converted MJD(TT).
pub fn utc_mjd_slice_to_tt(mjds_utc: &[f64]) -> Vec<f64> {
    mjds_utc
        .iter()
        .map(|&m| Epoch::from_mjd_utc(m).to_mjd_tt_days())
        .collect()
}

/// Map a slice of JD(UTC) into a newly allocated Vec of MJD(TT).
///
/// Arguments
/// -----------------
/// * `jds_utc`: Slice of JD expressed on the UTC time scale.
///
/// Return
/// ----------
/// * A `Vec<f64>` containing the converted MJD(TT).
pub fn utc_jd_slice_to_tt_mjd(jds_utc: &[f64]) -> Vec<f64> {
    jds_utc
        .iter()
        .map(|&m| Epoch::from_jde_utc(m).to_mjd_tt_days())
        .collect()
}

/// Transformation from date in the format YYYY MM DD.FFFFF UTC frame to modified julian date (MJD) TT frame
///
/// Argument
/// --------
/// * `date_str`: a string representing the date in the format YYYY MM DD.FFFFF in the UTC frame
///
/// Return
/// ------
/// * a float representing the input date in modified julian date (MJD) in the TT frame
pub fn frac_date_to_mjd(date_str: &str) -> Result<f64, String> {
    let parts: Vec<&str> = date_str.split_whitespace().collect();
    if parts.len() != 3 {
        return Err("Invalid format, expected: YYYY MM DD.FFFFF".to_string());
    }

    // Extract values
    let year = i32::from_str(parts[0]).map_err(|_| "invalid year")?;
    let month = u8::from_str(parts[1]).map_err(|_| "invalid month")?;
    let day_fraction = f64::from_str(parts[2]).map_err(|_| "invalid frac day")?;

    // Separation of day and fraction day
    let day = day_fraction.trunc() as u8;
    let fraction = day_fraction - day as f64;

    let hour = (fraction * 24.0).trunc() as u8;
    let minute = ((fraction * 24.0 - hour as f64) * 60.0).trunc() as u8;
    let second = (((fraction * 24.0 - hour as f64) * 60.0 - minute as f64) * 60.0) as u8;
    let nano = ((((fraction * 24.0 - hour as f64) * 60.0 - minute as f64) * 60.0 - second as f64)
        * 1e9) as u32;

    // Creation of epoch
    let epoch = Epoch::from_gregorian(year, month, day, hour, minute, second, nano, TimeScale::UTC);

    // Convert to MJD
    Ok(epoch.to_mjd_tt_days())
}

/// Compute the Greenwich Mean Sidereal Time (GMST) in radians
/// for a given Modified Julian Date (UT1 time scale).
///
/// This function implements the IAU 1982/2000 polynomial formula
/// for the mean sidereal time at 0h UT1, plus the fractional-day
/// correction term due to Earth's rotation rate.
///
/// # Arguments
/// * `tjm` - Modified Julian Date (MJD, UT1 time scale)
///
/// # Returns
/// * GMST angle in radians, normalized to the interval [0, 2π).
///
/// # Details
/// The GMST is computed in two steps:
/// 1. Use a cubic polynomial (coefficients C0–C3) to get GMST at 0h UT1
///    in seconds for the given date.
/// 2. Add the contribution of Earth's rotation during the fractional day
///    using the factor `RAP`, which converts solar days to sidereal days.
///
/// # References
/// * IAU 1982, IERS Conventions 1996/2000.
/// * Explanatory Supplement to the Astronomical Almanac (1992).
pub fn gmst(tjm: f64) -> f64 {
    // Polynomial coefficients for GMST at 0h UT1 (in seconds)
    const C0: f64 = 24110.54841;
    const C1: f64 = 8640184.812866;
    const C2: f64 = 9.3104e-2;
    const C3: f64 = -6.2e-6;

    // Ratio of sidereal day to solar day
    const RAP: f64 = 1.00273790934;

    // Extract the integer MJD (0h UT1) and compute centuries since J2000.0
    let itjm = tjm.floor();
    let t = (itjm - T2000) / 36525.0;

    // Step 1: GMST at 0h UT1 using the polynomial expression
    let mut gmst0 = ((C3 * t + C2) * t + C1) * t + C0;

    // Convert GMST from seconds to radians (86400 seconds per day)
    gmst0 *= DPI / 86400.0;

    // Step 2: Add the contribution from the fraction of the day
    // tjm.fract() is the fraction of the current day (0 to 1)
    // Multiplied by 2π (DPI) to convert into radians of a solar day,
    // and then scaled by RAP to account for the faster rotation of sidereal time.
    let h = tjm.fract() * DPI;
    let mut gmst = gmst0 + h * RAP;

    // Normalize GMST to the [0, 2π) range
    let mut i: i64 = (gmst / DPI).floor() as i64;
    if gmst < 0.0 {
        i -= 1;
    }
    gmst -= i as f64 * DPI;

    gmst
}

#[cfg(test)]
mod time_test {
    use super::*;

    #[test]
    fn test_date_to_mjd() {
        let date = vec!["2021-01-01T00:00:00", "2021-01-02T00:00:00"];
        let mjd = date_to_mjd(&date);
        assert_eq!(mjd, vec![59215.0, 59216.0]);
    }

    #[test]
    fn test_mjd_to_jd() {
        let mjd = vec![59215.0, 59216.0];
        let jd = mjd_to_jd(&mjd);
        assert_eq!(jd, vec![2459215.5, 2459216.5]);
    }

    #[test]
    fn test_jd_to_mjd() {
        let jd = vec![2459215.5, 2459216.5];
        let mjd = jd_to_mjd(&jd);
        assert_eq!(mjd, vec![59215.0, 59216.0]);
    }

    #[test]
    fn test_frac_date_to_mjd() {
        let mjd = frac_date_to_mjd("2021 1 1.0").unwrap();
        assert_eq!(mjd, 59215.00080074074);

        let mjd = frac_date_to_mjd("2021 1 1.5").unwrap();
        assert_eq!(mjd, 59215.50080074074);

        let mjd = frac_date_to_mjd("2021 1 1.75").unwrap();
        assert_eq!(mjd, 59215.75080074074);

        let mjd = frac_date_to_mjd("2021 1 1.875").unwrap();
        assert_eq!(mjd, 59215.87580074074);

        let mjd = frac_date_to_mjd("2021 1 1.999").unwrap();
        assert_eq!(mjd, 59215.99980074074);

        let mjd = frac_date_to_mjd("2021 1 1.9999").unwrap();
        assert_eq!(mjd, 59216.00070074073);

        let mjd = frac_date_to_mjd("1976 09 20.93878").unwrap();
        assert_eq!(mjd, 43041.93932611111);
    }

    #[test]
    fn test_gmst() {
        let tut = 57028.478514610404;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.851925725092499);

        let tut = T2000;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.894961212789145);
    }
}
