//! # MPC 80-Column Observation Reader
//!
//! Utilities to parse **MPC 80-column** astrometric observations and turn them into
//! [`Observation`] values usable by the orbit-determination pipeline.
//!
//! ## Overview
//! -----------------
//! This module provides:
//! - A small error type [`ParseObsError`] describing MPC parsing failures.
//! - A crate-internal line parser (`from_80col`) that converts a single 80-col line
//!   into an [`Observation`] with angles in **radians** and time in **MJD (TT)**.
//! - A crate-visible batch routine \[`extract_80col`\] that reads an entire file,
//!   returns all parsed [`Observation`]s, and extracts the **object number** from
//!   the header line.
//!
//! The implementation enforces **MPC CCD** observations (rejects non-CCD lines) and
//! converts RA/Dec strings using robust helpers (`parse_ra_to_deg`, `parse_dec_to_deg`),
//! then applies uncertainty handling that accounts for `cos δ` on the RA component.
//!
//! ## Units & Conventions
//! -----------------
//! - **Input format:** MPC 80-column fixed-width ASCII lines.
//! - **Angles:** RA/Dec are parsed from strings into **degrees**, then converted to
//!   **radians**.
//! - **Uncertainties:** parsed in input units (hour/deg for RA/Dec) and mapped to
//!   **radians**; RA uncertainties are divided by `cos δ`.
//! - **Time scale:** Observing time is parsed from fractional date and converted
//!   to **MJD (TT)** via [`frac_date_to_mjd`].
//! - **Observer site:** Extracted from columns 77–80 (MPC code), mapped to the
//!   compact `u16` site id via \[`Outfit::uint16_from_mpc_code`\].
//!
//! ## File-Level Object Number
//! -----------------
//! \[`extract_80col`\] retrieves an **object number** from the first line using the
//! conventional MPC header locations (columns `0..5` or fallback to `5..12`), trims
//! leading zeros, and returns it as an [`ObjectNumber::String`].
//!
//! ## Error Handling
//! -----------------
//! Parser failures are wrapped into [`OutfitError::Parsing80ColumnFileError`] with a
//! [`ParseObsError`] payload for precise diagnostics (e.g. *line too short*,
//! *invalid RA*, *invalid date*). Non-CCD lines are filtered out by default.
//!
//! ## See also
//! ------------
//! * [`parse_ra_to_deg`] – RA string → degrees with uncertainty.
//! * [`parse_dec_to_deg`] – Dec string → degrees with uncertainty.
//! * [`frac_date_to_mjd`] – Fractional date → MJD (TT).
//! * [`Observation`] – Internal astrometric sample type.
//! * [`ObjectNumber`] – Logical object identifier.
//! * [`Outfit`] – Global state (site registry, time scales, etc.).
use std::ops::Range;

use camino::Utf8Path;
use thiserror::Error;

use crate::{
    constants::Observations,
    conversion::{parse_dec_to_deg, parse_ra_to_deg},
    observations::Observation,
    time::frac_date_to_mjd,
    ObjectNumber, Outfit, OutfitError, RADH, RADSEC,
};

/// Line-level parsing errors for MPC 80-column observations.
///
/// Variants
/// -----------------
/// * `TooShortLine` – The line does not reach 80 characters.
/// * `NotCCDObs` – The line is not flagged as a CCD observation (column 15 is `'s'`).
/// * `InvalidRA` – Failed to parse RA field (`line[32..44]`); payload carries the offending slice.
/// * `InvalidDec` – Failed to parse Dec field (`line[44..56]`); payload carries the offending slice.
/// * `InvalidDate` – Failed to parse fractional date (`line[15..32]`); payload carries the offending slice.
///
/// See also
/// ------------
/// * [`parse_ra_to_deg`] – Robust RA parser with uncertainty extraction.
/// * [`parse_dec_to_deg`] – Robust Dec parser with uncertainty extraction.
/// * [`frac_date_to_mjd`] – Converts fractional date to MJD (TT).
#[derive(Error, Debug, PartialEq)]
pub enum ParseObsError {
    #[error("The line is too short")]
    TooShortLine,
    #[error("The line is not a CCD observation")]
    NotCCDObs,
    #[error("Error parsing RA: {0}")]
    InvalidRA(String),
    #[error("Invalid Dec value: {0}")]
    InvalidDec(String),
    #[error("Invalid date: {0}")]
    InvalidDate(String),
}

/// Parse a single **MPC 80-column** line into an [`Observation`] (crate-private helper).
///
/// This routine enforces **CCD** observations, parses RA/Dec/time fields, maps the
/// MPC site code to the compact observer id, and builds an [`Observation`] in **radians**
/// with a **MJD (TT)** epoch. RA uncertainties are divided by `cos δ` to reflect
/// the geometry on the sphere.
///
/// Arguments
/// -----------------
/// * `env_state` – Global state used to resolve MPC site codes and build observations.
/// * `line` – A single 80-column ASCII line.
///
/// Return
/// ----------
/// * A parsed [`Observation`] or an [`OutfitError::Parsing80ColumnFileError`] on failure.
///
/// Panics
/// ----------
/// * Never panics directly; errors are surfaced as [`OutfitError`]. Bounds use fixed slices.
///
/// Field Layout (MPC 80-col subset used here)
/// -----------------
/// * `15..32` – Fractional date (UTC-like string expected by `frac_date_to_mjd`).  
/// * `32..44` – Right ascension string (parsed by `parse_ra_to_deg`).  
/// * `44..56` – Declination string (parsed by `parse_dec_to_deg`).  
/// * `77..80` – MPC site code.
///
/// See also
/// ------------
/// * [`parse_ra_to_deg`] – RA parsing (degrees + uncertainty).
/// * [`parse_dec_to_deg`] – Dec parsing (degrees + uncertainty).
/// * [`frac_date_to_mjd`] – Fractional date → MJD (TT).
fn from_80col(env_state: &mut Outfit, line: &str) -> Result<Observation, OutfitError> {
    if line.len() < 80 {
        return Err(OutfitError::Parsing80ColumnFileError(
            ParseObsError::TooShortLine,
        ));
    }

    if line.chars().nth(14) == Some('s') {
        return Err(OutfitError::Parsing80ColumnFileError(
            ParseObsError::NotCCDObs,
        ));
    }

    let (ra, error_ra) = parse_ra_to_deg(line[32..44].trim()).ok_or_else(|| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidRA(
            line[32..44].trim().to_string(),
        ))
    })?;

    let (dec, error_dec) = parse_dec_to_deg(line[44..56].trim()).ok_or_else(|| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidDec(
            line[44..56].trim().to_string(),
        ))
    })?;

    let time = frac_date_to_mjd(line[15..32].trim()).map_err(|_| {
        OutfitError::Parsing80ColumnFileError(ParseObsError::InvalidDate(
            line[15..32].trim().to_string(),
        ))
    })?;

    let observer_id = env_state.uint16_from_mpc_code(&line[77..80].trim().into());
    let observer = env_state.get_observer_from_uint16(observer_id);

    let max_rms = |observation_error: f64, observer_error: f64, factor: f64| {
        f64::max(observation_error, observer_error * factor)
    };

    let dec_radians = dec.to_radians();
    let dec_rad_cos = dec_radians.cos();

    let ra_error = max_rms(
        (error_ra * RADH) / dec_rad_cos,
        observer.ra_accuracy.map(|v| v.into_inner()).unwrap_or(0.0),
        RADSEC / dec_rad_cos,
    );

    let dec_error = max_rms(
        error_dec.to_radians(),
        observer.dec_accuracy.map(|v| v.into_inner()).unwrap_or(0.0),
        RADSEC,
    );

    let observation = Observation::new(
        env_state,
        observer_id,
        ra.to_radians(),
        ra_error,
        dec_radians,
        dec_error,
        time,
    )?;
    Ok(observation)
}

/// Read a full **MPC 80-column** file, returning parsed observations and the object number.
///
/// Lines that are not CCD observations are **silently skipped**. Any other parsing error
/// triggers a panic with context (the current strategy is fail-fast for corrupted inputs).
/// The object number is extracted from the first line (`0..5`, fallback `5..12`), trimming
/// leading zeros.
///
/// Arguments
/// -----------------
/// * `env_state` – Global state used to resolve MPC site codes and build observations.
/// * `colfile` – Path to the MPC 80-column file.
///
/// Return
/// ----------
/// * A tuple `(observations, object_number)` where:
///   - `observations` is a `Vec<Observation>` with angles in **radians** and epochs in **MJD (TT)**,
///   - `object_number` is the header-derived [`ObjectNumber::String`].
///
/// Panics
/// ----------
/// * Panics if the file cannot be read or if a non-CCD line fails to parse for reasons
///   other than the expected `NotCCDObs` (fail-fast behavior).
///
/// See also
/// ------------
/// * [`Observation`] – Parsed astrometric sample.
/// * [`ObjectNumber`] – Identifier extracted from header columns.
/// * [`parse_ra_to_deg`], [`parse_dec_to_deg`], [`frac_date_to_mjd`] – Parsing helpers.
pub(crate) fn extract_80col(
    env_state: &mut Outfit,
    colfile: &Utf8Path,
) -> Result<(Observations, ObjectNumber), OutfitError> {
    let file_content = std::fs::read_to_string(colfile)
        .unwrap_or_else(|_| panic!("Could not read file {}", colfile.as_str()));

    let first_line = file_content
        .lines()
        .next()
        .unwrap_or_else(|| panic!("Could not read first line of file {}", colfile.as_str()));

    fn get_object_number(line: &str, range: Range<usize>) -> String {
        line[range].trim_start_matches('0').trim().to_string()
    }

    let mut object_number = get_object_number(first_line, 0..5);
    if object_number.is_empty() {
        object_number = get_object_number(first_line, 5..12);
    }

    Ok((
        file_content
            .lines()
            .filter_map(|line| match from_80col(env_state, line) {
                Ok(obs) => Some(obs),
                Err(OutfitError::Parsing80ColumnFileError(ParseObsError::NotCCDObs)) => None,
                Err(e) => panic!("Error parsing line: {e:?}"),
            })
            .collect(),
        ObjectNumber::String(object_number),
    ))
}

#[cfg(test)]
#[cfg(feature = "jpl-download")]
mod mpc_80col_test {
    use super::*;

    #[test]
    fn test_from_80col_valid_line() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C2009 09 15.23433 22 52 22.62 -14 47 03.2          20.8 Vr~097wG96";

        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(result.is_ok());
        let obs = result.unwrap();

        assert_eq!(
            obs,
            Observation {
                observer: 0,
                ra: 5.988124307160555,
                error_ra: 1.2535340843609459e-6,
                dec: -0.25803335512429054,
                error_dec: 1.0181086985431635e-6,
                time: 55089.23509601851,
                observer_earth_position: [
                    3.0499942822953885e-5,
                    -8.594304778250371e-6,
                    2.8491013919142154e-5
                ]
                .into(),
                observer_helio_position: [
                    0.9968138444702415,
                    -0.12221921296802639,
                    -0.05295724448160355
                ]
                .into(),
            }
        );
    }

    #[test]
    fn test_from_80col_too_short_line() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line = "short line";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::TooShortLine
            ))
        ));
    }

    #[test]
    fn test_from_80col_invalid_date() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C20xx 09 15.23433 22 52 22.62 -14 47 03.2          20.8 Vr~097wG96";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::InvalidDate(_)
            ))
        ));
    }

    #[test]
    fn test_from_80col_invalid_ra_dec() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let line =
            "     K09R05F  C2009 09 15.23433 XX YY ZZ.ZZ -AA BB CC.C          20.8 Vr~097wG96";
        let mut env_state = OUTFIT_HORIZON_TEST.0.clone();
        let result = from_80col(&mut env_state, line);

        assert!(matches!(
            result,
            Err(OutfitError::Parsing80ColumnFileError(
                ParseObsError::InvalidRA(_)
            ))
        ));
    }
}
