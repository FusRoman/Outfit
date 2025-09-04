//! NAIF/JPL ephemeris version enumeration and helpers.
//!
//! This module defines the [`NaifVersion`] enum listing a few widely used
//! JPL ephemerides (DE430..DE442), plus utilities to map a version string
//! to the enum and to retrieve the canonical BSP filename.
//!
//! # Examples
//! ```rust, no_run
//! use std::str::FromStr;
//! use outfit::jpl_ephem::naif::naif_version::NaifVersion; // adjust path if needed
//!
//! let v: NaifVersion = "DE440".parse().unwrap();
//! assert_eq!(v.get_filename(), "de440.bsp");
//! ```
//!
//! # Notes
//! * Parsing is **case‑sensitive** and expects exact tokens like `"DE440"`,
//!   `"DE441_part-1"`, `"DE441_part-2"`, etc.
//! * Filenames returned by [`NaifVersion::get_filename`] are the canonical
//!   NAIF BSP kernel names used by the JPL distribution.
//!
//! # See also
//! ------------
//! * `download_jpl_file` component in your crate — where these filenames are used.
//! * `NaifData` loader — consumes a BSP file and exposes interpolation.

use std::str::FromStr;

/// Ephemeris version identifiers supported by this crate.
///
/// See also
/// ------------
/// * `NaifVersion::get_filename` – Returns the canonical BSP filename.
/// * `FromStr` for `NaifVersion` – Parses strings like `"DE440"` or `"DE441_part-1"`.
#[derive(Debug, Clone)]
pub enum NaifVersion {
    DE430,
    DE431p1,
    DE431p2,
    DE432,
    DE435,
    DE438,
    DE440,
    DE440s,
    DE441p1,
    DE441p2,
    DE442,
}

impl NaifVersion {
    /// Return the canonical BSP filename associated with this version.
    ///
    /// Arguments
    /// -----------------
    /// * `&self`: The ephemeris version enum.
    ///
    /// Return
    /// ----------
    /// * `&str` – Canonical NAIF filename (e.g., `"de440.bsp"`).
    ///
    /// See also
    /// ------------
    /// * `FromStr` for `NaifVersion` – Convert `"DE440"` into `NaifVersion::DE440`.
    pub fn get_filename(&self) -> &str {
        match self {
            NaifVersion::DE430 => "de430.bsp",
            NaifVersion::DE431p1 => "de431_part-1.bsp",
            NaifVersion::DE431p2 => "de431_part-2.bsp",
            NaifVersion::DE432 => "de432.bsp",
            NaifVersion::DE435 => "de435.bsp",
            NaifVersion::DE438 => "de438.bsp",
            NaifVersion::DE440 => "de440.bsp",
            NaifVersion::DE440s => "de440s.bsp",
            NaifVersion::DE441p1 => "de441_part-1.bsp",
            NaifVersion::DE441p2 => "de441_part-2.bsp",
            NaifVersion::DE442 => "de442.bsp",
        }
    }

    /// Internal helper to parse a version string into a [`NaifVersion`].
    ///
    /// Arguments
    /// -----------------
    /// * `s`: Input string, e.g. `"DE440"`, `"DE441_part-1"`, `"DE441_part-2"`.
    ///
    /// Return
    /// ----------
    /// * `Option<Self>` – `Some(variant)` if recognized; `None` otherwise.
    ///
    /// See also
    /// ------------
    /// * `FromStr` for `NaifVersion` – Public, error-carrying parsing API.
    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE430" => Some(NaifVersion::DE430),
            "DE431_part-1" => Some(NaifVersion::DE431p1),
            "DE431_part-2" => Some(NaifVersion::DE431p2),
            "DE432" => Some(NaifVersion::DE432),
            "DE435" => Some(NaifVersion::DE435),
            "DE438" => Some(NaifVersion::DE438),
            "DE440" => Some(NaifVersion::DE440),
            "DE440s" => Some(NaifVersion::DE440s),
            "DE441_part-1" => Some(NaifVersion::DE441p1),
            "DE441_part-2" => Some(NaifVersion::DE441p2),
            "DE442" => Some(NaifVersion::DE442),
            _ => None,
        }
    }
}

impl FromStr for NaifVersion {
    type Err = String;

    /// Parse a `NaifVersion` from its textual identifier.
    ///
    /// Arguments
    /// -----------------
    /// * `s`: Version string, e.g. `"DE440"`, `"DE441_part-1"`.
    ///
    /// Return
    /// ----------
    /// * `Result<NaifVersion, String>` – The parsed variant, or an error message
    ///   of the form `"Invalid NAIF version: <input>"`.
    ///
    /// See also
    /// ------------
    /// * [`NaifVersion::get_filename`] – Map the version to its BSP filename.
    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        NaifVersion::from_str(s).ok_or_else(|| format!("Invalid NAIF version: {s}"))
    }
}
