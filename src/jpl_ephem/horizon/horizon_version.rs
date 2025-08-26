//! Version helpers for JPL Horizons DE ephemerides.
//!
//! This module defines [`JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion), an enum covering commonly used
//! JPL DE solutions (e.g., DE430, DE440, DE441) and provides utilities to map
//! those labels to:
//! - the **legacy Horizons Linux distribution path fragments**
//!   (as served under `.../eph/planets/Linux/`), via
//!   [`JPLHorizonVersion::get_filename`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion::get_filename),
//! - the **canonical NAIF SPK kernel filenames** (`*.bsp`), via
//!   [`JPLHorizonVersion::to_filename`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion::to_filename) and
//!   [`JPLHorizonVersion::from_filename`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion::from_filename).
//!
//! Typical use
//! -----------------
//! ```rust
//! use std::str::FromStr;
//! use outfit::jpl_ephem::horizon::horizon_version::JPLHorizonVersion;
//!
//! let v = JPLHorizonVersion::from_str("DE440").unwrap();
//! assert_eq!(v.get_filename(), "de440/linux_p1550p2650.440");
//! assert_eq!(v.to_filename(), "DE440.bsp");
//! assert_eq!(JPLHorizonVersion::from_filename("DE440.bsp"), Some(JPLHorizonVersion::DE440));
//! ```
//!
//! See also
//! ------------
//! * [`crate::jpl_ephem::horizon::horizon_data::HorizonData`] – parsing & interpolation entry point.
//! * JPL Horizons FTP layout (legacy `planets/Linux/`) if you need the base URL handling.
use std::str::FromStr;

/// Enumerates supported JPL DE ephemeris solutions.
///
/// Variants correspond to public JPL solutions distributed through Horizons
/// and/or NAIF. Some have a `t` suffix (e.g. `DE430t`, `DE440t`), which
/// denotes a **trimmed/variant** distribution provided by JPL; choose the one
/// required by your pipeline or data source.
///
/// See also
/// ------------
/// * [`JPLHorizonVersion::get_filename`] – legacy Horizons path fragment.
/// * [`JPLHorizonVersion::to_filename`] – NAIF SPK kernel filename.
/// * [`JPLHorizonVersion::from_filename`] – parse an SPK filename back to a version.
#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum JPLHorizonVersion {
    DE102,
    DE200,
    DE202,
    DE403,
    DE405,
    DE406,
    DE410,
    DE413,
    DE414,
    DE418,
    DE421,
    DE422,
    DE423,
    DE430,
    DE430t,
    DE431,
    DE440,
    DE440t,
    DE441,
}

impl JPLHorizonVersion {
    /// Return the **legacy Horizons path fragment** for this DE version.
    ///
    /// This is the relative path used under the Horizons Linux directory
    /// tree (e.g., `de440/linux_p1550p2650.440`). It is not a full URL or
    /// filesystem path; prepend your base location accordingly.
    ///
    /// Return
    /// ----------
    /// * A `&str` path fragment suitable for building download or lookup paths.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// use outfit::jpl_ephem::horizon::horizon_version::JPLHorizonVersion;
    /// assert_eq!(JPLHorizonVersion::DE430.get_filename(), "de430/linux_p1550p2650.430");
    /// ```
    ///
    /// See also
    /// ------------
    /// * [`JPLHorizonVersion::to_filename`] – SPK filename if you are using NAIF kernels.
    pub fn get_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "de102/lnxm1410p3002.102",
            JPLHorizonVersion::DE200 => "de200/lnxm1600p2170.200",
            JPLHorizonVersion::DE202 => "de202/lnxp1900p2050.202",
            JPLHorizonVersion::DE403 => "de403/lnxp1600p2200.403",
            JPLHorizonVersion::DE405 => "de405/lnxp1600p2200.405",
            JPLHorizonVersion::DE406 => "de406/lnxm3000p3000.406",
            JPLHorizonVersion::DE410 => "de410/lnxp1960p2020.410",
            JPLHorizonVersion::DE413 => "de413/lnxp1900p2050.413",
            JPLHorizonVersion::DE414 => "de414/lnxp1600p2200.414",
            JPLHorizonVersion::DE418 => "de418/lnxp1900p2050.418",
            JPLHorizonVersion::DE421 => "de421/lnxp1900p2053.421",
            JPLHorizonVersion::DE422 => "de422/lnxm3000p3000.422",
            JPLHorizonVersion::DE423 => "de423/lnxp1800p2200.423",
            JPLHorizonVersion::DE430 => "de430/linux_p1550p2650.430",
            JPLHorizonVersion::DE430t => "de430t/linux_p1550p2650.430t",
            JPLHorizonVersion::DE431 => "de431/lnxm13000p17000.431",
            JPLHorizonVersion::DE440 => "de440/linux_p1550p2650.440",
            JPLHorizonVersion::DE440t => "de440t/linux_p1550p2650.440t",
            JPLHorizonVersion::DE441 => "de441/linux_m13000p17000.441",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE102" => Some(JPLHorizonVersion::DE102),
            "DE200" => Some(JPLHorizonVersion::DE200),
            "DE202" => Some(JPLHorizonVersion::DE202),
            "DE403" => Some(JPLHorizonVersion::DE403),
            "DE405" => Some(JPLHorizonVersion::DE405),
            "DE406" => Some(JPLHorizonVersion::DE406),
            "DE410" => Some(JPLHorizonVersion::DE410),
            "DE413" => Some(JPLHorizonVersion::DE413),
            "DE414" => Some(JPLHorizonVersion::DE414),
            "DE418" => Some(JPLHorizonVersion::DE418),
            "DE421" => Some(JPLHorizonVersion::DE421),
            "DE422" => Some(JPLHorizonVersion::DE422),
            "DE423" => Some(JPLHorizonVersion::DE423),
            "DE430" => Some(JPLHorizonVersion::DE430),
            "DE430t" => Some(JPLHorizonVersion::DE430t),
            "DE431" => Some(JPLHorizonVersion::DE431),
            "DE440" => Some(JPLHorizonVersion::DE440),
            "DE440t" => Some(JPLHorizonVersion::DE440t),
            "DE441" => Some(JPLHorizonVersion::DE441),
            _ => None,
        }
    }

    /// Return the canonical **NAIF SPK kernel filename** for this DE version.
    ///
    /// This is the standard filename (e.g., `DE440.bsp`) as used by NAIF
    /// SPICE kernels. It is useful when integrating with an SPK/DAF backend.
    ///
    /// Return
    /// ----------
    /// * A `&str` filename ending in `.bsp`.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// use outfit::jpl_ephem::horizon::horizon_version::JPLHorizonVersion;
    /// assert_eq!(JPLHorizonVersion::DE441.to_filename(), "DE441.bsp");
    /// ```
    ///
    /// See also
    /// ------------
    /// * [`JPLHorizonVersion::from_filename`] – parse back from an SPK filename.
    pub fn to_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "DE102.bsp",
            JPLHorizonVersion::DE200 => "DE200.bsp",
            JPLHorizonVersion::DE202 => "DE202.bsp",
            JPLHorizonVersion::DE403 => "DE403.bsp",
            JPLHorizonVersion::DE405 => "DE405.bsp",
            JPLHorizonVersion::DE406 => "DE406.bsp",
            JPLHorizonVersion::DE410 => "DE410.bsp",
            JPLHorizonVersion::DE413 => "DE413.bsp",
            JPLHorizonVersion::DE414 => "DE414.bsp",
            JPLHorizonVersion::DE418 => "DE418.bsp",
            JPLHorizonVersion::DE421 => "DE421.bsp",
            JPLHorizonVersion::DE422 => "DE422.bsp",
            JPLHorizonVersion::DE423 => "DE423.bsp",
            JPLHorizonVersion::DE430 => "DE430.bsp",
            JPLHorizonVersion::DE430t => "DE430t.bsp",
            JPLHorizonVersion::DE431 => "DE431.bsp",
            JPLHorizonVersion::DE440 => "DE440.bsp",
            JPLHorizonVersion::DE440t => "DE440t.bsp",
            JPLHorizonVersion::DE441 => "DE441.bsp",
        }
    }

    /// Parse a `JPLHorizonVersion` from a canonical SPK filename.
    ///
    /// Arguments
    /// -----------------
    /// * `filename` — A kernel filename like `DE440.bsp`.
    ///
    /// Return
    /// ----------
    /// * `Some(JPLHorizonVersion)` if the filename matches a known kernel;
    ///   otherwise `None`.
    ///
    /// Examples
    /// ----------
    /// ```rust, no_run
    /// use outfit::jpl_ephem::horizon::horizon_version::JPLHorizonVersion;
    /// assert_eq!(JPLHorizonVersion::from_filename("DE430t.bsp"), Some(JPLHorizonVersion::DE430t));
    /// assert!(JPLHorizonVersion::from_filename("UNKNOWN.bsp").is_none());
    /// ```
    pub fn from_filename(filename: &str) -> Option<Self> {
        match filename {
            "DE102.bsp" => Some(JPLHorizonVersion::DE102),
            "DE200.bsp" => Some(JPLHorizonVersion::DE200),
            "DE202.bsp" => Some(JPLHorizonVersion::DE202),
            "DE403.bsp" => Some(JPLHorizonVersion::DE403),
            "DE405.bsp" => Some(JPLHorizonVersion::DE405),
            "DE406.bsp" => Some(JPLHorizonVersion::DE406),
            "DE410.bsp" => Some(JPLHorizonVersion::DE410),
            "DE413.bsp" => Some(JPLHorizonVersion::DE413),
            "DE414.bsp" => Some(JPLHorizonVersion::DE414),
            "DE418.bsp" => Some(JPLHorizonVersion::DE418),
            "DE421.bsp" => Some(JPLHorizonVersion::DE421),
            "DE422.bsp" => Some(JPLHorizonVersion::DE422),
            "DE423.bsp" => Some(JPLHorizonVersion::DE423),
            "DE430.bsp" => Some(JPLHorizonVersion::DE430),
            "DE430t.bsp" => Some(JPLHorizonVersion::DE430t),
            "DE431.bsp" => Some(JPLHorizonVersion::DE431),
            "DE440.bsp" => Some(JPLHorizonVersion::DE440),
            "DE440t.bsp" => Some(JPLHorizonVersion::DE440t),
            "DE441.bsp" => Some(JPLHorizonVersion::DE441),
            _ => None,
        }
    }
}

/// Parse a [`JPLHorizonVersion`] from its textual label (e.g., `"DE440"`).
///
/// Implements standard `FromStr` so you can do:
/// ```rust, no_run
/// use std::str::FromStr;
/// use outfit::jpl_ephem::horizon::horizon_version::JPLHorizonVersion;
/// let v = JPLHorizonVersion::from_str("DE431").unwrap();
/// assert_eq!(v.to_filename(), "DE431.bsp");
/// ```
///
/// Errors
/// ----------
/// * Returns an error string if the label does not match a known version.
///
/// See also
/// ------------
/// * [`JPLHorizonVersion::from_filename`] – parse from an SPK filename.
/// * [`JPLHorizonVersion::get_filename`] – legacy Horizons path fragment.
impl FromStr for JPLHorizonVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        JPLHorizonVersion::from_str(s).ok_or_else(|| format!("Invalid JPL Horizon version: {s}"))
    }
}
