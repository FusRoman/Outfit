//! # Tabular display for astrometric observations
//!
//! This module provides zero-allocation, table-style renderers for an
//! [`Observations`] collection (a `SmallVec<[Observation; 6]>`).
//!
//! ## Overview
//!
//! The main entry point is the **display adaptor** [`ObservationsDisplay`], which borrows
//! the underlying observation slice and prints a formatted table when used with
//! Rust’s formatting machinery (`{}`).
//!
//! Three **layouts** are supported:
//!
//! - **Default** (compact): `# | Site | MJD (TT) | RA[hms] ±σ["] | DEC[dms] ±σ["]`
//! - **Wide** (diagnostic): adds `JD(TT) | RA[rad] | DEC[rad] | |r_geo| AU | |r_hel| AU`
//! - **ISO** (timestamp-centric): `ISO (TT)` and `ISO (UTC)` replace the MJD/JD columns
//!
//! Precision controls:
//!
//! - `sec_prec` — fractional digits for sexagesimal seconds **and** ISO seconds
//! - `dist_prec` — fixed-point digits for AU distances (wide mode)
//!
//! ## Units & conventions
//!
//! - **Time**: MJD and JD are **TT** in the tables (ISO mode prints **TT** and **UTC**).
//! - **Angles**: RA/DEC are rendered in sexagesimal (RA in **hours**, DEC in **degrees**).
//! - **Uncertainties**: shown in **arcseconds**; converted from radians using [`RAD2ARC`].
//! - **Positions**: vector norms are in **AU** (equatorial mean J2000 frame by convention).
//!
//! ## Examples
//!
//! ```rust,ignore
//! use outfit::observations::display::{ObservationsDisplayExt};
//!
//! // Compact table
//! println!("{}", observations.show());
//!
//! // Wide table (adds JD, radians, |r| in AU)
//! println!("{}", observations.show().wide(true));
//!
//! // ISO table (ISO TT + ISO UTC)
//! println!("{}", observations.table_iso());
//!
//! // Customize precisions
//! println!("{}", observations.table_iso().with_seconds_precision(4));
//! println!("{}", observations.table_wide().with_distance_precision(7));
//!
//! // Get an owned string (compact mode)
//! let s = observations.show_string();
//! ```
//!
//! ## See also
//!
//! - [`Observation`](crate::observations::Observation) — single-observation pretty-printer and helpers
//! - [`crate::conversion::ra_hms_prec`] / [`crate::conversion::dec_sdms_prec`] — sexagesimal decomposition with carry
//! - [`crate::time::fmt_ss`] — seconds string `"SS.sss"` with 2-digit integer part
//! - [`crate::time::iso_tt_from_epoch`] / [`crate::time::iso_utc_from_epoch`] — ISO renderers via `hifitime`
//!
//! ## Notes
//!
//! - The adaptor **does not allocate** for the table itself; only transient
//!   small strings are created per row.
//! - Header widths are chosen for readability; they are not intended to be parsed.
//! - ISO rendering relies on `hifitime`’s leap-second tables for UTC.

use std::fmt;

use hifitime::{Epoch, TimeScale};

use crate::constants::{JDTOMJD, RAD2ARC};
use crate::conversion::{dec_sdms_prec, ra_hms_prec};
use crate::time::{fmt_ss, iso_tt_from_epoch, iso_utc_from_epoch};
use crate::Observations;

/// Internal layout selector for the table renderer.
///
/// This enum is crate-internal and selects the columns printed by
/// [`ObservationsDisplay`]. It is not exposed publicly.
///
/// Variants
/// -----------------
/// * `Default` — Compact columns: MJD(TT), RA±σ, DEC±σ.
/// * `Wide` — Adds JD(TT), RA/DEC in radians, and AU vector norms.
/// * `Iso` — Replaces MJD/JD with `ISO (TT)` and `ISO (UTC)`.
enum TableMode {
    Default, // MJD(TT) + RA/DEC ±σ
    Wide,    // + JD(TT), RA/DEC [rad], |r_geo|, |r_hel|
    Iso,     // ISO TT + ISO UTC instead of MJD/JD
}

/// Display adaptor to render an [`Observations`] collection as a **table**.
///
/// Render modes
/// -----------------
/// * **Default** (via [`ObservationsDisplayExt::show`]):  
///   columns `# | Site | MJD (TT) | RA[hms] ±σ["] | DEC[dms] ±σ["]`.
/// * **Wide** (via [`ObservationsDisplayExt::table_wide`]):  
///   adds `JD (TT) | RA [rad] | DEC [rad] | |r_geo| AU | |r_hel| AU`.
/// * **ISO** (via [`ObservationsDisplayExt::table_iso`]):  
///   replaces MJD/JD with `ISO (TT)` and `ISO (UTC)` timestamps.
///
/// Precision
/// -----------------
/// * `sec_prec` controls the number of fractional digits for sexagesimal seconds
///   **and** ISO seconds.
/// * `dist_prec` controls fixed-point digits for AU distances (wide mode).
///
/// Notes
/// ----------
/// * This adaptor **borrows** the observation slice; it does not own or copy it.
/// * RA/DEC are decomposed with carry-safe helpers (see *See also*).
/// * Uncertainties are printed in **arcseconds**, using [`RAD2ARC`] for conversion.
///
/// See also
/// ------------
/// * [`ObservationsDisplayExt`] – Ergonomic builders for each mode.
/// * [`Observation`](crate::observations::Observation) – Per-row semantics used to derive the columns.
/// * [`crate::conversion::ra_hms_prec`] / [`crate::conversion::dec_sdms_prec`].
pub struct ObservationsDisplay<'a> {
    /// Borrowed collection to render. No allocation or copying occurs.
    obs: &'a Observations,
    /// Column layout selector (default / wide / iso).
    mode: TableMode,
    /// Fractional digits for sexagesimal and ISO seconds (default = 3).
    sec_prec: usize,
    /// Fixed-point digits for AU distances in wide mode (default = 6).
    dist_prec: usize,
}

impl<'a> ObservationsDisplay<'a> {
    /// Build a new table adaptor (default: **compact** columns).
    ///
    /// Arguments
    /// -----------------
    /// * `obs` – Borrowed observation container to render.
    ///
    /// Return
    /// ----------
    /// * An `ObservationsDisplay` configured for the **Default** mode.
    ///
    /// See also
    /// ------------
    /// * [`Self::wide`] – Enable the wide layout.
    /// * [`Self::iso`] – Enable the ISO layout.
    pub fn new(obs: &'a Observations) -> Self {
        Self {
            obs,
            mode: TableMode::Default,
            sec_prec: 3,
            dist_prec: 6,
        }
    }

    /// Switch to **wide** mode (adds JD, radians, vector norms).
    ///
    /// Arguments
    /// -----------------
    /// * `yes` – If `true`, selects the **Wide** layout; otherwise resets to **Default**.
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style), allowing chained configuration.
    ///
    /// See also
    /// ------------
    /// * [`ObservationsDisplayExt::table_wide`]
    pub fn wide(mut self, yes: bool) -> Self {
        self.mode = if yes {
            TableMode::Wide
        } else {
            TableMode::Default
        };
        self
    }

    /// Switch to **ISO** mode (replace MJD/JD with ISO TT and ISO UTC columns).
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style), allowing chained configuration.
    ///
    /// Notes
    /// ----------
    /// * ISO strings are produced via `hifitime`, using leap-second tables for UTC.
    ///
    /// See also
    /// ------------
    /// * [`ObservationsDisplayExt::table_iso`]
    pub fn iso(mut self) -> Self {
        self.mode = TableMode::Iso;
        self
    }

    /// Set seconds precision for **sexagesimal** and **ISO** seconds.
    ///
    /// Arguments
    /// -----------------
    /// * `p` – Number of fractional digits to render (typ. `0..=9`).
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style).
    ///
    /// Notes
    /// ----------
    /// * Affects both RA/DEC seconds and ISO seconds.
    pub fn with_seconds_precision(mut self, p: usize) -> Self {
        self.sec_prec = p;
        self
    }

    /// Set decimal precision for **AU** distances (wide mode only).
    ///
    /// Arguments
    /// -----------------
    /// * `p` – Fixed-point fractional digits for `|r_geo|` and `|r_hel|`.
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style).
    pub fn with_distance_precision(mut self, p: usize) -> Self {
        self.dist_prec = p;
        self
    }
}

/// Ergonomic extension to create table adaptors from an [`Observations`] collection.
///
/// Provided builders
/// -----------------
/// * [`ObservationsDisplayExt::show`] – Default compact table.
/// * [`ObservationsDisplayExt::table_wide`] – Wide table with diagnostics.
/// * [`ObservationsDisplayExt::table_iso`] – ISO-centric table (TT + UTC).
///
/// Examples
/// ----------
/// ```rust,ignore
/// println!("{}", observations.show());          // Default
/// println!("{}", observations.table_wide());    // Wide
/// println!("{}", observations.table_iso());     // ISO
/// println!("{}", observations.show().with_seconds_precision(4));
/// ```
pub trait ObservationsDisplayExt {
    /// Wide table (adds JD, RA/DEC in radians, vector norms in AU).
    ///
    /// Return
    /// ----------
    /// * A configured [`ObservationsDisplay`] in **Wide** mode.
    fn table_wide(&self) -> ObservationsDisplay<'_>;

    /// ISO table (replaces MJD/JD with `ISO (TT)` and `ISO (UTC)`).
    ///
    /// Return
    /// ----------
    /// * A configured [`ObservationsDisplay`] in **ISO** mode.
    fn table_iso(&self) -> ObservationsDisplay<'_>;

    /// Create a zero-allocation display adaptor (Default/compact mode).
    ///
    /// Examples
    /// ----------
    /// ```rust,ignore
    /// // Compact table
    /// println!("{}", observations.show());
    ///
    /// // Derive other modes from it (builder style)
    /// println!("{}", observations.show().wide(true)); // Wide
    /// println!("{}", observations.show().iso());      // ISO
    /// ```
    fn show(&self) -> ObservationsDisplay<'_>;

    /// Convenience: return a formatted `String` in **compact** mode.
    ///
    /// Return
    /// ----------
    /// * An owned `String` containing the compact table.
    fn show_string(&self) -> String {
        format!("{}", self.show())
    }
}

impl ObservationsDisplayExt for Observations {
    fn table_wide(&self) -> ObservationsDisplay<'_> {
        ObservationsDisplay::new(self).wide(true)
    }
    fn table_iso(&self) -> ObservationsDisplay<'_> {
        ObservationsDisplay::new(self).iso()
    }
    fn show(&self) -> ObservationsDisplay<'_> {
        ObservationsDisplay::new(self)
    }
}

impl fmt::Display for ObservationsDisplay<'_> {
    /// Render the table according to the selected mode.
    ///
    /// Notes
    /// ----------
    /// * Rows are printed in the **current order** of the underlying slice.
    /// * Numeric fields are right-aligned; headers have fixed widths for readability.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.obs.len();
        writeln!(f, "Observations (n={n})")?;
        writeln!(f, "-------------------")?;

        // Headers per mode
        match self.mode {
            TableMode::Default => {
                writeln!(
                    f,
                    "{:>3}  {:>5}  {:>14}  {:>20}  {:>20}",
                    "#", "Site", "MJD (TT)", "RA ±σ[arcsec]", "DEC ±σ[arcsec]"
                )?;
            }
            TableMode::Wide => {
                writeln!(
                    f,
                    "{:>3}  {:>5}  {:>14}  {:>14}  {:>20}  {:>11}  {:>20}  {:>11}  {:>12}  {:>12}",
                    "#",
                    "Site",
                    "MJD (TT)",
                    "JD (TT)",
                    "RA ±σ[arcsec]",
                    "RA [rad]",
                    "DEC ±σ[arcsec]",
                    "DEC [rad]",
                    "|r_geo| AU",
                    "|r_hel| AU"
                )?;
            }
            TableMode::Iso => {
                writeln!(
                    f,
                    "{:>3}  {:>5}  {:>26}  {:>26}  {:>20}  {:>20}",
                    "#", "Site", "ISO (TT)", "ISO (UTC)", "RA ±σ[arcsec]", "DEC ±σ[arcsec]"
                )?;
            }
        }

        let sp = self.sec_prec;
        let dp = self.dist_prec;

        for (i, o) in self.obs.iter().enumerate() {
            let site = o.observer;

            // Angles & uncertainties
            let ra_rad: f64 = o.ra;
            let dec_rad: f64 = o.dec;
            let (hh, mm, ss) = ra_hms_prec(ra_rad, sp);
            let (sgn, dd, dm, ds) = dec_sdms_prec(dec_rad, sp);
            let ss_s = fmt_ss(ss, sp);
            let ds_s = fmt_ss(ds, sp);
            let sra_as: f64 = o.error_ra * RAD2ARC;
            let sdec_as: f64 = o.error_dec * RAD2ARC;
            let ra_str = format!("{hh:02}h{mm:02}m{ss_s}s ± {sra_as:.3}\"");
            let dec_str = format!("{sgn}{dd:02}°{dm:02}'{ds_s}\" ± {sdec_as:.3}\"");

            match self.mode {
                TableMode::Default => {
                    let mjd_tt: f64 = o.time;
                    writeln!(
                        f,
                        "{i:>3}  {site:>5}  {mjd_tt:>14.6}  {ra_str:>20}  {dec_str:>20}"
                    )?;
                }
                TableMode::Wide => {
                    let mjd_tt: f64 = o.time;
                    let jd_tt = mjd_tt + JDTOMJD;
                    let r_geo = o.observer_earth_position.norm();
                    let r_hel = o.observer_helio_position.norm();
                    writeln!(
                        f,
                        "{i:>3}  {site:>5}  {mjd_tt:>14.6}  {jd_tt:>14.6}  {ra_str:>20}  {ra_rad:>11.7}  {dec_str:>20}  {dec_rad:>11.7}  {r_geo:>12.dp$}  {r_hel:>12.dp$}"
                    )?;
                }
                TableMode::Iso => {
                    // Build ISO strings from MJD(TT) via hifitime helpers.
                    let mjd_tt: f64 = o.time;
                    let epoch_tt = Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT);
                    let iso_tt = iso_tt_from_epoch(epoch_tt, sp);
                    let iso_utc = iso_utc_from_epoch(epoch_tt, sp);
                    writeln!(
                        f,
                        "{i:>3}  {site:>5}  {iso_tt:>26}  {iso_utc:>26}  {ra_str:>20}  {dec_str:>20}"
                    )?;
                }
            }
        }
        Ok(())
    }
}
