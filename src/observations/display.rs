//! # Tabular display for astrometric observations
//!
//! This module provides zero-copy, table-style renderers for an
//! [`Observations`] collection (a `SmallVec<[Observation; 6]>`).
//!
//! ## Overview
//!
//! The main entry point is the **display adaptor** [`ObservationsDisplay`], which borrows
//! the underlying observation slice and prints a formatted table via Rust’s
//! formatting machinery (`{}`).
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
//! **Sorting option**
//! -----------------
//! Use [`ObservationsDisplay::sorted`] to print rows **sorted by epoch** (MJD TT).
//! The `#` column always shows the **original index** (pre-sort) for traceability.
//!
//! ## Units & conventions
//!
//! - **Time**: MJD and JD columns are on the **TT** scale (ISO mode prints **TT** and **UTC**).
//! - **Angles**: RA/DEC are rendered in sexagesimal (RA in **hours**, DEC in **degrees**).
//! - **Uncertainties**: shown in **arcseconds**; converted from radians using [`RAD2ARC`].
//! - **Positions**: vector norms are in **AU** (equatorial mean J2000 frame by convention).
//!
//! ## Performance
//!
//! - The adaptor does **not** move or clone observations; it iterates indices and formats rows.
//! - Per-row small `String`s (sexagesimal and ISO fields) are allocated transiently.
//! - With `sorted(true)`, only a temporary **index vector** is sorted; the data order is preserved.
//!
//! ## Examples
//!
//! ```rust,ignore
//! use outfit::observations::display::ObservationsDisplayExt;
//!
//! // Compact table, sorted by epoch
//! println!("{}", observations.show().sorted(true));
//!
////! // Wide table (adds JD, radians, |r| in AU), sorted
//! println!("{}", observations.table_wide().sorted(true));
//!
//! // ISO table (ISO TT + ISO UTC), custom seconds precision, sorted
//! println!("{}", observations.table_iso().with_seconds_precision(4).sorted(true));
//!
//! // Get an owned string (compact, unsorted)
//! let s = observations.show_string();
//! ```
//!
//! ## See also
//!
//! - [`Observation`] — single-observation pretty-printer and helpers
//! - [`crate::conversion::ra_hms_prec`] / [`crate::conversion::dec_sdms_prec`] — sexagesimal decomposition with carry
//! - [`crate::time::fmt_ss`] — seconds string `"SS.sss"` with 2-digit integer part
//! - [`crate::time::iso_tt_from_epoch`] / [`crate::time::iso_utc_from_epoch`] — ISO renderers via `hifitime`
//!
//! ## Notes
//!
//! - Header widths are chosen for readability; they are not intended to be parsed.
//! - ISO rendering relies on `hifitime`’s leap-second tables for UTC.

use std::fmt;

use hifitime::{Epoch, TimeScale};

use crate::constants::{JDTOMJD, RAD2ARC};
use crate::conversion::{dec_sdms_prec, ra_hms_prec};
use crate::observations::Observation;
use crate::time::{fmt_ss, iso_tt_from_epoch, iso_utc_from_epoch};
use crate::Observations;

/// Internal layout selector for the table renderer.
///
/// This enum is crate-internal and selects the columns printed by
/// [`ObservationsDisplay`]. It is not part of the public API.
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
/// Sorting
/// -----------------
/// * Call [`Self::sorted`] to display rows **sorted by MJD (TT)** in ascending order.
/// * The `#` column always shows the **original index** (pre-sort) for traceability.
/// * Ties (identical epochs) keep a stable order by original index.
///
/// See also
/// ------------
/// * [`ObservationsDisplayExt`] – Ergonomic builders for each mode.
/// * [`Observation`] – Per-row semantics used to derive the columns.
pub struct ObservationsDisplay<'a> {
    /// Borrowed collection to render. No allocation or copying occurs.
    obs: &'a Observations,
    /// Column layout selector (default / wide / iso).
    mode: TableMode,
    /// Fractional digits for sexagesimal and ISO seconds (default = 3).
    sec_prec: usize,
    /// Fixed-point digits for AU distances in wide mode (default = 6).
    dist_prec: usize,
    /// If `true`, rows are printed **sorted by epoch** (MJD TT) ascending.
    sorted: bool,
}

/// Pre-computed, per-row fields shared across all modes.
///
/// Notes
/// ----------
/// * This structure is internal to avoid recomputing formatting across modes.
/// * Optional fields are only populated in the relevant mode (e.g., `jd_tt` in **Wide**).
struct RowFields {
    i: usize,
    site: u16,
    // Base time & angles
    mjd_tt: f64,
    ra_rad: f64,
    dec_rad: f64,
    // Rendered sexagesimal with uncertainties
    ra_str: String,
    dec_str: String,
    // Optional extras by mode
    jd_tt: Option<f64>,
    r_geo: Option<f64>,
    r_hel: Option<f64>,
    iso_tt: Option<String>,
    iso_utc: Option<String>,
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
            sorted: false,
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
    /// * `p` – Fixed-point fractional digits for the `|r_geo|` and `|r_hel|` columns.
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style).
    pub fn with_distance_precision(mut self, p: usize) -> Self {
        self.dist_prec = p;
        self
    }

    /// Enable/disable **time-sorted** display.
    ///
    /// Arguments
    /// -----------------
    /// * `yes` – If `true`, rows are sorted by `Observation::time` (MJD TT) ascending.
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style), allowing chained configuration.
    ///
    /// Notes
    /// ----------
    /// * Sorting uses a **stable** index order (no reordering or cloning of observations).
    /// * The `#` column prints the **original index** (pre-sort).
    pub fn sorted(mut self, yes: bool) -> Self {
        self.sorted = yes;
        self
    }

    /// Write the table header according to the selected mode.
    fn write_header(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.mode {
            TableMode::Default => {
                writeln!(
                    f,
                    "{:>3}  {:>5}  {:>14}  {:>20}  {:>20}",
                    "#", "Site", "MJD (TT)", "RA ±σ[arcsec]", "DEC ±σ[arcsec]"
                )
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
                )
            }
            TableMode::Iso => {
                writeln!(
                    f,
                    "{:>3}  {:>5}  {:>26}  {:>26}  {:>20}  {:>20}",
                    "#", "Site", "ISO (TT)", "ISO (UTC)", "RA ±σ[arcsec]", "DEC ±σ[arcsec]"
                )
            }
        }
    }

    /// Build an iterator over `(original_index, &Observation)` honoring the `sorted` flag.
    ///
    /// Return
    /// ----------
    /// * A boxed iterator over `(index_before_sort, &Observation)`.
    ///
    /// Notes
    /// ----------
    /// * When `sorted == true`, indices are ordered by MJD(TT) ascending, with a stable
    ///   tie-break on the original index.
    fn row_iter(&self) -> Box<dyn Iterator<Item = (usize, &Observation)> + '_> {
        if self.sorted {
            use std::cmp::Ordering;
            let mut order: Vec<usize> = (0..self.obs.len()).collect();
            order.sort_by(|&a, &b| {
                let ta: f64 = self.obs[a].time;
                let tb: f64 = self.obs[b].time;
                match ta.partial_cmp(&tb) {
                    Some(ord) => ord,
                    None => Ordering::Equal, // NaN-safe
                }
                .then_with(|| a.cmp(&b))
            });
            Box::new(order.into_iter().map(|i| (i, &self.obs[i])))
        } else {
            Box::new(self.obs.iter().enumerate())
        }
    }

    /// Compute common formatted values once for a given row.
    ///
    /// Arguments
    /// -----------------
    /// * `i` – Original index of the observation (pre-sort).
    /// * `o` – Borrowed [`Observation`] for this row.
    ///
    /// Return
    /// ----------
    /// * A populated [`RowFields`] struct reused by the row writer.
    fn format_row_fields(&self, i: usize, o: &Observation) -> RowFields {
        let sp = self.sec_prec;

        // Sexagesimal decomposition
        let ra_rad: f64 = o.ra;
        let dec_rad: f64 = o.dec;
        let (hh, mm, ss) = ra_hms_prec(ra_rad, sp);
        let (sgn, dd, dm, ds) = dec_sdms_prec(dec_rad, sp);
        let ss_s = fmt_ss(ss, sp);
        let ds_s = fmt_ss(ds, sp);

        // Uncertainties [arcsec]
        let sra_as = o.error_ra * RAD2ARC;
        let sdec_as = o.error_dec * RAD2ARC;

        // Common formatted strings
        let ra_str = format!("{hh:02}h{mm:02}m{ss_s}s ± {sra_as:.3}\"");
        let dec_str = format!("{sgn}{dd:02}°{dm:02}'{ds_s}\" ± {sdec_as:.3}\"");

        // Base time
        let mjd_tt: f64 = o.time;

        // Mode-dependent extras (computed lazily below)
        let (jd_tt, r_geo, r_hel, iso_tt, iso_utc) = match self.mode {
            TableMode::Default => (None, None, None, None, None),
            TableMode::Wide => {
                let jd = mjd_tt + JDTOMJD;
                let g = o.observer_earth_position.norm();
                let h = o.observer_helio_position.norm();
                (Some(jd), Some(g), Some(h), None, None)
            }
            TableMode::Iso => {
                let epoch_tt = Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT);
                let tt_str = iso_tt_from_epoch(epoch_tt, sp);
                let utc_str = iso_utc_from_epoch(epoch_tt, sp);
                (None, None, None, Some(tt_str), Some(utc_str))
            }
        };

        RowFields {
            i,
            site: o.observer,
            mjd_tt,
            ra_rad,
            dec_rad,
            ra_str,
            dec_str,
            jd_tt,
            r_geo,
            r_hel,
            iso_tt,
            iso_utc,
        }
    }

    /// Write a single table row using pre-computed [`RowFields`].
    ///
    /// Arguments
    /// -----------------
    /// * `f` – Destination formatter.
    /// * `r` – Pre-formatted row fields.
    fn write_row(&self, f: &mut fmt::Formatter<'_>, r: &RowFields) -> fmt::Result {
        match self.mode {
            TableMode::Default => {
                writeln!(
                    f,
                    "{i:>3}  {site:>5}  {mjd:>14.6}  {ra:>20}  {dec:>20}",
                    i = r.i,
                    site = r.site,
                    mjd = r.mjd_tt,
                    ra = r.ra_str,
                    dec = r.dec_str
                )
            }
            TableMode::Wide => {
                let dp = self.dist_prec;
                writeln!(
                    f,
                    "{i:>3}  {site:>5}  {mjd:>14.6}  {jd:>14.6}  {ra:>20}  {ra_rad:>11.7}  {dec:>20}  {dec_rad:>11.7}  {rgeo:>12.dp$}  {rhel:>12.dp$}",
                    i = r.i,
                    site = r.site,
                    mjd = r.mjd_tt,
                    jd = r.jd_tt.unwrap_or_default(),
                    ra = r.ra_str,
                    ra_rad = r.ra_rad,
                    dec = r.dec_str,
                    dec_rad = r.dec_rad,
                    rgeo = r.r_geo.unwrap_or_default(),
                    rhel = r.r_hel.unwrap_or_default(),
                    dp = dp
                )
            }
            TableMode::Iso => {
                writeln!(
                    f,
                    "{i:>3}  {site:>5}  {iso_tt:>26}  {iso_utc:>26}  {ra:>20}  {dec:>20}",
                    i = r.i,
                    site = r.site,
                    iso_tt = r.iso_tt.as_deref().unwrap_or(""),
                    iso_utc = r.iso_utc.as_deref().unwrap_or(""),
                    ra = r.ra_str,
                    dec = r.dec_str
                )
            }
        }
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
    /// * Rows are printed in the **current order** unless `sorted(true)` is set.
    /// * When sorted, ordering is by MJD(TT) ascending, with ties broken by original index.
    /// * Numeric fields are right-aligned; headers have fixed widths for readability.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let n = self.obs.len();
        writeln!(f, "Observations (n={n})")?;
        writeln!(f, "-------------------")?;
        self.write_header(f)?;

        for (i, o) in self.row_iter() {
            let row = self.format_row_fields(i, o);
            self.write_row(f, &row)?;
        }
        Ok(())
    }
}
