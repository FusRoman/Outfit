//! # Tabular display for astrometric observations
//!
//! Pretty, zero-copy renderers to print an [`Observations`] collection
//! (a `SmallVec<[Observation; 6]>`) as a **table**.
//!
//! ## Overview
//!
//! The main entry point is the display adaptor [`ObservationsDisplay`]. It **borrows**
//! your observations and renders a formatted table when used with Rust formatting
//! (`{}` or `{:#}`), without cloning or moving data.
//!
//! Three layouts are available:
//!
//! - **Default** (compact, fixed-width):  
//!   `# | Site | MJD (TT) | RA[hms] ±σ["] | DEC[dms] ±σ["]`
//! - **Wide** (diagnostic, uses `comfy-table`):  
//!   adds `JD (TT) | RA [rad] | DEC [rad] | |r_geo| AU | |r_hel| AU`
//! - **ISO** (timestamp-centric, uses `comfy-table`):  
//!   replaces MJD/JD with `ISO (TT)` and `ISO (UTC)`
//!
//! ## Units & Conventions
//!
//! - **Time**: MJD/JD columns are on the **TT** scale. In ISO mode, both **TT** and **UTC**
//!   timestamps are shown (UTC includes leap seconds via `hifitime`).
//! - **Angles**: RA/DEC are formatted in sexagesimal (RA in **hours**, DEC in **degrees**).
//! - **Uncertainties**: printed in **arcseconds**, converted from radians with [`RAD2ARC`].
//! - **Positions**: vector norms (wide mode) are in **AU**, conventional **equatorial mean J2000**.
//!
//! ## Precision & Sorting
//!
//! - `with_seconds_precision(p)` — controls fractional digits for sexagesimal and ISO seconds.
//! - `with_distance_precision(p)` — controls fixed-point digits for AU distances (wide mode).
//! - `sorted()` — prints rows **sorted by epoch** (MJD TT ascending). The first column `#`
//!   always shows the **original index** (pre-sort) for traceability.
//!
//! ## Observer names
//!
//! If you pass an [`Outfit`] with [`ObservationsDisplay::with_env`], site labels render
//! as `"Name (#id)"` when available; otherwise the numeric **site id** is shown.
//!
//! ## Performance
//!
//! - The adaptor never clones/moves `Observation`s. It sorts a **vector of indices** when
//!   `sorted()` is used, and builds small transient strings per row (sexagesimal & ISO).
//! - **Default** layout writes fixed-width lines directly.  
//!   **Wide** and **ISO** layouts use [`comfy-table`] to build the table; this is still fast
//!   but implies allocating the table representation before printing.
//!
//! ## Quick examples
//!
//! ```rust,ignore
//! use outfit::observations::display::ObservationsDisplayExt;
//!
//! // 1) Compact table (fixed-width), sorted by epoch
//! println!("{}", observations.show().sorted());
//!
//! // 2) Wide table (adds JD, radians, |r| in AU), with custom precisions
//! println!("{}", observations
//!     .table_wide()
//!     .with_seconds_precision(4)
//!     .with_distance_precision(8)
//!     .sorted());
//!
//! // 3) ISO table (ISO TT + ISO UTC), resolve site names via Outfit
//! println!("{}", observations
//!     .table_iso()
//!     .with_env(&env)
//!     .with_seconds_precision(4)
//!     .sorted());
//!
//! // 4) Owned string (compact, unsorted)
//! let s = observations.show_string();
//! ```
//!
//! ## See also
//!
//! - [`Observation`] — single-observation pretty-printer and helpers.
//! - [`crate::conversion::ra_hms_prec`] / [`crate::conversion::dec_sdms_prec`]
//!   — sexagesimal decomposition with carry.
//! - [`crate::time::fmt_ss`] — seconds string `"SS.sss"` with 2-digit integer part.
//! - [`crate::time::iso_tt_from_epoch`] / [`crate::time::iso_utc_from_epoch`]
//!   — ISO renderers via `hifitime`.
//!
//! [`comfy-table`]: https://crates.io/crates/comfy-table
use std::fmt;

use hifitime::{Epoch, TimeScale};

use crate::constants::{JDTOMJD, RAD2ARC};
use crate::conversion::{dec_sdms_prec, ra_hms_prec};
use crate::observations::Observation;
use crate::time::{fmt_ss, iso_tt_from_epoch, iso_utc_from_epoch};
use crate::{Observations, Outfit};

use comfy_table::{presets::UTF8_FULL, Cell, CellAlignment, ContentArrangement, Row, Table};

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
    /// Optional Outfit environment to resolve observer names.
    env: Option<&'a Outfit>,
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
    site_label: String,
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
            env: None,
            mode: TableMode::Default,
            sec_prec: 3,
            dist_prec: 6,
            sorted: false,
        }
    }

    /// Switch to **wide** mode (adds JD, radians, vector norms).
    ///
    /// Adds the following columns:
    ///     - `JD (TT)`
    ///     - `RA [rad]`, `DEC [rad]`
    ///     - `|r_geo| AU`, `|r_hel| AU`
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
    /// Replaces the time columns with:
    ///     - `ISO (TT)` — Gregorian breakdown on TT,
    ///     - `ISO (UTC)` — TT converted to UTC (leap seconds handled by `hifitime`).
    ///
    /// Return
    /// ----------
    /// * `Self` (builder style), allowing chained configuration.
    ///
    /// Notes
    /// ----------
    /// * ISO strings are produced via `hifitime`, using leap-second tables for UTC.
    /// * This mode uses [`comfy-table`](https://docs.rs/comfy-table/latest/comfy_table/).
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
    pub fn sorted(mut self) -> Self {
        self.sorted = true;
        self
    }

    /// Attach an [`Outfit`] to resolve **observer names** in the `Site` column.
    ///
    /// If a name is available, rows show `"Name (#id)"`. Otherwise the numeric
    /// site id is displayed.
    ///
    /// Example
    /// -------
    /// ```rust,no_run
    /// println!("{}", observations.table_iso().with_env(&env).sorted());
    /// ```
    ///
    /// Arguments
    /// -----------------
    /// * `env` – The [`Outfit`] environment to use for resolving observer names.
    pub fn with_env(mut self, env: &'a Outfit) -> Self {
        self.env = Some(env);
        self
    }

    /// Generate the label for the observer site of a given observation.
    ///
    /// Arguments
    /// -----------------
    /// * `i` – Original index of the observation (pre-sort).
    ///
    /// Return
    /// -----------------
    /// * A string label for the site, either `"Name (#id)"` or `"#id"`.
    fn site_label(&self, i: usize) -> String {
        if let Some(env) = self.env {
            let site_id = self.obs[i].observer;
            let site = env.get_observer_from_uint16(site_id);
            if let Some(name) = site.name.as_deref() {
                if !name.is_empty() {
                    return format!("{name} (#{site_id})");
                }
            }
            format!("Site ID #{site_id}")
        } else {
            // No environment: keep a compact ID for compatibility
            format!("{}", self.obs[i].observer)
        }
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
                    "{:>3}  {:>5}  {:>14}  {:>14}  {:>26}  {:>11}  {:>26}  {:>11}  {:>12}  {:>12}",
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
                    "{:>3}  {:>5}  {:>26}  {:>26}  {:>26}  {:>26}",
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
            site_label: self.site_label(i),
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

    /// Render the WIDE table using comfy-table.
    fn render_wide_comfy(&self) -> String {
        let mut table = Table::new();
        table
            .load_preset(UTF8_FULL)
            .set_content_arrangement(ContentArrangement::Dynamic);

        // Header
        table.set_header(vec![
            Cell::new("#"),
            Cell::new("Site"),
            Cell::new("MJD (TT)"),
            Cell::new("JD (TT)"),
            Cell::new("RA ±σ[arcsec]"),
            Cell::new("RA [rad]"),
            Cell::new("DEC ±σ[arcsec]"),
            Cell::new("DEC [rad]"),
            Cell::new("|r_geo| AU"),
            Cell::new("|r_hel| AU"),
        ]);

        // Rows
        for (i, o) in self.row_iter() {
            let r = self.format_row_fields(i, o);
            let dp = self.dist_prec;

            table.add_row(Row::from(vec![
                Cell::new(r.i).set_alignment(CellAlignment::Right),
                Cell::new(r.site_label).set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.6}", r.mjd_tt)).set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.6}", r.jd_tt.unwrap_or_default()))
                    .set_alignment(CellAlignment::Right),
                Cell::new(r.ra_str.clone()).set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.7}", r.ra_rad)).set_alignment(CellAlignment::Right),
                Cell::new(r.dec_str.clone()).set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.7}", r.dec_rad)).set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.*}", dp, r.r_geo.unwrap_or_default()))
                    .set_alignment(CellAlignment::Right),
                Cell::new(format!("{:.*}", dp, r.r_hel.unwrap_or_default()))
                    .set_alignment(CellAlignment::Right),
            ]));
        }

        table.to_string()
    }

    /// Render the ISO table using comfy-table.
    fn render_iso_comfy(&self) -> String {
        let mut table = Table::new();
        table
            .load_preset(UTF8_FULL)
            .set_content_arrangement(ContentArrangement::Dynamic);

        // Header
        table.set_header(vec![
            Cell::new("#"),
            Cell::new("Site"),
            Cell::new("ISO (TT)"),
            Cell::new("ISO (UTC)"),
            Cell::new("RA ±σ[arcsec]"),
            Cell::new("DEC ±σ[arcsec]"),
        ]);

        // Rows
        for (i, o) in self.row_iter() {
            let r = self.format_row_fields(i, o);
            table.add_row(Row::from(vec![
                Cell::new(r.i).set_alignment(CellAlignment::Right),
                Cell::new(r.site_label).set_alignment(CellAlignment::Right),
                Cell::new(r.iso_tt.as_deref().unwrap_or("")).set_alignment(CellAlignment::Right),
                Cell::new(r.iso_utc.as_deref().unwrap_or("")).set_alignment(CellAlignment::Right),
                Cell::new(r.ra_str.clone()).set_alignment(CellAlignment::Right),
                Cell::new(r.dec_str.clone()).set_alignment(CellAlignment::Right),
            ]));
        }

        table.to_string()
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
                    site = r.site_label,
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
                    site = r.site_label,
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
                    site = r.site_label,
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

        match self.mode {
            TableMode::Wide => {
                // comfy-table rendering
                let out = self.render_wide_comfy();
                f.write_str(&out)?;
            }
            TableMode::Iso => {
                // comfy-table rendering
                let out = self.render_iso_comfy();
                f.write_str(&out)?;
            }
            TableMode::Default => {
                // Legacy compact mode: keep your existing fixed-width rendering
                self.write_header(f)?;
                for (i, o) in self.row_iter() {
                    let row = self.format_row_fields(i, o);
                    self.write_row(f, &row)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod observation_display_tests {
    use super::*;
    use crate::observations::display::ObservationsDisplayExt;
    use crate::Observations;
    use nalgebra::Vector3;

    /// Build a minimal Observation for tests.
    /// Angles are provided in degrees; uncertainties in arcseconds; time is MJD (TT).
    fn make_obs(
        site: u16,
        ra_deg: f64,
        dec_deg: f64,
        err_arcsec: f64,
        mjd_tt: f64,
        rgeo: (f64, f64, f64),
        rhel: (f64, f64, f64),
    ) -> Observation {
        // Convert helpers
        let ra_rad = ra_deg.to_radians();
        let dec_rad = dec_deg.to_radians();
        let err_rad = (err_arcsec / 3600.0).to_radians();

        Observation {
            observer: site,
            ra: ra_rad,         // Radian from f64
            error_ra: err_rad,  // Radian from f64
            dec: dec_rad,       // Radian from f64
            error_dec: err_rad, // Radian from f64
            time: mjd_tt,       // MJD from f64
            observer_earth_position: Vector3::new(rgeo.0, rgeo.1, rgeo.2),
            observer_helio_position: Vector3::new(rhel.0, rhel.1, rhel.2),
        }
    }

    /// Build a small, heterogeneous set of observations for table tests.
    fn sample_observations() -> Observations {
        let mut obs: Observations = Observations::default();
        // idx = 0 (later than #1), zero vectors to simplify wide-mode distance checks
        obs.push(make_obs(
            809,
            0.0,
            0.0,
            1.0,
            60000.123456,
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ));
        // idx = 1 (earliest), negative DEC to verify sign & DMS formatting
        obs.push(make_obs(
            2,
            0.0,
            -10.0,
            0.5,
            59000.0,
            (0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0),
        ));
        // idx = 2 (latest), non-zero vectors to exercise wide-mode distances (1 and 2 AU)
        obs.push(make_obs(
            500,
            180.0,
            20.0,
            1.2,
            60001.0,
            (1.0, 0.0, 0.0),
            (0.0, 2.0, 0.0),
        ));
        obs
    }

    #[test]
    fn default_headers_and_basic_format() {
        let obs = sample_observations();
        let s = format!("{}", obs.show()); // Default, unsorted

        // Headers present
        assert!(s.contains("MJD (TT)"));
        assert!(s.contains("RA ±σ[arcsec]"));
        assert!(s.contains("DEC ±σ[arcsec]"));

        // RA/DEC sexagesimal + uncertainties for idx 0 (RA=0h, DEC=+0°)
        // Cells do not include 'RA=' / 'DEC=' prefixes in the table.
        assert!(s.contains("00h00m00.000s ± 1.000\""));
        assert!(s.contains("+00°00'00.000\" ± 1.000\""));
    }

    #[test]
    fn dec_negative_sign_is_preserved_in_table() {
        let obs = sample_observations();
        let s = format!("{}", obs.show()); // Default, unsorted
                                           // idx 1 has DEC = -10° (no 'DEC=' prefix in table cells)
        assert!(
            s.contains("-10°00'00.000\""),
            "Negative DEC sign or DMS formatting incorrect: {s}"
        );
        // And its sigma is 0.500"
        assert!(s.contains("± 0.500\""));
    }

    #[test]
    fn sorted_orders_by_time_and_keeps_original_index() {
        let obs = sample_observations();
        let s = format!("{}", obs.show().sorted());

        // Split lines: title, underline, header, then data rows...
        let mut lines = s.lines();
        let _title = lines.next().unwrap_or_default();
        let _rule = lines.next().unwrap_or_default();
        let _hdr = lines.next().unwrap_or_default();

        // First data row should correspond to the earliest epoch (idx = 1)
        let first_row = lines.next().unwrap_or_default();
        assert!(
            first_row.trim_start().starts_with("1"),
            "Expected first printed row to be original index 1, got: {first_row}"
        );
        // A later row should include index 0 (natural order changed by sort)
        let rest = lines.collect::<Vec<_>>().join("\n");
        assert!(
            rest.contains("\n  0") || rest.starts_with("  0"),
            "Index 0 should also appear."
        );
    }

    #[test]
    fn wide_mode_headers_and_radians_and_distances() {
        let obs = sample_observations();
        let s = format!("{}", obs.table_wide()); // Wide, unsorted

        // Headers present
        assert!(s.contains("JD (TT)"));
        assert!(s.contains("RA [rad]"));
        assert!(s.contains("DEC [rad]"));
        assert!(s.contains("|r_geo| AU"));
        assert!(s.contains("|r_hel| AU"));

        // idx 1 has DEC = -10° => about -0.1745329 rad (7 decimals printed)
        assert!(
            s.contains("-0.1745329"),
            "Expected DEC in radians around -0.1745329 rad in wide mode: {s}"
        );

        // idx 0 vectors are zeros => distances should include 0.000000
        assert!(
            s.contains("  0.000000"),
            "Expected at least one zero distance in wide mode for idx 0: {s}"
        );

        // idx 2 vectors norms are 1 AU and 2 AU
        assert!(
            s.contains("  1.000000") && s.contains("  2.000000"),
            "Expected distances 1.000000 and 2.000000 AU for idx 2: {s}"
        );
    }

    #[test]
    fn iso_mode_headers_and_suffixes() {
        let obs = sample_observations();
        let s = format!("{}", obs.table_iso()); // ISO, unsorted

        // Headers present
        assert!(s.contains("ISO (TT)"));
        assert!(s.contains("ISO (UTC)"));

        // Body should contain ' TT' and 'Z' suffixes (at least once)
        assert!(s.contains(" TT"), "TT suffix missing in ISO TT column: {s}");
        assert!(s.contains('Z'), "Z suffix missing in ISO UTC column: {s}");
    }

    #[test]
    fn seconds_and_distance_precision_knobs() {
        let obs = sample_observations();

        // Seconds precision = 4: "00.0000" for RA seconds at 0h (no 'RA=' prefix in table cells)
        let s_iso = format!("{}", obs.table_iso().with_seconds_precision(4));
        assert!(
            s_iso.contains("00h00m00.0000s"),
            "Seconds precision not applied: {s_iso}"
        );

        // Distance precision = 4: look for "  0.0000" in wide mode (idx 0 has zero vectors)
        let s_wide = format!("{}", obs.table_wide().with_distance_precision(4));
        assert!(
            s_wide.contains("  0.0000"),
            "Distance precision not applied (expected 4 decimals): {s_wide}"
        );
    }

    #[test]
    fn show_string_matches_display_default() {
        let obs = sample_observations();
        let s1 = obs.show_string();
        let s2 = format!("{}", obs.show());
        assert_eq!(s1, s2, "show_string() must match Display in default mode");
    }
}
