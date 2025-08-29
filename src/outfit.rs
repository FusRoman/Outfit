//! # Outfit: environment, ephemerides, and observatory registry
//!
//! This module defines the [`Outfit`](crate::outfit::Outfit) struct, the central façade that wires together:
//!
//! 1. **Environment state** ([`OutfitEnv`](crate::env_state::OutfitEnv)) — providers and configuration (e.g., UT1).
//! 2. **JPL ephemerides access** — lazy, cached handle over a chosen source
//!    ([`EphemFileSource`](crate::jpl_ephem::download_jpl_file::EphemFileSource) → [`JPLEphem`](crate::jpl_ephem::JPLEphem)).
//! 3. **Observatory registry** — MPC Observatory Codes parsed into [`Observer`](crate::observers::Observer) instances,
//!    with stable integer IDs for compact indexing and storage.
//! 4. **Astrometric error models** — per-site bias/RMS lookup for RA/DEC accuracies.
//!
//! The design emphasizes *lazy initialization* and *idempotent caching*:
//! - The ephemeris file is opened on first use via [`OnceCell`](once_cell::sync::OnceCell), then reused.
//! - The MPC observatory table is fetched and parsed once, then retained.
//!
//! ## Key responsibilities
//!
//! - Single source of truth for **JPL ephemerides** (HORIZONS/NAIF) through [`get_jpl_ephem`](crate::outfit::Outfit::get_jpl_ephem)
//! - Access to **UT1 provider** for Earth-rotation dependent calculations
//! - **MPC observatory code → Observer** resolution and the inverse (**Observer → u16 index**)
//! - Enrichment of observers with **bias/RMS** angular accuracies from the configured
//!   [`ErrorModel`](crate::error_models::ErrorModel) (e.g., *FCCT14*)
//!
//! ## Typical usage
//!
//! ```rust, no_run
//! use outfit::outfit::Outfit;
//! use outfit::error_models::ErrorModel;
//!
//! // Instantiate the context with a JPL source and an error model
//! let outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
//!
//! // On-demand: the ephemeris is opened only once and cached
//! let jpl = outfit.get_jpl_ephem().unwrap();
//!
//! // Resolve an observer by MPC code
//! let haleakala = outfit.get_observer_from_mpc_code(&"F51".into());
//! ```
//!
//! ## Notes
//!
//! - The MPC table is pulled from:
//!   `https://www.minorplanetcenter.net/iau/lists/ObsCodes.html`
//!   and minimally parsed from its `<pre>` text block.
//! - Per-site **RA/DEC accuracies** (bias+RMS) are looked up with [`get_bias_rms`](crate::error_models::get_bias_rms),
//!   currently assuming catalog code `"c"` unless indicated otherwise (TODO).
//!
//! ## See also
//! ------------
//! * [`JPLEphem`](crate::jpl_ephem::JPLEphem) – Ephemerides access layer.
//! * [`Observer`](crate::observers::Observer) – Geodetic/parallax observer representation with optional RA/DEC accuracies.
//! * [`ErrorModel`](crate::error_models::ErrorModel) / [`get_bias_rms`](crate::error_models::get_bias_rms) – Site accuracy enrichment.
//! * [`OutfitEnv`](crate::env_state::OutfitEnv) – Providers (e.g., UT1) and environment state.
//! * [`EphemFileSource`](crate::jpl_ephem::download_jpl_file::EphemFileSource) – Source selection for JPL files (HORIZONS/NAIF).
//!
//! ## Panics & errors
//!
//! - Functions that *must* find an MPC code will `panic!` if the code is unknown.
//!   Prefer adding a fallible variant if you need graceful handling.
//! - I/O and parsing failures are surfaced as [`OutfitError`](crate::outfit_errors::OutfitError) where applicable.

use std::{collections::HashMap, sync::Arc};

use nalgebra::Matrix3;
use once_cell::sync::OnceCell;

use crate::{
    constants::{Degree, Kilometer, MpcCode, MpcCodeObs},
    env_state::OutfitEnv,
    error_models::{get_bias_rms, ErrorModel, ErrorModelData},
    jpl_ephem::download_jpl_file::EphemFileSource,
    observers::{observatories::Observatories, Observer},
    outfit_errors::OutfitError,
    ref_system::{rotpn, RefEpoch, RefSystem},
};

use crate::jpl_ephem::JPLEphem;

#[derive(Debug, Clone)]
pub struct Outfit {
    env_state: OutfitEnv,
    observatories: Observatories,
    jpl_source: EphemFileSource,
    jpl_ephem: OnceCell<JPLEphem>,
    error_model: ErrorModelData,
    rot_equmj2000_to_eclmj2000: Matrix3<f64>,
    rot_eclmj2000_to_equmj2000: Matrix3<f64>,
}

impl Outfit {
    /// Construct a new [`Outfit`] context.
    ///
    /// Initializes the environment, sets the JPL ephemerides source, and loads the configured
    /// error model from disk. The ephemeris file itself is **not** opened yet; it is lazily
    /// initialized the first time [`get_jpl_ephem`](crate::outfit::Outfit::get_jpl_ephem) is called.
    ///
    /// Arguments
    /// -----------------
    /// * `jpl_file`: A source descriptor resolvable into an [`EphemFileSource`]
    ///   (e.g., `"horizon:DE440"` or a NAIF path).
    /// * `error_model`: The site accuracy model to load (e.g., [`ErrorModel::FCCT14`]).
    ///
    /// Return
    /// ----------
    /// * A new [`Outfit`] instance or an [`OutfitError`] if the error model cannot be read.
    ///
    /// See also
    /// ------------
    /// * [`get_jpl_ephem`](crate::outfit::Outfit::get_jpl_ephem) – Lazy initialization and access to the ephemeris handle.
    /// * [`ErrorModel::read_error_model_file`](crate::error_models::ErrorModel::read_error_model_file) – Underlying loader for the model data.
    pub fn new(jpl_file: &str, error_model: ErrorModel) -> Result<Self, OutfitError> {
        let rot1 = rotpn(
            &RefSystem::Equm(RefEpoch::J2000),
            &RefSystem::Eclm(RefEpoch::J2000),
        )?;

        let rot2 = rotpn(
            &RefSystem::Eclm(RefEpoch::J2000),
            &RefSystem::Equm(RefEpoch::J2000),
        )?;

        Ok(Outfit {
            env_state: OutfitEnv::new(),
            observatories: Observatories::new(),
            jpl_source: jpl_file.try_into()?,
            jpl_ephem: OnceCell::new(),
            error_model: error_model.read_error_model_file()?,
            rot_equmj2000_to_eclmj2000: rot1,
            rot_eclmj2000_to_equmj2000: rot2,
        })
    }

    /// Get the rotation matrix from equatorial J2000 to ecliptic J2000.
    /// This matrix is used to transform coordinates from the equatorial frame to the ecliptic frame.
    pub fn get_rot_equmj2000_to_eclmj2000(&self) -> &Matrix3<f64> {
        &self.rot_equmj2000_to_eclmj2000
    }

    pub fn get_rot_eclmj2000_to_equmj2000(&self) -> &Matrix3<f64> {
        &self.rot_eclmj2000_to_equmj2000
    }

    /// Get the lazily-initialized JPL ephemerides handle.
    ///
    /// If this is the first call, the ephemeris is opened and cached in an internal [`OnceCell`].
    /// Subsequent calls return the same reference.
    ///
    /// Arguments
    /// -----------------
    /// *None*
    ///
    /// Return
    /// ----------
    /// * `&JPLEphem` on success, or an [`OutfitError`] if the source cannot be opened.
    ///
    /// See also
    /// ------------
    /// * [`EphemFileSource`] – Source configuration.
    /// * [`OnceCell::get_or_try_init`] – Lazy initialization helper.
    pub fn get_jpl_ephem(&self) -> Result<&JPLEphem, OutfitError> {
        self.jpl_ephem
            .get_or_try_init(|| JPLEphem::new(&self.jpl_source))
    }

    /// Access the UT1 provider from the environment.
    ///
    /// This is useful for Earth-rotation dependent calculations (e.g., GMST, sidereal time).
    ///
    /// Arguments
    /// -----------------
    /// *None*
    ///
    /// Return
    /// ----------
    /// * A reference to the [`hifitime::ut1::Ut1Provider`].
    ///
    /// See also
    /// ------------
    /// * [`OutfitEnv`] – Environment state and providers.
    pub fn get_ut1_provider(&self) -> &hifitime::ut1::Ut1Provider {
        &self.env_state.ut1_provider
    }

    /// Get the lazily built MPC observatory map (MPC code → [`Observer`]).
    ///
    /// The map is fetched and parsed from the MPC HTML table on first use, then cached.
    ///
    /// Return
    /// ----------
    /// * A reference to the shared map: [`MpcCodeObs`] = `HashMap<MpcCode, Arc<Observer>>`.
    ///
    /// See also
    /// ------------
    /// * [`init_observatories`](crate::outfit::Outfit::init_observatories) – Builder invoked on first access.
    /// * [`get_observer_from_mpc_code`](crate::outfit::Outfit::get_observer_from_mpc_code) – Convenience accessor for one site.
    pub(crate) fn get_observatories(&self) -> &MpcCodeObs {
        self.observatories
            .mpc_code_obs
            .get_or_init(|| self.init_observatories())
    }

    /// Resolve an [`Observer`] from a given MPC observatory code.
    ///
    /// This accessor panics if the code is unknown. Use it when unknown codes are exceptional.
    ///
    /// Arguments
    /// -----------------
    /// * `mpc_code`: The MPC observatory code (e.g., `"F51"`).
    ///
    /// Return
    /// ----------
    /// * An `Arc<Observer>` for the requested site.
    pub fn get_observer_from_mpc_code(&self, mpc_code: &MpcCode) -> Arc<Observer> {
        self.get_observatories()
            .get(mpc_code)
            .unwrap_or_else(|| panic!("MPC code not found: {mpc_code}"))
            .clone()
    }

    /// Build the MPC observatory registry by fetching and parsing the MPC list.
    ///
    /// For each row, the routine extracts:
    /// - Longitude (deg), ρ·cosφ, ρ·sinφ (parallax factors),
    /// - Human-readable name,
    /// - Optional RA/DEC accuracies derived from the loaded [`ErrorModelData`]
    ///   via [`get_bias_rms`] (currently using catalog code `"c"`, TODO).
    ///
    /// Return
    /// ----------
    /// * A freshly constructed [`MpcCodeObs`] map.
    ///
    /// See also
    /// ------------
    /// * [`get_observatories`](crate::outfit::Outfit::get_observatories) – Lazy wrapper that caches this map.
    /// * [`get_bias_rms`] – Site accuracy lookup by (mpc_code, catalog_code).
    pub(crate) fn init_observatories(&self) -> MpcCodeObs {
        let mut observatories: MpcCodeObs = HashMap::new();

        let mpc_code_response = self
            .env_state
            .get_from_url("https://www.minorplanetcenter.net/iau/lists/ObsCodes.html");

        let mpc_code_csv = mpc_code_response
            .trim()
            .strip_prefix("<pre>")
            .and_then(|s| s.strip_suffix("</pre>"))
            .expect("Failed to strip pre tags");

        for lines in mpc_code_csv.lines().skip(2) {
            let line = lines.trim();

            if let Some((code, remain)) = line.split_at_checked(3) {
                let remain = remain.trim_end();

                let (longitude, cos, sin, name) = parse_remain(remain, code);

                // TODO: support per-site catalog codes (not always "c")
                let bias_rms = get_bias_rms(&self.error_model, code.to_string(), "c".to_string());

                let observer = Observer::from_parallax(
                    longitude as f64,
                    cos as f64,
                    sin as f64,
                    Some(name),
                    bias_rms.map(|(ra, _)| ra as f64),
                    bias_rms.map(|(_, dec)| dec as f64),
                )
                .expect("Failed to create observer");
                observatories.insert(code.to_string(), Arc::new(observer));
            };
        }
        observatories
    }

    /// Convert an MPC code to its stable 16-bit observatory index.
    ///
    /// Useful for compact storage of observer references in catalogs, measurements,
    /// and ephemeris products.
    ///
    /// Arguments
    /// -----------------
    /// * `mpc_code`: The MPC observatory code.
    ///
    /// Return
    /// ----------
    /// * The `u16` index associated with the given observer.
    ///
    /// See also
    /// ------------
    /// * [`get_observer_from_mpc_code`](crate::outfit::Outfit::get_observer_from_mpc_code) – Resolve the observer first (panic on unknown).
    /// * [`uint16_from_observer`](crate::outfit::Outfit::uint16_from_observer) – Indexing for arbitrary/new observers.
    pub(crate) fn uint16_from_mpc_code(&mut self, mpc_code: &MpcCode) -> u16 {
        let observer = self.get_observer_from_mpc_code(mpc_code);
        self.observatories.uint16_from_observer(observer)
    }

    /// Convert an [`Observer`] handle to its stable 16-bit index.
    ///
    /// Arguments
    /// -----------------
    /// * `observer`: The observer to be indexed.
    ///
    /// Return
    /// ----------
    /// * The `u16` index associated with this observer (inserting if new).
    ///
    /// See also
    /// ------------
    /// * [`get_observer_from_uint16`](crate::outfit::Outfit::get_observer_from_uint16) – Recover a reference from an index.
    /// * [`new_observer`](crate::outfit::Outfit::new_observer) – Create and register a new custom observer.
    pub(crate) fn uint16_from_observer(&mut self, observer: Arc<Observer>) -> u16 {
        self.observatories.uint16_from_observer(observer)
    }

    /// Recover an [`Observer`] reference from a 16-bit index.
    ///
    /// Arguments
    /// -----------------
    /// * `observer_idx`: The previously assigned index.
    ///
    /// Return
    /// ----------
    /// * A reference to the corresponding [`Observer`].
    ///
    /// See also
    /// ------------
    /// * [`uint16_from_observer`](crate::outfit::Outfit::uint16_from_observer) – Assign/lookup indices.
    /// * [`get_observer_from_mpc_code`](crate::outfit::Outfit::get_observer_from_mpc_code) – Resolve by MPC code instead.
    pub(crate) fn get_observer_from_uint16(&self, observer_idx: u16) -> &Observer {
        self.observatories.get_observer_from_uint16(observer_idx)
    }

    /// Create and register a new **custom** observer.
    ///
    /// This helper converts geodetic inputs to the internal parallax representation
    /// (ρ·cosφ, ρ·sinφ) and stores the new [`Observer`] with an optional display name.
    ///
    /// Arguments
    /// -----------------
    /// * `longitude`: Geodetic longitude in **degrees** (east-positive).
    /// * `latitude`:  Geodetic latitude in **degrees**.
    /// * `elevation`: Elevation in **kilometers** above the ellipsoid/geoid (model-dependent).
    /// * `name`:      Optional human-readable name for the site.
    ///
    /// Return
    /// ----------
    /// * An `Arc<Observer>` handle to the newly created observer.
    pub fn new_observer(
        &mut self,
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
    ) -> Arc<Observer> {
        self.observatories
            .create_observer(longitude, latitude, elevation, name)
    }

    pub(crate) fn add_observer_internal(&mut self, observer: Arc<Observer>) -> u16 {
        self.observatories.add_observer(observer)
    }

    pub fn add_observer(&mut self, observer: Arc<Observer>) {
        self.add_observer_internal(observer);
    }
}

/// Parse a fixed-width float slice from an MPC observatory line.
///
/// Arguments
/// -----------------
/// * `s`:    Full line (trailing part after the 3-char MPC code).
/// * `slice`: Byte range selecting the numeric field.
/// * `code`:  MPC code (for diagnostics).
///
/// Return
/// ----------
/// * `Ok(f32)` parsed value, or the parsing error.
///
/// See also
/// ------------
/// * [`parse_remain`] – Higher-level field extraction for one line.
fn parse_f32(
    s: &str,
    slice: std::ops::Range<usize>,
    code: &str,
) -> Result<f32, std::num::ParseFloatError> {
    s.get(slice)
        .unwrap_or_else(|| panic!("Failed to parse float for observer code: {code}"))
        .trim()
        .parse()
}

/// Extract longitude, ρ·cosφ, ρ·sinφ, and name from a fixed-width MPC row.
///
/// This helper returns partial values (with zeros) if any field fails to parse,
/// allowing the caller to still record the site name while signaling missing data
/// implicitly through zeros.
///
/// Arguments
/// -----------------
/// * `remain`: Fixed-width tail of the line (after the 3-char MPC code).
/// * `code`:   MPC code (for diagnostics).
///
/// Return
/// ----------
/// * `(longitude_deg, rho_cos_phi, rho_sin_phi, name)`
fn parse_remain(remain: &str, code: &str) -> (f32, f32, f32, String) {
    let name = remain
        .get(27..)
        .unwrap_or_else(|| panic!("Failed to parse name value for code: {code}"));

    let Some(longitude) = parse_f32(remain, 1..10, code).ok() else {
        return (0.0, 0.0, 0.0, name.to_string());
    };

    let Some(cos) = parse_f32(remain, 10..18, code).ok() else {
        return (longitude, 0.0, 0.0, name.to_string());
    };

    let Some(sin) = parse_f32(remain, 18..27, code).ok() else {
        return (longitude, cos, 0.0, name.to_string());
    };
    (longitude, cos, sin, name.to_string())
}
