//! Unified access to JPL ephemerides (legacy DE binaries and NAIF SPK/DAF).
//!
//! This module abstracts two **binary** JPL ephemeris formats behind a single API:
//!
//! - **Legacy JPL DE binaries** (TTL/CNAM/IPT layout) handled by the [`horizon`](crate::jpl_ephem::horizon) backend:
//!   parsed into Chebyshev segments per body and interpolated on demand.
//! - **NAIF SPK/DAF kernels** handled by the [`naif`](crate::jpl_ephem::naif) backend: parsed via a DAF
//!   header + summary/directory records into Chebyshev segments (SPK types).
//!
//! The enum [`JPLEphem`](crate::jpl_ephem::JPLEphem) wraps either backend and exposes a common entry point to
//! get Earth state vectors via [`JPLEphem::earth_ephemeris`](crate::jpl_ephem::JPLEphem::earth_ephemeris).
//!
//! # Frames, centers, and units
//!
//! - **Legacy DE backend (`horizon`)**
//!   - **Query:** `Earth` relative to `Sun`.
//!   - **Internal units:** kilometers; velocity in kilometers **per day** (JD spacing).
//!   - After calling `.to_au()`: position in **AU**, velocity in **AU/day**.
//!   - Earth geocenter is derived from EMB and Moon using the Earth–Moon mass ratio.
//!
//! - **NAIF SPK/DAF backend (`naif`)**
//!   - **Query:** `Earth–Moon Barycenter (EMB)` relative to `Solar System Barycenter (SSB)`.
//!   - **Internal units:** kilometers; velocity in kilometers **per second**.
//!   - After calling `.to_au()`: position in **AU**, velocity in **AU/s**.
//!
//! > ⚠️ **Caveat:** the velocity unit differs between backends. If your pipeline needs
//! > a single unit (e.g. AU/day), convert NAIF velocities by multiplying by `86400.0`.
//! > If you require **geocentric Earth** when using SPK/DAF, convert EMB → Earth geocenter
//! > (or query a segment providing geocenter directly, if available).
//!
//! # Time scales
//! - `horizon`: input time is converted to **MJD(TT)** for interpolation on JD intervals.
//! - `naif`: input time is converted to **ET seconds** (SPICE TDB‑like) for SPK evaluation.
//!
//! # Submodules overview
//!
//! - [`download_jpl_file`](crate::jpl_ephem::download_jpl_file) — Ephemeris file resolution and retrieval (local/cache/remote).
//! - [`horizon`](crate::jpl_ephem::horizon) — Legacy DE reader and interpolator:
//!   * [`horizon::HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData) — top‑level loader and query interface,
//!   * [`horizon::horizon_records::HorizonRecord`](crate::jpl_ephem::horizon::horizon_records::HorizonRecord) — per‑interval Chebyshev coefficients,
//!   * [`horizon::horizon_ids::HorizonID`](crate::jpl_ephem::horizon::horizon_ids::HorizonID) — body/center identifiers,
//!   * [`horizon::horizon_version::JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion) — DE version to filename mapping,
//!   * [`horizon::interpolation_result::InterpResult`](crate::jpl_ephem::horizon::interpolation_result::InterpResult) — (pos, vel?, acc?) with `.to_au()`.
//! - [`naif`](crate::jpl_ephem::naif) — SPK/DAF reader and interpolator:
//!   * [`naif::NaifData`](crate::jpl_ephem::naif::naif_data::NaifData) — top‑level loader and query interface,
//!   * [`naif::naif_ids`](crate::jpl_ephem::naif::naif_ids) — NAIF ID enums (planetary barycenters, SSB, …),
//!   * [`naif::naif_version::NaifVersion`](crate::jpl_ephem::naif::naif_version::NaifVersion) — SPK file names by official DE label.
//!
//! # Example
//! ```no_run
//! use outfit::jpl::{download_jpl_file::EphemFileSource, JPLEphem};
//! use hifitime::Epoch;
//!
//! let eph = JPLEphem::new(&EphemFileSource::Auto)?;
//! let t = Epoch::from_tai_seconds(1_700_000_000.0);
//!
//! // Position always in AU. Velocity is AU/day (legacy DE) or AU/s (SPK/DAF).
//! let (r_au, v_opt) = eph.earth_ephemeris(&t, true);
//!
//! // Normalize NAIF velocities to AU/day if your pipeline expects it:
//! let v_au_per_day = v_opt.map(|v| v * 86400.0);
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! # See also
//! * [`horizon::HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris) — Earth (geocenter) w.r.t. Sun, Chebyshev on JD.
//! * [`naif::NaifData::ephemeris`](crate::jpl_ephem::naif::naif_data::NaifData::ephemeris) — EMB w.r.t. SSB, Chebyshev on ET seconds.
//! * [`hifitime::Epoch`] — conversions to MJD(TT) and ET seconds.

use download_jpl_file::{EphemFilePath, EphemFileSource};
use hifitime::Epoch;
use horizon::{horizon_data::HorizonData, horizon_ids::HorizonID};
use naif::{
    naif_data::NaifData,
    naif_ids::{planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds},
};
use nalgebra::Vector3;

use crate::outfit_errors::OutfitError;

pub mod download_jpl_file;
pub mod horizon;
pub mod naif;

/// Runtime‑selected ephemeris backend (legacy DE vs NAIF SPK/DAF).
///
/// # See also
/// * [`download_jpl_file::EphemFileSource`] – Source resolution policy.
/// * [`horizon::HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData) – Legacy DE binary reader/interpolator.
/// * [`naif::NaifData`](crate::jpl_ephem::naif::naif_data::NaifData) – SPK/DAF high‑level accessor.
#[derive(Debug, Clone)]
pub enum JPLEphem {
    HorizonFile(HorizonData),
    NaifFile(NaifData),
}

impl JPLEphem {
    /// Construct a [`JPLEphem`] from an ephemeris source policy.
    ///
    /// Resolves the concrete file path (download if requested), detects the
    /// format (legacy DE vs SPK/DAF), then loads the appropriate backend.
    ///
    /// # Errors
    /// Returns [`OutfitError`] if the ephemeris file cannot be resolved or parsed.
    ///
    /// # See also
    /// * [`download_jpl_file::EphemFilePath::get_ephemeris_file`]
    /// * [`horizon::HorizonData::read_horizon_file`](crate::jpl_ephem::horizon::horizon_data::HorizonData::read_horizon_file)
    /// * [`naif::NaifData::read_naif_file`](crate::jpl_ephem::naif::naif_data::NaifData::read_naif_file)
    pub fn new(file_source: &EphemFileSource) -> Result<Self, OutfitError> {
        let file_path = EphemFilePath::get_ephemeris_file(file_source)?;
        match file_path {
            EphemFilePath::JPLHorizon(..) => {
                let horizon_data = HorizonData::read_horizon_file(&file_path);
                Ok(JPLEphem::HorizonFile(horizon_data))
            }
            EphemFilePath::Naif(..) => {
                let naif_data = NaifData::read_naif_file(&file_path);
                Ok(JPLEphem::NaifFile(naif_data))
            }
        }
    }

    /// Return Earth state vectors at `ephem_time`.
    ///
    /// - **Legacy DE (`horizon`)**: `Earth` − `Sun`; after `.to_au()`: **AU** and **AU/day**.
    /// - **NAIF SPK/DAF (`naif`)**: `EMB` − `SSB`; after `.to_au()`: **AU** and **AU/s**.
    ///
    /// If a uniform velocity unit is required, convert NAIF velocities by `* 86400.0`.
    ///
    /// # Parameters
    /// * `ephem_time` — Observation epoch (`hifitime::Epoch`).
    /// * `compute_velocity` — Whether to compute and return the velocity.
    ///
    /// # Returns
    /// `(position_au, velocity_opt)` with velocity present only if requested.
    ///
    /// # See also
    /// * [`horizon::HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris)
    /// * [`naif::NaifData::ephemeris`](crate::jpl_ephem::naif::naif_data::NaifData::ephemeris)
    pub fn earth_ephemeris(
        &self,
        ephem_time: &Epoch,
        compute_velocity: bool,
    ) -> (Vector3<f64>, Option<Vector3<f64>>) {
        match self {
            JPLEphem::HorizonFile(horizon_data) => {
                let ephem_res = horizon_data
                    .ephemeris(
                        HorizonID::Earth,
                        HorizonID::Sun,
                        ephem_time.to_mjd_tt_days(),
                        compute_velocity,
                        false,
                    )
                    .to_au();
                (ephem_res.position, ephem_res.velocity)
            }
            JPLEphem::NaifFile(naif_data) => {
                let ephem_res = naif_data
                    .ephemeris(
                        NaifIds::PB(PlanetaryBary::EarthMoon),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        ephem_time.to_et_seconds(),
                    )
                    .to_au();
                (ephem_res.position, ephem_res.velocity)
            }
        }
    }
}
