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
//! ```rust, no_run
//! use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem, horizon::horizon_version::JPLHorizonVersion};
//! use hifitime::Epoch;
//!
//! let eph = JPLEphem::new(&EphemFileSource::JPLHorizon(JPLHorizonVersion::DE440))?;
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
    naif_ids::{
        planet_bary::PlanetaryBary, satellite_mass::SatelliteMassCenter,
        solar_system_bary::SolarSystemBary, NaifIds,
    },
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
    pub fn new(source: impl Into<EphemFileSource>) -> Result<Self, OutfitError> {
        let file_path = EphemFilePath::get_ephemeris_file(&source.into())?;
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
    ///
    /// # Parameters
    /// * `ephem_time` — Observation epoch (`hifitime::Epoch`).
    /// * `compute_velocity` — Whether to compute and return the velocity.
    ///
    /// # Returns
    /// `(position_au, velocity_opt)` with velocity present only if requested.
    ///     - `position_au`: Earth position in **AU** (geocenter for `horizon`, EMB for `naif`).
    ///     - `velocity_opt`: Earth velocity in **AU/day**.
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
                (ephem_res.position, ephem_res.velocity.map(|v| v / 86400.0)) // Convert from AU/s to AU/day
            }
        }
    }

    pub fn try_into_horizon(self) -> Result<HorizonData, OutfitError> {
        match self {
            JPLEphem::HorizonFile(horizon_data) => Ok(horizon_data),
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "Expected a JPL Horizon source".to_string(),
            )),
        }
    }

    pub fn try_into_naif(self) -> Result<NaifData, OutfitError> {
        match self {
            JPLEphem::NaifFile(naif_data) => Ok(naif_data),
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "Expected a NAIF source".to_string(),
            )),
        }
    }

    /// Return heliocentric position and velocity (AU, AU/day) of `body` at `epoch`.
    ///
    /// Works for both backends:
    /// - **Horizon**: maps `NaifIds` to the corresponding `HorizonID` and queries relative to `Sun`.
    /// - **NAIF**: queries `(body, SSB)` then subtracts `(Sun, SSB)` to obtain heliocentric state.
    ///
    /// # Errors
    /// Returns [`OutfitError::EphemerisBodyNotSupported`] if the body cannot be mapped to
    /// the active backend.
    pub fn body_ephemeris(
        &self,
        body: NaifIds,
        epoch: &Epoch,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        match self {
            JPLEphem::HorizonFile(horizon_data) => {
                let horizon_id = naif_to_horizon_id(body)?;
                let ephem_res = horizon_data
                    .ephemeris(
                        horizon_id,
                        HorizonID::Sun,
                        epoch.to_mjd_tt_days(),
                        true,
                        false,
                    )
                    .to_au();
                let vel = ephem_res.velocity.unwrap_or_else(Vector3::zeros) * 86400.0; // AU/s → AU/day
                Ok((ephem_res.position, vel))
            }
            JPLEphem::NaifFile(naif_data) => {
                let et = epoch.to_et_seconds();
                // Query body w.r.t. SSB
                let body_res = naif_data
                    .ephemeris(body, NaifIds::SSB(SolarSystemBary::SSB), et)
                    .to_au();
                // Query Sun w.r.t. SSB
                let sun_res = naif_data
                    .ephemeris(
                        NaifIds::SSB(SolarSystemBary::Sun),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        et,
                    )
                    .to_au();
                // Heliocentric = body_ssb - sun_ssb; convert AU/s → AU/day
                let pos = body_res.position - sun_res.position;
                let vel = (body_res.velocity.unwrap_or_else(Vector3::zeros)
                    - sun_res.velocity.unwrap_or_else(Vector3::zeros))
                    / 86400.0;
                Ok((pos, vel))
            }
        }
    }
}

impl TryFrom<&str> for JPLEphem {
    type Error = OutfitError;
    fn try_from(s: &str) -> Result<Self, Self::Error> {
        let source = EphemFileSource::try_from(s)?;
        JPLEphem::new(source)
    }
}

impl TryFrom<String> for JPLEphem {
    type Error = OutfitError;
    fn try_from(s: String) -> Result<Self, Self::Error> {
        JPLEphem::try_from(s.as_str())
    }
}

/// Map a [`NaifIds`] to the corresponding [`HorizonID`] for the Horizon backend.
///
/// Only bodies that are stored in JPL Horizon DE files are supported.
/// Returns [`OutfitError::EphemerisBodyNotSupported`] for anything else.
fn naif_to_horizon_id(body: NaifIds) -> Result<HorizonID, OutfitError> {
    match body {
        NaifIds::SSB(SolarSystemBary::Sun) => Ok(HorizonID::Sun),
        NaifIds::SSB(SolarSystemBary::SSB) => Err(OutfitError::EphemerisBodyNotSupported(
            "Solar System Barycenter is not a physical body".to_string(),
        )),
        NaifIds::PB(PlanetaryBary::Mercury) => Ok(HorizonID::Mercury),
        NaifIds::PB(PlanetaryBary::Venus) => Ok(HorizonID::Venus),
        NaifIds::PB(PlanetaryBary::EarthMoon) => Ok(HorizonID::Earth),
        NaifIds::PB(PlanetaryBary::Mars) => Ok(HorizonID::Mars),
        NaifIds::PB(PlanetaryBary::Jupiter) => Ok(HorizonID::Jupiter),
        NaifIds::PB(PlanetaryBary::Saturn) => Ok(HorizonID::Saturn),
        NaifIds::PB(PlanetaryBary::Uranus) => Ok(HorizonID::Uranus),
        NaifIds::PB(PlanetaryBary::Neptune) => Ok(HorizonID::Neptune),
        NaifIds::PB(PlanetaryBary::Pluto) => Ok(HorizonID::Pluto),
        NaifIds::SMC(SatelliteMassCenter::Moon) => Ok(HorizonID::Moon),
        other => Err(OutfitError::EphemerisBodyNotSupported(format!(
            "{other} is not available in the Horizon backend"
        ))),
    }
}
