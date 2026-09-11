//! Unified access to JPL ephemerides behind a single API.
//!
//! The enum [`JPLEphem`](crate::jpl_ephem::JPLEphem) wraps a planetary ephemeris
//! backend and exposes a common entry point to get Earth state vectors via
//! [`JPLEphem::earth_ephemeris`](crate::jpl_ephem::JPLEphem::earth_ephemeris) and
//! arbitrary body state vectors via
//! [`JPLEphem::body_ephemeris`](crate::jpl_ephem::JPLEphem::body_ephemeris).
//!
//! # Backends and feature flags
//!
//! The backend is chosen at compile time by Cargo feature. At least one must be
//! enabled; building with neither is a compile error.
//!
//! - **`ephem-builtin`** (in the `default` set) — the in-house reader:
//!   - **Legacy JPL DE binaries** (TTL/CNAM/IPT layout) via the
//!     [`horizon`](crate::jpl_ephem::horizon) backend, parsed into Chebyshev
//!     segments per body and interpolated on demand.
//!   - **NAIF SPK/DAF kernels** via the [`naif`](crate::jpl_ephem::naif) backend,
//!     parsed via a DAF header + summary/directory records into Chebyshev
//!     segments.
//!   Exposed as [`JPLEphem::HorizonFile`] and [`JPLEphem::NaifFile`].
//! - **`ephem-anise`** — the ANISE toolkit, a Rust reimplementation of the NAIF
//!   SPICE toolkit. It reads NAIF SPK kernels only, so it requires a `naif:DE###`
//!   source; a `horizon:` source is rejected. Exposed as [`JPLEphem::Anise`].
//!
//! [`JPLEphem::new`] resolves the source and builds the backend for the enabled
//! features — the ANISE backend when both are compiled. Use
//! [`JPLEphem::from_builtin`] / [`JPLEphem::from_anise`] to pick explicitly.
//!
//! Which bodies a backend can resolve also varies —
//! [`naif::naif_ids::NaifIds`](crate::jpl_ephem::naif::naif_ids::NaifIds)'s
//! "Ephemeris backend support" section lists exactly which variant works with
//! which backend and what happens (always a clean error, never a panic) when
//! it doesn't.
//!
//! The lightweight identifier and version modules
//! ([`naif::naif_ids`](crate::jpl_ephem::naif::naif_ids),
//! [`naif::naif_version`](crate::jpl_ephem::naif::naif_version),
//! [`horizon::horizon_version`](crate::jpl_ephem::horizon::horizon_version)) and
//! the ephemeris-file resolver
//! ([`download_jpl_file`](crate::jpl_ephem::download_jpl_file)) are available for
//! every feature combination.
//!
//! # Frames, centers, and units
//!
//! Every backend delivers states in the **equatorial mean J2000 (ICRF)** frame.
//! [`JPLEphem::earth_ephemeris`] and [`JPLEphem::body_ephemeris`] take an
//! [`EphemerisFrame`] argument: [`EphemerisFrame::Equatorial`] returns the state
//! unchanged, [`EphemerisFrame::Ecliptic`] rotates it to ecliptic mean J2000.
//! Both entry points return AU and AU/day, heliocentric (relative to the Sun),
//! whatever the backend.
//!
//! The Earth query returns the Earth **geocenter** relative to the Sun on the
//! legacy DE backend and on the ANISE backend; the built-in SPK reader returns
//! the Earth–Moon barycenter relative to the Solar System Barycenter. The
//! built-in DE backend derives the geocenter from the Earth–Moon barycenter and
//! the Moon using the Earth–Moon mass ratio.
//!
//! # Time scales
//!
//! The public entry points take a [`hifitime::Epoch`]; each backend converts it
//! internally — MJD(TT) for the legacy DE reader, ET/TDB seconds for the SPK
//! readers (built-in and ANISE).
//!
//! # Submodules overview
//!
//! - [`download_jpl_file`](crate::jpl_ephem::download_jpl_file) — Ephemeris file resolution and retrieval (local cache / download); every backend.
//! - [`horizon`](crate::jpl_ephem::horizon) — legacy DE reader (feature `ephem-builtin`); its `horizon_version` submodule is always compiled.
//! - [`naif`](crate::jpl_ephem::naif) — built-in SPK/DAF reader (feature `ephem-builtin`); its `naif_ids` and `naif_version` submodules are always compiled.
//! - [`anise_backend`](crate::jpl_ephem::anise_backend) — ANISE-backed SPK accessor (feature `ephem-anise`).
//!
//! # Example
//! ```rust, no_run
//! use outfit::jpl_ephem::{EphemerisFrame, JPLEphem};
//! use hifitime::Epoch;
//!
//! // `"naif:DE440"` works with every backend; the built-in reader also accepts
//! // `"horizon:DE440"`.
//! let eph: JPLEphem = "naif:DE440".try_into()?;
//! let t = Epoch::from_tai_seconds(1_700_000_000.0);
//!
//! // Position in AU, velocity in AU/day, equatorial mean J2000.
//! let (r_au, v_opt) = eph.earth_ephemeris(&t, EphemerisFrame::Equatorial, true);
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! # See also
//! * [`horizon::HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris) — Earth (geocenter) w.r.t. Sun, Chebyshev on JD.
//! * [`naif::NaifData::ephemeris`](crate::jpl_ephem::naif::naif_data::NaifData::ephemeris) — EMB w.r.t. SSB, Chebyshev on ET seconds.
//! * [`hifitime::Epoch`] — conversions to MJD(TT) and ET seconds.

#[cfg(not(any(feature = "ephem-builtin", feature = "ephem-anise")))]
compile_error!(
    "Outfit needs a planetary ephemeris backend: enable feature `ephem-builtin` (default) or `ephem-anise`."
);

use download_jpl_file::{EphemFilePath, EphemFileSource};
use hifitime::Epoch;
#[cfg(feature = "ephem-builtin")]
use horizon::{horizon_data::HorizonData, horizon_ids::HorizonID};
#[cfg(feature = "ephem-builtin")]
use naif::naif_data::NaifData;
use naif::naif_ids::NaifIds;
#[cfg(feature = "ephem-builtin")]
use naif::naif_ids::{
    planet_bary::PlanetaryBary, satellite_mass::SatelliteMassCenter,
    solar_system_bary::SolarSystemBary,
};
use nalgebra::Vector3;

use crate::constants::ROT_EQUMJ2000_TO_ECLMJ2000;
use crate::outfit_errors::OutfitError;

#[cfg(feature = "ephem-anise")]
pub mod anise_backend;
pub mod download_jpl_file;
pub mod horizon;
pub mod naif;

/// J2000 frame orientation in which an ephemeris state is expressed.
///
/// JPL DE and NAIF SPK data are stored in the equatorial ICRF/J2000 frame.
/// Selecting [`EphemerisFrame::Ecliptic`] rotates the returned state by the mean
/// obliquity of the ecliptic at J2000.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EphemerisFrame {
    /// Equatorial mean J2000 (ICRF), returned without rotation.
    Equatorial,
    /// Ecliptic mean J2000, obtained by rotating the equatorial state around the
    /// X-axis by the mean obliquity of the ecliptic at J2000.
    Ecliptic,
}

/// Rotate an equatorial mean J2000 `(position, velocity)` pair into ecliptic
/// mean J2000.
///
/// Applies [`ROT_EQUMJ2000_TO_ECLMJ2000`](crate::constants::ROT_EQUMJ2000_TO_ECLMJ2000),
/// a rotation of `-ε` around the X-axis where `ε` is the obliquity of the
/// ecliptic at J2000.
///
/// # Arguments
///
/// * `position_equ` – Position in equatorial mean J2000. Units: AU.
/// * `velocity_equ` – Velocity in equatorial mean J2000. Units: AU/day.
///
/// # Returns
///
/// `(position_ecl, velocity_ecl)` in ecliptic mean J2000, same units.
#[inline]
pub(crate) fn equ_state_to_ecl(
    position_equ: Vector3<f64>,
    velocity_equ: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    (
        ROT_EQUMJ2000_TO_ECLMJ2000 * position_equ,
        ROT_EQUMJ2000_TO_ECLMJ2000 * velocity_equ,
    )
}

/// Express an equatorial mean J2000 Earth state in the requested output frame.
///
/// The velocity is carried through untouched (still `None` when it was not
/// requested), so the "velocity present only if asked for" contract of
/// [`JPLEphem::earth_ephemeris`] is preserved.
///
/// # Arguments
///
/// * `frame` – requested output frame orientation.
/// * `position_equ` – Earth position in equatorial mean J2000. Units: AU.
/// * `velocity_equ` – optional Earth velocity in equatorial mean J2000. Units: AU/day.
///
/// # Returns
///
/// `(position, velocity)` in the requested frame: returned unchanged for
/// [`EphemerisFrame::Equatorial`], rotated by [`equ_state_to_ecl`] on both
/// components for [`EphemerisFrame::Ecliptic`].
#[inline]
pub(crate) fn select_earth_output(
    frame: EphemerisFrame,
    position_equ: Vector3<f64>,
    velocity_equ: Option<Vector3<f64>>,
) -> (Vector3<f64>, Option<Vector3<f64>>) {
    match frame {
        EphemerisFrame::Equatorial => (position_equ, velocity_equ),
        EphemerisFrame::Ecliptic => (
            ROT_EQUMJ2000_TO_ECLMJ2000 * position_equ,
            velocity_equ.map(|v| ROT_EQUMJ2000_TO_ECLMJ2000 * v),
        ),
    }
}

/// Express an equatorial mean J2000 body state in the requested output frame.
///
/// # Arguments
///
/// * `frame` – requested output frame orientation.
/// * `position_equ` – body position in equatorial mean J2000. Units: AU.
/// * `velocity_equ` – body velocity in equatorial mean J2000. Units: AU/day.
///
/// # Returns
///
/// `(position, velocity)` in the requested frame: returned unchanged for
/// [`EphemerisFrame::Equatorial`], rotated by [`equ_state_to_ecl`] for
/// [`EphemerisFrame::Ecliptic`].
#[inline]
pub(crate) fn select_body_output(
    frame: EphemerisFrame,
    position_equ: Vector3<f64>,
    velocity_equ: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    match frame {
        EphemerisFrame::Equatorial => (position_equ, velocity_equ),
        EphemerisFrame::Ecliptic => equ_state_to_ecl(position_equ, velocity_equ),
    }
}

/// Planetary ephemeris handle, holding whichever backend was selected at compile
/// time.
///
/// # Backends and feature flags
///
/// Which variants exist is chosen by Cargo feature:
///
/// * `ephem-builtin` (default) – the in-house reader, exposed as
///   [`JPLEphem::HorizonFile`] (legacy JPL DE binary) and [`JPLEphem::NaifFile`]
///   (NAIF SPK/DAF kernel).
/// * `ephem-anise` – the ANISE toolkit, exposed as [`JPLEphem::Anise`]; it reads
///   NAIF SPK kernels only, so a `horizon:` source is rejected.
///
/// Enabling both features compiles every variant; [`JPLEphem::new`] then builds
/// the ANISE backend, while [`JPLEphem::from_builtin`] still forces the in-house
/// one. Building with neither feature is a compile error.
///
/// Regardless of the backend, [`JPLEphem::earth_ephemeris`] and
/// [`JPLEphem::body_ephemeris`] return heliocentric states in AU and AU/day,
/// equatorial mean J2000 unless [`EphemerisFrame::Ecliptic`] is requested.
#[derive(Debug, Clone)]
pub enum JPLEphem {
    /// Legacy JPL DE binary reader (feature `ephem-builtin`).
    #[cfg(feature = "ephem-builtin")]
    HorizonFile(HorizonData),
    /// NAIF SPK/DAF kernel reader (feature `ephem-builtin`).
    #[cfg(feature = "ephem-builtin")]
    NaifFile(NaifData),
    /// ANISE-backed SPK kernel accessor (feature `ephem-anise`).
    #[cfg(feature = "ephem-anise")]
    Anise(anise_backend::AniseEphem),
}

impl JPLEphem {
    /// Construct a [`JPLEphem`] using the compile-time-selected backend.
    ///
    /// With a single backend feature enabled, that backend is used. With both
    /// `ephem-builtin` and `ephem-anise` enabled, the ANISE backend is selected;
    /// use [`JPLEphem::from_builtin`] to force the in-house reader.
    ///
    /// # Arguments
    ///
    /// * `source` – ephemeris source policy (backend token plus DE version), for
    ///   example the string form `"horizon:DE440"` or `"naif:DE440"`.
    ///
    /// # Returns
    ///
    /// A ready-to-query [`JPLEphem`]: the ephemeris file is resolved (downloaded
    /// into the local cache if missing) and parsed.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if the ephemeris file cannot be resolved,
    /// downloaded, or parsed, or if the `source` is not compatible with the
    /// active backend.
    pub fn new(source: impl Into<EphemFileSource>) -> Result<Self, OutfitError> {
        #[cfg(feature = "ephem-anise")]
        {
            Self::from_anise(source)
        }
        #[cfg(all(feature = "ephem-builtin", not(feature = "ephem-anise")))]
        {
            Self::from_builtin(source)
        }
    }

    /// Construct a [`JPLEphem`] backed by the in-house JPL DE / NAIF SPK reader.
    ///
    /// Resolves the concrete file path (downloading into the local cache when
    /// missing), then loads the reader matching the file kind: a legacy DE
    /// binary for a `horizon:` source, a NAIF SPK/DAF kernel for a `naif:` source.
    ///
    /// # Arguments
    ///
    /// * `source` – ephemeris source policy; both `horizon:` and `naif:` tokens
    ///   are accepted.
    ///
    /// # Returns
    ///
    /// A [`JPLEphem`] wrapping the parsed in-house backend.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if the ephemeris file cannot be resolved,
    /// downloaded, or parsed.
    #[cfg(feature = "ephem-builtin")]
    pub fn from_builtin(source: impl Into<EphemFileSource>) -> Result<Self, OutfitError> {
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
            #[cfg(feature = "ephem-anise")]
            EphemFilePath::MainBeltAsteroids(..) => Err(OutfitError::InvalidJPLEphemFileSource(
                "the in-house reader cannot read the main-belt asteroid supplementary kernel; \
                 use the ephem-anise backend"
                    .to_string(),
            )),
        }
    }

    /// Construct a [`JPLEphem`] backed by the ANISE toolkit.
    ///
    /// Resolves the concrete file path (downloading into the local cache when
    /// missing) and loads it as a NAIF SPK kernel.
    ///
    /// # Arguments
    ///
    /// * `source` – ephemeris source policy; only the `naif:` token is supported.
    ///   A `horizon:` (legacy DE binary) source is rejected.
    ///
    /// # Returns
    ///
    /// A [`JPLEphem`] wrapping an ANISE almanac loaded from the resolved kernel.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::InvalidJPLEphemFileSource`] for a legacy-DE source,
    /// or [`OutfitError`] if the kernel cannot be resolved, downloaded, or parsed.
    #[cfg(feature = "ephem-anise")]
    pub fn from_anise(source: impl Into<EphemFileSource>) -> Result<Self, OutfitError> {
        let file_path = EphemFilePath::get_ephemeris_file(&source.into())?;
        Ok(JPLEphem::Anise(anise_backend::AniseEphem::from_file_path(
            &file_path,
        )?))
    }

    /// Load the main-belt asteroid supplementary kernel into this handle, in
    /// place.
    ///
    /// Adds a second SPK kernel (`codes_300ast_20100725.bsp`, downloaded into
    /// the local cache on first use) covering 300 numbered main-belt
    /// asteroids, on top of whichever primary kernel this handle already has
    /// loaded. Once this call returns, `NaifIds::AST(_)` bodies (see
    /// [`crate::propagator::planet_gm::known_main_belt_asteroids`]) resolve
    /// through [`JPLEphem::body_ephemeris`] like any other body.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::InvalidJPLEphemFileSource`] if this handle is
    /// not backed by the ANISE toolkit (see [`JPLEphem::from_anise`]) — the
    /// in-house reader has no code path for supplementary kernels — or
    /// [`OutfitError`] if the supplementary kernel cannot be resolved,
    /// downloaded, or parsed.
    #[cfg(feature = "ephem-anise")]
    pub fn with_main_belt_asteroids(&mut self) -> Result<(), OutfitError> {
        match self {
            JPLEphem::Anise(anise) => {
                let asteroids_path =
                    EphemFilePath::get_ephemeris_file(&EphemFileSource::MainBeltAsteroids)?;
                anise.with_supplementary_kernel(&asteroids_path)
            }
            #[cfg(feature = "ephem-builtin")]
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "the main-belt asteroid supplementary kernel requires the ephem-anise backend"
                    .to_string(),
            )),
        }
    }

    /// Return Earth's heliocentric state at `ephem_time`.
    ///
    /// The state is heliocentric (Earth geocenter relative to the Sun for the
    /// legacy DE and the ANISE backends; Earth–Moon barycenter relative to the
    /// Solar System Barycenter for the built-in SPK reader), in AU and AU/day. It
    /// is returned in equatorial mean J2000 for [`EphemerisFrame::Equatorial`]
    /// and rotated to ecliptic mean J2000 for [`EphemerisFrame::Ecliptic`].
    ///
    /// # Arguments
    ///
    /// * `ephem_time` – evaluation epoch.
    /// * `frame` – frame orientation of the returned state.
    /// * `compute_velocity` – whether to also compute and return the velocity.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` with `position` in AU and `velocity` in AU/day,
    /// heliocentric, in the requested frame; `velocity` is `Some` only when
    /// `compute_velocity` is `true`.
    ///
    /// # Panics
    ///
    /// Every loaded kernel is expected to cover Earth for any epoch a real
    /// ephemeris file spans, so this call is treated as infallible; it panics
    /// if the active backend cannot resolve Earth at `ephem_time` (for example
    /// an epoch outside the kernel's coverage). Use [`JPLEphem::body_ephemeris`]
    /// for a fallible query.
    pub fn earth_ephemeris(
        &self,
        ephem_time: &Epoch,
        frame: EphemerisFrame,
        compute_velocity: bool,
    ) -> (Vector3<f64>, Option<Vector3<f64>>) {
        let (position_equ, velocity_equ) =
            self.earth_ephemeris_equatorial(ephem_time, compute_velocity);
        select_earth_output(frame, position_equ, velocity_equ)
    }

    /// Heliocentric equatorial mean J2000 Earth state, in AU and AU/day.
    ///
    /// Backend-specific worker shared by every [`EphemerisFrame`] branch of
    /// [`JPLEphem::earth_ephemeris`].
    ///
    /// # Arguments
    ///
    /// * `ephem_time` – evaluation epoch.
    /// * `compute_velocity` – whether to also compute and return the velocity.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` heliocentric, equatorial mean J2000; `velocity` is
    /// `Some` only when `compute_velocity` is `true`.
    fn earth_ephemeris_equatorial(
        &self,
        ephem_time: &Epoch,
        compute_velocity: bool,
    ) -> (Vector3<f64>, Option<Vector3<f64>>) {
        match self {
            #[cfg(feature = "ephem-builtin")]
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
            #[cfg(feature = "ephem-builtin")]
            JPLEphem::NaifFile(naif_data) => {
                let ephem_res = naif_data
                    .ephemeris(
                        NaifIds::PB(PlanetaryBary::EarthMoon),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        ephem_time.to_et_seconds(),
                    )
                    .unwrap_or_else(|err| panic!("NAIF Earth ephemeris lookup failed: {err}"))
                    .to_au();
                (ephem_res.position, ephem_res.velocity.map(|v| v * 86400.0)) // Convert from AU/s to AU/day
            }
            #[cfg(feature = "ephem-anise")]
            JPLEphem::Anise(anise) => {
                anise.earth_ephemeris_equatorial(ephem_time, compute_velocity)
            }
        }
    }

    /// Consume this handle and return the inner legacy DE binary reader.
    ///
    /// # Returns
    ///
    /// The wrapped reader when this handle holds a legacy DE binary backend.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::InvalidJPLEphemFileSource`] when this handle holds a
    /// different backend.
    #[cfg(feature = "ephem-builtin")]
    pub fn try_into_horizon(self) -> Result<HorizonData, OutfitError> {
        match self {
            JPLEphem::HorizonFile(horizon_data) => Ok(horizon_data),
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "Expected a JPL Horizon source".to_string(),
            )),
        }
    }

    /// Consume this handle and return the inner built-in NAIF SPK reader.
    ///
    /// # Returns
    ///
    /// The wrapped reader when this handle holds a built-in NAIF SPK backend.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::InvalidJPLEphemFileSource`] when this handle holds a
    /// different backend.
    #[cfg(feature = "ephem-builtin")]
    pub fn try_into_naif(self) -> Result<NaifData, OutfitError> {
        match self {
            JPLEphem::NaifFile(naif_data) => Ok(naif_data),
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "Expected a NAIF source".to_string(),
            )),
        }
    }

    /// Consume this handle and return the inner ANISE backend.
    ///
    /// # Returns
    ///
    /// The wrapped [`anise_backend::AniseEphem`] when this handle holds the ANISE
    /// backend.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::InvalidJPLEphemFileSource`] when this handle holds a
    /// different backend.
    #[cfg(feature = "ephem-anise")]
    pub fn try_into_anise(self) -> Result<anise_backend::AniseEphem, OutfitError> {
        match self {
            JPLEphem::Anise(anise) => Ok(anise),
            #[cfg(feature = "ephem-builtin")]
            _ => Err(OutfitError::InvalidJPLEphemFileSource(
                "Expected an ANISE source".to_string(),
            )),
        }
    }

    /// Return the heliocentric position and velocity of `body` at `epoch`.
    ///
    /// The result is normalised to a single unit system regardless of backend:
    /// position in AU and velocity in AU/day, relative to the Sun. It is returned
    /// in equatorial mean J2000 for [`EphemerisFrame::Equatorial`] and rotated to
    /// ecliptic mean J2000 for [`EphemerisFrame::Ecliptic`].
    ///
    /// # Arguments
    ///
    /// * `body` – perturbing body to look up.
    /// * `epoch` – evaluation epoch.
    /// * `frame` – frame orientation of the returned state.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` with `position` in AU and `velocity` in AU/day,
    /// heliocentric, in the requested frame.
    ///
    /// # Errors
    ///
    /// - [`OutfitError::EphemerisBodyNotSupported`] if `body` cannot be resolved
    ///   by the active backend (see [`NaifIds`]'s "Ephemeris backend support"
    ///   section for which variants each backend accepts) — this never panics,
    ///   including for a body absent from the loaded kernel.
    /// - Any error propagated by the underlying ephemeris query.
    pub fn body_ephemeris(
        &self,
        body: NaifIds,
        epoch: &Epoch,
        frame: EphemerisFrame,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        let (position_equ, velocity_equ) = self.body_ephemeris_equatorial(body, epoch)?;
        Ok(select_body_output(frame, position_equ, velocity_equ))
    }

    /// Heliocentric equatorial mean J2000 `(position, velocity)` of `body`, in
    /// AU and AU/day.
    ///
    /// This is the backend-specific query shared by every [`EphemerisFrame`]
    /// branch of [`JPLEphem::body_ephemeris`].
    ///
    /// # Arguments
    ///
    /// * `body` – perturbing body to look up.
    /// * `epoch` – evaluation epoch.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` heliocentric, equatorial mean J2000.
    ///
    /// # Errors
    ///
    /// - [`OutfitError::EphemerisBodyNotSupported`] if `body` cannot be resolved
    ///   by the active backend.
    fn body_ephemeris_equatorial(
        &self,
        body: NaifIds,
        epoch: &Epoch,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        match self {
            #[cfg(feature = "ephem-anise")]
            JPLEphem::Anise(anise) => anise.body_ephemeris_equatorial(body, epoch),
            #[cfg(feature = "ephem-builtin")]
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
                // `to_au()` already yields km/day → AU/day for the Horizon
                // backend; there is no per-second conversion here (unlike the
                // NAIF branch below, whose raw velocity is AU/s).
                let velocity = ephem_res.velocity.unwrap_or_else(Vector3::zeros);
                Ok((ephem_res.position, velocity))
            }
            #[cfg(feature = "ephem-builtin")]
            JPLEphem::NaifFile(naif_data) => {
                let et = epoch.to_et_seconds();
                // Query body w.r.t. SSB. A body this reader's kernel does not
                // carry (e.g. a numbered asteroid) surfaces as
                // `EphemerisBodyNotSupported` here rather than panicking.
                let body_res = naif_data
                    .ephemeris(body, NaifIds::SSB(SolarSystemBary::SSB), et)?
                    .to_au();
                // Query Sun w.r.t. SSB
                let sun_res = naif_data
                    .ephemeris(
                        NaifIds::SSB(SolarSystemBary::Sun),
                        NaifIds::SSB(SolarSystemBary::SSB),
                        et,
                    )?
                    .to_au();
                // Heliocentric = body_ssb - sun_ssb; convert AU/s → AU/day
                let pos = body_res.position - sun_res.position;
                let vel = (body_res.velocity.unwrap_or_else(Vector3::zeros)
                    - sun_res.velocity.unwrap_or_else(Vector3::zeros))
                    * 86400.0;
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
#[cfg(feature = "ephem-builtin")]
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

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(all(test, feature = "ephem-builtin"))]
mod jpl_ephem_tests {
    use super::*;
    use crate::test_fixture::{JPL_EPHEM_HORIZON, JPL_EPHEM_NAIF};
    use approx::assert_abs_diff_eq;
    use hifitime::TimeScale;
    use proptest::prelude::*;

    /// MJD-TT epoch comfortably inside the DE440 coverage used by the fixtures.
    fn epoch_tt(mjd_tt: f64) -> Epoch {
        Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT)
    }

    /// Central finite-difference estimate of a body's heliocentric velocity
    /// (AU/day) from two [`JPLEphem::body_ephemeris`] position samples.
    ///
    /// # Arguments
    ///
    /// * `jpl` – ephemeris backend to query.
    /// * `body` – perturbing body.
    /// * `mjd_tt` – central epoch (MJD TT).
    /// * `h_days` – half-step in days; must be strictly positive.
    ///
    /// # Returns
    ///
    /// `(r(t + h) − r(t − h)) / (2 h)` in AU/day.
    fn velocity_by_central_difference(
        jpl: &JPLEphem,
        body: NaifIds,
        mjd_tt: f64,
        h_days: f64,
    ) -> Vector3<f64> {
        let r_plus = jpl
            .body_ephemeris(body, &epoch_tt(mjd_tt + h_days), EphemerisFrame::Equatorial)
            .unwrap()
            .0;
        let r_minus = jpl
            .body_ephemeris(body, &epoch_tt(mjd_tt - h_days), EphemerisFrame::Equatorial)
            .unwrap()
            .0;
        (r_plus - r_minus) / (2.0 * h_days)
    }

    /// On the Horizon backend `body_ephemeris(EarthMoon)` and `earth_ephemeris`
    /// resolve the exact same query (`naif_to_horizon_id(EarthMoon) = Earth`),
    /// so both position and velocity must match bit-for-bit. This locks the
    /// AU/day contract: a stray `* 86400.0` would blow the velocity check up by
    /// almost five orders of magnitude.
    #[test]
    fn horizon_body_ephemeris_velocity_agrees_with_earth_ephemeris() {
        let jpl = &*JPL_EPHEM_HORIZON;
        let epoch = epoch_tt(59_000.0);

        let (body_pos, body_vel) = jpl
            .body_ephemeris(
                NaifIds::PB(PlanetaryBary::EarthMoon),
                &epoch,
                EphemerisFrame::Equatorial,
            )
            .unwrap();
        let (earth_pos, earth_vel) = jpl.earth_ephemeris(&epoch, EphemerisFrame::Equatorial, true);
        let earth_vel = earth_vel.expect("compute_velocity = true must return a velocity");

        assert_abs_diff_eq!(body_pos, earth_pos, epsilon = 1e-15);
        assert_abs_diff_eq!(body_vel, earth_vel, epsilon = 1e-15);
    }

    /// On the Horizon backend the returned velocity must be the time derivative
    /// of the returned position, i.e. genuinely in AU/day.
    #[test]
    fn horizon_body_ephemeris_velocity_is_the_position_derivative() {
        let bodies = [
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::PB(PlanetaryBary::Mars),
            NaifIds::PB(PlanetaryBary::Jupiter),
            NaifIds::PB(PlanetaryBary::Saturn),
        ];
        let h = 0.25_f64;
        let mjd_tt = 59_500.0_f64;

        for &body in &bodies {
            let v = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch_tt(mjd_tt), EphemerisFrame::Equatorial)
                .unwrap()
                .1;
            let v_fd = velocity_by_central_difference(&JPL_EPHEM_HORIZON, body, mjd_tt, h);
            let rel_err = (v - v_fd).norm() / v.norm();
            assert!(
                rel_err < 1e-5,
                "{body:?}: velocity {v:?} vs finite difference {v_fd:?}, rel err {rel_err:e}"
            );
        }
    }

    /// Position agrees between the two DE440 backends (they read the same
    /// ephemeris; only the binary layout / interpolation granularity differ).
    #[test]
    fn body_ephemeris_position_agrees_between_horizon_and_naif() {
        let epoch = epoch_tt(59_777.0);
        for body in [
            NaifIds::PB(PlanetaryBary::Mars),
            NaifIds::PB(PlanetaryBary::Jupiter),
            NaifIds::PB(PlanetaryBary::Saturn),
        ] {
            let h_pos = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                .unwrap()
                .0;
            let n_pos = JPL_EPHEM_NAIF
                .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                .unwrap()
                .0;
            assert!(
                (h_pos - n_pos).norm() < 1e-6,
                "{body:?}: position disagreement {:e} AU",
                (h_pos - n_pos).norm()
            );
        }
    }

    /// On the NAIF backend the returned velocity must be the time derivative
    /// of the returned position, i.e. genuinely in AU/day — same contract as
    /// [`horizon_body_ephemeris_velocity_is_the_position_derivative`].
    ///
    /// Regression test for a former bug: `EphemerisRecord::interpolate` scaled
    /// the Chebyshev derivative by `2.0 / radius` instead of the correct
    /// chain-rule factor `1.0 / radius` (velocities 2× too fast), and both
    /// `body_ephemeris` / `earth_ephemeris` then *divided* by 86400 where they
    /// must *multiply* (AU/s → AU/day). Net effect: NAIF `body_ephemeris`
    /// velocity was `2 / 86400²` of the true value — a factor of about
    /// 3.7 billion.
    #[test]
    fn naif_body_ephemeris_velocity_is_the_position_derivative() {
        let bodies = [
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::PB(PlanetaryBary::Mars),
            NaifIds::PB(PlanetaryBary::Jupiter),
            NaifIds::PB(PlanetaryBary::Saturn),
        ];
        let h = 0.25_f64;
        let mjd_tt = 59_500.0_f64;

        for &body in &bodies {
            let v = JPL_EPHEM_NAIF
                .body_ephemeris(body, &epoch_tt(mjd_tt), EphemerisFrame::Equatorial)
                .unwrap()
                .1;
            let v_fd = velocity_by_central_difference(&JPL_EPHEM_NAIF, body, mjd_tt, h);
            let rel_err = (v - v_fd).norm() / v.norm();
            assert!(
                rel_err < 1e-5,
                "{body:?}: velocity {v:?} vs finite difference {v_fd:?}, rel err {rel_err:e}"
            );
        }
    }

    proptest! {
        /// Property: over the fixture coverage window and for any slow-moving
        /// planet, the Horizon-backend velocity is the central-difference
        /// derivative of the position (AU/day contract).
        #[test]
        fn horizon_velocity_matches_central_difference(
            mjd_tt in 55_000.0_f64..62_000.0,
            body_idx in 0usize..4,
        ) {
            let body = [
                NaifIds::PB(PlanetaryBary::Mars),
                NaifIds::PB(PlanetaryBary::Jupiter),
                NaifIds::PB(PlanetaryBary::Saturn),
                NaifIds::PB(PlanetaryBary::Uranus),
            ][body_idx];

            let v = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch_tt(mjd_tt), EphemerisFrame::Equatorial)
                .unwrap()
                .1;
            let v_fd = velocity_by_central_difference(&JPL_EPHEM_HORIZON, body, mjd_tt, 0.25);
            prop_assert!((v - v_fd).norm() / v.norm() < 1e-5);
        }

        /// Same property as [`horizon_velocity_matches_central_difference`],
        /// on the NAIF backend.
        #[test]
        fn naif_velocity_matches_central_difference(
            mjd_tt in 55_000.0_f64..62_000.0,
            body_idx in 0usize..4,
        ) {
            let body = [
                NaifIds::PB(PlanetaryBary::Mars),
                NaifIds::PB(PlanetaryBary::Jupiter),
                NaifIds::PB(PlanetaryBary::Saturn),
                NaifIds::PB(PlanetaryBary::Uranus),
            ][body_idx];

            let v = JPL_EPHEM_NAIF
                .body_ephemeris(body, &epoch_tt(mjd_tt), EphemerisFrame::Equatorial)
                .unwrap()
                .1;
            let v_fd = velocity_by_central_difference(&JPL_EPHEM_NAIF, body, mjd_tt, 0.25);
            prop_assert!((v - v_fd).norm() / v.norm() < 1e-5);
        }
    }

    // ── Reference-frame selection ────────────────────────────────────────────

    /// All non-Sun perturbing bodies used by the frame tests.
    const FRAME_TEST_BODIES: [NaifIds; 4] = [
        NaifIds::PB(PlanetaryBary::EarthMoon),
        NaifIds::PB(PlanetaryBary::Mars),
        NaifIds::PB(PlanetaryBary::Jupiter),
        NaifIds::PB(PlanetaryBary::Saturn),
    ];

    /// `equ_state_to_ecl` is exactly the `ROT_EQUMJ2000_TO_ECLMJ2000` product on
    /// both components.
    #[test]
    fn equ_state_to_ecl_matches_reference() {
        let position = Vector3::new(1.3, -0.7, 0.42);
        let velocity = Vector3::new(-0.011, 0.008, 0.003);
        let (pos_ecl, vel_ecl) = equ_state_to_ecl(position, velocity);
        assert_abs_diff_eq!(
            pos_ecl,
            ROT_EQUMJ2000_TO_ECLMJ2000 * position,
            epsilon = 1e-15
        );
        assert_abs_diff_eq!(
            vel_ecl,
            ROT_EQUMJ2000_TO_ECLMJ2000 * velocity,
            epsilon = 1e-15
        );
    }

    /// The heliocentric position of the Earth–Moon barycenter defines the
    /// ecliptic plane: sampled over a full year, its `z` component never leaves a
    /// ~10⁻⁴ AU band in ecliptic mean J2000, whereas in equatorial mean J2000 it
    /// swings up to ~0.4 AU (it also crosses zero twice a year, near the
    /// equinoxes). This discriminates the two frames on both backends.
    #[test]
    fn earth_moon_barycenter_ecliptic_z_is_small() {
        let body = NaifIds::PB(PlanetaryBary::EarthMoon);
        for (label, jpl) in [("horizon", &*JPL_EPHEM_HORIZON), ("naif", &*JPL_EPHEM_NAIF)] {
            let mut max_z_ecl = 0.0_f64;
            let mut max_z_equ = 0.0_f64;
            for step in 0..24 {
                let epoch = epoch_tt(58_900.0 + step as f64 * 15.5);
                max_z_ecl = max_z_ecl.max(
                    jpl.body_ephemeris(body, &epoch, EphemerisFrame::Ecliptic)
                        .unwrap()
                        .0[2]
                        .abs(),
                );
                max_z_equ = max_z_equ.max(
                    jpl.body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                        .unwrap()
                        .0[2]
                        .abs(),
                );
            }
            assert!(
                max_z_ecl < 3e-3,
                "{label}: max ecliptic |z| = {max_z_ecl} AU over a year (expected ~1e-4)"
            );
            assert!(
                max_z_equ > 0.35,
                "{label}: max equatorial |z| = {max_z_equ} AU over a year (expected ~0.4)"
            );
        }
    }

    /// On the Horizon backend, `earth_ephemeris` and `body_ephemeris(EarthMoon)`
    /// resolve the same query, so they must agree in every frame.
    #[test]
    fn earth_ephemeris_frame_is_consistent_with_body_ephemeris() {
        let jpl = &*JPL_EPHEM_HORIZON;
        let epoch = epoch_tt(59_321.0);
        for frame in [EphemerisFrame::Equatorial, EphemerisFrame::Ecliptic] {
            let earth_pos = jpl.earth_ephemeris(&epoch, frame, true).0;
            let body_pos = jpl
                .body_ephemeris(NaifIds::PB(PlanetaryBary::EarthMoon), &epoch, frame)
                .unwrap()
                .0;
            assert_abs_diff_eq!(earth_pos, body_pos, epsilon = 1e-15);
        }
    }

    proptest! {
        /// The ecliptic state is the equatorial state left-multiplied by
        /// `ROT_EQUMJ2000_TO_ECLMJ2000`, for position and velocity alike.
        #[test]
        fn body_ephemeris_ecliptic_is_equatorial_rotated(
            mjd_tt in 55_000.0_f64..62_000.0,
            body_idx in 0usize..4,
        ) {
            let body = FRAME_TEST_BODIES[body_idx];
            let epoch = epoch_tt(mjd_tt);
            let (pos_equ, vel_equ) = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                .unwrap();
            let (pos_ecl, vel_ecl) = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch, EphemerisFrame::Ecliptic)
                .unwrap();
            prop_assert!((pos_ecl - ROT_EQUMJ2000_TO_ECLMJ2000 * pos_equ).norm() < 1e-13);
            prop_assert!((vel_ecl - ROT_EQUMJ2000_TO_ECLMJ2000 * vel_equ).norm() < 1e-13);
        }

        /// The frame rotation is about the X-axis, so it preserves the `x`
        /// component and the vector norm.
        #[test]
        fn frame_choice_preserves_x_and_norm(
            mjd_tt in 55_000.0_f64..62_000.0,
            body_idx in 0usize..4,
        ) {
            let body = FRAME_TEST_BODIES[body_idx];
            let epoch = epoch_tt(mjd_tt);
            let pos_equ = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                .unwrap()
                .0;
            let pos_ecl = JPL_EPHEM_HORIZON
                .body_ephemeris(body, &epoch, EphemerisFrame::Ecliptic)
                .unwrap()
                .0;
            prop_assert!((pos_ecl[0] - pos_equ[0]).abs() < 1e-12);
            prop_assert!((pos_ecl.norm() - pos_equ.norm()).abs() < 1e-12 * pos_equ.norm());
        }
    }

    /// `with_main_belt_asteroids` requires an ANISE-backed handle; called on
    /// the in-house reader it must return a clean error, never panic. Only
    /// reachable when both backends are compiled in — with `ephem-anise`
    /// alone, `JPLEphem` has no other variant to call it on.
    #[test]
    #[cfg(all(feature = "ephem-builtin", feature = "ephem-anise"))]
    fn with_main_belt_asteroids_rejects_the_builtin_backend() {
        let source: EphemFileSource = "naif:DE440"
            .try_into()
            .expect("failed to parse JPL ephemeris source");
        let mut jpl =
            JPLEphem::from_builtin(source).expect("failed to load the in-house NAIF reader");
        let err = jpl.with_main_belt_asteroids().unwrap_err();
        assert!(matches!(err, OutfitError::InvalidJPLEphemFileSource(_)));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Frame-selection helper tests (backend-independent)
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod frame_selection_tests {
    use super::*;
    use proptest::prelude::*;

    /// `Equatorial` leaves both components untouched.
    #[test]
    fn equatorial_is_identity() {
        let p = Vector3::new(1.3, -0.7, 0.42);
        let v = Vector3::new(-0.011, 0.008, 0.003);

        let (pe, ve) = select_body_output(EphemerisFrame::Equatorial, p, v);
        assert_eq!(pe, p);
        assert_eq!(ve, v);

        let (pe, ve) = select_earth_output(EphemerisFrame::Equatorial, p, Some(v));
        assert_eq!(pe, p);
        assert_eq!(ve, Some(v));

        let (pe, ve) = select_earth_output(EphemerisFrame::Equatorial, p, None);
        assert_eq!(pe, p);
        assert_eq!(ve, None);
    }

    /// `Ecliptic` applies exactly [`equ_state_to_ecl`] to both components.
    #[test]
    fn ecliptic_matches_reference_rotation() {
        let p = Vector3::new(0.9, 0.31, -0.57);
        let v = Vector3::new(0.004, -0.017, 0.006);
        let (rp, rv) = equ_state_to_ecl(p, v);

        let (pe, ve) = select_body_output(EphemerisFrame::Ecliptic, p, v);
        assert_eq!(pe, rp);
        assert_eq!(ve, rv);

        let (pe, ve) = select_earth_output(EphemerisFrame::Ecliptic, p, Some(v));
        assert_eq!(pe, rp);
        assert_eq!(ve, Some(rv));
    }

    /// The `compute_velocity` contract survives frame selection: a `None`
    /// velocity stays `None`, a `Some` stays `Some`.
    #[test]
    fn earth_velocity_option_is_preserved() {
        let p = Vector3::new(1.0, 2.0, 3.0);
        assert!(select_earth_output(EphemerisFrame::Ecliptic, p, None)
            .1
            .is_none());
        assert!(select_earth_output(EphemerisFrame::Equatorial, p, None)
            .1
            .is_none());
        let v = Vector3::new(0.1, 0.2, 0.3);
        assert!(select_earth_output(EphemerisFrame::Ecliptic, p, Some(v))
            .1
            .is_some());
    }

    proptest! {
        /// The ecliptic rotation is about the X axis: it preserves the `x`
        /// component and the Euclidean norm of the position.
        #[test]
        fn ecliptic_preserves_x_and_norm(
            x in -50.0_f64..50.0,
            y in -50.0_f64..50.0,
            z in -50.0_f64..50.0,
        ) {
            let p = Vector3::new(x, y, z);
            let v = Vector3::new(0.0, 0.0, 0.0);
            let (pe, _) = select_body_output(EphemerisFrame::Ecliptic, p, v);
            prop_assert!((pe[0] - p[0]).abs() < 1e-12);
            prop_assert!((pe.norm() - p.norm()).abs() < 1e-12 * (1.0 + p.norm()));
        }

        /// `select_earth_output` and `select_body_output` agree on the position
        /// and on the velocity whenever the earth velocity is present.
        #[test]
        fn earth_and_body_helpers_agree(
            x in -10.0_f64..10.0,
            y in -10.0_f64..10.0,
            z in -10.0_f64..10.0,
            vx in -1.0_f64..1.0,
            vy in -1.0_f64..1.0,
            vz in -1.0_f64..1.0,
        ) {
            let p = Vector3::new(x, y, z);
            let v = Vector3::new(vx, vy, vz);
            for frame in [EphemerisFrame::Equatorial, EphemerisFrame::Ecliptic] {
                let (bp, bv) = select_body_output(frame, p, v);
                let (ep, ev) = select_earth_output(frame, p, Some(v));
                prop_assert_eq!(bp, ep);
                prop_assert_eq!(Some(bv), ev);
            }
        }
    }
}
