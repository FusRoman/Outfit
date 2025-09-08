//! # Observer & Site Geometry (top-level module)
//!
//! This module gathers **observer/site handling** and associated geometry utilities used in
//! orbit determination. It provides:
//!
//! - A robust [`Observer`](crate::observers::Observer) type storing **geocentric parallax coordinates** (ρ·cosφ, ρ·sinφ),
//!   geodetic longitude, optional astrometric accuracies, and **precomputed body-fixed** position
//!   and velocity vectors.
//! - High-level routines to compute the observer’s **geocentric PV** in the ecliptic J2000 frame
//!   ([`Observer::pvobs`](crate::observers::Observer::pvobs)) and its **heliocentric equatorial position** ([`Observer::helio_position`](crate::observers::Observer::helio_position)).
//! - Helpers to convert geodetic latitude/elevation to normalized parallax coordinates
//!   ([`geodetic_to_parallax`](crate::observers::geodetic_to_parallax)) and to lift optional floats into NaN-safe values
//!   ([`to_opt_notnan`](crate::observers::to_opt_notnan)).
//! - A convenience function to compute **three observers’ heliocentric positions** at three epochs
//!   in one call ([`helio_obs_pos`](crate::observers::helio_obs_pos)) — useful for Gauss/Vaisala IOD triplets.
//!
//! ## Frames & conventions
//!
//! - **Earth-fixed (body-fixed)**: the frame in which site coordinates and rotation rate are defined.
//! - **Ecliptic mean J2000**: default geocentric output of [`Observer::pvobs`](crate::observers::Observer::pvobs) (ICRS ecliptic plane).
//! - **Equatorial mean J2000**: default heliocentric output of [`Observer::helio_position`](crate::observers::Observer::helio_position) and
//!   [`helio_obs_pos`](crate::observers::helio_obs_pos) (ICRS-aligned).
//!
//! ```text
//! Body-fixed  --(Earth rotation)-->  Earth-equatorial   --(precession+nutation)-->  Ecliptic J2000
//!                                                                                 \-> Equatorial J2000
//! ```
//!
//! Internally, time-dependent Earth orientation uses **GMST** and the **equation of the equinoxes**,
//! and frame changes are performed via [`rotpn`](crate::ref_system::rotpn) / [`rotmt`](crate::ref_system::rotmt) utilities.
//!
//! ## Units
//!
//! - Longitudes: **degrees** (east positive).
//! - Geocentric parallax (ρ·cosφ, ρ·sinφ): **Earth radii** (dimensionless scaling of geocentric distance).
//! - Positions: **AU**.
//! - Velocities: **AU/day** (from `ω × r`, with `ω = (0, 0, 2π·1.00273790934)` rad/day).
//! - `ra_accuracy`, `dec_accuracy`: **radians**.
//!
//! ## Data flow (typical IOD usage)
//!
//! 1. Build an [`Observer`](crate::observers::Observer) from geodetic inputs (`longitude`, `latitude`, `elevation`) via
//!    [`Observer::new`](crate::observers::Observer::new) (internally calls [`geodetic_to_parallax`](crate::observers::geodetic_to_parallax)) **or** from known (ρ·cosφ, ρ·sinφ)
//!    via [`Observer::from_parallax`](crate::observers::Observer::from_parallax).
//! 2. Compute geocentric PV in **ecliptic J2000** with [`Observer::pvobs`](crate::observers::Observer::pvobs) (needs UT1 provider).
//! 3. Obtain Earth heliocentric state from JPL ephemerides and sum to get the observer’s
//!    **heliocentric equatorial** position with [`Observer::helio_position`](crate::observers::Observer::helio_position).
//! 4. For triplets, call [`helio_obs_pos`](crate::observers::helio_obs_pos) to get the 3×3 matrix of heliocentric positions.
//!
//! ## Quick start
//!
//! ```rust,no_run
//! use hifitime::{Epoch, TimeScale};
//! use nalgebra::{Vector3, Matrix3};
//! use outfit::outfit::Outfit;
//! use outfit::error_models::ErrorModel;
//! use outfit::observers::{Observer, helio_obs_pos};
//!
//! // 1) Environment (JPL ephem + UT1) and site
//! let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//! let site = Observer::new(203.74409, 20.707233557, 3067.694, Some("Pan-STARRS 1".into()), None, None)?;
//!
////! // 2) Geocentric PV (ecliptic J2000)
//! let t = Epoch::from_mjd_in_time_scale(57028.479297592596, TimeScale::TT);
//! let (_x_ecl, _v_ecl) = site.pvobs(&t, state.get_ut1_provider())?;
//!
//! // 3) Heliocentric position (equatorial J2000)
//! let r_helio_eq = site.helio_position(&state, &t, &Vector3::identity())?;
//!
//! // 4) Batch (3 observers, 3 epochs)
//! let tmjd = Vector3::new(57028.479297592596, 57049.24514759259, 57063.97711759259);
//! let R: Matrix3<f64> = helio_obs_pos([&site, &site, &site], &tmjd, &state)?;
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## Design & invariants
//!
//! - [`Observer`](crate::observers::Observer) stores **precomputed body-fixed** position and velocity to avoid recomputing
//!   `ω × r` and trigonometric terms at every call. This is beneficial in tight IOD loops.
//! - `NotNan<f64>` is used for fields where **NaN must be forbidden** (e.g., site geometry).
//!   Use [`to_opt_notnan`](crate::observers::to_opt_notnan) for optional measurement accuracies.
//! - The geodetic-to-parallax conversion accounts for Earth oblateness via
//!   [`EARTH_MAJOR_AXIS`](crate::constants::EARTH_MAJOR_AXIS) / [`EARTH_MINOR_AXIS`](crate::constants::EARTH_MINOR_AXIS).
//!
//! ## Errors
//!
//! - Constructors and helpers return [`OutfitError`](crate::outfit_errors::OutfitError) when NaNs are encountered or a frame
//!   conversion fails; [`to_opt_notnan`](crate::observers::to_opt_notnan) returns `ordered_float::FloatIsNan` if given `Some(NaN)`.
//!
//! ## Testing
//!
//! The module includes unit tests for site construction, body-fixed coordinates,
//! geocentric PV against known values, and multi-epoch heliocentric positions (behind the
//! `jpl-download` feature).
//!
//! ## See also
//! ------------
//! * [`Observer`](crate::observers::Observer) – Site container with precomputed body-fixed state.
//! * [`Observer::pvobs`](crate::observers::Observer::pvobs) – Geocentric PV in **ecliptic J2000**.
//! * [`Observer::helio_position`](crate::observers::Observer::helio_position) – Heliocentric **equatorial J2000** position.
//! * [`helio_obs_pos`](crate::observers::helio_obs_pos) – Batch heliocentric positions for triplets.
//! * [`geodetic_to_parallax`](crate::observers::geodetic_to_parallax) – Geodetic latitude/elevation → (ρ·cosφ, ρ·sinφ).
//! * [`rotpn`](crate::ref_system::rotpn), [`rotmt`](crate::ref_system::rotmt) – Reference-frame transformations.
//! * [`gmst`](crate::time::gmst), [`equequ`](crate::earth_orientation::equequ) – Earth orientation (sidereal time & equation of equinoxes).
//! * [`Outfit`](crate::outfit::Outfit) – Access to JPL ephemerides and UT1 provider.

pub mod bimap;
pub mod observatories;

use hifitime::ut1::Ut1Provider;
use hifitime::Epoch;
use nalgebra::{Matrix3, Vector3};
use ordered_float::NotNan;

use crate::constants::{Degree, Meter, EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS, MJD};
use crate::constants::{DPI, ERAU};
use crate::earth_orientation::equequ;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use crate::ref_system::{rotmt, rotpn, RefEpoch, RefSystem};
use crate::time::gmst;

/// Convert an `Option<f64>` into an `Option<NotNan<f64>>`, propagating `NaN` as an error.
///
/// This helper lifts a possibly missing floating-point value into a `NotNan` container
/// while keeping the outer `Option`. If the inner value is `Some(x)` and `x.is_nan()`,
/// the function returns `Err(FloatIsNan)`. If it is `None`, the result is `Ok(None)`.
///
/// Arguments
/// -----------------
/// * `x`: The optional floating-point value to wrap.
///
/// Return
/// ----------
/// * A `Result` containing `Some(NotNan<f64>)` when `x` is finite, `Ok(None)` when `x` is `None`,
///   or an error if `x` is `NaN`.
///
/// Errors
/// ----------
/// * `ordered_float::FloatIsNan` if `x` is `Some(NaN)`.
///
/// See also
/// ------------
/// * [`ordered_float::NotNan`] – NaN-forbidding wrapper used across the crate.
pub fn to_opt_notnan(x: Option<f64>) -> Result<Option<NotNan<f64>>, ordered_float::FloatIsNan> {
    x.map(NotNan::new).transpose()
}

/// Observer geocentric parameters and precomputed body-fixed state.
///
/// This struct stores:
/// - The observer's **geocentric parallax coordinates** (ρ·cosφ, ρ·sinφ), where ρ is the
///   geocentric distance in **Earth radii** and φ is the **geocentric** latitude.
/// - The **geodetic longitude** (degrees, east of Greenwich).
/// - Optional **astrometric accuracies** for right ascension and declination (radians).
/// - Precomputed **body-fixed** position and velocity vectors used in orbit determination.
///
/// Units
/// -----
/// * `longitude`: degrees (east positive).
/// * `rho_cos_phi`, `rho_sin_phi`: Earth radii (dimensionless scale factor ρ times trig of φ).
/// * `observer_fixed_coord`: astronomical units (AU).
/// * `observer_velocity`: AU/day (from Earth rotation cross product).
/// * `ra_accuracy`, `dec_accuracy`: radians.
///
/// Notes
/// -----
/// The precomputed body-fixed vectors assume a constant Earth rotation rate
/// ω = (0, 0, 2π·1.00273790934) rad/day. Position is scaled by `ERAU` (Earth radius in AU),
/// hence the resulting velocity is in AU/day.
///
/// See also
/// ------------
/// * [`geodetic_to_parallax`] – Converts geodetic latitude/elevation to (ρ·cosφ, ρ·sinφ).
/// * [`Observer::new`] – Construct from geodetic longitude/latitude/elevation.
/// * [`Observer::from_parallax`] – Construct directly from (ρ·cosφ, ρ·sinφ).
/// * [`crate::ref_system::rotpn`] – Reference frame rotations used elsewhere in the pipeline.
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Observer {
    /// Geodetic longitude in **degrees** east of Greenwich.
    pub longitude: NotNan<f64>,

    /// ρ·cosφ (geocentric latitude φ), in **Earth radii** (dimensionless scale).
    pub rho_cos_phi: NotNan<f64>,

    /// ρ·sinφ (geocentric latitude φ), in **Earth radii** (dimensionless scale).
    pub rho_sin_phi: NotNan<f64>,

    /// Optional human-readable site name.
    pub name: Option<String>,

    /// Right ascension measurement accuracy, in **radians** (optional).
    pub ra_accuracy: Option<NotNan<f64>>,

    /// Declination measurement accuracy, in **radians** (optional).
    pub dec_accuracy: Option<NotNan<f64>>,

    /// Precomputed **body-fixed** position of the observer in **AU**.
    observer_fixed_coord: Vector3<NotNan<f64>>,

    /// Precomputed **body-fixed** velocity of the observer in **AU/day**.
    observer_velocity: Vector3<NotNan<f64>>,
}

impl Observer {
    /// Create a new observer from geodetic coordinates.
    ///
    /// This constructor converts `(latitude, elevation)` into geocentric parallax
    /// coordinates `(ρ·cosφ, ρ·sinφ)` using [`geodetic_to_parallax`], builds the
    /// body-fixed position vector in **AU** (scaled by `ERAU`), and computes the
    /// body-fixed velocity as `ω × r` with `ω = (0, 0, 2π·1.00273790934)` in rad/day,
    /// yielding **AU/day**.
    ///
    /// Arguments
    /// -----------------
    /// * `longitude`: Geodetic longitude in **degrees** (east positive).
    /// * `latitude`: Geodetic latitude in **degrees**.
    /// * `elevation`: Height above the reference ellipsoid in **meters**.
    /// * `name`: Optional site name.
    /// * `ra_accuracy`: Optional RA accuracy in **radians**.
    /// * `dec_accuracy`: Optional DEC accuracy in **radians**.
    ///
    /// Return
    /// ----------
    /// * A constructed [`Observer`] with precomputed body-fixed state.
    ///
    /// Errors
    /// ----------
    /// * [`OutfitError`] if the inputs cannot be represented as `NotNan` (e.g., NaN encountered).
    ///
    /// See also
    /// ------------
    /// * [`geodetic_to_parallax`] – Geodetic-to-geocentric parallax conversion.
    /// * [`Observer::from_parallax`] – Build directly from (ρ·cosφ, ρ·sinφ).
    /// * [`crate::ref_system::rotpn`] – Frame rotation utilities used later in the pipeline.
    pub fn new(
        longitude: Degree,
        latitude: Degree,
        elevation: Meter,
        name: Option<String>,
        ra_accuracy: Option<f64>,
        dec_accuracy: Option<f64>,
    ) -> Result<Observer, OutfitError> {
        let (rho_cos_phi, rho_sin_phi) = geodetic_to_parallax(latitude, elevation);

        // Angular velocity of Earth rotation (rad/day) on the z-axis.
        let omega: Vector3<NotNan<f64>> = Vector3::new(
            NotNan::new(0.0)?,
            NotNan::new(0.0)?,
            NotNan::new(DPI * 1.00273790934)?,
        );

        // Body-fixed position in AU from (ρ·cosφ, ρ·sinφ) scaled by Earth radius (AU).
        let lon_radians = longitude.to_radians();
        let body_fixed_coord: Vector3<NotNan<f64>> = Vector3::new(
            NotNan::new(ERAU * rho_cos_phi * lon_radians.cos())?,
            NotNan::new(ERAU * rho_cos_phi * lon_radians.sin())?,
            NotNan::new(ERAU * rho_sin_phi)?,
        );

        // Body-fixed velocity from Earth rotation.
        let dvbf = omega.cross(&body_fixed_coord);

        Ok(Observer {
            longitude: NotNan::try_from(longitude)?,
            rho_cos_phi: NotNan::try_from(rho_cos_phi)?,
            rho_sin_phi: NotNan::try_from(rho_sin_phi)?,
            name,
            ra_accuracy: to_opt_notnan(ra_accuracy)?,
            dec_accuracy: to_opt_notnan(dec_accuracy)?,
            observer_fixed_coord: body_fixed_coord,
            observer_velocity: dvbf,
        })
    }

    /// Create a new observer from geocentric parallax coordinates.
    ///
    /// This constructor skips the geodetic-to-parallax conversion and directly uses
    /// `(ρ·cosφ, ρ·sinφ)` to build the body-fixed position in **AU** (scaled by `ERAU`),
    /// and the body-fixed velocity in **AU/day** as `ω × r` with
    /// `ω = (0, 0, 2π·1.00273790934)` rad/day.
    ///
    /// Arguments
    /// -----------------
    /// * `longitude`: Geodetic longitude in **degrees** (east positive).
    /// * `rho_cos_phi`: ρ·cosφ in **Earth radii** (dimensionless).
    /// * `rho_sin_phi`: ρ·sinφ in **Earth radii** (dimensionless).
    /// * `name`: Optional site name.
    /// * `ra_accuracy`: Optional RA accuracy in **radians**.
    /// * `dec_accuracy`: Optional DEC accuracy in **radians**.
    ///
    /// Return
    /// ----------
    /// * A constructed [`Observer`] with precomputed body-fixed state.
    ///
    /// Errors
    /// ----------
    /// * [`OutfitError`] if inputs cannot be represented as `NotNan` (e.g., NaN encountered).
    ///
    /// See also
    /// ------------
    /// * [`Observer::new`] – Build from geodetic latitude and elevation.
    /// * [`geodetic_to_parallax`] – For the forward conversion when geodetic inputs are available.
    pub fn from_parallax(
        longitude: Degree,
        rho_cos_phi: f64,
        rho_sin_phi: f64,
        name: Option<String>,
        ra_accuracy: Option<f64>,
        dec_accuracy: Option<f64>,
    ) -> Result<Observer, OutfitError> {
        // Angular velocity of Earth rotation (rad/day) on the z-axis.
        let omega: Vector3<NotNan<f64>> = Vector3::new(
            NotNan::new(0.0)?,
            NotNan::new(0.0)?,
            NotNan::new(DPI * 1.00273790934)?,
        );

        // Body-fixed position in AU from (ρ·cosφ, ρ·sinφ) scaled by Earth radius (AU).
        let lon_radians = longitude.to_radians();
        let body_fixed_coord: Vector3<NotNan<f64>> = Vector3::new(
            NotNan::new(ERAU * rho_cos_phi * lon_radians.cos())?,
            NotNan::new(ERAU * rho_cos_phi * lon_radians.sin())?,
            NotNan::new(ERAU * rho_sin_phi)?,
        );

        // Body-fixed velocity from Earth rotation.
        let dvbf = omega.cross(&body_fixed_coord);

        Ok(Observer {
            longitude: NotNan::try_from(longitude)?,
            rho_cos_phi: NotNan::try_from(rho_cos_phi)?,
            rho_sin_phi: NotNan::try_from(rho_sin_phi)?,
            name,
            ra_accuracy: to_opt_notnan(ra_accuracy)?,
            dec_accuracy: to_opt_notnan(dec_accuracy)?,
            observer_fixed_coord: body_fixed_coord,
            observer_velocity: dvbf,
        })
    }

    /// Get the fixed position of an observatory using its geographic coordinates
    ///
    /// Argument
    /// --------
    /// * longitude: observer longitude in Degree
    /// * latitude: observer latitude in Degree
    /// * height: observer height in Degree
    ///
    /// Return
    /// ------
    /// * observer fixed coordinates vector on the Earth (not corrected from Earth motion)
    /// * units is AU
    pub fn body_fixed_coord(&self) -> Vector3<f64> {
        let lon_radians = self.longitude.to_radians();

        Vector3::new(
            ERAU * self.rho_cos_phi.into_inner() * lon_radians.cos(),
            ERAU * self.rho_cos_phi.into_inner() * lon_radians.sin(),
            ERAU * self.rho_sin_phi.into_inner(),
        )
    }

    /// Compute the observer’s geocentric position and velocity in the ecliptic J2000 frame.
    ///
    /// This function calculates the position and velocity of a ground-based observer relative to the Earth's
    /// center of mass, accounting for Earth rotation (via GMST), nutation, and the observer’s geographic location.
    /// The result is expressed in the ecliptic mean J2000 frame, suitable for use in orbital initial determination.
    ///
    /// Arguments
    /// ---------
    /// * `observer`: a reference to an [`Observer`] containing the site longitude and parallax parameters.
    /// * `tmjd`: observation epoch as a [`hifitime::Epoch`] in TT.
    /// * `ut1_provider`: a reference to a [`hifitime::ut1::Ut1Provider`] for accurate UT1 conversion.
    ///
    /// Returns
    /// --------
    /// * `(dx, dv)` – Tuple of:
    ///     - `dx`: observer geocentric position vector in ecliptic mean J2000 frame \[AU\].
    ///     - `dv`: observer velocity vector due to Earth's rotation, in the same frame \[AU/day\].
    ///
    /// Remarks
    /// -------
    /// * Internally, this function:
    ///     1. get the body-fixed coordinates of the observer.
    ///     2. get its rotational velocity: `v = ω × r`.
    ///     3. Applies Earth orientation corrections using:
    ///         - Greenwich Mean Sidereal Time (GMST),
    ///         - Equation of the equinoxes,
    ///         - Precession and nutation transformation (`rotpn`).
    ///     4. Returns position and velocity in the J2000 ecliptic frame (used in classical orbital mechanics).
    ///
    /// # See also
    /// * [`Observer::body_fixed_coord`] – observer's base vector in Earth-fixed frame
    /// * [`rotpn`] – rotation between reference frames
    /// * [`gmst`], [`equequ`] – time-dependent Earth orientation
    pub fn pvobs(
        &self,
        tmjd: &Epoch,
        ut1_provider: &Ut1Provider,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        // Get observer position and velocity in the Earth-fixed frame
        let dxbf = self.observer_fixed_coord.map(|x| x.into_inner());
        let dvbf = self.observer_velocity.map(|x| x.into_inner());

        // deviation from Orbfit, use of another conversion from MJD UTC (ET scale) to UT1 scale
        // based on the hifitime crate
        let mjd_ut1 = tmjd.to_ut1(ut1_provider);
        let tut = mjd_ut1.to_mjd_tai_days();

        // Compute the Greenwich sideral apparent time
        let gast = gmst(tut) + equequ(tmjd.to_mjd_tt_days());

        // Earth rotation matrix
        let rot = rotmt(-gast, 2);

        // Compute the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let rer_sys1 = RefSystem::Equt(RefEpoch::Epoch(tmjd.to_mjd_tt_days()));
        let rer_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let rot1 = rotpn(&rer_sys1, &rer_sys2)?;

        let rot1_mat = rot1.transpose();
        let rot_mat = rot.transpose();

        let rotmat = rot1_mat * rot_mat;

        // Apply transformation to the observer position and velocity
        let dx = rotmat * dxbf;
        let dv = rotmat * dvbf;

        Ok((dx, dv))
    }

    /// Compute the observer’s heliocentric position in the **equatorial mean J2000** frame.
    ///
    /// This method forms the full heliocentric position of the observing site by combining:
    /// - the site **geocentric** position vector at `epoch`, and
    /// - the Earth’s **heliocentric** position from the JPL ephemerides.
    ///
    /// The input geocentric vector is assumed to be expressed in the **ecliptic mean J2000** frame
    /// (AU). It is rotated to **equatorial mean J2000**, then added to Earth’s heliocentric
    /// position (also in equatorial mean J2000).
    ///
    /// Arguments
    /// -----------------
    /// * `state` – [`Outfit`] environment providing JPL ephemerides and frame rotations.
    /// * `epoch` – Observation epoch in the **TT** time scale.
    /// * `observer_geocentric_position` – Geocentric site position **in ecliptic mean J2000** (AU).
    ///
    /// Return
    /// ----------
    /// * `Result<Vector3<f64>, OutfitError>` – Observer’s **heliocentric** position at `epoch`,
    ///   in **AU**, expressed in **equatorial mean J2000**.
    ///
    /// Remarks
    /// -------------
    /// * If your geocentric site vector is already in **equatorial** J2000, rotate it to
    ///   **ecliptic** before calling this method, or adapt the rotation accordingly.
    /// * This routine is typically used internally when constructing per-observation geometry
    ///   (e.g., within `Observation::new`), ensuring consistent frames for Gauss IOD.
    ///
    /// See also
    /// ------------
    /// * [`Observer::pvobs`] – Geocentric position (and velocity) of the site at `epoch`.
    /// * [`Outfit::get_jpl_ephem`] – Access Earth’s heliocentric state from JPL ephemerides.
    /// * [`Outfit::get_rot_eclmj2000_to_equmj2000`] – Rotation between ecliptic and equatorial J2000.
    pub fn helio_position(
        &self,
        state: &Outfit,
        epoch: &Epoch,
        observer_geocentric_position: &Vector3<f64>,
    ) -> Result<Vector3<f64>, OutfitError> {
        let jpl = state.get_jpl_ephem().unwrap();

        // Earth's heliocentric position
        let earth_pos = jpl.earth_ephemeris(epoch, false).0;

        // Transform observer position from ecliptic to equatorial J2000
        let rot_matrix = state.get_rot_eclmj2000_to_equmj2000().transpose();

        Ok(earth_pos + rot_matrix * observer_geocentric_position)
    }

    /// Recover geodetic latitude and ellipsoidal height (WGS-84) from parallax constants.
    ///
    /// Inverts the stored parallax coordinates `(ρ·cosφ, ρ·sinφ)` to the **geodetic**
    /// latitude `φ` (degrees) and the ellipsoidal height `h` (meters) above the
    /// WGS-84 reference ellipsoid. The computation uses **Bowring’s closed-form**
    /// formula (no iteration), which is usually sufficient for double-precision
    /// accuracy at the centimeter level or better.
    ///
    /// Units & model
    /// -------------
    /// * Inputs: `ρ·cosφ` and `ρ·sinφ` are dimensionless, expressed in **Earth radii**
    ///   (normalized by the equatorial radius). They are scaled internally by `a`
    ///   (the equatorial radius) to meters.
    /// * Output: latitude in **degrees**, height in **meters** (ellipsoidal height, not geoid/orthometric).
    /// * Ellipsoid: WGS-84 radii from `constants.rs` (`EARTH_MAJOR_AXIS` = `a`, `EARTH_MINOR_AXIS` = `b`).
    ///   If you prefer exact GRS-80 reproduction, use consistent `b` there; the difference vs WGS-84 is sub-millimetric.
    ///
    /// Notes
    /// -----
    /// * This routine **does not** compute the geodetic longitude; it only returns `(lat, h)`.
    ///   Your `Observer` already stores the geodetic longitude independently.
    /// * Numerical robustness is good across latitudes, including near the poles.
    /// * If you require bit-for-bit parity with an external reference using a different ellipsoid,
    ///   ensure `a`/`b` match that reference.
    ///
    /// Arguments
    /// -----------------
    /// * None.
    ///
    /// Return
    /// ----------
    /// * `(geodetic_latitude_deg, height_meters)` — latitude in degrees, ellipsoidal height in meters.
    ///
    /// See also
    /// ------------
    /// * [`geodetic_to_parallax`] – Forward conversion used at construction.
    /// * [`Observer::from_parallax`] – Builds an observer from `(ρ·cosφ, ρ·sinφ)`.
    /// * `constants::EARTH_MAJOR_AXIS` / `constants::EARTH_MINOR_AXIS` – Ellipsoid radii used here.
    pub fn geodetic_lat_height_wgs84(&self) -> (f64, f64) {
        let a = EARTH_MAJOR_AXIS;
        let b = EARTH_MINOR_AXIS;
        let e2 = 1.0 - (b * b) / (a * a);
        let ep2 = (a * a) / (b * b) - 1.0;

        let p = self.rho_cos_phi.into_inner() * a; // distance in equatorial plane [m]
        let z = self.rho_sin_phi.into_inner() * a; // z [m]

        // Bowring’s formula:
        let theta = (z * a).atan2(p * b);
        let st = theta.sin();
        let ct = theta.cos();
        let phi = (z + ep2 * b * st.powi(3)).atan2(p - e2 * a * ct.powi(3));

        let s = phi.sin();
        let n = a / (1.0 - e2 * s * s).sqrt();
        let h = p / phi.cos() - n;

        (phi.to_degrees(), h)
    }
}

/// Convert geodetic latitude and height into normalized parallax coordinates
/// on the Earth.
///
/// This transformation is used to compute the observer's position on Earth
/// in a way that accounts for the Earth's oblateness. The resulting values
/// are dimensionless and are expressed in units of the Earth's equatorial
/// radius (`EARTH_MAJOR_AXIS`).
///
/// Arguments
/// ---------
/// * `lat` - Geodetic latitude of the observer in **radians**.
/// * `height` - Observer's altitude above the reference ellipsoid in **meters**.
///
/// Returns
/// -------
/// A tuple `(rho_cos_phi, rho_sin_phi)`:
/// * `rho_cos_phi`: normalized distance of the observer projected on
///   the Earth's equatorial plane.
/// * `rho_sin_phi`: normalized distance of the observer projected on
///   the Earth's rotation (polar) axis.
///
/// Details
/// -------
/// The computation uses the reference ellipsoid defined by:
/// * `EARTH_MAJOR_AXIS`: Equatorial radius (m),
/// * `EARTH_MINOR_AXIS`: Polar radius (m).
///
/// The formula comes from standard geodetic to geocentric conversion:
///
/// ```text
/// u = atan( (sin φ * (b/a)) / cos φ )
/// ρ_sinφ = (b/a) * sin u + (h/a) * sin φ
/// ρ_cosφ = cos u + (h/a) * cos φ
/// ```
///
/// where `a` and `b` are the Earth's semi-major and semi-minor axes,
/// and `h` is the height above the ellipsoid.
///
/// See also
/// --------
/// * [`Observer::body_fixed_coord`] – Uses this function to compute
///   the observer's fixed position in Earth-centered coordinates.
pub fn lat_alt_to_parallax(lat: f64, height: f64) -> (f64, f64) {
    // Ratio of the Earth's minor to major axis (flattening factor)
    let axis_ratio = EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS;

    // Compute the auxiliary angle u (parametric latitude)
    // This corrects for the Earth's oblateness.
    let u = (lat.sin() * axis_ratio).atan2(lat.cos());

    // Compute the normalized distance along the polar axis
    let rho_sin_phi = axis_ratio * u.sin() + (height / EARTH_MAJOR_AXIS) * lat.sin();

    // Compute the normalized distance along the equatorial plane
    let rho_cos_phi = u.cos() + (height / EARTH_MAJOR_AXIS) * lat.cos();

    (rho_cos_phi, rho_sin_phi)
}

/// Convert geodetic latitude (in degrees) and height (in meters)
/// into normalized parallax coordinates.
///
/// This is a convenience wrapper around [`lat_alt_to_parallax`] that
/// performs the degrees-to-radians conversion before calling the main
/// function.
///
/// Arguments
/// ---------
/// * `lat` - Geodetic latitude of the observer in **degrees**.
/// * `height` - Observer's altitude above the reference ellipsoid in **meters**.
///
/// Returns
/// -------
/// A tuple `(rho_cos_phi, rho_sin_phi)`:
/// * `rho_cos_phi`: normalized distance of the observer projected on
///   the Earth's equatorial plane.
/// * `rho_sin_phi`: normalized distance of the observer projected on
///   the Earth's rotation (polar) axis.
///
/// Details
/// -------
/// This function simply converts `lat` to radians and delegates the
/// computation to [`lat_alt_to_parallax`].
///
/// See also
/// --------
/// * [`lat_alt_to_parallax`] – Performs the actual computation given latitude in radians.
pub fn geodetic_to_parallax(lat: f64, height: f64) -> (f64, f64) {
    // Convert latitude from degrees to radians
    let latitude_rad = lat.to_radians();

    // Call the main routine that works with radians
    let (rho_cos_phi, rho_sin_phi) = lat_alt_to_parallax(latitude_rad, height);

    (rho_cos_phi, rho_sin_phi)
}

/// Compute the heliocentric positions of three observers at their respective epochs,
/// expressed in the **equatorial mean J2000** frame (ICRS-aligned).
///
/// Overview
/// -----------------
/// This routine builds a `3×3` matrix of observer positions by combining:
/// - the **geocentric site position** of each observer (from [`Observer::pvobs`]),
/// - the Earth’s **heliocentric barycentric position** from JPL ephemerides,
/// - a frame transformation from **ecliptic mean J2000** (site positions) to
///   **equatorial mean J2000** (final output).
///
/// The result is a compact representation where each column corresponds to one
/// observer/epoch pair:  
/// `observers[0] ↔ mjd_tt.x`,  
/// `observers[1] ↔ mjd_tt.y`,  
/// `observers[2] ↔ mjd_tt.z`.
///
/// Arguments
/// -----------------
/// * `observers` – Array of three [`Observer`] references, each encoding the site geometry
///   (longitude, normalized geocentric radius components, etc.).
/// * `mjd_tt` – [`Vector3<MJD>`] of observation epochs in Terrestrial Time (TT), one per observer.
/// * `state` – [`Outfit`] environment providing:
///   - JPL planetary ephemerides (via [`Outfit::get_jpl_ephem`]),
///   - UT1 provider for Earth rotation/orientation (via [`Outfit::get_ut1_provider`]).
///
/// Return
/// ----------
/// * `Result<Matrix3<f64>, OutfitError>` – A 3×3 matrix of observer heliocentric positions, with:  
///   - **Columns** = `[r₁, r₂, r₃]`, one per observer/epoch,  
///   - **Units** = astronomical units (AU),  
///   - **Frame** = equatorial mean J2000 (ICRS-aligned).
///
/// Remarks
/// -------------
/// * For each observer/time pair:
///   1. The site’s **geocentric position** is computed via [`Observer::pvobs`] (AU, ecliptic J2000).
///   2. Earth’s heliocentric position is retrieved from the JPL ephemeris.
///   3. The site position is rotated into **equatorial mean J2000** using the frame rotation.
///   4. The Earth + rotated site vectors give the full heliocentric observer position.
/// * This function is mainly used during **Gauss IOD** preparation to populate the
///   observer position matrix stored in [`GaussObs`](crate::initial_orbit_determination::gauss::GaussObs).
///
/// See also
/// ------------
/// * [`Observer::pvobs`] – Geocentric observer position at a given epoch (ecliptic J2000).
/// * [`Observer::helio_position`] – Per-observer heliocentric position (equatorial J2000).
/// * [`Outfit::get_jpl_ephem`] – Access to planetary ephemerides (Earth state).
/// * [`Outfit::get_ut1_provider`] – Provides Earth orientation parameters (ΔUT1).
pub fn helio_obs_pos(
    observers: [&Observer; 3],
    mjd_tt: &Vector3<MJD>,
    state: &Outfit,
) -> Result<Matrix3<f64>, OutfitError> {
    let epochs = [
        Epoch::from_mjd_in_time_scale(mjd_tt.x, hifitime::TimeScale::TT),
        Epoch::from_mjd_in_time_scale(mjd_tt.y, hifitime::TimeScale::TT),
        Epoch::from_mjd_in_time_scale(mjd_tt.z, hifitime::TimeScale::TT),
    ];

    let pvobs1 = observers[0].pvobs(&epochs[0], state.get_ut1_provider())?;
    let pvobs2 = observers[1].pvobs(&epochs[1], state.get_ut1_provider())?;
    let pvobs3 = observers[2].pvobs(&epochs[2], state.get_ut1_provider())?;

    let positions = [
        observers[0].helio_position(state, &epochs[0], &pvobs1.0)?,
        observers[1].helio_position(state, &epochs[1], &pvobs2.0)?,
        observers[2].helio_position(state, &epochs[2], &pvobs3.0)?,
    ];

    Ok(Matrix3::from_columns(&positions))
}

#[cfg(test)]
mod observer_test {

    use crate::{error_models::ErrorModel, outfit::Outfit};

    use super::*;

    #[test]
    fn test_observer_constructor() {
        let observer = Observer::new(0.0, 0.0, 0.0, None, None, None).unwrap();
        assert_eq!(observer.longitude, 0.0);
        assert_eq!(observer.rho_cos_phi, 1.0);
        assert_eq!(observer.rho_sin_phi, 0.0);

        let observer = Observer::new(
            289.25058,
            -30.2446,
            2647.,
            Some("Rubin Observatory".to_string()),
            Some(0.0001),
            Some(0.0001),
        )
        .unwrap();

        assert_eq!(observer.longitude, 289.25058);
        assert_eq!(observer.rho_cos_phi, 0.8649760504617418);
        assert_eq!(observer.rho_sin_phi, -0.5009551027512434);
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, None, None, None).unwrap();
        assert_eq!(
            pan_starrs.body_fixed_coord(),
            Vector3::new(
                -0.00003653799439776371,
                -0.00001607260397528885,
                0.000014988110430544328
            )
        );

        assert_eq!(
            pan_starrs.observer_fixed_coord,
            Vector3::new(
                NotNan::new(-0.00003653799439776371).unwrap(),
                NotNan::new(-0.00001607260397528885).unwrap(),
                NotNan::new(0.000014988110430544328).unwrap()
            )
        )
    }

    #[test]
    fn pvobs_test() {
        let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let tmjd = 57028.479297592596;
        let epoch = Epoch::from_mjd_in_time_scale(tmjd, hifitime::TimeScale::TT);
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs =
            Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None).unwrap();

        let (observer_position, observer_velocity) =
            &pan_starrs.pvobs(&epoch, state.get_ut1_provider()).unwrap();

        assert_eq!(
            observer_position.as_slice(),
            [
                -2.086211182493635e-5,
                3.718476815327979e-5,
                2.4978996447997476e-7
            ]
        );
        assert_eq!(
            observer_velocity.as_slice(),
            [
                -0.0002143246535691577,
                -0.00012059801691431748,
                5.262184624215718e-5
            ]
        );
    }

    #[test]
    fn geodetic_to_parallax_test() {
        // latitude and height of Pan-STARRS 1, Haleakala
        let (pxy1, pz1) = geodetic_to_parallax(20.707233557, 3067.694);
        assert_eq!(pxy1, 0.9362410003211518);
        assert_eq!(pz1, 0.35154299856304305);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_helio_pos_obs() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let tmjd = Vector3::new(
            57028.479297592596,
            57_049.245_147_592_59,
            57_063.977_117_592_59,
        );

        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs =
            Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None).unwrap();

        // Now we need a Vector3<Observer> with three identical copies
        let observers = [&pan_starrs, &pan_starrs, &pan_starrs];

        let helio_pos = helio_obs_pos(observers, &tmjd, &OUTFIT_HORIZON_TEST.0).unwrap();

        assert_eq!(
            helio_pos.as_slice(),
            [
                -0.2645666171464416,
                0.8689351643701766,
                0.3766996211107864,
                -0.5891631852137064,
                0.7238872516824697,
                0.3138186516540669,
                -0.7743280306286537,
                0.5612532665812755,
                0.24333415479994636
            ]
        );
    }

    #[cfg(test)]
    mod geodetic_inverse_tests {
        use super::*;
        use crate::constants::Degree;
        use approx::assert_abs_diff_eq;

        /// Round-trip a single site through (lon, lat, h) -> parallax -> inverse
        /// and check that we recover the original geodetic latitude & height.
        ///
        /// Notes
        /// -----
        /// * `Observer::new` is given `h_m` in meters (as per current API usage).
        /// * `geodetic_lat_height_wgs84()` returns height in **meters**; we convert to meters.
        fn roundtrip_site(name: &str, lon_deg: Degree, lat_deg: Degree, h_m: f64) {
            // Build observer (forward: geodetic -> parallax is done inside `Observer::new`)
            let obs = Observer::new(lon_deg, lat_deg, h_m, Some(name.to_string()), None, None)
                .expect("Failed to create observer");

            // Inverse: parallax -> geodetic (WGS-84)
            let (lat_rec_deg, h_rec_m) = obs.geodetic_lat_height_wgs84();

            // Tolerances:
            // - Latitude: 1e-6 deg (~3.6 mas) – tight but should pass for double precision Bowring + 0–1 Newton step
            // - Height:   1e-2 m
            let tol_lat_deg = 1e-6;
            let tol_h_m = 1e-2;

            assert_abs_diff_eq!(lat_rec_deg, lat_deg, epsilon = tol_lat_deg);
            assert_abs_diff_eq!(h_rec_m, h_m, epsilon = tol_h_m);
        }

        /// See also
        /// ------------
        /// * [`Observer::new`] – Forward geodetic->parallax conversion under test by round-trip.
        /// * [`Observer::from_parallax`] – Alternative constructor, if you want to inject ρ·cosφ/ρ·sinφ.
        /// * `geodetic_to_parallax` – The forward routine used internally by `Observer::new`.

        #[test]
        fn geodetic_roundtrip_known_observatories_wgs84() {
            // NOTE:
            // The heights below are commonly quoted "above sea level" (orthometric).
            // For pure algorithmic round-trip testing, that's acceptable because we feed
            // the same height into forward and inverse. If you want strict ellipsoidal
            // (h) values, substitute official WGS-84 heights here.
            let sites: &[(&str, Degree, Degree, f64)] = &[
                // name,                      lon_deg (E+), lat_deg (N+), height_m
                ("Haleakala (PS1 I41)", -156.2575, 20.7075, 3055.0),
                ("Mauna Kea (CFHT)", -155.4700, 19.8261, 4205.0),
                ("ESO Paranal", -70.4025, -24.6252, 2635.0),
                ("Cerro Pachon (Rubin)", -70.7366, -30.2407, 2663.0),
                ("La Silla", -70.7346, -29.2613, 2400.0),
                ("Kitt Peak", -111.5967, 31.9583, 2096.0),
                ("Roque de los Muchachos", -17.8947, 28.7606, 2396.0),
            ];

            for (name, lon, lat, h_m) in sites.iter().copied() {
                roundtrip_site(name, lon, lat, h_m);
            }
        }

        #[test]
        fn geodetic_roundtrip_extremes_equator_and_pole() {
            // Equator, sea level
            roundtrip_site("Equator (0°, 0 m)", 0.0, 0.0, 0.0);

            // Near-North-Pole and Near-South-Pole with modest height
            roundtrip_site("Near North Pole", 0.0, 89.999, 1000.0);
            roundtrip_site("Near South Pole", 0.0, -89.999, 1000.0);
        }

        #[test]
        fn geodetic_roundtrip_high_altitude_and_negative() {
            // Very high site (simulate balloon/aircraft)
            roundtrip_site("High Alt 10 m", 10.0, 45.0, 10_000.0);

            // Negative height (below ellipsoid, synthetic but tests robustness)
            roundtrip_site("Below ellipsoid -50 m", -30.0, -10.0, -50.0);
        }
    }
}
