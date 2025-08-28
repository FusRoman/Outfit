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
//! let r_helio_eq = site.helio_position(&t, &state)?;
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

use crate::constants::{Degree, Kilometer, EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS, MJD};
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
    /// * `elevation`: Height above the reference ellipsoid in **kilometers**.
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
        elevation: Kilometer,
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

    /// Compute the heliocentric position of this observer in the **equatorial J2000 frame**.
    ///
    /// This method combines the observer’s **geocentric position** (computed from Earth rotation,
    /// orientation, and the observer’s fixed geographic coordinates) with the Earth's **heliocentric
    /// barycentric position** from the JPL ephemerides, to obtain the observer’s full heliocentric
    /// position vector.
    ///
    /// # Reference frame
    /// The output vector is expressed in the **mean equatorial J2000 frame (ICRS-aligned)**,
    /// with units in astronomical units (AU).
    ///
    /// Arguments
    /// ---------
    /// * `epoch` – Epoch of observation (TT time scale).
    /// * `state` – [`Outfit`] providing:
    ///   * A JPL planetary ephemeris via [`Outfit::get_jpl_ephem`],
    ///   * A UT1 provider via [`Outfit::get_ut1_provider`].
    ///
    /// Returns
    /// -------
    /// * `Vector3<f64>` – The heliocentric position of this observer at the specified epoch,
    ///   in AU and in the equatorial J2000 frame.
    ///
    /// Algorithm
    /// ---------
    /// 1. Compute the observer’s **geocentric position** in the ecliptic J2000 frame using [`Observer::pvobs`].
    /// 2. Query the JPL ephemerides for the Earth's heliocentric position.
    /// 3. Transform the geocentric position from **ecliptic J2000** to **equatorial J2000** using [`rotpn`].
    /// 4. Add the Earth’s heliocentric position and the transformed geocentric vector to obtain the
    ///    observer’s heliocentric position.
    ///
    /// See also
    /// --------
    /// * [`Observer::pvobs`] – Geocentric position of the observer
    /// * [`rotpn`] – Frame transformation between ecliptic and equatorial J2000
    /// * [`Outfit::get_jpl_ephem`] – Access to planetary ephemerides
    pub fn helio_position(
        &self,
        epoch: &Epoch,
        state: &Outfit,
    ) -> Result<Vector3<f64>, OutfitError> {
        let ut1_provider = state.get_ut1_provider();
        let jpl = state.get_jpl_ephem().unwrap();

        // Geocentric position in ecliptic J2000
        let obs_pos_ecl = self.pvobs(epoch, ut1_provider)?.0;

        // Earth's heliocentric position
        let earth_pos = jpl.earth_ephemeris(epoch, false).0;

        // Transform observer position from ecliptic to equatorial J2000
        let rot_matrix = state.get_rot_eclmj2000_to_equmj2000().transpose();

        Ok(earth_pos + rot_matrix * obs_pos_ecl)
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
/// * `height` - Observer's altitude above the reference ellipsoid in **kilometers**.
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
/// * `EARTH_MAJOR_AXIS`: Equatorial radius (km),
/// * `EARTH_MINOR_AXIS`: Polar radius (km).
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

/// Convert geodetic latitude (in degrees) and height (in kilometers)
/// into normalized parallax coordinates.
///
/// This is a convenience wrapper around [`lat_alt_to_parallax`] that
/// performs the degrees-to-radians conversion before calling the main
/// function.
///
/// Arguments
/// ---------
/// * `lat` - Geodetic latitude of the observer in **degrees**.
/// * `height` - Observer's altitude above the reference ellipsoid in **kilometers**.
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

/// Compute the heliocentric equatorial positions of three different observers at their respective observation times.
///
/// This function returns the heliocentric positions of three ground-based observers at three distinct epochs,
/// expressed in the equatorial mean J2000 frame (ICRS-aligned), by combining:
/// - the barycentric position of the Earth (from a JPL planetary ephemeris),
/// - the geocentric position of each observer (accounting for Earth orientation and rotation).
///
/// Each observer is paired with a corresponding time: `observers[0]` ↔ `mjd_tt.x`, `observers[1]` ↔ `mjd_tt.y`, `observers[2]` ↔ `mjd_tt.z`.
///
/// Arguments
/// ---------
/// * `observers`: an array of references `[&Observer; 3]`, each containing the geodetic parameters
///   (longitude, normalized radius components) of one observer.
/// * `mjd_tt`: a [`Vector3<MJD>`] of observation epochs in Terrestrial Time (TT), one per observer.
/// * `state`: an [`Outfit`] providing access to:
///     - the JPL planetary ephemeris (`get_jpl_ephem()`),
///     - the UT1 provider (`get_ut1_provider()`).
///
/// Returns
/// --------
/// * `Matrix3<f64>`: a 3×3 matrix where each column is the heliocentric position vector of the corresponding observer:
///     - Columns: `[r₁, r₂, r₃]` for observers at epochs `mjd_tt.x`, `mjd_tt.y`, `mjd_tt.z`
///     - Units: astronomical units (AU)
///     - Frame: equatorial mean J2000 (ICRS-aligned)
///
/// Remarks
/// -------
/// * This function performs the following steps for each observer/time pair:
///     1. Computes the observer’s geocentric position in the ecliptic mean J2000 frame via [`Observer::pvobs`].
///     2. Retrieves the heliocentric Earth position from the JPL ephemeris.
///     3. Transforms the observer’s position from the ecliptic mean J2000 frame to the equatorial mean J2000 frame using [`rotpn`].
///     4. Computes the heliocentric observer position by summing the Earth and transformed observer vectors.
///
///
/// # See also
/// * [`Observer::pvobs`] – computes the observer’s geocentric position in the ecliptic mean J2000 frame
/// * [`rotpn`] – transforms vectors between celestial reference frames (e.g., Eclm to Equm)
/// * [`Outfit::get_jpl_ephem`] – accesses planetary ephemerides for Earth position
/// * [`Outfit::get_ut1_provider`] – provides Earth orientation parameters (e.g., ΔUT1)
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

    let positions = [
        observers[0].helio_position(&epochs[0], state)?,
        observers[1].helio_position(&epochs[1], state)?,
        observers[2].helio_position(&epochs[2], state)?,
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
}
