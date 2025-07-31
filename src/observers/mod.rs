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

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Observer {
    // in degrees east of Greenwich
    pub longitude: ordered_float::NotNan<f64>,
    // phi' is the geocentric latitude and rho is the geocentric distance in earth radii
    // rho cos phi in degrees
    pub rho_cos_phi: ordered_float::NotNan<f64>,
    // rho sin phi in kilometers
    pub rho_sin_phi: ordered_float::NotNan<f64>,
    pub name: Option<String>,
    // Accuracy of the right ascension and declination
    pub ra_accuracy: Option<ordered_float::NotNan<f64>>,
    pub dec_accuracy: Option<ordered_float::NotNan<f64>>,
}

impl Observer {
    /// Create a new observer
    ///
    /// Argument
    /// --------
    /// * `longitude`: observer longitude
    /// * `latitude`: observer latitude
    /// * `elevation`: observer elevation
    /// * `name`: observer name
    ///
    /// Return
    /// ------
    /// * Observer object
    pub(crate) fn new(
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
        ra_accuracy: Option<ordered_float::NotNan<f64>>,
        dec_accuracy: Option<ordered_float::NotNan<f64>>,
    ) -> Observer {
        let (rho_cos_phi, rho_sin_phi) = geodetic_to_parallax(latitude, elevation);
        Observer {
            longitude: NotNan::try_from(longitude).expect("Longitude cannot be NaN"),
            rho_cos_phi: NotNan::try_from(rho_cos_phi).expect("Longitude cannot be NaN"),
            rho_sin_phi: NotNan::try_from(rho_sin_phi).expect("Longitude cannot be NaN"),
            name,
            ra_accuracy,
            dec_accuracy,
        }
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
    pub(crate) fn body_fixed_coord(&self) -> Vector3<f64> {
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
    ///     1. Computes the body-fixed coordinates of the observer.
    ///     2. Derives its rotational velocity: `v = ω × r`.
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
    pub(crate) fn pvobs(
        &self,
        tmjd: &Epoch,
        ut1_provider: &Ut1Provider,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        // Initialisation
        let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

        // Get the coordinates of the observer on Earth
        let dxbf = self.body_fixed_coord();

        // Get the observer velocity due to Earth rotation
        let dvbf = omega.cross(&dxbf);

        // deviation from Orbfit, use of another conversion from MJD UTC (ET scale) to UT1 scale
        // based on the hifitime crate
        let mjd_ut1 = tmjd.to_ut1(ut1_provider.to_owned());
        let tut = mjd_ut1.to_mjd_tai_days();

        // Compute the Greenwich sideral apparent time
        let gast = gmst(tut) + equequ(tmjd.to_mjd_tt_days());

        // Earth rotation matrix
        let rot = rotmt(-gast, 2);

        // Compute the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let rer_sys1 = RefSystem::Equt(RefEpoch::Epoch(tmjd.to_mjd_tt_days()));
        let rer_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let rot1 = rotpn(&rer_sys1, &rer_sys2)?;

        let rot1_mat = Matrix3::from(rot1).transpose();
        let rot_mat = Matrix3::from(rot).transpose();

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
    /// # Arguments
    /// * `epoch` – Epoch of observation (TT time scale).
    /// * `state` – [`Outfit`] providing:
    ///   * A JPL planetary ephemeris via [`Outfit::get_jpl_ephem`],
    ///   * A UT1 provider via [`Outfit::get_ut1_provider`].
    ///
    /// # Returns
    /// * `Vector3<f64>` – The heliocentric position of this observer at the specified epoch,
    ///   in AU and in the equatorial J2000 frame.
    ///
    /// # Algorithm
    /// 1. Compute the observer’s **geocentric position** in the ecliptic J2000 frame using [`Observer::pvobs`].
    /// 2. Query the JPL ephemerides for the Earth's heliocentric position.
    /// 3. Transform the geocentric position from **ecliptic J2000** to **equatorial J2000** using [`rotpn`].
    /// 4. Add the Earth’s heliocentric position and the transformed geocentric vector to obtain the
    ///    observer’s heliocentric position.
    ///
    /// # See also
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
        let ref_sys1 = RefSystem::Eclm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let rot = rotpn(&ref_sys1, &ref_sys2)?;
        let rot_matrix = Matrix3::from(rot).transpose();

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
/// # Arguments
/// * `lat` - Geodetic latitude of the observer in **radians**.
/// * `height` - Observer's altitude above the reference ellipsoid in **kilometers**.
///
/// # Returns
/// A tuple `(rho_cos_phi, rho_sin_phi)`:
/// * `rho_cos_phi`: normalized distance of the observer projected on
///   the Earth's equatorial plane.
/// * `rho_sin_phi`: normalized distance of the observer projected on
///   the Earth's rotation (polar) axis.
///
/// # Details
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
/// # See also
/// * [`Observer::body_fixed_coord`] – Uses this function to compute
///   the observer's fixed position in Earth-centered coordinates.
fn lat_alt_to_parallax(lat: f64, height: f64) -> (f64, f64) {
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
/// # Arguments
/// * `lat` - Geodetic latitude of the observer in **degrees**.
/// * `height` - Observer's altitude above the reference ellipsoid in **kilometers**.
///
/// # Returns
/// A tuple `(rho_cos_phi, rho_sin_phi)`:
/// * `rho_cos_phi`: normalized distance of the observer projected on
///   the Earth's equatorial plane.
/// * `rho_sin_phi`: normalized distance of the observer projected on
///   the Earth's rotation (polar) axis.
///
/// # Details
/// This function simply converts `lat` to radians and delegates the
/// computation to [`lat_alt_to_parallax`].
///
/// # See also
/// * [`lat_alt_to_parallax`] – Performs the actual computation given latitude in radians.
pub(crate) fn geodetic_to_parallax(lat: f64, height: f64) -> (f64, f64) {
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
        let observer = Observer::new(0.0, 0.0, 0.0, None, None, None);
        assert_eq!(observer.longitude, 0.0);
        assert_eq!(observer.rho_cos_phi, 1.0);
        assert_eq!(observer.rho_sin_phi, 0.0);

        let observer = Observer::new(
            289.25058,
            -30.2446,
            2647.,
            Some("Rubin Observatory".to_string()),
            Some(NotNan::new(0.0001).unwrap()),
            Some(NotNan::new(0.0001).unwrap()),
        );

        assert_eq!(observer.longitude, 289.25058);
        assert_eq!(observer.rho_cos_phi, 0.8649760504617418);
        assert_eq!(observer.rho_sin_phi, -0.5009551027512434);
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, None, None, None);
        let obs_fixed_vector = pan_starrs.body_fixed_coord();
        assert_eq!(
            obs_fixed_vector,
            Vector3::new(
                -0.00003653799439776371,
                -0.00001607260397528885,
                0.000014988110430544328
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
        let pan_starrs = Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None);

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
        let pan_starrs = Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None);

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
