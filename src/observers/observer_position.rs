use crate::constants::MJD;
use crate::outfit::Outfit;

use super::super::observers::observers::Observer;
use nalgebra::{Matrix3, Vector3};

use super::super::constants::{DPI, EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS, RADSEC, T2000};
use super::super::ref_system::{nutn80, obleq, rotmt, rotpn};
use hifitime::prelude::Epoch;
use hifitime::ut1::Ut1Provider;

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
///                (longitude, normalized radius components) of one observer.
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
///     1. Computes the observer’s geocentric position in the ecliptic mean J2000 frame via [`pvobs`].
///     2. Retrieves the heliocentric Earth position from the JPL ephemeris.
///     3. Transforms the observer’s position from the ecliptic mean J2000 frame to the equatorial mean J2000 frame using [`rotpn`].
///     4. Computes the heliocentric observer position by summing the Earth and transformed observer vectors.
///
/// # Validation
/// This function is tested in [`observer_pos_tests::test_helio_pos_obs`] and validated against known heliocentric
/// positions from authoritative sources such as JPL Horizons and OrbFit for real observatories.
///
/// # See also
/// * [`pvobs`] – computes the observer’s geocentric position in the ecliptic mean J2000 frame
/// * [`rotpn`] – transforms vectors between celestial reference frames (e.g., ECLM to EQUM)
/// * [`Outfit::get_jpl_ephem`] – accesses planetary ephemerides for Earth position
/// * [`Outfit::get_ut1_provider`] – provides Earth orientation parameters (e.g., ΔUT1)
pub fn helio_obs_pos(
    observers: [&Observer; 3],
    mjd_tt: &Vector3<MJD>,
    state: &Outfit,
) -> Matrix3<f64> {
    let ut1_provider = state.get_ut1_provider();
    let jpl = state.get_jpl_ephem().unwrap();

    let epochs = [
        Epoch::from_mjd_in_time_scale(mjd_tt.x, hifitime::TimeScale::TT),
        Epoch::from_mjd_in_time_scale(mjd_tt.y, hifitime::TimeScale::TT),
        Epoch::from_mjd_in_time_scale(mjd_tt.z, hifitime::TimeScale::TT),
    ];

    let pos_obs = [
        pvobs(observers[0], &epochs[0], &ut1_provider).0,
        pvobs(observers[1], &epochs[1], &ut1_provider).0,
        pvobs(observers[2], &epochs[2], &ut1_provider).0,
    ];
    let pos_obs_matrix = Matrix3::from_columns(&pos_obs);

    let earth_pos = [
        jpl.earth_ephemeris(&epochs[0], false).0,
        jpl.earth_ephemeris(&epochs[1], false).0,
        jpl.earth_ephemeris(&epochs[2], false).0,
    ];
    let earth_pos_matrix = Matrix3::from_columns(&earth_pos);

    let mut rot = [[0.0; 3]; 3];
    rotpn(&mut rot, "ECLM", "J2000", 0.0, "EQUM", "J2000", 0.0);
    let rot_matrix = Matrix3::from(rot).transpose();

    earth_pos_matrix + rot_matrix * pos_obs_matrix
}

pub(crate) fn geo_obs_pos(
    observer: &Observer,
    tmjd: &Epoch,
    ut1_provider: &Ut1Provider,
) -> (Vector3<f64>, Vector3<f64>) {
    // Initialisation
    let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

    // Get the coordinates of the observer on Earth
    let dxbf = observer.body_fixed_coord();

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
    let rot_mat = Matrix3::from(rot).transpose();

    let obs_pos = rot_mat * dxbf;
    let obs_vel = rot_mat * dvbf;

    // Transformation in the ecliptic mean J2000
    let mut rot1 = [[0.; 3]; 3];
    rotpn(
        &mut rot1,
        "EQUT",
        "OFDATE",
        tmjd.to_mjd_tt_days(),
        "ECLM",
        "J2000",
        0.,
    );
    let rot1_mat = Matrix3::from(rot1).transpose();
    let obs_pos = rot1_mat * obs_pos;
    let obs_vel = rot1_mat * obs_vel;
    (obs_pos, obs_vel)
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
///     - `dx`: observer geocentric position vector in ecliptic mean J2000 frame [AU].
///     - `dv`: observer velocity vector due to Earth's rotation, in the same frame [AU/day].
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
    observer: &Observer,
    tmjd: &Epoch,
    ut1_provider: &Ut1Provider,
) -> (Vector3<f64>, Vector3<f64>) {
    // Initialisation
    let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

    // Get the coordinates of the observer on Earth
    let dxbf = observer.body_fixed_coord();

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

    // Transformation in the ecliptic mean J2000
    let mut rot1 = [[0.; 3]; 3];
    rotpn(
        &mut rot1,
        "EQUT",
        "OFDATE",
        tmjd.to_mjd_tt_days(),
        "ECLM",
        "J2000",
        0.,
    );

    let rot1_mat = Matrix3::from(rot1).transpose();
    let rot_mat = Matrix3::from(rot).transpose();

    let rotmat = rot1_mat * rot_mat;

    // Apply transformation to the observer position and velocity
    let dx = rotmat * dxbf;
    let dv = rotmat * dvbf;

    (dx, dv)
}

/// Compute the Greenwich Mean Sidereal Time (GMST)
/// in radians, for a modified julian date (UT1).
///
/// Arguments
/// ---------
/// * `tjm` - Modified Julian Date (MJD)
///
/// Retour
/// ------
/// GMST in radians, in the [0, 2π) interval.
fn gmst(tjm: f64) -> f64 {
    /// Coefficients du polynôme pour GMST à 0h UT1 (en secondes)
    const C0: f64 = 24110.54841;
    const C1: f64 = 8640184.812866;
    const C2: f64 = 9.3104e-2;
    const C3: f64 = -6.2e-6;
    /// Rapport entre jour sidéral et jour solaire
    const RAP: f64 = 1.00273790934;

    let itjm = tjm.floor();
    let t = (itjm - T2000) / 36525.0;

    // Calcul du temps sidéral moyen à 0h UT1
    let mut gmst0 = ((C3 * t + C2) * t + C1) * t + C0;

    gmst0 *= DPI / 86400.0;

    // Incrément de GMST à partir de 0h
    // let h = (57028.476562500000 - 57028.0) * DPI;
    let h = tjm.fract() * DPI;
    let mut gmst = gmst0 + h * RAP;

    // Ajustement pour rester dans [0, 2π]
    let mut i: i64 = (gmst / DPI).floor() as i64;
    if gmst < 0.0 {
        i = i - 1;
    }
    gmst -= i as f64 * DPI;

    gmst
}

/// Compute the equinox equation
///
/// Arguments
/// ---------
/// * `tjm`: Modified Julian Date (MJD)
///
/// Retour
/// ------
/// * Equinox equation in radians
fn equequ(tjm: f64) -> f64 {
    let oblm = obleq(tjm);
    let (dpsi, _deps) = nutn80(tjm);
    RADSEC * dpsi * oblm.cos()
}

/// Convert latitude and height in parallax coordinates on the Earth
///
/// Argument
/// --------
/// * `lat`: observer latitude in radians
/// * `height`: observer height in kilometer
///
/// Return
/// ------
/// * rho_cos_phi: normalized radius of the observer projected on the equatorial plane
/// * rho_sin_phi: normalized radius of the observer projected on the polar axis.
fn lat_alt_to_parallax(lat: f64, height: f64) -> (f64, f64) {
    let axis_ratio = EARTH_MINOR_AXIS / EARTH_MAJOR_AXIS;
    let u = (lat.sin() * axis_ratio).atan2(lat.cos());

    let rho_sin_phi = axis_ratio * u.sin() + (height / EARTH_MAJOR_AXIS) * lat.sin();
    let rho_cos_phi = u.cos() + (height / EARTH_MAJOR_AXIS) * lat.cos();

    (rho_cos_phi, rho_sin_phi)
}

/// Convert latitude in degree and height in parallax coordinate
///
/// Argument
/// --------
/// * `lat`: observer latitude in degree
/// * `height`: observer height in kilometer
///
/// Return
/// ------
/// * `rho_cos_phi`: normalized radius of the observer projected on the equatorial plane
/// * `rho_sin_phi`: normalized radius of the observer projected on the polar axis.
pub(crate) fn geodetic_to_parallax(lat: f64, height: f64) -> (f64, f64) {
    let latitude_rad = lat.to_radians();

    let (rho_cos_phi, rho_sin_phi) = lat_alt_to_parallax(latitude_rad, height);

    (rho_cos_phi, rho_sin_phi)
}

#[cfg(test)]
mod observer_pos_tests {

    use super::*;
    use crate::{error_models::ErrorModel, outfit::Outfit};

    #[test]
    fn geodetic_to_parallax_test() {
        // latitude and height of Pan-STARRS 1, Haleakala
        let (pxy1, pz1) = geodetic_to_parallax(20.707233557, 3067.694);
        assert_eq!(pxy1, 0.9362410003211518);
        assert_eq!(pz1, 0.35154299856304305);
    }

    #[test]
    fn test_gmst() {
        let tut = 57028.478514610404;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.851925725092499);

        let tut = T2000;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.894961212789145);
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
            pvobs(&pan_starrs, &epoch, state.get_ut1_provider());

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
    #[cfg(feature = "jpl-download")]
    fn test_helio_pos_obs() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let tmjd = Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593);

        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None);

        // Now we need a Vector3<Observer> with three identical copies
        let observers = [&pan_starrs, &pan_starrs, &pan_starrs];

        let helio_pos = helio_obs_pos(observers, &tmjd, &OUTFIT_HORIZON_TEST.0);

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
