use nalgebra::{Matrix3, Vector3};
use super::super::constants::{DPI, EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS, ERAU, RADSEC, T2000};
use super::super::env_state::OutfitState;
use super::super::ref_system::{nutn80, obleq, rotmt, rotpn};
use super::earth_pos::get_earth_position;
use hifitime::prelude::Epoch;
use hifitime::ut1::Ut1Provider;


/// Get the matrix of the heliocentric position vector at the time of three observations
/// 
/// Argument
/// --------
/// * mjd_tt: vector of observation time in modified julian date (MJD)
/// * longitude: observer longitude on Earth in degree
/// * latitude: observer latitude on Earth in degree
/// * height: observer height on Earth in degree
/// * state: need the http_client and the ut1_provider
/// 
/// Return
/// ------
/// * a 3x3 matrix containing the x,y,z coordinates of the observer at the time of the three 
///     observations (reference frame: Equatorial mean J2000, units: AU)
pub async fn helio_obs_pos(
    mjd_tt: &Vector3<f64>,
    longitude: f64,
    latitude: f64,
    height: f64,
    state: &OutfitState,
) -> Matrix3<f64> {
    let position_obs_time = mjd_tt
        .iter()
        .map(|mjd_el| pvobs(*mjd_el, longitude, latitude, height, &state.ut1_provider).0)
        .collect::<Vec<Vector3<f64>>>();

    let pos_obs_matrix = Matrix3::from_columns(&position_obs_time);

    let time_vec = vec![mjd_tt.x, mjd_tt.y, mjd_tt.z];
    let earth_helio_position = get_earth_position(&time_vec, &state.http_client).await;

    let earth_pos_vec = earth_helio_position
        .iter()
        .map(|pos_record| pos_record.pos_vector())
        .collect::<Vec<Vector3<f64>>>();

    let earth_pos_matrix = Matrix3::from_columns(&earth_pos_vec);

    let mut rot = [[0.; 3]; 3];
    rotpn(&mut rot, "ECLM", "J2000", 0., "EQUM", "J2000", 0.);
    let rot_matrix = Matrix3::from(rot).transpose();
    let dx = rot_matrix * pos_obs_matrix;
    earth_pos_matrix + dx
}

/// Get the observer position and velocity on the Earth
///
/// Argument
/// --------
/// * tmjd: time of the observation in modified julian date (MJD)
/// * longitude: observer longitude on Earth in degree
/// * latitude: observer latitude on Earth in degree
/// * height: observer height on Earth in degree
/// * ut1_provider: the ut1 provider from hifitime containing the delta time in second between TAI and UT1 from the JPL
///
/// Return
/// ------
/// * dx: corrected observer position with respect to the center of mass of Earth (in ecliptic J2000)
/// * dy: corrected observer velocity with respect to the center of mass of Earth (in ecliptic J2000)
fn pvobs(
    tmjd: f64,
    longitude: f64,
    latitude: f64,
    height: f64,
    ut1_provider: &Ut1Provider,
) -> (Vector3<f64>, Vector3<f64>) {
    // Initialisation
    let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

    // Get the coordinates of the observer on Earth
    let dxbf = body_fixed_coord(longitude, latitude, height);

    // Get the observer velocity due to Earth rotation
    let dvbf = omega.cross(&dxbf);

    // deviation from Orbfit, use of another conversion from MJD UTC (ET scale) to UT1 scale
    // based on the hifitime crate
    let epoch_mjd = Epoch::from_mjd_utc(tmjd);
    let mjd_ut1 = epoch_mjd.to_ut1(ut1_provider.to_owned());
    let tut = mjd_ut1.to_mjd_utc_days();

    // Compute the Greenwich sideral apparent time
    let gast = gmst(tut) + equequ(tmjd);

    // Earth rotation matrix
    let rot = rotmt(-gast, 2);

    // Transformation in the ecliptic mean J2000
    let mut rot1 = [[0.; 3]; 3];
    rotpn(&mut rot1, "EQUT", "OFDATE", tmjd, "ECLM", "J2000", 0.);

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
/// * lat: observer latitude in radians
/// * height: observer height in kilometer
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
fn geodetic_to_parallax(lat: f64, height: f64) -> (f64, f64) {
    let latitude_rad = lat.to_radians();

    let (rho_cos_phi, rho_sin_phi) = lat_alt_to_parallax(latitude_rad, height);

    (rho_cos_phi, rho_sin_phi)
}

/// Get the fixed position of an observatory using its geographic coordinates
///
/// Argument
/// --------
/// * longitude: observer longitude in degree
/// * latitude: observer latitude in degree
/// * height: observer height in degree
///
/// Return
/// ------
/// * observer fixed coordinates vector on the Earth (not corrected from Earth motion)
/// * units is AU
fn body_fixed_coord(longitude: f64, latitude: f64, height: f64) -> Vector3<f64> {
    let (pxy1, pz1) = geodetic_to_parallax(latitude, height);
    let lon_radians = longitude.to_radians();

    Vector3::new(
        ERAU * pxy1 * lon_radians.cos(),
        ERAU * pxy1 * lon_radians.sin(),
        ERAU * pz1,
    )
}

#[cfg(test)]
mod observer_pos_tests {

    use super::super::super::env_state::OutfitState;
    use super::*;

    #[test]
    fn geodetic_to_parallax_test() {
        // latitude and height of Pan-STARRS 1, Haleakala
        let (pxy1, pz1) = geodetic_to_parallax(20.707233557, 3067.694);
        assert_eq!(pxy1, 0.9362410003211518);
        assert_eq!(pz1, 0.35154299856304305);
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let obs_fixed_vector = body_fixed_coord(lon, lat, h);
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
    fn test_gmst() {
        let tut = 57028.478514610404;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.851925725092499);

        let tut = T2000;
        let res_gmst = gmst(tut);
        assert_eq!(res_gmst, 4.894961212789145);
    }

    #[tokio::test]
    async fn pvobs_test() {
        let state = OutfitState::new().await;
        let tmjd = 57028.479297592596;
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);

        let (observer_position, observer_velocity) = pvobs(tmjd, lon, lat, h, &state.ut1_provider);

        assert_eq!(
            observer_position.as_slice(),
            [
                -2.1029664445055886e-5,
                3.7089965349631534e-5,
                2.911548164794497e-7
            ]
        );
        assert_eq!(
            observer_velocity.as_slice(),
            [
                -0.00021367298085517918,
                -0.00012156695591212987,
                5.304083328775301e-5
            ]
        );
    }

    #[tokio::test]
    async fn test_helio_pos_obs() {
        let state = OutfitState::new().await;
        let tmjd = Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593);

        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);

        let helio_pos = helio_obs_pos(&tmjd, lon, lat, h, &state).await;

        assert_eq!(
            helio_pos.as_slice(),
            [
                -0.26456678469889994,
                0.8689350609363788,
                0.3766996213519332,
                -0.5891633520067185,
                0.7238873564700253,
                0.31381865190488956,
                -0.7743279501475471,
                0.5612534464162001,
                0.24333415468923875
            ]
        )
    }
}
