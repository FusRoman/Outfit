use nalgebra::{Matrix3, Vector3};

use super::super::constants::{
    DPI, EARTH_MAJOR_AXIS, EARTH_MINOR_AXIS, ERAU, RADSEC, SECONDS_PER_DAY, T2000,
};
use super::super::ref_system::{nutn80, obleq, rotmt, rotpn};

// Structure pour stocker la position et la vitesse
#[derive(Debug)]
struct ObserverState {
    position: Vector3<f64>,
    velocity: Vector3<f64>,
}

// Fonction principale
fn pvobs(tmjd: f64, longitude: f64, latitude: f64, height: f64) -> ObserverState {
    // Initialisation
    let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);
    let mut dx = Vector3::zeros();
    let mut dv = Vector3::zeros();

    // Récupération des coordonnées de l'observatoire
    let dxbf = body_fixed_coord(longitude, latitude, height);

    // Calcul de la vitesse due à la rotation terrestre
    let dvbf = omega.cross(&dxbf);

    println!("dxbf: {dxbf}");
    println!("dvbf: {dvbf}");

    // TODO: conversion from ET to UT1 not working here
    // using the crate https://docs.rs/hifitime/latest/hifitime/index.html
    // could resolve this problem

    // see line 335 in reference_system.f90 
    // and subroutine cnvtim in time_scale.f90 line 564
    // in orbfit
    let (mjd1, sec1) = (tmjd.floor(), (tmjd - tmjd.floor()) * SECONDS_PER_DAY);
    let tut = (mjd1 + sec1 / SECONDS_PER_DAY) as f64;

    // Calcul du temps sidéral apparent de Greenwich
    let gast = gmst(tut) + equequ(tmjd);

    // Matrice de rotation de la Terre
    let rot = rotmt(-gast, 2);

    // Transformation vers le référentiel J2000
    let mut rot1 = [[0.; 3]; 3];
    rotpn(&mut rot1, "EQUT", "OFDATE", tmjd, "ECLM", "J2000", 0.);
    let rotmat = Matrix3::from(rot1) * Matrix3::from(rot);

    // Application des transformations
    dx = rotmat * dxbf;
    dv = rotmat * dvbf;

    println!("t: {tmjd}");
    println!("omega: {omega}");
    println!("dxbf: {dxbf}");
    println!("dvbf: {dvbf}");
    println!("mjd1: {mjd1}");
    println!("sec1: {sec1}");
    println!("tut: {tut}");
    println!("gast: {gast}");
    println!("rot: {rot:?}");
    println!("rot1: {rot1:?}");
    println!("rotmat: {rotmat}");

    ObserverState {
        position: dx,
        velocity: dv,
    }
}

/// Compute the Greenwich Mean Sidereal Time (GMST)
/// in radians, for a modified julian date (UT1).
///
/// # Arguments
/// * `tjm` - Modified Julian Date (MJD)
///
/// # Retour
/// GMST in radians, in the [0, 2π) interval.
fn gmst(tjm: f64) -> f64 {
    /// Coefficients du polynôme pour GMST à 0h UT1 (en secondes)
    const C0: f64 = 24110.54841;
    const C1: f64 = 8640184.812866;
    const C2: f64 = 9.3104e-2;
    const C3: f64 = -6.2e-6;
    /// Rapport entre jour sidéral et jour solaire
    const RAP: f64 = 1.00273790934;

    // Partie entière de tjm correspondant à 0h UT1
    let itjm = tjm.floor();
    // Nombre de siècles julien écoulés depuis J2000
    let t = (itjm - T2000) / 36525.0;

    // GMST à 0h UT1 en secondes via le polynôme, converti en radians
    let gmst0 = (((C3 * t + C2) * t + C1) * t + C0) * DPI / 86400.0;

    // Incrément de GMST dû à la fraction de jour (convertie en radians)
    let h = (tjm - itjm) * DPI;

    // Calcul brut du GMST avec correction UT1
    let raw_gmst = gmst0 + h * RAP;

    // Normalisation dans l'intervalle [0, 2π)
    raw_gmst.rem_euclid(DPI)
}

/// Compute the equinox equation
///
/// # Arguments
/// * `tjm`: Modified Julian Date (MJD)
///
/// # Retour
/// * Equinox equation in radians
fn equequ(tjm: f64) -> f64 {
    let oblm = obleq(tjm);
    let (dpsi, _deps) = nutn80(tjm);
    RADSEC * dpsi * oblm.cos()
}

// Fonction de génération d'une matrice de rotation
fn rotation_matrix(angle: f64, axis: usize) -> Matrix3<f64> {
    let cos_a = angle.to_radians().cos();
    let sin_a = angle.to_radians().sin();

    match axis {
        3 => Matrix3::new(cos_a, sin_a, 0.0, -sin_a, cos_a, 0.0, 0.0, 0.0, 1.0),
        _ => Matrix3::identity(),
    }
}

// Fonction de transformation de référentiel (simplifiée)
fn transformation_matrix(
    _from: &str,
    _to: &str,
    _t: f64,
    _frame1: &str,
    _frame2: &str,
) -> Matrix3<f64> {
    Matrix3::identity() // Placeholder
}

/// Convert latitude and height in parallax coordinates on the the Earth
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
/// * lat: observer latitude in degree
/// * height: observer height in kilometer
///
/// Return
/// ------
/// * rho_cos_phi: normalized radius of the observer projected on the equatorial plane
/// * rho_sin_phi: normalized radius of the observer projected on the polar axis.
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

    use super::*;

    #[test]
    fn geodetic_to_parallax_test() {
        /// latitude and height of Pan-STARRS 1, Haleakala
        let (pxy1, pz1) = geodetic_to_parallax(20.707233557, 3067.694);
        assert_eq!(pxy1, 0.9362410003211518);
        assert_eq!(pz1, 0.35154299856304305);
    }

    #[test]
    fn body_fixed_coord_test() {
        /// longitude, latitude and height of Pan-STARRS 1, Haleakala
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
    fn pvobs_test() {
        let tmjd = 57028.479297592596;
        /// longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);

        let observer_state = pvobs(tmjd, lon, lat, h);

        println!("Position (J2000 ecliptic): {:?}", observer_state.position);
        println!("Velocity (J2000 ecliptic): {:?}", observer_state.velocity);
    }
}
