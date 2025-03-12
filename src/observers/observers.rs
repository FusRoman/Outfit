use nalgebra::Vector3;

use crate::constants::ERAU;

use super::super::jpl_request::observer_pos::geodetic_to_parallax;

pub struct Observer {
    // in degrees east of Greenwich
    pub longitude: f64,
    // phi' is the geocentric latitude and rho is the geocentric distance in earth radii
    // rho cos phi
    pub rho_cos_phi: f64,
    // rho sin phi
    pub rho_sin_phi: f64,
    pub name: Option<String>,
}

impl Observer {
    pub fn new(longitude: f64, latitude: f64, elevation: f64, name: Option<String>) -> Observer {
        let (rho_cos_phi, rho_sin_phi) = geodetic_to_parallax(latitude, elevation);
        Observer {
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            name,
        }
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
    pub(crate) fn body_fixed_coord(&self) -> Vector3<f64> {
        let lon_radians = self.longitude.to_radians();

        Vector3::new(
            ERAU * self.rho_cos_phi * lon_radians.cos(),
            ERAU * self.rho_cos_phi * lon_radians.sin(),
            ERAU * self.rho_sin_phi,
        )
    }
}

#[cfg(test)]
mod observer_test {

    use super::*;

    #[test]
    fn test_observer_constructor() {
        let observer = Observer::new(0.0, 0.0, 0.0, None);
        assert_eq!(observer.longitude, 0.0);
        assert_eq!(observer.rho_cos_phi, 1.0);
        assert_eq!(observer.rho_sin_phi, 0.0);

        let observer = Observer::new(
            289.25058,
            -30.2446,
            2647.,
            Some("Rubin Observatory".to_string()),
        );

        assert_eq!(observer.longitude, 289.25058);
        assert_eq!(observer.rho_cos_phi, 0.8649760504617418);
        assert_eq!(observer.rho_sin_phi, -0.5009551027512434);
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, None);
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
}
