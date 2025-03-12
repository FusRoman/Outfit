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
}

#[cfg(test)]
mod observer_test {

    #[test]
    fn test_observer_constructor() {
        use super::Observer;
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
}
