use nalgebra::Vector3;
use std::collections::HashMap;

use crate::{constants::ERAU, env_state::OutfitState};

use super::super::jpl_request::observer_pos::geodetic_to_parallax;
use crate::constants::{Degree, Kilometer, MpcCodeObs};
use serde::Deserialize;

#[derive(Debug, PartialEq, Deserialize)]
pub struct Observer {
    // in degrees east of Greenwich
    pub longitude: Degree,
    // phi' is the geocentric latitude and rho is the geocentric distance in earth radii
    // rho cos phi
    pub rho_cos_phi: Degree,
    // rho sin phi
    pub rho_sin_phi: Kilometer,
    pub name: Option<String>,
}

fn parse_f32(
    s: &str,
    slice: std::ops::Range<usize>,
    code: &str,
) -> Result<f32, std::num::ParseFloatError> {
    s.get(slice)
        .expect(format!("Failed to parse float for observer code: {code}").as_str())
        .trim()
        .parse()
}

fn parse_remain(remain: &str, code: &str) -> (f32, f32, f32, String) {
    let name = remain
        .get(27..)
        .expect(format!("Failed to parse name value for code: {code}").as_str());

    let Some(longitude) = parse_f32(remain, 1..10, code).ok() else {
        return (0.0, 0.0, 0.0, name.to_string());
    };

    let Some(cos) = parse_f32(remain, 10..18, code).ok() else {
        return (longitude, 0.0, 0.0, name.to_string());
    };

    let Some(sin) = parse_f32(remain, 18..27, code).ok() else {
        return (longitude, cos, 0.0, name.to_string());
    };
    (longitude, cos, sin, name.to_string())
}

pub(crate) async fn init_observatories(state: &OutfitState) -> MpcCodeObs {
    let mut observatories: HashMap<String, Observer> = HashMap::new();

    let mpc_code_response = state
        .http_client
        .get("https://minorplanetcenter.net/iau/lists/ObsCodes.html")
        .send()
        .await
        .expect("Request to get mpc code observatory failed")
        .text()
        .await
        .expect("Failed to get text from mpc code request");

    let mpc_code_csv = mpc_code_response
        .trim()
        .strip_prefix("<pre>")
        .and_then(|s| s.strip_suffix("</pre>"))
        .expect("Failed to strip pre tags");

    for lines in mpc_code_csv.lines().skip(2) {
        let line = lines.trim();

        if let Some((code, remain)) = line.split_at_checked(3) {
            let remain = remain.trim_end();

            let (longitude, cos, sin, name) = parse_remain(remain, code);

            observatories.insert(
                code.to_string(),
                Observer {
                    longitude: longitude as f64,
                    rho_cos_phi: cos as f64,
                    rho_sin_phi: sin as f64,
                    name: Some(name),
                },
            );
        };
    }
    observatories
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
    pub fn new(
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
    ) -> Observer {
        let (rho_cos_phi, rho_sin_phi) = geodetic_to_parallax(latitude, elevation);
        Observer {
            longitude,
            rho_cos_phi,
            rho_sin_phi,
            name,
        }
    }

    pub fn from_mpc_code<'a>(mpc_code: &str, state: &'a OutfitState) -> &'a Observer {
        let observatories = state.get_observatories();
        observatories
            .get(mpc_code)
            .expect("Failed to get observatory from mpc code")
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

    #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
    async fn test_observer_from_mpc_code() {
        let state = OutfitState::new().await;

        let observer = Observer::from_mpc_code("000", &state);
        let test = Observer {
            longitude: 0.,
            rho_cos_phi: 0.6241099834442139,
            rho_sin_phi: 0.7787299752235413,
            name: Some("Greenwich".to_string()),
        };
        assert_eq!(observer, &test);

        let observer = Observer::from_mpc_code("C51", &state);
        let test = Observer {
            longitude: 0.,
            rho_cos_phi: 0.,
            rho_sin_phi: 0.,
            name: Some("WISE".to_string()),
        };
        assert_eq!(observer, &test);

        let observer = Observer::from_mpc_code("Z50", &state);
        let test = Observer {
            longitude: 355.2843017578125,
            rho_cos_phi: 0.7440530061721802,
            rho_sin_phi: 0.666068971157074,
            name: Some("Mazariegos".to_string()),
        };
        assert_eq!(observer, &test);
    }
}
