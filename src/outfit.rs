use std::{collections::HashMap, sync::Arc};

use ordered_float::NotNan;
use ureq::http::{self, Uri};

use crate::{
    constants::{Degree, Kilometer, MpcCode, MpcCodeObs},
    env_state::OutfitEnv,
    observers::{observatories::Observatories, observers::Observer},
};

pub struct Outfit {
    env_state: OutfitEnv,
    observatories: Observatories,
}

impl Outfit {
    pub fn new() -> Self {
        Outfit {
            env_state: OutfitEnv::new(),
            observatories: Observatories::new(),
        }
    }

    pub(crate) fn get_ut1_provider(&self) -> &hifitime::ut1::Ut1Provider {
        &self.env_state.ut1_provider
    }

    pub(crate) fn post_url<T, I, K, V>(&self, url: T, form: I) -> String
    where
        Uri: TryFrom<T>,
        <Uri as TryFrom<T>>::Error: Into<http::Error>,
        I: IntoIterator<Item = (K, V)>,
        K: AsRef<str>,
        V: AsRef<str>,
    {
        self.env_state.post_from_url(url, form)
    }

    /// Get the observatories from the Minor Planet Center
    ///
    /// Return
    /// ------
    /// * A map of observatories from the Minor Planet Center
    ///    The key is the MPC code and the value is the observer
    pub(crate) fn get_observatories(&self) -> &MpcCodeObs {
        self.observatories
            .mpc_code_obs
            .get_or_init(|| self.init_observatories())
    }

    /// Get an observer from an MPC code
    ///
    /// Arguments
    /// ---------
    /// * `mpc_code`: the MPC code
    ///
    /// Return
    /// ------
    /// * The observer
    pub fn get_observer_from_mpc_code(&self, mpc_code: &MpcCode) -> Arc<Observer> {
        self.get_observatories()
            .get(mpc_code)
            .expect(format!("MPC code not found: {}", mpc_code).as_str())
            .clone()
    }

    /// Initialize the observatories map from the Minor Planet Center
    /// The map is a lazy map and is only loaded once
    ///
    /// Return
    /// ------
    /// * A map of observatories from the Minor Planet Center
    ///    The key is the MPC code and the value is the observer
    pub(crate) fn init_observatories(&self) -> MpcCodeObs {
        let mut observatories: MpcCodeObs = HashMap::new();

        let mpc_code_response = self
            .env_state
            .get_from_url("https://www.minorplanetcenter.net/iau/lists/ObsCodes.html");

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

                let observer = Observer {
                    longitude: NotNan::try_from(longitude as f64).expect("Longitude cannot be NaN"),
                    rho_cos_phi: NotNan::try_from(cos as f64).expect("Longitude cannot be NaN"),
                    rho_sin_phi: NotNan::try_from(sin as f64).expect("Longitude cannot be NaN"),
                    name: Some(name),
                };
                observatories.insert(code.to_string(), Arc::new(observer));
            };
        }
        observatories
    }

    /// Get an observatory index from an MPC code
    ///
    /// Argument
    /// --------
    /// * `mpc_code`: MPC code
    ///
    /// Return
    /// ------
    /// * The observatory index
    pub(crate) fn uint16_from_mpc_code(&mut self, mpc_code: &MpcCode) -> u16 {
        let observer = self.get_observer_from_mpc_code(mpc_code);
        self.observatories.uint16_from_observer(observer)
    }

    /// Get an observer from an observer index
    ///
    /// Arguments
    /// ---------
    /// * `observer_idx`: the observer index
    ///
    /// Return
    /// ------
    /// * The observer
    pub(crate) fn uint16_from_observer(&mut self, observer: Arc<Observer>) -> u16 {
        self.observatories.uint16_from_observer(observer)
    }

    /// Get an observer from an observer index
    ///
    /// Arguments
    /// ---------
    /// * `observer_idx`: the observer index
    ///
    /// Return
    /// ------
    /// * The observer
    pub(crate) fn get_observer_from_uint16(&self, observer_idx: u16) -> &Observer {
        self.observatories.get_observer_from_uint16(observer_idx)
    }

    pub fn new_observer(
        &mut self,
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
    ) -> Arc<Observer> {
        self.observatories
            .add_observer(longitude, latitude, elevation, name)
    }
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

#[cfg(test)]
mod outfit_struct_test {
    use ordered_float::NotNan;

    use crate::{observers::observers::Observer, outfit::Outfit};

    #[tokio::test(flavor = "multi_thread", worker_threads = 2)]
    async fn test_observer_from_mpc_code() {
        let outfit = Outfit::new();

        let observer = outfit.get_observer_from_mpc_code(&"000".into());
        let test = Observer {
            longitude: NotNan::try_from(0.).unwrap(),
            rho_cos_phi: NotNan::try_from(0.6241099834442139).unwrap(),
            rho_sin_phi: NotNan::try_from(0.7787299752235413).unwrap(),
            name: Some("Greenwich".to_string()),
        };
        assert_eq!(observer, test.into());

        let observer = outfit.get_observer_from_mpc_code(&"C51".into());
        let test = Observer {
            longitude: NotNan::try_from(0.).unwrap(),
            rho_cos_phi: NotNan::try_from(0.).unwrap(),
            rho_sin_phi: NotNan::try_from(0.).unwrap(),
            name: Some("WISE".to_string()),
        };
        assert_eq!(observer, test.into());

        let observer = outfit.get_observer_from_mpc_code(&"Z50".into());
        let test = Observer {
            longitude: NotNan::try_from(355.2843017578125).unwrap(),
            rho_cos_phi: NotNan::try_from(0.7440530061721802).unwrap(),
            rho_sin_phi: NotNan::try_from(0.666068971157074).unwrap(),
            name: Some("Mazariegos".to_string()),
        };
        assert_eq!(observer, test.into());
    }

    #[test]
    fn test_add_observer() {
        let mut outfit = Outfit::new();
        let obs = outfit.new_observer(1.0, 2.0, 3.0, Some("Test".to_string()));
        assert_eq!(obs.longitude, 1.0);
        assert_eq!(obs.rho_cos_phi, 0.999395371426802);
        assert_eq!(obs.rho_sin_phi, 0.0346660237964843);
        assert_eq!(obs.name, Some("Test".to_string()));
        assert_eq!(outfit.observatories.uint16_from_observer(obs), 0);

        let obs2 = outfit.new_observer(4.0, 5.0, 6.0, Some("Test2".to_string()));
        assert_eq!(outfit.observatories.uint16_from_observer(obs2), 1);
    }
}
