use std::{collections::HashMap, sync::Arc};

use once_cell::sync::OnceCell;

use ordered_float::NotNan;

use crate::{
    constants::{Degree, Kilometer, MpcCode, MpcCodeObs},
    env_state::OutfitEnv,
    error_models::{get_bias_rms, ErrorModel, ErrorModelData},
    jpl_ephem::download_jpl_file::EphemFileSource,
    observers::{observatories::Observatories, observers::Observer},
    outfit_errors::OutfitError,
};

use crate::jpl_ephem::JPLEphem;

#[derive(Debug, Clone)]
pub struct Outfit {
    env_state: OutfitEnv,
    observatories: Observatories,
    jpl_source: EphemFileSource,
    jpl_ephem: OnceCell<JPLEphem>,
    error_model: ErrorModelData,
}

impl Outfit {
    pub fn new(jpl_file: &str, error_model: ErrorModel) -> Result<Self, OutfitError> {
        Ok(Outfit {
            env_state: OutfitEnv::new(),
            observatories: Observatories::new(),
            jpl_source: jpl_file.try_into()?,
            jpl_ephem: OnceCell::new(),
            error_model: error_model.read_error_model_file()?,
        })
    }

    pub fn get_jpl_ephem(&self) -> Result<&JPLEphem, OutfitError> {
        self.jpl_ephem
            .get_or_try_init(|| JPLEphem::new(&self.jpl_source))
    }

    pub(crate) fn get_ut1_provider(&self) -> &hifitime::ut1::Ut1Provider {
        &self.env_state.ut1_provider
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
            .unwrap_or_else(|| panic!("MPC code not found: {}", mpc_code))
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

                // TODO: need to handle the catalog code
                // For now, we assume the catalog code is "c" (for "catalog")
                // This is a simplification, as the catalog code can vary
                // depending on the observer and the error model
                let bias_rms = get_bias_rms(&self.error_model, code.to_string(), "c".to_string());

                let observer = Observer {
                    longitude: NotNan::try_from(longitude as f64).expect("Longitude cannot be NaN"),
                    rho_cos_phi: NotNan::try_from(cos as f64).expect("Longitude cannot be NaN"),
                    rho_sin_phi: NotNan::try_from(sin as f64).expect("Longitude cannot be NaN"),
                    name: Some(name),
                    ra_accuracy: bias_rms.map(|(ra, _)| NotNan::try_from(ra as f64).unwrap()),
                    dec_accuracy: bias_rms.map(|(_, dec)| NotNan::try_from(dec as f64).unwrap()),
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
        .unwrap_or_else(|| panic!("Failed to parse float for observer code: {code}"))
        .trim()
        .parse()
}

fn parse_remain(remain: &str, code: &str) -> (f32, f32, f32, String) {
    let name = remain
        .get(27..)
        .unwrap_or_else(|| panic!("Failed to parse name value for code: {code}"));

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
    use super::*;
    use ordered_float::NotNan;

    use crate::{observers::observers::Observer, outfit::Outfit};

    #[test]
    fn test_observer_from_mpc_code() {
        let outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

        let observer = outfit.get_observer_from_mpc_code(&"000".into());
        let test = Observer {
            longitude: NotNan::try_from(0.).unwrap(),
            rho_cos_phi: NotNan::try_from(0.6241099834442139).unwrap(),
            rho_sin_phi: NotNan::try_from(0.7787299752235413).unwrap(),
            name: Some("Greenwich".to_string()),
            ra_accuracy: Some(NotNan::try_from(0.5099999904632568).unwrap()),
            dec_accuracy: Some(NotNan::try_from(0.4000000059604645).unwrap()),
        };
        assert_eq!(observer, test.into());

        let observer = outfit.get_observer_from_mpc_code(&"C51".into());
        let test = Observer {
            longitude: NotNan::try_from(0.).unwrap(),
            rho_cos_phi: NotNan::try_from(0.).unwrap(),
            rho_sin_phi: NotNan::try_from(0.).unwrap(),
            name: Some("WISE".to_string()),
            ra_accuracy: Some(NotNan::try_from(0.5099999904632568).unwrap()),
            dec_accuracy: Some(NotNan::try_from(0.4000000059604645).unwrap()),
        };
        assert_eq!(observer, test.into());

        let observer = outfit.get_observer_from_mpc_code(&"Z50".into());
        let test = Observer {
            longitude: NotNan::try_from(355.2843017578125).unwrap(),
            rho_cos_phi: NotNan::try_from(0.7440530061721802).unwrap(),
            rho_sin_phi: NotNan::try_from(0.666068971157074).unwrap(),
            name: Some("Mazariegos".to_string()),
            ra_accuracy: Some(NotNan::try_from(0.5099999904632568).unwrap()),
            dec_accuracy: Some(NotNan::try_from(0.4000000059604645).unwrap()),
        };
        assert_eq!(observer, test.into());
    }

    #[test]
    fn test_add_observer() {
        let mut outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let obs = outfit.new_observer(1.0, 2.0, 3.0, Some("Test".to_string()));
        assert_eq!(obs.longitude, 1.0);
        assert_eq!(obs.rho_cos_phi, 0.999395371426802);
        assert_eq!(obs.rho_sin_phi, 0.0346660237964843);
        assert_eq!(obs.name, Some("Test".to_string()));
        assert_eq!(outfit.observatories.uint16_from_observer(obs), 0);

        let obs2 = outfit.new_observer(4.0, 5.0, 6.0, Some("Test2".to_string()));
        assert_eq!(outfit.observatories.uint16_from_observer(obs2), 1);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_get_jpl_ephem_from_horizon() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;
        let jpl_ephem = OUTFIT_HORIZON_TEST.0.get_jpl_ephem();
        assert!(jpl_ephem.is_ok(), "Failed to get JPL ephemeris file");
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_get_jpl_ephem_from_naif() {
        use crate::unit_test_global::OUTFIT_NAIF_TEST;
        let jpl_ephem = OUTFIT_NAIF_TEST.get_jpl_ephem();
        assert!(jpl_ephem.is_ok(), "Failed to get JPL ephemeris file");
    }
}
