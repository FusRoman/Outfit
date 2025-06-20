use super::bimap::BiMap;
use super::observers::Observer;
use crate::constants::{Degree, Kilometer, MpcCodeObs};
use std::sync::{Arc, OnceLock};

#[derive(Debug)]
pub(crate) struct Observatories {
    pub(crate) mpc_code_obs: OnceLock<MpcCodeObs>,
    obs_to_uint16: BiMap<Arc<Observer>, u16>,
}

impl Observatories {
    pub(crate) fn new() -> Self {
        Observatories {
            mpc_code_obs: OnceLock::new(),
            obs_to_uint16: BiMap::new(),
        }
    }

    pub(crate) fn add_observer(
        &mut self,
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
    ) -> Arc<Observer> {
        let obs = Observer::new(longitude, latitude, elevation, name.clone(), None, None);
        let arc_observer = Arc::new(obs);
        self.obs_to_uint16
            .entry_or_insert_by_key(arc_observer.clone(), self.obs_to_uint16.len() as u16);
        arc_observer
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
        self.obs_to_uint16
            .get_by_value(&observer_idx)
            .expect(format!("Observer index not found: {}", observer_idx).as_str())
    }

    /// Get an observer index from an observer
    /// If the observer is not already in the bimap, it is added
    ///
    /// Arguments
    /// ---------
    /// * `observer`: the observer
    ///
    /// Return
    /// ------
    /// * The observer index
    pub(crate) fn uint16_from_observer(&mut self, observer: Arc<Observer>) -> u16 {
        let obs_idx = self.obs_to_uint16.len() as u16;
        self.obs_to_uint16
            .entry_or_insert_by_key(observer, obs_idx)
            .clone()
    }
}

#[cfg(test)]
mod observatories_test {
    use super::*;

    #[test]
    fn test_observatories() {
        let mut observatories = Observatories::new();
        let obs = observatories.add_observer(1.0, 2.0, 3.0, Some("Test".to_string()));
        assert_eq!(obs.longitude, 1.0);
        assert_eq!(obs.rho_cos_phi, 0.999395371426802);
        assert_eq!(obs.rho_sin_phi, 0.0346660237964843);
        assert_eq!(obs.name, Some("Test".to_string()));
        assert_eq!(observatories.obs_to_uint16.len(), 1);
        let observer = observatories.get_observer_from_uint16(0);
        assert_eq!(observer.name, Some("Test".to_string()));

        observatories.add_observer(4.0, 5.0, 6.0, Some("Test2".to_string()));
        assert_eq!(observatories.obs_to_uint16.len(), 2);
        let observer = observatories.get_observer_from_uint16(1);
        assert_eq!(observer.name, Some("Test2".to_string()));
    }
}
