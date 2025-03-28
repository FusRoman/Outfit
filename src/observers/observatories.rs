use super::bimap::BiMap;
use super::observers::Observer;
use crate::constants::MpcCodeObs;
use once_cell::sync::OnceCell;
use std::sync::Arc;

#[derive(Debug)]
pub(crate) struct Observatories {
    pub (crate) mpc_code_obs: OnceCell<MpcCodeObs>,
    obs_to_uint16: BiMap<Arc<Observer>, u16>,
}

impl Observatories {
    pub(crate) fn new() -> Self {
        Observatories {
            mpc_code_obs: OnceCell::new(),
            obs_to_uint16: BiMap::new(),
        }
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
