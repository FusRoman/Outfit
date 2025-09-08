use super::bimap::BiMap;
use super::Observer;
use crate::constants::{Degree, Kilometer, MpcCodeObs};
use std::{
    fmt,
    sync::{Arc, OnceLock},
};

#[derive(Debug, Clone)]
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

    pub(crate) fn create_observer(
        &mut self,
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
    ) -> Arc<Observer> {
        let obs = Observer::new(longitude, latitude, elevation, name.clone(), None, None)
            .expect("Failed to create observer");
        let arc_observer = Arc::new(obs);
        self.obs_to_uint16
            .entry_or_insert_by_key(arc_observer.clone(), self.obs_to_uint16.len() as u16);
        arc_observer
    }

    pub(crate) fn add_observer(&mut self, observer: Arc<Observer>) -> u16 {
        let obs_idx = self.obs_to_uint16.len() as u16;
        *self.obs_to_uint16.entry_or_insert_by_key(observer, obs_idx)
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
            .unwrap_or_else(|| panic!("Observer index not found: {observer_idx}"))
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
        *self.obs_to_uint16.entry_or_insert_by_key(observer, obs_idx)
    }
}

impl fmt::Display for Observatories {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "User-defined observers:")?;
        for obs in self.obs_to_uint16.keys() {
            let (lat, height) = obs.geodetic_lat_height_wgs84();

            writeln!(
                f,
                "  {} (lon: {:.6}°, lat: {:.6}°, elev: {:.2} km)",
                obs.name.clone().unwrap_or_else(|| "Unnamed".to_string()),
                obs.longitude,
                lat,
                height
            )?;
        }

        if let Some(mpc_code_obs) = self.mpc_code_obs.get() {
            writeln!(f, "MPC observers:")?;
            for (code, obs) in mpc_code_obs.iter() {
                let (lat, height) = obs.geodetic_lat_height_wgs84();

                writeln!(
                    f,
                    "  {} [{}] (lon: {:.6}°, lat: {:.6}°, elev: {:.2} km)",
                    obs.name.clone().unwrap_or_else(|| "Unnamed".to_string()),
                    code,
                    obs.longitude,
                    lat,
                    height
                )?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod observatories_test {
    use super::*;

    #[test]
    fn test_observatories() {
        let mut observatories = Observatories::new();
        let obs = observatories.create_observer(1.0, 2.0, 3.0, Some("Test".to_string()));
        assert_eq!(obs.longitude, 1.0);
        assert_eq!(obs.rho_cos_phi, 0.999395371426802);
        assert_eq!(obs.rho_sin_phi, 0.0346660237964843);
        assert_eq!(obs.name, Some("Test".to_string()));
        assert_eq!(observatories.obs_to_uint16.len(), 1);
        let observer = observatories.get_observer_from_uint16(0);
        assert_eq!(observer.name, Some("Test".to_string()));

        observatories.create_observer(4.0, 5.0, 6.0, Some("Test2".to_string()));
        assert_eq!(observatories.obs_to_uint16.len(), 2);
        let observer = observatories.get_observer_from_uint16(1);
        assert_eq!(observer.name, Some("Test2".to_string()));
    }

    #[cfg(test)]
    mod observatories_display_tests {
        use super::*;

        /// Ensure the "user-defined" section is printed and includes both user observers.
        ///
        /// Notes
        /// -----
        /// * We don't assume any iteration order (HashMap-backed bi-map).
        /// * We check for the header and the presence of each observer line fragment.
        #[test]
        fn display_user_defined_only() {
            let mut obs = Observatories::new();

            // Build two user-defined observers (elevation in kilometers).
            obs.create_observer(10.0, 0.0, 0.0, Some("UserA".to_string()));
            obs.create_observer(20.0, 45.0, 2.0, Some("UserB".to_string()));

            let s = format!("{obs}");

            // Header must be present
            assert!(
                s.starts_with("User-defined observers:\n"),
                "Missing 'User-defined observers:' header. Got:\n{s}"
            );

            // Each user observer should be listed with their name and longitude fragment
            assert!(
                s.contains("UserA (lon: 10.000000°"),
                "Missing formatted line for UserA. Got:\n{s}"
            );
            assert!(
                s.contains("UserB (lon: 20.000000°"),
                "Missing formatted line for UserB. Got:\n{s}"
            );

            // The MPC section should not appear if not initialized
            assert!(
                !s.contains("MPC observers:"),
                "Unexpected 'MPC observers:' section when OnceLock is unset. Got:\n{s}"
            );
        }

        /// If the MPC table is initialized, ensure the "MPC observers" section appears.
        ///
        /// Notes
        /// -----
        /// * We set the OnceLock<MpcCodeObs> with a single entry.
        /// * We only check for presence of the section and the MPC code tag.
        #[test]
        fn display_includes_mpc_section_when_set() {
            let mut obs = Observatories::new();

            // One user-defined observer so the first section is non-empty.
            obs.create_observer(0.0, 0.0, 0.0, Some("UserOnly".to_string()));

            // Prepare an MPC observer entry.
            let mpc_site = Observer::new(
                -156.2575,
                20.7075,
                3.055,
                Some("Haleakala".to_string()),
                None,
                None,
            )
            .expect("Failed to create MPC observer");

            // Build an MpcCodeObs map with a single code.
            // If your `MpcCodeObs` is a type alias, this should compile as-is.
            // Example: `pub type MpcCodeObs = std::collections::HashMap<String, Observer>` (or Arc<Observer>).
            let mut mpc_table: MpcCodeObs = Default::default();
            // Adjust Arc<Observer> vs Observer depending on your alias:
            // If it is `HashMap<String, Arc<Observer>>`, wrap with `Arc::new(mpc_site)`.
            use std::sync::Arc;
            mpc_table.insert("I41".to_string(), Arc::new(mpc_site));

            // Initialize the OnceLock (only once)
            obs.mpc_code_obs
                .set(mpc_table)
                .expect("OnceLock<MpcCodeObs> was already initialized");

            let s = format!("{obs}");

            // MPC section header must be present now
            assert!(
                s.contains("MPC observers:"),
                "Missing 'MPC observers:' header after setting OnceLock. Got:\n{s}"
            );

            // The code tag should appear in the MPC line
            assert!(
                s.contains("[I41]"),
                "Missing MPC code tag '[I41]' in output. Got:\n{s}"
            );
        }
    }
}
