pub mod observer_centric_cache;
pub mod observer_fixed_cache;

use hifitime::ut1::Ut1Provider;
use photom::{observation_dataset::ObsDataset, observer::dataset::ObserverId, ObsIndex};

use crate::{
    cache::{
        observer_centric_cache::{
            build_centric_observer_cache, CentricObserverCache, ObserverCentricCache,
            ObserverGeocentricPosition, ObserverGeocentricVelocity, ObserverHeliocentricPosition,
        },
        observer_fixed_cache::{
            build_fixed_observer_cache, BodyFixedObserverCache, ObserverFixedCache,
        },
    },
    JPLEphem, OutfitError,
};

/// Precomputed observer positions for all observations in a dataset.
///
/// Built once before any trajectory fitting. Each entry is keyed by the
/// observation's [`ObsIndex`], which is stable for the lifetime of the
/// [`ObsDataset`].
pub struct OutfitCache {
    /// len == number of observations in the dataset. Indexed by [`ObsIndex`].
    observer_centric: CentricObserverCache,
    /// len == number of observer in the dataset. Indexed by observer ID.
    observer_fixed: BodyFixedObserverCache,
}

impl OutfitCache {
    pub fn get_centric(&self, idx: ObsIndex) -> &ObserverCentricCache {
        &self.observer_centric[idx]
    }

    pub fn get_fixed(&self, observer_id: ObserverId) -> Option<&ObserverFixedCache> {
        self.observer_fixed.get(&observer_id)
    }

    /// Build the cache for every observation in `obs_dataset`.
    pub fn build(
        obs_dataset: &ObsDataset,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
    ) -> Result<Self, OutfitError> {
        let observer_iter = obs_dataset.iter_observer()?;

        let observer_fixed_cache = build_fixed_observer_cache(observer_iter)?;

        let observer_centric_cache =
            build_centric_observer_cache(jpl, ut1_provider, obs_dataset, &observer_fixed_cache)?;

        Ok(Self {
            observer_centric: observer_centric_cache,
            observer_fixed: observer_fixed_cache,
        })
    }

    pub fn get_observer_fixed_cache(&self, observer_id: ObserverId) -> Option<&ObserverFixedCache> {
        self.observer_fixed.get(&observer_id)
    }

    /// Accessor for the precomputed geocentric position of an observer.
    pub fn get_observer_geocentric_position(&self, idx: ObsIndex) -> &ObserverGeocentricPosition {
        &self.get_centric(idx).geo_position
    }

    pub fn get_observer_geocentric_velocity(&self, idx: ObsIndex) -> &ObserverGeocentricVelocity {
        &self.get_centric(idx).geo_velocity
    }

    /// Accessor for the precomputed heliocentric position of an observer.
    pub fn get_helio_position(&self, idx: ObsIndex) -> &ObserverHeliocentricPosition {
        &self.get_centric(idx).helio_position
    }
}
