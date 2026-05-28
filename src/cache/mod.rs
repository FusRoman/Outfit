//! Observer position cache for trajectory fitting.
//!
//! This module is the public entry point for the two-level observer cache used by
//! Outfit before any trajectory fitting or residual computation. The cache is built
//! **once** from a [`photom::observation_dataset::ObsDataset`] and then queried by index or observer ID
//! throughout the fitting pipeline, avoiding repeated ephemeris lookups.
//!
//! # Two-level design
//!
//! The cache is split into two complementary layers:
//!
//! ## 1. Body-fixed layer — [`observer_fixed_cache`](crate::cache::observer_fixed_cache)
//!
//! Stores the **time-independent** Earth-fixed (ECEF-like) position and velocity
//! of every observer. Indexed by `ObserverId`. Built first because the centric
//! layer depends on it.
//!
//! ## 2. Observer-centric layer — [`observer_centric_cache`](crate::cache::observer_centric_cache)
//!
//! Stores the **epoch-dependent** geocentric and heliocentric position (and
//! geocentric velocity) of the observer at the precise time of each observation.
//! Indexed by `ObsIndex` (i.e., by observation order in the dataset).
//!
//! # Build order
//!
//! ```text
//! ObsDataset
//!   │
//!   ├─ iter_observer() ──► build_fixed_observer_cache()  ──► BodyFixedObserverCache
//!   │                                                               │
//!   └─ iter_observations() + BodyFixedObserverCache ──► build_centric_observer_cache()
//!                                                               │
//!                                                        CentricObserverCache
//! ```
//!
//! Both caches are then wrapped in [`OutfitCache`](crate::cache::OutfitCache) and exposed through typed accessors.
//!
//! # Usage
//!
//! ```rust,ignore
//! let cache = OutfitCache::build(&obs_dataset, &jpl, &ut1_provider)?;
//!
//! // Per-observation accessors (indexed by ObsIndex):
//! let geo_pos  = cache.get_observer_geocentric_position(idx);
//! let geo_vel  = cache.get_observer_geocentric_velocity(idx);
//! let helio    = cache.get_helio_position(idx);
//!
//! // Per-observer accessor (indexed by ObserverId):
//! let fixed    = cache.get_observer_fixed_cache(observer_id);
//! ```

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
/// [`OutfitCache`] is the top-level cache built **once** before any trajectory
/// fitting, from a fully loaded [`ObsDataset`]. It encapsulates:
///
/// - a per-observation centric cache ([`CentricObserverCache`]) holding geocentric
///   and heliocentric positions at each observation epoch;
/// - a per-observer body-fixed cache ([`BodyFixedObserverCache`]) holding the
///   time-independent Earth-fixed state of each observer.
///
/// All positions are in the **ecliptic mean J2000** frame (AU / AU·day⁻¹).
///
/// # Build
///
/// Use [`OutfitCache::build`] to construct the cache. Accessor methods then
/// provide O(1) lookups keyed by [`ObsIndex`] or [`ObserverId`].
#[derive(Debug)]
pub struct OutfitCache {
    /// Per-observation centric cache. Length equals the number of observations
    /// in the dataset. Indexed by [`ObsIndex`].
    observer_centric: CentricObserverCache,
    /// Per-observer body-fixed cache. Indexed by [`ObserverId`].
    observer_fixed: BodyFixedObserverCache,
}

impl OutfitCache {
    /// Returns the full [`ObserverCentricCache`] entry for the observation at `idx`.
    ///
    /// The returned reference gives access to the geocentric position, geocentric
    /// velocity, and heliocentric position at the epoch of that observation.
    ///
    /// # Panics
    ///
    /// Panics if `idx` is out of bounds (i.e., ≥ number of observations in the
    /// dataset used to build this cache).
    pub fn get_centric(&self, idx: ObsIndex) -> &ObserverCentricCache {
        &self.observer_centric[idx]
    }

    /// Returns the [`ObserverFixedCache`] for the given observer, if present.
    ///
    /// Returns `None` if `observer_id` is not found in the body-fixed cache
    /// (e.g., the observer was not part of the dataset used to build this cache).
    pub fn get_fixed(&self, observer_id: ObserverId) -> Option<&ObserverFixedCache> {
        self.observer_fixed.get(&observer_id)
    }

    /// Builds the cache for every observation in `obs_dataset`.
    ///
    /// This is the primary constructor. It proceeds in two steps:
    ///
    /// 1. Build the [`BodyFixedObserverCache`] from the dataset's observer list.
    /// 2. Build the [`CentricObserverCache`] by computing epoch-dependent positions
    ///    for each observation using the JPL ephemeris and UT1 provider.
    ///
    /// # Arguments
    ///
    /// - `obs_dataset` — the observation dataset to cache; all observers must have
    ///   valid geodetic coordinates and associated observer IDs.
    /// - `jpl` — JPL planetary ephemeris (DE440 or compatible) for heliocentric
    ///   Earth positions.
    /// - `ut1_provider` — UT1 time scale data for computing Earth's sidereal angle
    ///   at each observation epoch.
    /// - `cache_velocity` — whether to compute and cache the velocity components in the resulting `ObserverCentricCache`. If `false`, the velocity fields will be set to `None` to save computation time and memory.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if:
    /// - the dataset's observer iterator fails, or
    /// - any observer's body-fixed position cannot be computed, or
    /// - any observation lacks an associated observer ID, or
    /// - the JPL ephemeris or UT1 provider cannot be evaluated at a given epoch.
    pub fn build(
        obs_dataset: &ObsDataset,
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        cache_velocity: bool,
    ) -> Result<Self, OutfitError> {
        let observer_iter = obs_dataset.iter_observer()?;

        let observer_fixed_cache = build_fixed_observer_cache(observer_iter)?;

        let observer_centric_cache = build_centric_observer_cache(
            jpl,
            ut1_provider,
            obs_dataset,
            &observer_fixed_cache,
            cache_velocity,
        )?;

        Ok(Self {
            observer_centric: observer_centric_cache,
            observer_fixed: observer_fixed_cache,
        })
    }

    /// Returns the [`ObserverFixedCache`] for the given observer, if present.
    ///
    /// Alias for [`OutfitCache::get_fixed`], provided for explicitness when the
    /// caller already has an [`ObserverId`].
    pub fn get_observer_fixed_cache(&self, observer_id: ObserverId) -> Option<&ObserverFixedCache> {
        self.observer_fixed.get(&observer_id)
    }

    /// Returns the precomputed geocentric position of the observer at observation `idx`.
    ///
    /// The position is in the **ecliptic mean J2000** frame, in **AU**.
    ///
    /// # Panics
    ///
    /// Panics if `idx` is out of bounds.
    pub fn get_observer_geocentric_position(&self, idx: ObsIndex) -> &ObserverGeocentricPosition {
        &self.get_centric(idx).geo_position
    }

    /// Returns the precomputed geocentric velocity of the observer at observation `idx`.
    ///
    /// The velocity is in the **ecliptic mean J2000** frame, in **AU/day**.
    ///
    /// # Panics
    ///
    /// Panics if `idx` is out of bounds.
    pub fn get_observer_geocentric_velocity(
        &self,
        idx: ObsIndex,
    ) -> &Option<ObserverGeocentricVelocity> {
        &self.get_centric(idx).geo_velocity
    }

    /// Returns the precomputed heliocentric position of the observer at observation `idx`.
    ///
    /// The position is in the **ecliptic mean J2000** frame, in **AU**.
    ///
    /// # Panics
    ///
    /// Panics if `idx` is out of bounds.
    pub fn get_helio_position(&self, idx: ObsIndex) -> &ObserverHeliocentricPosition {
        &self.get_centric(idx).helio_position
    }
}
