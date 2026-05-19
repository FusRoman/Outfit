//! Precomputed body-fixed observer positions and velocities.
//!
//! This module computes and caches the **Earth-fixed (body-fixed) position and
//! velocity** of each observer present in an observation dataset. These quantities
//! are time-independent — they depend only on the observer's geographic coordinates
//! (longitude, latitude, height) — so they are computed once at cache-build time
//! and reused for every observation associated with the same observer.
//!
//! # Coordinate system
//!
//! All vectors are expressed in the **geocentric Earth-fixed frame** (ECEF-like),
//! with components given in **astronomical units (AU)** for positions and
//! **AU/day** for velocities.
//!
//! The body-fixed velocity is derived from Earth's sidereal rotation:
//!
//! ```text
//! v_fixed = ω_earth × r_fixed
//! ```
//!
//! where `ω_earth` is the Earth rotation vector (see [`crate::constants::EARTH_ROTATION`]).
//!
//! # Organisation
//!
//! - [`ObserverFixedPosition`] / [`ObserverFixedVelocity`] — type aliases for 3-vectors.
//! - [`ObserverFixedCache`] — holds the fixed position and velocity for one observer.
//! - [`BodyFixedObserverCache`] — map from [`ObserverId`] to [`ObserverFixedCache`].
//! - [`build_fixed_observer_cache`] — constructs the map from an iterator of observers.

use ahash::AHashMap;
use nalgebra::Vector3;
use ordered_float::NotNan;
use photom::observer::{dataset::ObserverId, Observer};

use crate::{
    constants::EARTH_ROTATION, conversion::ToNotNan, observer_extension::ResolvedObserver,
    OutfitError,
};

/// Precomputed **body-fixed** position of the observer in **AU**.
///
/// The vector is expressed in the geocentric Earth-fixed frame. Its components
/// are stored as [`NotNan<f64>`] to guarantee the absence of NaN values at
/// construction time.
pub type ObserverFixedPosition = Vector3<NotNan<f64>>;

/// Precomputed **body-fixed** velocity of the observer in **AU/day**.
///
/// Derived from the cross product of Earth's rotation vector with the observer's
/// body-fixed position: `v = ω × r`. Stored as [`NotNan<f64>`] for the same
/// NaN-safety guarantee as [`ObserverFixedPosition`].
pub type ObserverFixedVelocity = Vector3<NotNan<f64>>;

/// Body-fixed position and velocity for a single ground-based observer.
///
/// This cache entry is time-independent: it is built once from the observer's
/// geographic coordinates and reused across all observations made by that observer.
///
/// # Fields
///
/// Both fields are in the geocentric Earth-fixed frame:
///
/// - position in **AU**
/// - velocity in **AU/day** (from Earth rotation)
#[derive(Debug)]
pub struct ObserverFixedCache {
    /// Geocentric Earth-fixed position of the observer, in AU.
    observer_fixed_positions: ObserverFixedPosition,
    /// Geocentric Earth-fixed velocity of the observer due to Earth's rotation, in AU/day.
    observer_fixed_velocities: ObserverFixedVelocity,
}

impl ObserverFixedCache {
    /// Constructs a new [`ObserverFixedCache`] from a ground-based [`Observer`].
    ///
    /// Computes the body-fixed position from the observer's geodetic coordinates
    /// (longitude, latitude, height above ellipsoid), then derives the body-fixed
    /// velocity using Earth's sidereal rotation vector:
    ///
    /// ```text
    /// v_fixed = ω_earth × r_fixed
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if:
    /// - the observer's Earth-fixed position cannot be computed
    ///   (e.g., invalid geodetic coordinates), or
    /// - a NaN is encountered when converting to [`NotNan<f64>`].
    pub fn new(observer: &Observer) -> Result<Self, OutfitError> {
        // Body-fixed position in AU from (ρ·cosφ, ρ·sinφ) scaled by Earth radius (AU).
        let body_fixed_pos = observer.earth_fixed_position()?;

        // Body-fixed velocity from Earth rotation.
        let body_fixed_vel: Vector3<NotNan<f64>> =
            EARTH_ROTATION.to_notnan()?.cross(&body_fixed_pos);

        Ok(Self {
            observer_fixed_positions: body_fixed_pos,
            observer_fixed_velocities: body_fixed_vel,
        })
    }

    /// Returns the precomputed body-fixed position of the observer, in AU.
    pub fn position(&self) -> &ObserverFixedPosition {
        &self.observer_fixed_positions
    }

    /// Returns the precomputed body-fixed velocity of the observer, in AU/day.
    ///
    /// This is the velocity due to Earth's sidereal rotation: `ω × r`.
    pub fn velocity(&self) -> &ObserverFixedVelocity {
        &self.observer_fixed_velocities
    }
}

impl TryFrom<&Observer> for ObserverFixedCache {
    type Error = OutfitError;

    /// Converts an [`Observer`] reference into an [`ObserverFixedCache`].
    ///
    /// Delegates to [`ObserverFixedCache::new`].
    ///
    /// # Errors
    ///
    /// Propagates any error from [`ObserverFixedCache::new`].
    fn try_from(resolved: &Observer) -> Result<Self, Self::Error> {
        Self::new(resolved)
    }
}

/// Cache mapping observer IDs to their precomputed body-fixed positions and velocities.
///
/// This hash map is built once before any trajectory fitting (see
/// [`build_fixed_observer_cache`]) and looked up by [`ObserverId`] for every
/// observation in the dataset. Using [`AHashMap`] provides fast, non-cryptographic
/// hashing suited for integer-keyed lookups.
pub type BodyFixedObserverCache = AHashMap<ObserverId, ObserverFixedCache>;

/// Builds the [`BodyFixedObserverCache`] from an iterator of `(ObserverId, &Observer)` pairs.
///
/// For each observer, computes the body-fixed position and velocity and stores
/// them in the map. This function is typically called once before the main
/// cache-building step in [`crate::cache::OutfitCache::build`].
///
/// # Arguments
///
/// - `observers` — an iterator yielding `(ObserverId, &Observer)` pairs, typically
///   obtained from [`photom::observation_dataset::ObsDataset::iter_observer`].
///
/// # Errors
///
/// Returns [`OutfitError`] if [`ObserverFixedCache::new`] fails for any observer
/// in the iterator (e.g., invalid geodetic coordinates or NaN conversion).
///
/// # Examples
///
/// ```rust,ignore
/// let cache = build_fixed_observer_cache(obs_dataset.iter_observer()?)?;
/// ```
pub fn build_fixed_observer_cache<'a>(
    observers: impl Iterator<Item = (ObserverId, &'a Observer)>,
) -> Result<BodyFixedObserverCache, OutfitError> {
    observers
        .map(|(id, obs)| -> Result<_, OutfitError> { Ok((id, obs.try_into()?)) })
        .collect::<Result<BodyFixedObserverCache, _>>()
}
