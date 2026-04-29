use ahash::AHashMap;
use nalgebra::Vector3;
use ordered_float::NotNan;
use photom::observer::{dataset::ObserverId, Observer};

use crate::{
    constants::EARTH_ROTATION, conversion::ToNotNan, observer_extension::ResolvedObserver,
    OutfitError,
};

/// Precomputed **body-fixed** position of the observer in **AU**.
pub type ObserverFixedPosition = Vector3<NotNan<f64>>;
/// Precomputed **body-fixed** velocity of the observer in **AU/day**.
pub type ObserverFixedVelocity = Vector3<NotNan<f64>>;

pub struct ObserverFixedCache {
    observer_fixed_positions: ObserverFixedPosition,
    observer_fixed_velocities: ObserverFixedVelocity,
}

impl ObserverFixedCache {
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

    pub fn position(&self) -> &ObserverFixedPosition {
        &self.observer_fixed_positions
    }

    pub fn velocity(&self) -> &ObserverFixedVelocity {
        &self.observer_fixed_velocities
    }
}

impl TryFrom<&Observer> for ObserverFixedCache {
    type Error = OutfitError;

    fn try_from(resolved: &Observer) -> Result<Self, Self::Error> {
        Self::new(resolved)
    }
}

/// Cache mapping observer IDs to their precomputed body-fixed positions and velocities.
pub type BodyFixedObserverCache = AHashMap<ObserverId, ObserverFixedCache>;

pub fn build_fixed_observer_cache<'a>(
    observers: impl Iterator<Item = (ObserverId, &'a Observer)>,
) -> Result<BodyFixedObserverCache, OutfitError> {
    observers
        .map(|(id, obs)| -> Result<_, OutfitError> { Ok((id, obs.try_into()?)) })
        .collect::<Result<BodyFixedObserverCache, _>>()
}
