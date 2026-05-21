//! Precomputed observer-centric positions for every observation epoch.
//!
//! This module computes and caches the **geocentric** and **heliocentric** position
//! (and velocity) of each observer at the precise epoch of each observation. Unlike
//! the body-fixed quantities in [`crate::cache::observer_fixed_cache`], these values
//! depend on the observation time: Earth's rotation and orbital motion must be
//! integrated at each epoch using the JPL ephemeris and the UT1 time scale.
//!
//! # Workflow
//!
//! 1. For each observation, look up the pre-built [`ObserverFixedCache`] by observer ID.
//! 2. Call [`photom::observer::Observer::pvobs`] with the observation epoch and the
//!    body-fixed cache to obtain the **geocentric position and velocity** in the
//!    mean ecliptic J2000 frame.
//! 3. Call [`photom::observer::Observer::helio_position`] with the JPL ephemeris to
//!    add the geocentric Earth position and obtain the **heliocentric position**.
//! 4. Store the results in [`ObserverCentricCache`].
//!
//! # Coordinate system
//!
//! All output vectors are expressed in the **ecliptic mean J2000** reference frame:
//!
//! - positions in **AU**
//! - velocities in **AU/day**
//!
//! # Organisation
//!
//! - [`ObserverGeocentricPosition`] / [`ObserverGeocentricVelocity`] /
//!   [`ObserverHeliocentricPosition`] — NaN-safe 3-vector type aliases.
//! - [`ObserverCentricCache`] — epoch-dependent position/velocity for one observation.
//! - [`CentricObserverCache`] — `Vec` of [`ObserverCentricCache`], indexed by `ObsIndex`.
//! - [`build_centric_observer_cache`] — constructs the full vector for a dataset.

use hifitime::{ut1::Ut1Provider, Epoch};
use nalgebra::Vector3;
use ordered_float::NotNan;
use photom::{observation_dataset::ObsDataset, observer::Observer, MJDTT};

use crate::{
    cache::observer_fixed_cache::{BodyFixedObserverCache, ObserverFixedCache},
    observer_extension::ResolvedObserver,
    JPLEphem, OutfitError,
};

/// Geocentric position of the observer at the epoch of an observation.
///
/// Expressed in the **ecliptic mean J2000** frame, in **AU**.
/// Uses [`NotNan<f64>`] components to enforce NaN-safety at construction time.
pub type ObserverGeocentricPosition = Vector3<NotNan<f64>>;

/// Geocentric velocity of the observer at the epoch of an observation.
///
/// Expressed in the **ecliptic mean J2000** frame, in **AU/day**.
/// Uses [`NotNan<f64>`] components to enforce NaN-safety at construction time.
pub type ObserverGeocentricVelocity = Vector3<NotNan<f64>>;

/// Heliocentric position of the observer at the epoch of an observation.
///
/// This is the sum of the geocentric observer position and the geocentric
/// position of the Earth's centre, both in the **ecliptic mean J2000** frame,
/// in **AU**.
/// Uses [`NotNan<f64>`] components to enforce NaN-safety at construction time.
pub type ObserverHeliocentricPosition = Vector3<NotNan<f64>>;

/// Geocentric and heliocentric observer state for a single observation epoch.
///
/// Built for each observation in the dataset by [`build_centric_observer_cache`].
/// The fields are indexed positionally: the *i*-th element of
/// [`CentricObserverCache`] corresponds to the *i*-th observation in the dataset
/// (i.e., at `ObsIndex` *i*).
///
/// All vectors are in the **ecliptic mean J2000** frame.
#[derive(Debug)]
pub struct ObserverCentricCache {
    /// Geocentric position of the observer at the observation epoch, in AU.
    pub geo_position: ObserverGeocentricPosition,
    /// Geocentric velocity of the observer at the observation epoch, in AU/day.
    pub geo_velocity: ObserverGeocentricVelocity,
    /// Heliocentric position of the observer at the observation epoch, in AU.
    pub helio_position: ObserverHeliocentricPosition,
}

impl ObserverCentricCache {
    /// Computes the geocentric and heliocentric observer state at a given observation epoch.
    ///
    /// # Arguments
    ///
    /// - `jpl` — JPL planetary ephemeris used to obtain the geocentric Earth position.
    /// - `ut1_provider` — UT1 time scale data required to compute Earth's sidereal angle
    ///   at the observation epoch.
    /// - `obs_time` — observation epoch as a Modified Julian Date in the TT time scale
    ///   (see [`MJDTT`]).
    /// - `observer_fixed_cache` — precomputed body-fixed position and velocity of the
    ///   observer (see [`ObserverFixedCache`]).
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if:
    /// - [`photom::observer::Observer::pvobs`] fails (e.g., missing UT1 data for the epoch), or
    /// - [`photom::observer::Observer::helio_position`] fails (e.g., epoch out of range
    ///   for the JPL ephemeris).
    pub fn new(
        jpl: &JPLEphem,
        ut1_provider: &Ut1Provider,
        obs_time: MJDTT,
        observer_fixed_cache: &ObserverFixedCache,
    ) -> Result<Self, OutfitError> {
        let obs_mjd = Epoch::from_mjd_in_time_scale(obs_time, hifitime::TimeScale::TT);
        let (geocentric_pos, geocentric_vel) =
            Observer::pvobs(&obs_mjd, ut1_provider, observer_fixed_cache)?;

        let heliocentric_pos = Observer::helio_position(jpl, &obs_mjd, &geocentric_pos)?;

        Ok(Self {
            geo_position: geocentric_pos,
            geo_velocity: geocentric_vel,
            helio_position: heliocentric_pos,
        })
    }
}

/// Full observer-centric cache for an entire observation dataset.
///
/// A contiguous `Vec` where the element at index *i* holds the precomputed
/// geocentric and heliocentric state for the *i*-th observation
/// (i.e., at [`photom::ObsIndex`] *i*).
pub type CentricObserverCache = Vec<ObserverCentricCache>;

/// Builds the [`CentricObserverCache`] for every observation in a dataset.
///
/// Iterates over all observations in `obs_dataset`, looks up the pre-built
/// body-fixed cache entry for the corresponding observer, and computes the
/// epoch-dependent geocentric and heliocentric state via
/// [`ObserverCentricCache::new`].
///
/// # Arguments
///
/// - `jpl` — JPL planetary ephemeris used to compute Earth's heliocentric position.
/// - `ut1_provider` — UT1 time scale data for Earth rotation at each epoch.
/// - `obs_dataset` — the full observation dataset to process.
/// - `observer_fixed_cache` — pre-built map from [`photom::observer::dataset::ObserverId`]
///   to body-fixed observer state (see [`BodyFixedObserverCache`]).
///
/// # Errors
///
/// Returns [`OutfitError`] if:
/// - any observation has no associated observer ID
///   ([`OutfitError::ObserverIdIsNone`]), or
/// - the observer ID is not found in `observer_fixed_cache`, or
/// - [`ObserverCentricCache::new`] fails for any observation epoch.
///
/// # Examples
///
/// ```rust,ignore
/// let centric_cache = build_centric_observer_cache(
///     &jpl,
///     &ut1_provider,
///     &obs_dataset,
///     &fixed_cache,
/// )?;
/// ```
pub fn build_centric_observer_cache(
    jpl: &JPLEphem,
    ut1_provider: &Ut1Provider,
    obs_dataset: &ObsDataset,
    observer_fixed_cache: &BodyFixedObserverCache,
) -> Result<CentricObserverCache, OutfitError> {
    #[cfg(not(feature = "parallel"))]
    let iter = obs_dataset.iter_observations();

    #[cfg(feature = "parallel")]
    use rayon::iter::{IndexedParallelIterator, ParallelIterator};
    #[cfg(feature = "parallel")]
    let iter = obs_dataset.par_iter_observations();

    iter.enumerate()
        .map(|(idx, obs)| {
            let observer_id = obs
                .observer_id()
                .ok_or_else(|| OutfitError::ObserverIdIsNone(idx as u64))?;

            let fixed_cache = observer_fixed_cache
                .get(observer_id)
                .ok_or_else(|| OutfitError::ObserverIdIsNone(idx as u64))?;

            ObserverCentricCache::new(jpl, ut1_provider, obs.mjd_tt(), fixed_cache)
        })
        .collect()
}

#[cfg(test)]
mod observer_test {

    use photom::{Meters, Radians};

    use crate::{
        cache::observer_centric_cache::ObserverCentricCache,
        test_fixture::{JPL_EPHEM_HORIZON, UT1_PROVIDER},
    };

    use super::*;

    fn to_observer(
        longitude: Radians,
        latitude: Radians,
        height: Meters,
        name: Option<String>,
        ra_accuracy: Option<f64>,
        dec_accuracy: Option<f64>,
    ) -> Observer {
        Observer::new(longitude, latitude, height, name, ra_accuracy, dec_accuracy)
            .expect("Failed to create Observer")
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (
            203.744090000_f64.to_radians(),
            20.707233557_f64.to_radians(),
            3067.694,
        );
        let pan_starrs = to_observer(lon, lat, h, None, None, None);
        assert_eq!(
            pan_starrs
                .earth_fixed_position()
                .unwrap()
                .map(|x| x.into_inner()),
            Vector3::new(
                -0.00003653799439776371,
                -0.00001607260397528885,
                0.000014988110430544328
            )
        );
    }

    #[test]
    fn pvobs_test() {
        let tmjd = 57028.479297592596;
        let epoch = Epoch::from_mjd_in_time_scale(tmjd, hifitime::TimeScale::TT);
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (
            203.744090000_f64.to_radians(),
            20.707233557_f64.to_radians(),
            3067.694,
        );

        let pan_starrs = to_observer(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None);

        let observer_fixed_cache: ObserverFixedCache = (&pan_starrs).try_into().unwrap();

        let (observer_position, observer_velocity) =
            Observer::pvobs(&epoch, &UT1_PROVIDER, &observer_fixed_cache).unwrap();

        assert_eq!(
            observer_position.as_slice(),
            [
                -2.086211182493635e-5,
                3.718476815327979e-5,
                2.4978996447997476e-7
            ]
        );
        assert_eq!(
            observer_velocity.as_slice(),
            [
                -0.0002143246535691577,
                -0.00012059801691431748,
                5.262184624215718e-5
            ]
        );
    }

    #[test]
    fn test_helio_pos_obs() {
        let (lon, lat, h) = (203.744090000_f64, 20.707233557_f64, 3067.694_f64);
        let pan_starrs = to_observer(
            lon.to_radians(),
            lat.to_radians(),
            h,
            Some("Pan-STARRS 1".to_string()),
            None,
            None,
        );
        let observer_fixed_cache: ObserverFixedCache = (&pan_starrs).try_into().unwrap();

        let cases = [
            (
                57_028.479_297_592_596,
                [-0.2645666171464416, 0.8689351643701766, 0.3766996211107864],
            ),
            (
                57_049.245_147_592_59,
                [-0.5891631852137064, 0.7238872516824697, 0.3138186516540669],
            ),
            (
                57_063.977_117_592_59,
                [-0.7743280306286537, 0.5612532665812755, 0.24333415479994636],
            ),
        ];

        for (tmjd, expected) in cases {
            let obs = ObserverCentricCache::new(
                &JPL_EPHEM_HORIZON,
                &UT1_PROVIDER,
                tmjd,
                &observer_fixed_cache,
            )
            .unwrap();

            assert_eq!(obs.helio_position.as_slice(), expected, "tmjd = {tmjd}");
        }
    }
}
