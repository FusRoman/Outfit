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

/// Heliocentric velocity of the observer at the epoch of an observation.
///
/// This is the sum of the geocentric observer velocity and the geocentric
/// velocity of the Earth's centre, both in the **ecliptic mean J2000** frame,
/// in **AU/day**.
/// Uses [`NotNan<f64>`] components to enforce NaN-safety at construction time.
pub type ObserverHeliocentricVelocity = Vector3<NotNan<f64>>;

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
    pub geo_velocity: Option<ObserverGeocentricVelocity>,
    /// Heliocentric position of the observer at the observation epoch, in AU.
    pub helio_position: ObserverHeliocentricPosition,
    /// Heliocentric velocity of the observer at the observation epoch, in AU/day.
    pub helio_velocity: Option<ObserverHeliocentricVelocity>,
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
    /// - `cache_velocity` — whether to compute and cache the velocity components. If `false`, the velocity fields will be set to `None` to save computation time and memory.
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
        cache_velocity: bool,
    ) -> Result<Self, OutfitError> {
        let obs_mjd = Epoch::from_mjd_in_time_scale(obs_time, hifitime::TimeScale::TT);
        let (geocentric_pos, geocentric_vel) =
            Observer::pvobs(&obs_mjd, ut1_provider, observer_fixed_cache, cache_velocity)?;

        let heliocentric_pos = Observer::helio_position(jpl, &obs_mjd, &geocentric_pos)?;
        let heliocentric_vel = if cache_velocity {
            Some(Observer::helio_velocity(jpl, &obs_mjd, &geocentric_vel)?)
        } else {
            None
        };

        Ok(Self {
            geo_position: geocentric_pos,
            geo_velocity: if cache_velocity {
                Some(geocentric_vel)
            } else {
                None
            },
            helio_position: heliocentric_pos,
            helio_velocity: heliocentric_vel,
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
/// - `cache_velocity` — whether to compute and cache the velocity components in the resulting `ObserverCentricCache`. If `false`, the velocity fields will be set to `None` to save computation time and memory.
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
    cache_velocity: bool,
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

            ObserverCentricCache::new(jpl, ut1_provider, obs.mjd_tt(), fixed_cache, cache_velocity)
        })
        .collect()
}

#[cfg(test)]
mod observer_test {

    use approx::assert_relative_eq;
    use photom::{Meters, Radians};

    use crate::{
        cache::observer_centric_cache::ObserverCentricCache,
        conversion::ToNotNan,
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
            Observer::pvobs(&epoch, &UT1_PROVIDER, &observer_fixed_cache, true).unwrap();

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

        let (_, observer_velocity) =
            Observer::pvobs(&epoch, &UT1_PROVIDER, &observer_fixed_cache, false).unwrap();

        assert_eq!(observer_velocity.as_slice(), [0.0, 0.0, 0.0]);
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
                true,
            )
            .unwrap();

            assert_eq!(obs.helio_position.as_slice(), expected, "tmjd = {tmjd}");
        }
    }

    fn v3(x: f64, y: f64, z: f64) -> Vector3<NotNan<f64>> {
        Vector3::new(x, y, z).to_notnan().unwrap()
    }

    fn to_f64(v: &Vector3<NotNan<f64>>) -> Vector3<f64> {
        v.map(|x| x.into_inner())
    }

    fn assert_v3_eq(actual: &Vector3<NotNan<f64>>, expected: &Vector3<NotNan<f64>>, eps: f64) {
        assert_relative_eq!(to_f64(actual), to_f64(expected), epsilon = eps);
    }

    #[test]
    fn test_helio_pos_vel_geocenter() {
        let geocenter =
            Observer::from_parallax(0.0, 0.0, 0.0, Some("Geocenter".to_string()), None, None)
                .unwrap();

        let observer_fixed_cache: ObserverFixedCache = (&geocenter).try_into().unwrap();

        let obs = ObserverCentricCache::new(
            &JPL_EPHEM_HORIZON,
            &UT1_PROVIDER,
            59000.0,
            &observer_fixed_cache,
            true,
        )
        .unwrap();

        assert_v3_eq(&obs.geo_position, &v3(0.0, 0.0, 0.0), 1e-10);
        assert_v3_eq(
            &obs.helio_position,
            &v3(
                -0.35112872984703947,
                -0.8726911829575209,
                -0.37831199013326505,
            ),
            1e-10,
        );
        assert_v3_eq(
            obs.helio_velocity.as_ref().unwrap(),
            &v3(
                0.015860197805364396,
                -0.005519387867661577,
                -0.002392757495907968,
            ),
            1e-10,
        );

        assert_v3_eq(&obs.geo_position, &v3(0.0, 0.0, 0.0), 1e-10);
        assert_v3_eq(
            obs.geo_velocity.as_ref().unwrap(),
            &v3(0.0, 0.0, 0.0),
            1e-10,
        );
    }

    #[test]
    fn test_helio_pos_vel_mauna_kea() {
        // MPC code 568 — Mauna Kea (W. M. Keck Observatory)
        // lon = 204.5284°, lat = 19.8260°, h = 4160 m
        let mauna_kea = Observer::from_parallax(
            204.5278_f64.to_radians(),
            0.94171,
            0.33725,
            Some("Maunakea".to_string()),
            None,
            None,
        )
        .unwrap();

        let observer_fixed_cache: ObserverFixedCache = (&mauna_kea).try_into().unwrap();

        let obs = ObserverCentricCache::new(
            &JPL_EPHEM_HORIZON,
            &UT1_PROVIDER,
            59000.0,
            &observer_fixed_cache,
            true,
        )
        .unwrap();

        assert_v3_eq(
            &obs.helio_position,
            &v3(
                -3.511307549159519e-01,
                -8.726510855672746e-01,
                -3.782976072020051e-01,
            ),
            1e-9,
        );
        assert_v3_eq(
            obs.helio_velocity.as_ref().unwrap(),
            &v3(
                1.560756863717671e-02,
                -5.532323168433832e-03,
                -2.392265222947331e-03,
            ),
            1e-8,
        );
        assert_v3_eq(
            &obs.geo_position,
            &v3(
                -2.025068912418855e-06,
                4.250983777758508e-05,
                -2.753744421400818e-06,
            ),
            1e-8,
        );
        assert_v3_eq(
            obs.geo_velocity.as_ref().unwrap(),
            &v3(
                -2.526291681876792e-04,
                -1.167209148779009e-05,
                5.597018763337186e-06,
            ),
            1e-8,
        );
    }
}
