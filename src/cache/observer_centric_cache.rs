use hifitime::{ut1::Ut1Provider, Epoch};
use nalgebra::Vector3;
use ordered_float::NotNan;
use photom::{observation_dataset::ObsDataset, observer::Observer, MJDTT};

use crate::{
    cache::observer_fixed_cache::{BodyFixedObserverCache, ObserverFixedCache},
    observer_extension::ResolvedObserver,
    JPLEphem, OutfitError,
};

/// Geocentric position of the observer at `time` of an observation (AU, ecliptic mean J2000).
pub type ObserverGeocentricPosition = Vector3<NotNan<f64>>;
/// Geocentric velocity of the observer at `time` of an observation (AU/day, ecliptic mean J2000).
pub type ObserverGeocentricVelocity = Vector3<NotNan<f64>>;
/// Heliocentric position of the observer at `time` of an observation (AU, ecliptic mean J2000).
pub type ObserverHeliocentricPosition = Vector3<NotNan<f64>>;

/// Geocentric and heliocentric observer positions for a single observation epoch.
#[derive(Debug)]
pub struct ObserverCentricCache {
    pub geo_position: ObserverGeocentricPosition,
    pub geo_velocity: ObserverGeocentricVelocity,
    pub helio_position: ObserverHeliocentricPosition,
}

impl ObserverCentricCache {
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

pub type CentricObserverCache = Vec<ObserverCentricCache>;

pub fn build_centric_observer_cache(
    jpl: &JPLEphem,
    ut1_provider: &Ut1Provider,
    obs_dataset: &ObsDataset,
    observer_fixed_cache: &BodyFixedObserverCache,
) -> Result<CentricObserverCache, OutfitError> {
    obs_dataset
        .iter_observations()
        .enumerate()
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
