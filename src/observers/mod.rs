pub mod bimap;
pub mod observatories;
pub(crate) mod observer_position;

use hifitime::ut1::Ut1Provider;
use hifitime::Epoch;
use nalgebra::{Matrix3, Vector3};
use ordered_float::NotNan;

use crate::constants::{Degree, Kilometer};
use crate::constants::{DPI, ERAU};
use crate::earth_orientation::equequ;
use crate::ref_system::{rotmt, rotpn, RefEpoch, RefSystem};
use crate::time::gmst;

use observer_position::geodetic_to_parallax;

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
pub struct Observer {
    // in degrees east of Greenwich
    pub longitude: ordered_float::NotNan<f64>,
    // phi' is the geocentric latitude and rho is the geocentric distance in earth radii
    // rho cos phi in degrees
    pub rho_cos_phi: ordered_float::NotNan<f64>,
    // rho sin phi in kilometers
    pub rho_sin_phi: ordered_float::NotNan<f64>,
    pub name: Option<String>,
    // Accuracy of the right ascension and declination
    pub ra_accuracy: Option<ordered_float::NotNan<f64>>,
    pub dec_accuracy: Option<ordered_float::NotNan<f64>>,
}

impl Observer {
    /// Create a new observer
    ///
    /// Argument
    /// --------
    /// * `longitude`: observer longitude
    /// * `latitude`: observer latitude
    /// * `elevation`: observer elevation
    /// * `name`: observer name
    ///
    /// Return
    /// ------
    /// * Observer object
    pub(crate) fn new(
        longitude: Degree,
        latitude: Degree,
        elevation: Kilometer,
        name: Option<String>,
        ra_accuracy: Option<ordered_float::NotNan<f64>>,
        dec_accuracy: Option<ordered_float::NotNan<f64>>,
    ) -> Observer {
        let (rho_cos_phi, rho_sin_phi) = geodetic_to_parallax(latitude, elevation);
        Observer {
            longitude: NotNan::try_from(longitude).expect("Longitude cannot be NaN"),
            rho_cos_phi: NotNan::try_from(rho_cos_phi).expect("Longitude cannot be NaN"),
            rho_sin_phi: NotNan::try_from(rho_sin_phi).expect("Longitude cannot be NaN"),
            name,
            ra_accuracy,
            dec_accuracy,
        }
    }

    /// Get the fixed position of an observatory using its geographic coordinates
    ///
    /// Argument
    /// --------
    /// * longitude: observer longitude in Degree
    /// * latitude: observer latitude in Degree
    /// * height: observer height in Degree
    ///
    /// Return
    /// ------
    /// * observer fixed coordinates vector on the Earth (not corrected from Earth motion)
    /// * units is AU
    pub(crate) fn body_fixed_coord(&self) -> Vector3<f64> {
        let lon_radians = self.longitude.to_radians();

        Vector3::new(
            ERAU * self.rho_cos_phi.into_inner() * lon_radians.cos(),
            ERAU * self.rho_cos_phi.into_inner() * lon_radians.sin(),
            ERAU * self.rho_sin_phi.into_inner(),
        )
    }

    /// Compute the observer’s geocentric position and velocity in the ecliptic J2000 frame.
    ///
    /// This function calculates the position and velocity of a ground-based observer relative to the Earth's
    /// center of mass, accounting for Earth rotation (via GMST), nutation, and the observer’s geographic location.
    /// The result is expressed in the ecliptic mean J2000 frame, suitable for use in orbital initial determination.
    ///
    /// Arguments
    /// ---------
    /// * `observer`: a reference to an [`Observer`] containing the site longitude and parallax parameters.
    /// * `tmjd`: observation epoch as a [`hifitime::Epoch`] in TT.
    /// * `ut1_provider`: a reference to a [`hifitime::ut1::Ut1Provider`] for accurate UT1 conversion.
    ///
    /// Returns
    /// --------
    /// * `(dx, dv)` – Tuple of:
    ///     - `dx`: observer geocentric position vector in ecliptic mean J2000 frame [AU].
    ///     - `dv`: observer velocity vector due to Earth's rotation, in the same frame [AU/day].
    ///
    /// Remarks
    /// -------
    /// * Internally, this function:
    ///     1. Computes the body-fixed coordinates of the observer.
    ///     2. Derives its rotational velocity: `v = ω × r`.
    ///     3. Applies Earth orientation corrections using:
    ///         - Greenwich Mean Sidereal Time (GMST),
    ///         - Equation of the equinoxes,
    ///         - Precession and nutation transformation (`rotpn`).
    ///     4. Returns position and velocity in the J2000 ecliptic frame (used in classical orbital mechanics).
    ///
    /// # See also
    /// * [`Observer::body_fixed_coord`] – observer's base vector in Earth-fixed frame
    /// * [`rotpn`] – rotation between reference frames
    /// * [`gmst`], [`equequ`] – time-dependent Earth orientation
    pub(crate) fn pvobs(
        &self,
        tmjd: &Epoch,
        ut1_provider: &Ut1Provider,
    ) -> (Vector3<f64>, Vector3<f64>) {
        // Initialisation
        let omega = Vector3::new(0.0, 0.0, DPI * 1.00273790934);

        // Get the coordinates of the observer on Earth
        let dxbf = self.body_fixed_coord();

        // Get the observer velocity due to Earth rotation
        let dvbf = omega.cross(&dxbf);

        // deviation from Orbfit, use of another conversion from MJD UTC (ET scale) to UT1 scale
        // based on the hifitime crate
        let mjd_ut1 = tmjd.to_ut1(ut1_provider.to_owned());
        let tut = mjd_ut1.to_mjd_tai_days();

        // Compute the Greenwich sideral apparent time
        let gast = gmst(tut) + equequ(tmjd.to_mjd_tt_days());

        // Earth rotation matrix
        let rot = rotmt(-gast, 2);

        // Compute the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let rer_sys1 = RefSystem::Equt(RefEpoch::Epoch(tmjd.to_mjd_tt_days()));
        let rer_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let rot1 = rotpn(&rer_sys1, &rer_sys2);

        let rot1_mat = Matrix3::from(rot1).transpose();
        let rot_mat = Matrix3::from(rot).transpose();

        let rotmat = rot1_mat * rot_mat;

        // Apply transformation to the observer position and velocity
        let dx = rotmat * dxbf;
        let dv = rotmat * dvbf;

        (dx, dv)
    }
}

#[cfg(test)]
mod observer_test {

    use crate::{error_models::ErrorModel, outfit::Outfit};

    use super::*;

    #[test]
    fn test_observer_constructor() {
        let observer = Observer::new(0.0, 0.0, 0.0, None, None, None);
        assert_eq!(observer.longitude, 0.0);
        assert_eq!(observer.rho_cos_phi, 1.0);
        assert_eq!(observer.rho_sin_phi, 0.0);

        let observer = Observer::new(
            289.25058,
            -30.2446,
            2647.,
            Some("Rubin Observatory".to_string()),
            Some(NotNan::new(0.0001).unwrap()),
            Some(NotNan::new(0.0001).unwrap()),
        );

        assert_eq!(observer.longitude, 289.25058);
        assert_eq!(observer.rho_cos_phi, 0.8649760504617418);
        assert_eq!(observer.rho_sin_phi, -0.5009551027512434);
    }

    #[test]
    fn body_fixed_coord_test() {
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, None, None, None);
        let obs_fixed_vector = pan_starrs.body_fixed_coord();
        assert_eq!(
            obs_fixed_vector,
            Vector3::new(
                -0.00003653799439776371,
                -0.00001607260397528885,
                0.000014988110430544328
            )
        )
    }

    #[test]
    fn pvobs_test() {
        let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let tmjd = 57028.479297592596;
        let epoch = Epoch::from_mjd_in_time_scale(tmjd, hifitime::TimeScale::TT);
        // longitude, latitude and height of Pan-STARRS 1, Haleakala
        let (lon, lat, h) = (203.744090000, 20.707233557, 3067.694);
        let pan_starrs = Observer::new(lon, lat, h, Some("Pan-STARRS 1".to_string()), None, None);

        let (observer_position, observer_velocity) =
            &pan_starrs.pvobs(&epoch, state.get_ut1_provider());

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
}
