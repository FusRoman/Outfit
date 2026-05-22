//! Ground-observer geometry: body-fixed and heliocentric position routines.
//!
//! This module provides the [`ResolvedObserver`](crate::observer_extension::ResolvedObserver) trait, which extends
//! [`photom::observer::Observer`] with the geometric computations required
//! before any orbit determination step:
//!
//! - [`ResolvedObserver::earth_fixed_position`](crate::observer_extension::ResolvedObserver::earth_fixed_position) — observer position in the
//!   geocentric Earth-fixed frame (AU).
//! - [`ResolvedObserver::earth_fixed_velocity`](crate::observer_extension::ResolvedObserver::earth_fixed_velocity) — observer velocity due to
//!   Earth's sidereal rotation (AU/day).
//! - [`ResolvedObserver::pvobs`](crate::observer_extension::ResolvedObserver::pvobs) — geocentric position and velocity in the
//!   **ecliptic mean J2000** frame at a given epoch, accounting for Earth
//!   rotation, nutation, and precession.
//! - [`ResolvedObserver::helio_position`](crate::observer_extension::ResolvedObserver::helio_position) — heliocentric observer position in
//!   the **equatorial mean J2000** frame (AU), formed by adding the Earth's
//!   JPL-ephemeris position to the geocentric site vector.

use hifitime::{ut1::Ut1Provider, Epoch};
use nalgebra::Vector3;
use ordered_float::NotNan;
use photom::{constants::ERAU, observer::Observer};

use crate::{
    cache::{
        observer_centric_cache::{
            ObserverGeocentricPosition, ObserverGeocentricVelocity, ObserverHeliocentricPosition,
            ObserverHeliocentricVelocity,
        },
        observer_fixed_cache::{ObserverFixedCache, ObserverFixedPosition, ObserverFixedVelocity},
    },
    constants::{EARTH_ROTATION, ROT_ECLMJ2000_TO_EQUMJ2000},
    conversion::ToNotNan,
    earth_orientation::equequ,
    ref_system::{rotmt, rotpn, RefEpoch, RefSystem},
    time::gmst,
    JPLEphem, OutfitError,
};

pub trait ResolvedObserver {
    /// Get the fixed position of an observatory using its geographic coordinates
    ///
    /// Return
    /// ------
    /// * observer fixed coordinates vector on the Earth (not corrected from Earth motion)
    /// * units is AU
    fn earth_fixed_position(&self) -> Result<ObserverFixedPosition, OutfitError>;

    /// Get the fixed velocity of an observatory due to Earth rotation, using its geographic coordinates
    ///
    /// Return
    /// ------
    /// * observer fixed velocity vector due to Earth rotation, in the Earth-fixed frame (not corrected from Earth motion)
    /// * units is AU/day
    fn earth_fixed_velocity(&self) -> Result<ObserverFixedVelocity, OutfitError>;

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
    /// * `compute_velocity`: whether to compute and return the observer's velocity due to Earth's rotation (true) or return a zero vector (false).
    ///
    /// Returns
    /// --------
    /// * `(dx, dv)` – Tuple of:
    ///     - `dx`: observer geocentric position vector in ecliptic mean J2000 frame \[AU\].
    ///     - `dv`: observer velocity vector due to Earth's rotation, in the same frame \[AU/day\].
    ///
    /// Remarks
    /// -------
    /// * Internally, this function:
    ///     1. get the body-fixed coordinates of the observer.
    ///     2. get its rotational velocity: `v = ω × r`.
    ///     3. Applies Earth orientation corrections using:
    ///         - Greenwich Mean Sidereal Time (GMST),
    ///         - Equation of the equinoxes,
    ///         - Precession and nutation transformation (`rotpn`).
    ///     4. Returns position and velocity in the J2000 ecliptic frame (used in classical orbital mechanics).
    ///
    /// # See also
    /// * [`ResolvedObserver::earth_fixed_position`] – observer's base vector in Earth-fixed frame
    /// * [`rotpn`] – rotation between reference frames
    /// * [`gmst`], [`equequ`] – time-dependent Earth orientation
    fn pvobs(
        tmjd: &Epoch,
        ut1_provider: &Ut1Provider,
        observer_fixed_vectors: &ObserverFixedCache,
        compute_velocity: bool,
    ) -> Result<(ObserverGeocentricPosition, ObserverGeocentricVelocity), OutfitError>;

    /// Compute the observer’s heliocentric position in the **equatorial mean J2000** frame.
    ///
    /// This method forms the full heliocentric position of the observing site by combining:
    /// - the site **geocentric** position vector at `epoch`, and
    /// - the Earth’s **heliocentric** position from the JPL ephemerides.
    ///
    /// The input geocentric vector is assumed to be expressed in the **ecliptic mean J2000** frame
    /// (AU). It is rotated to **equatorial mean J2000**, then added to Earth’s heliocentric
    /// position (also in equatorial mean J2000).
    ///
    /// Arguments
    /// -----------------
    /// * `jpl` – [`JPLEphem`] providing Earth's heliocentric state.
    /// * `epoch` – Observation epoch in the **TT** time scale.
    /// * `observer_geocentric_position` – Geocentric site position **in ecliptic mean J2000** (AU).
    ///
    /// Return
    /// ----------
    /// * `Result<ObserverHeliocentricPosition, OutfitError>` – Observer’s **heliocentric** position at `epoch`,
    ///   in **AU**, expressed in **equatorial mean J2000**.
    ///
    /// Remarks
    /// -------------
    /// * If your geocentric site vector is already in **equatorial** J2000, rotate it to
    ///   **ecliptic** before calling this method, or adapt the rotation accordingly.
    /// * This routine is typically used internally when constructing per-observation geometry
    ///   (e.g., within `Observation::new`), ensuring consistent frames for Gauss IOD.
    ///
    /// See also
    /// ------------
    /// * [`ResolvedObserver::pvobs`] – Geocentric position (and velocity) of the site at `epoch`.
    /// * [`JPLEphem`] – Access Earth's heliocentric state from JPL ephemerides.
    /// * [`crate::constants::ROT_ECLMJ2000_TO_EQUMJ2000`] – Rotation between ecliptic and equatorial J2000.
    fn helio_position(
        jpl: &JPLEphem,
        epoch: &Epoch,
        observer_geocentric_position: &ObserverGeocentricPosition,
    ) -> Result<ObserverHeliocentricPosition, OutfitError>;

    /// Compute the observer’s heliocentric velocity in the **equatorial mean J2000** frame.
    ///
    /// This method forms the full heliocentric velocity of the observing site by combining:
    /// - the site **geocentric** velocity vector at `epoch`, and
    /// - the Earth’s **heliocentric** velocity from the JPL ephemerides.
    ///
    /// # Arguments
    ///
    /// * `jpl` – [`JPLEphem`] providing Earth's heliocentric state.
    /// * `epoch` – Observation epoch in the **TT** time scale.
    /// * `observer_geocentric_velocity` – Geocentric site velocity **in ecliptic mean J2000** (AU/day).
    ///
    /// Return
    ///
    /// * `Result<ObserverHeliocentricVelocity, OutfitError>` – Observer’s **heliocentric** velocity at `epoch` in **AU/day**, expressed in **equatorial mean J2000**.
    fn helio_velocity(
        jpl: &JPLEphem,
        epoch: &Epoch,
        observer_geocentric_velocity: &ObserverGeocentricVelocity,
    ) -> Result<ObserverHeliocentricVelocity, OutfitError>;
}

impl ResolvedObserver for Observer {
    fn earth_fixed_position(&self) -> Result<ObserverFixedPosition, OutfitError> {
        let (sin_lon, cos_lon): (NotNan<f64>, NotNan<f64>) = {
            let (s, c) = self.longitude.sin_cos();
            (s.to_notnan()?, c.to_notnan()?)
        };
        let erau_not_nan = ERAU.to_notnan()?;

        Ok(Vector3::new(
            erau_not_nan * self.rho_cos_phi * cos_lon,
            erau_not_nan * self.rho_cos_phi * sin_lon,
            erau_not_nan * self.rho_sin_phi,
        ))
    }

    #[inline]
    fn earth_fixed_velocity(&self) -> Result<ObserverFixedVelocity, OutfitError> {
        Ok(EARTH_ROTATION
            .to_notnan()?
            .cross(&self.earth_fixed_position()?))
    }

    fn pvobs(
        tmjd: &Epoch,
        ut1_provider: &Ut1Provider,
        observer_fixed_vectors: &ObserverFixedCache,
        compute_velocity: bool,
    ) -> Result<(ObserverGeocentricPosition, ObserverGeocentricVelocity), OutfitError> {
        // Get observer position and velocity in the Earth-fixed frame
        let dxbf = observer_fixed_vectors.position();

        // deviation from Orbfit, use of another conversion from MJD UTC (ET scale) to UT1 scale
        // based on the hifitime crate
        let mjd_ut1 = tmjd.to_ut1(ut1_provider);
        let tut = mjd_ut1.to_mjd_tai_days();

        // Compute the Greenwich sideral apparent time
        let gast = gmst(tut) + equequ(tmjd.to_mjd_tt_days());

        // Earth rotation matrix
        let rot = rotmt(-gast, 2);

        // Compute the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let rer_sys1 = RefSystem::Equt(RefEpoch::Epoch(tmjd.to_mjd_tt_days()));
        let rer_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let rot1 = rotpn(&rer_sys1, &rer_sys2)?;

        let rot1_mat = rot1.transpose().to_notnan()?;
        let rot_mat = rot.transpose().to_notnan()?;

        let rotmat = rot1_mat * rot_mat;

        // Apply transformation to the observer position and velocity
        let dx = rotmat * dxbf;

        let dv = if compute_velocity {
            let dvbf = observer_fixed_vectors.velocity();
            rotmat * dvbf
        } else {
            Vector3::zeros()
        };

        Ok((dx, dv))
    }

    fn helio_position(
        jpl: &JPLEphem,
        epoch: &Epoch,
        observer_geocentric_position: &ObserverGeocentricPosition,
    ) -> Result<ObserverHeliocentricPosition, OutfitError> {
        // Earth's heliocentric position
        let earth_pos = jpl.earth_ephemeris(epoch, false).0.to_notnan()?;

        // Transform observer position from ecliptic to equatorial J2000
        let rot_matrix = ROT_ECLMJ2000_TO_EQUMJ2000.to_notnan()?.transpose();

        let helio_pos = earth_pos + rot_matrix * observer_geocentric_position;

        Ok(helio_pos)
    }

    fn helio_velocity(
        jpl: &JPLEphem,
        epoch: &Epoch,
        observer_geocentric_velocity: &ObserverGeocentricVelocity,
    ) -> Result<ObserverHeliocentricVelocity, OutfitError> {
        // Earth's heliocentric velocity — already in ecliptic J2000, AU/day
        let earth_vel = jpl
            .earth_ephemeris(epoch, true)
            .1
            .expect("Velocity is always available, this should not happen")
            .to_notnan()?;

        // geo_velocity is in equatorial J2000 → rotate to ecliptic (same as helio_position)
        let rot_matrix = ROT_ECLMJ2000_TO_EQUMJ2000.to_notnan()?.transpose();
        let helio_vel = earth_vel + rot_matrix * observer_geocentric_velocity;
        Ok(helio_vel)
    }
}
