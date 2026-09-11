//! Planetary ephemeris backend backed by a NAIF SPK kernel loaded through the
//! ANISE toolkit (feature `ephem-anise`).
//!
//! # What this backend provides
//!
//! A single type, [`AniseEphem`], that answers the same two queries as the
//! in-house reader:
//!
//! * Earth's heliocentric state (geocenter relative to the Sun),
//! * any supported body's heliocentric state (relative to the Sun),
//!
//! both in astronomical units and AU per day, expressed in equatorial mean J2000.
//! Frame selection (equatorial vs ecliptic) is applied by the caller in
//! [`crate::jpl_ephem`], not here.
//!
//! # Structure
//!
//! The numeric transforms are free, side-effect-free functions
//! (`resolve_translation_target`, `cartesian_km_s_to_au_day`) that are unit- and
//! property-tested without any kernel. The only method that reads the loaded
//! kernel is `AniseEphem::heliocentric_state_km`; it is deliberately tiny so that
//! all arithmetic stays in the pure functions.
//!
//! # Units and frames
//!
//! The kernel is queried for a geometric translation (no aberration correction),
//! yielding kilometres and kilometres per second in equatorial mean J2000. This
//! module converts kilometres to AU with [`crate::constants::AU`] and kilometres
//! per second to AU per day with the same factor times the number of seconds in a
//! day.

use std::sync::Arc;

use anise::prelude::{Almanac, Frame};
use camino::Utf8Path;
use hifitime::Epoch;
use nalgebra::Vector3;

use crate::constants::{AU, ROT_ECLMJ2000_TO_EQUMJ2000};
use crate::jpl_ephem::download_jpl_file::EphemFilePath;
use crate::jpl_ephem::naif::naif_ids::{solar_system_bary::SolarSystemBary, NaifIds};
use crate::outfit_errors::OutfitError;

/// NAIF integer identifier of the Sun, used as the heliocentric observer.
const SUN_NAIF_ID: i32 = 10;

/// NAIF integer identifier of the Earth geocenter.
const EARTH_GEOCENTER_NAIF_ID: i32 = 399;

/// Kilometres to astronomical units.
const KM_TO_AU: f64 = 1.0 / AU;

/// Kilometres per second to astronomical units per day.
const KM_S_TO_AU_DAY: f64 = 86_400.0 / AU;

/// Resolve the NAIF integer identifier used to translate `body` against the Sun.
///
/// # Arguments
///
/// * `body` – solar-system body to look up.
///
/// # Returns
///
/// The NAIF integer identifier of `body`.
///
/// # Errors
///
/// Returns [`OutfitError::EphemerisBodyNotSupported`] when `body` is the Solar
/// System Barycenter, which is a reference point rather than a physical body and
/// therefore cannot be requested as a translation target.
fn resolve_translation_target(body: NaifIds) -> Result<i32, OutfitError> {
    match body {
        NaifIds::SSB(SolarSystemBary::SSB) => Err(OutfitError::EphemerisBodyNotSupported(
            "Solar System Barycenter is not a physical body".to_string(),
        )),
        other => Ok(i32::from(other)),
    }
}

/// Convert a raw kilometre / kilometre-per-second state to AU and AU per day.
///
/// # Arguments
///
/// * `radius_km` – position components in kilometres.
/// * `velocity_km_s` – velocity components in kilometres per second.
///
/// # Returns
///
/// `(position, velocity)` with `position` in AU and `velocity` in AU per day, in
/// the same axes as the input.
fn cartesian_km_s_to_au_day(
    radius_km: [f64; 3],
    velocity_km_s: [f64; 3],
) -> (Vector3<f64>, Vector3<f64>) {
    let position_au = Vector3::new(
        radius_km[0] * KM_TO_AU,
        radius_km[1] * KM_TO_AU,
        radius_km[2] * KM_TO_AU,
    );
    let velocity_au_day = Vector3::new(
        velocity_km_s[0] * KM_S_TO_AU_DAY,
        velocity_km_s[1] * KM_S_TO_AU_DAY,
        velocity_km_s[2] * KM_S_TO_AU_DAY,
    );
    (position_au, velocity_au_day)
}

/// True for every body served by the main-belt asteroid supplementary kernel
/// (`codes_300ast_20100725.bsp`, NAIF ids `2_000_000..=2_999_999`).
///
/// That kernel stores its states in the ECLIPJ2000_DE405 frame (ecliptic mean
/// J2000 at the DE405 obliquity), not equatorial mean J2000/ICRF. ANISE's
/// translation-only API returns whatever frame the segment is stored in
/// without rotating it, so bodies from this kernel need the correction in
/// [`correct_ecliptic_de405_frame`] before they can be treated as equatorial.
///
/// # Arguments
///
/// * `target_naif_id` – NAIF integer identifier being queried.
///
/// # Returns
///
/// `true` if `target_naif_id` falls in the main-belt asteroid id range.
fn needs_ecliptic_de405_correction(target_naif_id: i32) -> bool {
    (2_000_000..3_000_000).contains(&target_naif_id)
}

/// Rotate a heliocentric state from ECLIPJ2000_DE405 to equatorial mean J2000,
/// when the body's source kernel requires it.
///
/// The DE405 obliquity (84381.412″) that defines ECLIPJ2000_DE405 is, to the
/// precision both are stated at, the same obliquity as
/// [`ROT_ECLMJ2000_TO_EQUMJ2000`](crate::constants::ROT_ECLMJ2000_TO_EQUMJ2000),
/// so that existing rotation is reused rather than deriving a new one.
///
/// # Arguments
///
/// * `target_naif_id` – NAIF integer identifier the state belongs to.
/// * `position_au` – position, AU, as returned by [`cartesian_km_s_to_au_day`].
/// * `velocity_au_day` – velocity, AU per day, same source.
///
/// # Returns
///
/// `(position, velocity)` in equatorial mean J2000: rotated when
/// [`needs_ecliptic_de405_correction`] is `true` for `target_naif_id`,
/// returned unchanged otherwise.
fn correct_ecliptic_de405_frame(
    target_naif_id: i32,
    position_au: Vector3<f64>,
    velocity_au_day: Vector3<f64>,
) -> (Vector3<f64>, Vector3<f64>) {
    if needs_ecliptic_de405_correction(target_naif_id) {
        (
            ROT_ECLMJ2000_TO_EQUMJ2000 * position_au,
            ROT_ECLMJ2000_TO_EQUMJ2000 * velocity_au_day,
        )
    } else {
        (position_au, velocity_au_day)
    }
}

/// Extract the on-disk path of an [`EphemFilePath`] that denotes a NAIF SPK
/// kernel, rejecting the legacy-DE variant this backend cannot read.
///
/// # Arguments
///
/// * `path` – resolved ephemeris file.
///
/// # Returns
///
/// The kernel's on-disk path.
///
/// # Errors
///
/// Returns [`OutfitError::InvalidJPLEphemFileSource`] if `path` denotes a
/// legacy DE binary (`EphemFilePath::JPLHorizon`).
fn spk_path(path: &EphemFilePath) -> Result<&Utf8Path, OutfitError> {
    match path {
        EphemFilePath::Naif(kernel_path, _) => Ok(kernel_path),
        EphemFilePath::MainBeltAsteroids(kernel_path) => Ok(kernel_path),
        EphemFilePath::JPLHorizon(..) => Err(OutfitError::InvalidJPLEphemFileSource(
            "the ANISE backend reads NAIF SPK kernels only; use a \"naif:DE###\" source"
                .to_string(),
        )),
    }
}

/// Planetary ephemeris backend holding a NAIF SPK kernel parsed by ANISE.
///
/// The kernel is reference-counted so that cloning the handle is cheap; the
/// parsed kernel itself is never mutated after construction.
#[derive(Clone)]
pub struct AniseEphem {
    /// Parsed kernel context, shared between clones.
    almanac: Arc<Almanac>,
}

impl core::fmt::Debug for AniseEphem {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "AniseEphem({})", self.almanac)
    }
}

impl AniseEphem {
    /// Build the backend from an already-resolved on-disk ephemeris path.
    ///
    /// # Arguments
    ///
    /// * `path` – resolved ephemeris file; it must denote a NAIF SPK kernel.
    ///
    /// # Returns
    ///
    /// A ready-to-query backend holding the parsed kernel.
    ///
    /// # Errors
    ///
    /// * [`OutfitError::InvalidJPLEphemFileSource`] – the path denotes a legacy DE
    ///   binary, which this backend cannot read.
    /// * [`OutfitError::AniseEphemerisError`] – the kernel could not be parsed.
    pub fn from_file_path(path: &EphemFilePath) -> Result<Self, OutfitError> {
        let almanac = Almanac::new(spk_path(path)?.as_str())
            .map_err(|err| OutfitError::AniseEphemerisError(err.to_string()))?;
        Ok(Self {
            almanac: Arc::new(almanac),
        })
    }

    /// Load an additional SPK kernel into this backend (e.g. the main-belt
    /// asteroid supplementary kernel), in place.
    ///
    /// After this call, this handle can resolve every body covered by either
    /// the kernel(s) it already had loaded or `path`.
    ///
    /// # Arguments
    ///
    /// * `path` – resolved ephemeris file to add; it must denote a NAIF SPK
    ///   kernel.
    ///
    /// # Errors
    ///
    /// * [`OutfitError::InvalidJPLEphemFileSource`] – the path denotes a
    ///   legacy DE binary, which this backend cannot read.
    /// * [`OutfitError::AniseEphemerisError`] – the kernel could not be parsed.
    pub fn with_supplementary_kernel(&mut self, path: &EphemFilePath) -> Result<(), OutfitError> {
        let almanac = (*self.almanac)
            .clone()
            .load(spk_path(path)?.as_str())
            .map_err(|err| OutfitError::AniseEphemerisError(err.to_string()))?;
        self.almanac = Arc::new(almanac);
        Ok(())
    }

    /// Heliocentric equatorial-mean-J2000 state of a NAIF body relative to the
    /// Sun, as raw kilometre / kilometre-per-second component arrays.
    ///
    /// This is the only method that reads the loaded kernel.
    ///
    /// # Arguments
    ///
    /// * `target_naif_id` – NAIF integer identifier of the body to look up.
    /// * `epoch` – evaluation epoch.
    ///
    /// # Returns
    ///
    /// `(radius_km, velocity_km_s)` component arrays, heliocentric, equatorial
    /// mean J2000.
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError::AniseEphemerisError`] when the kernel does not cover
    /// the epoch or the body chain cannot be resolved against the Sun.
    fn heliocentric_state_km(
        &self,
        target_naif_id: i32,
        epoch: &Epoch,
    ) -> Result<([f64; 3], [f64; 3]), OutfitError> {
        let state = self
            .almanac
            .translate_geometric(
                Frame::from_ephem_j2000(target_naif_id),
                Frame::from_ephem_j2000(SUN_NAIF_ID),
                *epoch,
            )
            .map_err(|err| OutfitError::AniseEphemerisError(err.to_string()))?;
        Ok((
            [state.radius_km[0], state.radius_km[1], state.radius_km[2]],
            [
                state.velocity_km_s[0],
                state.velocity_km_s[1],
                state.velocity_km_s[2],
            ],
        ))
    }

    /// Earth's heliocentric equatorial-mean-J2000 state, in AU and AU per day.
    ///
    /// The Earth geocenter is translated against the Sun.
    ///
    /// # Arguments
    ///
    /// * `epoch` – evaluation epoch.
    /// * `compute_velocity` – whether to also return the velocity.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` with `position` in AU and `velocity` in AU per day;
    /// `velocity` is `Some` only when `compute_velocity` is `true`.
    ///
    /// # Panics
    ///
    /// Panics when the kernel lookup fails (for example an out-of-coverage
    /// epoch); this matches the infallible contract of
    /// [`crate::jpl_ephem::JPLEphem::earth_ephemeris`]. Use
    /// [`AniseEphem::body_ephemeris_equatorial`] for a fallible query.
    pub(crate) fn earth_ephemeris_equatorial(
        &self,
        epoch: &Epoch,
        compute_velocity: bool,
    ) -> (Vector3<f64>, Option<Vector3<f64>>) {
        let (radius_km, velocity_km_s) = self
            .heliocentric_state_km(EARTH_GEOCENTER_NAIF_ID, epoch)
            .unwrap_or_else(|err| panic!("ANISE Earth ephemeris lookup failed: {err}"));
        let (position_au, velocity_au_day) = cartesian_km_s_to_au_day(radius_km, velocity_km_s);
        (position_au, compute_velocity.then_some(velocity_au_day))
    }

    /// Heliocentric equatorial-mean-J2000 state of `body`, in AU and AU per day.
    ///
    /// # Arguments
    ///
    /// * `body` – body to look up.
    /// * `epoch` – evaluation epoch.
    ///
    /// # Returns
    ///
    /// `(position, velocity)` with `position` in AU and `velocity` in AU per day,
    /// heliocentric, equatorial mean J2000.
    ///
    /// # Errors
    ///
    /// * [`OutfitError::EphemerisBodyNotSupported`] – `body` is the Solar System
    ///   Barycenter.
    /// * [`OutfitError::AniseEphemerisError`] – the kernel does not cover the
    ///   epoch or the body chain cannot be resolved.
    pub(crate) fn body_ephemeris_equatorial(
        &self,
        body: NaifIds,
        epoch: &Epoch,
    ) -> Result<(Vector3<f64>, Vector3<f64>), OutfitError> {
        let target_naif_id = resolve_translation_target(body)?;
        let (radius_km, velocity_km_s) = self.heliocentric_state_km(target_naif_id, epoch)?;
        let (position_au, velocity_au_day) = cartesian_km_s_to_au_day(radius_km, velocity_km_s);
        Ok(correct_ecliptic_de405_frame(
            target_naif_id,
            position_au,
            velocity_au_day,
        ))
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod pure_function_tests {
    use super::*;
    use crate::jpl_ephem::naif::naif_ids::{
        planet_bary::PlanetaryBary, satellite_mass::SatelliteMassCenter,
    };
    use proptest::prelude::*;

    #[test]
    fn resolve_translation_target_maps_known_bodies() {
        assert_eq!(
            resolve_translation_target(NaifIds::PB(PlanetaryBary::Mars)).unwrap(),
            4
        );
        assert_eq!(
            resolve_translation_target(NaifIds::SMC(SatelliteMassCenter::Moon)).unwrap(),
            301
        );
        assert_eq!(
            resolve_translation_target(NaifIds::SSB(SolarSystemBary::Sun)).unwrap(),
            10
        );
    }

    #[test]
    fn resolve_translation_target_rejects_the_barycenter() {
        let err = resolve_translation_target(NaifIds::SSB(SolarSystemBary::SSB)).unwrap_err();
        assert!(matches!(err, OutfitError::EphemerisBodyNotSupported(_)));
    }

    #[test]
    fn cartesian_conversion_maps_one_au() {
        let (position_au, velocity_au_day) =
            cartesian_km_s_to_au_day([AU, 0.0, 0.0], [AU / 86_400.0, 0.0, 0.0]);
        assert!((position_au[0] - 1.0).abs() < 1e-12);
        assert!(position_au[1].abs() < 1e-12 && position_au[2].abs() < 1e-12);
        assert!((velocity_au_day[0] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn cartesian_conversion_maps_zero_to_zero() {
        let (position_au, velocity_au_day) =
            cartesian_km_s_to_au_day([0.0, 0.0, 0.0], [0.0, 0.0, 0.0]);
        assert_eq!(position_au, Vector3::zeros());
        assert_eq!(velocity_au_day, Vector3::zeros());
    }

    #[test]
    fn ecliptic_correction_applies_only_to_asteroid_ids() {
        assert!(needs_ecliptic_de405_correction(2_000_001)); // Ceres
        assert!(needs_ecliptic_de405_correction(2_999_999));
        assert!(!needs_ecliptic_de405_correction(10)); // Sun
        assert!(!needs_ecliptic_de405_correction(399)); // Earth
        assert!(!needs_ecliptic_de405_correction(1_999_999));
        assert!(!needs_ecliptic_de405_correction(3_000_000));
    }

    #[test]
    fn ecliptic_correction_is_identity_for_a_planet() {
        let p = Vector3::new(1.3, -0.7, 0.42);
        let v = Vector3::new(-0.011, 0.008, 0.003);
        let (rp, rv) = correct_ecliptic_de405_frame(4, p, v); // Mars barycenter id
        assert_eq!(rp, p);
        assert_eq!(rv, v);
    }

    #[test]
    fn ecliptic_correction_rotates_an_asteroid_state() {
        let p = Vector3::new(1.3, -0.7, 0.42);
        let v = Vector3::new(-0.011, 0.008, 0.003);
        let (rp, rv) = correct_ecliptic_de405_frame(2_000_001, p, v); // Ceres
        assert_eq!(rp, ROT_ECLMJ2000_TO_EQUMJ2000 * p);
        assert_eq!(rv, ROT_ECLMJ2000_TO_EQUMJ2000 * v);
        assert_ne!(rp, p, "the correction must actually rotate the state");
    }

    proptest! {
        /// For every non-barycenter body the target id is exactly the raw NAIF
        /// integer code.
        #[test]
        fn resolve_translation_target_agrees_with_raw_code(
            selector in 0usize..11,
        ) {
            let bodies = [
                NaifIds::SSB(SolarSystemBary::Sun),
                NaifIds::PB(PlanetaryBary::Mercury),
                NaifIds::PB(PlanetaryBary::Venus),
                NaifIds::PB(PlanetaryBary::EarthMoon),
                NaifIds::PB(PlanetaryBary::Mars),
                NaifIds::PB(PlanetaryBary::Jupiter),
                NaifIds::PB(PlanetaryBary::Saturn),
                NaifIds::PB(PlanetaryBary::Uranus),
                NaifIds::PB(PlanetaryBary::Neptune),
                NaifIds::PB(PlanetaryBary::Pluto),
                NaifIds::SMC(SatelliteMassCenter::Moon),
            ];
            let body = bodies[selector];
            prop_assert_eq!(resolve_translation_target(body).unwrap(), i32::from(body));
        }

        /// Scaling the kilometre state by a positive factor scales the AU state
        /// by the same factor.
        #[test]
        fn cartesian_conversion_is_linear(
            x in -1e9_f64..1e9,
            y in -1e9_f64..1e9,
            z in -1e9_f64..1e9,
            vx in -1e3_f64..1e3,
            vy in -1e3_f64..1e3,
            vz in -1e3_f64..1e3,
            scale in 0.25_f64..4.0,
        ) {
            let (p1, v1) = cartesian_km_s_to_au_day([x, y, z], [vx, vy, vz]);
            let (p2, v2) = cartesian_km_s_to_au_day(
                [x * scale, y * scale, z * scale],
                [vx * scale, vy * scale, vz * scale],
            );
            prop_assert!((p2 - p1 * scale).norm() < 1e-6 * (1.0 + p1.norm()));
            prop_assert!((v2 - v1 * scale).norm() < 1e-9 * (1.0 + v1.norm()));
        }

        /// A finite kilometre state always maps to a finite AU state.
        #[test]
        fn cartesian_conversion_stays_finite(
            x in -1e12_f64..1e12,
            y in -1e12_f64..1e12,
            z in -1e12_f64..1e12,
        ) {
            let (p, _) = cartesian_km_s_to_au_day([x, y, z], [0.0, 0.0, 0.0]);
            prop_assert!(p.iter().all(|c| c.is_finite()));
        }

        /// A rotation preserves vector norms, whether or not the correction
        /// actually fires for the given id.
        #[test]
        fn ecliptic_correction_preserves_norm(
            target_naif_id in prop_oneof![
                Just(10), Just(399), Just(2_000_001), Just(2_001_467),
            ],
            x in -50.0_f64..50.0,
            y in -50.0_f64..50.0,
            z in -50.0_f64..50.0,
        ) {
            let p = Vector3::new(x, y, z);
            let v = Vector3::zeros();
            let (rp, _) = correct_ecliptic_de405_frame(target_naif_id, p, v);
            prop_assert!((rp.norm() - p.norm()).abs() < 1e-12 * (1.0 + p.norm()));
        }
    }
}

#[cfg(all(test, feature = "ephem-anise"))]
mod kernel_backed_tests {
    use super::*;
    use crate::jpl_ephem::naif::naif_ids::planet_bary::PlanetaryBary;
    use crate::jpl_ephem::EphemerisFrame;
    use crate::test_fixture::JPL_EPHEM_ANISE;
    use hifitime::TimeScale;

    /// MJD-TT epoch comfortably inside DE440 coverage.
    fn epoch_tt(mjd_tt: f64) -> Epoch {
        Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT)
    }

    /// The ANISE velocity is the time derivative of the ANISE position: this is
    /// the check the built-in NAIF backend currently fails.
    #[test]
    fn body_velocity_matches_central_difference() {
        let jpl = &*JPL_EPHEM_ANISE;
        let body = NaifIds::PB(PlanetaryBary::Mars);
        let mjd_tt = 59_500.0_f64;
        let h = 0.25_f64;

        let v = jpl
            .body_ephemeris(body, &epoch_tt(mjd_tt), EphemerisFrame::Equatorial)
            .unwrap()
            .1;
        let r_plus = jpl
            .body_ephemeris(body, &epoch_tt(mjd_tt + h), EphemerisFrame::Equatorial)
            .unwrap()
            .0;
        let r_minus = jpl
            .body_ephemeris(body, &epoch_tt(mjd_tt - h), EphemerisFrame::Equatorial)
            .unwrap()
            .0;
        let v_fd = (r_plus - r_minus) / (2.0 * h);

        let rel_err = (v - v_fd).norm() / v.norm();
        assert!(rel_err < 1e-6, "relative velocity error {rel_err:e}");
    }

    /// The Solar System Barycenter is rejected as a body.
    #[test]
    fn barycenter_is_rejected() {
        let err = JPL_EPHEM_ANISE
            .body_ephemeris(
                NaifIds::SSB(SolarSystemBary::SSB),
                &epoch_tt(59_000.0),
                EphemerisFrame::Equatorial,
            )
            .unwrap_err();
        assert!(matches!(err, OutfitError::EphemerisBodyNotSupported(_)));
    }

    /// `earth_ephemeris` returns a plausible ~1 AU heliocentric position and a
    /// velocity only when requested.
    #[test]
    fn earth_state_is_about_one_au() {
        let (position, velocity) =
            JPL_EPHEM_ANISE.earth_ephemeris(&epoch_tt(59_000.0), EphemerisFrame::Equatorial, false);
        assert!((position.norm() - 1.0).abs() < 0.05);
        assert!(velocity.is_none());

        let (_, velocity) =
            JPL_EPHEM_ANISE.earth_ephemeris(&epoch_tt(59_000.0), EphemerisFrame::Equatorial, true);
        let velocity = velocity.expect("velocity requested");
        // Earth's mean orbital speed is ~0.0172 AU/day.
        assert!((velocity.norm() - 0.0172).abs() < 0.002);
    }
}

#[cfg(all(test, feature = "ephem-anise", feature = "ephem-builtin"))]
mod parity_tests {
    use super::*;
    use crate::jpl_ephem::naif::naif_ids::planet_bary::PlanetaryBary;
    use crate::jpl_ephem::EphemerisFrame;
    use crate::test_fixture::{JPL_EPHEM_ANISE, JPL_EPHEM_HORIZON};
    use hifitime::TimeScale;

    fn epoch_tt(mjd_tt: f64) -> Epoch {
        Epoch::from_mjd_in_time_scale(mjd_tt, TimeScale::TT)
    }

    /// ANISE and the built-in Horizon reader agree on heliocentric planet
    /// positions to well below a metre-scale fraction of an AU.
    #[test]
    fn body_positions_agree_with_horizon() {
        let bodies = [
            NaifIds::PB(PlanetaryBary::Mars),
            NaifIds::PB(PlanetaryBary::Jupiter),
            NaifIds::PB(PlanetaryBary::Saturn),
        ];
        for mjd_tt in [58_900.0_f64, 59_500.0, 60_100.0] {
            let epoch = epoch_tt(mjd_tt);
            for &body in &bodies {
                let anise = JPL_EPHEM_ANISE
                    .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                    .unwrap()
                    .0;
                let horizon = JPL_EPHEM_HORIZON
                    .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
                    .unwrap()
                    .0;
                let diff = (anise - horizon).norm();
                assert!(diff < 1e-8, "{body:?} at {mjd_tt}: |Δp| = {diff:e} AU");
            }
        }
    }

    /// ANISE and the built-in Horizon reader agree on Earth's heliocentric
    /// position and velocity.
    #[test]
    fn earth_state_agrees_with_horizon() {
        for mjd_tt in [58_900.0_f64, 59_500.0, 60_100.0] {
            let epoch = epoch_tt(mjd_tt);
            let (anise_p, anise_v) =
                JPL_EPHEM_ANISE.earth_ephemeris(&epoch, EphemerisFrame::Equatorial, true);
            let (horizon_p, horizon_v) =
                JPL_EPHEM_HORIZON.earth_ephemeris(&epoch, EphemerisFrame::Equatorial, true);
            assert!((anise_p - horizon_p).norm() < 1e-9, "position at {mjd_tt}");
            let dv = (anise_v.unwrap() - horizon_v.unwrap()).norm();
            assert!(dv < 1e-10, "velocity at {mjd_tt}: |Δv| = {dv:e} AU/day");
        }
    }

    /// The ecliptic branch is the equatorial state rotated by the shared
    /// reference rotation, on both backends alike.
    #[test]
    fn ecliptic_selection_is_consistent() {
        let epoch = epoch_tt(59_321.0);
        let body = NaifIds::PB(PlanetaryBary::Jupiter);
        let (equ, _) = JPL_EPHEM_ANISE
            .body_ephemeris(body, &epoch, EphemerisFrame::Equatorial)
            .unwrap();
        let (ecl, _) = JPL_EPHEM_ANISE
            .body_ephemeris(body, &epoch, EphemerisFrame::Ecliptic)
            .unwrap();
        let (expected, _) = crate::jpl_ephem::equ_state_to_ecl(equ, equ);
        assert!((ecl - expected).norm() < 1e-13);
    }
}

#[cfg(all(test, feature = "ephem-anise"))]
mod main_belt_asteroid_tests {
    use super::*;
    use crate::jpl_ephem::naif::naif_ids::main_belt::AsteroidNumber;
    use crate::jpl_ephem::EphemerisFrame;
    use crate::test_fixture::JPL_EPHEM_ANISE_WITH_ASTEROIDS;
    use hifitime::TimeScale;

    /// Heliocentric equatorial mean J2000 (ICRF) state of Ceres at
    /// 2021-Sep-19 00:00 TDB (JD 2459476.5), from JPL Horizons
    /// (`COMMAND='1;' CENTER='500@10' REF_PLANE=FRAME OUT_UNITS=AU-D`,
    /// queried 2026-09-11): position in AU, velocity in AU/day.
    const CERES_HORIZONS_POSITION_AU: [f64; 3] = [
        1.782_341_075_522_114,
        2.073_176_490_439_472,
        0.614_852_658_891_694_7,
    ];
    const CERES_HORIZONS_VELOCITY_AU_DAY: [f64; 3] = [
        -0.008_136_059_326_451_917,
        0.004_769_609_050_730_651,
        0.003_905_915_532_942_483,
    ];

    /// The epoch the constants above were sampled at (MJD 59476.0). Tagged TT
    /// rather than TDB: the two differ by a couple of milliseconds, far below
    /// the tolerance used to compare against Horizons below.
    fn ceres_reference_epoch() -> Epoch {
        Epoch::from_mjd_in_time_scale(59_476.0, TimeScale::TT)
    }

    /// Querying a main-belt asteroid works once the supplementary kernel is
    /// loaded, and the resolved position is a plausible main-belt distance
    /// (Ceres: semi-major axis ~2.77 AU, eccentricity ~0.08).
    #[test]
    fn ceres_position_is_a_plausible_main_belt_distance() {
        let (position, _) = JPL_EPHEM_ANISE_WITH_ASTEROIDS
            .body_ephemeris(
                NaifIds::AST(AsteroidNumber::CERES),
                &ceres_reference_epoch(),
                EphemerisFrame::Equatorial,
            )
            .unwrap();
        assert!(
            (2.5..3.0).contains(&position.norm()),
            "|r| = {}",
            position.norm()
        );
    }

    /// The resolved state agrees, component by component, with an
    /// independent JPL Horizons state vector at the same epoch. Comparing
    /// components (not just the norm) is essential here: a state that was
    /// left in the kernel's native ECLIPJ2000_DE405 frame instead of being
    /// rotated to equatorial mean J2000 would still have the right distance
    /// from the Sun, just pointed in the wrong direction.
    #[test]
    fn ceres_state_agrees_with_horizons() {
        let (position, velocity) = JPL_EPHEM_ANISE_WITH_ASTEROIDS
            .body_ephemeris(
                NaifIds::AST(AsteroidNumber::CERES),
                &ceres_reference_epoch(),
                EphemerisFrame::Equatorial,
            )
            .unwrap();

        let expected_position = Vector3::from(CERES_HORIZONS_POSITION_AU);
        let expected_velocity = Vector3::from(CERES_HORIZONS_VELOCITY_AU_DAY);

        // `codes_300ast_20100725.bsp` is an independent, ~2010-era orbit fit
        // rather than a DE440-consistent one, so the tolerance is looser than
        // the sub-1e-8 AU parity checks against the built-in reader elsewhere
        // in this file; observed agreement is ~7e-7 AU / ~2e-9 AU/day, so
        // 1e-4 / 1e-6 leaves headroom while still catching a missing or
        // doubled frame rotation outright (that would be off by ~1 AU or a
        // large relative angle, not a small fraction of a percent).
        let dp = (position - expected_position).norm();
        let dv = (velocity - expected_velocity).norm();
        assert!(dp < 1e-4, "|Δposition| = {dp:e} AU");
        assert!(dv < 1e-6, "|Δvelocity| = {dv:e} AU/day");
    }

    /// The Solar System Barycenter still cannot be requested as a body, even
    /// once the asteroid kernel is loaded.
    #[test]
    fn barycenter_is_still_rejected() {
        let err = JPL_EPHEM_ANISE_WITH_ASTEROIDS
            .body_ephemeris(
                NaifIds::SSB(SolarSystemBary::SSB),
                &ceres_reference_epoch(),
                EphemerisFrame::Equatorial,
            )
            .unwrap_err();
        assert!(matches!(err, OutfitError::EphemerisBodyNotSupported(_)));
    }

    /// Loading the supplementary kernel does not break access to bodies from
    /// the primary DE440 kernel.
    #[test]
    fn planets_still_resolve_after_loading_asteroids() {
        let (position, _) = JPL_EPHEM_ANISE_WITH_ASTEROIDS.earth_ephemeris(
            &ceres_reference_epoch(),
            EphemerisFrame::Equatorial,
            false,
        );
        assert!((position.norm() - 1.0).abs() < 0.05);
    }
}
