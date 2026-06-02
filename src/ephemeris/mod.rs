//! Public façade for ephemeris computation.
//!
//! This module exposes [`ApparentPosition`], [`BodyGeometry`],
//! [`EphemerisConfig`], and [`AberrationOrder`], plus the request/result
//! system built around [`EphemerisRequest`] and [`EphemerisResult`].
//!
//! # Quick-start
//!
//! Build a request, add one or more `(observer, mode)` pairs, then call
//! [`OrbitalElements::compute`]:
//!
//! ```rust,ignore
//! use outfit::{
//!     Combined, EphemerisConfig, EphemerisMode, EphemerisRequest,
//! };
//! use hifitime::{Epoch, Duration};
//!
//! let result = elements.compute(
//!     &EphemerisRequest::<Combined>::new(EphemerisConfig::default())
//!         .add(observer_a, EphemerisMode::Range {
//!             start: Epoch::from_mjd_tt(60310.0),
//!             end:   Epoch::from_mjd_tt(60340.0),
//!             step:  Duration::from_days(1.0),
//!         })
//!         .add(observer_b, EphemerisMode::At(vec![t1, t2, t3]))
//!         .add(observer_c, EphemerisMode::Single(t)),
//!     &jpl,
//!     &ut1,
//! );
//!
//! for entry in result.successes() {
//!     let (pos, geo) = entry.result.as_ref().unwrap();
//!     println!("{}: RA={:.4} phase={:.4}", entry.epoch, pos.coord.ra, geo.phase_angle);
//! }
//! ```
//!
//! # Output kinds
//!
//! The type parameter on [`EphemerisRequest`] selects what is computed:
//!
//! | Marker | Output per epoch |
//! |---|---|
//! | [`Position`]  | [`ApparentPosition`] |
//! | [`Geometry`]  | [`BodyGeometry`] |
//! | [`Combined`]  | `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` |
//!
//! # Generation modes
//!
//! Each observer in the request is paired with an [`EphemerisMode`]:
//!
//! | Variant | Epochs |
//! |---|---|
//! | [`EphemerisMode::Single`] | Exactly one epoch |
//! | [`EphemerisMode::Range`]  | Uniform grid |
//! | [`EphemerisMode::At`]     | Arbitrary list |
//!
//! # Coordinate conventions
//!
//! - Positions in **AU**, velocities in **AU/day**, angles in **radians**.
//! - Intermediate frames: **ecliptic mean J2000**.
//! - Final output: **equatorial mean J2000** (RA ∈ \[0, 2π), Dec ∈ (−π/2, π/2)).
//! - Time input: any [`hifitime::Epoch`]; converted to **MJD TT** internally.
//!
//! # Pipeline overview
//!
//! ```text
//! OrbitalElements
//!       │
//!       ▼  to_equinoctial()
//! EquinoctialElements
//!       │
//!       ▼  propagate  (TwoBody or NBody)
//! (heliocentric position, velocity) [ecliptic J2000]
//!       │
//!       ▼  ROT_ECLMJ2000_TO_EQUMJ2000
//! (heliocentric position, velocity) [equatorial J2000]
//!       │
//!       ├──────────────────────────────────────────┐
//!       ▼                                          ▼
//! topocentric vector                        distances
//! (body − observer)                   (geocentric, heliocentric)
//!       │
//!       ▼  correct_aberration_{first,second}_order()
//! aberration-corrected line-of-sight
//!       │
//!       ├──────────────────────────────────────────────────────────┐
//!       ▼                                                          ▼
//! ApparentPosition                                          BodyGeometry
//! { coord, geocentric_dist, heliocentric_dist }   { phase_angle, solar_elongation,
//!                                                   radial_velocity, d_ra_dt, d_dec_dt }
//! ```

pub(crate) mod aberration;
pub(crate) mod apparent_position;
pub mod batch;
pub(crate) mod geometry;
pub(crate) mod observation_ephemeris;
pub mod request;
pub mod result;
pub use aberration::AberrationOrder;
pub use apparent_position::ApparentPosition;
pub use batch::FullOrbitResultExt;
pub use geometry::BodyGeometry;
pub use request::{
    Combined, EphemerisMode, EphemerisOutputKind, EphemerisRequest, Geometry, ObserverRequest,
    Position,
};
pub use result::{EphemerisEntry, EphemerisResult};

use hifitime::ut1::Ut1Provider;

use crate::{
    cache::observer_fixed_cache::ObserverFixedCache,
    ephemeris::observation_ephemeris::check_elliptical_orbit, propagator::PropagatorKind, JPLEphem,
    OrbitalElements, OutfitError,
};

// ---------------------------------------------------------------------------
// EphemerisConfig
// ---------------------------------------------------------------------------

/// Configuration for ephemeris computation.
///
/// Controls which propagation strategy and aberration correction are applied
/// when computing a predicted apparent position.  The default uses the
/// analytic two-body (Keplerian) propagator and the first-order aberration
/// correction, which are fast and sufficient for most targets.
#[derive(Debug, Clone, Default)]
pub struct EphemerisConfig {
    /// Propagator to use for computing predicted positions.
    ///
    /// - [`PropagatorKind::TwoBody`] (default): analytic Keplerian propagation.
    /// - [`PropagatorKind::NBody`]: numerical DOP853 N-body integration with
    ///   user-specified perturbing bodies.
    pub propagator: PropagatorKind,

    /// Aberration correction order.
    ///
    /// - [`AberrationOrder::First`] (default): linear light-travel-time shift.
    /// - [`AberrationOrder::Second`]: two-step Keplerian back-propagation.
    pub aberration: AberrationOrder,
}

// ---------------------------------------------------------------------------
// OrbitalElements — ephemeris entry point
// ---------------------------------------------------------------------------

impl OrbitalElements {
    /// Execute an [`EphemerisRequest`] and return an [`EphemerisResult`].
    ///
    /// Iterates over every `(observer, mode)` pair in `request`, expands the
    /// mode into a concrete list of epochs, and computes the ephemeris for
    /// each `(epoch, observer)` combination using the output kind `O`.
    ///
    /// Errors at individual epochs are recorded in the result rather than
    /// aborting the whole computation.
    ///
    /// # Arguments
    ///
    /// - `request` – Typed request carrying observers, modes, and config.
    /// - `jpl`     – JPL ephemeris.
    /// - `ut1`     – UT1 provider.
    ///
    /// # Returns
    ///
    /// An [`EphemerisResult<O::Output>`] whose entries are ordered: all epochs
    /// for the first observer, then all epochs for the second observer, etc.
    ///
    /// # Errors (per-entry)
    ///
    /// - Hyperbolic or parabolic orbit (eccentricity ≥ 1).
    /// - Propagation failure.
    /// - JPL ephemeris data unavailable for the requested epoch.
    /// - Conversion to equinoctial elements fails.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let result = elements.compute(
    ///     &EphemerisRequest::<Combined>::new(EphemerisConfig::default())
    ///         .add(observer, EphemerisMode::Range { start, end, step }),
    ///     &jpl, &ut1,
    /// );
    /// for entry in result.successes() {
    ///     let (pos, geo) = entry.result.as_ref().unwrap();
    ///     println!("{}: RA={:.4}", entry.epoch, pos.coord.ra);
    /// }
    /// ```
    pub fn compute<O: EphemerisOutputKind>(
        &self,
        request: &EphemerisRequest<O>,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
    ) -> EphemerisResult<O::Output> {
        // Convert to equinoctial once; propagation failures are per-entry.
        let equi = match self.to_equinoctial_for_ephemeris() {
            Ok(e) => e,
            Err(err) => {
                // If conversion fails, every entry is an error.
                let total: usize = request
                    .observers
                    .iter()
                    .map(|r| r.mode.epochs().len())
                    .sum();
                let mut result = EphemerisResult::with_capacity(total);
                for obs_req in &request.observers {
                    for epoch in obs_req.mode.epochs() {
                        result.push(
                            epoch,
                            obs_req.observer.clone(),
                            Err(OutfitError::InvalidConversion(err.to_string())),
                        );
                    }
                }
                return result;
            }
        };

        // Optim 3 — check eccentricity once before any epoch loop.
        // Avoids recomputing sqrt(h²+k²) for every (epoch, observer) pair.
        // If the orbit is hyperbolic/parabolic every entry would fail anyway,
        // so we short-circuit immediately with a uniform error result.
        if let Err(err) = check_elliptical_orbit(&equi) {
            let total: usize = request
                .observers
                .iter()
                .map(|r| r.mode.epochs().len())
                .sum();
            let mut result = EphemerisResult::with_capacity(total);
            for obs_req in &request.observers {
                for epoch in obs_req.mode.epochs() {
                    result.push(
                        epoch,
                        obs_req.observer.clone(),
                        Err(OutfitError::InvalidConversion(err.to_string())),
                    );
                }
            }
            return result;
        }

        let total: usize = request
            .observers
            .iter()
            .map(|r| r.mode.epochs().len())
            .sum();
        let mut result = EphemerisResult::with_capacity(total);

        for obs_req in &request.observers {
            // Optim — build ObserverFixedCache once per observer slot.
            //fixed_cache
            // ObserverFixedCache holds the body-fixed position and velocity of
            // the observing site (sin/cos of longitude, ρ factors, cross-product
            // for sidereal rotation).  These quantities are epoch-invariant: they
            // depend only on the observer's geodetic coordinates.  Building the
            // cache here — once per observer, outside the epoch loop — avoids
            // repeating the trig conversion for every epoch.
            let fixed_cache = match ObserverFixedCache::try_from(&obs_req.observer) {
                Ok(c) => c,
                Err(err) => {
                    // If the observer geometry is invalid, mark every epoch of
                    // this observer as failed and move on to the next observer.
                    for epoch in obs_req.mode.epochs() {
                        result.push(
                            epoch,
                            obs_req.observer.clone(),
                            Err(OutfitError::InvalidConversion(err.to_string())),
                        );
                    }
                    continue;
                }
            };

            for epoch in obs_req.mode.epochs() {
                let obs_time_mjd = epoch.to_mjd_tt_days();
                let value = O::compute_one(
                    &equi,
                    obs_time_mjd,
                    &obs_req.observer,
                    &fixed_cache,
                    jpl,
                    ut1,
                    &request.config,
                );
                result.push(epoch, obs_req.observer.clone(), value);
            }
        }

        result
    }

    /// Convert `self` to [`crate::EquinoctialElements`] for use in ephemeris
    /// computation.
    fn to_equinoctial_for_ephemeris(&self) -> Result<crate::EquinoctialElements, OutfitError> {
        self.to_equinoctial()
    }
}
