//! Public façade for ephemeris computation.
//!
//! This module exposes [`ApparentPosition`], [`BodyGeometry`],
//! [`EphemerisTable`], [`AberrationOrder`], and [`EphemerisConfig`], together
//! with the three concrete table aliases [`ApparentPositionTable`],
//! [`GeometryTable`], and [`ApparentPositionAndGeometryTable`].
//!
//! Eight methods on [`OrbitalElements`] cover all combinations of single-epoch
//! vs. bulk and position-only vs. geometry-only vs. combined:
//!
//! | Method | Returns | Propagations |
//! |---|---|---|
//! | [`apparent_position`](OrbitalElements::apparent_position) | [`ApparentPosition`] | 1 |
//! | [`body_geometry`](OrbitalElements::body_geometry) | [`BodyGeometry`] | 1 |
//! | [`apparent_position_and_geometry`](OrbitalElements::apparent_position_and_geometry) | `(`[`ApparentPosition`]`,` [`BodyGeometry`]`)` | 1 |
//! | [`apparent_position_range`](OrbitalElements::apparent_position_range) | [`ApparentPositionTable`] | 1 / epoch |
//! | [`apparent_position_at`](OrbitalElements::apparent_position_at) | [`ApparentPositionTable`] | 1 / epoch |
//! | [`body_geometry_range`](OrbitalElements::body_geometry_range) | [`GeometryTable`] | 1 / epoch |
//! | [`body_geometry_at`](OrbitalElements::body_geometry_at) | [`GeometryTable`] | 1 / epoch |
//! | [`apparent_position_and_geometry_range`](OrbitalElements::apparent_position_and_geometry_range) | [`ApparentPositionAndGeometryTable`] | 1 / epoch |
//! | [`apparent_position_and_geometry_at`](OrbitalElements::apparent_position_and_geometry_at) | [`ApparentPositionAndGeometryTable`] | 1 / epoch |
//!
//! # Time input
//!
//! All entry points accept any type that implements `Into<`[`hifitime::Epoch`]`>`:
//!
//! ```rust,ignore
//! use hifitime::{Epoch, Duration};
//! use core::str::FromStr;
//!
//! // From an Epoch directly
//! elements.apparent_position(epoch, &observer, &jpl, &ut1, &config)?;
//!
//! // From an explicit MJD TT value
//! elements.apparent_position(Epoch::from_mjd_tt(60310.5), &observer, &jpl, &ut1, &config)?;
//!
//! // From an ISO 8601 / RFC 3339 string (parse error handled by the caller)
//! let t = Epoch::from_str("2025-01-01T00:00:00 UTC")?;
//! elements.apparent_position(t, &observer, &jpl, &ut1, &config)?;
//! ```
//!
//! Internally, epochs are converted to MJD TT before entering the pipeline.
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
pub(crate) mod geometry;
pub(crate) mod observation_ephemeris;
mod table;

pub use aberration::AberrationOrder;
pub use apparent_position::ApparentPosition;
pub use geometry::BodyGeometry;
pub use table::{
    ApparentPositionAndGeometryTable, ApparentPositionTable, EphemerisTable, GeometryTable,
};

use hifitime::{ut1::Ut1Provider, Duration, Epoch};
use photom::observer::Observer;

use crate::{propagator::PropagatorKind, JPLEphem, OrbitalElements, OutfitError};

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
// OrbitalElements — ephemeris entry points
// ---------------------------------------------------------------------------

impl OrbitalElements {
    /// Compute the apparent equatorial position of the body at a single epoch.
    ///
    /// This is the primary user-facing entry point for single-epoch ephemeris
    /// computation.  Internally, it converts `self` to
    /// [`crate::EquinoctialElements`], propagates the orbit to `obs_time`,
    /// resolves the observer's heliocentric position, applies the aberration
    /// correction, and converts the resulting line-of-sight vector to
    /// equatorial coordinates.
    ///
    /// If you also need geometric quantities (phase angle, elongation, …) at
    /// the same epoch, use [`apparent_position_and_geometry`](Self::apparent_position_and_geometry)
    /// instead — it runs the same single propagation and returns both.
    ///
    /// See the [module-level documentation](self) for the full pipeline diagram
    /// and time-input examples.
    ///
    /// # Arguments
    ///
    /// - `obs_time`  – Observation epoch.  Accepts any type implementing
    ///   `Into<`[`hifitime::Epoch`]`>`.
    /// - `observer`  – Observing site (geodetic coordinates or parallax
    ///   constants).
    /// - `jpl`       – JPL ephemeris providing Earth's heliocentric state
    ///   and (for N-body) planetary states.
    /// - `ut1`       – UT1 provider for Greenwich Sidereal Time and observer
    ///   geocentric position and velocity.
    /// - `config`    – Propagator selection and aberration correction order.
    ///
    /// # Returns
    ///
    /// An [`ApparentPosition`] containing:
    ///
    /// - `coord` – predicted equatorial coordinates (RA, Dec in radians);
    ///   error fields are `0.0` (prediction, not measurement).
    /// - `geocentric_dist` – Earth-centre to body distance \[AU\].
    /// - `heliocentric_dist` – Sun to body distance \[AU\].
    ///
    /// # Errors
    ///
    /// Returns [`OutfitError`] if:
    /// - The orbit is hyperbolic or parabolic (eccentricity ≥ 1).
    /// - Orbit propagation fails.
    /// - The JPL ephemeris data is unavailable for the requested epoch.
    /// - `self` cannot be converted to equinoctial elements.
    pub fn apparent_position(
        &self,
        obs_time: impl Into<Epoch>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<ApparentPosition, OutfitError> {
        let obs_time_mjd = obs_time.into().to_mjd_tt_days();
        let equi = self.to_equinoctial_for_ephemeris()?;
        apparent_position::compute(&equi, obs_time_mjd, observer, jpl, ut1, config)
    }

    /// Compute the geometric quantities of the body at a single epoch.
    ///
    /// Returns a [`BodyGeometry`] containing phase angle, solar elongation,
    /// radial velocity and apparent angular rates — without computing the full
    /// sky-coordinate conversion in [`ApparentPosition`].
    ///
    /// If you also need `(RA, Dec)` at the same epoch, use
    /// [`apparent_position_and_geometry`](Self::apparent_position_and_geometry)
    /// to avoid a second propagation.
    ///
    /// # Arguments
    ///
    /// Same as [`apparent_position`](Self::apparent_position).
    ///
    /// # Errors
    ///
    /// Same conditions as [`apparent_position`](Self::apparent_position).
    pub fn body_geometry(
        &self,
        obs_time: impl Into<Epoch>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<BodyGeometry, OutfitError> {
        let obs_time_mjd = obs_time.into().to_mjd_tt_days();
        let equi = self.to_equinoctial_for_ephemeris()?;
        let state = apparent_position::propagate(&equi, obs_time_mjd, observer, jpl, ut1, config)?;
        geometry::compute_geometry(&state, &equi, &config.aberration)
    }

    /// Compute both the apparent position and geometric quantities at a single
    /// epoch using a **single orbit propagation**.
    ///
    /// Equivalent to calling [`apparent_position`](Self::apparent_position) and
    /// [`body_geometry`](Self::body_geometry) separately, but avoids propagating
    /// the orbit twice.  Prefer this method when both types of output are needed.
    ///
    /// # Arguments
    ///
    /// Same as [`apparent_position`](Self::apparent_position).
    ///
    /// # Returns
    ///
    /// `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)`.
    ///
    /// # Errors
    ///
    /// Same conditions as [`apparent_position`](Self::apparent_position).
    pub fn apparent_position_and_geometry(
        &self,
        obs_time: impl Into<Epoch>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<(ApparentPosition, BodyGeometry), OutfitError> {
        let obs_time_mjd = obs_time.into().to_mjd_tt_days();
        let equi = self.to_equinoctial_for_ephemeris()?;
        apparent_position::compute_with_geometry(&equi, obs_time_mjd, observer, jpl, ut1, config)
    }

    /// Compute apparent positions over a uniformly-spaced time range.
    ///
    /// Generates epochs `start, start + step, start + 2·step, …` for as long
    /// as the current epoch is ≤ `end`, then computes the apparent position at
    /// each one.  Errors are collected per-epoch — a failure at one instant
    /// does not abort the rest of the computation.
    ///
    /// If `step` is zero or negative, or if `start > end`, the returned table
    /// is empty.
    ///
    /// # Arguments
    ///
    /// - `start`    – First observation epoch (inclusive).
    /// - `end`      – Last observation epoch (inclusive).
    /// - `step`     – Time step between consecutive epochs.
    /// - `observer` – Observing site.
    /// - `jpl`      – JPL ephemeris.
    /// - `ut1`      – UT1 provider.
    /// - `config`   – Computation configuration.
    ///
    /// # Returns
    ///
    /// An [`ApparentPositionTable`] whose entries are in chronological order.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::{Epoch, Duration};
    /// use core::str::FromStr;
    ///
    /// let start = Epoch::from_str("2025-01-01T00:00:00 UTC").unwrap();
    /// let end   = Epoch::from_str("2025-01-31T00:00:00 UTC").unwrap();
    /// let step  = Duration::from_days(1.0);
    ///
    /// let table = elements.apparent_position_range(
    ///     start, end, step, &observer, &jpl, &ut1, &config,
    /// );
    /// for (epoch, result) in table {
    ///     match result {
    ///         Ok(pos) => println!("{epoch}: RA = {:.6} rad", pos.coord.ra),
    ///         Err(e)  => eprintln!("{epoch}: error — {e}"),
    ///     }
    /// }
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn apparent_position_range(
        &self,
        start: impl Into<Epoch>,
        end: impl Into<Epoch>,
        step: Duration,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> ApparentPositionTable {
        let start = start.into();
        let end = end.into();

        if step <= Duration::ZERO || start > end {
            return ApparentPositionTable::with_capacity(0);
        }

        let approx_n = ((end - start).to_seconds() / step.to_seconds()).ceil() as usize + 1;
        let mut table = ApparentPositionTable::with_capacity(approx_n);

        let mut current = start;
        while current <= end {
            let result = self.apparent_position(current, observer, jpl, ut1, config);
            table.push(current, result);
            current += step;
        }

        table
    }

    /// Compute apparent positions at an arbitrary set of epochs.
    ///
    /// Accepts any value implementing `IntoIterator` whose items convert to
    /// [`hifitime::Epoch`] — for example a `Vec<Epoch>`, a slice `&[Epoch]`,
    /// or any lazy iterator of epochs.  Errors are collected per-epoch.
    ///
    /// # Arguments
    ///
    /// - `times`    – Any iterable of epochs (or values convertible to
    ///   [`hifitime::Epoch`]).
    /// - `observer` – Observing site.
    /// - `jpl`      – JPL ephemeris.
    /// - `ut1`      – UT1 provider.
    /// - `config`   – Computation configuration.
    ///
    /// # Returns
    ///
    /// An [`ApparentPositionTable`] in the same order as the input iterator.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::Epoch;
    /// use core::str::FromStr;
    ///
    /// let times = vec![
    ///     Epoch::from_str("2025-01-01T00:00:00 UTC").unwrap(),
    ///     Epoch::from_str("2025-06-15T12:00:00 UTC").unwrap(),
    ///     Epoch::from_mjd_tt(60500.0),
    /// ];
    ///
    /// let table = elements.apparent_position_at(
    ///     times, &observer, &jpl, &ut1, &config,
    /// );
    /// println!("{} successes, {} errors", table.success_count(), table.error_count());
    /// ```
    pub fn apparent_position_at(
        &self,
        times: impl IntoIterator<Item: Into<Epoch>>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> ApparentPositionTable {
        let iter = times.into_iter();
        let (lower, _) = iter.size_hint();
        let mut table = ApparentPositionTable::with_capacity(lower);

        for t in iter {
            let epoch = t.into();
            let result = self.apparent_position(epoch, observer, jpl, ut1, config);
            table.push(epoch, result);
        }

        table
    }

    /// Compute geometric quantities over a uniformly-spaced time range.
    ///
    /// Same epoch-generation logic as [`apparent_position_range`](Self::apparent_position_range),
    /// but computes only [`BodyGeometry`] at each epoch (phase angle, solar
    /// elongation, radial velocity, apparent angular rates) without performing
    /// the full sky-coordinate conversion.
    ///
    /// If you also need `(RA, Dec)` for the same time range, use
    /// [`apparent_position_and_geometry_range`](Self::apparent_position_and_geometry_range)
    /// instead — it avoids a second propagation.
    ///
    /// If `step` is zero or negative, or if `start > end`, the returned table
    /// is empty.
    ///
    /// # Returns
    ///
    /// A [`GeometryTable`] whose entries are in chronological order.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::{Epoch, Duration};
    /// use core::str::FromStr;
    ///
    /// let start = Epoch::from_str("2025-01-01T00:00:00 UTC").unwrap();
    /// let end   = Epoch::from_str("2025-01-31T00:00:00 UTC").unwrap();
    /// let step  = Duration::from_days(1.0);
    ///
    /// let table = elements.body_geometry_range(
    ///     start, end, step, &observer, &jpl, &ut1, &config,
    /// );
    /// for (epoch, result) in table {
    ///     if let Ok(geo) = result {
    ///         println!("{epoch}: phase = {:.4} rad", geo.phase_angle);
    ///     }
    /// }
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn body_geometry_range(
        &self,
        start: impl Into<Epoch>,
        end: impl Into<Epoch>,
        step: Duration,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> GeometryTable {
        let start = start.into();
        let end = end.into();

        if step <= Duration::ZERO || start > end {
            return GeometryTable::with_capacity(0);
        }

        let approx_n = ((end - start).to_seconds() / step.to_seconds()).ceil() as usize + 1;
        let mut table = GeometryTable::with_capacity(approx_n);

        let mut current = start;
        while current <= end {
            let result = self.body_geometry(current, observer, jpl, ut1, config);
            table.push(current, result);
            current += step;
        }

        table
    }

    /// Compute geometric quantities at an arbitrary set of epochs.
    ///
    /// Same as [`body_geometry_range`](Self::body_geometry_range) but accepts
    /// any iterable of epochs rather than a uniform time range.  Errors are
    /// collected per-epoch.
    ///
    /// If you also need `(RA, Dec)` for the same epochs, use
    /// [`apparent_position_and_geometry_at`](Self::apparent_position_and_geometry_at)
    /// instead — it avoids a second propagation.
    ///
    /// # Returns
    ///
    /// A [`GeometryTable`] in the same order as the input iterator.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::Epoch;
    ///
    /// let times = vec![Epoch::from_mjd_tt(60310.0), Epoch::from_mjd_tt(60320.0)];
    /// let table = elements.body_geometry_at(times, &observer, &jpl, &ut1, &config);
    /// println!("{} successes", table.success_count());
    /// ```
    pub fn body_geometry_at(
        &self,
        times: impl IntoIterator<Item: Into<Epoch>>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> GeometryTable {
        let iter = times.into_iter();
        let (lower, _) = iter.size_hint();
        let mut table = GeometryTable::with_capacity(lower);

        for t in iter {
            let epoch = t.into();
            let result = self.body_geometry(epoch, observer, jpl, ut1, config);
            table.push(epoch, result);
        }

        table
    }

    /// Compute both apparent position and geometric quantities over a
    /// uniformly-spaced time range, using a **single propagation per epoch**.
    ///
    /// Equivalent to calling [`apparent_position_range`](Self::apparent_position_range)
    /// and [`body_geometry_range`](Self::body_geometry_range) separately, but
    /// avoids propagating the orbit twice at each epoch.  Prefer this method
    /// when both output types are needed for a time range.
    ///
    /// If `step` is zero or negative, or if `start > end`, the returned table
    /// is empty.
    ///
    /// # Returns
    ///
    /// An [`ApparentPositionAndGeometryTable`] whose entries are in
    /// chronological order.  Each successful entry is a
    /// `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` tuple.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::{Epoch, Duration};
    /// use core::str::FromStr;
    ///
    /// let start = Epoch::from_str("2025-01-01T00:00:00 UTC").unwrap();
    /// let end   = Epoch::from_str("2025-01-31T00:00:00 UTC").unwrap();
    /// let step  = Duration::from_days(1.0);
    ///
    /// let table = elements.apparent_position_and_geometry_range(
    ///     start, end, step, &observer, &jpl, &ut1, &config,
    /// );
    /// for (epoch, result) in table {
    ///     if let Ok((pos, geo)) = result {
    ///         println!(
    ///             "{epoch}: RA = {:.4} rad, phase = {:.4} rad",
    ///             pos.coord.ra, geo.phase_angle,
    ///         );
    ///     }
    /// }
    /// ```
    #[allow(clippy::too_many_arguments)]
    pub fn apparent_position_and_geometry_range(
        &self,
        start: impl Into<Epoch>,
        end: impl Into<Epoch>,
        step: Duration,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> ApparentPositionAndGeometryTable {
        let start = start.into();
        let end = end.into();

        if step <= Duration::ZERO || start > end {
            return ApparentPositionAndGeometryTable::with_capacity(0);
        }

        let approx_n = ((end - start).to_seconds() / step.to_seconds()).ceil() as usize + 1;
        let mut table = ApparentPositionAndGeometryTable::with_capacity(approx_n);

        let mut current = start;
        while current <= end {
            let result = self.apparent_position_and_geometry(current, observer, jpl, ut1, config);
            table.push(current, result);
            current += step;
        }

        table
    }

    /// Compute both apparent position and geometric quantities at an arbitrary
    /// set of epochs, using a **single propagation per epoch**.
    ///
    /// Same as [`apparent_position_and_geometry_range`](Self::apparent_position_and_geometry_range)
    /// but accepts any iterable of epochs rather than a uniform time range.
    /// Errors are collected per-epoch.
    ///
    /// # Returns
    ///
    /// An [`ApparentPositionAndGeometryTable`] in the same order as the input
    /// iterator.  Each successful entry is a
    /// `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` tuple.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// use hifitime::Epoch;
    ///
    /// let times = vec![Epoch::from_mjd_tt(60310.0), Epoch::from_mjd_tt(60320.0)];
    /// let table = elements.apparent_position_and_geometry_at(
    ///     times, &observer, &jpl, &ut1, &config,
    /// );
    /// for (epoch, result) in table.successes() {
    ///     let (pos, geo) = result;
    ///     println!("{epoch}: elongation = {:.4} rad", geo.solar_elongation);
    /// }
    /// ```
    pub fn apparent_position_and_geometry_at(
        &self,
        times: impl IntoIterator<Item: Into<Epoch>>,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> ApparentPositionAndGeometryTable {
        let iter = times.into_iter();
        let (lower, _) = iter.size_hint();
        let mut table = ApparentPositionAndGeometryTable::with_capacity(lower);

        for t in iter {
            let epoch = t.into();
            let result = self.apparent_position_and_geometry(epoch, observer, jpl, ut1, config);
            table.push(epoch, result);
        }

        table
    }

    /// Convert `self` to [`crate::EquinoctialElements`] for use in ephemeris
    /// computation.
    ///
    /// Delegates to [`OrbitalElements::to_equinoctial`], which handles all
    /// supported element variants and returns an appropriate error for
    /// unsupported orbit types (e.g. cometary elements or hyperbolic orbits).
    fn to_equinoctial_for_ephemeris(&self) -> Result<crate::EquinoctialElements, OutfitError> {
        self.to_equinoctial()
    }
}
