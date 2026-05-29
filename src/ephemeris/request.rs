//! Ephemeris request builder and output-kind trait.
//!
//! The central entry point for ephemeris generation is [`EphemerisRequest`],
//! a typed builder that collects any number of `(observer, mode)` pairs and
//! carries the desired output kind as a type parameter.
//!
//! # Output kinds
//!
//! Three zero-sized marker types select what is computed at each epoch:
//!
//! | Marker | Output type per epoch |
//! |---|---|
//! | [`Position`]  | [`ApparentPosition`] |
//! | [`Geometry`]  | [`BodyGeometry`] |
//! | [`Combined`]  | `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` |
//!
//! # Generation modes
//!
//! Each observer is paired with an [`EphemerisMode`] that describes *when* to
//! compute:
//!
//! | Variant | Epochs generated |
//! |---|---|
//! | [`EphemerisMode::Single`] | Exactly one epoch |
//! | [`EphemerisMode::Range`]  | Uniform grid from `start` to `end` |
//! | [`EphemerisMode::At`]     | Arbitrary list of epochs |
//!
//! Different observers may use different modes in the same request.
//!
//! # Example
//!
//! ```rust,ignore
//! use outfit::{
//!     EphemerisConfig, EphemerisRequest, EphemerisMode,
//!     ephemeris::request::Combined,
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
//!         .add(observer_b, EphemerisMode::At(vec![t1, t2, t3])),
//!     &jpl,
//!     &ut1,
//! );
//! ```

use hifitime::{Duration, Epoch};
use photom::observer::Observer;
use std::marker::PhantomData;

use super::{ApparentPosition, BodyGeometry, EphemerisConfig};
use crate::{EquinoctialElements, JPLEphem, OutfitError};
use hifitime::ut1::Ut1Provider;

// ---------------------------------------------------------------------------
// Output-kind sealed trait
// ---------------------------------------------------------------------------

mod sealed {
    pub trait Sealed {}
}

/// Trait implemented by [`Position`], [`Geometry`], and [`Combined`].
///
/// It is sealed so that no external crate can implement it.  Use one of the
/// three provided marker types to select what an [`EphemerisRequest`] computes.
pub trait EphemerisOutputKind: sealed::Sealed {
    /// The value produced for each epoch.
    type Output;

    /// Compute the output for a single `(epoch, observer)` pair.
    #[doc(hidden)]
    fn compute_one(
        equi: &EquinoctialElements,
        obs_time_mjd: f64,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<Self::Output, OutfitError>;
}

// ---------------------------------------------------------------------------
// Marker types
// ---------------------------------------------------------------------------

/// Output-kind marker: compute only the apparent sky position.
///
/// Use this when you need `(RA, Dec)` and distances but not the geometric
/// quantities (phase angle, elongation, …).
///
/// Produces [`ApparentPosition`] per epoch.
pub struct Position;

/// Output-kind marker: compute only the geometric quantities.
///
/// Use this when you need the phase angle, solar elongation, radial velocity
/// and apparent angular rates but not the full sky-coordinate conversion.
///
/// Produces [`BodyGeometry`] per epoch.
pub struct Geometry;

/// Output-kind marker: compute both apparent position and geometric quantities
/// from a **single orbit propagation** per epoch.
///
/// More efficient than computing [`Position`] and [`Geometry`] separately.
///
/// Produces `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` per epoch.
pub struct Combined;

impl sealed::Sealed for Position {}
impl sealed::Sealed for Geometry {}
impl sealed::Sealed for Combined {}

impl EphemerisOutputKind for Position {
    type Output = ApparentPosition;

    fn compute_one(
        equi: &EquinoctialElements,
        obs_time_mjd: f64,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<Self::Output, OutfitError> {
        super::apparent_position::compute(equi, obs_time_mjd, observer, jpl, ut1, config)
    }
}

impl EphemerisOutputKind for Geometry {
    type Output = BodyGeometry;

    fn compute_one(
        equi: &EquinoctialElements,
        obs_time_mjd: f64,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<Self::Output, OutfitError> {
        let state =
            super::apparent_position::propagate(equi, obs_time_mjd, observer, jpl, ut1, config)?;
        super::geometry::compute_geometry(&state, equi, &config.aberration)
    }
}

impl EphemerisOutputKind for Combined {
    type Output = (ApparentPosition, BodyGeometry);

    fn compute_one(
        equi: &EquinoctialElements,
        obs_time_mjd: f64,
        observer: &Observer,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
        config: &EphemerisConfig,
    ) -> Result<Self::Output, OutfitError> {
        super::apparent_position::compute_with_geometry(
            equi,
            obs_time_mjd,
            observer,
            jpl,
            ut1,
            config,
        )
    }
}

// ---------------------------------------------------------------------------
// EphemerisMode
// ---------------------------------------------------------------------------

/// Describes *when* to compute ephemerides for one observer.
///
/// Different observers in the same [`EphemerisRequest`] may use different
/// modes.
#[derive(Debug, Clone)]
pub enum EphemerisMode {
    /// Compute at exactly one epoch.
    Single(Epoch),

    /// Compute over a uniformly-spaced time range.
    ///
    /// Generates epochs `start, start + step, start + 2·step, …` for as long
    /// as the current epoch is ≤ `end`.  If `step ≤ 0` or `start > end` the
    /// observer produces no entries.
    Range {
        /// First epoch (inclusive).
        start: Epoch,
        /// Last epoch (inclusive if reachable by an exact multiple of `step`).
        end: Epoch,
        /// Time step between consecutive epochs.
        step: Duration,
    },

    /// Compute at an arbitrary list of epochs.
    At(Vec<Epoch>),
}

impl EphemerisMode {
    /// Expand this mode into the concrete list of epochs it covers.
    pub(crate) fn epochs(&self) -> Vec<Epoch> {
        match self {
            EphemerisMode::Single(t) => vec![*t],
            EphemerisMode::Range { start, end, step } => {
                if *step <= Duration::ZERO || start > end {
                    return vec![];
                }
                let approx_n =
                    (((*end - *start).to_seconds()) / step.to_seconds()).ceil() as usize + 1;
                let mut epochs = Vec::with_capacity(approx_n);
                let mut current = *start;
                while current <= *end {
                    epochs.push(current);
                    current += *step;
                }
                epochs
            }
            EphemerisMode::At(times) => times.clone(),
        }
    }
}

// ---------------------------------------------------------------------------
// ObserverRequest
// ---------------------------------------------------------------------------

/// A single `(observer, mode)` slot inside an [`EphemerisRequest`].
#[derive(Debug, Clone)]
pub struct ObserverRequest {
    /// The observing site.
    pub observer: Observer,
    /// When to compute ephemerides for this observer.
    pub mode: EphemerisMode,
}

// ---------------------------------------------------------------------------
// EphemerisRequest
// ---------------------------------------------------------------------------

/// A typed ephemeris generation request.
///
/// Collects any number of `(`[`Observer`]`, `[`EphemerisMode`]`)` pairs and
/// carries the desired output kind `O` as a compile-time type parameter.
///
/// Build the request with [`new`](Self::new) and chain [`add`](Self::add)
/// calls, then pass it to [`crate::OrbitalElements::compute`].
///
/// # Type parameter
///
/// `O` must be one of [`Position`], [`Geometry`], or [`Combined`].
///
/// # Example
///
/// ```rust,ignore
/// let request = EphemerisRequest::<Combined>::new(EphemerisConfig::default())
///     .add(observer_a, EphemerisMode::Range { start, end, step })
///     .add(observer_b, EphemerisMode::At(vec![t1, t2, t3]))
///     .add(observer_c, EphemerisMode::Single(t));
///
/// let result = elements.compute(&request, &jpl, &ut1);
/// ```
#[derive(Debug, Clone)]
pub struct EphemerisRequest<O: EphemerisOutputKind> {
    pub(crate) observers: Vec<ObserverRequest>,
    pub(crate) config: EphemerisConfig,
    _output: PhantomData<O>,
}

impl<O: EphemerisOutputKind> EphemerisRequest<O> {
    /// Create a new empty request with the given [`EphemerisConfig`].
    pub fn new(config: EphemerisConfig) -> Self {
        Self {
            observers: Vec::new(),
            config,
            _output: PhantomData,
        }
    }

    /// Add an `(observer, mode)` pair to the request.
    ///
    /// Returns `self` to allow builder-style chaining.
    pub fn add(mut self, observer: Observer, mode: EphemerisMode) -> Self {
        self.observers.push(ObserverRequest { observer, mode });
        self
    }

    /// Number of `(observer, mode)` pairs currently in the request.
    pub fn len(&self) -> usize {
        self.observers.len()
    }

    /// Returns `true` if no observers have been added yet.
    pub fn is_empty(&self) -> bool {
        self.observers.is_empty()
    }
}
