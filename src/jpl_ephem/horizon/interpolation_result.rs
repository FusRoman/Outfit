//! Interpolation result type for Horizons state evaluation.
//!
//! Overview
//! -----------------
//! `InterpResult` is the unified container returned by interpolation routines
//! in the Horizons backend. It always carries a **position** vector, and
//! may optionally include **velocity** and **acceleration** depending on the
//! caller’s request.
//!
//! Units
//! -----------------
//! * `position`: kilometers (km)  
//! * `velocity`: kilometers per day (km/day)  
//! * `acceleration`: kilometers per day² (km/day²)  
//! Use [`InterpResult::to_au`](crate::jpl_ephem::horizon::interpolation_result::InterpResult::to_au) to convert all present fields to AU-based units
//! (AU, AU/day, AU/day²).
//!
//! Arithmetic semantics
//! -----------------
//! Addition/subtraction are defined component-wise. Optional fields propagate
//! **only when present on both operands**; otherwise they are dropped (remain
//! `None`) to avoid silently mixing partial information.
//!
//! Typical workflow
//! -----------------
//! 1. Query ephemerides via
//!    [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris).
//! 2. Optionally convert to AU with [`InterpResult::to_au`](crate::jpl_ephem::horizon::interpolation_result::InterpResult::to_au).
//! 3. Combine results with `+`/`-` if needed (e.g., frame changes).
//!
//! See also
//! -----------------
//! * [`crate::jpl_ephem::horizon::horizon_records::HorizonRecord::interpolate`]
//!   – low-level interpolator producing an `InterpResult`.
//! * [`crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris`]
//!   – high-level query returning an `InterpResult`.

use crate::constants::AU;
use nalgebra::Vector3;
use std::ops::{Add, Div, Sub};

/// Interpolation result for a celestial body state.
///
/// Holds position and, optionally, velocity and acceleration for a given
/// evaluation time. Returned by high-level queries and low-level record
/// interpolation.
///
/// Fields
/// -----------------
/// * `position` — Cartesian position (km).
/// * `velocity` — Optional Cartesian velocity (km/day).
/// * `acceleration` — Optional Cartesian acceleration (km/day²).
///
/// See also
/// -----------------
/// * [`Self::to_au`] – convert present fields to AU units.
/// * [`crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris`] – high-level source.
/// * [`crate::jpl_ephem::horizon::horizon_records::HorizonRecord::interpolate`] – low-level source.
#[derive(Debug, PartialEq, Clone)]
pub struct InterpResult {
    pub position: Vector3<f64>,
    pub velocity: Option<Vector3<f64>>,
    pub acceleration: Option<Vector3<f64>>,
}

impl InterpResult {
    /// Convert the result to AU-based units.
    ///
    /// Scales all present vectors by `1 / AU`:
    /// * `position`: km → AU,
    /// * `velocity`: km/day → AU/day (if present),
    /// * `acceleration`: km/day² → AU/day² (if present).
    ///
    /// Return
    /// -----------------
    /// * A new `InterpResult` expressed in AU units.
    ///
    /// See also
    /// -----------------
    /// * [`AU`] – astronomical unit constant used for scaling.
    #[must_use = "`.to_au()` returns a new InterpResult; assign or use it"]
    pub fn to_au(&self) -> Self {
        self / AU
    }
}

impl Add for InterpResult {
    type Output = Self;

    /// Component-wise addition of two `InterpResult`s.
    ///
    /// Semantics
    /// -----------------
    /// * `position` is always added.
    /// * `velocity`/`acceleration` are added **only if present on both**;
    ///   otherwise they remain `None` in the result.
    ///
    /// Return
    /// -----------------
    /// * A new `InterpResult` combining both operands.
    ///
    /// Examples
    /// -----------------
    /// ```rust
    /// use nalgebra::Vector3;
    /// use outfit::jpl_ephem::horizon::interpolation_result::InterpResult;
    ///
    /// let a = InterpResult { position: Vector3::new(1.0, 0.0, 0.0), velocity: None, acceleration: None };
    /// let b = InterpResult { position: Vector3::new(2.0, 0.0, 0.0), velocity: None, acceleration: None };
    /// let c = a + b;
    /// assert_eq!(c.position.x, 3.0);
    /// assert!(c.velocity.is_none());
    /// ```
    fn add(self, other: Self) -> Self::Output {
        InterpResult {
            position: self.position + other.position,
            velocity: match (self.velocity, other.velocity) {
                (Some(v1), Some(v2)) => Some(v1 + v2),
                _ => None,
            },
            acceleration: match (self.acceleration, other.acceleration) {
                (Some(a1), Some(a2)) => Some(a1 + a2),
                _ => None,
            },
        }
    }
}

impl Add for &InterpResult {
    type Output = InterpResult;

    /// Component-wise addition for borrowed operands.
    ///
    /// Semantics identical to [`InterpResult::add`], but does not consume the inputs.
    fn add(self, other: Self) -> Self::Output {
        InterpResult {
            position: self.position + other.position,
            velocity: match (self.velocity, other.velocity) {
                (Some(v1), Some(v2)) => Some(v1 + v2),
                _ => None,
            },
            acceleration: match (self.acceleration, other.acceleration) {
                (Some(a1), Some(a2)) => Some(a1 + a2),
                _ => None,
            },
        }
    }
}

impl Sub for InterpResult {
    type Output = Self;

    /// Component-wise subtraction of two `InterpResult`s.
    ///
    /// Semantics
    /// -----------------
    /// * `position` is always subtracted.
    /// * `velocity`/`acceleration` are subtracted **only if present on both**;
    ///   otherwise they remain `None` in the result.
    fn sub(self, other: Self) -> Self::Output {
        InterpResult {
            position: self.position - other.position,
            velocity: match (self.velocity, other.velocity) {
                (Some(v1), Some(v2)) => Some(v1 - v2),
                _ => None,
            },
            acceleration: match (self.acceleration, other.acceleration) {
                (Some(a1), Some(a2)) => Some(a1 - a2),
                _ => None,
            },
        }
    }
}

impl Sub for &InterpResult {
    type Output = InterpResult;

    /// Component-wise subtraction for borrowed operands.
    ///
    /// Semantics identical to [`InterpResult::sub`], but does not consume the inputs.
    fn sub(self, other: Self) -> Self::Output {
        InterpResult {
            position: self.position - other.position,
            velocity: match (self.velocity, other.velocity) {
                (Some(v1), Some(v2)) => Some(v1 - v2),
                _ => None,
            },
            acceleration: match (self.acceleration, other.acceleration) {
                (Some(a1), Some(a2)) => Some(a1 - a2),
                _ => None,
            },
        }
    }
}

impl Div<f64> for InterpResult {
    type Output = Self;

    /// Divide all present fields by a scalar.
    ///
    /// Semantics
    /// -----------------
    /// * Scales `position`, and scales `velocity`/`acceleration` if present.
    /// * Leaves absent optional fields as `None`.
    fn div(self, rhs: f64) -> Self::Output {
        InterpResult {
            position: self.position / rhs,
            velocity: self.velocity.map(|v| v / rhs),
            acceleration: self.acceleration.map(|a| a / rhs),
        }
    }
}

impl Div<f64> for &InterpResult {
    type Output = InterpResult;

    /// Divide all present fields by a scalar (borrowed).
    ///
    /// Semantics identical to [`InterpResult::div`], but does not consume `self`.
    fn div(self, rhs: f64) -> Self::Output {
        InterpResult {
            position: self.position / rhs,
            velocity: self.velocity.map(|v| v / rhs),
            acceleration: self.acceleration.map(|a| a / rhs),
        }
    }
}
