use std::ops::{Add, Div, Sub};

use nalgebra::Vector3;

use crate::constants::AU;

/// Interplation result for celestial objects.
/// This struct contains the interpolated position, velocity, and acceleration
/// of a celestial object at a given time.
/// The position is always present, while velocity and acceleration are
/// optional, depending on the user request.
#[derive(Debug, PartialEq, Clone)]
pub struct InterpResult {
    pub position: Vector3<f64>,
    pub velocity: Option<Vector3<f64>>,
    pub acceleration: Option<Vector3<f64>>,
}

impl InterpResult {
    pub fn to_au(&self) -> Self {
        self / AU
    }
}

impl Add for InterpResult {
    type Output = Self;

    /// Adds two InterpResult instances together.
    /// The position is added directly, while velocity and acceleration are
    /// added only if both instances have them.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
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

impl<'a> Add for &'a InterpResult {
    type Output = InterpResult;

    /// Adds a reference to an InterpResult instance to another InterpResult.
    /// The position is added directly, while velocity and acceleration are
    /// added only if both instances have them.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
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

    /// Subtracts two InterpResult instances.
    /// The position is subtracted directly, while velocity and acceleration are
    /// subtracted only if both instances have them.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
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

impl<'a> Sub for &'a InterpResult {
    type Output = InterpResult;

    /// Subtracts a reference to an InterpResult instance from another InterpResult.
    /// The position is subtracted directly, while velocity and acceleration are
    /// subtracted only if both instances have them.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
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

    /// Divides the InterpResult instance by a scalar value.
    /// The position is divided directly, while velocity and acceleration are
    /// divided only if they are present.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
    fn div(self, rhs: f64) -> Self::Output {
        InterpResult {
            position: self.position.div(rhs),
            velocity: self.velocity.map(|v| v.div(rhs)),
            acceleration: self.acceleration.map(|a| a.div(rhs)),
        }
    }
}

impl<'a> Div<f64> for &'a InterpResult {
    type Output = InterpResult;

    /// Divides a reference to an InterpResult instance by a scalar value.
    /// The position is divided directly, while velocity and acceleration are
    /// divided only if they are present.
    /// If either instance does not have velocity or acceleration, the result
    /// will not have them either.
    fn div(self, rhs: f64) -> Self::Output {
        InterpResult {
            position: self.position.div(rhs),
            velocity: self.velocity.map(|v| v.div(rhs)),
            acceleration: self.acceleration.map(|a| a.div(rhs)),
        }
    }
}
