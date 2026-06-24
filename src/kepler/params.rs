//! Input parameters bundle for the universal-variable Kepler solver.

use super::orbit_type::OrbitType;

/// Parameters required to solve the universal Kepler equation.
///
/// This struct bundles together the state and constants needed to propagate an
/// orbit using the universal variable formulation:
///
/// * `dt`: Propagation time interval Δt = t - t₀ (in days).
/// * `r0`: Heliocentric distance at the reference epoch t₀ (in AU).
/// * `sig0`: Radial velocity component at t₀ (in AU/day).
/// * `mu`: Standard gravitational parameter μ = GM (in AU³/day²).
/// * `alpha`: Twice the specific orbital energy (2E).
/// * `e0`: Orbital eccentricity (unitless).
///
/// The associated method [`orbit_type`](UniversalKeplerParams::orbit_type)
/// classifies the orbit into elliptical, parabolic, or hyperbolic regimes
/// based on `alpha`.
#[derive(Debug, Clone, Copy)]
pub struct UniversalKeplerParams {
    /// Propagation time interval Δt = t - t₀, in days.
    pub dt: f64,
    /// Heliocentric distance at the reference epoch t₀, in AU.
    pub r0: f64,
    /// Radial velocity component at t₀ (= r₀·v₀ / r₀), in AU/day.
    pub sig0: f64,
    /// Standard gravitational parameter μ = GM, in AU³/day².
    pub mu: f64,
    /// Twice the specific orbital energy (alpha = 2E).
    pub alpha: f64,
    /// Orbital eccentricity (unitless).
    pub e0: f64,
}

impl UniversalKeplerParams {
    /// Returns the [`OrbitType`](crate::kepler::OrbitType) (elliptic, parabolic, or hyperbolic)
    /// corresponding to the value of `alpha`.
    pub fn orbit_type(&self) -> OrbitType {
        OrbitType::from_alpha(self.alpha)
    }
}
