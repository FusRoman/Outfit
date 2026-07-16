//! Classification of two-body orbital regimes from the energy parameter `alpha`.

/// Classifies the orbital regime based on the sign of the energy parameter `alpha`.
///
/// The parameter `alpha` is defined as:
///
/// ```text
/// alpha = 2 * specific_orbital_energy = 2 * (v²/2 - μ/r)
/// ```
///
/// Depending on its value:
/// - **Elliptic** (`alpha < 0`) – Closed orbit, bounded motion.
/// - **Parabolic** (`alpha = 0`) – Escape trajectory with zero excess velocity (marginally unbound).
/// - **Hyperbolic** (`alpha > 0`) – Open orbit, unbounded motion.
///
/// This classification is used to select the correct branch when solving the
/// universal Kepler equation and related initial approximations.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum OrbitType {
    /// Bound, closed orbit (alpha < 0).
    Elliptic,
    /// Open, unbound orbit (alpha > 0).
    Hyperbolic,
    /// Parabolic trajectory (alpha = 0).
    Parabolic,
}

impl OrbitType {
    /// Determine the [`OrbitType`](crate::kepler::OrbitType) from a given value of `alpha`.
    ///
    /// # Arguments
    /// * `alpha` – Twice the specific orbital energy.
    ///
    /// # Returns
    /// * `OrbitType::Elliptic` if `alpha < 0`
    /// * `OrbitType::Hyperbolic` if `alpha > 0`
    /// * `OrbitType::Parabolic` if `alpha == 0`
    pub fn from_alpha(alpha: f64) -> Self {
        match alpha {
            energy if energy < 0.0 => OrbitType::Elliptic,
            energy if energy > 0.0 => OrbitType::Hyperbolic,
            _ => OrbitType::Parabolic,
        }
    }
}
