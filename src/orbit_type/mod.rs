//! Orbital element representations and uncertainty propagation
//!
//! This module defines multiple **canonical orbital element sets** with full support for
//! **uncertainty propagation** during conversions between representations. It provides
//! three distinct parameterizations of Keplerian orbits, each with advantages for different
//! orbital regimes:
//!
//! - [`keplerian_element`](crate::orbit_type::keplerian_element) — Classical Keplerian elements $(a, e, i, \Omega, \omega, M)$,
//!   widely used for elliptic and hyperbolic orbits but singular for circular and equatorial cases.
//! - [`equinoctial_element`](crate::orbit_type::equinoctial_element) — Equinoctial elements $(a, h, k, p, q, \lambda)$,
//!   a **non-singular formulation** ideal for orbit determination and propagation near
//!   zero eccentricity or inclination.
//! - [`cometary_element`](crate::orbit_type::cometary_element) — Perihelion-based representation $(q, e, i, \Omega, \omega, \nu)$,
//!   natural for parabolic and hyperbolic orbits where perihelion distance $q$ is better
//!   defined than semi-major axis $a$.
//!
//! The [`OrbitalElements`] enum acts as a **type-erased wrapper** that can hold any of these
//! three representations along with optional **uncertainty** and **covariance** information,
//! providing uniform constructors and rigorous conversion methods with uncertainty propagation.
//!
//! ## Choosing an orbital element representation
//!
//! Each representation has optimal use cases:
//!
//! **Keplerian elements** $(a, e, i, \Omega, \omega, M)$:
//! - Intuitive physical interpretation: size, shape, orientation, and position in orbit
//! - Standard in classical celestial mechanics and literature
//! - Singular when $e \to 0$ (circular) or $i \to 0$ (equatorial): derivatives become undefined
//! - Best for: communication, visualization, moderate eccentricity/inclination orbits
//!
//! **Equinoctial elements** $(a, h, k, p, q, \lambda)$:
//! - Non-singular for all $e < 1$ and $0 \leq i < \pi$
//! - Smooth, well-defined derivatives enable robust orbit fitting and least-squares optimization
//! - Ideal for numerical propagation and uncertainty analysis
//! - Best for: orbit determination, propagation, covariance propagation, near-circular or near-equatorial orbits
//!
//! **Cometary elements** $(q, e, i, \Omega, \omega, \nu)$:
//! - Uses perihelion distance $q = a(1 - e)$ instead of semi-major axis
//! - Well-behaved for $e \geq 1$ (parabolic and hyperbolic orbits)
//! - Natural for cometary and hyperbolic trajectories
//! - Best for: comets, interstellar objects, hyperbolic encounters
//!
//! ## Uncertainty representation and propagation
//!
//! Orbital uncertainties arise from measurement errors, numerical approximations, and model limitations.
//! This module represents uncertainties in two complementary ways:
//!
//! 1. **Standard deviations**: individual $1\sigma$ uncertainties on each element
//!    (see [`KeplerianUncertainty`](crate::orbit_type::uncertainty::KeplerianUncertainty), [`EquinoctialUncertainty`](crate::orbit_type::uncertainty::EquinoctialUncertainty), [`CometaryUncertainty`](crate::orbit_type::uncertainty::CometaryUncertainty))
//!
//! 2. **Covariance matrices**: full $6 \times 6$ symmetric positive semi-definite matrices capturing
//!    element correlations (see [`OrbitalCovariance`](crate::orbit_type::uncertainty::OrbitalCovariance))
//!
//! When converting between representations, **covariance matrices are transformed** using the Jacobian
//! of the conversion:
//!
//! $$
//! \Sigma_y = J \, \Sigma_x \, J^\top
//! $$
//!
//! where $\Sigma_x$ is the covariance in the source representation, $\Sigma_y$ is the covariance in the
//! target representation, and $J = \partial \mathbf{y} / \partial \mathbf{x}$ is the Jacobian matrix of
//! partial derivatives evaluated at the nominal element values.
//!
//! This **linear covariance propagation** preserves the statistical properties of uncertainties under
//! first-order approximation, which is accurate for small uncertainties relative to the element values.
//!
//! ## Mathematical background: Jacobian matrices
//!
//! Each element representation provides Jacobian methods for conversions:
//!
//! ### Keplerian ↔ Equinoctial
//!
//! The transformation between Keplerian $(a, e, i, \Omega, \omega, M)$ and equinoctial
//! $(a, h, k, p, q, \lambda)$ is defined by:
//!
//! $$
//! \begin{aligned}
//! h &= e \sin(\Omega + \omega) \\
//! k &= e \cos(\Omega + \omega) \\
//! p &= \tan(i/2) \sin(\Omega) \\
//! q &= \tan(i/2) \cos(\Omega) \\
//! \lambda &= \Omega + \omega + M
//! \end{aligned}
//! $$
//!
//! The Jacobian $J_{\text{Kep} \to \text{Eq}}$ is computed analytically in
//! [`KeplerianElements::jacobian_to_equinoctial`], and the inverse in
//! [`EquinoctialElements::jacobian_to_keplerian`].
//!
//! **Singular cases**: When $e \approx 0$, the longitude of periapsis $\varpi = \Omega + \omega$
//! is undefined, and derivatives involving $\partial \varpi / \partial h$ and
//! $\partial \varpi / \partial k$ become singular. These are handled by setting the derivatives
//! to zero when $e < \epsilon$ (typically $\epsilon = 10^{-12}$). Similarly, when $i \approx 0$,
//! $\Omega$ is undefined and related derivatives are zeroed.
//!
//! ### Cometary ↔ Keplerian
//!
//! The cometary representation $(q, e, i, \Omega, \omega, \nu)$ relates to Keplerian elements through:
//!
//! $$
//! a = \frac{q}{1 - e}, \quad M = M(e, \nu)
//! $$
//!
//! where $M(e, \nu)$ is the mean anomaly computed from the true anomaly $\nu$ via the eccentric
//! anomaly $E$ (for $e < 1$) or hyperbolic anomaly $H$ (for $e > 1$). The Jacobian
//! $J_{\text{Com} \to \text{Kep}}$ is computed in [`CometaryElements::jacobian_to_keplerian`].
//!
//! For conversions to equinoctial, the **chain rule** is applied:
//!
//! $$
//! J_{\text{Com} \to \text{Eq}} = J_{\text{Kep} \to \text{Eq}} \cdot J_{\text{Com} \to \text{Kep}}
//! $$
//!
//! This is implemented in [`CometaryElements::jacobian_to_equinoctial`].
//!
//! ## Conversion methods with uncertainty propagation
//!
//! The [`OrbitalElements`] enum provides conversion methods that automatically propagate covariance:
//!
//! - [`OrbitalElements::to_keplerian`] — Convert to Keplerian, propagating covariance if present
//! - [`OrbitalElements::to_equinoctial`] — Convert to equinoctial, propagating covariance if present
//!
//! If a covariance matrix is attached to the source elements, it is transformed using the appropriate
//! Jacobian and attached to the result. The 1-σ uncertainties are then recomputed from the diagonal
//! of the transformed covariance.
//!
//! **Error handling**: Conversions that are mathematically undefined (e.g., parabolic cometary → Keplerian,
//! where $a = \infty$) return `Err(OutfitError::InvalidConversion)`.
//!
//! ## Units and conventions
//!
//! - **Lengths**: astronomical units (AU)
//! - **Angles**: radians
//! - **Time**: Modified Julian Date (MJD) in TDB scale for epochs
//! - **Velocities**: AU/day
//! - **Reference frame**: heliocentric ecliptic J2000
//!
//! ## Typical workflow
//!
//! ```rust, no_run
//! use nalgebra::Vector3;
//! use outfit::orbit_type::OrbitalElements;
//!
//! // State vector (heliocentric J2000)
//! let r = Vector3::new(1.0, 0.0, 0.0);
//! let v = Vector3::new(0.0, 1.0, 0.0);
//!
//! // Build canonical orbital elements from state
//! let elems = OrbitalElements::from_orbital_state(&r, &v, 2460000.5);
//!
//! // Convert to Keplerian form if possible (uncertainty propagated if present)
//! if let Ok(kep) = elems.to_keplerian() {
//!     if let OrbitalElements::Keplerian { elements, uncertainty, .. } = kep {
//!         println!("semi-major axis = {}", elements.semi_major_axis);
//!         if let Some(unc) = uncertainty {
//!             println!("  ± {} AU", unc.semi_major_axis);
//!         }
//!     }
//! }
//! ```
//!
//! ## See also
//!
//! - [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements) — Classical elements with Jacobian methods
//! - [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements) — Non-singular elements with propagation
//! - [`CometaryElements`](crate::orbit_type::cometary_element::CometaryElements) — Perihelion-based elements for hyperbolic orbits
//! - [`uncertainty`](crate::orbit_type::uncertainty) — Uncertainty structures and covariance propagation
//!
//! ## References
//!
//! - Milani & Gronchi, *Theory of Orbit Determination* (2010), Chapter 5
//! - Walker et al., "A Set of Modified Equinoctial Orbital Elements", *Celestial Mechanics* 36 (1985)
use nalgebra::Vector3;

use crate::{
    orb_elem::ccek1,
    orbit_type::{
        cometary_element::CometaryElements,
        equinoctial_element::EquinoctialElements,
        keplerian_element::KeplerianElements,
        uncertainty::{
            CometaryUncertainty, EquinoctialUncertainty, KeplerianUncertainty, OrbitalCovariance,
        },
    },
    OutfitError,
};

/// Equinoctial orbital elements and related conversions.
pub mod equinoctial_element;

/// Classical Keplerian elements structure and utilities.
pub mod keplerian_element;

/// Cometary (parabolic/hyperbolic) orbital elements and related conversions.
pub mod cometary_element;

/// Uncertainty structures for orbital elements.
pub mod uncertainty;

/// Canonical orbital elements in multiple representations.
///
/// This enum acts as a sum type over several orbital-element parameterizations.
/// It lets callers request or carry elements without committing to a single
/// representation at the type level.
///
/// Variants
/// --------
/// * `Keplerian`  — Classical elements `(a, e, i, Ω, ω, M)`; best for elliptic orbits.
/// * `Equinoctial` — Non-singular elements `(a, h, k, p, q, λ)`; robust near e≈0 and i≈0.
/// * `Cometary`   — Perihelion form `(q, e, i, Ω, ω, ν)`; convenient for e ≥ 1.
///
/// See also
/// --------
/// * [`KeplerianElements`] – Classical Keplerian elements.
/// * [`EquinoctialElements`] – Non-singular elements (equinoctial).
/// * [`CometaryElements`] – Perihelion distance representation for e≥1.
/// * [`crate::orbit_type::OrbitalElements::to_equinoctial`] – Lossless conversion for elliptic orbits.
/// * [`crate::orbit_type::OrbitalElements::to_keplerian`] – Conversion with domain checks.
#[derive(Debug, Clone, PartialEq)]
pub enum OrbitalElements {
    Keplerian {
        elements: KeplerianElements,
        uncertainty: Option<KeplerianUncertainty>,
        covariance: Option<OrbitalCovariance>,
    },
    Equinoctial {
        elements: EquinoctialElements,
        uncertainty: Option<EquinoctialUncertainty>,
        covariance: Option<OrbitalCovariance>,
    },
    Cometary {
        elements: CometaryElements,
        uncertainty: Option<CometaryUncertainty>,
        covariance: Option<OrbitalCovariance>,
    },
}

impl OrbitalElements {
    /// Build orbital elements from a heliocentric Cartesian state vector.
    ///
    /// This constructor converts a position–velocity pair `[r, v]` in the J2000
    /// equatorial frame into the appropriate orbital element representation.
    ///
    /// Depending on the reciprocal semi-major axis `1/a`, the method selects the
    /// correct formulation:
    ///
    /// * **Elliptical orbits** (`1/a > 0`): returns [`OrbitalElements::Keplerian`]
    ///   containing a [`KeplerianElements`] struct.  
    /// * **Parabolic or hyperbolic orbits** (`1/a ≤ 0`): returns
    ///   [`OrbitalElements::Cometary`] containing a [`CometaryElements`] struct.  
    ///
    /// Arguments
    /// -----------------
    /// * `position`: Position vector `[x, y, z]` in AU, heliocentric J2000.  
    /// * `velocity`: Velocity vector `[vx, vy, vz]` in AU/day, heliocentric J2000.  
    /// * `reference_epoch`: Epoch of the state vector in Julian Date (TDB).  
    ///
    /// Return
    /// ----------
    /// * An [`OrbitalElements`] enum instance, wrapping either:
    ///   - [`KeplerianElements`] (elliptical), or  
    ///   - [`CometaryElements`] (parabolic/hyperbolic).  
    ///
    /// Note
    /// ----------
    /// * No planetary perturbations or relativistic corrections are applied.  
    ///
    /// See also
    /// ------------
    /// * [`KeplerianElements`] – Orbital elements for elliptical orbits.  
    /// * [`CometaryElements`] – Orbital elements for parabolic and hyperbolic solutions.
    pub fn from_orbital_state(
        position: &Vector3<f64>,
        velocity: &Vector3<f64>,
        reference_epoch: f64,
    ) -> Self {
        ccek1(position, velocity, reference_epoch)
    }

    /// Convert to Keplerian representation with full uncertainty propagation
    ///
    /// This method converts the orbital elements to Keplerian form $(a, e, i, \Omega, \omega, M)$
    /// and, if a covariance matrix is present, propagates it through the transformation using
    /// the appropriate Jacobian matrix.
    ///
    /// ## Uncertainty propagation
    ///
    /// When the source elements have an attached covariance matrix $\Sigma_x$, this method:
    ///
    /// 1. Computes the Jacobian $J = \partial \mathbf{y}_{\text{Kep}} / \partial \mathbf{x}_{\text{src}}$
    ///    where $\mathbf{x}_{\text{src}}$ is the source parameterization and
    ///    $\mathbf{y}_{\text{Kep}} = [a, e, i, \Omega, \omega, M]^\top$
    ///
    /// 2. Transforms the covariance using **linear covariance propagation**:
    ///    $$\Sigma_{\text{Kep}} = J \, \Sigma_x \, J^\top$$
    ///
    /// 3. Extracts 1-σ uncertainties from the diagonal: $\sigma_i = \sqrt{\Sigma_{ii}}$
    ///
    /// This transformation preserves statistical properties under first-order approximation,
    /// which is accurate when uncertainties are small relative to element values.
    ///
    /// ## Conversions
    ///
    /// - **From Keplerian**: returns a clone (no transformation needed)
    /// - **From Equinoctial**: uses [`EquinoctialElements::jacobian_to_keplerian`]
    /// - **From Cometary**: uses [`CometaryElements::jacobian_to_keplerian`]
    ///
    /// ## Return
    ///
    /// * `Ok(OrbitalElements::Keplerian)` – Converted elements with propagated
    ///   uncertainty and covariance when available
    /// * `Err(OutfitError)` – If the conversion is mathematically undefined
    ///   (e.g., parabolic cometary elements with $e = 1$, where $a = \infty$)
    pub fn to_keplerian(&self) -> Result<OrbitalElements, OutfitError> {
        match self {
            OrbitalElements::Keplerian { .. } => Ok(self.clone()),

            OrbitalElements::Equinoctial {
                elements,
                covariance,
                ..
            } => {
                let kep = KeplerianElements::from(elements);
                let jacobian = elements.jacobian_to_keplerian();

                let new_cov = covariance.as_ref().map(|c| c.propagate(&jacobian));
                let new_unc = new_cov.as_ref().map(KeplerianUncertainty::from_covariance);

                Ok(OrbitalElements::Keplerian {
                    elements: kep,
                    uncertainty: new_unc,
                    covariance: new_cov,
                })
            }

            OrbitalElements::Cometary {
                elements,
                covariance,
                ..
            } => {
                let kep = KeplerianElements::try_from(elements)?;
                let jacobian = elements.jacobian_to_keplerian();

                let new_cov = covariance.as_ref().map(|c| c.propagate(&jacobian));
                let new_unc = new_cov.as_ref().map(KeplerianUncertainty::from_covariance);

                Ok(OrbitalElements::Keplerian {
                    elements: kep,
                    uncertainty: new_unc,
                    covariance: new_cov,
                })
            }
        }
    }

    /// Convert to equinoctial representation with full uncertainty propagation
    ///
    /// This method converts the orbital elements to equinoctial form $(a, h, k, p, q, \lambda)$
    /// and, if a covariance matrix is present, propagates it through the transformation using
    /// the appropriate Jacobian matrix.
    ///
    /// ## Uncertainty propagation
    ///
    /// When the source elements have an attached covariance matrix $\Sigma_x$, this method:
    ///
    /// 1. Computes the Jacobian $J = \partial \mathbf{y}_{\text{Eq}} / \partial \mathbf{x}_{\text{src}}$
    ///    where $\mathbf{x}_{\text{src}}$ is the source parameterization and
    ///    $\mathbf{y}_{\text{Eq}} = [a, h, k, p, q, \lambda]^\top$
    ///
    /// 2. Transforms the covariance using **linear covariance propagation**:
    ///    $$\Sigma_{\text{Eq}} = J \, \Sigma_x \, J^\top$$
    ///
    /// 3. Extracts 1-σ uncertainties from the diagonal: $\sigma_i = \sqrt{\Sigma_{ii}}$
    ///
    /// The equinoctial representation is **non-singular** for $e < 1$ and $0 \leq i < \pi$,
    /// making it ideal for uncertainty analysis and propagation when Keplerian elements
    /// would be near-singular.
    ///
    /// ## Conversions
    ///
    /// - **From Equinoctial**: returns a clone (no transformation needed)
    /// - **From Keplerian**: uses [`KeplerianElements::jacobian_to_equinoctial`]
    /// - **From Cometary**: uses [`CometaryElements::jacobian_to_equinoctial`] (chain rule through Keplerian)
    ///
    /// ## Return
    ///
    /// * `Ok(OrbitalElements::Equinoctial)` – Converted elements with propagated
    ///   uncertainty and covariance when available
    /// * `Err(OutfitError)` – If the conversion is mathematically undefined
    pub fn to_equinoctial(&self) -> Result<OrbitalElements, OutfitError> {
        match self {
            OrbitalElements::Equinoctial { .. } => Ok(self.clone()),

            OrbitalElements::Keplerian {
                elements,
                covariance,
                ..
            } => {
                let eq = EquinoctialElements::from(elements);
                let jacobian = elements.jacobian_to_equinoctial();

                let new_cov = covariance.as_ref().map(|c| c.propagate(&jacobian));
                let new_unc = new_cov
                    .as_ref()
                    .map(EquinoctialUncertainty::from_covariance);

                Ok(OrbitalElements::Equinoctial {
                    elements: eq,
                    uncertainty: new_unc,
                    covariance: new_cov,
                })
            }

            OrbitalElements::Cometary {
                elements,
                covariance,
                ..
            } => {
                let eq = EquinoctialElements::try_from(elements)?;
                let jacobian = elements.jacobian_to_equinoctial()?;

                let new_cov = covariance.as_ref().map(|c| c.propagate(&jacobian));
                let new_unc = new_cov
                    .as_ref()
                    .map(EquinoctialUncertainty::from_covariance);

                Ok(OrbitalElements::Equinoctial {
                    elements: eq,
                    uncertainty: new_unc,
                    covariance: new_cov,
                })
            }
        }
    }

    /// Get a reference to the underlying [`KeplerianElements`] if this is `Keplerian`.
    pub fn as_keplerian_ref(&self) -> Option<&KeplerianElements> {
        if let OrbitalElements::Keplerian { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Get the owned underlying [`KeplerianElements`] if this is `Keplerian`.
    pub fn as_keplerian(self) -> Option<KeplerianElements> {
        if let OrbitalElements::Keplerian { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Get a reference to the underlying [`EquinoctialElements`] if this is `Equinoctial`.
    pub fn as_equinoctial_ref(&self) -> Option<&EquinoctialElements> {
        if let OrbitalElements::Equinoctial { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Get the owned underlying [`EquinoctialElements`] if this is `Equinoctial`.
    pub fn as_equinoctial(self) -> Option<EquinoctialElements> {
        if let OrbitalElements::Equinoctial { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Get a reference to the underlying [`CometaryElements`] if this is `Cometary`.
    pub fn as_cometary_ref(&self) -> Option<&CometaryElements> {
        if let OrbitalElements::Cometary { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Get the owned underlying [`CometaryElements`] if this is `Cometary`.
    pub fn as_cometary(self) -> Option<CometaryElements> {
        if let OrbitalElements::Cometary { elements, .. } = self {
            Some(elements)
        } else {
            None
        }
    }

    /// Convert to [`KeplerianElements`], propagating covariance if present.
    ///
    /// Shorthand for `.to_keplerian()?.as_keplerian()`.
    ///
    /// Return
    /// ------
    /// * `Ok(KeplerianElements)` – Converted elements.
    /// * `Err(OutfitError)` – If the conversion is not defined for the current
    ///   element set (e.g. parabolic cometary elements).
    pub fn into_keplerian(self) -> Result<KeplerianElements, OutfitError> {
        self.to_keplerian()?
            .as_keplerian()
            .ok_or(OutfitError::InvalidConversion(
                "Conversion to Keplerian elements failed".to_string(),
            ))
    }

    /// Convert to [`EquinoctialElements`], propagating covariance if present.
    ///
    /// Shorthand for `.to_equinoctial()?.as_equinoctial()`.
    ///
    /// Return
    /// ------
    /// * `Ok(EquinoctialElements)` – Converted elements.
    /// * `Err(OutfitError)` – If the conversion is not defined for the current
    ///   element set (e.g. hyperbolic cometary elements).
    pub fn into_equinoctial(self) -> Result<EquinoctialElements, OutfitError> {
        self.to_equinoctial()?
            .as_equinoctial()
            .ok_or(OutfitError::InvalidConversion(
                "Conversion to equinoctial elements failed".to_string(),
            ))
    }
}

use std::fmt;

impl fmt::Display for OrbitalElements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrbitalElements::Keplerian { elements, .. } => {
                writeln!(f, "[Keplerian]")?;
                write!(f, "{elements}")
            }
            OrbitalElements::Equinoctial { elements, .. } => {
                writeln!(f, "[Equinoctial]")?;
                write!(f, "{elements}")
            }
            OrbitalElements::Cometary { elements, .. } => {
                writeln!(f, "[Cometary]")?;
                write!(f, "{elements}")
            }
        }
    }
}

#[cfg(test)]
pub(crate) mod orbit_type_test {
    use super::*;
    use approx::{abs_diff_eq, assert_abs_diff_eq, assert_relative_eq};
    use std::f64::consts::PI;

    #[allow(dead_code)]
    pub(crate) fn approx_equal(
        current: &OrbitalElements,
        other: &OrbitalElements,
        tol: f64,
    ) -> bool {
        match (current, other) {
            (
                OrbitalElements::Keplerian { elements: ke1, .. },
                OrbitalElements::Keplerian { elements: ke2, .. },
            ) => {
                abs_diff_eq!(ke1.semi_major_axis, ke2.semi_major_axis, epsilon = tol)
                    && abs_diff_eq!(ke1.eccentricity, ke2.eccentricity, epsilon = tol)
                    && abs_diff_eq!(ke1.inclination, ke2.inclination, epsilon = tol)
                    && abs_diff_eq!(
                        ke1.ascending_node_longitude,
                        ke2.ascending_node_longitude,
                        epsilon = tol
                    )
                    && abs_diff_eq!(
                        ke1.periapsis_argument,
                        ke2.periapsis_argument,
                        epsilon = tol
                    )
                    && abs_diff_eq!(ke1.mean_anomaly, ke2.mean_anomaly, epsilon = tol)
            }
            (
                OrbitalElements::Equinoctial { elements: ee1, .. },
                OrbitalElements::Equinoctial { elements: ee2, .. },
            ) => {
                abs_diff_eq!(ee1.semi_major_axis, 0.0, epsilon = tol)
                    && abs_diff_eq!(
                        ee1.eccentricity_sin_lon,
                        ee2.eccentricity_sin_lon,
                        epsilon = tol
                    )
                    && abs_diff_eq!(
                        ee1.eccentricity_cos_lon,
                        ee2.eccentricity_cos_lon,
                        epsilon = tol
                    )
                    && abs_diff_eq!(
                        ee1.tan_half_incl_sin_node,
                        ee2.tan_half_incl_sin_node,
                        epsilon = tol
                    )
                    && abs_diff_eq!(
                        ee1.tan_half_incl_cos_node,
                        ee2.tan_half_incl_cos_node,
                        epsilon = tol
                    )
                    && abs_diff_eq!(ee1.mean_longitude, ee2.mean_longitude, epsilon = tol)
            }
            (
                OrbitalElements::Cometary { elements: ce1, .. },
                OrbitalElements::Cometary { elements: ce2, .. },
            ) => {
                abs_diff_eq!(
                    ce1.perihelion_distance,
                    ce2.perihelion_distance,
                    epsilon = tol
                ) && abs_diff_eq!(ce1.eccentricity, ce2.eccentricity, epsilon = tol)
                    && abs_diff_eq!(ce1.inclination, ce2.inclination, epsilon = tol)
                    && abs_diff_eq!(
                        ce1.ascending_node_longitude,
                        ce2.ascending_node_longitude,
                        epsilon = tol
                    )
                    && abs_diff_eq!(
                        ce1.periapsis_argument,
                        ce2.periapsis_argument,
                        epsilon = tol
                    )
                    && abs_diff_eq!(ce1.true_anomaly, ce2.true_anomaly, epsilon = tol)
            }
            _ => false, // Different types cannot be equal
        }
    }

    /// Degrees to radians helper.
    fn deg(x: f64) -> f64 {
        x * PI / 180.0
    }

    /// Compare two angles modulo 2π with an absolute epsilon.
    fn assert_angle_eq(a: f64, b: f64, eps: f64) {
        fn wrap_to_pi(x: f64) -> f64 {
            // Shift into [-π, π] using rem_euclid for clarity and robustness
            let two_pi = 2.0 * PI;
            (x + PI).rem_euclid(two_pi) - PI
        }

        let d = wrap_to_pi(a - b);
        assert_abs_diff_eq!(d, 0.0, epsilon = eps);
    }

    // ---------- construction from state (elliptic vs hyperbolic) ----------

    #[test]
    fn from_state_returns_keplerian_for_elliptic_orbit() {
        // Example taken from user's earlier vector dump (AU and AU/day).
        let r = Vector3::new(
            -0.632_363_815_650_731_f64,
            1.130_641_635_350_301_7_f64,
            0.491_670_066_794_699_5_f64,
        );
        let v = Vector3::new(
            -0.015_515_672_728_876_466_f64,
            -0.004_456_064_593_471_24_f64,
            -0.001_879_709_408_428_156_1_f64,
        );
        let jd_tdb = 2460000.5_f64;

        let elems = OrbitalElements::from_orbital_state(&r, &v, jd_tdb);

        match elems {
            OrbitalElements::Keplerian { elements: ke, .. } => {
                assert!(ke.semi_major_axis > 0.0);
                // Loose target for a sanity check, not a golden number.
                assert_abs_diff_eq!(ke.semi_major_axis, 1.8155, epsilon = 5e-3);
                assert!(ke.inclination >= 0.0 && ke.inclination <= PI);
                assert!(ke.eccentricity >= 0.0);
            }
            _ => panic!("Expected Keplerian elements for elliptic state"),
        }
    }

    #[test]
    fn from_state_returns_cometary_for_hyperbolic_orbit() {
        // Position at 1 AU on x-axis, velocity set > escape speed on y-axis.
        let r = Vector3::new(1.0_f64, 0.0, 0.0);

        // Use Gauss constant to compute v_escape in AU/day:
        let k = 0.017_202_098_95_f64; // Gauss gravitational constant [AU^(3/2)/day]
        let v_circ = (k * k).sqrt(); // circular speed at 1 AU
        let v_esc = (2.0_f64).sqrt() * v_circ;
        let v = Vector3::new(0.0, 1.05 * v_esc, 0.0);

        let jd_tdb = 2460000.5_f64;
        let elems = OrbitalElements::from_orbital_state(&r, &v, jd_tdb);

        match elems {
            OrbitalElements::Cometary { elements: ce, .. } => {
                assert!(ce.eccentricity >= 1.0);
                assert!(ce.perihelion_distance > 0.0);
            }
            _ => panic!("Expected Cometary elements for hyperbolic state"),
        }
    }

    // ---------- conversions to Keplerian / Equinoctial ----------

    #[test]
    fn to_keplerian_from_keplerian_is_identity_like() {
        let ke = KeplerianElements {
            reference_epoch: 2460000.5,
            semi_major_axis: 1.2,
            eccentricity: 0.1,
            inclination: deg(10.0),
            ascending_node_longitude: deg(20.0),
            periapsis_argument: deg(30.0),
            mean_anomaly: deg(40.0),
        };
        let oe = OrbitalElements::Keplerian {
            elements: ke.clone(),
            uncertainty: None,
            covariance: None,
        };

        let back = oe
            .to_keplerian()
            .expect("Keplerian -> Keplerian should succeed")
            .as_keplerian()
            .expect("Failed to convert to Keplerian");

        assert_abs_diff_eq!(back.semi_major_axis, ke.semi_major_axis, epsilon = 1e-14);
        assert_abs_diff_eq!(back.eccentricity, ke.eccentricity, epsilon = 1e-14);
        assert_abs_diff_eq!(back.inclination, ke.inclination, epsilon = 1e-14);
        assert_abs_diff_eq!(
            back.ascending_node_longitude,
            ke.ascending_node_longitude,
            epsilon = 1e-14
        );
        assert_abs_diff_eq!(
            back.periapsis_argument,
            ke.periapsis_argument,
            epsilon = 1e-14
        );
        assert_angle_eq(back.mean_anomaly, ke.mean_anomaly, 1e-12);
    }

    #[test]
    fn to_equinoctial_from_keplerian_round_trips_reasonably() {
        let ke = KeplerianElements {
            reference_epoch: 2460000.5,
            semi_major_axis: 2.0,
            eccentricity: 0.3,
            inclination: deg(15.0),
            ascending_node_longitude: deg(25.0),
            periapsis_argument: deg(35.0),
            mean_anomaly: deg(45.0),
        };
        let oe = OrbitalElements::Keplerian {
            elements: ke.clone(),
            uncertainty: None,
            covariance: None,
        };

        let eq = oe
            .to_equinoctial()
            .expect("Keplerian -> Equinoctial should succeed")
            .as_equinoctial()
            .expect("Failed to convert to Equinoctial");
        let ke_back = KeplerianElements::from(&eq);

        // Use a mix of absolute and relative checks to be robust to scaling.
        assert_abs_diff_eq!(ke_back.semi_major_axis, ke.semi_major_axis, epsilon = 1e-12);
        assert_relative_eq!(ke_back.eccentricity, ke.eccentricity, max_relative = 1e-10);
        assert_abs_diff_eq!(ke_back.inclination, ke.inclination, epsilon = 1e-12);
        assert_abs_diff_eq!(
            ke_back.ascending_node_longitude,
            ke.ascending_node_longitude,
            epsilon = 1e-12
        );
        assert_abs_diff_eq!(
            ke_back.periapsis_argument,
            ke.periapsis_argument,
            epsilon = 1e-12
        );
        assert_angle_eq(ke_back.mean_anomaly, ke.mean_anomaly, 1e-10);
    }

    #[test]
    fn to_keplerian_from_equinoctial_matches_direct_from() {
        let ke = KeplerianElements {
            reference_epoch: 2460000.5,
            semi_major_axis: 1.7,
            eccentricity: 0.05,
            inclination: deg(5.0),
            ascending_node_longitude: deg(40.0),
            periapsis_argument: deg(80.0),
            mean_anomaly: deg(10.0),
        };
        let eq = EquinoctialElements::from(&ke);
        let oe = OrbitalElements::Equinoctial {
            elements: eq.clone(),
            uncertainty: None,
            covariance: None,
        };

        let back_via_enum = oe
            .to_keplerian()
            .expect("Eq -> Kep should succeed")
            .as_keplerian()
            .expect("Failed to convert to Keplerian");
        let back_via_direct = KeplerianElements::from(&eq);

        assert_abs_diff_eq!(
            back_via_enum.semi_major_axis,
            back_via_direct.semi_major_axis,
            epsilon = 1e-13
        );
        assert_relative_eq!(
            back_via_enum.eccentricity,
            back_via_direct.eccentricity,
            max_relative = 1e-11
        );
        assert_abs_diff_eq!(
            back_via_enum.inclination,
            back_via_direct.inclination,
            epsilon = 1e-13
        );
        assert_angle_eq(
            back_via_enum.ascending_node_longitude,
            back_via_direct.ascending_node_longitude,
            1e-12,
        );
        assert_angle_eq(
            back_via_enum.periapsis_argument,
            back_via_direct.periapsis_argument,
            1e-12,
        );
        assert_angle_eq(
            back_via_enum.mean_anomaly,
            back_via_direct.mean_anomaly,
            1e-10,
        );
    }

    #[test]
    fn cometary_to_keplerian_hyperbolic_succeeds_and_is_consistent() {
        use std::f64::consts::PI;
        let deg = |x: f64| x * PI / 180.0;

        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.9,
            eccentricity: 1.1,
            inclination: deg(12.0),
            ascending_node_longitude: deg(33.0),
            periapsis_argument: deg(45.0),
            true_anomaly: deg(5.0),
        };
        let oe = OrbitalElements::Cometary {
            elements: ce,
            uncertainty: None,
            covariance: None,
        };

        let ke = oe
            .to_keplerian()
            .expect("Cometary(e>1) -> Keplerian should succeed for hyperbolic orbits")
            .as_keplerian()
            .expect("Failed to convert to Keplerian");

        assert!(ke.eccentricity >= 1.0, "expected e >= 1 for hyperbola");
        assert!(ke.semi_major_axis < 0.0, "expected a < 0 for hyperbola");
    }

    #[test]
    fn cometary_to_equinoctial_hyperbolic_maybe_succeeds_and_if_so_is_consistent() {
        use std::f64::consts::PI;
        let deg = |x: f64| x * PI / 180.0;

        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.7,
            eccentricity: 1.3,
            inclination: deg(10.0),
            ascending_node_longitude: deg(20.0),
            periapsis_argument: deg(30.0),
            true_anomaly: deg(15.0),
        };
        let oe = OrbitalElements::Cometary {
            elements: ce,
            uncertainty: None,
            covariance: None,
        };

        if let Ok(eq) = oe.to_equinoctial() {
            let eq = eq
                .as_equinoctial()
                .expect("Failed to convert to Equinoctial after successful conversion");
            assert!(
                eq.semi_major_axis < 0.0,
                "equinoctial a should be < 0 for hyperbolic orbits"
            );
        }
    }

    // ---------- as_* accessors ----------

    #[test]
    fn as_accessors_return_some_only_for_matching_variant() {
        let ke = KeplerianElements {
            reference_epoch: 2460000.5,
            semi_major_axis: 1.0,
            eccentricity: 0.0,
            inclination: 0.0,
            ascending_node_longitude: 0.0,
            periapsis_argument: 0.0,
            mean_anomaly: 0.0,
        };
        let eq = EquinoctialElements::from(&ke);
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 1.0,
            eccentricity: 1.0,
            inclination: 0.0,
            ascending_node_longitude: 0.0,
            periapsis_argument: 0.0,
            true_anomaly: 0.0,
        };

        let oe_k = OrbitalElements::Keplerian {
            elements: ke.clone(),
            uncertainty: None,
            covariance: None,
        };
        assert!(oe_k.as_keplerian_ref().is_some());
        assert!(oe_k.as_equinoctial_ref().is_none());
        assert!(oe_k.as_cometary_ref().is_none());

        let oe_e = OrbitalElements::Equinoctial {
            elements: eq,
            uncertainty: None,
            covariance: None,
        };
        assert!(oe_e.as_keplerian_ref().is_none());
        assert!(oe_e.as_equinoctial_ref().is_some());
        assert!(oe_e.as_cometary_ref().is_none());

        let oe_c = OrbitalElements::Cometary {
            elements: ce,
            uncertainty: None,
            covariance: None,
        };
        assert!(oe_c.as_keplerian_ref().is_none());
        assert!(oe_c.as_equinoctial_ref().is_none());
        assert!(oe_c.as_cometary_ref().is_some());
    }

    // ---------- Display formatting ----------

    #[test]
    fn display_prefix_matches_variant() {
        let ke = KeplerianElements {
            reference_epoch: 2460000.5,
            semi_major_axis: 1.0,
            eccentricity: 0.01,
            inclination: deg(1.0),
            ascending_node_longitude: deg(2.0),
            periapsis_argument: deg(3.0),
            mean_anomaly: deg(4.0),
        };
        let eq = EquinoctialElements::from(&ke);
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.5,
            eccentricity: 1.2,
            inclination: deg(10.0),
            ascending_node_longitude: deg(20.0),
            periapsis_argument: deg(30.0),
            true_anomaly: deg(0.0),
        };

        let s_k = format!(
            "{}",
            OrbitalElements::Keplerian {
                elements: ke,
                uncertainty: None,
                covariance: None,
            }
        );
        assert!(s_k.starts_with("[Keplerian]"));

        let s_e = format!(
            "{}",
            OrbitalElements::Equinoctial {
                elements: eq,
                uncertainty: None,
                covariance: None,
            }
        );
        assert!(s_e.starts_with("[Equinoctial]"));

        let s_c = format!(
            "{}",
            OrbitalElements::Cometary {
                elements: ce,
                uncertainty: None,
                covariance: None,
            }
        );
        assert!(s_c.starts_with("[Cometary]"));
    }
}
