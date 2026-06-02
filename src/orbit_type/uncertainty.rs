//! Orbital uncertainty representation and covariance propagation
//!
//! This module defines uncertainty structures for orbital elements and provides mechanisms
//! for propagating uncertainties through coordinate transformations.
//!
//! ## Overview
//!
//! Orbital uncertainties arise from multiple sources:
//! - **Measurement errors** in astrometric observations (position, timing)
//! - **Numerical approximations** in orbit fitting and propagation
//! - **Model limitations** (unmodeled perturbations, approximations in dynamics)
//!
//! These uncertainties are represented in two complementary forms:
//!
//! 1. **Standard deviations** ($1\sigma$ uncertainties) — individual uncertainties on each element,
//!    extracted from the diagonal of the covariance matrix
//! 2. **Covariance matrices** — full $6 \times 6$ symmetric matrices capturing correlations between elements
//!
//! ## Uncertainty structures
//!
//! Three uncertainty structures correspond to the three orbital element representations:
//!
//! - [`KeplerianUncertainty`] — for classical Keplerian elements $(a, e, i, \Omega, \omega, M)$
//! - [`EquinoctialUncertainty`] — for equinoctial elements $(a, h, k, p, q, \lambda)$
//! - [`CometaryUncertainty`] — for cometary elements $(q, e, i, \Omega, \omega, \nu)$
//!
//! Each structure stores the standard deviation $\sigma_i$ for each element $x_i$, extracted from
//! the diagonal of the covariance matrix: $\sigma_i = \sqrt{\Sigma_{ii}}$.
//!
//! All units match those of the parent element struct (AU for distances, radians for angles).
//!
//! ## Covariance matrix representation
//!
//! The [`OrbitalCovariance`] structure holds a full $6 \times 6$ **symmetric positive semi-definite**
//! covariance matrix $\Sigma$. The matrix element $\Sigma_{ij}$ represents the covariance between
//! orbital elements $x_i$ and $x_j$:
//!
//! $$
//! \Sigma_{ij} = \text{Cov}(x_i, x_j) = E[(x_i - \mu_i)(x_j - \mu_j)]
//! $$
//!
//! where $\mu_i$ is the expected value (nominal element value), with $\mu_i = \mathbb{E}(x_i)$ and $\mathbb{E}(\cdot)$ denoting expectation.
//!
//! **Diagonal entries** $\Sigma_{ii} = \sigma_i^2$ are the variances of each element.
//!
//! **Off-diagonal entries** $\Sigma_{ij}$ ($i \neq j$) capture correlations. The correlation coefficient is:
//!
//! $$
//! \rho_{ij} = \frac{\Sigma_{ij}}{\sigma_i \sigma_j} \in [-1, 1]
//! $$
//!
//! ## Linear covariance propagation
//!
//! When transforming orbital elements from one representation to another, the covariance matrix
//! is propagated using **first-order linear approximation**:
//!
//! $$
//! \Sigma_y = J \, \Sigma_x \, J^\top
//! $$
//!
//! where:
//! - $\Sigma_x$ is the covariance in the source representation $\mathbf{x}$
//! - $\Sigma_y$ is the covariance in the target representation $\mathbf{y}$
//! - $J = \frac{\partial \mathbf{y}}{\partial \mathbf{x}}$ is the $6 \times 6$ Jacobian matrix
//!   evaluated at the nominal element values
//!
//! This transformation preserves the **statistical properties** of the uncertainty distribution
//! under the assumption that:
//! 1. The transformation $\mathbf{y} = f(\mathbf{x})$ is smooth and differentiable
//! 2. Uncertainties are small relative to element values (linear approximation valid)
//! 3. The distribution is approximately Gaussian (covariance fully characterizes uncertainty)
//!
//! ## Jacobian computation
//!
//! Each orbital element representation provides methods to compute analytical Jacobians:
//!
//! - [`KeplerianElements::jacobian_to_equinoctial`](crate::orbit_type::keplerian_element::KeplerianElements::jacobian_to_equinoctial)
//!   — $J_{\text{Kep} \to \text{Eq}} = \partial(a,h,k,p,q,\lambda) / \partial(a,e,i,\Omega,\omega,M)$
//!
//! - [`EquinoctialElements::jacobian_to_keplerian`](crate::orbit_type::equinoctial_element::EquinoctialElements::jacobian_to_keplerian)
//!   — $J_{\text{Eq} \to \text{Kep}} = \partial(a,e,i,\Omega,\omega,M) / \partial(a,h,k,p,q,\lambda)$
//!
//! - [`CometaryElements::jacobian_to_keplerian`](crate::orbit_type::cometary_element::CometaryElements::jacobian_to_keplerian)
//!   — $J_{\text{Com} \to \text{Kep}} = \partial(a,e,i,\Omega,\omega,M) / \partial(q,e,i,\Omega,\omega,\nu)$
//!
//! - [`CometaryElements::jacobian_to_equinoctial`](crate::orbit_type::cometary_element::CometaryElements::jacobian_to_equinoctial)
//!   — $J_{\text{Com} \to \text{Eq}}$ computed via chain rule
//!
//! These Jacobians are computed analytically for numerical accuracy and efficiency.
//!
//! ## Handling singularities
//!
//! Some orbital element representations have **singularities** where derivatives are undefined:
//!
//! - **Keplerian elements**: singular when $e \to 0$ (circular) or $i \to 0$ (equatorial)
//! - **Equinoctial elements**: non-singular for $e < 1$ and $0 \leq i < \pi$
//!
//! When computing Jacobians near singular points, derivatives involving undefined quantities
//! (e.g., $\partial \varpi / \partial h$ when $e = 0$) are set to zero. This ensures numerical
//! stability but may underestimate uncertainties in degenerate configurations.
//!
//! For orbits near circular or equatorial, **equinoctial elements** should be preferred as the
//! primary representation to avoid singularities in both the transformation and its Jacobian.
//!
//! ## Usage
//!
//! Uncertainties are typically created by orbit determination codes and attached to
//! [`OrbitalElements`](crate::orbit_type::OrbitalElements) variants. When converting between
//! representations, the covariance is automatically propagated:
//!
//! ```rust
//! # use outfit::orbit_type::{OrbitalElements, keplerian_element::KeplerianElements};
//! # use outfit::orbit_type::uncertainty::{KeplerianUncertainty, OrbitalCovariance};
//! # use nalgebra::Matrix6;
//! // Assume we have Keplerian elements with uncertainty
//! let kep = KeplerianElements {
//!     reference_epoch: 60000.0,
//!     semi_major_axis: 2.5,
//!     eccentricity: 0.1,
//!     inclination: 0.2,
//!     ascending_node_longitude: 1.0,
//!     periapsis_argument: 0.5,
//!     mean_anomaly: 2.0,
//! };
//!
//! // With an associated covariance (example: diagonal matrix)
//! let cov = OrbitalCovariance {
//!     matrix: Matrix6::from_diagonal_element(1e-6),
//! };
//!
//! let oe = OrbitalElements::Keplerian {
//!     elements: kep,
//!     uncertainty: Some(KeplerianUncertainty::from_covariance(&cov)),
//!     covariance: Some(cov),
//! };
//!
//! // Convert to equinoctial — covariance automatically propagated
//! let oe_eq = oe.to_equinoctial().unwrap();
//!
//! // Uncertainty in equinoctial representation now available
//! if let OrbitalElements::Equinoctial { uncertainty, .. } = oe_eq {
//!     if let Some(unc) = uncertainty {
//!         println!("Uncertainty in h: {}", unc.eccentricity_sin_lon);
//!     }
//! }
//! ```
//!
//! ## See also
//!
//! - [`OrbitalElements`](crate::orbit_type::OrbitalElements) — Container with uncertainty
//! - [`OrbitalCovariance::propagate`] — Covariance transformation method
//! - Module-level documentation in [`crate::orbit_type`] for conversion details
//!
//! ## References
//!
//! - Tapley, Schutz, & Born, *Statistical Orbit Determination* (2004), Chapter 4
//! - Milani & Gronchi, *Theory of Orbit Determination* (2010), Chapter 5
use nalgebra::Matrix6;
use photom::Radians;

/// One-sigma uncertainties on Keplerian elements.
///
/// Units
/// -----
/// * `semi_major_axis`:            AU.
/// * `eccentricity`:               unitless.
/// * `inclination`:                radians.
/// * `ascending_node_longitude`:   radians.
/// * `periapsis_argument`:         radians.
/// * `mean_anomaly`:               radians.
///
/// Notes
/// -----
/// The reference epoch is treated as exact (no uncertainty propagated here).
///
/// See also
/// --------
/// * [`crate::orbit_type::keplerian_element::KeplerianElements`]  – Associated element struct.
/// * [`OrbitalCovariance`]  – Full 6×6 covariance when available.
#[derive(Debug, Clone, PartialEq)]
pub struct KeplerianUncertainty {
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: Radians,
    pub ascending_node_longitude: Radians,
    pub periapsis_argument: Radians,
    pub mean_anomaly: Radians,
}

/// One-sigma uncertainties on equinoctial elements.
///
/// Units
/// -----
/// * `semi_major_axis`:          AU.
/// * `eccentricity_sin_lon`:     unitless.
/// * `eccentricity_cos_lon`:     unitless.
/// * `tan_half_incl_sin_node`:   unitless.
/// * `tan_half_incl_cos_node`:   unitless.
/// * `mean_longitude`:           radians.
///
/// See also
/// --------
/// * [`crate::orbit_type::equinoctial_element::EquinoctialElements`] – Associated element struct.
/// * [`OrbitalCovariance`]   – Full 6×6 covariance when available.
#[derive(Debug, Clone, PartialEq)]
pub struct EquinoctialUncertainty {
    pub semi_major_axis: f64,
    pub eccentricity_sin_lon: f64,
    pub eccentricity_cos_lon: f64,
    pub tan_half_incl_sin_node: f64,
    pub tan_half_incl_cos_node: f64,
    pub mean_longitude: Radians,
}

/// One-sigma uncertainties on cometary elements.
///
/// Units
/// -----
/// * `perihelion_distance`:        AU.
/// * `eccentricity`:               unitless.
/// * `inclination`:                radians.
/// * `ascending_node_longitude`:   radians.
/// * `periapsis_argument`:         radians.
/// * `true_anomaly`:               radians.
///
/// See also
/// --------
/// * [`crate::orbit_type::cometary_element::CometaryElements`]  – Associated element struct.
/// * [`OrbitalCovariance`] – Full 6×6 covariance when available.
#[derive(Debug, Clone, PartialEq)]
pub struct CometaryUncertainty {
    pub perihelion_distance: f64,
    pub eccentricity: f64,
    pub inclination: Radians,
    pub ascending_node_longitude: Radians,
    pub periapsis_argument: Radians,
    pub true_anomaly: Radians,
}

impl KeplerianUncertainty {
    /// Extract 1-σ standard deviations from the diagonal of a covariance matrix.
    ///
    /// Element ordering: $[a, e, i, \Omega, \omega, M]$.
    pub fn from_covariance(cov: &OrbitalCovariance) -> Self {
        let v = cov.variances();
        KeplerianUncertainty {
            semi_major_axis: v[0].sqrt(),
            eccentricity: v[1].sqrt(),
            inclination: v[2].sqrt(),
            ascending_node_longitude: v[3].sqrt(),
            periapsis_argument: v[4].sqrt(),
            mean_anomaly: v[5].sqrt(),
        }
    }
}

impl EquinoctialUncertainty {
    /// Extract 1-σ standard deviations from the diagonal of a covariance matrix.
    ///
    /// Element ordering: $[a, h, k, p, q, \lambda]$.
    pub fn from_covariance(cov: &OrbitalCovariance) -> Self {
        let v = cov.variances();
        EquinoctialUncertainty {
            semi_major_axis: v[0].sqrt(),
            eccentricity_sin_lon: v[1].sqrt(),
            eccentricity_cos_lon: v[2].sqrt(),
            tan_half_incl_sin_node: v[3].sqrt(),
            tan_half_incl_cos_node: v[4].sqrt(),
            mean_longitude: v[5].sqrt(),
        }
    }
}

impl CometaryUncertainty {
    /// Extract 1-σ standard deviations from the diagonal of a covariance matrix.
    ///
    /// Element ordering: $[q, e, i, \Omega, \omega, \nu]$.
    pub fn from_covariance(cov: &OrbitalCovariance) -> Self {
        let v = cov.variances();
        CometaryUncertainty {
            perihelion_distance: v[0].sqrt(),
            eccentricity: v[1].sqrt(),
            inclination: v[2].sqrt(),
            ascending_node_longitude: v[3].sqrt(),
            periapsis_argument: v[4].sqrt(),
            true_anomaly: v[5].sqrt(),
        }
    }
}

/// Full 6×6 covariance matrix for a set of orbital elements.
///
/// The matrix is stored as a flat array of 36 elements in **row-major** order.
/// The parameterization (Keplerian, equinoctial, or cometary) is determined by
/// the variant it accompanies inside [`crate::orbit_type::OrbitalElements`].
///
/// The matrix is assumed to be **symmetric positive semi-definite**.
/// Only the upper triangle is guaranteed to be written by solvers; the lower
/// triangle is kept consistent by convention.
///
/// Element ordering
/// ----------------
/// For Keplerian elements: `[a, e, i, Ω, ω, M]`.  
/// For equinoctial elements: `[a, h, k, p, q, λ]`.  
/// For cometary elements: `[q, e, i, Ω, ω, ν]`.
///
/// See also
/// --------
/// * [`crate::orbit_type::OrbitalElements`] – Container that holds this matrix.
/// * [`KeplerianUncertainty`], [`EquinoctialUncertainty`], [`CometaryUncertainty`] – Diagonal 1-σ summaries.
#[derive(Debug, Clone, PartialEq)]
pub struct OrbitalCovariance {
    /// Row-major 6×6 covariance matrix.
    pub matrix: Matrix6<f64>,
}

impl OrbitalCovariance {
    /// Returns the diagonal entries, i.e., the variances of each element.
    ///
    /// The square roots of the returned values equal the 1-σ standard
    /// deviations stored in the matching uncertainty struct
    /// ([`KeplerianUncertainty`], [`EquinoctialUncertainty`], or [`CometaryUncertainty`]).
    #[inline]
    pub fn variances(&self) -> [f64; 6] {
        [
            self.matrix[(0, 0)],
            self.matrix[(1, 1)],
            self.matrix[(2, 2)],
            self.matrix[(3, 3)],
            self.matrix[(4, 4)],
            self.matrix[(5, 5)],
        ]
    }

    /// Propagate this covariance matrix through a coordinate transformation
    ///
    /// This method implements **linear covariance propagation**, transforming uncertainties
    /// from one orbital element representation to another using the Jacobian matrix of the
    /// transformation.
    ///
    /// ## Mathematical formulation
    ///
    /// Given a transformation $\mathbf{y} = f(\mathbf{x})$ where $\mathbf{x}$ and $\mathbf{y}$
    /// are 6-element orbital parameter vectors, the covariance in the target space is:
    ///
    /// $$
    /// \Sigma_y = J \, \Sigma_x \, J^\top
    /// $$
    ///
    /// where:
    /// - $\Sigma_x$ is the $6 \times 6$ covariance matrix in the source representation (this matrix)
    /// - $\Sigma_y$ is the $6 \times 6$ covariance matrix in the target representation (returned)
    /// - $J = \frac{\partial \mathbf{y}}{\partial \mathbf{x}} \bigg|_{\mathbf{x}_0}$
    ///   is the Jacobian matrix evaluated at the nominal element values $\mathbf{x}_0$
    ///
    /// ## Derivation
    ///
    /// Under a first-order Taylor expansion around the nominal point $\mathbf{x}_0$:
    ///
    /// $$
    /// \mathbf{y} \approx f(\mathbf{x}_0) + J(\mathbf{x} - \mathbf{x}_0)
    /// $$
    ///
    /// The covariance of $\mathbf{y}$ is:
    ///
    /// $$
    /// \begin{aligned}
    /// \Sigma_y &= E[(\mathbf{y} - E[\mathbf{y}])(\mathbf{y} - E[\mathbf{y}])^\top] \\
    /// &= E[J(\mathbf{x} - \mathbf{x}_0)(\mathbf{x} - \mathbf{x}_0)^\top J^\top] \\
    /// &= J \, E[(\mathbf{x} - \mathbf{x}_0)(\mathbf{x} - \mathbf{x}_0)^\top] \, J^\top \\
    /// &= J \, \Sigma_x \, J^\top
    /// </end{aligned}
    /// $$
    ///
    /// ## Validity conditions
    ///
    /// This linear approximation is accurate when:
    /// 1. **Small uncertainties**: $\|\Sigma_x\|$ is small relative to $\|\mathbf{x}_0\|^2$,
    ///    ensuring the linear term dominates the Taylor expansion
    /// 2. **Smooth transformation**: $f$ is twice continuously differentiable in the uncertainty region
    /// 3. **Gaussian distribution**: The underlying probability distribution is approximately Gaussian,
    ///    so the covariance matrix fully characterizes the uncertainty
    ///
    /// For large uncertainties or highly nonlinear transformations, higher-order methods
    /// (unscented transform, Monte Carlo) may be more appropriate.
    ///
    /// ## Symmetry preservation
    ///
    /// The result $\Sigma_y$ is guaranteed to be symmetric because:
    ///
    /// $$
    /// \Sigma_y^\top = (J \Sigma_x J^\top)^\top = J \Sigma_x^\top J^\top = J \Sigma_x J^\top = \Sigma_y
    /// $$
    ///
    /// provided $\Sigma_x$ is symmetric (as required for a valid covariance matrix).
    ///
    /// ## Arguments
    ///
    /// * `jacobian` – The $6 \times 6$ Jacobian matrix $J = \partial\mathbf{y}/\partial\mathbf{x}$
    ///   evaluated at the nominal element values
    ///
    /// ## Return
    ///
    /// * `OrbitalCovariance` – The propagated covariance matrix $\Sigma_y$ in the target
    ///   parameterization
    ///
    /// ## See also
    ///
    /// - [`KeplerianElements::jacobian_to_equinoctial`](crate::orbit_type::keplerian_element::KeplerianElements::jacobian_to_equinoctial)
    /// - [`EquinoctialElements::jacobian_to_keplerian`](crate::orbit_type::equinoctial_element::EquinoctialElements::jacobian_to_keplerian)
    /// - [`CometaryElements::jacobian_to_keplerian`](crate::orbit_type::cometary_element::CometaryElements::jacobian_to_keplerian)
    pub fn propagate(&self, jacobian: &Matrix6<f64>) -> OrbitalCovariance {
        OrbitalCovariance {
            matrix: jacobian * self.matrix * jacobian.transpose(),
        }
    }
}
