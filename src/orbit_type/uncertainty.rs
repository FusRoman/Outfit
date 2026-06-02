/// Uncertainty (1-σ) associated with a set of orbital elements.
///
/// Each field corresponds to a diagonal entry of the covariance matrix,
/// i.e., the standard deviation of the matching element.  Off-diagonal
/// correlations are carried by [`OrbitalCovariance`] when a full matrix
/// is available.
///
/// All units match those of the parent element struct.
///
/// See also
/// --------
/// * [`KeplerianUncertainty`]    – Uncertainties for [`crate::orbit_type::keplerian_element::KeplerianElements`].
/// * [`EquinoctialUncertainty`]  – Uncertainties for [`crate::orbit_type::equinoctial_element::EquinoctialElements`].
/// * [`CometaryUncertainty`]     – Uncertainties for [`crate::orbit_type::cometary_element::CometaryElements`].
/// * [`OrbitalCovariance`]       – Full covariance matrix (optional).
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

    /// Propagate this covariance matrix through a linear map.
    ///
    /// Computes $\Sigma_y = J \, \Sigma_x \, J^\top$ where $J$ is the
    /// Jacobian of the transformation evaluated at the nominal element set.
    ///
    /// Arguments
    /// ---------
    /// * `jacobian` – $6 \times 6$ Jacobian matrix $\partial\mathbf{y}/\partial\mathbf{x}$.
    ///
    /// Return
    /// ------
    /// * `OrbitalCovariance` – Propagated covariance in the target parameterization.
    pub fn propagate(&self, jacobian: &Matrix6<f64>) -> OrbitalCovariance {
        OrbitalCovariance {
            matrix: jacobian * self.matrix * jacobian.transpose(),
        }
    }
}
