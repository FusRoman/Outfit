//! # Keplerian orbital elements
//!
//! This module defines the [`crate::orbit_type::keplerian_element::KeplerianElements`] struct
//! and its associated conversion routines, providing the **classical orbital element
//! representation** widely used in celestial mechanics.
//!
//! ## What are Keplerian elements?
//!
//! The six Keplerian elements are:
//!
//! 1. **a** – Semi-major axis (AU)  
//! 2. **e** – Eccentricity (unitless)  
//! 3. **i** – Inclination (radians)  
//! 4. **Ω** – Longitude of ascending node (radians)  
//! 5. **ω** – Argument of periapsis (radians)  
//! 6. **M** – Mean anomaly at epoch (radians)  
//!
//! Together with the epoch of reference (usually in Modified Julian Date, MJD),
//! these parameters fully describe an orbit under the two-body approximation.
//!
//! ## Provided functionality
//!
//! - **Conversion** between [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements) and
//!   [`crate::orbit_type::equinoctial_element::EquinoctialElements`].
//!   * [`KeplerianElements::from_equinoctial_internal`](crate::orbit_type::keplerian_element::KeplerianElements::from_equinoctial_internal) – inverse mapping from equinoctial form.
//!   * [`From<KeplerianElements>`] and [`From<&KeplerianElements>`] – forward mapping.
//!
//! - **Normalization** of angular parameters using [`crate::kepler::principal_angle`],
//!   ensuring all angles are restricted to the range `[0, 2π)`.
//!
//! ## Units
//!
//! - Lengths: **AU**
//! - Angles: **radians**
//! - Time: **days** (epoch usually in **MJD**, TDB/TT scale)
//!
//! ## Degeneracies
//!
//! Classical Keplerian elements suffer from singularities:
//!
//! - **Circular orbits (`e → 0`)**: periapsis argument ω becomes undefined.  
//!   → conventionally set to `0.0` during conversion.
//!
//! - **Equatorial orbits (`i → 0`)**: ascending node Ω becomes undefined.  
//!   → conventionally set to `0.0` during conversion.
//!
//! For robust numerical work, the [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements) representation is recommended.
//!
//! ## Example
//!
//! ```rust, no_run
//! use outfit::orbit_type::keplerian_element::KeplerianElements;
//! use outfit::orbit_type::equinoctial_element::EquinoctialElements;
//!
//! // Define Keplerian elements
//! let kep = KeplerianElements {
//!     reference_epoch: 59000.0,
//!     semi_major_axis: 1.0,
//!     eccentricity: 0.0167,
//!     inclination: 0.4091,
//!     ascending_node_longitude: 0.0,
//!     periapsis_argument: 1.7966,
//!     mean_anomaly: 0.0,
//! };
//!
//! // Convert to equinoctial elements
//! let equ: EquinoctialElements = (&kep).into();
//!
//! // Convert back to Keplerian
//! let kep2 = KeplerianElements::from_equinoctial_internal(
//!     equ.reference_epoch,
//!     equ.semi_major_axis,
//!     equ.eccentricity_sin_lon,
//!     equ.eccentricity_cos_lon,
//!     equ.tan_half_incl_sin_node,
//!     equ.tan_half_incl_cos_node,
//!     equ.mean_longitude,
//! );
//! ```
//!
//! ## See also
//!
//! - [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements) – regularized, non-singular form.
//! - [`principal_angle`](crate::kepler::principal_angle) – helper to normalize angular elements.
//! - Milani & Gronchi, *Theory of Orbit Determination* (2010).

use crate::{kepler::principal_angle, orbit_type::equinoctial_element::EquinoctialElements};
use std::fmt;

/// Keplerian orbital elements (osculating, two-body).
///
/// Units
/// -----
/// * `reference_epoch`: MJD (Modified Julian Date).
/// * `semi_major_axis`: Astronomical Units (AU).
/// * `eccentricity`: unitless.
/// * `inclination`: radians.
/// * `ascending_node_longitude`: radians (Ω).
/// * `periapsis_argument`: radians (ω).
/// * `mean_anomaly`: radians (M).
///
/// Notes
/// -----
/// This struct represents the classical Keplerian element set:
/// (a, e, i, Ω, ω, M) at a given epoch. It is primarily intended
/// as a user-facing representation; internal propagation routines
/// can work in equinoctial form for numerical robustness.
///
/// See also
/// --------
/// * [`EquinoctialElements`] – Regularized element set used internally.
/// * [`principal_angle`] – Angle normalization helper to [0, 2π).
#[derive(Debug, PartialEq, Clone)]
pub struct KeplerianElements {
    pub reference_epoch: f64,
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub ascending_node_longitude: f64,
    pub periapsis_argument: f64,
    pub mean_anomaly: f64,
}

impl KeplerianElements {
    /// Convert equinoctial elements (internal representation) to Keplerian elements.
    ///
    /// This routine implements the inverse mapping from equinoctial parameters
    /// `(a, h ≡ e·sin(ϖ), k ≡ e·cos(ϖ), p ≡ tan(i/2)·sin(Ω), q ≡ tan(i/2)·cos(Ω), λ ≡ M+ϖ)`
    /// to classical Keplerian elements `(a, e, i, Ω, ω, M)`, where ϖ = ω + Ω is the
    /// longitude of periapsis and λ the mean longitude.
    ///
    /// Degenerate cases are handled as follows:
    /// - If `e ≈ 0`, we set `ϖ = 0` (thus ω = −Ω mod 2π) to avoid undefined direction of periapsis.
    /// - If `tan(i/2) ≈ 0`, we set `Ω = 0` (node undefined on the reference plane).
    ///   All final angles are normalized to `[0, 2π)` using [`principal_angle`].
    ///
    /// Arguments
    /// ---------
    /// * `reference_epoch` – Epoch of validity (MJD).
    /// * `semi_major_axis` – Semi-major axis `a` (AU).
    /// * `eccentricity_sin_lon` – `h = e·sin(ϖ)` with ϖ = ω + Ω.
    /// * `eccentricity_cos_lon` – `k = e·cos(ϖ)` with ϖ = ω + Ω.
    /// * `tan_half_incl_sin_node` – `p = tan(i/2)·sin(Ω)`.
    /// * `tan_half_incl_cos_node` – `q = tan(i/2)·cos(Ω)`.
    /// * `mean_longitude` – `λ = M + ϖ` (radians).
    ///
    /// Return
    /// ------
    /// * A new [`KeplerianElements`] instance `(a, e, i, Ω, ω, M)` at `reference_epoch`.
    ///
    /// See also
    /// --------
    /// * [`principal_angle`] – Angle normalization helper.
    pub fn from_equinoctial_internal(
        reference_epoch: f64,
        semi_major_axis: f64,
        eccentricity_sin_lon: f64,
        eccentricity_cos_lon: f64,
        tan_half_incl_sin_node: f64,
        tan_half_incl_cos_node: f64,
        mean_longitude: f64,
    ) -> Self {
        let eps = 1.0e-12; // thresholds for near-circular / near-equatorial tests
        let a = semi_major_axis;

        // e = sqrt(h^2 + k^2)
        let ecc = (eccentricity_sin_lon.powi(2) + eccentricity_cos_lon.powi(2)).sqrt();

        // ϖ = atan2(h, k), undefined when e ~ 0
        let dig = if ecc < eps {
            0.0
        } else {
            eccentricity_sin_lon.atan2(eccentricity_cos_lon)
        };

        // t = tan(i/2) = sqrt(p^2 + q^2)
        let tgi2 = (tan_half_incl_sin_node.powi(2) + tan_half_incl_cos_node.powi(2)).sqrt();

        // Ω = atan2(p, q), undefined when t ~ 0
        let omega_node = if tgi2 < eps {
            0.0
        } else {
            tan_half_incl_sin_node.atan2(tan_half_incl_cos_node)
        };

        // i = 2 atan(t)
        let inclination = 2.0 * tgi2.atan();

        // ω = ϖ − Ω ,  M = λ − ϖ
        let periapsis_arg = principal_angle(dig - omega_node);
        let mean_anomaly = principal_angle(mean_longitude - dig);

        Self {
            reference_epoch,
            semi_major_axis: a,
            eccentricity: ecc,
            inclination,
            ascending_node_longitude: omega_node,
            periapsis_argument: periapsis_arg,
            mean_anomaly,
        }
    }
}

impl From<KeplerianElements> for EquinoctialElements {
    /// Forward conversion to equinoctial elements (by value).
    ///
    /// Arguments
    /// ---------
    /// * `k` – Keplerian elements `(a, e, i, Ω, ω, M)`.
    ///
    /// Return
    /// ------
    /// * [`EquinoctialElements`] with components `(a, h, k, p, q, λ)`.
    fn from(k: KeplerianElements) -> Self {
        EquinoctialElements::from_kepler_internal(
            k.reference_epoch,
            k.semi_major_axis,
            k.eccentricity,
            k.inclination,
            k.ascending_node_longitude,
            k.periapsis_argument,
            k.mean_anomaly,
        )
    }
}

impl From<&KeplerianElements> for EquinoctialElements {
    /// Forward conversion to equinoctial elements (by reference).
    ///
    /// Arguments
    /// ---------
    /// * `k` – Reference to Keplerian elements `(a, e, i, Ω, ω, M)`.
    ///
    /// Return
    /// ------
    /// * [`EquinoctialElements`] with components `(a, h, k, p, q, λ)`.
    fn from(k: &KeplerianElements) -> Self {
        EquinoctialElements::from_kepler_internal(
            k.reference_epoch,
            k.semi_major_axis,
            k.eccentricity,
            k.inclination,
            k.ascending_node_longitude,
            k.periapsis_argument,
            k.mean_anomaly,
        )
    }
}

/// Keplerian orbital elements (osculating, two-body).
impl fmt::Display for KeplerianElements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        writeln!(
            f,
            "Keplerian Elements @ epoch (MJD): {:.6}",
            self.reference_epoch
        )?;
        writeln!(f, "-------------------------------------------")?;
        writeln!(
            f,
            "  a   (semi-major axis)       = {:.6} AU",
            self.semi_major_axis
        )?;
        writeln!(
            f,
            "  e   (eccentricity)          = {:.6}",
            self.eccentricity
        )?;
        writeln!(
            f,
            "  i   (inclination)           = {:.6} rad ({:.6}°)",
            self.inclination,
            self.inclination * rad_to_deg
        )?;
        writeln!(
            f,
            "  Ω   (longitude of node)     = {:.6} rad ({:.6}°)",
            self.ascending_node_longitude,
            self.ascending_node_longitude * rad_to_deg
        )?;
        writeln!(
            f,
            "  ω   (argument of periapsis) = {:.6} rad ({:.6}°)",
            self.periapsis_argument,
            self.periapsis_argument * rad_to_deg
        )?;
        writeln!(
            f,
            "  M   (mean anomaly)          = {:.6} rad ({:.6}°)",
            self.mean_anomaly,
            self.mean_anomaly * rad_to_deg
        )
    }
}

#[cfg(test)]
pub(crate) mod test_keplerian_element {
    use super::*;
    use crate::orbit_type::equinoctial_element::EquinoctialElements;

    #[test]
    fn test_keplerian_conversion() {
        let kepler = KeplerianElements {
            reference_epoch: 0.0,
            semi_major_axis: 1.8017360713,
            eccentricity: 0.2835591457,
            inclination: 0.2026738329,
            ascending_node_longitude: 0.0079559790,
            periapsis_argument: 1.2451951388,
            mean_anomaly: 0.4405458902,
        };

        let equinoctial: EquinoctialElements = kepler.into();

        assert_eq!(
            equinoctial,
            EquinoctialElements {
                reference_epoch: 0.0,
                semi_major_axis: 1.8017360713,
                eccentricity_sin_lon: 0.2693736809404963,
                eccentricity_cos_lon: 0.08856415260522467,
                tan_half_incl_sin_node: 0.0008089970142830734,
                tan_half_incl_cos_node: 0.10168201110394352,
                mean_longitude: 1.693697008,
            }
        );
    }
}
