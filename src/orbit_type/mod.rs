use crate::orbit_type::{
    cometary_element::CometaryElements, equinoctial_element::EquinoctialElements,
    keplerian_element::KeplerianElements,
};

/// Equinoctial orbital elements and related conversions.
pub mod equinoctial_element;

/// Classical Keplerian elements structure and utilities.
pub mod keplerian_element;

/// Cometary (parabolic/hyperbolic) orbital elements and related conversions.
pub mod cometary_element;

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
/// * [`OrbitalElements::to_equinoctial`] – Lossless conversion for elliptic orbits.
/// * [`OrbitalElements::to_keplerian`] – Conversion with domain checks.
#[derive(Debug, Clone, PartialEq)]
pub enum OrbitalElements {
    Keplerian(KeplerianElements),
    Equinoctial(EquinoctialElements),
    Cometary(CometaryElements),
}

impl OrbitalElements {
    /// Convert to Keplerian elements, if possible.
    ///
    /// This conversion may fail if the current representation is not suitable
    /// for Keplerian elements. In particular, `Cometary` elements with
    /// eccentricity e < 1 cannot be converted to Keplerian form.
    ///
    /// # Errors
    ///
    /// Returns an `OutfitError::InvalidOrbit` if the conversion is not possible.
    pub fn to_keplerian(&self) -> Result<KeplerianElements, crate::outfit_errors::OutfitError> {
        match self {
            OrbitalElements::Keplerian(ke) => Ok(ke.clone()),
            OrbitalElements::Equinoctial(ee) => Ok(KeplerianElements::from(ee)),
            OrbitalElements::Cometary(ce) => KeplerianElements::try_from(ce),
        }
    }

    /// Convert to Equinoctial elements, if possible.
    ///
    /// This conversion may fail if the current representation is not suitable
    /// for Equinoctial elements. In particular, `Cometary` elements with
    /// eccentricity e < 1 cannot be converted to Equinoctial form.
    ///
    /// # Errors
    ///
    /// Returns an `OutfitError::InvalidOrbit` if the conversion is not possible.
    pub fn to_equinoctial(&self) -> Result<EquinoctialElements, crate::outfit_errors::OutfitError> {
        match self {
            OrbitalElements::Keplerian(ke) => Ok(EquinoctialElements::from(ke)),
            OrbitalElements::Equinoctial(ee) => Ok(ee.clone()),
            OrbitalElements::Cometary(ce) => EquinoctialElements::try_from(ce),
        }
    }

    pub fn as_keplerian(&self) -> Option<&KeplerianElements> {
        if let OrbitalElements::Keplerian(ref k) = self {
            Some(k)
        } else {
            None
        }
    }

    pub fn as_equinoctial(&self) -> Option<&EquinoctialElements> {
        if let OrbitalElements::Equinoctial(ref e) = self {
            Some(e)
        } else {
            None
        }
    }

    pub fn as_cometary(&self) -> Option<&CometaryElements> {
        if let OrbitalElements::Cometary(ref c) = self {
            Some(c)
        } else {
            None
        }
    }
}

use std::fmt;

impl fmt::Display for OrbitalElements {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OrbitalElements::Keplerian(k) => {
                writeln!(f, "[Keplerian representation]")?;
                write!(f, "{k}")
            }
            OrbitalElements::Equinoctial(e) => {
                writeln!(f, "[Equinoctial representation]")?;
                write!(f, "{e}")
            }
            OrbitalElements::Cometary(c) => {
                writeln!(f, "[Cometary representation]")?;
                write!(f, "{c}")
            }
        }
    }
}

#[cfg(test)]
pub(crate) mod orbit_type_test {
    use super::*;
    use approx::abs_diff_eq;

    pub fn approx_equal(current: &OrbitalElements, other: &OrbitalElements, tol: f64) -> bool {
        match (current, other) {
            (OrbitalElements::Keplerian(ke1), OrbitalElements::Keplerian(ke2)) => {
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
            (OrbitalElements::Equinoctial(ee1), OrbitalElements::Equinoctial(ee2)) => {
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
            (OrbitalElements::Cometary(ce1), OrbitalElements::Cometary(ce2)) => {
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
}
