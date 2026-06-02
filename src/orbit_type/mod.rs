//! # Orbital element representations
//!
//! This module defines multiple **canonical orbital element sets** and the
//! associated conversions between them:
//!
//! - [`keplerian_element`](crate::orbit_type::keplerian_element) — Classical Keplerian elements `(a, e, i, Ω, ω, M)`,
//!   valid for elliptic and hyperbolic orbits.
//! - [`equinoctial_element`](crate::orbit_type::equinoctial_element) — Equinoctial elements `(a, h, k, p, q, λ)`,
//!   a **non-singular formulation** well suited for orbit determination near
//!   zero eccentricity or inclination.
//! - [`cometary_element`](crate::orbit_type::cometary_element) — Perihelion-based representation `(q, e, i, Ω, ω, ν)`,
//!   convenient for parabolic and hyperbolic orbits.
//!
//! The [`OrbitalElements`](crate::orbit_type::OrbitalElements) enum acts as a **type-erased wrapper** that can hold
//! any of these three representations, while providing uniform constructors and
//! conversion methods.
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
//! // Convert to Keplerian form if possible
//! if let Ok(kep) = elems.to_keplerian() {
//!     if let OrbitalElements::Keplerian { elements, .. } = kep {
//!         println!("semi-major axis = {}", elements.semi_major_axis);
//!     }
//! }
//! ```
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

    /// Convert to Keplerian representation, propagating covariance if present.
    ///
    /// Return
    /// ------
    /// * `Ok(OrbitalElements::Keplerian)` – Converted elements with propagated
    ///   uncertainty and covariance when available.
    /// * `Err(OutfitError)` – If the conversion is not defined for the current
    ///   element set (e.g. parabolic cometary elements).
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

    /// Convert to equinoctial representation, propagating covariance if present.
    ///
    /// Return
    /// ------
    /// * `Ok(OrbitalElements::Equinoctial)` – Converted elements with propagated
    ///   uncertainty and covariance when available.
    /// * `Err(OutfitError)` – If the conversion is not defined for the current
    ///   element set (e.g. hyperbolic cometary elements).
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
