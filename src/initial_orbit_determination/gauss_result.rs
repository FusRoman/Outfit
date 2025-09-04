//! # Gauss orbit determination result
//!
//! This module defines [`GaussResult`], an enum representing the outcome of the
//! **Gauss initial orbit determination method** applied to a triplet of astrometric
//! observations.  
//! It encapsulates the orbital elements computed from the Gauss method and distinguishes
//! between the *preliminary* and *corrected* stages of the algorithm.
//!
//! ## Variants
//!
//! - **`PrelimOrbit`**  
//!   A preliminary orbit computed directly from the root of the 8th-degree polynomial
//!   and the Gibbs method, without any iterative refinement.
//!
//! - **`CorrectedOrbit`**  
//!   A corrected orbit obtained after at least one iteration of refinement, typically
//!   using velocity correction or Lagrange coefficient adjustments.
//!
//! Both variants wrap an [`OrbitalElements`] value, which can hold Keplerian,
//! Equinoctial, or Cometary elements depending on the solution domain.
//!
//! ## Features
//!
//! - Query methods:
//!   * [`GaussResult::is_prelim`] – check if the result is a preliminary orbit.
//!   * [`GaussResult::is_corrected`] – check if the result includes refinement.
//! - Accessors:
//!   * [`GaussResult::get_orbit`] – borrow the inner [`OrbitalElements`].
//!   * [`GaussResult::as_inner`] – alias to `get_orbit`, for generic contexts.
//!   * [`GaussResult::into_inner`] – consume the enum and return the inner [`OrbitalElements`].
//!
//! Additionally, `GaussResult` implements [`Display`](std::fmt::Display) to provide a formatted,
//! human-readable view of the result and its orbital elements.
//!
//! ## Usage
//!
//! This type is returned by functions such as
//! [`GaussObs::prelim_orbit`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit)
//! or [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit).
//!
//! ```rust,no_run
//! use outfit::initial_orbit_determination::gauss_result::GaussResult;
//!
//! fn handle_result(result: GaussResult) {
//!     if result.is_corrected() {
//!         println!("Corrected orbit:\n{result}");
//!     } else {
//!         println!("Preliminary orbit:\n{result}");
//!     }
//! }
//! ```
//!
//! ## Notes
//!
//! - If you need a specific representation, use [`OrbitalElements::to_keplerian`] or
//!   [`OrbitalElements::to_equinoctial`] for explicit conversion.
//!
//! ## See also
//!
//! - [`OrbitalElements`] – Canonical orbital element sum type.
//! - [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements), [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements), [`CometaryElements`](crate::orbit_type::cometary_element::CometaryElements) – The three supported element forms.
//! - [`GaussObs`](crate::initial_orbit_determination::gauss::GaussObs)
//! - Milani & Gronchi (2010), *Theory of Orbit Determination*

use crate::orbit_type::OrbitalElements;
use std::fmt;

/// Result of the Gauss initial orbit determination method.
///
/// This enum represents the possible outcomes of a Gauss-based orbit estimation
/// from a triplet of astrometric observations. It distinguishes between:
///
/// - A preliminary orbit derived directly from the polynomial root and position/velocity reconstruction,
/// - A corrected orbit obtained after iterative refinement of the velocity and positions.
///
/// Variants
/// ---------
/// * `PrelimOrbit` – Orbit computed without correction; based solely on the root of the 8th-degree polynomial and
///   the initial Gibbs velocity estimate.
/// * `CorrectedOrbit` – Orbit obtained after convergence of the iterative refinement loop (typically involving
///   Lagrange coefficient fitting or velocity optimization).
///
/// Notes
/// -------
/// * Both variants wrap an [`OrbitalElements`] value, which may itself contain:
///   - [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements) (classical orbital elements),
///   - [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements) (non-singular formulation), or
///   - [`CometaryElements`](crate::orbit_type::cometary_element::CometaryElements) (perihelion form for e ≥ 1).
/// * Use [`GaussResult::get_orbit`] or [`GaussResult::as_inner`] to borrow the inner
///   [`OrbitalElements`] without matching on the enum.
/// * If a specific representation is needed, use conversion helpers on [`OrbitalElements`]
///   such as [`OrbitalElements::to_keplerian`] or [`OrbitalElements::to_equinoctial`].
///
/// # See also
/// * [`OrbitalElements`] – Sum type for canonical orbital elements.
/// * [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements), [`EquinoctialElements`](crate::orbit_type::equinoctial_element::EquinoctialElements), [`CometaryElements`](crate::orbit_type::cometary_element::CometaryElements) – Supported element forms.
/// * [`GaussObs::prelim_orbit`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit) – Main entry point that returns a `GaussResult`.
#[derive(PartialEq, Clone, Debug)]
pub enum GaussResult {
    PrelimOrbit(OrbitalElements),
    CorrectedOrbit(OrbitalElements),
}

impl GaussResult {
    /// Check whether the result is a preliminary (uncorrected) orbit.
    ///
    /// Arguments
    /// -----------------
    /// * None
    ///
    /// Return
    /// ----------
    /// * `true` if the variant is [`GaussResult::PrelimOrbit`], `false` otherwise.
    pub fn is_prelim(&self) -> bool {
        matches!(self, GaussResult::PrelimOrbit(_))
    }

    /// Check whether the result is a corrected orbit after refinement.
    ///
    /// Arguments
    /// -----------------
    /// * None
    ///
    /// Return
    /// ----------
    /// * `true` if the variant is [`GaussResult::CorrectedOrbit`], `false` otherwise.
    pub fn is_corrected(&self) -> bool {
        matches!(self, GaussResult::CorrectedOrbit(_))
    }

    /// Borrow the orbital elements associated with this Gauss result.
    ///
    /// The returned elements are kept in their native representation as produced
    /// by the solver (i.e., [`OrbitalElements::Keplerian`], [`OrbitalElements::Equinoctial`],
    /// or [`OrbitalElements::Cometary`]). Use conversion helpers on
    /// [`OrbitalElements`] if you need a specific parameterization.
    ///
    /// Arguments
    /// -----------------
    /// * None
    ///
    /// Return
    /// ----------
    /// * `&OrbitalElements` – A shared reference to the inner orbital elements.
    ///
    /// See also
    /// ------------
    /// * [`OrbitalElements::to_keplerian`] – Convert to classical Keplerian elements.
    /// * [`OrbitalElements::to_equinoctial`] – Convert to non-singular equinoctial elements.
    pub fn get_orbit(&self) -> &OrbitalElements {
        match self {
            GaussResult::PrelimOrbit(orbit) => orbit,
            GaussResult::CorrectedOrbit(orbit) => orbit,
        }
    }

    /// Alias for [`GaussResult::get_orbit`].
    ///
    /// Provided for naming consistency in generic contexts where `as_inner`
    /// better communicates the intent of borrowing the wrapped value.
    ///
    /// Arguments
    /// -----------------
    /// * None
    ///
    /// Return
    /// ----------
    /// * `&OrbitalElements` – A shared reference to the inner orbital elements.
    ///
    /// See also
    /// ------------
    /// * [`GaussResult::get_orbit`] – Primary accessor.
    pub fn as_inner(&self) -> &OrbitalElements {
        self.get_orbit()
    }

    /// Consume the enum and return the orbital elements by value.
    ///
    /// This method transfers ownership of the contained [`OrbitalElements`],
    /// allowing you to store or transform it without borrowing constraints.
    ///
    /// Arguments
    /// -----------------
    /// * `self` – The [`GaussResult`] to be consumed.
    ///
    /// Return
    /// ----------
    /// * `OrbitalElements` – The inner orbital elements, moved out of `self`.
    ///
    /// See also
    /// ------------
    /// * [`GaussResult::get_orbit`] – Borrow instead of moving.
    /// * [`OrbitalElements::to_keplerian`] – Convert to Keplerian if needed.
    /// * [`OrbitalElements::to_equinoctial`] – Convert to equinoctial if needed.
    pub fn into_inner(self) -> OrbitalElements {
        match self {
            GaussResult::PrelimOrbit(orbit) => orbit,
            GaussResult::CorrectedOrbit(orbit) => orbit,
        }
    }
}

impl fmt::Display for GaussResult {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            GaussResult::PrelimOrbit(orbit) => {
                writeln!(f, "Gauss IOD Result: Preliminary Orbit")?;
                write!(f, "{orbit}")
            }
            GaussResult::CorrectedOrbit(orbit) => {
                writeln!(f, "Gauss IOD Result: Corrected Orbit")?;
                write!(f, "{orbit}")
            }
        }
    }
}

#[cfg(test)]
mod gauss_results_tests {
    use crate::orbit_type::keplerian_element::KeplerianElements;

    use super::*;
    use std::format;

    fn dummy_orbit() -> OrbitalElements {
        OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 59000.0,
            semi_major_axis: 2.5,
            eccentricity: 0.12,
            inclination: 0.1,
            ascending_node_longitude: 1.2,
            periapsis_argument: 0.8,
            mean_anomaly: 0.3,
        })
    }

    #[test]
    fn test_variant_identification() {
        let prelim = GaussResult::PrelimOrbit(dummy_orbit());
        let corrected = GaussResult::CorrectedOrbit(dummy_orbit());

        assert!(prelim.is_prelim());
        assert!(!prelim.is_corrected());

        assert!(corrected.is_corrected());
        assert!(!corrected.is_prelim());
    }

    #[test]
    fn test_get_orbit_consistency() {
        let orbit = dummy_orbit();
        let result = GaussResult::CorrectedOrbit(orbit.clone());
        assert_eq!(result.get_orbit(), &orbit);
        assert_eq!(result.as_inner(), &orbit);
    }

    #[test]
    fn test_into_inner_extracts_correct_value() {
        let orbit = dummy_orbit();
        let result = GaussResult::PrelimOrbit(orbit.clone());
        let extracted = result.into_inner();
        assert_eq!(extracted, orbit);
    }

    #[test]
    fn test_display_format_summary() {
        let orbit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 59001.5,
            semi_major_axis: 1.234567,
            eccentricity: 0.1,
            inclination: 0.2,
            ascending_node_longitude: 0.3,
            periapsis_argument: 0.4,
            mean_anomaly: 0.5,
        });

        let result = GaussResult::CorrectedOrbit(orbit);

        let output = format!("{result}");

        // Vérifie que le type est bien annoncé
        assert!(output.starts_with("Gauss IOD Result: Corrected Orbit"));
        // Vérifie quelques champs clés (avec la nouvelle mise en forme)
        assert!(output.contains("a   (semi-major axis)       = 1.234567 AU"));
        assert!(output.contains("e   (eccentricity)          = 0.100000"));
        assert!(output.contains("i   (inclination)           = 0.200000 rad (11.459156°)"));
        assert!(output.contains("Keplerian Elements @ epoch (MJD): 59001.500000"));
    }
}
