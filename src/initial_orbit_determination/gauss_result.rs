use crate::keplerian_element::KeplerianElements;
use std::fmt;
use std::ops::Deref;

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
/// * Both variants wrap a [`KeplerianElements`] object, which contains the classical orbital parameters.
/// * This type implements [`Deref`] to `KeplerianElements`, allowing direct access to orbital fields
///   using dot notation (e.g., `result.mean_anomaly`). While idiomatic in Rust, this behavior may be
///   surprising in some contexts. Use [`get_orbit`] if you prefer to make the access explicit.
///
/// # See also
/// * [`KeplerianElements`](crate::keplerian_element::KeplerianElements) – Struct representing the orbital elements.
/// * [`GaussObs::prelim_orbit`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit) – Main entry point that returns a `GaussResult`.
/// * [`Deref`](https://doc.rust-lang.org/std/ops/trait.Deref.html) – trait allowing ergonomic field access.
#[derive(PartialEq, Clone, Debug)]
pub enum GaussResult {
    PrelimOrbit(KeplerianElements),
    CorrectedOrbit(KeplerianElements),
}

impl GaussResult {
    /// Check whether the result is a preliminary (uncorrected) orbit.
    ///
    /// Returns `true` if the orbit was computed directly from the root of the 8th-degree polynomial
    /// without any iterative correction.
    ///
    /// Returns
    /// --------
    /// * `true` if the variant is [`GaussResult::PrelimOrbit`],
    /// * `false` otherwise.
    pub fn is_prelim(&self) -> bool {
        matches!(self, GaussResult::PrelimOrbit(_))
    }

    /// Check whether the result is a corrected orbit after refinement.
    ///
    /// Returns `true` if the orbit underwent at least one successful iteration of correction
    /// (typically via velocity or position refinement).
    ///
    /// Returns
    /// --------
    /// * `true` if the variant is [`GaussResult::CorrectedOrbit`],
    /// * `false` otherwise.
    pub fn is_corrected(&self) -> bool {
        matches!(self, GaussResult::CorrectedOrbit(_))
    }

    /// Retrieve the Keplerian orbital elements associated with the result.
    ///
    /// Returns a reference to the [`KeplerianElements`] struct regardless of whether the orbit is preliminary or corrected.
    ///
    /// Returns
    /// --------
    /// * `&KeplerianElements` – orbital elements associated with the result.
    ///
    /// Notes
    /// ------
    /// This method provides a unified interface for accessing the orbital solution without needing
    /// to match explicitly on the enum variant.
    pub fn get_orbit(&self) -> &KeplerianElements {
        match self {
            GaussResult::PrelimOrbit(orbit) => orbit,
            GaussResult::CorrectedOrbit(orbit) => orbit,
        }
    }

    /// Borrow the inner [`KeplerianElements`] struct immutably.
    ///
    /// This is equivalent to [`get_orbit`], and provided for naming consistency
    /// when used in generic contexts or with other enums that implement `as_inner`.
    ///
    /// Returns
    /// --------
    /// * `&KeplerianElements` – a shared reference to the orbital elements.
    pub fn as_inner(&self) -> &KeplerianElements {
        self.get_orbit()
    }

    /// Consume the enum and extract the inner [`KeplerianElements`] by value.
    ///
    /// This method moves the contents out of the enum, allowing ownership
    /// of the orbital solution.
    ///
    /// Returns
    /// --------
    /// * `KeplerianElements` – the inner orbital elements, consuming `self`.
    pub fn into_inner(self) -> KeplerianElements {
        match self {
            GaussResult::PrelimOrbit(orbit) => orbit,
            GaussResult::CorrectedOrbit(orbit) => orbit,
        }
    }
}

impl Deref for GaussResult {
    type Target = KeplerianElements;

    /// Deref implementation for `GaussResult`.
    ///
    /// This allows direct field access to the underlying [`KeplerianElements`] fields
    /// using dot notation, e.g. `gauss_result.mean_anomaly`.
    ///
    /// Returns
    /// --------
    /// * A shared reference to the inner `KeplerianElements` regardless of variant.
    fn deref(&self) -> &Self::Target {
        self.get_orbit()
    }
}

impl fmt::Display for GaussResult {
    /// Pretty-print the orbital elements contained in a GaussResult.
    ///
    /// Outputs all six Keplerian elements with angles converted to degrees,
    /// plus the reference epoch (MJD). The output is formatted for CLI or logs.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let label = match self {
            GaussResult::PrelimOrbit(_) => "[Prelim orbit]",
            GaussResult::CorrectedOrbit(_) => "[Corrected orbit]",
        };

        let orb = self.get_orbit();

        // Convert all angular quantities to degrees
        let i_deg = orb.inclination.to_degrees();
        let omega_deg = orb.ascending_node_longitude.to_degrees();
        let argp_deg = orb.periapsis_argument.to_degrees();
        let m_deg = orb.mean_anomaly.to_degrees();

        writeln!(f, "{}", label)?;
        writeln!(f, "  Epoch (MJD): {:.5}", orb.reference_epoch)?;
        writeln!(f, "  a   (AU)   : {:.6}", orb.semi_major_axis)?;
        writeln!(f, "  e           : {:.6}", orb.eccentricity)?;
        writeln!(f, "  i   (deg)  : {:.6}", i_deg)?;
        writeln!(f, "  Ω   (deg)  : {:.6}", omega_deg)?;
        writeln!(f, "  ω   (deg)  : {:.6}", argp_deg)?;
        write!(f, "  M   (deg)  : {:.6}", m_deg)
    }
}

#[cfg(test)]
mod gauss_results_tests {
    use super::*;
    use std::format;

    fn dummy_orbit() -> KeplerianElements {
        KeplerianElements {
            reference_epoch: 59000.0,
            semi_major_axis: 2.5,
            eccentricity: 0.12,
            inclination: 0.1,
            ascending_node_longitude: 1.2,
            periapsis_argument: 0.8,
            mean_anomaly: 0.3,
        }
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
        assert_eq!(*result, orbit); // test Deref
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
        let orbit = KeplerianElements {
            reference_epoch: 59001.5,
            semi_major_axis: 1.234567,
            eccentricity: 0.1,
            inclination: 0.2,
            ascending_node_longitude: 0.3,
            periapsis_argument: 0.4,
            mean_anomaly: 0.5,
        };

        let result = GaussResult::CorrectedOrbit(orbit);

        let output = format!("{}", result);

        assert!(output.starts_with("[Corrected orbit]"));
        assert!(output.contains("a   (AU)   : 1.234567"));
        assert!(output.contains("e           : 0.100000"));
        assert!(output.contains("i   (deg)  : 11.459156"));
        assert!(output.contains("Epoch (MJD): 59001.50000"));
    }
}
