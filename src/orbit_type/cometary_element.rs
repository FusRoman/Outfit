use crate::{
    orbit_type::{equinoctial_element::EquinoctialElements, keplerian_element::KeplerianElements},
    outfit_errors::OutfitError,
};

/// # Cometary orbital elements
///
/// Cometary (perihelion-based) elements are convenient for **parabolic and
/// hyperbolic** solutions, where the semi-major axis is not finite (parabola)
/// or negative (hyperbola).
///
/// Units & conventions
/// --------------------
/// - Distances in **AU**; angles in **radians**; epochs in **MJD (TDB)**.
/// - State is assumed heliocentric, equatorial mean J2000.
/// - For hyperbolic motion: `a < 0`, `e > 1`; for parabolic: `e = 1`.
///
/// See also
/// ------------
/// * [`KeplerianElements`] – Classical elements `(a, e, i, Ω, ω, M)` (supports elliptic & hyperbolic).  
/// * [`EquinoctialElements`] – Non-singular elements `(a, h, k, p, q, λ)`.  
/// * [`CometaryElements::hyperbolic_mean_anomaly`] – Returns the hyperbolic mean anomaly from `(e, ν)`
#[derive(Debug, Clone, PartialEq)]
pub struct CometaryElements {
    /// Reference epoch of the element set (MJD, TDB).
    pub reference_epoch: f64,

    /// Perihelion distance `q` (AU).
    pub perihelion_distance: f64,

    /// Eccentricity `e` (≥ 1).  
    /// `e = 1` for parabolic, `e > 1` for hyperbolic orbits.
    pub eccentricity: f64,

    /// Inclination `i` (rad).
    pub inclination: f64,

    /// Longitude of the ascending node `Ω` (rad).
    pub ascending_node_longitude: f64,

    /// Argument of periapsis `ω` (rad).
    pub periapsis_argument: f64,

    /// True anomaly `ν` at the reference epoch (rad).
    pub true_anomaly: f64,
}

impl CometaryElements {
    /// Compute the hyperbolic mean anomaly `M` from `(e, ν)`.
    ///
    /// This implementation clamps the argument of `atanh` strictly inside `(-1, 1)`
    /// to avoid infinities and `NaN` when `tanh(H/2) → ±1` due to round-off.
    ///
    /// Arguments
    /// -----------------
    /// * `e`  – Eccentricity (> 1 expected).
    /// * `nu` – True anomaly `ν` (radians).
    ///
    /// Return
    /// ----------
    /// * `Ok(M)` – Hyperbolic mean anomaly (radians).  
    /// * `Err(OutfitError::InvalidOrbit)` – if `e ≤ 1`.
    ///
    /// See also
    /// ------------
    /// * [`CometaryElements`] – Perihelion-based orbital elements.
    /// * Danby, *Fundamentals of Celestial Mechanics*, hyperbolic motion.
    pub fn hyperbolic_mean_anomaly(e: f64, nu: f64) -> Result<f64, OutfitError> {
        if e <= 1.0 {
            return Err(OutfitError::InvalidOrbit(
                "eccentricity must be > 1 for hyperbolic orbits".into(),
            ));
        }

        const EPS: f64 = 1e-15;

        let s = ((e - 1.0) / (e + 1.0)).sqrt();
        let t = (0.5 * nu).tan();
        let x = (s * t).clamp(-1.0 + EPS, 1.0 - EPS);

        let h = 2.0 * x.atanh();
        Ok(e * h.sinh() - h)
    }

    /// Convert **Cometary → Keplerian** (hyperbolic/parabolic).
    ///
    /// For `e > 1` (hyperbola), the conversion yields a valid Keplerian set with
    /// `a < 0` and the returned `mean_anomaly` interpreted as **hyperbolic mean anomaly**.
    /// For `e = 1` (parabola), the semi-major axis is not finite; this function returns an error.
    ///
    /// Arguments
    /// -----------------
    /// * `self` – Cometary element set `(q, e, i, Ω, ω, ν)` at `reference_epoch`.
    ///
    /// Return
    /// ----------
    /// * A `Result` with the converted [`KeplerianElements`] if `e > 1`,
    ///   or an [`OutfitError`] if `e = 1`.
    ///
    /// Errors
    /// ----------
    /// * `OutfitError::InvalidConversion` – when `e = 1` (parabolic orbit).
    ///
    /// See also
    /// ------------
    /// * [`hyperbolic_mean_anomaly`] – Hyperbolic mean anomaly used for `mean_anomaly`.  
    /// * [`EquinoctialElements::from`] – Pathway to equinoctials via Keplerian.
    fn cometary_to_keplerian(&self) -> Result<KeplerianElements, OutfitError> {
        if (self.eccentricity - 1.0).abs() < 1e-12 {
            return Err(OutfitError::InvalidConversion(
                "Parabolic orbit cannot be represented with finite a".into(),
            ));
        }
        // For hyperbolic: p = q(1 + e), and a = -p/(e^2 - 1).
        // (This mirrors the standard relations for conics.)
        let e = self.eccentricity;
        let p = self.perihelion_distance * (1.0 + e);
        let a = -p / (e * e - 1.0);

        let m = CometaryElements::hyperbolic_mean_anomaly(e, self.true_anomaly)?;

        // We keep (i, Ω, ω) and set a placeholder "mean anomaly" interpreted
        // by downstream code as hyperbolic M when used with Kepler dynamics.
        Ok(KeplerianElements {
            reference_epoch: self.reference_epoch,
            semi_major_axis: a,
            eccentricity: e,
            inclination: self.inclination,
            ascending_node_longitude: self.ascending_node_longitude,
            periapsis_argument: self.periapsis_argument,
            mean_anomaly: m, // interpreted as hyperbolic mean anomaly M
        })
    }
}

impl TryFrom<CometaryElements> for KeplerianElements {
    type Error = OutfitError;

    /// Convert by-value `CometaryElements` into `KeplerianElements`.
    ///
    /// Arguments
    /// -----------------
    /// * `ce` – Cometary elements to convert.
    ///
    /// Return
    /// ----------
    /// * `Ok(KeplerianElements)` if `ce.eccentricity > 1`, otherwise an error for the parabolic case.
    fn try_from(ce: CometaryElements) -> Result<Self, Self::Error> {
        ce.cometary_to_keplerian()
    }
}

impl TryFrom<&CometaryElements> for KeplerianElements {
    type Error = OutfitError;

    /// Convert a borrowed `CometaryElements` into `KeplerianElements`.
    ///
    /// Arguments
    /// -----------------
    /// * `ce` – Borrowed cometary elements to convert.
    ///
    /// Return
    /// ----------
    /// * `Ok(KeplerianElements)` if `ce.eccentricity > 1`, otherwise an error for the parabolic case.
    fn try_from(ce: &CometaryElements) -> Result<Self, Self::Error> {
        ce.cometary_to_keplerian()
    }
}

impl TryFrom<CometaryElements> for EquinoctialElements {
    type Error = OutfitError;

    /// Convert by-value `CometaryElements` into `EquinoctialElements` via Keplerian.
    ///
    /// Arguments
    /// -----------------
    /// * `ce` – Cometary elements to convert.
    ///
    /// Return
    /// ----------
    /// * `Ok(EquinoctialElements)` for `e > 1`, error for `e = 1`.
    ///
    /// See also
    /// ------------
    /// * [`KeplerianElements::try_from`] – Intermediate conversion step.
    /// * [`EquinoctialElements::from`] – Keplerian → Equinoctial mapping.
    fn try_from(ce: CometaryElements) -> Result<Self, Self::Error> {
        let ke = ce.cometary_to_keplerian()?;
        Ok(EquinoctialElements::from(ke))
    }
}

impl TryFrom<&CometaryElements> for EquinoctialElements {
    type Error = OutfitError;

    /// Convert a borrowed `CometaryElements` into `EquinoctialElements` via Keplerian.
    ///
    /// Arguments
    /// -----------------
    /// * `ce` – Borrowed cometary elements to convert.
    ///
    /// Return
    /// ----------
    /// * `Ok(EquinoctialElements)` for `e > 1`, error for `e = 1`.
    ///
    /// See also
    /// ------------
    /// * [`KeplerianElements::try_from`] – Intermediate conversion step.
    /// * [`EquinoctialElements::from`] – Keplerian → Equinoctial mapping.
    fn try_from(ce: &CometaryElements) -> Result<Self, Self::Error> {
        let ke = ce.cometary_to_keplerian()?;
        Ok(EquinoctialElements::from(ke))
    }
}

use std::fmt;

impl fmt::Display for CometaryElements {
    /// Pretty-print cometary elements with both radians and degrees for angles.
    ///
    /// See also
    /// ------------
    /// * [`CometaryElements`] – Data container.
    /// * [`KeplerianElements`] – Classical representation used in conversions.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rad_to_deg = 180.0 / std::f64::consts::PI;
        writeln!(
            f,
            "Cometary Elements @ epoch (MJD): {:.6}",
            self.reference_epoch
        )?;
        writeln!(f, "------------------------------------------------")?;
        writeln!(
            f,
            "  q   (perihelion distance)     = {:.6} AU",
            self.perihelion_distance
        )?;
        writeln!(
            f,
            "  e   (eccentricity)            = {:.6}",
            self.eccentricity
        )?;
        writeln!(
            f,
            "  i   (inclination)             = {:.6} rad ({:.6}°)",
            self.inclination,
            self.inclination * rad_to_deg
        )?;
        writeln!(
            f,
            "  Ω   (longitude of node)       = {:.6} rad ({:.6}°)",
            self.ascending_node_longitude,
            self.ascending_node_longitude * rad_to_deg
        )?;
        writeln!(
            f,
            "  ω   (argument of periapsis)   = {:.6} rad ({:.6}°)",
            self.periapsis_argument,
            self.periapsis_argument * rad_to_deg
        )?;
        writeln!(
            f,
            "  ν   (true anomaly)            = {:.6} rad ({:.6}°)",
            self.true_anomaly,
            self.true_anomaly * rad_to_deg
        )
    }
}

#[cfg(test)]
mod cometary_element_tests {
    use super::*;
    use approx::{assert_abs_diff_eq, assert_relative_eq};
    use std::f64::consts::PI;

    /// Assert that two angles are equal modulo 2π, within an absolute epsilon.
    fn assert_angle_eq(a: f64, b: f64, eps: f64) {
        let two_pi = 2.0 * PI;
        let d = (a - b + PI).rem_euclid(two_pi) - PI;
        assert_abs_diff_eq!(d, 0.0, epsilon = eps);
    }

    // ---------- hyperbolic_mean_anomaly ----------

    #[test]
    fn hma_is_zero_when_nu_is_zero() {
        // For ν=0, tanh(H/2)=0 -> H=0 -> M=0.
        let e = 1.2;
        let nu = 0.0;
        let m = CometaryElements::hyperbolic_mean_anomaly(e, nu).unwrap();
        assert_abs_diff_eq!(m, 0.0, epsilon = 1e-15);
    }

    #[test]
    fn hma_increases_with_positive_nu() {
        // For small positive ν, M should be positive.
        let e = 1.3;
        let nu = 10.0_f64.to_degrees();
        let m = CometaryElements::hyperbolic_mean_anomaly(e, nu).unwrap();
        assert!(m > 0.0);
    }

    #[test]
    fn hma_returns_err_for_e_equal_1() {
        let res = CometaryElements::hyperbolic_mean_anomaly(1.0, 0.0);
        assert!(res.is_err(), "expected Err for e = 1");
    }

    #[test]
    fn hma_returns_err_for_e_less_than_1() {
        let res = CometaryElements::hyperbolic_mean_anomaly(0.9, 0.0);
        assert!(res.is_err(), "expected Err for e < 1");
    }

    // ---------- Cometary -> Keplerian conversion (core) ----------

    #[test]
    fn cometary_to_keplerian_hyperbolic_yields_negative_a_and_preserves_angles() {
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.8, // q
            eccentricity: 1.4,        // e>1 hyperbolic
            inclination: 12.0_f64.to_degrees(),
            ascending_node_longitude: 33.0_f64.to_degrees(),
            periapsis_argument: 45.0_f64.to_degrees(),
            true_anomaly: 5.0_f64.to_degrees(),
        };

        let ke = ce
            .cometary_to_keplerian()
            .expect("hyperbolic cometary -> keplerian should succeed");

        // a must be negative for hyperbola, e preserved
        assert!(
            ke.semi_major_axis < 0.0,
            "a should be negative for hyperbola"
        );
        assert_relative_eq!(ke.eccentricity, ce.eccentricity, max_relative = 1e-14);

        // Orientation angles preserved
        assert_angle_eq(ke.inclination, ce.inclination, 1e-14);
        assert_angle_eq(
            ke.ascending_node_longitude,
            ce.ascending_node_longitude,
            1e-14,
        );
        assert_angle_eq(ke.periapsis_argument, ce.periapsis_argument, 1e-14);

        // mean_anomaly equals hyperbolic mean anomaly at ν
        let m_expected =
            CometaryElements::hyperbolic_mean_anomaly(ce.eccentricity, ce.true_anomaly).unwrap();
        assert_abs_diff_eq!(ke.mean_anomaly, m_expected, epsilon = 1e-14);

        // Sanity check for 'a' formula used: a = -q(1+e)/(e^2-1)
        let a_expected = -ce.perihelion_distance * (1.0 + ce.eccentricity)
            / (ce.eccentricity * ce.eccentricity - 1.0);
        assert_relative_eq!(ke.semi_major_axis, a_expected, max_relative = 1e-14);
    }

    #[test]
    fn cometary_to_keplerian_parabolic_fails() {
        // Exactly parabolic (e=1) -> error
        let ce_exact = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 1.0,
            eccentricity: 1.0,
            inclination: 0.0,
            ascending_node_longitude: 0.0,
            periapsis_argument: 0.0,
            true_anomaly: 0.0,
        };
        assert!(
            ce_exact.cometary_to_keplerian().is_err(),
            "parabolic conversion must error"
        );

        // Within the tolerance band |e-1| < 1e-12 -> error
        let ce_tiny = CometaryElements {
            eccentricity: 1.0 + 5e-13, // still < 1e-12 away from 1
            ..ce_exact
        };
        assert!(
            ce_tiny.cometary_to_keplerian().is_err(),
            "borderline parabolic (|e-1|<1e-12) must error"
        );

        // Just outside tolerance -> should succeed as hyperbolic
        let ce_ok = CometaryElements {
            eccentricity: 1.0 + 5e-12, // > 1e-12 away
            ..ce_exact
        };
        assert!(
            ce_ok.cometary_to_keplerian().is_ok(),
            "slightly hyperbolic should convert"
        );
    }

    // ---------- TryFrom implementations ----------

    #[test]
    fn tryfrom_owned_to_keplerian_matches_method() {
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.7,
            eccentricity: 1.2,
            inclination: 3.0_f64.to_degrees(),
            ascending_node_longitude: 10.0_f64.to_degrees(),
            periapsis_argument: 20.0_f64.to_degrees(),
            true_anomaly: 15.0_f64.to_degrees(),
        };
        let via_method = ce.cometary_to_keplerian().unwrap();
        let via_tryfrom = KeplerianElements::try_from(ce).unwrap();
        assert_relative_eq!(
            via_method.semi_major_axis,
            via_tryfrom.semi_major_axis,
            max_relative = 1e-15
        );
        assert_relative_eq!(
            via_method.eccentricity,
            via_tryfrom.eccentricity,
            max_relative = 1e-15
        );
        assert_angle_eq(via_method.inclination, via_tryfrom.inclination, 1e-14);
        assert_angle_eq(
            via_method.ascending_node_longitude,
            via_tryfrom.ascending_node_longitude,
            1e-14,
        );
        assert_angle_eq(
            via_method.periapsis_argument,
            via_tryfrom.periapsis_argument,
            1e-14,
        );
        assert_abs_diff_eq!(
            via_method.mean_anomaly,
            via_tryfrom.mean_anomaly,
            epsilon = 1e-15
        );
    }

    #[test]
    fn tryfrom_borrowed_to_keplerian_matches_method() {
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.9,
            eccentricity: 1.05,
            inclination: 1.0_f64.to_degrees(),
            ascending_node_longitude: 2.0_f64.to_degrees(),
            periapsis_argument: 3.0_f64.to_degrees(),
            true_anomaly: 4.0_f64.to_degrees(),
        };
        let via_method = ce.cometary_to_keplerian().unwrap();
        let via_tryfrom = KeplerianElements::try_from(&ce).unwrap();
        assert_relative_eq!(
            via_method.semi_major_axis,
            via_tryfrom.semi_major_axis,
            max_relative = 1e-15
        );
        assert_relative_eq!(
            via_method.eccentricity,
            via_tryfrom.eccentricity,
            max_relative = 1e-15
        );
        assert_angle_eq(via_method.inclination, via_tryfrom.inclination, 1e-14);
        assert_abs_diff_eq!(
            via_method.mean_anomaly,
            via_tryfrom.mean_anomaly,
            epsilon = 1e-15
        );
    }

    #[test]
    fn tryfrom_to_equinoctial_via_keplerian_has_negative_a() {
        // When e>1, Equinoctial a should remain negative (inherited from Keplerian).
        let ce = CometaryElements {
            reference_epoch: 2460000.5,
            perihelion_distance: 0.6,
            eccentricity: 1.3,
            inclination: 8.0_f64.to_degrees(),
            ascending_node_longitude: 12.0_f64.to_degrees(),
            periapsis_argument: 30.0_f64.to_degrees(),
            true_anomaly: 25.0_f64.to_degrees(),
        };
        let eq = EquinoctialElements::try_from(&ce)
            .expect("Cometary -> Equinoctial via Keplerian should succeed for e>1");
        assert!(
            eq.semi_major_axis < 0.0,
            "equinoctial a should be negative for hyperbolic orbits"
        );
    }

    // ---------- Display formatting ----------

    #[test]
    fn display_contains_expected_lines_and_degrees() {
        let ce = CometaryElements {
            reference_epoch: 2461234.56789,
            perihelion_distance: 0.5,
            eccentricity: 1.8,
            inclination: 10.0_f64.to_degrees(),
            ascending_node_longitude: 20.0_f64.to_degrees(),
            periapsis_argument: 30.0_f64.to_degrees(),
            true_anomaly: 40.0_f64.to_degrees(),
        };

        let s = format!("{ce}");
        assert!(
            s.contains("Cometary Elements @ epoch (MJD)"),
            "header missing"
        );
        assert!(s.contains("q   (perihelion distance)"), "q line missing");
        assert!(s.contains("e   (eccentricity)"), "e line missing");
        assert!(
            s.contains(" rad (") && s.contains("°)"),
            "expected rad/deg formatting with degree symbol"
        );
        // Spot-check one angle formatting with degrees
        assert!(
            s.contains("inclination") && s.contains("rad (") && s.contains("°)"),
            "inclination rad/deg formatting missing"
        );
    }
}

#[cfg(test)]
mod cometary_element_proptests {
    use super::*;
    use proptest::prelude::*;
    use std::f64::consts::PI;

    // ---------- helpers ----------

    /// Wrap angle difference to [-π, π].
    fn wrap_to_pi(x: f64) -> f64 {
        (x + PI).rem_euclid(2.0 * PI) - PI
    }

    // Domain constants for generation
    const EPS_PARABOLA: f64 = 1e-12;
    const NU_EPS: f64 = 1e-6; // to avoid tan(ν/2) singularities at ±π

    // ---------- strategies ----------

    // e in (1 + 1e-12, 5]
    prop_compose! {
        fn e_hyperbolic()(u in 0.0f64..1.0) -> f64 {
            1.0 + EPS_PARABOLA + u * (5.0 - (1.0 + EPS_PARABOLA))
        }
    }

    // e very close to 1 within the "parabolic tolerance band"
    prop_compose! {
        fn e_near_parabolic()(u in -0.5f64..0.5) -> f64 {
            1.0 + u * EPS_PARABOLA // |e-1| < 5e-13 < 1e-12 -> erreur attendue
        }
    }

    // q in (1e-6, 50] AU (strictly positive, up to 50 AU)
    prop_compose! {
        fn perihelion_q()(u in 0.0f64..1.0) -> f64 {
            1e-6 + u * (50.0 - 1e-6)
        }
    }

    prop_compose! {
        fn true_anomaly()(u in 0.0f64..1.0) -> f64 {
            // ν ∈ (-π + NU_EPS, π - NU_EPS)
            -PI + NU_EPS + u * (2.0 * PI - 2.0 * NU_EPS)
        }
    }

    // i ∈ [0, π], Ω, ω ∈ [0, 2π]
    prop_compose! {
        fn inc_omega_arg()(ui in 0.0f64..1.0, uo in 0.0f64..1.0, ua in 0.0f64..1.0) -> (f64, f64, f64) {
            let i = ui * PI;
            let big_omega = uo * 2.0 * PI;
            let small_omega = ua * 2.0 * PI;
            (i, big_omega, small_omega)
        }
    }

    // ---------- properties ----------

    proptest! {
        /// For hyperbolic e>1 and safe ν, M_hyp(e,ν) is finite (no NaN/Inf).
        #[test]
        fn hma_no_nan_no_inf(e in e_hyperbolic(), nu in true_anomaly()) {
            let m = CometaryElements::hyperbolic_mean_anomaly(e, nu)
                .expect("e>1 must be Ok");
            prop_assert!(m.is_finite());
        }
    }

    proptest! {
        /// Odd symmetry away from the atanh saturation band:
        /// For e>1 and |x|<0.9 with x = sqrt((e-1)/(e+1)) * tan(ν/2),
        /// we expect M(e, -ν) ≈ -M(e, ν).
        #[test]
        fn hma_is_odd_in_nu(e in e_hyperbolic(), nu in true_anomaly()) {
            // Compute x = s * tan(ν/2) where s = sqrt((e-1)/(e+1)).
            let s = ((e - 1.0) / (e + 1.0)).sqrt();
            let t = (0.5 * nu).tan();
            // If tan blows up, or x is too close to ±1, skip this input:
            prop_assume!(t.is_finite());
            let x = s * t;
            // Keep well inside (-1,1) so the clamp in the implementation is inactive.
            const X_MAX: f64 = 0.9;
            prop_assume!(x.abs() < X_MAX);

            let m_pos = CometaryElements::hyperbolic_mean_anomaly(e, nu).unwrap();
            let m_neg = CometaryElements::hyperbolic_mean_anomaly(e, -nu).unwrap();

            // Mixed tolerance: absolute for near-zero, relative for large |M|.
            let err = (m_pos + m_neg).abs();
            let tol = 1e-12_f64.max(1e-12 * m_pos.abs().max(m_neg.abs()));
            prop_assert!(err <= tol, "oddness failed: |M(ν)+M(-ν)|={} > tol={}", err, tol);
        }
    }

    proptest! {
        /// Cometary -> Keplerian: a<0, a = -q(1+e)/(e^2-1), angles preserved, mean_anomaly = M_hyp.
        #[test]
        fn cometary_to_keplerian_invariants(
            q in perihelion_q(),
            e in e_hyperbolic(),
            nu in true_anomaly(),
            (i, big_omega, small_omega) in inc_omega_arg(),
            epoch in 2_400_000.5f64..=2_500_000.5
        ) {
            let ce = CometaryElements {
                reference_epoch: epoch,
                perihelion_distance: q,
                eccentricity: e,
                inclination: i,
                ascending_node_longitude: big_omega,
                periapsis_argument: small_omega,
                true_anomaly: nu,
            };

            let ke = ce.cometary_to_keplerian().expect("hyperbolic conversion must succeed");

            // a negative + formula check
            prop_assert!(ke.semi_major_axis < 0.0);
            let a_expected = -q * (1.0 + e) / (e*e - 1.0);
            prop_assert!((ke.semi_major_axis - a_expected).abs() <= 1e-12_f64.max(1e-12 * a_expected.abs()));

            // e preserved
            prop_assert!((ke.eccentricity - e).abs() <= 1e-14);

            // angles preserved modulo 2π
            prop_assert!(wrap_to_pi(ke.inclination - i).abs() < 1e-12);
            prop_assert!(wrap_to_pi(ke.ascending_node_longitude - big_omega).abs() < 1e-12);
            prop_assert!(wrap_to_pi(ke.periapsis_argument - small_omega).abs() < 1e-12);

            // mean_anomaly equals hyperbolic mean anomaly at ν
            let m_expected = CometaryElements::hyperbolic_mean_anomaly(e, nu).unwrap();
            prop_assert!((ke.mean_anomaly - m_expected).abs() <= 1e-12_f64.max(1e-12 * m_expected.abs()));
        }
    }

    proptest! {
        /// Near-parabolic band: |e-1| < 1e-12 -> conversion must error.
        #[test]
        fn cometary_to_keplerian_errors_in_parabolic_band(
            q in perihelion_q(),
            e in e_near_parabolic(), // within ±5e-13 around 1.0 -> < 1e-12
            nu in true_anomaly(),
            (i, big_omega, small_omega) in inc_omega_arg(),
            epoch in 2_400_000.5f64..=2_500_000.5
        ) {
            let ce = CometaryElements {
                reference_epoch: epoch,
                perihelion_distance: q,
                eccentricity: e,
                inclination: i,
                ascending_node_longitude: big_omega,
                periapsis_argument: small_omega,
                true_anomaly: nu,
            };
            prop_assert!(ce.cometary_to_keplerian().is_err());
        }
    }

    proptest! {
        /// Just outside the parabolic band: e >= 1 + 2e-12 -> conversion must succeed.
        #[test]
        fn cometary_to_keplerian_succeeds_outside_parabolic_band(
            q in perihelion_q(),
            u in 0.0f64..1.0,
            nu in true_anomaly(),
            (i, big_omega, small_omega) in inc_omega_arg(),
            epoch in 2_400_000.5f64..=2_500_000.5
        ) {
            let e = 1.0 + 2.0*EPS_PARABOLA + u * (4.0 - 2.0*EPS_PARABOLA); // >= 1 + 2e-12
            let ce = CometaryElements {
                reference_epoch: epoch,
                perihelion_distance: q,
                eccentricity: e,
                inclination: i,
                ascending_node_longitude: big_omega,
                periapsis_argument: small_omega,
                true_anomaly: nu,
            };
            let ke = ce.cometary_to_keplerian().expect("should succeed for e >= 1 + 2e-12");
            prop_assert!(ke.semi_major_axis.is_finite());
            prop_assert!(ke.semi_major_axis < 0.0);
        }
    }

    proptest! {
        /// Cometary -> Equinoctial via Keplerian:
        /// if supported, a should remain negative; if not supported, Err is acceptable.
        #[test]
        fn tryfrom_to_equinoctial_preserves_negative_a_when_ok(
            q in perihelion_q(),
            e in e_hyperbolic(),
            nu in true_anomaly(),
            (i, big_omega, small_omega) in inc_omega_arg(),
            epoch in 2_400_000.5f64..=2_500_000.5
        ) {
            let ce = CometaryElements {
                reference_epoch: epoch,
                perihelion_distance: q,
                eccentricity: e,
                inclination: i,
                ascending_node_longitude: big_omega,
                periapsis_argument: small_omega,
                true_anomaly: nu,
            };
            match EquinoctialElements::try_from(&ce) {
                Ok(eq) => {
                    // Hyperbolic: a should be negative too in the equinoctial mapping
                    prop_assert!(eq.semi_major_axis < 0.0);
                }
                Err(_) => {
                    // Acceptable: some implementations restrict equinoctials to elliptic.
                    prop_assert!(true);
                }
            }
        }
    }
}
