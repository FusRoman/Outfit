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
    /// Compute the **hyperbolic mean anomaly** `M` from cometary elements (`e > 1`).
    ///
    /// This uses the standard relation:
    /// 1) Convert true anomaly `ν` to hyperbolic anomaly `H` via  
    ///    `tanh(H/2) = sqrt((e-1)/(e+1)) * tan(ν/2)`  
    /// 2) Then `M = e*sinh(H) - H`.
    ///
    /// Arguments
    /// -----------------
    /// * `e`  – Eccentricity (> 1).  
    /// * `nu` – True anomaly `ν` (radians).
    ///
    /// Return
    /// ----------
    /// * Hyperbolic mean anomaly `M` (radians).
    ///
    /// Panics
    /// ----------
    /// * If `e ≤ 1` (domain violation).
    ///
    /// See also
    /// ------------
    /// * [`CometaryElements`] – Perihelion-based orbital elements.  
    /// * [`KeplerianElements`] – Classical elements used as an intermediate in conversions.
    ///
    /// Reference
    /// ---------
    /// * Danby, *Fundamentals of Celestial Mechanics*, hyperbolic motion relations.
    pub fn hyperbolic_mean_anomaly(e: f64, nu: f64) -> f64 {
        assert!(e > 1.0, "eccentricity must be > 1 for hyperbolic orbits");

        // Compute hyperbolic anomaly H from true anomaly ν
        let tanh_half_h = ((e - 1.0) / (e + 1.0)).sqrt() * (nu / 2.0).tan();
        let h = 2.0 * tanh_half_h.atanh();

        // Hyperbolic mean anomaly M = e sinh H - H
        e * h.sinh() - h
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

        let m = CometaryElements::hyperbolic_mean_anomaly(e, self.true_anomaly);

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
