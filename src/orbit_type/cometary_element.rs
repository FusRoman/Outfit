use crate::{
    orbit_type::{equinoctial_element::EquinoctialElements, keplerian_element::KeplerianElements},
    outfit_errors::OutfitError,
};

/// Elements in "cometary" form used for parabolic/hyperbolic solutions.
/// Units: q in AU; angles in radians.
#[derive(Debug, Clone, PartialEq)]
pub struct CometaryElements {
    pub reference_epoch: f64,          // MJD
    pub perihelion_distance: f64,      // q (AU)
    pub eccentricity: f64,             // e (>= 1 for hyperbolic; =1 parabolic)
    pub inclination: f64,              // i (rad)
    pub ascending_node_longitude: f64, // Ω (rad)
    pub periapsis_argument: f64,       // ω (rad)
    pub true_anomaly: f64,             // ν (rad)
}

impl CometaryElements {
    /// Compute the hyperbolic mean anomaly from cometary elements (e > 1).
    ///
    /// Arguments
    /// -----------------
    /// * `q` – perihelion distance (AU)
    /// * `e` – eccentricity (> 1)
    /// * `nu` – true anomaly ν (radians)
    ///
    /// Return
    /// ----------
    /// * Hyperbolic mean anomaly `M` in radians.
    ///
    /// See also
    /// ------------
    /// * [`CometaryElements`] – Representation of hyperbolic/parabolic elements.
    /// * [`KeplerianElements`] – Classical orbital elements.
    ///
    /// Reference:
    /// - Danby, "Fundamentals of Celestial Mechanics"
    pub fn hyperbolic_mean_anomaly(e: f64, nu: f64) -> f64 {
        assert!(e > 1.0, "eccentricity must be > 1 for hyperbolic orbits");

        // Compute hyperbolic anomaly H from true anomaly ν
        let tanh_half_h = ((e - 1.0) / (e + 1.0)).sqrt() * (nu / 2.0).tan();
        let h = 2.0 * tanh_half_h.atanh();

        // Hyperbolic mean anomaly M = e sinh H - H
        e * h.sinh() - h
    }

    /// Convert Cometary → Keplerian (hyperbolic/parabolic).
    ///
    /// Note
    /// ----------
    /// * For e = 1 (parabolic), Keplerian "a" is formally infinite; this function
    ///   returns an error. For e > 1, `a` is negative and we interpret `M` as the
    ///   hyperbolic mean anomaly when/if needed by downstream code.
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

    fn try_from(ce: CometaryElements) -> Result<Self, Self::Error> {
        ce.cometary_to_keplerian()
    }
}

impl TryFrom<&CometaryElements> for KeplerianElements {
    type Error = OutfitError;

    fn try_from(ce: &CometaryElements) -> Result<Self, Self::Error> {
        ce.cometary_to_keplerian()
    }
}

impl TryFrom<CometaryElements> for EquinoctialElements {
    type Error = OutfitError;

    fn try_from(ce: CometaryElements) -> Result<Self, Self::Error> {
        let ke = ce.cometary_to_keplerian()?;
        Ok(EquinoctialElements::from(ke))
    }
}

impl TryFrom<&CometaryElements> for EquinoctialElements {
    type Error = OutfitError;

    fn try_from(ce: &CometaryElements) -> Result<Self, Self::Error> {
        let ke = ce.cometary_to_keplerian()?;
        Ok(EquinoctialElements::from(ke))
    }
}

use std::fmt;

impl fmt::Display for CometaryElements {
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
