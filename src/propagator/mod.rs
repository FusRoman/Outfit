//! Orbit propagation strategies for differential orbit correction.
//!
//! This module exposes two propagator kinds:
//!
//! - `PropagatorKind::TwoBody` – classic Keplerian propagation through the
//!   analytic two-body solution.  This is the default and requires no external
//!   ephemeris data beyond the standard solar GM.
//!
//! - `PropagatorKind::NBody` – numerical N-body propagation via an 8th-order
//!   Dormand-Prince (DOP853) integrator.  The target object is integrated
//!   together with its 6×6 state transition matrix (STM) under the influence of
//!   the Sun and any additional perturbing bodies selected in `NBodyConfig`.
//!   Planetary positions are looked up from the `JPLEphem` file supplied to the
//!   differential corrector at runtime.

pub mod nbody;
pub mod planet_gm;

use crate::jpl_ephem::naif::naif_ids::NaifIds;

/// Select which propagator is used during differential orbit correction.
#[derive(Debug, Clone, Default)]
pub enum PropagatorKind {
    /// Keplerian (analytic) two-body propagation.  Fast and self-contained.
    #[default]
    TwoBody,

    /// Numerical N-body propagation with user-specified perturbing bodies.
    NBody(NBodyConfig),
}

/// Configuration for the N-body propagator.
#[derive(Debug, Clone)]
pub struct NBodyConfig {
    /// Bodies whose gravitational attraction perturbs the target orbit.
    ///
    /// Defaults to `[Sun]`.  For each body the GM is taken from [`planet_gm`]
    /// and its ephemeris is obtained from the `JPLEphem` file passed to the
    /// integrator at runtime.
    pub perturbing_bodies: Vec<NaifIds>,

    /// Absolute tolerance for the DOP853 step-size control.
    ///
    /// Units: AU.  Defaults to `1e-12`.
    pub abs_tol: f64,

    /// Relative tolerance for the DOP853 step-size control.
    ///
    /// Defaults to `1e-12`.
    pub rel_tol: f64,
}

impl Default for NBodyConfig {
    fn default() -> Self {
        use crate::jpl_ephem::naif::naif_ids::solar_system_bary::SolarSystemBary;
        Self {
            perturbing_bodies: vec![NaifIds::SSB(SolarSystemBary::Sun)],
            abs_tol: 1e-12,
            rel_tol: 1e-12,
        }
    }
}
