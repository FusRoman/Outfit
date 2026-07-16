//! Preliminary (initial-guess) solvers for the universal anomaly `psi`.
//!
//! These modules do not iterate to full convergence: they only provide a
//! starting value for `psi`, meant to be refined afterwards by
//! [`solve_kepuni`](crate::kepler::solve_kepuni).

pub mod prelim_elliptic;
pub mod prelim_hyperbolic;
pub mod prelim_parabolic;
