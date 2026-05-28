//! N-body propagator with state transition matrix (DOP853).
//!
//! # Overview
//!
//! This module integrates the heliocentric equations of motion for a small
//! solar-system body together with the 6×6 state transition matrix (STM) Φ,
//! using the explicit 8th-order Dormand–Prince (DOP853) integrator provided by
//! the [`differential-equations`](differential_equations) crate.
//!
//! ## State layout
//! The augmented state vector has 42 components stored as `[f64; 42]`:
//!
//! ```text
//! y[0..3]   – heliocentric position   r  [AU]
//! y[3..6]   – heliocentric velocity   v  [AU/day]
//! y[6..42]  – STM Φ, stored column-major (nalgebra convention) as a flat 36-element slice
//! ```
//!
//! ## Acceleration model
//! The Newtonian heliocentric acceleration for perturber `i` is
//!
//! ```text
//! a += −GMᵢ/|d|³ · d  +  GMᵢ/|rᵢ|³ · rᵢ        (indirect term)
//! ```
//!
//! where `d = r − rᵢ` is the asteroid–perturber vector, and the second term
//! removes the Sun's perturbation on the integration centre (heliocentric
//! formulation).  When the Sun itself is a perturber the indirect term cancels
//! the direct term exactly, yielding the standard two-body acceleration.
//!
//! ## Variational equations
//! ```text
//! dΦ/dt = A(t) · Φ
//! ```
//! with `A = [[0, I], [∂a/∂r, 0]]` and
//!
//! ```text
//! ∂a/∂r = Σᵢ GMᵢ · (3 d dᵀ/|d|⁵ − I/|d|³)
//! ```
//!
//! ## Element Jacobian
//! The result of this module is `dpos_delem_ecl` (Matrix6x3, rows = elements,
//! cols = ecliptic Cartesian components) that is directly substitutable for the
//! matrix returned by `solve_two_body_problem`.  It is computed as
//!
//! ```text
//! J(t1) = [Φ_pp | Φ_pv] · J0     (top 3 rows of Φ(t1) · J0_full)
//! ```
//!
//! where `J0_full` is the 6×6 matrix formed by stacking `dpos_delem` and
//! `dvel_delem` (each 6×3) from a zero-propagation call to
//! `solve_two_body_problem`.

use differential_equations::ode::ODE;
use differential_equations::prelude::*;
use nalgebra::{Matrix3, Matrix6, Matrix6x3, Vector3};

use crate::{
    jpl_ephem::JPLEphem, orbit_type::equinoctial_element::EquinoctialElements,
    outfit_errors::OutfitError,
};

use super::{planet_gm::gm_au3_day2, NBodyConfig};

// ---------------------------------------------------------------------------
// ODE right-hand side
// ---------------------------------------------------------------------------

// The augmented state is [f64; 42]:
//   y[0..3]  = position (AU)
//   y[3..6]  = velocity (AU/day)
//   y[6..42] = STM Φ stored col-major

// Snapshot of a perturbing body: (heliocentric position [AU], GM [AU³/day²]).
// The snapshot is taken at t0 and held constant over the integration arc.
// This is accurate for short arcs (≲ 30 days) where planetary motion is slow.
struct PerturberSnapshot {
    pos: Vector3<f64>,
    gm: f64,
}

struct NBodyOde {
    perturbers: Vec<PerturberSnapshot>,
}

impl ODE<f64, [f64; 42]> for NBodyOde {
    fn diff(&self, _t: f64, y: &[f64; 42], dydt: &mut [f64; 42]) {
        let r = Vector3::new(y[0], y[1], y[2]);

        let mut acc = Vector3::<f64>::zeros();
        let mut da_dr = Matrix3::<f64>::zeros();

        for p in &self.perturbers {
            let d: Vector3<f64> = r - p.pos;
            let d_norm: f64 = d.norm();
            let d3: f64 = d_norm * d_norm * d_norm;
            let d5: f64 = d3 * d_norm * d_norm;

            // Direct acceleration: -GM/|d|³ · d
            acc -= p.gm / d3 * d;

            // Indirect term (heliocentric formulation): +GM/|r_body|³ · r_body
            let rb_norm: f64 = p.pos.norm();
            if rb_norm > 1e-10 {
                let rb3 = rb_norm * rb_norm * rb_norm;
                acc += p.gm / rb3 * p.pos;
            }

            // Gravity gradient: ∂a/∂r contribution
            // = −GM · (I/d³ − 3 d dᵀ/d⁵)
            let outer: Matrix3<f64> = d * d.transpose();
            da_dr -= p.gm * (Matrix3::<f64>::identity() * (1.0 / d3) - outer * (3.0 / d5));
        }

        // State derivatives
        dydt[0] = y[3];
        dydt[1] = y[4];
        dydt[2] = y[5];
        dydt[3] = acc[0];
        dydt[4] = acc[1];
        dydt[5] = acc[2];

        // Variational equations: dΦ/dt = A(t) · Φ
        // A = [[0, I], [da_dr, 0]] — 6×6
        let mut a_mat = Matrix6::<f64>::zeros();
        // Upper-right 3×3 block: identity
        for i in 0..3 {
            a_mat[(i, i + 3)] = 1.0;
        }
        // Lower-left 3×3 block: da_dr
        for row in 0..3 {
            for col in 0..3 {
                a_mat[(row + 3, col)] = da_dr[(row, col)];
            }
        }

        // Φ from y[6..42], col-major
        let phi = Matrix6::<f64>::from_column_slice(&y[6..42]);
        let dphi_dt = a_mat * phi;

        // Store dΦ/dt back into dydt[6..42]
        let dphi_slice = dphi_dt.as_slice();
        dydt[6..42].copy_from_slice(dphi_slice);
    }
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

/// Result of a single N-body propagation step.
pub struct NBodyResult {
    /// Heliocentric position at t1 (AU), ecliptic J2000.
    pub position: Vector3<f64>,
    /// Heliocentric velocity at t1 [AU/day], ecliptic J2000.
    pub velocity: Vector3<f64>,
    /// ∂pos(t1)/∂elem  (6×3, rows = elements, cols = ecliptic Cartesian).
    pub dpos_delem: Matrix6x3<f64>,
    /// ∂vel(t1)/∂elem  (6×3, rows = elements, cols = ecliptic Cartesian).
    pub dvel_delem: Matrix6x3<f64>,
}

/// Propagate `elements` from their reference epoch to `t1_mjd_tt` using
/// numerical N-body integration.
///
/// # Arguments
/// * `elements` – equinoctial orbital elements (reference epoch = t0).
/// * `t1_mjd_tt` – target epoch as MJD-TT (days).
/// * `jpl` – ephemeris file used to query perturber positions.
/// * `config` – N-body configuration (perturbers, tolerances).
///
/// # Errors
/// Returns [`OutfitError::NBodyPropagationFailed`] if the integrator fails,
/// or [`OutfitError::EphemerisBodyNotSupported`] if a perturber cannot be
/// looked up in the supplied ephemeris.
pub fn propagate_nbody(
    elements: &EquinoctialElements,
    t1_mjd_tt: f64,
    jpl: &JPLEphem,
    config: &NBodyConfig,
) -> Result<NBodyResult, OutfitError> {
    let t0_mjd_tt = elements.reference_epoch;
    let dt = t1_mjd_tt - t0_mjd_tt;

    // Initial cartesian state and element Jacobian J0 = ∂(x0,v0)/∂elem
    let (pos0, vel0, jacs0) = elements.solve_two_body_problem(0.0, 0.0, true)?;
    let (dpos_delem0, dvel_delem0) =
        jacs0.expect("solve_two_body_problem(0,0,true) must return jacobians");

    // If dt ≈ 0 skip integration
    if dt.abs() < 1e-14 {
        return Ok(NBodyResult {
            position: pos0,
            velocity: vel0,
            dpos_delem: dpos_delem0,
            dvel_delem: dvel_delem0,
        });
    }

    // Build augmented initial state: [pos, vel, vec(Φ=I)]
    let mut y0 = [0.0_f64; 42];
    y0[0] = pos0[0];
    y0[1] = pos0[1];
    y0[2] = pos0[2];
    y0[3] = vel0[0];
    y0[4] = vel0[1];
    y0[5] = vel0[2];
    // Φ₀ = I₆, col-major: diagonal elements at indices 6,13,20,27,34,41
    let phi0 = Matrix6::<f64>::identity();
    y0[6..42].copy_from_slice(phi0.as_slice());

    // Build perturber snapshot at t0
    let epoch_t0 = hifitime::Epoch::from_mjd_in_time_scale(t0_mjd_tt, hifitime::TimeScale::TT);
    let mut perturbers = Vec::with_capacity(config.perturbing_bodies.len());
    for &body in &config.perturbing_bodies {
        let gm = gm_au3_day2(body).ok_or_else(|| {
            OutfitError::EphemerisBodyNotSupported(format!(
                "No GM available for perturber {body:?}"
            ))
        })?;
        let (pos_body, _) = jpl.body_ephemeris(body, &epoch_t0)?;
        perturbers.push(PerturberSnapshot { pos: pos_body, gm });
    }

    let ode = NBodyOde { perturbers };

    let sol = IVP::ode(&ode, 0.0_f64, dt, y0)
        .method(
            ExplicitRungeKutta::dop853()
                .atol(config.abs_tol)
                .rtol(config.rel_tol),
        )
        .solve()
        .map_err(|e| OutfitError::NBodyPropagationFailed(format!("{e:?}")))?;

    let y1 = *sol.y.last().ok_or_else(|| {
        OutfitError::NBodyPropagationFailed("Integrator returned no steps".to_string())
    })?;

    let pos1 = Vector3::new(y1[0], y1[1], y1[2]);
    let vel1 = Vector3::new(y1[3], y1[4], y1[5]);
    let phi1 = Matrix6::<f64>::from_column_slice(&y1[6..42]);

    // Element Jacobian via chain rule
    //
    // J0_state_by_elem: 6×6 where
    //   rows 0..3 = dpos/delem (transposed from dpos_delem0 which is 6elem×3pos)
    //   rows 3..6 = dvel/delem (transposed from dvel_delem0)
    //   cols = element index
    let mut j0 = Matrix6::<f64>::zeros();
    for elem_j in 0..6 {
        for pos_i in 0..3 {
            j0[(pos_i, elem_j)] = dpos_delem0[(elem_j, pos_i)];
            j0[(pos_i + 3, elem_j)] = dvel_delem0[(elem_j, pos_i)];
        }
    }

    let jac1 = phi1 * j0; // 6×6: rows=state_t1 component, cols=element

    let mut dpos_delem1 = Matrix6x3::<f64>::zeros();
    let mut dvel_delem1 = Matrix6x3::<f64>::zeros();
    for elem_j in 0..6 {
        for pos_i in 0..3 {
            dpos_delem1[(elem_j, pos_i)] = jac1[(pos_i, elem_j)];
            dvel_delem1[(elem_j, pos_i)] = jac1[(pos_i + 3, elem_j)];
        }
    }

    Ok(NBodyResult {
        position: pos1,
        velocity: vel1,
        dpos_delem: dpos_delem1,
        dvel_delem: dvel_delem1,
    })
}
