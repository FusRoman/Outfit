//! N-body propagator with state transition matrix (DOP853).
//!
//! # Overview
//!
//! This module integrates the heliocentric equations of motion for a small
//! solar-system body together with the 6×6 **state transition matrix** (STM) Φ —
//! the matrix of partial derivatives ∂(r,v)(t1)/∂(r,v)(t0) that maps small
//! perturbations of the initial state to perturbations of the final state,
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

use crate::{jpl_ephem::JPLEphem, outfit_errors::OutfitError};

use super::{planet_gm::gm_au3_day2, NBodyConfig};

// ---------------------------------------------------------------------------
// ODE right-hand side
// ---------------------------------------------------------------------------

// The augmented state is [f64; 42]:
//   y[0..3]  = position (AU)
//   y[3..6]  = velocity (AU/day)
//   y[6..42] = STM Φ stored col-major

/// Snapshot of a perturbing body at a fixed epoch.
///
/// The snapshot is taken at t0 and held constant over the integration arc.
/// This is accurate for short arcs (≲ 30 days) where planetary motion is slow.
pub(crate) struct PerturberSnapshot {
    /// Heliocentric position of the perturber at t0, in the ecliptic J2000 frame.
    ///
    /// Units: AU.
    heliocentric_position: Vector3<f64>,

    /// Gravitational parameter GM of the perturber.
    ///
    /// Units: AU³/day².
    gravitational_parameter: f64,
}

/// ODE right-hand side for the augmented state (position, velocity, STM).
///
/// Implements [`ODE<f64, [f64; 42]>`] so that it can be driven by the DOP853
/// integrator.  The dynamics are frozen at t0: perturber positions are sampled
/// once and held constant throughout the integration arc.
pub(crate) struct NBodyOde {
    /// Perturber snapshots evaluated at t0.
    ///
    /// Each entry holds the heliocentric position and GM of one perturbing body.
    /// These values are used on every call to [`ODE::diff`] to compute the total
    /// heliocentric acceleration and the gravity-gradient matrix required by the
    /// variational equations.
    pub(crate) perturbers: Vec<PerturberSnapshot>,
}

// ---------------------------------------------------------------------------
// Acceleration sub-functions
// ---------------------------------------------------------------------------

/// Computes the direct gravitational acceleration exerted by a single perturber
/// on the small body:
///
/// ```text
/// a_direct = −GM / |asteroid_to_perturber|³ · asteroid_to_perturber
/// ```
///
/// # Arguments
///
/// * `asteroid_to_perturber` – Vector from the heliocentric position of the
///   small body to the heliocentric position of the perturber, i.e.
///   `r_perturber − r_asteroid`. Units: AU.
/// * `gravitational_parameter` – Gravitational parameter GM of the perturber.
///   Units: AU³/day².
///
/// # Returns
///
/// The direct acceleration vector (AU/day²) pointing from the small body
/// toward the perturber, scaled by GM / |d|³.
fn direct_acceleration(
    asteroid_to_perturber: Vector3<f64>,
    gravitational_parameter: f64,
) -> Vector3<f64> {
    let distance = asteroid_to_perturber.norm();
    let distance_cubed = distance.powi(3);
    -gravitational_parameter / distance_cubed * asteroid_to_perturber
}

/// Computes the indirect (heliocentric frame correction) acceleration from a
/// single perturber:
///
/// ```text
/// a_indirect = +GM / |perturber_heliocentric_pos|³ · perturber_heliocentric_pos
/// ```
///
/// Returns zero if the perturber is at (or very near) the origin (i.e. it is
/// the Sun itself, which cancels with its direct term).
///
/// # Arguments
///
/// * `perturber_heliocentric_position` – Heliocentric position of the perturber
///   in the ecliptic J2000 frame. Units: AU.
/// * `gravitational_parameter` – Gravitational parameter GM of the perturber.
///   Units: AU³/day².
///
/// # Returns
///
/// The indirect acceleration correction vector (AU/day²) that removes the
/// apparent acceleration of the heliocentric origin due to the perturber.
/// Returns [`Vector3::zeros`] when the perturber distance is ≤ 1 × 10⁻¹⁰ AU
/// (i.e. the perturber coincides with the Sun).
fn indirect_acceleration(
    perturber_heliocentric_position: Vector3<f64>,
    gravitational_parameter: f64,
) -> Vector3<f64> {
    let perturber_distance = perturber_heliocentric_position.norm();
    if perturber_distance > 1e-10 {
        let perturber_distance_cubed = perturber_distance.powi(3);
        gravitational_parameter / perturber_distance_cubed * perturber_heliocentric_position
    } else {
        Vector3::zeros()
    }
}

/// Computes the gravity-gradient matrix contribution (∂a/∂r) from a single
/// perturber:
///
/// ```text
/// ∂a/∂r += −GM · (I/|d|³ − 3 d dᵀ/|d|⁵)
/// ```
///
/// # Arguments
///
/// * `asteroid_to_perturber` – Vector from the heliocentric position of the
///   small body to that of the perturber (`r_perturber − r_asteroid`).
///   Units: AU.
/// * `gravitational_parameter` – Gravitational parameter GM of the perturber.
///   Units: AU³/day².
///
/// # Returns
///
/// A 3×3 matrix (units: day⁻²) representing this perturber's contribution to
/// ∂a/∂r, the partial derivative of the acceleration with respect to the
/// small-body position. Summing this matrix over all perturbers yields the
/// lower-left block of the linearised dynamics matrix A used in the variational
/// equations `dΦ/dt = A · Φ`.
fn gravity_gradient_contribution(
    asteroid_to_perturber: Vector3<f64>,
    gravitational_parameter: f64,
) -> Matrix3<f64> {
    let distance = asteroid_to_perturber.norm();
    let distance_cubed = distance.powi(3);
    let distance_fifth = distance.powi(5);
    let outer_product = asteroid_to_perturber * asteroid_to_perturber.transpose();
    -gravitational_parameter
        * (Matrix3::<f64>::identity() * (1.0 / distance_cubed)
            - outer_product * (3.0 / distance_fifth))
}

/// Accumulates total heliocentric acceleration and gravity-gradient matrix
/// over all perturbers.
///
/// Iterates over every [`PerturberSnapshot`] and sums the direct acceleration,
/// indirect acceleration correction, and gravity-gradient contribution from
/// each body.
///
/// # Arguments
///
/// * `asteroid_heliocentric_position` – Heliocentric position of the small body
///   in the ecliptic J2000 frame. Units: AU.
/// * `perturbers` – Slice of perturber snapshots evaluated at t0. Each entry
///   provides the perturber's heliocentric position (AU) and GM (AU³/day²).
///
/// # Returns
///
/// A tuple `(total_acceleration, da_dr)` where:
/// - `total_acceleration` – Combined heliocentric acceleration vector (AU/day²)
///   from all perturbers, including indirect terms.
/// - `da_dr` – 3×3 gravity-gradient matrix (day⁻²), i.e. the sum of each
///   perturber's ∂a/∂r contribution, used to build the linearised dynamics
///   matrix for the STM variational equations.
fn accumulate_perturber_effects(
    asteroid_heliocentric_position: Vector3<f64>,
    perturbers: &[PerturberSnapshot],
) -> (Vector3<f64>, Matrix3<f64>) {
    perturbers.iter().fold(
        (Vector3::zeros(), Matrix3::zeros()),
        |(acc_total, da_dr_total), perturber| {
            let asteroid_to_perturber =
                asteroid_heliocentric_position - perturber.heliocentric_position;

            let acc_direct =
                direct_acceleration(asteroid_to_perturber, perturber.gravitational_parameter);
            let acc_indirect = indirect_acceleration(
                perturber.heliocentric_position,
                perturber.gravitational_parameter,
            );
            let da_dr_contribution = gravity_gradient_contribution(
                asteroid_to_perturber,
                perturber.gravitational_parameter,
            );

            (
                acc_total + acc_direct + acc_indirect,
                da_dr_total + da_dr_contribution,
            )
        },
    )
}

/// Builds the 6×6 linearised dynamics matrix A from the gravity-gradient block:
///
/// ```text
/// A = | 0    I  |
///     | ∂a/∂r 0 |
/// ```
///
/// # Arguments
///
/// * `gravity_gradient` – 3×3 gravity-gradient matrix ∂a/∂r (day⁻²), as
///   returned by [`accumulate_perturber_effects`]. It occupies the lower-left
///   3×3 block of A.
///
/// # Returns
///
/// The 6×6 linearised dynamics matrix A (units: mixed — upper-right block is
/// dimensionless identity, lower-left block has units day⁻²).  This matrix
/// satisfies the variational equation `dΦ/dt = A · Φ`.
fn build_dynamics_matrix(gravity_gradient: Matrix3<f64>) -> Matrix6<f64> {
    let mut dynamics_matrix = Matrix6::<f64>::zeros();
    // Upper-right 3×3 block: identity (velocity → position coupling)
    for i in 0..3 {
        dynamics_matrix[(i, i + 3)] = 1.0;
    }
    // Lower-left 3×3 block: ∂a/∂r (acceleration → velocity coupling)
    for row in 0..3 {
        for col in 0..3 {
            dynamics_matrix[(row + 3, col)] = gravity_gradient[(row, col)];
        }
    }
    dynamics_matrix
}

/// Writes position and velocity derivatives into the first 6 components of
/// `state_derivative` from the current velocity and computed acceleration.
///
/// # Arguments
///
/// * `current_velocity` – Full 42-element augmented state vector. Indices 3–5
///   are used as the current velocity components (AU/day), which become the
///   time derivative of position.
/// * `acceleration` – Total heliocentric acceleration vector (AU/day²)
///   previously computed by [`accumulate_perturber_effects`].
/// * `state_derivative` – Mutable reference to the 42-element output derivative
///   array. On return, indices 0–2 contain `ṙ = v` (AU/day) and indices 3–5
///   contain `v̇ = a` (AU/day²).
fn write_position_velocity_derivatives(
    current_velocity: &[f64; 42],
    acceleration: Vector3<f64>,
    state_derivative: &mut [f64; 42],
) {
    state_derivative[0] = current_velocity[3];
    state_derivative[1] = current_velocity[4];
    state_derivative[2] = current_velocity[5];
    state_derivative[3] = acceleration[0];
    state_derivative[4] = acceleration[1];
    state_derivative[5] = acceleration[2];
}

/// Computes dΦ/dt = A · Φ and writes the result into `state_derivative[6..42]`.
///
/// # Arguments
///
/// * `dynamics_matrix` – 6×6 linearised dynamics matrix A built by
///   [`build_dynamics_matrix`]. Units: mixed (see that function's documentation).
/// * `augmented_state` – Full 42-element augmented state vector. Indices 6–41
///   hold the current STM Φ stored in column-major order.
/// * `state_derivative` – Mutable reference to the 42-element output derivative
///   array. On return, indices 6–41 contain the flattened (column-major) entries
///   of dΦ/dt = A · Φ.
fn write_stm_derivative(
    dynamics_matrix: Matrix6<f64>,
    augmented_state: &[f64; 42],
    state_derivative: &mut [f64; 42],
) {
    let phi = Matrix6::<f64>::from_column_slice(&augmented_state[6..42]);
    let dphi_dt = dynamics_matrix * phi;
    state_derivative[6..42].copy_from_slice(dphi_dt.as_slice());
}

impl ODE<f64, [f64; 42]> for NBodyOde {
    fn diff(&self, _time: f64, augmented_state: &[f64; 42], state_derivative: &mut [f64; 42]) {
        let asteroid_heliocentric_position =
            Vector3::new(augmented_state[0], augmented_state[1], augmented_state[2]);

        let (total_acceleration, gravity_gradient) =
            accumulate_perturber_effects(asteroid_heliocentric_position, &self.perturbers);

        write_position_velocity_derivatives(augmented_state, total_acceleration, state_derivative);

        let dynamics_matrix = build_dynamics_matrix(gravity_gradient);
        write_stm_derivative(dynamics_matrix, augmented_state, state_derivative);
    }
}

// ---------------------------------------------------------------------------
// Public interface
// ---------------------------------------------------------------------------

/// Result of a single N-body propagation step.
pub struct NBodyResult {
    /// Heliocentric position of the small body at t1, in the ecliptic J2000 frame.
    ///
    /// Units: AU.
    pub position: Vector3<f64>,

    /// Heliocentric velocity of the small body at t1, in the ecliptic J2000 frame.
    ///
    /// Units: AU/day.
    pub velocity: Vector3<f64>,

    /// Partial derivatives of the propagated position with respect to the six
    /// equinoctial orbital elements, evaluated at t1.
    ///
    /// Shape: 6×3 (rows = equinoctial elements `[a, h, k, p, q, λ]`,
    /// cols = ecliptic Cartesian components `[x, y, z]`).
    /// Units: AU / (element unit).
    pub dpos_delem: Matrix6x3<f64>,

    /// Partial derivatives of the propagated velocity with respect to the six
    /// equinoctial orbital elements, evaluated at t1.
    ///
    /// Shape: 6×3 (rows = equinoctial elements `[a, h, k, p, q, λ]`,
    /// cols = ecliptic Cartesian velocity components `[ẋ, ẏ, ż]`).
    /// Units: (AU/day) / (element unit).
    pub dvel_delem: Matrix6x3<f64>,
}

// ---------------------------------------------------------------------------
// propagate_nbody sub-functions
// ---------------------------------------------------------------------------

/// Builds the augmented initial state vector `[pos, vel, vec(Φ=I₆)]`.
///
/// Packs position and velocity into indices 0–5 and initialises the STM block
/// (indices 6–41) to the identity matrix I₆ stored in column-major order, so
/// that the integration starts with Φ(t0) = I.
///
/// # Arguments
///
/// * `initial_position` – Heliocentric position of the small body at t0.
///   Units: AU.
/// * `initial_velocity` – Heliocentric velocity of the small body at t0.
///   Units: AU/day.
///
/// # Returns
///
/// A 42-element array `[f64; 42]` with layout:
/// - indices 0–2: position components (AU),
/// - indices 3–5: velocity components (AU/day),
/// - indices 6–41: column-major entries of Φ₀ = I₆ (dimensionless).
pub(crate) fn build_augmented_initial_state(
    initial_position: Vector3<f64>,
    initial_velocity: Vector3<f64>,
) -> [f64; 42] {
    let mut augmented_state = [0.0_f64; 42];
    augmented_state[0] = initial_position[0];
    augmented_state[1] = initial_position[1];
    augmented_state[2] = initial_position[2];
    augmented_state[3] = initial_velocity[0];
    augmented_state[4] = initial_velocity[1];
    augmented_state[5] = initial_velocity[2];
    // Φ₀ = I₆ stored col-major
    augmented_state[6..42].copy_from_slice(Matrix6::<f64>::identity().as_slice());
    augmented_state
}

/// Queries the ephemeris and builds a perturber snapshot vector at the given
/// epoch.
///
/// For each body listed in [`NBodyConfig::perturbing_bodies`], this function
/// looks up the gravitational parameter from the static table in
/// [`planet_gm`](super::planet_gm) and queries the heliocentric position from
/// the supplied JPL ephemeris file.
///
/// # Arguments
///
/// * `config` – N-body configuration specifying which perturbing bodies to
///   include and the integrator tolerances.
/// * `jpl` – Opened JPL ephemeris file used to query the heliocentric position
///   of each perturbing body.
/// * `epoch` – Reference epoch at which perturber positions are sampled.
///   Passed directly to [`JPLEphem::body_ephemeris`].
///
/// # Returns
///
/// A [`Vec<PerturberSnapshot>`] with one entry per body listed in
/// `config.perturbing_bodies`, ordered identically. Each entry contains the
/// body's heliocentric position (AU) and GM (AU³/day²) at the given epoch.
///
/// # Errors
///
/// - Returns [`OutfitError::EphemerisBodyNotSupported`] if a perturbing body
///   has no GM entry in the static table or cannot be resolved by the JPL
///   ephemeris file.
pub(crate) fn build_perturber_snapshots(
    config: &NBodyConfig,
    jpl: &JPLEphem,
    epoch: &hifitime::Epoch,
) -> Result<Vec<PerturberSnapshot>, OutfitError> {
    config
        .perturbing_bodies
        .iter()
        .map(|&body| {
            let gravitational_parameter = gm_au3_day2(body).ok_or_else(|| {
                OutfitError::EphemerisBodyNotSupported(format!(
                    "No GM available for perturber {body:?}"
                ))
            })?;
            let (heliocentric_position, _velocity) = jpl.body_ephemeris(body, epoch)?;
            Ok(PerturberSnapshot {
                heliocentric_position,
                gravitational_parameter,
            })
        })
        .collect()
}

/// Runs the DOP853 integrator from t=0 to t=`time_span_days` and returns the
/// final augmented state vector.
///
/// Constructs a DOP853 initial-value problem from the provided ODE, integrates
/// from time 0 to `time_span_days` using the absolute and relative tolerances
/// specified in `config`, and returns the last computed augmented state.
///
/// # Arguments
///
/// * `ode` – Reference to the [`NBodyOde`] instance holding the frozen perturber
///   snapshots. Implements the ODE right-hand side.
/// * `augmented_initial_state` – 42-element initial augmented state vector at
///   t=0, as produced by [`build_augmented_initial_state`].
/// * `time_span_days` – Integration duration. Positive for forward propagation,
///   negative for backward propagation. Units: days.
/// * `config` – N-body configuration supplying `abs_tol` (AU) and `rel_tol`
///   (dimensionless) for the adaptive step-size control of DOP853.
///
/// # Returns
///
/// The 42-element augmented state at t=`time_span_days`:
/// - indices 0–2: propagated position (AU),
/// - indices 3–5: propagated velocity (AU/day),
/// - indices 6–41: column-major entries of Φ(t1) (dimensionless).
///
/// # Errors
///
/// - Returns [`OutfitError::NBodyPropagationFailed`] if the DOP853 solver
///   returns an error or if the solution contains no steps.
pub(crate) fn integrate_augmented_state(
    ode: &NBodyOde,
    augmented_initial_state: [f64; 42],
    time_span_days: f64,
    config: &NBodyConfig,
) -> Result<[f64; 42], OutfitError> {
    let solution = IVP::ode(ode, 0.0_f64, time_span_days, augmented_initial_state)
        .method(
            ExplicitRungeKutta::dop853()
                .atol(config.abs_tol)
                .rtol(config.rel_tol),
        )
        .solve()
        .map_err(|e| OutfitError::NBodyPropagationFailed(format!("{e:?}")))?;

    solution.y.last().copied().ok_or_else(|| {
        OutfitError::NBodyPropagationFailed("Integrator returned no steps".to_string())
    })
}

/// Converts the t0 element Jacobians `(dpos_delem0, dvel_delem0)` (each 6×3)
/// into a single 6×6 matrix whose rows are state components and whose columns
/// are element indices:
///
/// ```text
/// J0 = | (dpos_delem0)ᵀ |
///      | (dvel_delem0)ᵀ |
/// ```
///
/// # Arguments
///
/// * `dpos_delem_at_t0` – 6×3 Jacobian of heliocentric position with respect to
///   the six equinoctial elements at t0 (rows = elements, cols = Cartesian
///   components). Units: AU / (element unit).
/// * `dvel_delem_at_t0` – 6×3 Jacobian of heliocentric velocity with respect to
///   the six equinoctial elements at t0 (rows = elements, cols = Cartesian
///   velocity components). Units: (AU/day) / (element unit).
///
/// # Returns
///
/// A 6×6 matrix J0 where:
/// - rows 0–2 correspond to the three position Cartesian components,
/// - rows 3–5 correspond to the three velocity Cartesian components,
/// - column `j` corresponds to equinoctial element index `j`.
///
/// This layout is compatible with the STM Φ so that the product `Φ(t1) · J0`
/// propagates the element Jacobian to t1.
pub(crate) fn build_initial_state_jacobian(
    dpos_delem_at_t0: &Matrix6x3<f64>,
    dvel_delem_at_t0: &Matrix6x3<f64>,
) -> Matrix6<f64> {
    let mut initial_state_jacobian = Matrix6::<f64>::zeros();
    for element_index in 0..6 {
        for cartesian_component in 0..3 {
            initial_state_jacobian[(cartesian_component, element_index)] =
                dpos_delem_at_t0[(element_index, cartesian_component)];
            initial_state_jacobian[(cartesian_component + 3, element_index)] =
                dvel_delem_at_t0[(element_index, cartesian_component)];
        }
    }
    initial_state_jacobian
}

/// Splits the propagated 6×6 Jacobian `Φ(t1) · J0` back into the
/// `(dpos_delem, dvel_delem)` pair at t1 (each 6×3).
///
/// Performs the inverse of the layout established by [`build_initial_state_jacobian`]:
/// extracts the top 3 rows into `dpos_delem_at_t1` and the bottom 3 rows into
/// `dvel_delem_at_t1`, transposing the index ordering back to the 6×3 convention
/// (rows = elements, cols = Cartesian components).
///
/// # Arguments
///
/// * `propagated_state_jacobian` – 6×6 matrix `Φ(t1) · J0` resulting from
///   multiplying the STM at t1 by the initial-state Jacobian J0.
///   - rows 0–2: propagated position partials,
///   - rows 3–5: propagated velocity partials,
///   - column `j`: partial with respect to equinoctial element `j`.
///
/// # Returns
///
/// A tuple `(dpos_delem_at_t1, dvel_delem_at_t1)` where each matrix has shape
/// 6×3 (rows = equinoctial elements, cols = ecliptic Cartesian components):
/// - `dpos_delem_at_t1`: ∂pos(t1)/∂elem. Units: AU / (element unit).
/// - `dvel_delem_at_t1`: ∂vel(t1)/∂elem. Units: (AU/day) / (element unit).
pub(crate) fn split_propagated_jacobian(
    propagated_state_jacobian: Matrix6<f64>,
) -> (Matrix6x3<f64>, Matrix6x3<f64>) {
    let mut dpos_delem_at_t1 = Matrix6x3::<f64>::zeros();
    let mut dvel_delem_at_t1 = Matrix6x3::<f64>::zeros();
    for element_index in 0..6 {
        for cartesian_component in 0..3 {
            dpos_delem_at_t1[(element_index, cartesian_component)] =
                propagated_state_jacobian[(cartesian_component, element_index)];
            dvel_delem_at_t1[(element_index, cartesian_component)] =
                propagated_state_jacobian[(cartesian_component + 3, element_index)];
        }
    }
    (dpos_delem_at_t1, dvel_delem_at_t1)
}
