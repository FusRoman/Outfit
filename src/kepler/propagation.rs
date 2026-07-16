//! Full Cartesian state propagation via the universal-variable formulation.

use nalgebra::Vector3;

use crate::kepler::params::SolverType;
use crate::outfit_errors::OutfitError;
use crate::GAUSS_GRAV;

use super::params::UniversalKeplerParams;

/// Output of the universal-variable two-body propagator.
pub struct UniversalPropagResult {
    /// Propagated position vector (au).
    pub r1: Vector3<f64>,
    /// Propagated velocity vector (au/day).
    pub v1: Vector3<f64>,
    /// Lagrange coefficient $f$ such that $r_1 = f \cdot r_0 + g \cdot v_0$.
    pub f_lag: f64,
    /// Lagrange coefficient $g$ such that $r_1 = f \cdot r_0 + g \cdot v_0$.
    pub g_lag: f64,
    /// Time-derivative of $f$: $\dot{f} = -\frac{\mu}{r_0 r_1} S_1$.
    pub f_dot: f64,
    /// Time-derivative of $g$: $\dot{g} = 1 - \frac{\mu}{r_1} S_2$.
    pub g_dot: f64,
    /// Universal anomaly $\psi$ (rad, generalized units).
    ///
    /// The root of the universal Kepler equation:
    ///
    /// $$f(\psi) = r_0 \cdot s_1 + \sigma_0 \cdot s_2 + \mu \cdot s_3 - \Delta t = 0$$
    pub psy: f64,
}

/// Propagate a Cartesian state vector using the universal-variable formulation.
///
/// Goal
/// ----
/// Given an initial heliocentric position $\mathbf{r}_0$ and velocity
/// $\mathbf{v}_0$ at epoch `t0`, compute the state $(\mathbf{r}_1, \mathbf{v}_1)$
/// at epoch `t1` by solving the universal Kepler equation.
///
/// This formulation handles all conic regimes in a unified framework:
///
/// - **Elliptic** ($\alpha < 0$)
/// - **Hyperbolic** ($\alpha > 0$)
///
/// The parabolic case ($\alpha = 0$) is not supported and returns an error.
///
/// Scientific background
/// ---------------------
/// The propagation uses the Lagrange coefficients $f$, $g$, $\dot{f}$,
/// $\dot{g}$, expressed in terms of the Stumpff-like auxiliary functions
/// $s_0, s_1, s_2, s_3$:
///
/// $$\mathbf{r}_1 = f \cdot \mathbf{r}_0 + g \cdot \mathbf{v}_0$$
/// $$\mathbf{v}_1 = \dot{f} \cdot \mathbf{r}_0 + \dot{g} \cdot \mathbf{v}_0$$
///
/// where:
///
/// $$f = 1 - \frac{\mu}{r_0} s_2, \quad g = r_0 s_1 + \sigma_0 s_2$$
/// $$\dot{f} = -\frac{\mu}{r_0 r_1} s_1, \quad \dot{g} = 1 - \frac{\mu}{r_1} s_2$$
///
/// The radius at $t_1$ is recovered from the universal variable relation:
///
/// $$r_1 = r_0 s_0 + \sigma_0 s_1 + \mu s_2$$
///
/// The normalized radial velocity proxy and energy parameter are:
///
/// $$\sigma_0 = \frac{\mathbf{r}_0 \cdot \mathbf{v}_0}{\sqrt{\mu}}, \quad
///   \alpha = v_0^2 - \frac{2\mu}{r_0}$$
///
/// The eccentricity is derived from the specific angular momentum
/// $h = \|\mathbf{r}_0 \times \mathbf{v}_0\|$:
///
/// $$e_0 = \sqrt{1 + \frac{\alpha}{\mu^2} h^2}$$
///
/// Arguments
/// ---------
/// * `position` – Initial heliocentric position vector (AU).
/// * `velocity` – Initial velocity vector (AU/day).
/// * `t0` – Initial epoch (days, MJD).
/// * `t1` – Target epoch (days, MJD).
/// * `convergency` – Optional solver tolerance for [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess).
///   Defaults to $100 \varepsilon$ if `None`.
/// * `psi_guess` – Optional warm-start for the universal anomaly $\psi$.
///   Pass `None` to use the heuristic from [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni).
///
/// Return
/// ------
/// * `Ok((r1, v1))` – Position (AU) and velocity (AU/day) at `t1`.
/// * `Err(OutfitError)` – If:
///   - the input position has zero norm,
///   - the orbit is parabolic ($\alpha \approx 0$),
///   - the universal Kepler solver does not converge (see [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess)),
///   - the propagated radius $r_1$ is degenerate.
///
/// Notes
/// -----
/// - The Lagrange identity $f \dot{g} - g \dot{f} = 1$ is preserved to
///   machine precision by construction.
/// - For nearly parabolic orbits ($|\alpha| \ll 1$), convergence may degrade.
///
/// See also
/// --------
/// * [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) – Universal Kepler solver used internally.
/// * [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams) – Input container for the solver.
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Lagrange-based velocity update from two positions.
pub fn propagate_universal(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    t0: f64,
    t1: f64,
    solver_type: SolverType,
) -> Result<UniversalPropagResult, OutfitError> {
    let gravitational_parameter = GAUSS_GRAV * GAUSS_GRAV;

    let initial_radius = position.norm();
    if initial_radius < f64::EPSILON {
        return Err(OutfitError::DegenerateState(format!(
            "initial position vector has zero norm ({initial_radius})"
        )));
    }

    let (radial_velocity_proxy, energy_parameter, eccentricity) =
        initial_orbital_state(position, velocity, initial_radius, gravitational_parameter);

    let time_of_flight = t1 - t0;

    let params = UniversalKeplerParams {
        r0: initial_radius,
        sig0: radial_velocity_proxy,
        mu: gravitational_parameter,
        alpha: energy_parameter,
        dt: time_of_flight,
        e0: eccentricity,
        solver_type,
    };

    let kepler_solution = params.solve()?;

    let (s0, s1, s2, _) = kepler_solution.as_raw_stumpff();

    let propagated_radius =
        initial_radius * s0 + radial_velocity_proxy * s1 + gravitational_parameter * s2;
    if propagated_radius < f64::EPSILON {
        return Err(OutfitError::DegenerateState(format!(
            "propagated radius r1 is zero or negative ({propagated_radius})"
        )));
    }

    let lagrange_f = 1.0 - (gravitational_parameter / initial_radius) * s2;
    let lagrange_g = initial_radius * s1 + radial_velocity_proxy * s2;
    let lagrange_f_dot = -(gravitational_parameter / (initial_radius * propagated_radius)) * s1;
    let lagrange_g_dot = 1.0 - (gravitational_parameter / propagated_radius) * s2;

    let propagated_position = lagrange_f * position + lagrange_g * velocity;
    let propagated_velocity = lagrange_f_dot * position + lagrange_g_dot * velocity;

    Ok(UniversalPropagResult {
        r1: propagated_position,
        v1: propagated_velocity,
        f_lag: lagrange_f,
        g_lag: lagrange_g,
        f_dot: lagrange_f_dot,
        g_dot: lagrange_g_dot,
        psy: kepler_solution.universal_anomaly,
    })
}

/// Derive the radial-velocity proxy, energy parameter, and eccentricity of
/// the initial state from the raw position/velocity vectors.
///
/// Returns `(radial_velocity_proxy, energy_parameter, eccentricity)` where:
/// * `radial_velocity_proxy = (r0 . v0) / sqrt(mu)`,
/// * `energy_parameter = alpha = |v0|^2 - 2*mu/r0`,
/// * `eccentricity = sqrt(1 + (alpha/mu) * |r0 x v0|^2 / mu)`, clamped at 0.
fn initial_orbital_state(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    initial_radius: f64,
    gravitational_parameter: f64,
) -> (f64, f64, f64) {
    let initial_speed_squared = velocity.norm_squared();
    let radial_velocity_proxy = position.dot(velocity) / gravitational_parameter.sqrt();
    let energy_parameter = initial_speed_squared - 2.0 * gravitational_parameter / initial_radius;
    let angular_momentum_squared = position.cross(velocity).norm_squared();
    let eccentricity = (1.0
        + (energy_parameter / gravitational_parameter) * angular_momentum_squared
            / gravitational_parameter)
        .sqrt()
        .max(0.0);

    (radial_velocity_proxy, energy_parameter, eccentricity)
}
