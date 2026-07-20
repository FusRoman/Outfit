//! Full Cartesian state propagation via the universal-variable formulation.

use nalgebra::Vector3;

use crate::kepler::params::SolverType;
use crate::outfit_errors::OutfitError;
use crate::GAUSS_GRAV;

use super::params::UniversalKeplerParams;

/// Output of the universal-variable two-body propagator.
#[derive(Debug)]
pub struct UniversalPropagResult {
    /// Propagated position vector (au).
    pub r1: Vector3<f64>,
    /// Propagated velocity vector (au/day).
    pub v1: Vector3<f64>,
    /// Lagrange coefficient $f$ such that $r_1 = f \cdot r_0 + g \cdot v_0$.
    pub f_lag: f64,
    /// Lagrange coefficient $g$ such that $r_1 = f \cdot r_0 + g \cdot v_0$.
    pub g_lag: f64,
    /// Time-derivative of $f$: $\dot{f} = -\frac{\sqrt{\mu}}{r_0 r_1} S_1$.
    pub f_dot: f64,
    /// Time-derivative of $g$: $\dot{g} = 1 - \frac{1}{r_1} S_2$.
    pub g_dot: f64,
    /// Universal anomaly $\psi$ (rad, generalized units).
    ///
    /// The root of the universal Kepler equation:
    ///
    /// $$f(\psi) = r_0 \cdot s_1 + \sigma_0 \cdot s_2 + s_3 - \sqrt{\mu} \cdot \Delta t = 0$$
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
/// $$f = 1 - \frac{s_2}{r_0}, \quad g = \frac{r_0 s_1 + \sigma_0 s_2}{\sqrt{\mu}}$$
/// $$\dot{f} = -\frac{\sqrt{\mu}}{r_0 r_1} s_1, \quad \dot{g} = 1 - \frac{s_2}{r_1}$$
///
/// The radius at $t_1$ is recovered from the universal variable relation:
///
/// $$r_1 = r_0 s_0 + \sigma_0 s_1 + s_2$$
///
/// The normalized radial velocity proxy and energy parameter are:
///
/// $$\sigma_0 = \frac{\mathbf{r}_0 \cdot \mathbf{v}_0}{\sqrt{\mu}}, \quad
///   \alpha = \frac{v_0^2 - \frac{2\mu}{r_0}}{\mu} = -\frac{1}{a}$$
///
/// `alpha` is the reciprocal semi-major axis convention (dimension
/// `1/length`) expected by [`s_funct`](crate::kepler::s_funct) and the rest
/// of the universal-variable machinery — *not* the raw vis-viva $2E$ (which
/// has dimension velocity$^2$). This division by $\mu$ is what keeps `alpha`
/// dimensionally consistent with the (mu-free) formulas above.
///
/// The eccentricity is derived from the specific angular momentum
/// $h = \|\mathbf{r}_0 \times \mathbf{v}_0\|$:
///
/// $$e_0 = \sqrt{1 + \frac{\alpha}{\mu} h^2}$$
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
    let sqrt_mu = gravitational_parameter.sqrt();

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

    let propagated_radius = initial_radius * s0 + radial_velocity_proxy * s1 + s2;
    if propagated_radius < f64::EPSILON {
        return Err(OutfitError::DegenerateState(format!(
            "propagated radius r1 is zero or negative ({propagated_radius})"
        )));
    }

    let lagrange_f = 1.0 - s2 / initial_radius;
    let lagrange_g = (initial_radius * s1 + radial_velocity_proxy * s2) / sqrt_mu;
    let lagrange_f_dot = -(sqrt_mu / (initial_radius * propagated_radius)) * s1;
    let lagrange_g_dot = 1.0 - s2 / propagated_radius;

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
/// * `energy_parameter = alpha = (|v0|^2 - 2*mu/r0) / mu = -1/a`,
/// * `eccentricity = sqrt(1 + alpha * |r0 x v0|^2 / mu)`, clamped at 0.
///
/// `alpha` is expressed in the `1/length` convention (reciprocal semi-major
/// axis) expected by [`s_funct`](crate::kepler::s_funct) and the rest of the
/// universal-variable machinery, not the raw vis-viva `2E = |v0|^2 - 2mu/r0`
/// (dimension velocity^2): dividing by `mu` here is what keeps `alpha`,
/// `s_funct`'s output, and the Lagrange coefficients below dimensionally
/// consistent.
fn initial_orbital_state(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    initial_radius: f64,
    gravitational_parameter: f64,
) -> (f64, f64, f64) {
    let initial_speed_squared = velocity.norm_squared();
    let radial_velocity_proxy = position.dot(velocity) / gravitational_parameter.sqrt();
    let energy_parameter = (initial_speed_squared - 2.0 * gravitational_parameter / initial_radius)
        / gravitational_parameter;
    let angular_momentum_squared = position.cross(velocity).norm_squared();
    let eccentricity = (1.0
        + energy_parameter * angular_momentum_squared / gravitational_parameter)
        .sqrt()
        .max(0.0);

    (radial_velocity_proxy, energy_parameter, eccentricity)
}

#[cfg(test)]
mod tests_propagation {

    use crate::kepler::{
        prelim_kepler::prelim_parabolic::ParabolicPrelimMethod, SolverKind, SolverParams,
    };

    use super::*;

    #[test]
    fn test_propag() {
        let position: Vector3<f64> = Vector3::new(
            -8.264_959_160_036_185e-1,
            3.919_660_608_486_096_3e-1,
            2.229_919_607_182_842_5e-2,
        );
        let velocity: Vector3<f64> = Vector3::new(
            -5.447_367_111_934_2e-3,
            -2.107_596_146_728_544e-2,
            1.560_811_152_125_889_6e-3,
        );
        let t0: f64 = 6.072_555_422_778_894e4;
        let t1: f64 = 6.072_754_670_468_881_5e4;
        let solver_type = SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: None,
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        };

        let result = propagate_universal(&position, &velocity, t0, t1, solver_type)
            .expect("propagation should converge on this real fink-fat Kalman-filter state");

        // Ground truth obtained by independently integrating the two-body
        // ODE (scipy `solve_ivp`, DOP853, rtol=1e-13) from the same initial
        // state, cross-checked against a classical (a, e, E) Kepler-equation
        // solve to 50 significant digits (mpmath). This is what exposed the
        // `alpha`/`mu` scaling bug this test was written to catch.
        let expected_r1 = Vector3::new(-0.83670766718652, 0.34968043043198, 0.02539102537652);
        let expected_v1 = Vector3::new(-0.00479883489139, -0.02136507308119, 0.00154221064858);

        assert!(
            (result.r1 - expected_r1).norm() < 1e-9,
            "propagated position mismatch: got {:?}, expected {expected_r1:?}",
            result.r1
        );
        assert!(
            (result.v1 - expected_v1).norm() < 1e-9,
            "propagated velocity mismatch: got {:?}, expected {expected_v1:?}",
            result.v1
        );
    }

    #[test]
    fn test_propag2() {
        let position = Vector3::new(
            -8.209_687_552_250_132e-1,
            3.782_813_412_927_746e-1,
            2.567_330_540_285_757_8e-2,
        );
        let velocity = Vector3::new(
            -5.290_803_826_727_631e-3,
            -2.120_754_244_524_938_2e-2,
            1.601_930_231_829_977e-3,
        );
        let t0 = 6.072_555_414_035_025e4;
        let t1 = 6.072_754_661_725_012_6e4;
        let solver_type = SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: None,
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        };

        let result = propagate_universal(&position, &velocity, t0, t1, solver_type)
            .expect("propagation should converge on this real fink-fat Kalman-filter state");

        // Ground truth obtained the same way as in `test_propag`: independent
        // two-body ODE integration (scipy `solve_ivp`, DOP853, rtol=1e-13),
        // cross-checked against a classical (a, e, E) Kepler-equation solve
        // to 50 significant digits (mpmath).
        let expected_r1 = Vector3::new(
            -0.8308499934162212,
            0.33573406780460846,
            0.028843689480680244,
        );
        let expected_v1 = Vector3::new(
            -0.004623556668660562,
            -0.021495885832796668,
            0.0015799033389438464,
        );

        assert!(
            (result.r1 - expected_r1).norm() < 1e-9,
            "propagated position mismatch: got {:?}, expected {expected_r1:?}",
            result.r1
        );
        assert!(
            (result.v1 - expected_v1).norm() < 1e-9,
            "propagated velocity mismatch: got {:?}, expected {expected_v1:?}",
            result.v1
        );
    }

    #[test]
    fn test_propag3() {
        let position = Vector3::new(
            -8.146_048_077_331_896e-1,
            3.625_248_181_551_134_5e-1,
            2.955_823_936_342_896e-2,
        );
        let velocity = Vector3::new(
            -5.110_839_457_442_879e-3,
            -2.135_829_675_942_633_3e-2,
            1.649_090_267_256_617_4e-3,
        );
        let t0 = 6.072_555_403_967_375e4;
        let t1 = 6.072_754_651_657_362_4e4;
        let solver_type = SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: None,
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        };

        let result = propagate_universal(&position, &velocity, t0, t1, solver_type)
            .expect("propagation should converge on this real fink-fat Kalman-filter state");

        // Ground truth obtained the same way as in `test_propag`: independent
        // two-body ODE integration (scipy `solve_ivp`, DOP853, rtol=1e-13),
        // cross-checked against a classical (a, e, E) Kepler-equation solve
        // to 50 significant digits (mpmath).
        let expected_r1 = Vector3::new(
            -0.8241054960270079,
            0.31967830644033735,
            0.03281843272600818,
        );
        let expected_v1 = Vector3::new(
            -0.004421449930078581,
            -0.02164520905453043,
            0.0016228438077301268,
        );

        assert!(
            (result.r1 - expected_r1).norm() < 1e-9,
            "propagated position mismatch: got {:?}, expected {expected_r1:?}",
            result.r1
        );
        assert!(
            (result.v1 - expected_v1).norm() < 1e-9,
            "propagated velocity mismatch: got {:?}, expected {expected_v1:?}",
            result.v1
        );
    }

    #[test]
    fn test_propag4() {
        let position: Vector3<f64> = Vector3::new(
            -7.955_255_488_725_416e-1,
            2.548_930_678_464_935_7e-1,
            4.862_868_244_426_771e-2,
        );
        let velocity: Vector3<f64> = Vector3::new(
            4.922_319_532_929_707e-3,
            -2.500_697_952_944_025_4e-2,
            1.270_734_650_756_455_2e-3,
        );
        let t0: f64 = 6.072_754_601_265_45e4;
        let t1: f64 = 6.072_754_948_485_463e4;
        let solver_type = SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: Some(4.117043559944161),
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        };

        let result = propagate_universal(&position, &velocity, t0, t1, solver_type)
            .expect("propagation should converge on this real fink-fat Kalman-filter state");

        // Ground truth obtained the same way as in `test_propag`: independent
        // two-body ODE integration (scipy `solve_ivp`, DOP853, rtol=1e-13),
        // cross-checked against a classical (a, e, E) Kepler-equation solve
        // to 50 significant digits (mpmath).
        let expected_r1 =
            Vector3::new(-0.795508455171958, 0.25480623783297346, 0.04863309454122529);
        let expected_v1 = Vector3::new(
            0.004923714681911317,
            -0.025007426475553637,
            0.0012706493636607933,
        );

        assert!(
            (result.r1 - expected_r1).norm() < 1e-9,
            "propagated position mismatch: got {:?}, expected {expected_r1:?}",
            result.r1
        );
        assert!(
            (result.v1 - expected_v1).norm() < 1e-9,
            "propagated velocity mismatch: got {:?}, expected {expected_v1:?}",
            result.v1
        );
    }
}

/// Deterministic edge-case tests for [`propagate_universal`].
///
/// Each orbital regime below has its own ground truth, obtained the same way
/// as the real-data cases in [`tests_propagation`]: independent two-body ODE
/// integration (scipy `solve_ivp`, DOP853, `rtol=1e-13`), verified against a
/// classical orbital-elements construction (the state vectors themselves are
/// built from `(q, e, i, Ω, ω, anomaly)` elements, so the elements and the
/// ODE integration are two independent checks on the same number).
#[cfg(test)]
mod tests_propagation_edge_cases {
    use super::*;

    use crate::kepler::{
        prelim_kepler::prelim_parabolic::ParabolicPrelimMethod, SolverKind, SolverParams,
    };

    /// Reference MJD epoch used by every test below (arbitrary; only `dt`
    /// matters for the two-body dynamics).
    const T0: f64 = 60000.0;

    fn default_solver_type() -> SolverType {
        SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: None,
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        }
    }

    fn assert_vec_close(actual: Vector3<f64>, expected: Vector3<f64>, tol: f64, label: &str) {
        let diff = (actual - expected).norm();
        assert!(
            diff < tol,
            "{label} mismatch: got {actual:?}, expected {expected:?}, diff={diff:e} (tol={tol:e})"
        );
    }

    #[test]
    fn test_quasi_circular_orbit() {
        // a=2.3 AU, e=1e-4: the near-circular MBA orbits a Kalman filter
        // converges to most often.
        let position = Vector3::new(0.21053202031103074, 2.2808559935794808, 0.20574336652469302);
        let velocity = Vector3::new(
            -0.011240896923812553,
            0.0009286862636099264,
            0.0012093812594121695,
        );

        let result = propagate_universal(&position, &velocity, T0, T0 + 1.0, default_solver_type())
            .expect("quasi-circular orbit should converge");

        let expected_r1 =
            Vector3::new(0.19928860805173557, 2.2817569317219633, 0.20695024021604705);
        let expected_v1 = Vector3::new(
            -0.011245882006125182,
            0.0008731863707498086,
            0.0012043612300275948,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_high_eccentricity_near_perihelion() {
        // a=2.5 AU, e=0.95, propagated 2 days from just past perihelion,
        // where s2/s3 vary fastest.
        let position = Vector3::new(
            -0.07334498751332834,
            -0.33895907499571465,
            -0.028392141622553987,
        );
        let velocity = Vector3::new(
            0.016602035346734694,
            -0.03512248571568977,
            -0.00855803655134319,
        );

        let result = propagate_universal(&position, &velocity, T0, T0 + 2.0, default_solver_type())
            .expect("high-eccentricity orbit should converge");

        let expected_r1 = Vector3::new(
            -0.03938955874240485,
            -0.4049151331571937,
            -0.045107962550348484,
        );
        let expected_v1 = Vector3::new(
            0.017239827060966013,
            -0.03104380931069492,
            -0.008159597609238184,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_near_parabolic_elliptic() {
        // alpha = -1e-6 (a = 1e6 AU), q = 1.5 AU: stresses the
        // `1/sqrt(|alpha|)` term in `prelim_elliptic` without alpha being
        // exactly zero.
        let position = Vector3::new(-39.52202041820489, -31.796429108199263, -8.899233495929744);
        let velocity = Vector3::new(
            -0.0021940078334573457,
            -0.0024735491628538825,
            -0.0007479533467381286,
        );

        let result = propagate_universal(&position, &velocity, T0, T0 + 5.0, default_solver_type())
            .expect("near-parabolic (elliptic side) orbit should converge");

        let expected_r1 = Vector3::new(-39.532989387326914, -31.808795993104074, -8.90297302170826);
        let expected_v1 = Vector3::new(
            -0.0021935798649369708,
            -0.002473204832571767,
            -0.0007478569735428885,
        );

        assert_vec_close(result.r1, expected_r1, 1e-8, "position");
        assert_vec_close(result.v1, expected_v1, 1e-8, "velocity");
    }

    #[test]
    fn test_near_parabolic_hyperbolic() {
        // alpha = +1e-6 (a = -1e6 AU), q = 1.5 AU: same stress as above, on
        // the hyperbolic side (`prelim_hyperbolic`).
        let position = Vector3::new(-39.52295104602465, -31.79684024886987, -8.899322047228939);
        let velocity = Vector3::new(
            0.002194072590815172,
            0.0024735663723493223,
            0.0007479554224776113,
        );

        let result = propagate_universal(&position, &velocity, T0, T0 + 5.0, default_solver_type())
            .expect("near-parabolic (hyperbolic side) orbit should converge");

        let expected_r1 = Vector3::new(-39.51197961256531, -31.784471555801662, -8.895582029083993);
        let expected_v1 = Vector3::new(
            0.0021945008425622334,
            0.0024739108884688134,
            0.0007480518443838213,
        );

        assert_vec_close(result.r1, expected_r1, 1e-8, "position");
        assert_vec_close(result.v1, expected_v1, 1e-8, "velocity");
    }

    #[test]
    fn test_hyperbolic_orbit() {
        // a=-1.5 AU, e=2.0: clearly unbound orbit (can occur if a Kalman
        // state estimate momentarily overshoots into an unbound regime).
        let position = Vector3::new(1.249244222099196, 0.79545767943562, 0.31874549259903184);
        let velocity = Vector3::new(
            0.011307076273210306,
            -0.019273992715882825,
            -0.00941294513552841,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 + 10.0, default_solver_type())
                .expect("hyperbolic orbit should converge");

        let expected_r1 = Vector3::new(1.356763056837531, 0.5995655052692187, 0.22337735587386795);
        let expected_v1 = Vector3::new(
            0.010175595301664717,
            -0.019879144421132384,
            -0.009648073244003378,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_dt_near_zero_returns_initial_state() {
        // dt = 1e-8 day (~0.9 ms): intra-night cadence limit. The propagated
        // state must be indistinguishable from the initial one.
        let position = Vector3::new(
            -8.264_959_160_036_185e-1,
            3.919_660_608_486_096_3e-1,
            2.229_919_607_182_842_5e-2,
        );
        let velocity = Vector3::new(
            -5.447_367_111_934_2e-3,
            -2.107_596_146_728_544e-2,
            1.560_811_152_125_889_6e-3,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 + 1e-8, default_solver_type())
                .expect("near-zero dt propagation should converge");

        assert_vec_close(result.r1, position, 1e-8, "position");
        assert_vec_close(result.v1, velocity, 1e-8, "velocity");
    }

    #[test]
    fn test_negative_dt_backward_propagation() {
        // a=1.8 AU, e=0.3, dt=-10 days: the Kalman filter also performs
        // backward (smoothing) corrections.
        let position = Vector3::new(
            -1.9971284992478424,
            -0.3712522565219666,
            0.07103537618726069,
        );
        let velocity = Vector3::new(
            -0.0011593869551458037,
            -0.010998451678232324,
            -0.0021125210607977427,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 - 10.0, default_solver_type())
                .expect("backward propagation should converge");

        let expected_r1 = Vector3::new(-1.981969219563329, -0.260669755669865, 0.09202082642052144);
        let expected_v1 = Vector3::new(
            -0.0018771502422668131,
            -0.011112287894550531,
            -0.0020830780553628306,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_gap_35_days_ztf_lsst_cadence() {
        // a=2.2 AU, e=0.15, dt=+35 days: the concrete "gap > 30 days" case
        // reported from real ZTF cadence, and the exact regime where the
        // `alpha`/`mu` scaling bug (fixed by `test_propag`) manifested.
        let position = Vector3::new(-0.5081279057987266, 1.9247110895634725, 0.35100377520249654);
        let velocity = Vector3::new(
            -0.01246619702904591,
            -0.0016710189605835082,
            0.00028431118257845765,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 + 35.0, default_solver_type())
                .expect("35-day gap propagation should converge");

        let expected_r1 = Vector3::new(-0.9305835567844181, 1.8256992698318637, 0.3534197183411793);
        let expected_v1 = Vector3::new(
            -0.011601678565060293,
            -0.0039348889411689745,
            -0.00014072224048158237,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_gap_45_days_negative() {
        // Same object as the 35-day case, dt=-45 days: a larger gap, backward.
        let position = Vector3::new(-0.5081279057987266, 1.9247110895634725, 0.35100377520249654);
        let velocity = Vector3::new(
            -0.01246619702904591,
            -0.0016710189605835082,
            0.00028431118257845765,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 - 45.0, default_solver_type())
                .expect("45-day backward gap propagation should converge");

        let expected_r1 = Vector3::new(0.06469927718927299, 1.9271232755936696, 0.3252726272245291);
        let expected_v1 = Vector3::new(
            -0.012836838087040878,
            0.0016190222436126842,
            0.0008615818662658044,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_gap_150_days_short_period_neo() {
        // a=1.0 AU, e=0.4 (period ~365 days), dt=150 days: recovery after a
        // long seasonal gap, still under one full revolution.
        let position = Vector3::new(1.2047460759903845, 0.6820351942449735, -0.09343868944527914);
        let velocity = Vector3::new(
            -0.006617740910480954,
            0.009290304236214792,
            0.0007113219381064047,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 + 150.0, default_solver_type())
                .expect("150-day gap propagation should converge");

        let expected_r1 = Vector3::new(
            -0.6120201957110859,
            0.10949611603813385,
            0.054394854746006055,
        );
        let expected_v1 = Vector3::new(
            -0.0007514307771635531,
            -0.02552809717488894,
            -0.00032308624673825666,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_gap_400_days_multi_revolution() {
        // Same object as the 150-day case, dt=400 days: exceeds the ~365-day
        // period, so `psi` must span more than one revolution — the failure
        // mode directly tied to the root cause fixed alongside `test_propag`.
        let position = Vector3::new(1.2047460759903845, 0.6820351942449735, -0.09343868944527914);
        let velocity = Vector3::new(
            -0.006617740910480954,
            0.009290304236214792,
            0.0007113219381064047,
        );

        let result =
            propagate_universal(&position, &velocity, T0, T0 + 400.0, default_solver_type())
                .expect("multi-revolution 400-day gap propagation should converge");

        let expected_r1 =
            Vector3::new(0.8970232454981946, 0.9499966243304765, -0.06285450365091588);
        let expected_v1 = Vector3::new(
            -0.011011765095008376,
            0.005846950052656893,
            0.001037596639746061,
        );

        assert_vec_close(result.r1, expected_r1, 1e-9, "position");
        assert_vec_close(result.v1, expected_v1, 1e-9, "velocity");
    }

    #[test]
    fn test_continuity_across_lunation_gap() {
        // Same object as the 35-day case, evaluated at dt=34.9 and dt=35.1
        // days: both must independently match their ground truth, which
        // rules out a discontinuous jump (e.g. an internal solver-branch or
        // revolution-count flip) between two very close gap lengths.
        let position = Vector3::new(-0.5081279057987266, 1.9247110895634725, 0.35100377520249654);
        let velocity = Vector3::new(
            -0.01246619702904591,
            -0.0016710189605835082,
            0.00028431118257845765,
        );

        let before =
            propagate_universal(&position, &velocity, T0, T0 + 34.9, default_solver_type())
                .expect("dt=34.9 propagation should converge");
        let after = propagate_universal(&position, &velocity, T0, T0 + 35.1, default_solver_type())
            .expect("dt=35.1 propagation should converge");

        assert_vec_close(
            before.r1,
            Vector3::new(-0.9294232358528046, 1.8260924582635572, 0.3534337324049776),
            1e-9,
            "position at dt=34.9",
        );
        assert_vec_close(
            after.r1,
            Vector3::new(-0.9317435714637801, 1.8253054805679922, 0.3534055879680079),
            1e-9,
            "position at dt=35.1",
        );

        // Sanity check on top of the ground-truth match: 0.2 days apart, the
        // two positions should differ by an amount consistent with the
        // local velocity (~|v|*0.2), not by a discontinuous jump.
        let expected_scale = before.v1.norm() * 0.2;
        let actual_gap = (after.r1 - before.r1).norm();
        assert!(
            actual_gap < 10.0 * expected_scale,
            "unexpectedly large jump between dt=34.9 and dt=35.1: {actual_gap} (expected order {expected_scale})"
        );
    }

    #[test]
    fn test_degenerate_position_returns_err() {
        let position = Vector3::new(0.0, 0.0, 0.0);
        let velocity = Vector3::new(0.01, 0.0, 0.0);

        let result = propagate_universal(&position, &velocity, T0, T0 + 1.0, default_solver_type());

        assert!(
            matches!(result, Err(OutfitError::DegenerateState(_))),
            "expected DegenerateState error, got {result:?}"
        );
    }

    #[test]
    fn test_warm_start_consistent_with_cold_start() {
        let position = Vector3::new(
            -8.264_959_160_036_185e-1,
            3.919_660_608_486_096_3e-1,
            2.229_919_607_182_842_5e-2,
        );
        let velocity = Vector3::new(
            -5.447_367_111_934_2e-3,
            -2.107_596_146_728_544e-2,
            1.560_811_152_125_889_6e-3,
        );

        let base = propagate_universal(&position, &velocity, T0, T0 + 10.0, default_solver_type())
            .expect("base propagation should converge");

        let mut warm_solver_type = default_solver_type();
        warm_solver_type.params.psi_guess = Some(base.psy);

        let warm = propagate_universal(&position, &velocity, T0, T0 + 10.5, warm_solver_type)
            .expect("warm-started propagation should converge");
        let cold = propagate_universal(&position, &velocity, T0, T0 + 10.5, default_solver_type())
            .expect("cold-started propagation should converge");

        assert_vec_close(warm.r1, cold.r1, 1e-9, "position");
        assert_vec_close(warm.v1, cold.v1, 1e-9, "velocity");
    }

    #[test]
    fn test_solver_kind_agreement() {
        let position = Vector3::new(
            -8.264_959_160_036_185e-1,
            3.919_660_608_486_096_3e-1,
            2.229_919_607_182_842_5e-2,
        );
        let velocity = Vector3::new(
            -5.447_367_111_934_2e-3,
            -2.107_596_146_728_544e-2,
            1.560_811_152_125_889_6e-3,
        );
        let t1 = T0 + 1.992476899875328;

        let mut newton_solver_type = default_solver_type();
        newton_solver_type.kind = SolverKind::NewtonRaphson;
        let mut brent_solver_type = default_solver_type();
        brent_solver_type.kind = SolverKind::BrentDecker;
        let mut auto_solver_type = default_solver_type();
        auto_solver_type.kind = SolverKind::Auto;

        let newton = propagate_universal(&position, &velocity, T0, t1, newton_solver_type)
            .expect("Newton-Raphson should converge");
        let brent = propagate_universal(&position, &velocity, T0, t1, brent_solver_type)
            .expect("Brent-Dekker should converge");
        let auto = propagate_universal(&position, &velocity, T0, t1, auto_solver_type)
            .expect("Auto should converge");

        assert_vec_close(newton.r1, brent.r1, 1e-9, "Newton vs Brent position");
        assert_vec_close(newton.v1, brent.v1, 1e-9, "Newton vs Brent velocity");
        assert_vec_close(newton.r1, auto.r1, 1e-9, "Newton vs Auto position");
        assert_vec_close(newton.v1, auto.v1, 1e-9, "Newton vs Auto velocity");
    }
}

/// Property-based invariants for [`propagate_universal`].
///
/// Unlike [`tests_propagation_edge_cases`], these tests require no external
/// ground truth: they check physical/mathematical invariants that must hold
/// for *any* valid two-body state, sampled broadly across the orbital
/// regimes and time-of-flight gaps a Kalman-filter-driven pipeline (fink-fat,
/// fed by survey cadences such as ZTF/LSST) actually encounters — including
/// gaps well beyond 30 days between observations, not just short arcs.
#[cfg(test)]
mod tests_propagation_invariants {
    use super::*;

    use crate::kepler::{
        prelim_kepler::prelim_parabolic::ParabolicPrelimMethod, SolverKind, SolverParams,
    };
    use nalgebra::Rotation3;
    use proptest::prelude::*;

    const T0: f64 = 60000.0;

    fn default_solver_type() -> SolverType {
        SolverType {
            kind: SolverKind::Auto,
            params: SolverParams {
                convergency: 2.220446049250313e-14,
                psi_guess: None,
                max_iter_prelim_kepuni: 20,
                parabolic_solving_method: ParabolicPrelimMethod::Cardano,
            },
        }
    }

    /// Independent reimplementation of the specific-orbital-energy formula
    /// (deliberately *not* reusing the private `initial_orbital_state`, so
    /// this genuinely checks conservation rather than tautologically
    /// re-deriving the same number the production code already computed).
    fn specific_energy_alpha(position: &Vector3<f64>, velocity: &Vector3<f64>) -> f64 {
        let mu = GAUSS_GRAV * GAUSS_GRAV;
        let r = position.norm();
        let v2 = velocity.norm_squared();
        (v2 - 2.0 * mu / r) / mu
    }

    /// Generate a physically plausible heliocentric `(position, velocity)`
    /// pair by sampling classical orbital elements `(q, e, true_anomaly, i,
    /// Ω, ω)` and converting to Cartesian state via the standard
    /// perifocal-frame conic formulas — valid uniformly across elliptic,
    /// near-parabolic, and hyperbolic regimes (no Kepler-equation solve
    /// needed, since the state is parametrised directly by true anomaly).
    ///
    /// `q` (periapsis distance) spans [0.05, 150] AU and `e` spans [0, 5),
    /// matching the plausibility bounds already used for IOD
    /// (`IODParams::r2_min_au`/`r2_max_au`/`max_ecc`). The resulting `r0` is
    /// filtered back into a comparable [0.02, 500] AU range, since a raw
    /// true-anomaly sweep on a near-parabolic/hyperbolic branch can
    /// otherwise wander arbitrarily far from the Sun.
    fn arb_heliocentric_state() -> impl Strategy<Value = (Vector3<f64>, Vector3<f64>)> {
        (
            0.05f64..150.0,                // periapsis distance q (AU)
            0.0f64..5.0,                   // target eccentricity e
            0.0f64..1.0,                   // true-anomaly fraction
            0.0f64..std::f64::consts::PI,  // inclination
            0.0f64..std::f64::consts::TAU, // ascending node longitude
            0.0f64..std::f64::consts::TAU, // argument of periapsis
        )
            .prop_map(|(q, e, nu_frac, inclination, node, periapsis_arg)| {
                let mu = GAUSS_GRAV * GAUSS_GRAV;

                // True anomaly: unrestricted for ellipses; clamped to 90% of
                // the asymptotic angle for parabola/hyperbola so `r` stays
                // finite.
                let true_anomaly = if e < 1.0 {
                    nu_frac * std::f64::consts::TAU - std::f64::consts::PI
                } else {
                    let nu_max = (-1.0 / e).acos() * 0.9;
                    (nu_frac * 2.0 - 1.0) * nu_max
                };

                let semi_latus_rectum = q * (1.0 + e);
                let h = (mu * semi_latus_rectum).sqrt();
                let radius = semi_latus_rectum / (1.0 + e * true_anomaly.cos());

                let position_perifocal = Vector3::new(
                    radius * true_anomaly.cos(),
                    radius * true_anomaly.sin(),
                    0.0,
                );
                let velocity_perifocal =
                    (mu / h) * Vector3::new(-true_anomaly.sin(), e + true_anomaly.cos(), 0.0);

                let rotation = Rotation3::from_axis_angle(&Vector3::z_axis(), node)
                    * Rotation3::from_axis_angle(&Vector3::x_axis(), inclination)
                    * Rotation3::from_axis_angle(&Vector3::z_axis(), periapsis_arg);

                (rotation * position_perifocal, rotation * velocity_perifocal)
            })
            .prop_filter(
                "r0 within a plausible heliocentric range",
                |(position, _)| (0.02..500.0).contains(&position.norm()),
            )
    }

    /// Time-of-flight strategy weighted to reflect real survey cadence
    /// (ZTF/LSST) rather than an arbitrary uniform range: intra-night
    /// timescales, inter-night arcs, inter-lunation gaps (the concrete
    /// "more than 30 days between detections" case), and inter-season/
    /// inter-opposition recoveries. Both signs are included, since the
    /// Kalman filter also performs backward (smoothing) corrections.
    fn arb_time_of_flight() -> impl Strategy<Value = f64> {
        prop_oneof![
            2 => -0.05f64..0.05,
            3 => -10.0f64..10.0,
            3 => -45.0f64..-20.0,
            3 => 20.0f64..45.0,
            2 => -400.0f64..-90.0,
            2 => 90.0f64..400.0,
        ]
    }

    proptest! {
        /// `propagate_universal` must never panic, regardless of orbital
        /// regime or gap length: it must always resolve to `Ok`/`Err`, and
        /// every field of an `Ok` result must be finite.
        #[test]
        fn prop_no_panic_ever(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            let result = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type());
            if let Ok(r) = result {
                prop_assert!(r.r1.iter().all(|x| x.is_finite()));
                prop_assert!(r.v1.iter().all(|x| x.is_finite()));
                prop_assert!(r.f_lag.is_finite());
                prop_assert!(r.g_lag.is_finite());
                prop_assert!(r.f_dot.is_finite());
                prop_assert!(r.g_dot.is_finite());
                prop_assert!(r.psy.is_finite());
            }
        }

        /// The Lagrange identity `f * g_dot - g * f_dot = 1` must hold to
        /// near machine precision for every converged propagation.
        #[test]
        fn prop_lagrange_identity(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            if let Ok(r) = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type()) {
                let identity = r.f_lag * r.g_dot - r.g_lag * r.f_dot;
                prop_assert!(
                    (identity - 1.0).abs() < 1e-8,
                    "Lagrange identity violated: f*gdot - g*fdot = {identity}"
                );
            }
        }

        /// Specific orbital energy (hence `alpha`) is conserved by two-body
        /// dynamics: recomputing it from `(r1, v1)` must match the value
        /// from `(r0, v0)`.
        #[test]
        fn prop_energy_conserved(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            if let Ok(r) = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type()) {
                let alpha_before = specific_energy_alpha(&position, &velocity);
                let alpha_after = specific_energy_alpha(&r.r1, &r.v1);
                let scale = 1.0 + alpha_before.abs();
                prop_assert!(
                    (alpha_after - alpha_before).abs() < 1e-8 * scale,
                    "energy not conserved: before={alpha_before}, after={alpha_after}"
                );
            }
        }

        /// Angular momentum norm `|r x v|` is conserved by two-body
        /// dynamics (planar, central-force motion).
        #[test]
        fn prop_angular_momentum_conserved(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            if let Ok(r) = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type()) {
                let h_before = position.cross(&velocity).norm();
                let h_after = r.r1.cross(&r.v1).norm();
                let scale = 1.0 + h_before;
                prop_assert!(
                    (h_after - h_before).abs() < 1e-8 * scale,
                    "angular momentum not conserved: before={h_before}, after={h_after}"
                );
            }
        }

        /// Propagating forward by `dt` then backward by the same `dt` must
        /// recover the original state. This is the strongest invariant here:
        /// it needs no external ground truth and directly catches unit/scale
        /// regressions like the one fixed alongside `test_propag`.
        #[test]
        fn prop_round_trip(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            if let Ok(fwd) = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type()) {
                if let Ok(back) = propagate_universal(&fwd.r1, &fwd.v1, T0 + dt, T0, default_solver_type()) {
                    let position_scale = 1.0 + position.norm();
                    let velocity_scale = 1.0 + velocity.norm();
                    prop_assert!(
                        (back.r1 - position).norm() < 1e-7 * position_scale,
                        "round-trip position mismatch: {:?} vs {:?}", back.r1, position
                    );
                    prop_assert!(
                        (back.v1 - velocity).norm() < 1e-7 * velocity_scale,
                        "round-trip velocity mismatch: {:?} vs {:?}", back.v1, velocity
                    );
                }
            }
        }

        /// When both solvers converge, `NewtonRaphson` and `BrentDecker`
        /// must agree on the propagated state.
        #[test]
        fn prop_solver_kind_agreement(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            let mut newton_solver_type = default_solver_type();
            newton_solver_type.kind = SolverKind::NewtonRaphson;
            let mut brent_solver_type = default_solver_type();
            brent_solver_type.kind = SolverKind::BrentDecker;

            let newton = propagate_universal(&position, &velocity, T0, T0 + dt, newton_solver_type);
            let brent = propagate_universal(&position, &velocity, T0, T0 + dt, brent_solver_type);

            if let (Ok(rn), Ok(rb)) = (newton, brent) {
                let scale = 1.0 + rn.r1.norm();
                prop_assert!((rn.r1 - rb.r1).norm() < 1e-7 * scale);
                prop_assert!((rn.v1 - rb.v1).norm() < 1e-7 * scale);
            }
        }

        /// Chaining several small, irregular propagation steps (mimicking a
        /// real sequence of survey detections) must land close to a single
        /// direct propagation over the same total time span.
        #[test]
        fn prop_chained_vs_single_step(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            let t1 = T0 + dt;
            // Irregular, non-uniform fractions of the total span (sums to 1.0).
            let fractions = [0.02, 0.28, 0.1, 0.6];
            let mut times = vec![T0];
            let mut acc = T0;
            for f in fractions {
                acc += dt * f;
                times.push(acc);
            }
            *times.last_mut().unwrap() = t1; // avoid float drift on the endpoint

            let mut pos = position;
            let mut vel = velocity;
            let mut chained_ok = true;
            for window in times.windows(2) {
                match propagate_universal(&pos, &vel, window[0], window[1], default_solver_type()) {
                    Ok(r) => {
                        pos = r.r1;
                        vel = r.v1;
                    }
                    Err(_) => {
                        chained_ok = false;
                        break;
                    }
                }
            }

            if chained_ok {
                if let Ok(direct) = propagate_universal(&position, &velocity, T0, t1, default_solver_type()) {
                    let position_scale = 1.0 + direct.r1.norm();
                    let velocity_scale = 1.0 + direct.v1.norm();
                    prop_assert!(
                        (pos - direct.r1).norm() < 1e-6 * position_scale,
                        "chained vs single-step position mismatch: {:?} vs {:?}", pos, direct.r1
                    );
                    prop_assert!(
                        (vel - direct.v1).norm() < 1e-6 * velocity_scale,
                        "chained vs single-step velocity mismatch: {:?} vs {:?}", vel, direct.v1
                    );
                }
            }
        }

        /// A tiny relative perturbation (`1e-8`) of the initial position
        /// must not produce a disproportionate change in the propagated
        /// position — guards against pathological ill-conditioning (e.g.
        /// near an internal solver-branch boundary).
        #[test]
        fn prop_continuity_position_perturbation(
            (position, velocity) in arb_heliocentric_state(),
            dt in arb_time_of_flight(),
        ) {
            const EPS: f64 = 1e-8;
            if let Ok(base) = propagate_universal(&position, &velocity, T0, T0 + dt, default_solver_type()) {
                let perturbed_position = position * (1.0 + EPS);
                if let Ok(perturbed) = propagate_universal(&perturbed_position, &velocity, T0, T0 + dt, default_solver_type()) {
                    let input_delta = (perturbed_position - position).norm();
                    let output_delta = (perturbed.r1 - base.r1).norm();
                    let output_scale = 1.0 + base.r1.norm();
                    prop_assert!(
                        output_delta < 1.0e6 * input_delta.max(1e-300) + 1e-9 * output_scale,
                        "position perturbation amplified disproportionately: input_delta={input_delta}, output_delta={output_delta}"
                    );
                }
            }
        }
    }
}
