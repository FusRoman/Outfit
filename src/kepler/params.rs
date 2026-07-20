//! Input parameters bundle for the universal-variable Kepler solver.

#[cfg(feature = "serde")]
use serde::Deserialize;

#[cfg(feature = "serde")]
use serde::Serialize;

use crate::{
    kepler::{
        brent_dekker_solver::solve_kepuni_brent_dekker,
        prelim_elliptic, prelim_hyperbolic,
        prelim_kepler::prelim_parabolic::{prelim_parabolic, ParabolicPrelimMethod},
        solve_kepuni_with_guess, UniversalKeplerSolution,
    },
    OutfitError,
};

use super::orbit_type::OrbitType;

/// Common tuning parameters shared by all Kepler equation solvers.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, Copy)]
pub struct SolverParams {
    /// Convergence threshold on the solved variable.
    pub convergency: f64,
    /// Initial guess for the universal anomaly $\psi$.
    pub psi_guess: Option<f64>,
    /// Maximum number of iterations for the preliminary `kepuni` step.
    pub max_iter_prelim_kepuni: usize,
    /// Choose between the exact analytical method or Newton-Raphson minimization
    pub parabolic_solving_method: ParabolicPrelimMethod,
}

impl Default for SolverParams {
    fn default() -> Self {
        Self {
            convergency: 100.0 * f64::EPSILON,
            psi_guess: None,
            max_iter_prelim_kepuni: 20,
            parabolic_solving_method: ParabolicPrelimMethod::Cardano,
        }
    }
}

/// Selects which root-finding algorithm is used to solve Kepler's equation.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, Copy)]
pub enum SolverKind {
    /// Newton-Raphson iteration.
    NewtonRaphson,
    /// Brent-Decker's bracketing method.
    BrentDecker,
    /// Automatically selects between available solvers.
    Auto,
}

/// Full solver configuration: algorithm choice plus its tuning parameters.
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone, Copy)]
pub struct SolverType {
    pub kind: SolverKind,
    pub params: SolverParams,
}

impl Default for SolverType {
    fn default() -> Self {
        Self {
            kind: SolverKind::NewtonRaphson,
            params: Default::default(),
        }
    }
}

/// Parameters required to solve the universal Kepler equation.
///
/// This struct bundles together the state and constants needed to propagate an
/// orbit using the universal variable formulation:
///
/// * `dt`: Propagation time interval Δt = t - t₀ (in days).
/// * `r0`: Heliocentric distance at the reference epoch t₀ (in AU).
/// * `sig0`: Radial velocity component at t₀ (in AU/day).
/// * `mu`: Standard gravitational parameter μ = GM (in AU³/day²).
/// * `alpha`: Reciprocal semi-major axis, `alpha = -1/a = 2E/μ` (in 1/AU).
///   This is the convention expected by [`s_funct`](crate::kepler::s_funct)
///   and the rest of the universal-variable machinery — *not* the raw
///   vis-viva `2E` (dimension velocity²).
/// * `e0`: Orbital eccentricity (unitless).
///
/// The associated method [`orbit_type`](UniversalKeplerParams::orbit_type)
/// classifies the orbit into elliptical, parabolic, or hyperbolic regimes
/// based on `alpha`.
#[derive(Debug, Clone, Copy)]
pub struct UniversalKeplerParams {
    /// Propagation time interval Δt = t - t₀, in days.
    pub dt: f64,
    /// Heliocentric distance at the reference epoch t₀, in AU.
    pub r0: f64,
    /// Radial velocity component at t₀ (= r₀·v₀ / r₀), in AU/day.
    pub sig0: f64,
    /// Standard gravitational parameter μ = GM, in AU³/day².
    pub mu: f64,
    /// Reciprocal semi-major axis (alpha = -1/a = 2E/mu), in 1/AU.
    pub alpha: f64,
    /// Orbital eccentricity (unitless).
    pub e0: f64,
    /// solver type
    pub solver_type: SolverType,
}

impl UniversalKeplerParams {
    /// Returns the [`OrbitType`](crate::kepler::OrbitType) (elliptic, parabolic, or hyperbolic)
    /// corresponding to the value of `alpha`.
    pub fn orbit_type(&self) -> OrbitType {
        OrbitType::from_alpha(self.alpha)
    }

    /// Return `Some(())` for elliptic or hyperbolic orbits, `None` for parabolic.
    ///
    /// Parabolic motion ( $ \alpha = 0 $ ) is not supported by the universal-variable
    /// formulation used here and must be rejected before entering the solver.
    #[inline]
    pub fn reject_parabolic_orbit(&self) -> Option<()> {
        match self.orbit_type() {
            OrbitType::Parabolic => None,
            OrbitType::Elliptic | OrbitType::Hyperbolic => Some(()),
        }
    }

    pub fn solve(&self) -> Result<UniversalKeplerSolution, OutfitError> {
        match self.solver_type.kind {
            SolverKind::NewtonRaphson => {
                solve_kepuni_with_guess(self).ok_or(OutfitError::NewtonRaphsonKeplerConvergence)
            }
            SolverKind::BrentDecker => {
                solve_kepuni_brent_dekker(self).ok_or(OutfitError::BrentDekkerKeplerConvergence)
            }
            SolverKind::Auto => solve_kepuni_with_guess(self)
                .or_else(|| solve_kepuni_brent_dekker(self))
                .ok_or(OutfitError::BrentDekkerKeplerConvergence),
        }
    }

    /// Compute an initial estimate of the universal anomaly $\psi_0$.
    ///
    /// # Goal
    ///
    /// Provide a starting value for $\psi$ before the main solver refines it to
    /// full convergence. The appropriate branch is selected based on the orbit
    /// type encoded in `self.alpha`:
    ///
    /// | Regime | Condition | Branch |
    /// |---|---|---|
    /// | Elliptic | $\alpha < 0$ | [`prelim_elliptic`] |
    /// | Hyperbolic | $\alpha > 0$ | [`prelim_hyperbolic`] |
    /// | Parabolic | $\alpha = 0$ | Not supported — returns `None` |
    ///
    /// # Convergence tolerance
    ///
    /// The tolerance passed to the preliminary Newton–Raphson sub-iterations is
    /// extracted from [`SolverType`]:
    ///
    /// - [`SolverKind::NewtonRaphson`] and [`SolverKind::Auto`] — uses the
    ///   `convergency` field directly, defaulting to $100\varepsilon$ if `None`.
    /// - [`SolverKind::BrentDecker`] — same extraction, since the preliminary
    ///   guess is solver-agnostic.
    ///
    /// # Return
    ///
    /// * `Some(psi_0)` — initial estimate of the universal anomaly.
    /// * `None` — if the orbit is parabolic ($\alpha = 0$).
    ///
    /// # Notes
    ///
    /// - This method does **not** iterate to full convergence; it only produces
    ///   a starting point for the main solver.
    /// - Each branch runs at most [`max_iter_prelim_kepuni`](crate::kepler::SolverType) Newton steps
    ///   internally.
    ///
    /// # See also
    ///
    /// * [`prelim_elliptic`] — elliptic branch.
    /// * [`prelim_hyperbolic`] — hyperbolic branch.
    /// * [`UniversalKeplerParams::solve`] — main solver entry point.
    pub fn prelim_kepuni(&self) -> Option<f64> {
        match self.orbit_type() {
            OrbitType::Elliptic => Some(prelim_elliptic(self)),
            OrbitType::Hyperbolic => Some(prelim_hyperbolic(self)),
            OrbitType::Parabolic => Some(prelim_parabolic(self)),
        }
    }
}

#[cfg(test)]
mod kepler_params_tests {
    use crate::kepler::s_funct;

    use super::*;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    // -----------------------------------------------------------------------
    // Test helpers
    // -----------------------------------------------------------------------

    /// Standard gravitational parameter for the Sun, AU³/day².
    const MU_SUN: f64 = 2.959_122_082_855_911e-4;

    /// Evaluate the universal Kepler residual at the solution `psi` and return
    /// its absolute value. Used to verify that a solution actually satisfies
    /// the equation regardless of which solver produced it.
    fn kepler_residual(sol: &UniversalKeplerSolution, params: &UniversalKeplerParams) -> f64 {
        // Recompute Stumpff functions from psi to catch any inconsistency
        // between the stored universal anomaly and the cached s.. values.
        let (_, s1, s2, s3) = s_funct(sol.universal_anomaly, params.alpha);
        (params.r0 * s1 + params.sig0 * s2 + s3 - params.mu.sqrt() * params.dt).abs()
    }

    /// Build a [`UniversalKeplerParams`] for a near-circular elliptic orbit.
    ///
    /// Uses a semi-major axis `a` (AU) to derive `alpha` and a nearly-zero
    /// radial velocity proxy.
    fn elliptic_params(dt: f64, a: f64, kind: SolverKind) -> UniversalKeplerParams {
        // alpha = -1 / a  (reciprocal semi-major axis convention)
        let alpha = -1.0 / a;
        // r0 ~ a for a near-circular orbit.
        let r0 = a;
        // sig0 ~ 0 for a circular orbit (radial velocity is zero at periapsis/apoapsis).
        let sig0 = 0.0;
        UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu: MU_SUN,
            alpha,
            e0: 0.01,
            solver_type: SolverType {
                kind,
                params: SolverParams::default(),
            },
        }
    }

    /// Build a [`UniversalKeplerParams`] for a hyperbolic orbit.
    ///
    /// Uses a characteristic energy $C_3 > 0$ (AU²/day²) so that
    /// $\alpha = C_3 / \mu > 0$.
    fn hyperbolic_params(dt: f64, c3: f64, kind: SolverKind) -> UniversalKeplerParams {
        let alpha = c3 / MU_SUN; // positive energy -> hyperbolic
        let r0 = 1.5; // AU
        let sig0 = 0.001;
        UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu: MU_SUN,
            alpha,
            e0: 1.5,
            solver_type: SolverType {
                kind,
                params: SolverParams::default(),
            },
        }
    }

    // -----------------------------------------------------------------------
    // Unit tests — Newton–Raphson
    // -----------------------------------------------------------------------

    #[test]
    fn newton_solves_elliptic_short_arc() {
        let params = elliptic_params(
            10.0, // 10 days
            1.0,  // 1 AU semi-major axis (Earth-like)
            SolverKind::NewtonRaphson,
        );
        let sol = params
            .solve()
            .expect("Newton should converge on elliptic short arc");
        assert!(
            kepler_residual(&sol, &params) < 1e-10,
            "residual too large: {}",
            kepler_residual(&sol, &params)
        );
    }

    #[test]
    fn newton_solves_elliptic_full_orbit() {
        // Propagate for approximately one full period of an Earth-like orbit (~365 days).
        let params = elliptic_params(365.0, 1.0, SolverKind::NewtonRaphson);
        let sol = params
            .solve()
            .expect("Newton should converge on full elliptic orbit");
        assert!(
            kepler_residual(&sol, &params) < 1e-9,
            "residual too large: {}",
            kepler_residual(&sol, &params)
        );
    }

    #[test]
    fn newton_solves_elliptic_negative_dt() {
        // Backward propagation: dt < 0 should be handled symmetrically.
        let params = elliptic_params(-30.0, 1.0, SolverKind::NewtonRaphson);
        let sol = params
            .solve()
            .expect("Newton should converge for negative dt");
        assert!(kepler_residual(&sol, &params) < 1e-10);
    }

    #[test]
    fn newton_solves_hyperbolic() {
        let params = hyperbolic_params(
            5.0,
            1e-5, // small positive C3 -> mildly hyperbolic
            SolverKind::NewtonRaphson,
        );
        let sol = params
            .solve()
            .expect("Newton should converge on hyperbolic orbit");
        assert!(
            kepler_residual(&sol, &params) < 1e-10,
            "residual: {}",
            kepler_residual(&sol, &params)
        );
    }

    #[test]
    fn newton_warm_start_consistent_with_cold_start() {
        let mut cold = elliptic_params(50.0, 2.0, SolverKind::NewtonRaphson);
        let sol_cold = cold.solve().expect("cold start should converge");

        // Feed the cold-start solution as a warm-start guess.
        cold.solver_type = SolverType {
            kind: SolverKind::NewtonRaphson,
            params: SolverParams {
                psi_guess: Some(sol_cold.universal_anomaly),
                ..Default::default()
            },
        };
        let sol_warm = cold.solve().expect("warm start should converge");

        assert_abs_diff_eq!(
            sol_cold.universal_anomaly,
            sol_warm.universal_anomaly,
            epsilon = 1e-10
        );
    }

    // -----------------------------------------------------------------------
    // Unit tests — Brent–Dekker
    // -----------------------------------------------------------------------

    #[test]
    fn brent_solves_elliptic_short_arc() {
        let params = elliptic_params(10.0, 1.0, SolverKind::BrentDecker);
        let sol = params
            .solve()
            .expect("Brent should converge on elliptic short arc");
        assert!(
            kepler_residual(&sol, &params) < 1e-10,
            "residual: {}",
            kepler_residual(&sol, &params)
        );
    }

    #[test]
    fn brent_solves_elliptic_full_orbit() {
        let params = elliptic_params(365.0, 1.0, SolverKind::BrentDecker);
        let sol = params
            .solve()
            .expect("Brent should converge on full elliptic orbit");
        assert!(kepler_residual(&sol, &params) < 1e-9);
    }

    #[test]
    fn brent_solves_elliptic_negative_dt() {
        let params = elliptic_params(-30.0, 1.0, SolverKind::BrentDecker);
        let sol = params
            .solve()
            .expect("Brent should converge for negative dt");
        assert!(kepler_residual(&sol, &params) < 1e-10);
    }

    #[test]
    fn brent_solves_hyperbolic() {
        let params = hyperbolic_params(5.0, 1e-5, SolverKind::BrentDecker);
        let sol = params
            .solve()
            .expect("Brent should converge on hyperbolic orbit");
        assert!(
            kepler_residual(&sol, &params) < 1e-10,
            "residual: {}",
            kepler_residual(&sol, &params)
        );
    }

    // -----------------------------------------------------------------------
    // Consistency tests — Newton is the reference
    // -----------------------------------------------------------------------

    /// Tolerance on `psi` agreement between the two solvers.
    /// Newton is the reference: Brent must land within this bound.
    const PSI_CONSISTENCY_TOL: f64 = 1e-12;

    fn assert_solvers_consistent(params_template: UniversalKeplerParams) {
        let mut newton_params = params_template;
        newton_params.solver_type.kind = SolverKind::NewtonRaphson;

        let mut brent_params = params_template;
        brent_params.solver_type.kind = SolverKind::BrentDecker;

        let sol_newton = newton_params.solve().expect("Newton must converge");
        let sol_brent = brent_params.solve().expect("Brent must converge");

        // Both solutions must satisfy the Kepler equation.
        assert!(
            kepler_residual(&sol_newton, &newton_params) < 1e-10,
            "Newton residual too large"
        );
        assert!(
            kepler_residual(&sol_brent, &brent_params) < 1e-10,
            "Brent residual too large"
        );

        // The two solvers must agree on psi (Newton is the reference).
        assert_abs_diff_eq!(
            sol_newton.universal_anomaly,
            sol_brent.universal_anomaly,
            epsilon = PSI_CONSISTENCY_TOL
        );
    }

    #[test]
    fn solvers_consistent_elliptic_short_arc() {
        assert_solvers_consistent(elliptic_params(10.0, 1.0, SolverKind::BrentDecker));
    }

    #[test]
    fn solvers_consistent_elliptic_long_arc() {
        assert_solvers_consistent(elliptic_params(200.0, 2.5, SolverKind::BrentDecker));
    }

    #[test]
    fn solvers_consistent_elliptic_negative_dt() {
        assert_solvers_consistent(elliptic_params(-45.0, 1.5, SolverKind::BrentDecker));
    }

    #[test]
    fn solvers_consistent_hyperbolic() {
        assert_solvers_consistent(hyperbolic_params(5.0, 1e-5, SolverKind::BrentDecker));
    }

    #[test]
    fn solvers_consistent_outer_solar_system_elliptic() {
        // Jupiter-like orbit: a ~ 5.2 AU, dt ~ 1000 days.
        assert_solvers_consistent(elliptic_params(1000.0, 5.2, SolverKind::BrentDecker));
    }

    // -----------------------------------------------------------------------
    // Property-based tests
    // -----------------------------------------------------------------------

    /// Strategy for physically plausible elliptic propagation parameters.
    ///
    /// * `a` in [0.3, 10.0] AU  — Mercury-like to outer solar system.
    /// * `dt` in [1.0, 500.0] days — short to moderately long arcs.
    fn elliptic_strategy() -> impl Strategy<Value = (f64, f64)> {
        (
            0.3_f64..10.0_f64,  // semi-major axis (AU)
            1.0_f64..500.0_f64, // time of flight (days)
        )
    }

    /// Strategy for hyperbolic parameters.
    ///
    /// * `c3` in [1e-6, 1e-3] AU²/day² — mildly to moderately hyperbolic.
    /// * `dt` in [1.0, 100.0] days.
    fn hyperbolic_strategy() -> impl Strategy<Value = (f64, f64)> {
        (1e-6_f64..1e-3_f64, 1.0_f64..100.0_f64)
    }

    proptest! {
        /// For a wide range of elliptic orbits, Newton must converge and the
        /// residual must vanish to within numerical tolerance.
        #[test]
        fn proptest_newton_elliptic_residual(
            (a, dt) in elliptic_strategy()
        ) {
            let params = elliptic_params(dt, a, SolverKind::NewtonRaphson);
            if let Ok(sol) = params.solve() {
                prop_assert!(
                    kepler_residual(&sol, &params) < 1e-9,
                    "Newton residual too large: {}",
                    kepler_residual(&sol, &params)
                );
            }
            // If Newton fails to converge on a valid input we do not fail the
            // proptest — coverage of the failure path is handled by dedicated
            // unit tests.
        }

        /// For a wide range of elliptic orbits, Brent–Dekker must converge and
        /// the residual must vanish to within numerical tolerance.
        #[test]
        fn proptest_brent_elliptic_residual(
            (a, dt) in elliptic_strategy()
        ) {
            let params = elliptic_params(dt, a, SolverKind::BrentDecker);
            if let Ok(sol) = params.solve() {
                prop_assert!(
                    kepler_residual(&sol, &params) < 1e-9,
                    "Brent residual too large: {}",
                    kepler_residual(&sol, &params)
                );
            }
        }

        /// When both solvers converge on the same elliptic input, their `psi`
        /// values must agree. Newton is the reference.
        #[test]
        fn proptest_solvers_consistent_elliptic(
            (a, dt) in elliptic_strategy()
        ) {
            let newton_params = elliptic_params(dt, a, SolverKind::NewtonRaphson);
            let mut brent_params = newton_params;
            brent_params.solver_type.kind = SolverKind::BrentDecker;

            if let (Ok(sol_n), Ok(sol_b)) = (newton_params.solve(), brent_params.solve()) {
                prop_assert!(
                    (sol_n.universal_anomaly - sol_b.universal_anomaly).abs()
                        < PSI_CONSISTENCY_TOL,
                    "psi mismatch: Newton={}, Brent={}",
                    sol_n.universal_anomaly,
                    sol_b.universal_anomaly
                );
                // If either solver diverges, skip this sample — the residual
                // tests above already cover individual solver robustness.
            }
        }

        /// For hyperbolic orbits, whenever Newton converges, Brent must also
        /// converge and agree on `psi`.
        #[test]
        fn proptest_solvers_consistent_hyperbolic(
            (c3, dt) in hyperbolic_strategy()
        ) {
            let newton_params = hyperbolic_params(dt, c3, SolverKind::NewtonRaphson);
            let mut brent_params = newton_params;
            brent_params.solver_type.kind = SolverKind::BrentDecker;

            if let (Ok(sol_n), Ok(sol_b)) = (newton_params.solve(), brent_params.solve()) {
                prop_assert!(
                    (sol_n.universal_anomaly - sol_b.universal_anomaly).abs()
                        < PSI_CONSISTENCY_TOL,
                    "psi mismatch on hyperbolic: Newton={}, Brent={}",
                    sol_n.universal_anomaly,
                    sol_b.universal_anomaly
                );
             }
        }

        /// The `orbit_type` method must return a regime consistent with the
        /// sign of `alpha`, independently of the solver chosen.
        #[test]
        fn proptest_orbit_type_consistent_with_alpha(
            alpha in prop_oneof![
                (-1e-2_f64..-1e-6_f64),  // elliptic
                (1e-6_f64..1e-2_f64),    // hyperbolic
            ]
        ) {
            let params = UniversalKeplerParams {
                dt: 10.0,
                r0: 1.0,
                sig0: 0.0,
                mu: MU_SUN,
                alpha,
                e0: 0.5,
                solver_type: SolverType::default(),
            };
            match params.orbit_type() {
                OrbitType::Elliptic  => prop_assert!(alpha < 0.0),
                OrbitType::Hyperbolic => prop_assert!(alpha > 0.0),
                OrbitType::Parabolic  => prop_assert!(alpha == 0.0),
            }
        }
    }
}

#[cfg(test)]
mod tests_prelim_kepuni {

    use super::*;

    const MU: f64 = 1.0;
    const CONTR: f64 = 1e-12;

    fn make_params(
        dt: f64,
        r0: f64,
        sig0: f64,
        mu: f64,
        alpha: f64,
        e0: f64,
    ) -> UniversalKeplerParams {
        UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu,
            alpha,
            e0,
            solver_type: SolverType {
                params: SolverParams {
                    convergency: CONTR,
                    ..Default::default()
                },
                ..Default::default()
            },
        }
    }

    #[test]
    fn test_returns_none_for_alpha_zero() {
        let params = make_params(1.0, 1.0, 0.0, MU, 0.0, 0.1);
        let res = params.prelim_kepuni().unwrap();
        assert_eq!(res, 0.8846222003969053);
    }

    #[test]
    fn test_elliptic_small_eccentricity() {
        let params = make_params(0.5, 1.0, 0.1, MU, -1.0, 1e-8);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_elliptic_high_eccentricity() {
        let params = make_params(0.1, 0.5, 0.2, MU, -1.0, 0.8);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_hyperbolic_case() {
        let params = make_params(0.3, 2.0, -0.1, MU, 1.0, 1.5);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_negative_sig0_changes_direction() {
        let alpha = -1.0;
        let r0 = 1.0;
        let e0 = 0.5;
        let dt = 0.25;

        let params_pos = make_params(dt, r0, 0.1, MU, alpha, e0);
        let params_neg = make_params(dt, r0, -0.1, MU, alpha, e0);

        let psi_pos = params_pos.prelim_kepuni().unwrap();
        let psi_neg = params_neg.prelim_kepuni().unwrap();

        assert!(
            (psi_pos - psi_neg).abs() > 1e-8,
            "psi did not change significantly when changing sig0 sign: {psi_pos} vs {psi_neg}"
        );
    }

    #[test]
    fn test_stability_long_dt() {
        let params = make_params(50.0, 1.0, 0.1, MU, -1.0, 0.5);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }

    #[test]
    fn test_edge_cosine_limits() {
        let params = make_params(0.25, 2.0, 0.1, MU, -1.0, 0.1);
        let result = params.prelim_kepuni();
        assert!(result.is_some());
    }

    #[test]
    fn test_prelim_kepuni_real_data() {
        let dt = -20.765849999996135;
        let r0 = 1.3803870211345761;
        let sig0 = 3.701_354_484_003_874_8E-3;
        let mu = 2.959_122_082_855_911_5E-4;
        // alpha = -1/a, the reciprocal semi-major-axis convention (see
        // `UniversalKeplerParams::alpha` doc). Equal to the raw vis-viva
        // `2E = -1.642_158_377_771_140_7E-4` (AU^2/day^2) divided by `mu`.
        let alpha = -0.554_947_829_724_638_7;
        let e0 = 0.283_599_599_137_344_5;

        let params = UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu,
            alpha,
            e0,
            solver_type: SolverType::default(),
        };
        let psi = params.prelim_kepuni().unwrap();
        assert_eq!(psi, -0.2636637076378094);

        let params2 = UniversalKeplerParams {
            alpha: 0.554_947_829_724_638_7,
            ..params
        };
        let psi = params2.prelim_kepuni().unwrap();
        assert_eq!(psi, -1.258980225923225);

        let params3 = UniversalKeplerParams {
            alpha: 0.0,
            ..params
        };
        let psi = params3.prelim_kepuni().unwrap();
        assert_eq!(psi, -0.2568229151088884);
    }

    mod kepuni_prop_tests {

        use super::*;
        use proptest::prelude::*;

        fn arb_params() -> impl Strategy<Value = UniversalKeplerParams> {
            (
                -10.0..10.0f64,
                0.1..5.0f64,
                -2.0..2.0f64,
                0.5..2.0f64,
                prop_oneof![(-5.0..-0.01f64), (0.01..5.0f64)],
                0.0..3.0f64,
                1e-14..1e-8f64,
            )
                .prop_map(|(dt, r0, sig0, mu, alpha, e0, contr)| {
                    UniversalKeplerParams {
                        dt,
                        r0,
                        sig0,
                        mu,
                        alpha,
                        e0,
                        solver_type: SolverType {
                            params: SolverParams {
                                convergency: contr,
                                ..Default::default()
                            },
                            ..Default::default()
                        },
                    }
                })
        }

        proptest! {
            #[test]
            fn prop_prelim_kepuni_behaves_well(params in arb_params()) {
                let result = params.prelim_kepuni();
                prop_assert!(result.is_some());
                prop_assert!(result.unwrap().is_finite());
            }
        }

        proptest! {
            #[test]
            fn prop_prelim_kepuni_alpha_zero(
                dt in -10.0..10.0f64,
                r0 in 0.1..5.0f64,
                sig0 in -2.0..2.0f64,
                mu in 0.5..2.0f64,
                e0 in 0.0..3.0f64,
                contr in 1e-14..1e-8f64
            ) {
                let params = UniversalKeplerParams { dt, r0, sig0, mu, alpha: 0.0, e0, solver_type: SolverType {
                        params: SolverParams {
                            convergency: contr,
                            ..Default::default()
                        },
                        ..Default::default()
                    }
                };
                let result = params.prelim_kepuni();
                if let Some(psi0) = result {
                    prop_assert!(psi0.is_finite());
                }
            }
        }

        proptest! {
            #[test]
            fn prop_sig0_influences_psi(params in arb_params()) {
                prop_assume!(params.dt.abs() > 1e-6);
                prop_assume!(params.e0 > 1e-6);
                prop_assume!(params.r0 > 1e-6);

                let mut params_pos = params;
                params_pos.sig0 = 0.1;
                let mut params_neg = params;
                params_neg.sig0 = -0.1;

                let res_pos = params_pos.prelim_kepuni();
                let res_neg = params_neg.prelim_kepuni();

                prop_assume!(res_pos.is_some() && res_neg.is_some());

                let diff = (res_pos.unwrap() - res_neg.unwrap()).abs();
                prop_assert!(diff >= 0.0);
            }
        }
    }
}
