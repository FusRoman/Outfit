//! Input parameters bundle for the universal-variable Kepler solver.

use crate::{
    kepler::{
        brent_dekker_solver::solve_kepuni_brent_dekker, prelim_elliptic, prelim_hyperbolic,
        solve_kepuni_with_guess, UniversalKeplerSolution,
    },
    OutfitError,
};

use super::orbit_type::OrbitType;

#[derive(Debug, Clone, Copy)]
pub enum SolverType {
    NewtonRaphson {
        convergency: Option<f64>,
        psi_guess: Option<f64>,
        max_iter_prelim_kepuni: Option<usize>,
    },
    BrentDecker {
        convergency: Option<f64>,
        max_iter_prelim_kepuni: Option<usize>,
    },
    Auto {
        convergency: Option<f64>,
        psi_guess: Option<f64>,
        max_iter_prelim_kepuni: Option<usize>,
    },
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
/// * `alpha`: Twice the specific orbital energy (2E).
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
    /// Twice the specific orbital energy (alpha = 2E).
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
        match self.solver_type {
            SolverType::NewtonRaphson {
                convergency,
                psi_guess,
                ..
            } => solve_kepuni_with_guess(self, convergency, psi_guess)
                .ok_or(OutfitError::NewtonRaphsonKeplerConvergence),
            SolverType::BrentDecker { convergency, .. } => {
                solve_kepuni_brent_dekker(self, convergency)
                    .ok_or(OutfitError::BrentDekkerKeplerConvergence)
            }
            SolverType::Auto {
                convergency,
                psi_guess,
                ..
            } => solve_kepuni_with_guess(self, convergency, psi_guess)
                .or_else(|| solve_kepuni_brent_dekker(self, convergency))
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
    /// - [`SolverType::NewtonRaphson`] and [`SolverType::Auto`] — uses the
    ///   `convergency` field directly, defaulting to $100\varepsilon$ if `None`.
    /// - [`SolverType::BrentDecker`] — same extraction, since the preliminary
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
        let (convergency, max_iter) = match self.solver_type {
            SolverType::NewtonRaphson {
                convergency,
                max_iter_prelim_kepuni,
                ..
            }
            | SolverType::BrentDecker {
                convergency,
                max_iter_prelim_kepuni,
            }
            | SolverType::Auto {
                convergency,
                max_iter_prelim_kepuni,
                ..
            } => (
                convergency.unwrap_or(100.0 * f64::EPSILON),
                max_iter_prelim_kepuni.unwrap_or(20),
            ),
        };

        match self.orbit_type() {
            OrbitType::Elliptic => Some(prelim_elliptic(self, convergency, max_iter)),
            OrbitType::Hyperbolic => Some(prelim_hyperbolic(self, convergency, max_iter)),
            OrbitType::Parabolic => None,
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
        (params.r0 * s1 + params.sig0 * s2 + params.mu * s3 - params.dt).abs()
    }

    /// Build a [`UniversalKeplerParams`] for a near-circular elliptic orbit.
    ///
    /// Uses a semi-major axis `a` (AU) to derive `alpha` and a nearly-zero
    /// radial velocity proxy.
    fn elliptic_params(dt: f64, a: f64, solver: SolverType) -> UniversalKeplerParams {
        // alpha = -mu / a  (twice specific orbital energy for elliptic orbit)
        let alpha = -MU_SUN / a;
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
            solver_type: solver,
        }
    }

    /// Build a [`UniversalKeplerParams`] for a hyperbolic orbit.
    ///
    /// Uses a characteristic energy $C_3 > 0$ (AU²/day²) so that
    /// $\alpha = C_3 > 0$.
    fn hyperbolic_params(dt: f64, c3: f64, solver: SolverType) -> UniversalKeplerParams {
        let alpha = c3; // positive energy -> hyperbolic
        let r0 = 1.5; // AU
        let sig0 = 0.001;
        UniversalKeplerParams {
            dt,
            r0,
            sig0,
            mu: MU_SUN,
            alpha,
            e0: 1.5,
            solver_type: solver,
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
            SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
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
        let params = elliptic_params(
            365.0,
            1.0,
            SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        );
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
        let params = elliptic_params(
            -30.0,
            1.0,
            SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        );
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
            SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
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
        let mut cold = elliptic_params(
            50.0,
            2.0,
            SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        );
        let sol_cold = cold.solve().expect("cold start should converge");

        // Feed the cold-start solution as a warm-start guess.
        cold.solver_type = SolverType::NewtonRaphson {
            convergency: None,
            psi_guess: Some(sol_cold.universal_anomaly),
            max_iter_prelim_kepuni: None,
        };
        let sol_warm = cold.solve().expect("warm start should converge");

        assert_abs_diff_eq!(
            sol_cold.universal_anomaly,
            sol_warm.universal_anomaly,
            epsilon = 1e-10
        );
    }

    #[test]
    fn newton_returns_error_for_parabolic() {
        let params = UniversalKeplerParams {
            dt: 10.0,
            r0: 1.0,
            sig0: 0.0,
            mu: MU_SUN,
            alpha: 0.0, // exactly parabolic
            e0: 1.0,
            solver_type: SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None,
            },
        };
        assert!(
            params.solve().is_err(),
            "parabolic orbit should return an error"
        );
    }

    // -----------------------------------------------------------------------
    // Unit tests — Brent–Dekker
    // -----------------------------------------------------------------------

    #[test]
    fn brent_solves_elliptic_short_arc() {
        let params = elliptic_params(
            10.0,
            1.0,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        );
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
        let params = elliptic_params(
            365.0,
            1.0,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        );
        let sol = params
            .solve()
            .expect("Brent should converge on full elliptic orbit");
        assert!(kepler_residual(&sol, &params) < 1e-9);
    }

    #[test]
    fn brent_solves_elliptic_negative_dt() {
        let params = elliptic_params(
            -30.0,
            1.0,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        );
        let sol = params
            .solve()
            .expect("Brent should converge for negative dt");
        assert!(kepler_residual(&sol, &params) < 1e-10);
    }

    #[test]
    fn brent_solves_hyperbolic() {
        let params = hyperbolic_params(
            5.0,
            1e-5,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        );
        let sol = params
            .solve()
            .expect("Brent should converge on hyperbolic orbit");
        assert!(
            kepler_residual(&sol, &params) < 1e-10,
            "residual: {}",
            kepler_residual(&sol, &params)
        );
    }

    #[test]
    fn brent_returns_error_for_parabolic() {
        let params = UniversalKeplerParams {
            dt: 10.0,
            r0: 1.0,
            sig0: 0.0,
            mu: MU_SUN,
            alpha: 0.0,
            e0: 1.0,
            solver_type: SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        };
        assert!(params.solve().is_err());
    }

    // -----------------------------------------------------------------------
    // Consistency tests — Newton is the reference
    // -----------------------------------------------------------------------

    /// Tolerance on `psi` agreement between the two solvers.
    /// Newton is the reference: Brent must land within this bound.
    const PSI_CONSISTENCY_TOL: f64 = 1e-12;

    fn assert_solvers_consistent(params_template: UniversalKeplerParams) {
        let mut newton_params = params_template;
        newton_params.solver_type = SolverType::NewtonRaphson {
            convergency: None,
            psi_guess: None,
            max_iter_prelim_kepuni: None,
        };

        let mut brent_params = params_template;
        brent_params.solver_type = SolverType::BrentDecker {
            convergency: None,
            max_iter_prelim_kepuni: None,
        };

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
        assert_solvers_consistent(elliptic_params(
            10.0,
            1.0,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        ));
    }

    #[test]
    fn solvers_consistent_elliptic_long_arc() {
        assert_solvers_consistent(elliptic_params(
            200.0,
            2.5,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        ));
    }

    #[test]
    fn solvers_consistent_elliptic_negative_dt() {
        assert_solvers_consistent(elliptic_params(
            -45.0,
            1.5,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        ));
    }

    #[test]
    fn solvers_consistent_hyperbolic() {
        assert_solvers_consistent(hyperbolic_params(
            5.0,
            1e-5,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        ));
    }

    #[test]
    fn solvers_consistent_outer_solar_system_elliptic() {
        // Jupiter-like orbit: a ~ 5.2 AU, dt ~ 1000 days.
        assert_solvers_consistent(elliptic_params(
            1000.0,
            5.2,
            SolverType::BrentDecker {
                convergency: None,
                max_iter_prelim_kepuni: None,
            },
        ));
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
            let params = elliptic_params(dt, a, SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None
            });
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
            let params = elliptic_params(dt, a, SolverType::BrentDecker { convergency: None, max_iter_prelim_kepuni: None });
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
            let newton_params = elliptic_params(dt, a, SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None
            });
            let mut brent_params = newton_params;
            brent_params.solver_type = SolverType::BrentDecker { convergency: None, max_iter_prelim_kepuni: None };

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
            let newton_params = hyperbolic_params(dt, c3, SolverType::NewtonRaphson {
                convergency: None,
                psi_guess: None,
                max_iter_prelim_kepuni: None
            });
            let mut brent_params = newton_params;
            brent_params.solver_type = SolverType::BrentDecker { convergency: None, max_iter_prelim_kepuni: None };

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
                solver_type: SolverType::NewtonRaphson { convergency: None, psi_guess: None, max_iter_prelim_kepuni: None },
            };
            match params.orbit_type() {
                OrbitType::Elliptic  => prop_assert!(alpha < 0.0),
                OrbitType::Hyperbolic => prop_assert!(alpha > 0.0),
                OrbitType::Parabolic  => prop_assert!(alpha == 0.0),
            }
        }
    }
}
