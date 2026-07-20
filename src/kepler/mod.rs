//! # Universal Kepler Propagation & f–g Velocity Update
//!
//! This module implements a **universal‐variable** formulation of the two‑body
//! problem to propagate orbital motion and to apply a **Lagrange f–g**
//! velocity correction between two position epochs.
//!
//! It provides:
//! - An **orbit classifier** [`OrbitType`](crate::kepler::OrbitType) from the energy sign `alpha = 2E`,
//! - A compact container of **universal Kepler inputs** [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams),
//! - A numerically robust evaluation of generalized **Stumpff‑like functions**
//!   via [`s_funct`](crate::kepler::s_funct),
//! - Practical **preliminary solvers** for the universal anomaly ψ in both
//!   elliptic and hyperbolic regimes ([`prelim_elliptic`](crate::kepler::prelim_elliptic), [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic),
//!   wrapped by [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni)),
//! - A Newton–Raphson **universal Kepler solver** [`solve_kepuni`](crate::kepler::solve_kepuni) that returns
//!   ψ and the Stumpff tuple (s0, s1, s2, s3),
//! - A high‑level **velocity correction** routine [`velocity_correction`](crate::kepler::velocity_correction) that
//!   computes `(v₂_corrected, f, g)` from two position vectors and an initial
//!   velocity.
//!
//! ## Units & Conventions
//! ----------------------
//! - Distance: **AU**
//! - Time: **days**
//! - Velocity: **AU/day**
//! - Gravitational parameter: `μ = GM = GAUSS_GRAV²` (AU³/day²)
//! - Angles: **radians** (normalized with [`principal_angle`](crate::kepler::principal_angle))
//! - Energy parameter: `alpha = 2 * (v²/2 − μ/r)`
//!
//! The module uses `nalgebra::Vector3<f64>` for state vectors.
//!
//! ## Mathematical Background
//! --------------------------
//! Universal Kepler’s equation is written with the generalized Stumpff
//! functions (here denoted `s0..s3`):
//!
//! ```text
//! r₀ s₁(ψ, α) + σ₀ s₂(ψ, α) + μ s₃(ψ, α) = Δt
//! with   σ₀ = r₀·ṙ₀ / r₀ = x₀·v₀
//! and    α = 2E  (twice specific orbital energy).
//! ```
//!
//! The invariants that must hold for any `(ψ, α)` are:
//!
//! ```text
//! s₀ = 1 + α s₂
//! s₁ = ψ + α s₃
//! ```
//!
//! Once ψ is known, the Lagrange coefficients for the f–g solution are:
//!
//! ```text
//! f = 1 − (μ/r₀) s₂   ,   g = Δt − μ s₃
//! ```
//!
//! which are used in [`velocity_correction`](crate::kepler::velocity_correction).
//!
//! ## Numerical Strategy
//! ---------------------
//! - [`s_funct`](crate::kepler::s_funct) evaluates `(s0..s3)` either by **power series** when
//!   `|β| = |α ψ²|` is small, or by **half‑angle reduction + duplication**
//!   when `|β|` is large, to avoid divergence.
//! - [`prelim_elliptic`](crate::kepler::prelim_elliptic) and [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) build robust initial guesses
//!   for ψ from the classical eccentric/hyperbolic anomalies (`u`, `F`) with
//!   Newton steps and sign handling from the radial motion (`sig0`).
//! - [`solve_kepuni`](crate::kepler::solve_kepuni) refines ψ with **Newton–Raphson**, includes sign‑change
//!   damping, and rejects obviously unstable configurations.
//! - **Parabolic (`alpha = 0`)** motion is **not supported** here and returns
//!   `None`; treat it with a dedicated routine if needed.
//!
//! ## API Overview
//! ---------------
//! - [`OrbitType`](crate::kepler::OrbitType): classify orbit from `alpha`.
//! - [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams): bundle of inputs `(dt, r0, sig0, mu, alpha, e0)`.
//! - [`s_funct`](crate::kepler::s_funct): compute `(s0, s1, s2, s3)` and satisfy `s0 = 1 + α s2`, `s1 = ψ + α s3`.
//! - [`principal_angle`](crate::kepler::principal_angle), [`angle_diff`](crate::kepler::angle_diff): angle normalization helpers.
//! - [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni): dispatch to [`prelim_elliptic`](crate::kepler::prelim_elliptic) / [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic) for ψ guess.
//! - [`solve_kepuni`](crate::kepler::solve_kepuni): universal Kepler solver → `(ψ, s0, s1, s2, s3)`.
//! - [`velocity_correction`](crate::kepler::velocity_correction): apply f–g to update velocity `(v₂_corrected, f, g)`.
//!
//! ## Error Handling
//! -----------------
//! - Parabolic case (`alpha = 0`) returns `None` from the preliminary/solver routines.
//! - [`velocity_correction`](crate::kepler::velocity_correction) returns `Result<… , OutfitError>` and can fail when:
//!   - the upstream [`eccentricity_control`](crate::orb_elem::eccentricity_control) rejects the orbit (e.g., zero angular momentum,
//!     perihelion/eccentricity bounds), or
//!   - [`solve_kepuni`](crate::kepler::solve_kepuni) does not converge.
//!
//! ## Testing & Invariants
//! -----------------------
//! The test suite covers:
//! - **Analytical limits:** `α = 0`, `ψ = 0`, symmetry in `ψ`,
//! - **Large |β|** branch consistency,
//! - **Residual checks:** `r₀ s₁ + σ₀ s₂ + μ s₃ − Δt ≈ 0`,
//! - **Property‑based tests** over wide parameter ranges,
//! - **“Real data”** vectors to guard against regressions.
//!
//! ## Performance Notes
//! --------------------
//! - The routines are allocation‑free and deterministic in `f64`.
//! - Tolerances default to `𝒪(ε)` with small safety factors; tune via the
//!   `contr`/`convergency` parameters when necessary.
//! - For extremely large time steps or highly hyperbolic cases, consider
//!   scaling states and/or tightening tolerances.
//!
//! ## Examples
//! -----------
//! ### Solve universal Kepler for ψ and Stumpff functions
//! ```rust, no_run
//! use nalgebra::Vector3;
//! use outfit::kepler::{UniversalKeplerParams, solve_kepuni, UniversalKeplerSolution, SolverType};
//!
//! let params = UniversalKeplerParams {
//!     dt: 0.25,          // days
//!     r0: 1.0,           // AU
//!     sig0: 0.0,         // AU/day
//!     mu: 0.00029591220828559115, // GAUSS_GRAV^2
//!     alpha: -1.0,       // elliptic
//!     e0: 0.1,
//!     solver_type: SolverType::default()
//! };
//!
//! let kepler_solution: UniversalKeplerSolution = solve_kepuni(&params)
//!     .expect("converged");
//! // Use (psi, s0..s3) to build f, g or propagate state.
//! ```
//!
//! ### Apply f–g velocity correction between two position epochs
//! ```rust
//! use nalgebra::Vector3;
//! use outfit::kepler::velocity_correction;
//!
//! let x1 = Vector3::new(1.0, 0.0, 0.0);   // r(t1) in AU
//! let x2 = Vector3::new(1.1, 0.0, 0.0);   // r(t2) in AU
//! let v2 = Vector3::new(0.0, 0.017, 0.0); // v(t2) in AU/day
//! let dt = 1.0;                           // t2 - t1 in days
//!
//! let (v2_corr, f, g) = velocity_correction(&x1, &x2, &v2, dt, 5.0, 0.9, 1e-12)?;
//! // v2_corr is the corrected velocity at t2 using the universal-variable f–g solution.
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## See also
//! ------------
//! * [`velocity_correction`](crate::kepler::velocity_correction) – Lagrange‑based velocity update.
//! * [`solve_kepuni`](crate::kepler::solve_kepuni) – Universal Kepler solver returning (ψ, s0..s3).
//! * [`prelim_elliptic`](crate::kepler::prelim_elliptic), [`prelim_hyperbolic`](crate::kepler::prelim_hyperbolic), [`prelim_kepuni`](crate::kepler::UniversalKeplerParams::prelim_kepuni) – Initial guesses for ψ.
//! * [`s_funct`](crate::kepler::s_funct) – Stumpff‑like functions and invariants.
//! * [`eccentricity_control`](crate::orb_elem::eccentricity_control) – Eccentricity and energy checks.
//! * [`GAUSS_GRAV`](crate::constants::GAUSS_GRAV), [`DPI`](crate::constants::DPI) – Constants.

mod angles;
mod brent_dekker_solver;
mod newton_solver;
mod orbit_type;
mod params;
mod prelim_kepler;
mod propagation;
mod stumpff;
mod universal_kepler_solution;
mod velocity;

pub use angles::{angle_diff, principal_angle};
pub use newton_solver::{solve_kepuni, solve_kepuni_with_guess};
pub use orbit_type::OrbitType;
pub use params::{SolverKind, SolverParams, SolverType, UniversalKeplerParams};
pub use prelim_kepler::prelim_elliptic::prelim_elliptic;
pub use prelim_kepler::prelim_hyperbolic::prelim_hyperbolic;
pub use prelim_kepler::prelim_parabolic::ParabolicPrelimMethod;
pub use propagation::{propagate_universal, UniversalPropagResult};
pub use stumpff::s_funct;
pub use universal_kepler_solution::UniversalKeplerSolution;
pub use velocity::{velocity_correction, velocity_correction_with_guess};
