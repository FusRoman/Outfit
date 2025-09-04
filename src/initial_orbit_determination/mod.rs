//! # Initial Orbit Determination (IOD) parameters
//!
//! This module defines the [`crate::initial_orbit_determination::IODParams`] configuration
//! struct and its builder, which control how the **Gauss method** selects observation triplets,
//! applies Monte Carlo perturbations, filters candidate solutions, and tunes numerical tolerances.
//!
//! ## Purpose
//!
//! The [`IODParams`](crate::initial_orbit_determination::IODParams) object centralizes all tunable parameters used by
//! [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit).
//! It allows you to:
//!
//! - Select and filter candidate observation triplets (time spans, downsampling, maximum counts),
//! - Apply Monte Carlo perturbations to simulate astrometric noise,
//! - Enforce physical plausibility constraints (eccentricity, perihelion, distance bounds),
//! - Adjust numerical tolerances for Newton–Raphson solvers and root filtering,
//! - Control time windows for RMS evaluation across the observation arc,
//! - **Control Gauss polynomial solving** (Aberth iterations/epsilon, real-root filter tolerance),
//! - **Constrain geometry** via a minimum admissible topocentric distance at the central epoch (`min_rho2_au`),
//! - **Cap solution scanning** with `max_tested_solutions` (how many candidate roots/orbits to keep).
//!
//! ## Pipeline overview
//!
//! 1. **Triplet generation**  
//!    Candidate 3-observation sets are generated, constrained by `dt_min`, `dt_max_triplet`,
//!    and `optimal_interval_time`. Oversized datasets are downsampled to `max_obs_for_triplets`
//!    before triplet selection.
//!
//! 2. **Monte Carlo perturbation**  
//!    Each triplet is expanded into `n_noise_realizations` noisy copies, drawn from Gaussian
//!    distributions scaled by `noise_scale`. This step increases robustness against measurement errors.
//!
//! 3. **Orbit computation & filtering**  
//!    Each (noisy) triplet is processed by the Gauss solver, generating a preliminary orbit.
//!    Orbits are filtered against physical bounds (`max_ecc`, `max_perihelion_au`, `r2_min_au`, `r2_max_au`,
//!    **`min_rho2_au`**) and validated numerically (`newton_eps`, `newton_max_it`, `root_imag_eps`,
//!    **`aberth_max_iter`**, **`aberth_eps`**, **`kepler_eps`**). Up to **`max_tested_solutions`**
//!    admissible solutions can be collected.
//!
//! 4. **RMS evaluation**  
//!    Candidate orbits are scored by computing RMS residuals over a time window around the triplet.
//!    The window is derived from `extf` × triplet span and clamped to at least `dtmax`.
//!    The best-fitting orbit (lowest RMS) is returned.
//!
//! ## Example
//!
//! ```rust,no_run
//! use outfit::initial_orbit_determination::IODParams;
//!
//! let params = IODParams::builder()
//!     .n_noise_realizations(10)
//!     .noise_scale(1.0)
//!     .dt_min(0.5)
//!     .dt_max_triplet(20.0)
//!     .max_ecc(3.0)
//!     .newton_eps(1e-12)
//!     // Optional fine control on Gauss/Aberth & geometry
//!     .min_rho2_au(0.02)
//!     .aberth_max_iter(100)
//!     .aberth_eps(5e-7)
//!     .kepler_eps(1e-12)
//!     .max_tested_solutions(3)
//!     .build()
//!     .unwrap();
//!
//! // Pass the configuration to the Gauss IOD entry point
//! // let (orbit, rms) = observations.estimate_best_orbit(&state, &error_model, &mut rng, &params)?;
//! ```
//!
//! ## See also
//!
//! * [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit) – main IOD entry point
//! * [`crate::initial_orbit_determination::gauss::GaussObs`] – triplet generation & Gauss solver
//! * [`crate::initial_orbit_determination::gauss_result::GaussResult`] – result type for preliminary/corrected orbits
//! * [`crate::initial_orbit_determination::gauss::GaussObs::solve_8poly`] – 8th-degree distance polynomial and root selection
use crate::outfit_errors::OutfitError;
use std::cmp::Ordering::{Equal, Greater, Less};
use std::fmt;

pub mod gauss;
pub mod gauss_result;

/// Configuration parameters controlling the behavior of
/// [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit).
///
/// This structure centralizes all tunable parameters used during
/// preliminary orbit determination with the Gauss method. It lets you adjust:
///
/// * how candidate triplets of observations are selected,
/// * how Monte Carlo perturbations are applied to simulate noise,
/// * how orbits are evaluated, filtered, and numerically refined.
///
/// Overview
/// -----------------
/// The **preliminary orbit determination** pipeline proceeds in stages:
///
/// 1) **Triplet generation** – all valid 3-observation combinations are formed,
///    filtered by `dt_min`, `dt_max_triplet`, and `optimal_interval_time`.
///
/// 2) **Noise perturbation (Monte Carlo)** – for each triplet, `n_noise_realizations`
///    noisy copies are produced using Gaussian perturbations scaled by `noise_scale`.
///
/// 3) **Orbit computation and evaluation** – each (possibly perturbed) triplet is
///    processed by the Gauss solver; the resulting orbit is scored against the
///    full observation arc. The orbit with the lowest RMS is returned.
///
/// This struct controls these steps as well as physical filtering and numerical tolerances.
///
/// Fields
/// -----------------
/// **Triplet generation / Monte Carlo**
/// * `n_noise_realizations` – number of Monte Carlo perturbations per original triplet
///   (higher → more robust, slower).
/// * `noise_scale` – scale factor applied to nominal RA/DEC uncertainties when
///   drawing Gaussian noise (1.0 uses nominal uncertainties).
/// * `extf` – extrapolation factor defining **how wide the RMS evaluation window is**.
///   After a candidate triplet is chosen, the time span covered by the triplet
///   (Δt between first and last obs) is multiplied by this factor:
///
///   ```text
///   dt_window = Δt_triplet × extf
///   ```
///
///   Window semantics:
///   - If `extf >= 0.0`, the window scales with the triplet span (e.g. Δt=5 d, `extf=1.5` → ~7.5 d on each side).
///   - If `extf < 0.0`, a broad fallback is used (e.g. a multiple of the total dataset span).
///     After this calculation, the window is clamped to be at least `dtmax`.
///
/// * `dtmax` – minimum window (days) when evaluating RMS over the observation arc
///   (acts as a floor for the `extf`-derived window).
/// * `dt_min` – minimum allowed span (days) between first and last obs in a triplet.
/// * `dt_max_triplet` – maximum allowed span (days) for any candidate triplet.
/// * `optimal_interval_time` – target spacing (days) within a triplet to favor well-separated obs.
/// * `max_obs_for_triplets` – cap on observation count used to build triplets; if exceeded,
///   the set is uniformly downsampled while keeping first/last.
/// * `max_triplets` – maximum number of triplets evaluated after filtering/sorting
///   (prevents combinatorial blow-up).
/// * `gap_max` – maximum allowed intra-batch time gap (days) for batch RMS corrections.
///
/// **Physical plausibility / filtering**
/// * `max_ecc` – maximum eccentricity accepted for preliminary/corrected solutions.
/// * `max_perihelion_au` – maximum perihelion distance (AU) accepted (filters extreme/far outliers).
/// * `min_rho2_au` – minimum admissible **topocentric** distance at the central epoch (AU)
///   used to reject spurious geometry (e.g., too-close solutions).
/// * `r2_min_au`, `r2_max_au` – plausibility bounds (AU) for the **central heliocentric distance** `r2`
///   used when selecting roots of the 8th-degree polynomial (helps avoid non-physical branches).
///
/// **Numerical tolerances / iterations**
/// * `newton_eps` – absolute tolerance for Newton–Raphson inner solves.
/// * `newton_max_it` – maximum iterations for Newton–Raphson inner solves.
/// * `root_imag_eps` – maximum allowed imaginary part magnitude when promoting a nearly-real
///   complex root to real (AU scale).
/// * `aberth_max_iter` – maximum iterations for the Aberth–Ehrlich polynomial solver.
/// * `aberth_eps` – convergence tolerance for the Aberth root finder.
/// * `kepler_eps` – tolerance used by the universal Kepler solver during velocity correction.
/// * `max_tested_solutions` – maximum number of admissible Gauss solutions to collect
///   (e.g., cap at 3 for stack-only smallvec paths).
///
/// Defaults
/// -----------------
/// The [`Default`] implementation provides a solid starting point:
///
/// ```rust,no_run
/// use outfit::initial_orbit_determination::IODParams;
/// let params = IODParams::default();
/// ```
///
/// Default values:
///
/// * `n_noise_realizations`: 20  
/// * `noise_scale`: 1.0  
/// * `extf`: -1.0 (use broad window)  
/// * `dtmax`: 30.0 d  
/// * `dt_min`: 0.03 d (~43 min)  
/// * `dt_max_triplet`: 150.0 d  
/// * `optimal_interval_time`: 20.0 d  
/// * `max_obs_for_triplets`: 100  
/// * `max_triplets`: 10  
/// * `gap_max`: 8.0 / 24.0 d (8 h)  
///
/// * `max_ecc`: 5.0  
/// * `max_perihelion_au`: 1.0e3 AU  
/// * `min_rho2_au`: 0.01 AU  
/// * `r2_min_au`: 0.05 AU  
/// * `r2_max_au`: 200.0 AU  
///
/// * `newton_eps`: 1.0e-10  
/// * `newton_max_it`: 50  
/// * `root_imag_eps`: 1.0e-6  
/// * `aberth_max_iter`: 50  
/// * `aberth_eps`: 1.0e-6  
/// * `kepler_eps`: 1e3 × `f64::EPSILON`  
/// * `max_tested_solutions`: 3  
///
/// Notes & Validation
/// -----------------
/// * Time parameters must be non-negative; `extf < 0.0` activates the broad fallback window.
/// * `0 < r2_min_au ≤ r2_max_au`, `max_ecc ≥ 0`, `max_perihelion_au ∈ (0, 1.0e3]`.
/// * `min_rho2_au > 0`, `aberth_max_iter ≥ 1`, `aberth_eps > 0`, `kepler_eps > 0`.
/// * `newton_eps > 0`, `newton_max_it ≥ 1`, `root_imag_eps ≥ 0`.
/// * `max_tested_solutions ≥ 1` (consider `≤ 3` when using stack-only smallvec paths).
///
/// See also
/// -----------------
/// * [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit) – main IOD entry point.
/// * [`crate::initial_orbit_determination::gauss::GaussObs`] – triplet generation & Gauss pipeline.
/// * [`crate::initial_orbit_determination::gauss_result::GaussResult`] – result types for preliminary/corrected solutions.
/// * [`crate::initial_orbit_determination::gauss::GaussObs::solve_8poly`] – 8th-degree distance polynomial and root selection.
#[derive(Debug, Clone)]
pub struct IODParams {
    // --- Triplet generation / Monte Carlo ---
    pub n_noise_realizations: usize,
    pub noise_scale: f64,
    pub extf: f64,
    pub dtmax: f64,
    pub dt_min: f64,
    pub dt_max_triplet: f64,
    pub optimal_interval_time: f64,
    pub max_obs_for_triplets: usize,
    pub max_triplets: u32,
    pub gap_max: f64,

    // --- Physical plausibility / filtering ---
    /// Maximum eccentricity accepted for preliminary/corrected solutions.
    pub max_ecc: f64,
    /// Maximum perihelion distance (AU) accepted for solutions (filters hyperbolic/far outliers).
    pub max_perihelion_au: f64,
    /// Minimum admissible topocentric distance at central epoch (AU).
    pub min_rho2_au: f64,
    // --- Gauss polynomial / solver controls ---
    /// Maximum iterations for the Aberth–Ehrlich polynomial solver.
    pub aberth_max_iter: u32,
    /// Convergence tolerance for the Aberth root finder.
    pub aberth_eps: f64,
    /// Tolerance used by the universal Kepler solver in velocity correction.
    pub kepler_eps: f64,
    /// Maximum number of admissible Gauss solutions to collect.
    pub max_tested_solutions: usize,

    /// Plausibility bounds for the central heliocentric distance r2 (AU).
    pub r2_min_au: f64,
    pub r2_max_au: f64,

    // --- Numerical tolerances / iterations ---
    /// Newton–Raphson absolute tolerance (non-dimensional) for inner solves.
    pub newton_eps: f64,
    /// Maximum iterations for Newton–Raphson inner solves.
    pub newton_max_it: usize,
    /// Maximum allowed imaginary part magnitude for complex roots promoted to real (AU).
    pub root_imag_eps: f64,
}

impl IODParams {
    /// Construct a new [`IODParams`] with sensible default values.
    ///
    /// This is equivalent to calling [`IODParams::default()`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Create a new [`IODParamsBuilder`] to configure custom parameters
    /// for Initial Orbit Determination (IOD).
    ///
    /// This is a **fluent builder API** for [`IODParams`], allowing you to
    /// override the default parameters step by step before running
    /// [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit).
    ///
    /// # Example
    ///
    /// ```rust,no_run
    /// use rand::{rngs::StdRng, SeedableRng};
    /// use outfit::initial_orbit_determination::IODParams;
    /// use outfit::constants::Observations;
    /// use outfit::observations::observations_ext::ObservationIOD;
    ///
    /// let params = IODParams::builder()
    ///     .n_noise_realizations(100)
    ///     .noise_scale(1.0)
    ///     .dtmax(30.0)
    ///     .max_triplets(50)
    ///     .build().unwrap();
    ///
    /// # let observations: Observations = unimplemented!();
    /// # let state = unimplemented!();
    /// # let error_model = unimplemented!();
    /// # let mut rng = StdRng::seed_from_u64(42);
    /// # let _ = observations.estimate_best_orbit(&state, &error_model, &mut rng, &params);
    /// ```
    ///
    /// # See also
    /// * [`IODParams`] – Holds all configuration parameters for IOD.
    /// * [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit) – Consumes these parameters to perform orbit determination.
    pub fn builder() -> IODParamsBuilder {
        IODParamsBuilder::new()
    }
}

impl Default for IODParams {
    fn default() -> Self {
        IODParams {
            // Triplet / MC
            n_noise_realizations: 20,
            noise_scale: 1.0,
            extf: -1.0,
            dtmax: 30.0,
            dt_min: 0.03,
            dt_max_triplet: 150.0,
            optimal_interval_time: 20.0,
            max_obs_for_triplets: 100,
            max_triplets: 10,
            gap_max: 8.0 / 24.0, // 8 hours in days

            // Physical filters
            max_ecc: 5.0,
            max_perihelion_au: 1.0e3,
            min_rho2_au: 0.01,

            // Gauss polynomial / solver
            aberth_max_iter: 50,
            aberth_eps: 1.0e-6,
            kepler_eps: 1e3 * f64::EPSILON,
            max_tested_solutions: 3,

            // Heliocentric r2 bounds
            r2_min_au: 0.05,
            r2_max_au: 200.0,

            // Numerics
            newton_eps: 1.0e-10,
            newton_max_it: 50,
            root_imag_eps: 1.0e-6,
        }
    }
}

/// Builder for [`IODParams`], with validation.
#[derive(Debug, Clone)]
pub struct IODParamsBuilder {
    params: IODParams,
}

impl Default for IODParamsBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl IODParamsBuilder {
    /// Create a new builder initialized with default values.
    pub fn new() -> Self {
        Self {
            params: IODParams::default(),
        }
    }

    // --- Triplet / MC ---
    pub fn n_noise_realizations(mut self, v: usize) -> Self {
        self.params.n_noise_realizations = v;
        self
    }
    pub fn noise_scale(mut self, v: f64) -> Self {
        self.params.noise_scale = v;
        self
    }
    pub fn extf(mut self, v: f64) -> Self {
        self.params.extf = v;
        self
    }
    pub fn dtmax(mut self, v: f64) -> Self {
        self.params.dtmax = v;
        self
    }
    pub fn dt_min(mut self, v: f64) -> Self {
        self.params.dt_min = v;
        self
    }
    pub fn dt_max_triplet(mut self, v: f64) -> Self {
        self.params.dt_max_triplet = v;
        self
    }
    pub fn optimal_interval_time(mut self, v: f64) -> Self {
        self.params.optimal_interval_time = v;
        self
    }
    pub fn max_obs_for_triplets(mut self, v: usize) -> Self {
        self.params.max_obs_for_triplets = v;
        self
    }
    pub fn max_triplets(mut self, v: u32) -> Self {
        self.params.max_triplets = v;
        self
    }
    pub fn gap_max(mut self, v: f64) -> Self {
        self.params.gap_max = v;
        self
    }

    // --- Physical filters ---
    pub fn max_ecc(mut self, v: f64) -> Self {
        self.params.max_ecc = v;
        self
    }
    pub fn max_perihelion_au(mut self, v: f64) -> Self {
        self.params.max_perihelion_au = v;
        self
    }
    pub fn min_rho2_au(mut self, v: f64) -> Self {
        self.params.min_rho2_au = v;
        self
    }
    pub fn r2_min_au(mut self, v: f64) -> Self {
        self.params.r2_min_au = v;
        self
    }
    pub fn r2_max_au(mut self, v: f64) -> Self {
        self.params.r2_max_au = v;
        self
    }

    // --- Gauss polynomial / solver ---
    pub fn aberth_max_iter(mut self, v: u32) -> Self {
        self.params.aberth_max_iter = v;
        self
    }
    pub fn aberth_eps(mut self, v: f64) -> Self {
        self.params.aberth_eps = v;
        self
    }
    pub fn kepler_eps(mut self, v: f64) -> Self {
        self.params.kepler_eps = v;
        self
    }
    pub fn max_tested_solutions(mut self, v: usize) -> Self {
        self.params.max_tested_solutions = v;
        self
    }

    // --- Numerics ---
    pub fn newton_eps(mut self, v: f64) -> Self {
        self.params.newton_eps = v;
        self
    }
    pub fn newton_max_it(mut self, v: usize) -> Self {
        self.params.newton_max_it = v;
        self
    }
    pub fn root_imag_eps(mut self, v: f64) -> Self {
        self.params.root_imag_eps = v;
        self
    }

    // ---- Numeric helpers for PartialOrd (handle NaN as invalid) ----

    /// Return true iff x > 0.0 and comparable (i.e., not NaN).
    #[inline]
    fn gt0(x: f64) -> bool {
        x.partial_cmp(&0.0) == Some(Greater)
    }

    /// Return true iff x >= 0.0 and comparable (i.e., not NaN).
    #[inline]
    fn ge0(x: f64) -> bool {
        matches!(x.partial_cmp(&0.0), Some(Greater) | Some(Equal))
    }

    /// Return true iff a <= b and comparable (i.e., not NaN).
    #[inline]
    fn le(a: f64, b: f64) -> bool {
        matches!(a.partial_cmp(&b), Some(Less) | Some(Equal))
    }

    /// Finalize the builder and produce an [`IODParams`] instance.
    ///
    /// This method applies a set of **validation rules** to ensure that the
    /// configured parameters are consistent and safe to use before constructing
    /// the final [`IODParams`] structure.
    ///
    /// Validation rules
    /// -----------------
    /// The following checks are performed:
    ///
    /// * `noise_scale >= 0.0` – Gaussian noise scaling factor cannot be negative.
    /// * `dt_min >= 0.0`, `dt_max_triplet >= 0.0`, `dtmax >= 0.0` – time spans must be non-negative.
    /// * `max_ecc >= 0.0` – eccentricity bound must be non-negative.
    /// * `0.0 < r2_min_au ≤ r2_max_au` – plausibility bounds for the central distance must be ordered and positive.
    /// * `0.0 < max_perihelion_au` – perihelion distance must be positive.
    /// * `min_rho2_au > 0.0` – topocentric floor must be strictly positive.
    /// * `aberth_max_iter ≥ 1`, `aberth_eps > 0.0` – Aberth solver parameters must be positive.
    /// * `kepler_eps > 0.0` – universal Kepler tolerance must be strictly positive.
    /// * `newton_eps > 0.0`, `newton_max_it ≥ 1`, `root_imag_eps ≥ 0.0`.
    /// * `max_tested_solutions ≥ 1`.
    ///
    /// Notes
    /// -----------------
    /// * If `extf < 0.0`, the algorithm ignores `extf` and uses a broad fallback window
    ///   for RMS evaluation (the final window is still clamped to be at least `dtmax`).
    ///
    /// Special cases
    /// -----------------
    /// **`n_noise_realizations = 0`**  
    /// * Allowed. In this case, the noise generation stage does not produce any
    ///   perturbed copies of a triplet. The orbit estimation uses **only the original triplet**
    ///   (so [`generate_noisy_realizations`](crate::initial_orbit_determination::gauss::GaussObs::generate_noisy_realizations)
    ///   returns a vector of size 1).
    ///
    /// **`max_obs_for_triplets < 3`**  
    /// * Accepted. The downsampling function (`downsample_uniform_with_edges_indices`)
    ///   still returns three indices (first, middle, last), ensuring at least one valid triplet.
    ///   Concretely:
    ///   - `max_obs_for_triplets = 1` or `2` behaves like `3`;
    ///   - if `max_obs_for_triplets ≥ n` (number of observations), no downsampling occurs.
    ///
    /// Returns
    /// -----------------
    /// * `Ok(IODParams)` if all values are valid and the parameters object can safely be used.
    /// * `Err(OutfitError)` if any validation rule fails.
    ///
    /// Examples
    /// -----------------
    /// ```rust,no_run
    /// use outfit::initial_orbit_determination::IODParams;
    ///
    /// // Use builder to customize parameters
    /// let params = IODParams::builder()
    ///     .n_noise_realizations(0) // disable noise clones
    ///     .max_obs_for_triplets(2) // behaves like 3: first, middle, last
    ///     .extf(2.0)
    ///     .build()
    ///     .unwrap();
    /// ```
    ///
    /// When `.build()` succeeds, the resulting [`IODParams`] can be passed
    /// to [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit).
    pub fn build(self) -> Result<IODParams, OutfitError> {
        let p = &self.params;

        // --- Basic non-negativity checks (accept zero) ---
        if !Self::ge0(p.noise_scale) {
            return Err(OutfitError::InvalidIODParameter(
                "noise_scale must be non-negative".into(),
            ));
        }
        if !Self::ge0(p.dt_min) || !Self::ge0(p.dt_max_triplet) || !Self::ge0(p.dtmax) {
            return Err(OutfitError::InvalidIODParameter(
                "time parameters must be non-negative".into(),
            ));
        }
        if !Self::ge0(p.max_ecc) {
            return Err(OutfitError::InvalidIODParameter(
                "max_ecc must be >= 0".into(),
            ));
        }
        if !Self::ge0(p.root_imag_eps) {
            return Err(OutfitError::InvalidIODParameter(
                "root_imag_eps must be >= 0".into(),
            ));
        }

        // --- Strictly positive checks (> 0) ---
        if !Self::gt0(p.max_perihelion_au) {
            return Err(OutfitError::InvalidIODParameter(
                "max_perihelion_au must be > 0".into(),
            ));
        }
        if !Self::gt0(p.min_rho2_au) {
            return Err(OutfitError::InvalidIODParameter(
                "min_rho2_au must be > 0".into(),
            ));
        }
        if !Self::gt0(p.aberth_eps) {
            return Err(OutfitError::InvalidIODParameter(
                "aberth_eps must be > 0".into(),
            ));
        }
        if !Self::gt0(p.kepler_eps) {
            return Err(OutfitError::InvalidIODParameter(
                "kepler_eps must be > 0".into(),
            ));
        }
        if !Self::gt0(p.newton_eps) {
            return Err(OutfitError::InvalidIODParameter(
                "newton_eps must be > 0".into(),
            ));
        }

        // --- Iteration counts (>= 1) ---
        if p.newton_max_it == 0 {
            return Err(OutfitError::InvalidIODParameter(
                "newton_max_it must be >= 1".into(),
            ));
        }
        if p.aberth_max_iter == 0 {
            return Err(OutfitError::InvalidIODParameter(
                "aberth_max_iter must be >= 1".into(),
            ));
        }
        if p.max_tested_solutions < 1 {
            return Err(OutfitError::InvalidIODParameter(
                "max_tested_solutions must be >= 1".into(),
            ));
        }

        // --- r2 bounds: 0 < r2_min_au <= r2_max_au ---
        let ok_r2min = Self::gt0(p.r2_min_au);
        let ok_r2max = Self::gt0(p.r2_max_au);
        let ok_order = Self::le(p.r2_min_au, p.r2_max_au);
        if !(ok_r2min && ok_r2max && ok_order) {
            return Err(OutfitError::InvalidIODParameter(
                "require 0 < r2_min_au <= r2_max_au".into(),
            ));
        }

        Ok(self.params)
    }
}

impl fmt::Display for IODParams {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            const PARAM_COL: usize = 50; // width reserved for "name = value"
            writeln!(f, "Initial Orbit Determination Parameters")?;
            writeln!(f, "-------------------------------------")?;

            macro_rules! line {
                ($fmt:expr, $val:expr, $comment:expr) => {{
                    let s = format!($fmt, $val);
                    let pad = if s.len() < PARAM_COL {
                        " ".repeat(PARAM_COL - s.len())
                    } else {
                        " ".to_string()
                    };
                    writeln!(f, "  {}{}# {}", s, pad, $comment)
                }};
            }

            // --- Triplet generation / Monte Carlo ---
            writeln!(f, "[Triplet generation / Monte Carlo]")?;
            line!(
                "n_noise_realizations = {}",
                self.n_noise_realizations,
                "Number of noisy copies per triplet"
            )?;
            line!(
                "noise_scale          = {:.3}",
                self.noise_scale,
                "Scale factor for RA/DEC uncertainties"
            )?;
            line!(
                "extf                 = {:.3}",
                self.extf,
                "Extrapolation factor for RMS window"
            )?;
            line!(
                "dtmax                = {:.3} d",
                self.dtmax,
                "Minimum RMS evaluation window"
            )?;
            line!(
                "dt_min               = {:.3} d",
                self.dt_min,
                "Minimum span between first/last obs"
            )?;
            line!(
                "dt_max_triplet       = {:.3} d",
                self.dt_max_triplet,
                "Maximum allowed triplet span"
            )?;
            line!(
                "optimal_interval_time= {:.3} d",
                self.optimal_interval_time,
                "Target spacing inside triplet"
            )?;
            line!(
                "max_obs_for_triplets = {}",
                self.max_obs_for_triplets,
                "Cap on obs used to build triplets"
            )?;
            line!(
                "max_triplets         = {}",
                self.max_triplets,
                "Maximum number of triplets evaluated"
            )?;
            line!(
                "gap_max              = {:.3} d",
                self.gap_max,
                "Max intra-batch time gap"
            )?;

            // --- Physical plausibility / filtering ---
            writeln!(f, "\n[Physical plausibility / filtering]")?;
            line!(
                "max_ecc              = {:.3}",
                self.max_ecc,
                "Maximum eccentricity accepted"
            )?;
            line!(
                "max_perihelion_au    = {:.3} AU",
                self.max_perihelion_au,
                "Maximum perihelion distance"
            )?;
            line!(
                "min_rho2_au          = {:.3} AU",
                self.min_rho2_au,
                "Minimum topocentric distance"
            )?;
            line!(
                "r2_min_au            = {:.3} AU",
                self.r2_min_au,
                "Minimum heliocentric distance"
            )?;
            line!(
                "r2_max_au            = {:.3} AU",
                self.r2_max_au,
                "Maximum heliocentric distance"
            )?;

            // --- Numerical tolerances / iterations ---
            writeln!(f, "\n[Numerical tolerances / iterations]")?;
            line!(
                "newton_eps           = {:.1e}",
                self.newton_eps,
                "Tolerance for Newton–Raphson"
            )?;
            line!(
                "newton_max_it        = {}",
                self.newton_max_it,
                "Max Newton–Raphson iterations"
            )?;
            line!(
                "root_imag_eps        = {:.1e}",
                self.root_imag_eps,
                "Max imaginary part for promoted roots"
            )?;
            line!(
                "aberth_max_iter      = {}",
                self.aberth_max_iter,
                "Max iterations for Aberth solver"
            )?;
            line!(
                "aberth_eps           = {:.1e}",
                self.aberth_eps,
                "Convergence tolerance for Aberth solver"
            )?;
            line!(
                "kepler_eps           = {:.1e}",
                self.kepler_eps,
                "Tolerance in Kepler solver"
            )?;
            line!(
                "max_tested_solutions = {}",
                self.max_tested_solutions,
                "Max Gauss solutions kept"
            )?;

            Ok(())
        } else {
            write!(
                f,
                "IODParams(n_noise_realizations={}, noise_scale={:.2}, extf={:.2}, dtmax={:.1}d, dt_min={:.2}d, dt_max_triplet={:.1}d, max_triplets={}, max_ecc={:.2}, perihelion≤{:.1}AU, r2∈[{:.2},{:.1}]AU)",
                self.n_noise_realizations,
                self.noise_scale,
                self.extf,
                self.dtmax,
                self.dt_min,
                self.dt_max_triplet,
                self.max_triplets,
                self.max_ecc,
                self.max_perihelion_au,
                self.r2_min_au,
                self.r2_max_au,
            )
        }
    }
}
