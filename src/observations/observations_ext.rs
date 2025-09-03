//! # Observation extensions for orbit determination
//!
//! This module extends the base [`Observations`] collection with methods
//! essential for **initial orbit determination (IOD)** and orbit quality
//! assessment.  
//! It provides two core traits:
//!
//! ## [`ObservationsExt`]
//!
//! High-level utilities for processing and analyzing sets of [`Observation`]s:
//!
//! - **Triplet generation**: [`compute_triplets`](ObservationsExt::compute_triplets)  
//!   Build optimized triplets of three observations for Gauss IOD.
//!
//! - **RMS evaluation**: [`select_rms_interval`](ObservationsExt::select_rms_interval),  
//!   [`rms_orbit_error`](ObservationsExt::rms_orbit_error)  
//!   Select subsets of observations around a triplet and compute RMS residuals
//!   for candidate orbits.
//!
//! - **Error handling**: [`apply_batch_rms_correction`](ObservationsExt::apply_batch_rms_correction)  
//!   Apply statistical corrections to observational errors based on temporal clustering.
//!
//! - **Uncertainty extraction**: [`extract_errors`](ObservationsExt::extract_errors)  
//!   Retrieve RA/DEC uncertainties for a given triplet.
//!
//! ## [`ObservationIOD`]
//!
//! A high-level trait encapsulating the full **initial orbit determination workflow**:
//!
//! 1. Apply uncertainty calibration with [`ObservationsExt::apply_batch_rms_correction`],  
//! 2. Generate candidate triplets with [`ObservationsExt::compute_triplets`],  
//! 3. Perform Monte Carlo perturbations to simulate astrometric noise,  
//! 4. Run the Gauss method on each triplet,  
//! 5. Select the orbit with the lowest RMS over the full arc using [`ObservationsExt::rms_orbit_error`].  
//!
//! This workflow is designed for **short arcs** (newly discovered asteroids, comets, NEOs),
//! where only a handful of observations are available.
//!
//! ## Typical usage
//!
//! ```rust, no_run
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::initial_orbit_determination::IODParams;
//! use outfit::constants::Observations;
//! use outfit::observations::observations_ext::ObservationIOD;
//!
//! let params = IODParams::builder()
//!     .n_noise_realizations(50)
//!     .noise_scale(1.0)
//!     .dtmax(30.0)
//!     .max_triplets(20)
//!     .build().unwrap();
//!
//! let observations: Observations = unimplemented!(); // Load observations here
//! let state = unimplemented!();                      // Outfit environment
//! let error_model = unimplemented!();                // Astrometric error model
//! let mut rng = StdRng::seed_from_u64(123);
//!
//! let (best_orbit, rms) = observations.estimate_best_orbit(
//!     &state, &error_model, &mut rng, &params
//! ).unwrap();
//!
//! if let Some(orbit) = best_orbit {
//!     println!("Best preliminary orbit RMS = {rms}");
//! }
//! ```
//!
//! ## References
//!
//! * Danby, J.M.A. (1992). *Fundamentals of Celestial Mechanics* (2nd ed.).  
//!   Willmann-Bell. Classic reference for Gauss, Laplace, and related IOD methods.
//!
//! * Milani, A., & Gronchi, G. F. (2009). *Theory of Orbit Determination*.  
//!   Cambridge University Press. [doi:10.1017/CBO9781139175371](https://doi.org/10.1017/CBO9781139175371)
//!
//! * Carpino, M., Milani, A., & Chesley, S. R. (2003). *OrbFit: Software for Preliminary Orbit Determination*.  
//!   [https://adams.dm.unipi.it/orbfit/](https://adams.dm.unipi.it/orbfit/)
//!
//! ## See also
//!
//! - [`Observation`] – Representation of a single astrometric measurement.
//! - [`GaussObs`] – Structure encoding a triplet of observations for Gauss IOD.
//! - [`GaussResult`] – Output of the Gauss preliminary orbit solver.
//! - [`IODParams`] – Parameters controlling IOD (triplet constraints, noise realizations).
use itertools::Itertools;
use nalgebra::Vector3;
use std::{collections::VecDeque, ops::ControlFlow};

use crate::{
    constants::{Observations, Radian},
    error_models::ErrorModel,
    initial_orbit_determination::{gauss::GaussObs, gauss_result::GaussResult, IODParams},
    observations::{triplets_iod::generate_triplets, Observation},
    orbit_type::equinoctial_element::EquinoctialElements,
    outfit::Outfit,
    outfit_errors::OutfitError,
};

/// Extension trait for [`Observations`] providing high-level operations
/// commonly used in orbit determination workflows.
///
/// This trait adds methods to process and analyze a collection of
/// astrometric [`Observation`] objects, including:
///
/// * Selection of observation triplets optimized for initial orbit determination (IOD).
/// * Selection of subsets of observations for root-mean-square (RMS) error calculation.
/// * Computation of orbit quality metrics (RMS of astrometric residuals).
/// * Statistical corrections of observational errors based on temporal clustering.
/// * Extraction of astrometric uncertainties for given observation indices.
///
/// # Typical usage
///
/// This trait is intended to be implemented on:
///
/// ```rust, ignore
/// pub type Observations = SmallVec<[Observation; 6]>;
/// ```
///
/// It provides functionality essential for the Gauss method and related
/// algorithms used in preliminary orbit determination.
///
/// # Provided methods
/// - [`compute_triplets`](ObservationsExt::compute_triplets):
///   Build time-filtered triplets of observations, sorted by weight.
/// - [`select_rms_interval`](ObservationsExt::select_rms_interval):
///   Given a triplet, determine which observations lie in an expanded time window.
/// - [`rms_orbit_error`](ObservationsExt::rms_orbit_error):
///   Evaluate the fit of an orbit by computing RMS residuals.
/// - [`apply_batch_rms_correction`](ObservationsExt::apply_batch_rms_correction):
///   Apply correlation-based scaling factors to astrometric errors based on temporal clustering.
/// - [`extract_errors`](ObservationsExt::extract_errors):
///   Retrieve the RA/DEC uncertainties for a given triplet of observations.
///
/// # See also
/// * [`GaussObs`] – Data structure used to represent a triplet of observations.
/// * [`Outfit`] – Global state providing ephemerides and reference frames.
/// * [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements) – Orbital elements used to propagate orbits.
/// * [`Observation`] – Individual observation data structure.
///
/// # References
///
/// * Danby, J.M.A. (1992). *Fundamentals of Celestial Mechanics* (2nd ed.).
///   Willmann-Bell.  
///   Classic reference for preliminary orbit determination methods
///   (Gauss, Laplace, Vaisala) and iterative improvement techniques.
///
/// * Milani, A., & Gronchi, G. F. (2009). *Theory of Orbit Determination*.
///   Cambridge University Press.  
///   Comprehensive modern treatment of statistical orbit determination,
///   including least-squares methods, weighting, and correlation handling.
///   [https://doi.org/10.1017/CBO9781139175371](https://doi.org/10.1017/CBO9781139175371)
///
/// * Carpino, M., Milani, A., & Chesley, S. R. (2003). *OrbFit: Software for Preliminary Orbit Determination*.  
///   Technical documentation distributed with the OrbFit package:
///   [https://adams.dm.unipi.it/orbfit/](https://adams.dm.unipi.it/orbfit/)
///
/// These references describe the algorithms used for:
/// - Building and filtering observation triplets (Gauss method)
/// - Propagating trial orbits and refining them via differential corrections
/// - Handling of astrometric uncertainties and RMS weighting.
///
/// This trait is central to orbit determination pipelines and is designed
/// to work with small batches of observations (often < 100 per object).
pub trait ObservationsExt {
    /// Compute **time-feasible, best-K** triplets of observations for Gauss IOD,
    /// leveraging a lazy **index stream** and a bounded **max-heap** on spacing weight.
    ///
    /// Overview
    /// -----------------
    /// This method is a convenience wrapper around [`generate_triplets`]. It operates
    /// directly on `self` (the current observation set) and returns up to `max_triplet`
    /// **best-scored** candidates for the Gauss preliminary solution. Internally it:
    ///
    /// 1) Uses a `TripletIndexGenerator` that:
    ///    - sorts epochs in place,
    ///    - downsamples to at most `max_obs_for_triplets` (uniform with edges),
    ///    - lazily **streams reduced indices** `(first, middle, last)` constrained by:
    ///      `dt_min ≤ t[last] − t[first] ≤ dt_max`.
    /// 2) Scores each feasible triplet with [`triplet_weight`](crate::observations::triplets_iod::triplet_weight) against `optimal_interval_time`.
    /// 3) Keeps only the **K** smallest weights in a bounded **max-heap** (best-K selection).
    /// 4) Materializes the survivors as [`GaussObs`] by (re)borrowing `self` immutably.
    ///
    /// Compared to brute-force `O(n³)`, the time-windowed enumeration drives the effective
    /// cost toward ~`O(n²)` in typical time distributions, plus `O(n log K)` for heap updates.
    ///
    /// Arguments
    /// -----------------
    /// * `dt_min` – Minimum allowed timespan `[same units as Observation::time]` between the first and last epoch of a triplet.
    /// * `dt_max` – Maximum allowed timespan between the first and last epoch of a triplet.
    /// * `optimal_interval_time` – Target per-gap spacing (e.g., days) used by [`triplet_weight`](crate::observations::triplets_iod::triplet_weight).
    /// * `max_obs_for_triplets` – Upper bound on observations kept after downsampling (uniform with endpoints).
    /// * `max_triplet` – Number `K` of best-scoring triplets to return.
    ///
    /// Return
    /// ----------
    /// * A `Vec<GaussObs>` of length `≤ max_triplet`, **sorted by increasing weight**
    ///   (best geometric spacing first), ready to be passed to `GaussObs::prelim_orbit`.
    ///
    /// Remarks
    /// -------------
    /// * Sorting is **in-place**; call sites should not rely on original ordering afterward.
    /// * The generator avoids overlapping borrows of `self`; only the final K triplets are materialized.
    /// * For robustness studies, each returned triplet can be expanded with
    ///   `GaussObs::realizations_iter` (lazy Monte-Carlo noise).
    ///
    /// Complexity
    /// -----------------
    /// * Enumeration: ~`O(n²)` (per-anchor time window).
    /// * Selection: `O(n log K)` (bounded max-heap).
    /// * Space: `O(1)` per yielded candidate; only K triplets are allocated at the end.
    ///
    /// See also
    /// ------------
    /// * [`generate_triplets`] – Low-level function performing the selection (index stream + heap + materialization).
    /// * [`TripletIndexGenerator`](crate::observations::triplets_generator::TripletIndexGenerator) – Lazy stream of reduced indices constrained by `(dt_min, dt_max)`.
    /// * [`triplet_weight`](crate::observations::triplets_iod::triplet_weight) – Spacing heuristic around `optimal_interval_time`.
    /// * [`GaussObs::realizations_iter`] – On-the-fly noisy realizations for a given triplet.
    fn compute_triplets(
        &mut self,
        dt_min: f64,
        dt_max: f64,
        optimal_interval_time: f64,
        max_obs_for_triplets: usize,
        max_triplet: u32,
    ) -> Vec<GaussObs>;

    /// Select the interval of observations for RMS calculation.
    ///
    /// This function selects the interval of observations for RMS calculation based on the provided triplet.
    /// It computes the maximum allowed interval and finds the start and end indices of the observations
    /// within that interval.
    ///
    /// Arguments
    /// ---------
    /// * `triplets`: A reference to a `GaussObs` representing the triplet of observations.
    /// * `extf`: A `f64` representing the external factor for the interval calculation.
    /// * `dtmax`: A `f64` representing the maximum allowed interval.
    ///
    /// Return
    /// ------
    /// * A `Result` containing a tuple of start and end indices of the observations within the interval,
    ///   or an `OutfitError` if an error occurs.
    fn select_rms_interval(
        &self,
        triplets: &GaussObs,
        extf: f64,
        dtmax: f64,
    ) -> Result<(usize, usize), OutfitError>;

    /// Evaluate the orbit quality by computing the RMS of normalized astrometric residuals
    /// over a time window centered on a Gauss triplet.
    ///
    /// Scientific context
    /// -------------------
    /// This function measures how well a preliminary orbit reproduces the observed
    /// astrometry (RA, DEC). It computes the **root-mean-square (RMS)** of the
    /// normalized residuals between predicted and observed positions, aggregated over
    /// a set of observations surrounding a Gauss triplet.
    ///
    /// Interval selection
    /// -------------------
    /// The observation arc is defined by:
    /// * `extf` – fractional extension factor applied around the triplet center,
    /// * `dtmax` – absolute maximum time span (days) allowed for the arc.
    ///
    /// The effective interval is determined by
    /// [`select_rms_interval`](Self::select_rms_interval), which returns the first
    /// and last indices of the observations to include.
    ///
    /// Computation
    /// ------------
    /// * Each observation contributes a squared normalized residual
    ///   from [`Observation::ephemeris_error`](crate::observations::Observation::ephemeris_error).
    /// * The final RMS is
    ///
    /// ```text
    /// RMS = √[ (1 / (2N)) · Σᵢ (ΔRAᵢ² + ΔDECᵢ²) ]
    /// ```
    ///
    /// where `N` is the number of observations in the selected interval.
    ///
    /// Pruning mode
    /// ------------
    /// If `prune_if_rms_ge` is set:
    /// * The summation stops early once the partial RMS reaches the threshold,
    ///   returning the pruning value directly.
    /// * If `prune_if_rms_ge = ∞`, no early exit occurs (equivalent to no pruning).
    ///
    /// Arguments
    /// ----------
    /// * `state` – Global context providing ephemerides, Earth orientation, and time conversion.
    /// * `triplets` – The Gauss triplet that defined the preliminary orbit.
    /// * `orbit_element` – The orbit (in equinoctial elements) to be tested against the arc.
    /// * `extf` – Fractional time extension of the interval around the triplet.
    /// * `dtmax` – Maximum arc duration (days).
    /// * `prune_if_rms_ge` – Optional RMS cutoff for early termination (see *Pruning mode*).
    ///
    /// Return
    /// -------
    /// * `Ok(rms)` – RMS of the normalized astrometric residuals (radians).
    /// * `Err(OutfitError)` – If interval selection fails or propagation/ephemeris lookup fails.
    ///
    /// Units
    /// -------
    /// * The returned RMS is dimensionless but expressed in **radians**.
    fn rms_orbit_error(
        &self,
        state: &Outfit,
        triplets: &GaussObs,
        orbit_element: &EquinoctialElements,
        extf: f64,
        dtmax: f64,
        prune_if_rms_ge: Option<f64>,
    ) -> Result<f64, OutfitError>;

    /// Apply RMS correction based on temporally clustered batches of observations.
    ///
    /// This method adjusts the astrometric uncertainties (`error_ra`, `error_dec`) of each observation
    /// based on the local density of observations in time and observer identity. Observations that are
    /// close in time (within 8 hours) and come from the same observer are grouped into batches, and a
    /// correction factor is applied to reflect statistical correlation or improvement due to redundancy.
    ///
    /// # Arguments
    /// ---------------
    /// * `error_model` - The error model to use when applying the batch correction. Supported values include:
    ///     - `"vfcc17"`: uses a reduced factor `√(n × 0.25)` if the batch has at least 5 observations,
    ///     - any other string: uses the standard `√n` factor.
    /// * `gap_max` - The maximum time gap (in days) to consider observations as part of the same batch.
    ///
    /// # Behavior
    /// ----------
    /// - Observations are grouped by `observer` and sorted in time.
    /// - A batch is formed when consecutive observations from the same observer are spaced by less than 8 hours.
    /// - Each observation in a batch of size `n` receives a correction:
    ///     - `√n` for standard models,
    ///     - `√(n × 0.25)` for `vfcc17` when `n ≥ 5`.
    /// - If `n < 5` with `vfcc17`, it falls back to `√n`.
    /// - Observations with fixed weights (`force_w`) are not affected (not yet implemented in this version).
    ///
    /// # Returns
    /// ----------
    /// * `()` - This function modifies the observations in-place; it does not return a value.
    ///
    /// # Computation Details
    /// ----------
    /// - The time comparison is based on Modified Julian Date (`MJD`), and the batch window is fixed at 8 hours (i.e., `8.0 / 24.0` days).
    /// - The error fields `error_ra` and `error_dec` are both scaled by the same batch correction factor.
    ///
    /// # Units
    /// ----------
    /// - Input and output uncertainties (`error_ra`, `error_dec`) are expressed in **radians**.
    fn apply_batch_rms_correction(&mut self, error_model: &ErrorModel, gap_max: f64);

    /// Extract astrometric uncertainties (RA and DEC) for a set of three observations.
    ///
    /// Given a triplet of observation indices, this function retrieves the corresponding
    /// astrometric errors in right ascension and declination from the observation set.
    ///
    /// # Arguments
    /// ---------------
    /// * `idx_obs` - A vector of three indices referring to the observations used in the triplet.
    ///
    /// # Returns
    /// ---------------
    /// * A tuple of two `Vector3<Radian>`:
    ///   - The first vector contains the RA uncertainties in radians.
    ///   - The second vector contains the DEC uncertainties in radians.
    ///
    /// # Panics
    /// ---------------
    /// This function will panic if any index in `idx_obs` is out of bounds of the observation set.
    fn extract_errors(&self, idx_obs: Vector3<usize>) -> (Vector3<Radian>, Vector3<Radian>);
}

/// Trait for performing Initial Orbit Determination (IOD) on a set of astrometric observations.
///
/// This trait encapsulates high-level algorithms to derive a **preliminary orbit**
/// from a time series of astrometric observations, using the **Gauss method**.
/// It focuses on searching for the best-fitting orbit over all valid triplets of observations.
///
/// ## Purpose
///
/// The main goal of this trait is to automate the process of:
/// 1. Building candidate triplets of observations (see [`compute_triplets`](ObservationsExt::compute_triplets)),
/// 2. Estimating a preliminary orbit for each triplet using the Gauss method,
/// 3. Perturbing triplets with Gaussian noise to simulate observational uncertainties,
/// 4. Selecting the orbit that minimizes the root-mean-square (RMS) of astrometric residuals
///    when compared to the entire set of observations.
///
/// This process is the standard entry point for orbit determination workflows.
/// It is designed for **short observational arcs**, such as those of newly discovered asteroids.
///
/// ## Typical usage
///
/// ```rust, no_run
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
/// let observations: Observations = unimplemented!(); // Your observations here
/// let state = unimplemented!(); // Your state here
/// let error_model = unimplemented!(); // Your error model here
/// let mut rng = StdRng::seed_from_u64(42);
///
/// let (best_orbit, rms) = observations.estimate_best_orbit(
///     &state, &error_model, &mut rng, &params).unwrap();
///
/// if let Some(orbit) = best_orbit {
///     println!("Best preliminary orbit RMS = {rms}");
/// }
/// ```
///
/// ## Algorithmic steps
///
/// 1. **Batch uncertainty correction:**  
///    Observations are first passed through [`ObservationsExt::apply_batch_rms_correction`].
///
/// 2. **Triplet generation:**  
///    Valid combinations of three observations are extracted with [`ObservationsExt::compute_triplets`],
///    using the configuration parameters from [`IODParams`] (e.g., `dt_min`, `dt_max_triplet`,
///    `optimal_interval_time`, `max_obs_for_triplets`, `max_triplets`).
///
/// 3. **Orbit estimation:**  
///    For each triplet, `n_noise_realizations` noisy variants are generated (Monte Carlo)
///    and processed by the Gauss method to obtain preliminary orbital elements.
///
/// 4. **Orbit evaluation:**  
///    Each preliminary orbit is propagated and compared to the full observation arc using
///    [`ObservationsExt::rms_orbit_error`]. The orbit with the smallest RMS is returned.
///
/// ## Performance considerations
///
/// * Typically applied to **short arcs with tens of observations**.
/// * The number of triplets can be limited via `params.max_triplets`.
/// * The Monte Carlo loop (`params.n_noise_realizations`) dominates runtime.
///
/// ## Returns
///
/// * `Ok((Some(best_orbit), best_rms))` – if a valid orbit was found,
/// * `Ok((None, f64::MAX))` – if no valid orbit could be estimated,
/// * `Err(OutfitError)` – if an error occurs during propagation or fitting.
///
/// ## Parameters
///
/// * [`IODParams`] – Controls the noise sampling, temporal constraints on triplets,
///   and the maximum number of triplets to evaluate.
///
/// ## See also
///
/// * [`ObservationsExt::compute_triplets`] – Generates candidate triplets from the observation set.
/// * [`ObservationsExt::rms_orbit_error`] – Evaluates the quality of an orbit over the full arc.
/// * [`GaussResult`] – Data structure holding the result of a single Gauss method run.
/// * [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements) – Orbital elements returned by successful preliminary orbit estimation.
/// * [`IODParams`] – Groups all configuration options for IOD.
pub trait ObservationIOD {
    /// Estimate the best-fitting preliminary orbit from a full set of astrometric observations.
    ///
    /// This method searches for the best preliminary orbit by evaluating a limited number of
    /// observation triplets generated from the dataset. The process includes:
    ///
    /// 1. **Error calibration**:
    ///    Observations are first preprocessed with [`ObservationsExt::apply_batch_rms_correction`] to account for
    ///    temporal clustering and observer-specific error models.
    ///
    /// 2. **Triplet generation**:
    ///    Candidate triplets are generated using [`ObservationsExt::compute_triplets`], which:
    ///      * Sorts observations by time,
    ///      * Optionally downsamples the dataset to at most `params.max_obs_for_triplets` points
    ///        (uniform in time, always keeping the first and last),
    ///      * Filters valid triplets according to `params.dt_min`, `params.dt_max_triplet`,
    ///        and `params.optimal_interval_time`.
    ///
    /// 3. **Monte Carlo noise sampling**:
    ///    For each triplet, `params.n_noise_realizations` perturbed versions are created using
    ///    Gaussian noise scaled by `params.noise_scale` times the nominal astrometric uncertainties.
    ///
    /// 4. **Orbit estimation and selection**:
    ///    For each (possibly perturbed) triplet, a preliminary orbit is computed with the Gauss method.
    ///    The resulting orbit is evaluated over the full set of observations using [`ObservationsExt::rms_orbit_error`].
    ///    The orbit with the smallest RMS is returned.
    ///
    /// # Arguments
    ///
    /// * `state` –
    ///   Global [`Outfit`] state, providing ephemerides and time conversions.
    /// * `error_model` –
    ///   The astrometric error model (typically per-band or per-observatory).
    /// * `rng` –
    ///   A random number generator used to draw Gaussian perturbations.
    /// * `params` –
    ///   Parameters controlling the initial orbit determination, including:
    ///     * `n_noise_realizations`: number of noisy triplet variants generated per original triplet,
    ///     * `noise_scale`: scaling factor for the noise,
    ///     * `extf`: extrapolation factor for RMS evaluation,
    ///     * `dtmax`: maximum time interval for RMS evaluation,
    ///     * `dt_min`, `dt_max_triplet`, `optimal_interval_time`: constraints on triplet spans,
    ///     * `max_obs_for_triplets`: maximum number of observations to keep when building triplets,
    ///     * `max_triplets`: maximum number of triplets to process,
    ///     * `gap_max`: maximum allowed time gap within a batch for RMS corrections.
    ///
    /// # Returns
    ///
    /// * `Ok((Some(best_orbit), best_rms))` – The best preliminary orbit found and its RMS.
    /// * `Ok((None, f64::MAX))` – No valid orbit could be estimated.
    /// * `Err(e)` – An error occurred during orbit estimation or RMS evaluation.
    ///
    /// # Notes
    ///
    /// - RMS values are computed with [`ObservationsExt::rms_orbit_error`], which accounts for
    ///   light-time correction and ephemeris propagation.
    /// - Each triplet can produce several preliminary orbit candidates due to
    ///   noise realizations.
    /// - The `max_obs_for_triplets` parameter is crucial for large datasets,
    ///   as it avoids the combinatorial explosion of triplets.
    ///
    /// # See also
    ///
    /// * [`ObservationsExt::compute_triplets`] – Selects triplets from the observation set.
    /// * [`GaussObs::generate_noisy_realizations`] – Creates perturbed triplets with Gaussian noise.
    /// * [`GaussObs::prelim_orbit`] – Computes a preliminary orbit from a single triplet.
    /// * [`ObservationsExt::rms_orbit_error`] – Measures the goodness-of-fit of an orbit against observations.
    /// * [`IODParams`] – Configuration options for the IOD process.
    fn estimate_best_orbit(
        &mut self,
        state: &Outfit,
        error_model: &ErrorModel,
        rng: &mut impl rand::Rng,
        params: &IODParams,
    ) -> Result<(Option<GaussResult>, f64), OutfitError>;
}

impl ObservationsExt for Observations {
    fn compute_triplets(
        &mut self,
        dt_min: f64,
        dt_max: f64,
        optimal_interval_time: f64,
        max_obs_for_triplets: usize,
        max_triplet: u32,
    ) -> Vec<GaussObs> {
        generate_triplets(
            self,
            dt_min,
            dt_max,
            optimal_interval_time,
            max_obs_for_triplets,
            max_triplet,
        )
    }

    fn select_rms_interval(
        &self,
        triplets: &GaussObs,
        extf: f64,
        dtmax: f64,
    ) -> Result<(usize, usize), OutfitError> {
        let nobs = self.len();

        let idx_obs1 = triplets.idx_obs[0];
        let obs1 = self
            .get(idx_obs1)
            .ok_or(OutfitError::ObservationNotFound(idx_obs1))?;

        let idx_obs3 = triplets.idx_obs[2];
        let obs3 = self
            .get(idx_obs3)
            .ok_or(OutfitError::ObservationNotFound(idx_obs3))?;

        let first_obs = self.first().ok_or(OutfitError::ObservationNotFound(0))?;
        let last_obs = self
            .last()
            .ok_or(OutfitError::ObservationNotFound(nobs - 1))?;
        // Step 1: Compute the maximum allowed interval
        let mut dt = if extf >= 0.0 {
            (obs3.time - obs1.time) * extf
        } else {
            10.0 * (last_obs.time - first_obs.time)
        };

        if dtmax >= 0.0 {
            dt = dt.max(dtmax);
        }

        let mut i_start = 0;

        for i in (0..=idx_obs1).rev() {
            if let Some(obs_i) = self.get(i) {
                if obs1.time - obs_i.time > dt {
                    break;
                }
                i_start = i;
            }
        }

        let mut i_end = nobs - 1;

        for i in idx_obs3..nobs {
            if let Some(obs_i) = self.get(i) {
                if obs_i.time - obs3.time > dt {
                    break;
                }
                i_end = i;
            }
        }

        Ok((i_start, i_end))
    }

    fn rms_orbit_error(
        &self,
        state: &Outfit,
        triplets: &GaussObs,
        orbit_element: &EquinoctialElements,
        extf: f64,
        dtmax: f64,
        prune_if_rms_ge: Option<f64>,
    ) -> Result<f64, OutfitError> {
        // Select the time interval [start_obs_rms, end_obs_rms] over which the RMS
        // error is evaluated. The interval depends on the triplet and on external
        // filtering parameters (extf, dtmax).
        let (start_obs_rms, end_obs_rms) = self.select_rms_interval(triplets, extf, dtmax)?;

        // Number of observations contributing to the RMS
        let n_obs = (end_obs_rms - start_obs_rms + 1) as f64;

        // Denominator of the RMS formula: here weighted by 2.0 for consistency
        // with the convention used elsewhere in the code.
        let denom = 2.0 * n_obs;

        // =========================================================================
        // Case 1: No pruning → behave like the "classical" RMS definition
        // =========================================================================
        if prune_if_rms_ge.is_none() {
            // Accumulate the squared ephemeris errors for each observation
            let sum = self[start_obs_rms..=end_obs_rms]
                .iter()
                .map(|obs| obs.ephemeris_error(state, orbit_element))
                // try_fold propagates errors from ephemeris_error while summing
                .try_fold(0.0, |acc, term| term.map(|v| acc + v))?;

            // Final RMS = sqrt( sum / denom )
            return Ok((sum / denom).sqrt());
        }

        // =========================================================================
        // Case 2: Pruning enabled → early stop if RMS exceeds a threshold
        // =========================================================================
        let prune = prune_if_rms_ge.unwrap();

        // Convert the RMS cutoff into a sum cutoff:
        // RMS² = sum / denom  →  stop if sum ≥ (prune² * denom).
        let sum_cutoff = if prune.is_finite() {
            prune * prune * denom
        } else {
            f64::INFINITY // "no real cutoff" if prune = ∞
        };

        // Iterate over observations and accumulate squared errors.
        // We use ControlFlow to allow early exit:
        //   - Continue(sum): keep summing,
        //   - Break(value):  stop early and return the pruning threshold.
        let folded: ControlFlow<f64, f64> = self[start_obs_rms..=end_obs_rms]
            .iter()
            .map(|obs| obs.ephemeris_error(state, orbit_element))
            .try_fold(0.0, |acc, term| match term {
                Ok(v) => {
                    let new_sum = acc + v;
                    if new_sum >= sum_cutoff {
                        // Early exit: threshold reached, return directly
                        ControlFlow::Break(prune)
                    } else {
                        ControlFlow::Continue(new_sum)
                    }
                }
                // In case of error in ephemeris_error, also exit with pruning value.
                Err(_) => ControlFlow::Break(prune),
            });

        // Final RMS depending on whether we exited early or not
        match folded {
            ControlFlow::Continue(sum) => Ok((sum / denom).sqrt()),
            ControlFlow::Break(rms) => Ok(rms),
        }
    }

    fn apply_batch_rms_correction(&mut self, error_model: &ErrorModel, gap_max: f64) {
        // Step 1: Sort in time
        self.sort_by(|a, b| a.time.partial_cmp(&b.time).unwrap());

        // Step 2: Group by observer
        for (_observer_id, group) in &self.into_iter().chunk_by(|obs| obs.observer) {
            // Step 3: Batch grouping using sliding window
            let mut batch: VecDeque<&mut Observation> = VecDeque::new();
            let mut iter = group.peekable();

            while let Some(obs) = iter.next() {
                batch.push_back(obs);

                // Extend batch while within gap_max
                while let Some(next) = iter.peek() {
                    let dt = next.time
                        - batch
                            .back()
                            .expect("in apply_batch_rms_correction: batch should not be empty")
                            .time;
                    if dt <= gap_max {
                        batch.push_back(iter.next().expect(
                            "in apply_batch_rms_correction: next in batch should not be None",
                        ));
                    } else {
                        break;
                    }
                }

                // Apply correction to current batch
                let n = batch.len();
                if n > 0 {
                    let factor = match error_model {
                        ErrorModel::VFCC17 if n >= 5 => (n as f64 * 0.25).sqrt(),
                        _ => (n as f64).sqrt(),
                    };

                    for obs in batch.drain(..) {
                        obs.error_ra *= factor;
                        obs.error_dec *= factor;
                    }
                }
            }
        }
    }

    fn extract_errors(&self, idx_obs: Vector3<usize>) -> (Vector3<Radian>, Vector3<Radian>) {
        let (errors_ra, errors_dec): (Vec<_>, Vec<_>) = idx_obs
            .into_iter()
            .map(|i| {
                let obs = &self[*i];
                (obs.error_ra, obs.error_dec)
            })
            .unzip();

        (
            Vector3::from_column_slice(&errors_ra),
            Vector3::from_column_slice(&errors_dec),
        )
    }
}

impl ObservationIOD for Observations {
    fn estimate_best_orbit(
        &mut self,
        state: &Outfit,
        error_model: &ErrorModel,
        rng: &mut impl rand::Rng,
        params: &IODParams,
    ) -> Result<(Option<GaussResult>, f64), OutfitError> {
        // --- Stage 1: Calibrate uncertainties once for the whole batch.
        // This aligns quoted per-obs errors with empirical RMS statistics.
        self.apply_batch_rms_correction(error_model, params.gap_max);

        // --- Stage 2: Enumerate candidate triplets under temporal constraints.
        let triplets = self.compute_triplets(
            params.dt_min,
            params.dt_max_triplet,
            params.optimal_interval_time,
            params.max_obs_for_triplets,
            params.max_triplets,
        );

        if triplets.is_empty() {
            let span = if self.is_empty() {
                0.0
            } else {
                self.last().unwrap().time - self.first().unwrap().time
            };
            return Err(OutfitError::NoFeasibleTriplets {
                span,
                n_obs: self.len(),
                dt_min: params.dt_min,
                dt_max: params.dt_max_triplet,
            });
        }

        // Current best (lowest) RMS and its orbit.
        // Using +∞ avoids Option branching in the hot path.
        let mut best_rms = f64::INFINITY;
        let mut best_orbit: Option<GaussResult> = None;

        // Keep the last encountered error so that we can report something meaningful if *all* fail.
        // We don't clone: we keep only the most recent error by moving it in.
        let mut last_error: Option<OutfitError> = None;

        // For diagnostics, count how many realizations we actually attempted.
        let mut n_attempts: usize = 0;

        // --- Stage 3: Explore each triplet.
        for triplet in triplets {
            // Extract 1-σ astrometric uncertainties for the three obs of this triplet.
            let (error_ra, error_dec) = self.extract_errors(triplet.idx_obs);

            // --- Stage 4: For each (lazy) noisy realization of this triplet...
            // The iterator yields the original triplet first, then noisy copies.
            for realization in triplet.realizations_iter(
                &error_ra,
                &error_dec,
                params.n_noise_realizations,
                params.noise_scale,
                rng,
            ) {
                n_attempts += 1;

                // 4.a) Preliminary Gauss solution on the current realization.
                let gauss_res = match realization.prelim_orbit(state, params) {
                    Ok(res) => res,
                    Err(e) => {
                        // Record the failure and continue exploring.
                        last_error = Some(e);
                        continue;
                    }
                };

                // 4.b) Convert to the element set required by the scorer.
                let equinoctial_elements = gauss_res.get_orbit().to_equinoctial()?;

                // 4.c) Score orbit vs. full observation set (RMS residual).
                let rms = match self.rms_orbit_error(
                    state,
                    &realization,
                    &equinoctial_elements,
                    params.extf,
                    params.dtmax,
                    Some(best_rms),
                ) {
                    Ok(v) => v,
                    Err(e) => {
                        last_error = Some(e);
                        continue;
                    }
                };

                // 4.d) Keep the best candidate so far.
                if rms < best_rms {
                    best_rms = rms;
                    best_orbit = Some(gauss_res);
                }
            }
        }

        // --- Stage 5: If at least one candidate succeeded, return the best; otherwise, propagate an error.
        if let Some(orbit) = best_orbit {
            Ok((Some(orbit), best_rms))
        } else {
            // If nothing succeeded, propagate a structured error with the last underlying cause.
            // Fallback to a domain-specific unit error if we never captured any (e.g., no attempts).
            let root_cause = match last_error {
                Some(e) => e,
                None => panic!("In estimate_best_orbit: no error captured but best_orbit is None, this should not happen"),
            };
            Err(OutfitError::NoViableOrbit {
                cause: Box::new(root_cause),
                attempts: n_attempts,
            })
        }
    }
}

#[cfg(test)]
mod test_obs_ext {

    use crate::error_models::ErrorModel;

    use super::*;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_select_rms_interval() {
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

        let traj_number = crate::constants::ObjectNumber::String("K09R05F".into());
        let traj_len = traj_set
            .get(&traj_number)
            .expect("Failed to get trajectory")
            .len();

        let traj = traj_set
            .get_mut(&traj_number)
            .expect("Failed to get trajectory");

        let triplets = traj.compute_triplets(0.03, 150.0, 20.0, traj_len, 10);
        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), -1., 30.)
            .unwrap();

        assert_eq!(u1, 0);
        assert_eq!(u2, 36);

        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), 10., 30.)
            .unwrap();

        assert_eq!(u1, 14);
        assert_eq!(u2, 36);

        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), 0.001, 3.)
            .unwrap();

        assert_eq!(u1, 17);
        assert_eq!(u2, 33);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_rms_trajectory() {
        use nalgebra::Matrix3;

        use crate::{
            orbit_type::keplerian_element::KeplerianElements, unit_test_global::OUTFIT_HORIZON_TEST,
        };

        let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

        let traj = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory");

        traj.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

        let triplets = GaussObs {
            idx_obs: Vector3::new(34, 35, 36),
            ra: [[
                1.789_797_623_341_267,
                1.789_865_909_348_251,
                1.7899347771316527,
            ]]
            .into(),
            dec: [[
                0.779_178_052_350_181,
                0.779_086_664_971_291_9,
                0.778_996_538_107_973_6,
            ]]
            .into(),
            time: [[
                57070.238017592594,
                57_070.250_007_592_59,
                57070.262067592594,
            ]]
            .into(),
            observer_helio_position: Matrix3::zeros(),
        };

        let kepler = KeplerianElements {
            reference_epoch: 57_049.242_334_573_75,
            semi_major_axis: 1.8017360713154256,
            eccentricity: 0.283_559_145_668_705_7,
            inclination: 0.20267383288689386,
            ascending_node_longitude: 7.955_979_023_693_781E-3,
            periapsis_argument: 1.2451951387589135,
            mean_anomaly: 0.44054589015887125,
        };

        let rms = traj
            .rms_orbit_error(
                &OUTFIT_HORIZON_TEST.0,
                &triplets,
                &kepler.into(),
                -1.0,
                30.,
                None,
            )
            .unwrap();

        assert_eq!(rms, 68.88650730830162);
    }

    mod test_batch_rms_correction {
        use crate::constants::MJD;
        use approx::assert_ulps_eq;
        use smallvec::smallvec;

        use super::*;

        fn obs(observer: u16, time: MJD) -> Observation {
            Observation {
                observer,
                ra: 1.0,
                error_ra: 1e-6,
                dec: 0.5,
                error_dec: 2e-6,
                time,
                observer_earth_position: Vector3::zeros(),
                observer_helio_position: Vector3::zeros(),
            }
        }

        #[test]
        fn test_single_batch_vfcc17_large() {
            let base_time = 59000.0;
            let mut obs: Observations = smallvec![
                obs(1, base_time),
                obs(1, base_time + 0.01),
                obs(1, base_time + 0.02),
                obs(1, base_time + 0.03),
                obs(1, base_time + 0.04), // n = 5
            ];

            obs.apply_batch_rms_correction(&ErrorModel::VFCC17, 8.0 / 24.0);

            let factor = (5.0_f64 * 0.25_f64).sqrt();
            for ob in &obs {
                assert_ulps_eq!(ob.error_ra, 1e-6 * factor, max_ulps = 2);
                assert_ulps_eq!(ob.error_dec, 2e-6 * factor, max_ulps = 2);
            }
        }

        #[test]
        fn test_single_batch_small_n() {
            let base_time = 59000.0;
            let mut obs: Observations = smallvec![
                obs(2, base_time),
                obs(2, base_time + 0.01), // n = 2
            ];

            obs.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

            let factor = (2.0f64).sqrt();
            for ob in &obs {
                assert_ulps_eq!(ob.error_ra, 1e-6 * factor, max_ulps = 2);
                assert_ulps_eq!(ob.error_dec, 2e-6 * factor, max_ulps = 2);
            }
        }

        #[test]
        fn test_multiple_batches_same_observer() {
            let base_time = 59000.0;
            let mut obs: Observations = smallvec![
                obs(3, base_time),
                obs(3, base_time + 0.01), // batch 1 (n = 2)
                obs(3, base_time + 1.0),  // isolated, batch 2 (n = 1)
            ];

            obs.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

            let factor1 = (2.0f64).sqrt();
            let factor2 = 1.0;

            assert_ulps_eq!(obs[0].error_ra, 1e-6 * factor1, max_ulps = 2);
            assert_ulps_eq!(obs[1].error_ra, 1e-6 * factor1, max_ulps = 2);
            assert_ulps_eq!(obs[2].error_ra, 1e-6 * factor2, max_ulps = 2);
        }

        #[test]
        fn test_different_observers_are_not_grouped() {
            let base_time = 59000.0;
            let mut obs: Observations = smallvec![
                obs(10, base_time),
                obs(11, base_time + 0.01),
                obs(12, base_time + 0.02),
            ];

            obs.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

            for ob in &obs {
                assert_ulps_eq!(ob.error_ra, 1e-6, max_ulps = 2);
                assert_ulps_eq!(ob.error_dec, 2e-6, max_ulps = 2);
            }
        }

        #[test]
        fn test_batch_gaps_exceed_gapmax() {
            let mut obs: Observations = smallvec![
                obs(5, 59000.0),
                obs(5, 59001.0), // > 8h => separate
            ];

            obs.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

            for ob in &obs {
                assert_ulps_eq!(ob.error_ra, 1e-6, max_ulps = 2);
                assert_ulps_eq!(ob.error_dec, 2e-6, max_ulps = 2);
            }
        }

        #[test]
        #[cfg(feature = "jpl-download")]
        fn test_batch_real_data() {
            use crate::unit_test_global::OUTFIT_HORIZON_TEST;

            let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

            let traj = traj_set
                .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
                .expect("Failed to get trajectory");

            traj.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

            assert_ulps_eq!(traj[0].error_ra, 2.507075226057322e-6, max_ulps = 2);
            assert_ulps_eq!(traj[0].error_dec, 2.036217397086327e-6, max_ulps = 2);

            assert_ulps_eq!(traj[1].error_ra, 2.5070681687218917e-6, max_ulps = 2);
            assert_ulps_eq!(traj[1].error_dec, 2.036217397086327e-6, max_ulps = 2);

            assert_ulps_eq!(traj[2].error_ra, 2.507_059_507_890_695_2E-6, max_ulps = 2);
            assert_ulps_eq!(traj[2].error_dec, 2.036217397086327e-6, max_ulps = 2);
        }
    }

    mod test_extract_errors {
        use super::*;
        use approx::assert_ulps_eq;
        use smallvec::smallvec;

        fn make_observations() -> Observations {
            smallvec![
                Observation {
                    observer: 0,
                    ra: 1.0,
                    dec: 0.5,
                    error_ra: 1e-6,
                    error_dec: 2e-6,
                    time: 59000.0,
                    observer_earth_position: Vector3::zeros(),
                    observer_helio_position: Vector3::zeros(),
                },
                Observation {
                    observer: 0,
                    ra: 1.1,
                    dec: 0.6,
                    error_ra: 3e-6,
                    error_dec: 4e-6,
                    time: 59000.1,
                    observer_earth_position: Vector3::zeros(),
                    observer_helio_position: Vector3::zeros(),
                },
                Observation {
                    observer: 0,
                    ra: 1.2,
                    dec: 0.7,
                    error_ra: 5e-6,
                    error_dec: 6e-6,
                    time: 59000.2,
                    observer_earth_position: Vector3::zeros(),
                    observer_helio_position: Vector3::zeros(),
                },
            ]
        }

        #[test]
        fn test_extract_errors_basic() {
            let obs = make_observations();
            let idx_obs = Vector3::new(0, 1, 2);

            let (ra_errors, dec_errors) = obs.extract_errors(idx_obs);

            assert_ulps_eq!(ra_errors[0], 1e-6, max_ulps = 2);
            assert_ulps_eq!(ra_errors[1], 3e-6, max_ulps = 2);
            assert_ulps_eq!(ra_errors[2], 5e-6, max_ulps = 2);

            assert_ulps_eq!(dec_errors[0], 2e-6, max_ulps = 2);
            assert_ulps_eq!(dec_errors[1], 4e-6, max_ulps = 2);
            assert_ulps_eq!(dec_errors[2], 6e-6, max_ulps = 2);
        }

        #[test]
        #[should_panic(expected = "index out of bounds")]
        fn test_extract_errors_out_of_bounds() {
            let obs = make_observations();
            let idx_obs = Vector3::new(0, 1, 10); // 10 is out of bounds
            let _ = obs.extract_errors(idx_obs);
        }
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_estimate_best_orbit() {
        use approx::assert_relative_eq;
        use rand::{rngs::StdRng, SeedableRng};

        use crate::{
            orbit_type::{
                keplerian_element::KeplerianElements, orbit_type_test::approx_equal,
                OrbitalElements,
            },
            unit_test_global::OUTFIT_HORIZON_TEST,
        };

        let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

        let traj_number = crate::constants::ObjectNumber::String("K09R05F".into());
        let traj_len = traj_set
            .get(&traj_number)
            .expect("Failed to get trajectory")
            .len();

        let traj = traj_set
            .get_mut(&traj_number)
            .expect("Failed to get trajectory");

        let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

        let gap_max = 8.0 / 24.0; // 8 hours in days

        let params = IODParams {
            n_noise_realizations: 5,
            max_obs_for_triplets: traj_len,
            gap_max,
            ..Default::default()
        };

        let (best_orbit, best_rms) = traj
            .estimate_best_orbit(
                &OUTFIT_HORIZON_TEST.0,
                &ErrorModel::FCCT14,
                &mut rng,
                &params,
            )
            .unwrap();

        let binding = best_orbit.unwrap();
        let orbit = binding.get_orbit();

        let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 57049.22904488294,
            semi_major_axis: 1.801748431600605,
            eccentricity: 0.283572284127787,
            inclination: 0.20266779609836036,
            ascending_node_longitude: 0.008022659889281067,
            periapsis_argument: 1.245060173584828,
            mean_anomaly: 0.44047943792316746,
        });

        assert!(approx_equal(orbit, &expected_orbit, 1e-14));
        assert_relative_eq!(best_rms, 55.14810894219461, epsilon = 1e-14);
    }
}
