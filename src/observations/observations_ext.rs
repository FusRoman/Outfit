use itertools::Itertools;
use nalgebra::Vector3;
use std::collections::VecDeque;

use crate::{
    constants::{Observations, Radian},
    equinoctial_element::EquinoctialElements,
    error_models::ErrorModel,
    initial_orbit_determination::{gauss::GaussObs, gauss_result::GaussResult},
    keplerian_element::KeplerianElements,
    observations::observations::Observation,
    outfit::Outfit,
    outfit_errors::OutfitError,
};

/// Calculate the weight of a triplet of observations based on the time difference.
///
/// Arguments
/// ---------
/// * `time`: A reference to a `Vector3<f64>` containing the time of the three observations.
/// * `dtw`: A `f64` representing an optimal interval time requested between the observations of the triplet.
///
/// Return
/// ------
/// * A `f64` representing the weight of the triplet.
fn triplet_weight(time1: f64, time2: f64, time3: f64, dtw: f64) -> f64 {
    fn s3dtw(dt: f64, dtw: f64) -> f64 {
        if dt <= dtw {
            dtw / dt
        } else {
            1.0 + dt / dtw
        }
    }

    let dt1 = time2 - time1;
    let dt2 = time3 - time2;

    s3dtw(dt1, dtw) + s3dtw(dt2, dtw)
}

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
/// ```rust
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
/// * [`KeplerianElements`] – Orbital elements used to propagate orbits.
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
trait ObservationsExt {
    /// Compute triplets of observations.
    ///
    /// This function computes triplets of observations from the given set of observations.
    /// It filters the observations based on the time difference between the first and last observation in the triplet,
    /// and sorts the triplets based on their weights.
    ///
    /// Arguments
    /// ---------
    /// * `dt_min`: An optional minimum time difference between the first and last observation in the triplet. (default: 0.03)
    /// * `dt_max`: An optional maximum time difference between the first and last observation in the triplet. (default: 150.0)
    /// * `optimal_interval_time`: An optional optimal interval time requested between the observations of the triplet. (default: 20.0)
    /// * `max_triplet`: An optional maximum number of triplets to return. (default: 10)
    ///
    /// Return
    /// ------
    /// * A vector of `GaussObs` representing the computed triplets of observations.
    fn compute_triplets(
        &mut self,
        state: &Outfit,
        dt_min: Option<f64>,
        dt_max: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplet: Option<u32>,
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

    /// Compute the RMS of normalized astrometric residuals over a selected interval of observations.
    ///
    /// This function evaluates the quality of a preliminary orbit by computing the root-mean-square (RMS)
    /// of normalized squared astrometric residuals (in RA and DEC) across a subset of observations selected
    /// from a `GaussObs` triplet. The subset is determined using a time interval expansion (`extf`) and a
    /// maximum allowed duration (`dtmax`), centered on the triplet.
    ///
    /// # Arguments
    /// ---------------
    /// * `state` - The current global state providing ephemerides, Earth orientation data, and time conversion.
    /// * `triplets` - A set of three observations used to generate the initial orbit (Gauss method).
    /// * `orbit_element` - The keplerian orbital elements to be tested against the observation arc.
    /// * `extf` - A fractional extension factor used to expand the time interval around the triplet center.
    /// * `dtmax` - A maximum time window (in days) allowed for including additional observations.
    ///
    /// # Returns
    /// ----------
    /// * `Result<f64, OutfitError>` - The RMS (root-mean-square) of the normalized residuals over the selected interval, in radians.
    ///
    /// # Computation Details
    /// ----------
    /// - The residuals per observation are computed using [`Observation::ephemeris_error`](crate::observations::observations::Observation::ephemeris_error), returning a weighted sum of squares (normalized).
    /// - The RMS is computed as:
    ///   RMS = √[ (1 / 2N) × Σᵢ (RAᵢ² + DECᵢ²) ]
    ///   where N is the number of observations in the selected interval.
    ///
    /// # Errors
    /// ----------
    /// Returns `OutfitError` if:
    /// - The interval selection fails (e.g., no valid observations in time range),
    /// - Orbit propagation or ephemeris lookup fails for any observation.
    ///
    /// # Units
    /// ----------
    /// - The final RMS is expressed in **radians**.
    fn rms_orbit_error(
        &self,
        state: &Outfit,
        triplets: &GaussObs,
        orbit_element: &KeplerianElements,
        extf: f64,
        dtmax: f64,
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
/// 1. Building candidate triplets of observations (see [`compute_triplets`]),
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
/// ```rust,ignore
/// let (best_orbit, rms) = observations
///     .estimate_best_orbit(
///         &state,
///         &error_model,
///         &mut rng,
///         50,        // noise realizations
///         1.0,       // noise scaling
///         1.5,       // extf factor
///         30.0,      // dtmax (days)
///         Some(0.03),
///         Some(150.0),
///         Some(20.0),
///         Some(10),
///         0.5,
///     )?;
/// if let Some(orbit) = best_orbit {
///     println!("Best preliminary orbit RMS = {rms}");
/// }
/// ```
///
/// ## Algorithmic steps
///
/// 1. **Batch uncertainty correction:**  
///    Observations are first passed through [`apply_batch_rms_correction`].
///
/// 2. **Triplet generation:**  
///    Valid combinations of three observations are extracted with [`compute_triplets`],
///    optionally constrained by `dt_min`, `dt_max_triplet`, and `optimal_interval_time`.
///
/// 3. **Orbit estimation:**  
///    For each triplet, noisy versions are generated (Monte Carlo) and processed by the
///    Gauss method to obtain preliminary orbital elements.
///
/// 4. **Orbit evaluation:**  
///    The preliminary orbits are tested on the full observation arc using [`rms_orbit_error`].
///    The orbit with the smallest RMS is returned.
///
/// ## Performance considerations
///
/// * Typically applied to **small arcs of tens of observations**.
/// * The number of tested triplets can be limited with `max_triplets`.
/// * The noise realization loop dominates runtime (Monte Carlo approach).
///
/// ## Returns
///
/// * `Ok((Some(best_orbit), best_rms))` – if a valid orbit was found,
/// * `Ok((None, f64::MAX))` – if no valid orbit could be estimated,
/// * `Err(OutfitError)` – if an error occurs during propagation or fitting.
///
/// ## See also
///
/// * [`compute_triplets`] – Generates candidate triplets from the observation set.
/// * [`rms_orbit_error`] – Evaluates the quality of an orbit over the full arc.
/// * [`GaussResult`] – Data structure holding the result of a single Gauss method run.
/// * [`KeplerianElements`] – Orbital elements returned by successful preliminary orbit estimation.
pub trait ObservationIOD {
    /// Estimate the best-fitting preliminary orbit from a full set of astrometric observations.
    ///
    /// This method evaluates all viable observation triplets over the entire observation set to
    /// find the best preliminary orbit using the Gauss method. Each triplet is perturbed multiple
    /// times using Gaussian noise to simulate observational uncertainty. The orbit with the smallest
    /// RMS astrometric residual is selected.
    ///
    /// # Arguments
    /// * `state` – Global Outfit state, providing access to ephemerides and time conversion.
    /// * `error_model` – The astrometric error model (typically per-band or per-observatory).
    /// * `rng` – A random number generator used for noise injection.
    /// * `n_noise_realizations` – Number of noisy variants to generate per triplet.
    /// * `noise_scale` – Scaling factor applied to the nominal RA/DEC uncertainties.
    /// * `extf` – Extrapolation fudge factor used in RMS error weighting.
    /// * `dtmax` – Maximum time interval over which to evaluate the RMS.
    /// * `dt_min` – Minimum time interval allowed between triplet endpoints (optional).
    /// * `dt_max_triplet` – Maximum allowed time span for any triplet (optional).
    /// * `optimal_interval_time` – Target time spacing between observations in triplets (optional).
    /// * `max_triplets` – Maximum number of triplets to test (optional).
    /// * `gap_max` – Maximum allowed time gap within a triplet (in days).
    ///
    /// # Returns
    /// * `Ok((Some(best_orbit), best_rms))` – The preliminary orbit and associated RMS if successful.
    /// * `Ok((None, f64::MAX))` – No valid orbit could be estimated.
    /// * `Err(e)` – An error occurred during RMS computation or orbit estimation.
    ///
    /// # Notes
    /// - RMS values are computed using [`rms_orbit_error`] which includes ephemeris propagation and light-time correction.
    /// - Each triplet may yield several preliminary orbits due to noise realizations.
    /// - Observations are preprocessed using [`apply_batch_rms_correction`] to calibrate uncertainties.
    ///
    /// # See also
    /// * [`compute_triplets`] – Extracts valid triplet combinations from the observation set.
    /// * [`generate_noisy_realizations`] – Creates perturbed triplets using RA/DEC uncertainties.
    /// * [`prelim_orbit`] – Estimates a preliminary orbit from a single triplet.
    /// * [`rms_orbit_error`] – Computes the goodness-of-fit for an orbit against all observations.
    fn estimate_best_orbit(
        &mut self,
        state: &Outfit,
        error_model: &ErrorModel,
        rng: &mut impl rand::Rng,
        n_noise_realizations: usize,
        noise_scale: f64,
        extf: f64,
        dtmax: f64,
        dt_min: Option<f64>,
        dt_max_triplet: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplets: Option<u32>,
        gap_max: f64,
    ) -> Result<(Option<GaussResult>, f64), OutfitError>;
}

impl ObservationsExt for Observations {
    fn compute_triplets(
        &mut self,
        state: &Outfit,
        dt_min: Option<f64>,
        dt_max: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplet: Option<u32>,
    ) -> Vec<GaussObs> {
        self.sort_by(|a, b| {
            a.time.partial_cmp(&b.time).expect(
                format!(
                    "Error sorting observations for object number: {:?} vs {:?}",
                    a.time, b.time
                )
                .as_str(),
            )
        });

        self.iter()
            .enumerate()
            .tuple_combinations::<(_, _, _)>()
            .filter_map(|(obs1, obs2, obs3)| {
                let dt13 = obs3.1.time - obs1.1.time;
                if dt13 < dt_min.unwrap_or(0.03) || dt13 > dt_max.unwrap_or(150.) {
                    return None;
                }

                Some((
                    (obs1.0, obs2.0, obs3.0),
                    triplet_weight(
                        obs1.1.time,
                        obs2.1.time,
                        obs3.1.time,
                        optimal_interval_time.unwrap_or(20.0),
                    ),
                ))
            })
            .sorted_by(|(_, w1), (_, w2)| {
                w1.partial_cmp(w2).expect("Error sorting triplet weights")
            })
            .take(max_triplet.unwrap_or(10) as usize)
            .map(|((idx1, idx2, idx3), _)| {
                let obs1 = &self[idx1];
                let obs2 = &self[idx2];
                let obs3 = &self[idx3];

                // Get the observer for each observation
                let observers = [
                    obs1.get_observer(state),
                    obs2.get_observer(state),
                    obs3.get_observer(state),
                ];

                // all coordinates and associated errors are in radians for the gauss preliminary orbit
                GaussObs::with_observer_position(
                    state,
                    Vector3::new(idx1, idx2, idx3),
                    Vector3::new(obs1.ra, obs2.ra, obs3.ra),
                    Vector3::new(obs1.dec, obs2.dec, obs3.dec),
                    Vector3::new(obs1.time, obs2.time, obs3.time),
                    observers,
                )
            })
            .collect::<Vec<GaussObs>>()
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
        orbit_element: &KeplerianElements,
        extf: f64,
        dtmax: f64,
    ) -> Result<f64, OutfitError> {
        let (start_obs_rms, end_obs_rms) = self.select_rms_interval(triplets, extf, dtmax)?;

        let equinoctial_elements: EquinoctialElements = orbit_element.into();

        let rms_all_obs = self[start_obs_rms..=end_obs_rms]
            .iter()
            .map(|obs| obs.ephemeris_error(state, &equinoctial_elements, obs.get_observer(state)))
            .try_fold(0.0, |acc, rms_obs| rms_obs.map(|rms| acc + rms))?;

        Ok((rms_all_obs / (2. * (end_obs_rms - start_obs_rms + 1) as f64)).sqrt())
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
        n_noise_realizations: usize,
        noise_scale: f64,
        extf: f64,
        dtmax: f64,
        dt_min: Option<f64>,
        dt_max_triplet: Option<f64>,
        optimal_interval_time: Option<f64>,
        max_triplets: Option<u32>,
        gap_max: f64,
    ) -> Result<(Option<GaussResult>, f64), OutfitError> {
        println!("\n\n ____ Estimating best orbit from observations...");
        println!("Estimating best orbit from {} observations", self.len());

        println!("Applying batch RMS correction with gap max: {}", gap_max);
        // Apply uncertainty calibration based on RMS statistics from the error model
        self.apply_batch_rms_correction(error_model, gap_max);

        println!("Computing triplets of observations...");

        // Generate candidate triplets (3-observation sets) based on temporal constraints
        let triplets = self.compute_triplets(
            state,
            dt_min,
            dt_max_triplet,
            optimal_interval_time,
            max_triplets,
        );

        println!("Found {} valid triplets\n", triplets.len());

        let mut best_rms = f64::MAX;
        let mut best_orbit = None;

        println!("Start triplet processing for orbit estimation...\n");
        // Iterate over each triplet to attempt orbit estimation
        for triplet in triplets {
            // Extract astrometric uncertainties for each observation in the triplet
            let (error_ra, error_dec) = self.extract_errors(triplet.idx_obs);

            println!(
                "Processing triplet: {:?} with RA errors: {:?}, DEC errors: {:?}",
                triplet.idx_obs, error_ra, error_dec
            );
            // Generate multiple noisy realizations of the triplet to simulate measurement noise
            let realizations = triplet.generate_noisy_realizations(
                &error_ra,
                &error_dec,
                n_noise_realizations,
                noise_scale,
                rng,
            )?;

            println!(
                "Generated {} noisy realizations for triplet: {:?}",
                realizations.len(),
                triplet.idx_obs
            );

            // For each noisy realization, attempt orbit determination
            for realization in realizations {
                match realization.prelim_orbit() {
                    Ok(orbit) => {
                        // Evaluate how well this orbit fits the full observation set
                        match self.rms_orbit_error(state, &realization, &orbit, extf, dtmax) {
                            Ok(rms) => {
                                // Keep orbit if it's the best (lowest RMS) so far
                                if rms < best_rms {
                                    best_rms = rms;
                                    best_orbit = Some(orbit);
                                }
                            }
                            Err(e) => {
                                eprintln!(
                                    "Failed to compute RMS for orbit from triplet {:?}: {}",
                                    triplet, e
                                );
                            }
                        }
                    }
                    Err(e) => {
                        eprintln!(
                            "Failed to compute preliminary orbit from triplet {:?}: {}",
                            triplet, e
                        );
                    }
                }
            }
        }

        println!("End of initial orbit determination ___ \n\n",);

        // Return best orbit found (if any) and corresponding RMS score
        Ok((best_orbit, best_rms))
    }
}

#[cfg(test)]
mod test_obs_ext {
    use approx::assert_relative_eq;
    use camino::Utf8Path;

    use crate::initial_orbit_determination::gauss::gauss_test::{
        assert_gauss_obs_approx_eq, assert_orbit_close,
    };
    use crate::{
        constants::TrajectorySet, error_models::ErrorModel,
        observations::trajectory_ext::TrajectoryExt, outfit::Outfit,
    };

    use super::*;

    #[test]
    fn test_compute_triplets() {
        let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let mut traj_set =
            TrajectorySet::new_from_80col(&mut env_state, &Utf8Path::new("tests/data/2015AB.obs"));

        let triplets = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory")
            .compute_triplets(&env_state, Some(0.03), Some(150.), None, Some(10));

        assert_eq!(
            triplets.len(),
            10,
            "Expected 10 triplets, got {}",
            triplets.len()
        );

        let expected_triplets = GaussObs {
            idx_obs: [[23, 24, 33]].into(),
            ra: [[1.6893715963476699, 1.689861452091063, 1.7527345385664372]].into(),
            dec: [[1.082468037385525, 0.9436790189346231, 0.8273762407899986]].into(),
            time: [[57028.479297592596, 57049.2318575926, 57063.97711759259]].into(),
            observer_position: [
                [-0.2645666171486676, 0.8689351643673471, 0.3766996211112465],
                [-0.5889735526502539, 0.7240117187952059, 0.3138734206791042],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        assert_gauss_obs_approx_eq(&triplets[0], &expected_triplets, 1e-12);

        let expected_triplet = GaussObs {
            idx_obs: [[21, 25, 33]].into(),
            ra: [[1.6894680985108947, 1.6898894500811472, 1.7527345385664372]].into(),
            dec: [[1.0825984522657437, 0.9435805047946215, 0.8273762407899986]].into(),
            time: [[57028.45404759259, 57049.245147592585, 57063.97711759259]].into(),
            observer_position: [
                [-0.26413563361674103, 0.8690466209095019, 0.3767466856686271],
                [-0.5891631852172257, 0.7238872516832191, 0.3138186516545291],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        assert_gauss_obs_approx_eq(&triplets[9], &expected_triplet, 1e-12);
    }

    #[test]
    fn test_select_rms_interval() {
        let mut env_state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let mut traj_set =
            TrajectorySet::new_from_80col(&mut env_state, &Utf8Path::new("tests/data/2015AB.obs"));

        let traj = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory");

        let triplets = traj.compute_triplets(&env_state, Some(0.03), Some(150.), None, Some(10));
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
        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

        let traj = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory");

        traj.apply_batch_rms_correction(&ErrorModel::FCCT14, 8.0 / 24.0);

        let triplets = GaussObs::new(
            Vector3::new(34, 35, 36),
            [[1.7897976233412669, 1.7898659093482510, 1.7899347771316527]].into(),
            [[
                0.77917805235018101,
                0.77908666497129186,
                0.77899653810797365,
            ]]
            .into(),
            [[57070.238017592594, 57070.250007592593, 57070.262067592594]].into(),
        );

        let kepler = KeplerianElements {
            reference_epoch: 57049.242334573748,
            semi_major_axis: 1.8017360713154256,
            eccentricity: 0.28355914566870571,
            inclination: 0.20267383288689386,
            ascending_node_longitude: 7.9559790236937815E-003,
            periapsis_argument: 1.2451951387589135,
            mean_anomaly: 0.44054589015887125,
        };

        let rms = traj
            .rms_orbit_error(&OUTFIT_HORIZON_TEST.0, &triplets, &kepler, -1.0, 30.)
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

            assert_ulps_eq!(traj[2].error_ra, 2.5070595078906952E-006, max_ulps = 2);
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
                },
                Observation {
                    observer: 0,
                    ra: 1.1,
                    dec: 0.6,
                    error_ra: 3e-6,
                    error_dec: 4e-6,
                    time: 59000.1,
                },
                Observation {
                    observer: 0,
                    ra: 1.2,
                    dec: 0.7,
                    error_ra: 5e-6,
                    error_dec: 6e-6,
                    time: 59000.2,
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

        use crate::unit_test_global::OUTFIT_HORIZON_TEST;

        let mut traj_set = OUTFIT_HORIZON_TEST.1.clone();

        let traj = traj_set
            .get_mut(&crate::constants::ObjectNumber::String("K09R05F".into()))
            .expect("Failed to get trajectory");

        let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

        let gap_max = 8.0 / 24.0; // 8 hours in days
        let (best_orbit, best_rms) = traj
            .estimate_best_orbit(
                &OUTFIT_HORIZON_TEST.0,
                &ErrorModel::FCCT14,
                &mut rng,
                5,
                1.0,
                -1.0,
                30.0,
                Some(0.03),
                Some(150.),
                None,
                Some(10),
                gap_max,
            )
            .unwrap();

        dbg!(&best_orbit, &best_rms);

        let expected_orbit = KeplerianElements {
            reference_epoch: 57049.22904488294,
            semi_major_axis: 1.801748431600605,
            eccentricity: 0.283572284127787,
            inclination: 0.20266779609836036,
            ascending_node_longitude: 0.008022659889281067,
            periapsis_argument: 1.245060173584828,
            mean_anomaly: 0.44047943792316746,
        };
        assert_orbit_close(&best_orbit.unwrap(), &expected_orbit, 1e-14);
        assert_relative_eq!(best_rms, 55.14810894219461, epsilon = 1e-14);
    }
}
