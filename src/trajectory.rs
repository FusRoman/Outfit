//! Core IOD pipeline trait over sorted observation slices.
//!
//! This module defines [`TrajectoryFit`], a `pub(crate)` trait implemented on
//! `Vec<&Observation>`. It encapsulates the sequential stages of the Gauss
//! Initial Orbit Determination pipeline:
//!
//! 1. **Triplet generation** — [`TrajectoryFit::compute_triplets`] selects the
//!    best-K observation triplets using a time-windowed enumeration and a
//!    bounded max-heap scorer.
//! 2. **RMS interval selection** — [`TrajectoryFit::select_rms_interval`]
//!    determines the observation arc over which the residual quality metric is
//!    evaluated.
//! 3. **RMS evaluation** — [`TrajectoryFit::rms_orbit_error`] measures how
//!    well a candidate orbit reproduces the arc's astrometry.
//! 4. **Orbit selection** — [`TrajectoryFit::estimate_best_orbit`] drives the
//!    full Monte-Carlo noise loop, combining triplets, Gauss solutions, and RMS
//!    scoring to return the best preliminary orbit.
//!
//! The public entry points are [`crate::obs_dataset::FitIOD::fit_iod`] and
//! [`crate::obs_dataset::FitIOD::fit_full_iod`], which call
//! [`TrajectoryFit::estimate_best_orbit`] after building the shared cache.

use std::ops::ControlFlow;

use nalgebra::Vector3;
use photom::{observation_dataset::observation::Observation, Radians};

use crate::{
    cache::OutfitCache,
    initial_orbit_determination::{gauss::GaussObs, triplet_generation::generate_triplets},
    observation_ephemeris::ObservationEphemeris,
    EquinoctialElements, GaussResult, IODParams, JPLEphem, OutfitError, IODRMS,
};

pub(crate) trait TrajectoryFit {
    /// Extract astrometric uncertainties (RA and DEC) for a set of three observations.
    ///
    /// Given a triplet of observation indices, this function retrieves the corresponding
    /// astrometric errors in right ascension and declination from the observation set.
    ///
    /// # Arguments
    ///
    /// - `idx_obs` - A vector of three indices referring to the observations used in the triplet.
    ///
    /// # Returns
    ///
    /// - A tuple of two `Vector3<Radian>`:
    ///   - The first vector contains the RA uncertainties in radians.
    ///   - The second vector contains the DEC uncertainties in radians.
    ///
    /// # Panics
    ///
    /// This function will panic if any index in `idx_obs` is out of bounds of the observation set.
    fn extract_errors(&self, idx_obs: Vector3<usize>) -> (Vector3<Radians>, Vector3<Radians>);

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
    fn compute_triplets(&self, cache: &OutfitCache, params: &IODParams) -> Vec<GaussObs>;

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
        params: &IODParams,
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
        cache: &OutfitCache,
        jpl: &JPLEphem,
        triplets: &GaussObs,
        orbit_element: &EquinoctialElements,
        params: &IODParams,
        prune_if_rms_ge: Option<f64>,
    ) -> Result<f64, OutfitError>;

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
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        params: &IODParams,
        rng: &mut impl rand::Rng,
    ) -> Result<(GaussResult, IODRMS), OutfitError>;
}

impl TrajectoryFit for Vec<&Observation> {
    fn extract_errors(&self, idx_obs: Vector3<usize>) -> (Vector3<Radians>, Vector3<Radians>) {
        let [i, j, k] = [idx_obs[0], idx_obs[1], idx_obs[2]];
        let [c1, c2, c3] = [
            self[i].equ_coord(),
            self[j].equ_coord(),
            self[k].equ_coord(),
        ];
        (
            Vector3::new(c1.ra_error, c2.ra_error, c3.ra_error),
            Vector3::new(c1.dec_error, c2.dec_error, c3.dec_error),
        )
    }

    fn compute_triplets(&self, cache: &OutfitCache, params: &IODParams) -> Vec<GaussObs> {
        generate_triplets(self, cache, params)
    }

    fn select_rms_interval(
        &self,
        triplets: &GaussObs,
        params: &IODParams,
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
        let mut dt = if params.extf >= 0.0 {
            (obs3.mjd_tt() - obs1.mjd_tt()) * params.extf
        } else {
            10.0 * (last_obs.mjd_tt() - first_obs.mjd_tt())
        };

        if params.dtmax >= 0.0 {
            dt = dt.max(params.dtmax);
        }

        let mut i_start = 0;

        for i in (0..=idx_obs1).rev() {
            if let Some(obs_i) = self.get(i) {
                if obs1.mjd_tt() - obs_i.mjd_tt() > dt {
                    break;
                }
                i_start = i;
            }
        }

        let mut i_end = nobs - 1;

        for i in idx_obs3..nobs {
            if let Some(obs_i) = self.get(i) {
                if obs_i.mjd_tt() - obs3.mjd_tt() > dt {
                    break;
                }
                i_end = i;
            }
        }

        Ok((i_start, i_end))
    }

    fn rms_orbit_error(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        triplets: &GaussObs,
        orbit_element: &EquinoctialElements,
        params: &IODParams,
        prune_if_rms_ge: Option<f64>,
    ) -> Result<f64, OutfitError> {
        // Select the time interval [start_obs_rms, end_obs_rms] over which the RMS
        // error is evaluated. The interval depends on the triplet and on external
        // filtering parameters (extf, dtmax).
        let (start_obs_rms, end_obs_rms) = self.select_rms_interval(triplets, params)?;

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
                .map(|obs| obs.ephemeris_error(cache, jpl, orbit_element))
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
            .map(|obs| obs.ephemeris_error(cache, jpl, orbit_element))
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

    fn estimate_best_orbit(
        &self,
        cache: &OutfitCache,
        jpl: &JPLEphem,
        params: &IODParams,
        rng: &mut impl rand::Rng,
    ) -> Result<(GaussResult, IODRMS), OutfitError> {
        // Compute candidate triplets for preliminary orbit estimation, sorted by geometric spacing.
        let triplets = self.compute_triplets(cache, params);

        if triplets.is_empty() {
            let span = if self.is_empty() {
                0.0
            } else {
                self.last().unwrap().mjd_tt() - self.first().unwrap().mjd_tt()
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

        for triplet in triplets {
            // Extract 1-σ astrometric uncertainties for the three obs of this triplet.
            let (error_ra, error_dec) = self.extract_errors(triplet.idx_obs);

            // For each (lazy) noisy realization of this triplet...
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
                let gauss_res = match realization.prelim_orbit(params) {
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
                    cache,
                    jpl,
                    &realization,
                    &equinoctial_elements,
                    params,
                    Some(best_rms),
                ) {
                    Ok(v) => {
                        if !v.is_finite() {
                            last_error = Some(OutfitError::NonFiniteScore(v));
                            continue;
                        } else {
                            v
                        }
                    }
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

        // If at least one candidate succeeded, return the best; otherwise, propagate an error.
        if let Some(orbit) = best_orbit {
            Ok((orbit, best_rms))
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
    use nalgebra::Matrix3;
    use photom::observer::error_model::{ModelCorrection, ObsErrorModel};

    use crate::{
        initial_orbit_determination::IODParamsBuilder,
        orbit_type::orbit_type_test::approx_equal,
        test_fixture::{DATASET_2015AB, JPL_EPHEM_HORIZON, UT1_PROVIDER},
        KeplerianElements, OrbitalElements,
    };

    use approx::assert_relative_eq;
    use rand::{rngs::StdRng, SeedableRng};

    use super::*;

    #[test]
    fn test_select_rms_interval() {
        let corrected_dataset = DATASET_2015AB
            .clone()
            .with_error_model(ObsErrorModel::FCCT14)
            .apply_model_errors();

        let traj = corrected_dataset
            .materialize_trajectory("K09R05F")
            .unwrap()
            .collect_into_vec();
        let traj_len = traj.len();

        let cache =
            OutfitCache::build(&corrected_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

        let iod_params = IODParams {
            dt_min: 0.03,
            dt_max_triplet: 150.0,
            optimal_interval_time: 20.0,
            max_obs_for_triplets: traj_len,
            max_triplets: 10,
            extf: -1.0,
            dtmax: 30.,
            ..Default::default()
        };

        let triplets = traj.compute_triplets(&cache, &iod_params);
        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), &iod_params)
            .unwrap();

        assert_eq!(u1, 0);
        assert_eq!(u2, 36);

        let new_params = IODParamsBuilder::from_params(iod_params)
            .extf(10.)
            .build()
            .unwrap();

        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), &new_params)
            .unwrap();

        assert_eq!(u1, 14);
        assert_eq!(u2, 36);

        let new_params = IODParamsBuilder::from_params(new_params)
            .extf(0.001)
            .dtmax(3.)
            .build()
            .unwrap();

        let (u1, u2) = traj
            .select_rms_interval(triplets.first().unwrap(), &new_params)
            .unwrap();

        assert_eq!(u1, 17);
        assert_eq!(u2, 33);
    }

    #[test]
    fn test_rms_trajectory() {
        let iod_params = IODParams {
            extf: -1.0,
            dtmax: 30.,
            ..Default::default()
        };

        let corrected_dataset = DATASET_2015AB
            .clone()
            .with_error_model(ObsErrorModel::FCCT14)
            .apply_batch_rms_correction(iod_params.gap_max);

        let traj = corrected_dataset
            .materialize_trajectory("K09R05F")
            .unwrap()
            .collect_into_vec();

        let cache =
            OutfitCache::build(&corrected_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

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
                &cache,
                &JPL_EPHEM_HORIZON,
                &triplets,
                &kepler.into(),
                &iod_params,
                None,
            )
            .unwrap();

        assert_eq!(rms, 153.84607281520138);
    }

    mod proptests_extract_errors {
        use super::*;
        use photom::{
            coordinates::equatorial::EquCoord,
            observation_dataset::{observation::ObservationInput, ObsDataset},
            observer::dataset::ObserverId,
            photometry::{Filter, Photometry},
        };
        use proptest::prelude::*;

        fn arb_equ_coord() -> impl Strategy<Value = EquCoord> {
            (
                0.0..(2.0 * std::f64::consts::PI),
                0.0..0.01f64,
                (-std::f64::consts::FRAC_PI_2)..std::f64::consts::FRAC_PI_2,
                0.0..0.01f64,
            )
                .prop_map(|(ra, ra_error, dec, dec_error)| EquCoord {
                    ra,
                    ra_error,
                    dec,
                    dec_error,
                })
        }

        fn make_obs_dataset_with_coords(coords: Vec<EquCoord>) -> ObsDataset {
            let inputs: Vec<ObservationInput> = coords
                .into_iter()
                .enumerate()
                .map(|(i, equ_coord)| {
                    ObservationInput::new(
                        i as u64,
                        equ_coord,
                        Photometry {
                            magnitude: 20.0,
                            error: 0.1,
                            filter: Filter::Int(0),
                        },
                        59000.0,
                        Some(ObserverId::MpcCode(*b"F51")),
                    )
                })
                .collect();

            ObsDataset::default().push_observation(inputs).unwrap().0
        }

        proptest! {
            /// `extract_errors` must return the `ra_error` and `dec_error` of each
            /// indexed observation, in the correct vector positions.
            #[test]
            fn proptest_extract_errors_returns_correct_components(
                coord0 in arb_equ_coord(),
                coord1 in arb_equ_coord(),
                coord2 in arb_equ_coord(),
            ) {
                let dataset = make_obs_dataset_with_coords(vec![coord0, coord1, coord2]);
                let obs: Vec<&Observation> = (0..3)
                    .map(|i| dataset.get_observation(i).unwrap())
                    .collect();

                let idx = Vector3::new(0usize, 1, 2);
                let (ra_errors, dec_errors) = obs.extract_errors(idx);

                prop_assert_eq!(ra_errors[0], coord0.ra_error);
                prop_assert_eq!(ra_errors[1], coord1.ra_error);
                prop_assert_eq!(ra_errors[2], coord2.ra_error);
                prop_assert_eq!(dec_errors[0], coord0.dec_error);
                prop_assert_eq!(dec_errors[1], coord1.dec_error);
                prop_assert_eq!(dec_errors[2], coord2.dec_error);
            }

            /// The index order passed to `extract_errors` must determine which
            /// observation's error ends up at which vector position.
            #[test]
            fn proptest_extract_errors_respects_index_order(
                coord0 in arb_equ_coord(),
                coord1 in arb_equ_coord(),
                coord2 in arb_equ_coord(),
            ) {
                let dataset = make_obs_dataset_with_coords(vec![coord0, coord1, coord2]);
                let obs: Vec<&Observation> = (0..3)
                    .map(|i| dataset.get_observation(i).unwrap())
                    .collect();

                // Permuted index: (2, 0, 1)
                let idx = Vector3::new(2usize, 0, 1);
                let (ra_errors, dec_errors) = obs.extract_errors(idx);

                prop_assert_eq!(ra_errors[0], coord2.ra_error);
                prop_assert_eq!(ra_errors[1], coord0.ra_error);
                prop_assert_eq!(ra_errors[2], coord1.ra_error);
                prop_assert_eq!(dec_errors[0], coord2.dec_error);
                prop_assert_eq!(dec_errors[1], coord0.dec_error);
                prop_assert_eq!(dec_errors[2], coord1.dec_error);
            }
        }
    }

    #[test]
    fn test_estimate_best_orbit() {
        let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

        let iod_params = IODParams::default();

        let corrected_dataset = DATASET_2015AB
            .clone()
            .with_error_model(ObsErrorModel::FCCT14)
            .apply_batch_rms_correction(iod_params.gap_max);

        let traj = corrected_dataset
            .materialize_trajectory("K09R05F")
            .unwrap()
            .collect_into_vec();

        let iod_params = IODParamsBuilder::from_params(iod_params)
            .n_noise_realizations(5)
            .max_obs_for_triplets(traj.len())
            .build()
            .unwrap();

        let cache =
            OutfitCache::build(&corrected_dataset, &JPL_EPHEM_HORIZON, &UT1_PROVIDER).unwrap();

        let (best_orbit, best_rms) = traj
            .estimate_best_orbit(&cache, &JPL_EPHEM_HORIZON, &iod_params, &mut rng)
            .unwrap();

        let binding = best_orbit;
        let orbit = binding.get_orbit();

        let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 57049.22904475403,
            semi_major_axis: 1.8017609974509807,
            eccentricity: 0.2835733667643381,
            inclination: 0.20267686119302475,
            ascending_node_longitude: 0.00799201841873464,
            periapsis_argument: 1.245034216916367,
            mean_anomaly: 0.4405089048961484,
        });

        assert!(approx_equal(orbit, &expected_orbit, 1e-14));
        assert_relative_eq!(best_rms, 222.16583195747745, epsilon = 1e-14);
    }
}
