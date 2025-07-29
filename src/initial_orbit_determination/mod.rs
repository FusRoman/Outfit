use crate::outfit_errors::OutfitError;

pub(crate) mod gauss;
pub mod gauss_result;

/// Configuration parameters controlling the behavior of [`ObservationIOD::estimate_best_orbit`].
///
/// This structure centralizes all the tunable parameters used during
/// preliminary orbit determination with the Gauss method. It allows
/// you to adjust:
///
/// * How candidate triplets of observations are selected,
/// * How Monte Carlo perturbations are applied to simulate noise,
/// * And how orbits are evaluated and filtered.
///
/// # Overview
///
/// The **preliminary orbit determination** process works in several stages:
///
/// 1. **Triplet generation**:
///    All valid combinations of three observations are generated,
///    filtered according to `dt_min`, `dt_max_triplet`, and
///    `optimal_interval_time`.
///
/// 2. **Noise perturbation (Monte Carlo)**:
///    For each triplet, `n_noise_realizations` noisy versions are created
///    using Gaussian perturbations scaled by `noise_scale`.
///
/// 3. **Orbit computation and evaluation**:
///    Each perturbed triplet is processed by the Gauss method, and the
///    resulting orbit is evaluated against the full observation arc.
///    The orbit with the lowest RMS is returned.
///
/// This struct controls these steps.
///
/// # Fields
///
/// * `n_noise_realizations` –
///   Number of Monte Carlo perturbations applied to each original triplet.
///   A higher value increases robustness but also runtime.
///
/// * `noise_scale` –
///   Scaling factor applied to the nominal RA/DEC uncertainties when
///   generating Gaussian noise. A value of `1.0` uses the nominal
///   uncertainties directly.
///
/// * `extf` –
///   Extrapolation factor that controls **how wide the RMS evaluation window is**.
///   
///   After a candidate triplet is chosen, the time span covered by that triplet
///   (Δt between the first and last observation of the triplet) is multiplied by
///   this factor:
///
///   ```text
///   dt_window = Δt_triplet × extf
///   ```
///
///   This window defines how far around the triplet the RMS residuals are computed:
///
///   * If `extf >= 0.0`, the window is proportional to the triplet span.
///     Example: if the triplet covers 5 days and `extf = 1.5`,
///     the RMS evaluation will include ~7.5 days before/after.
///   * If `extf < 0.0`, the algorithm ignores this factor and uses a broad
///     window covering the whole dataset (`10 × total_timespan`).
///
///   After this calculation, the window is clamped to ensure it is at least `dtmax`.
///
///   **Effect:**  
///   - A small value of `extf` means the RMS is computed close to the triplet  
///     → emphasizes local fit quality.  
///   - A large value means the RMS uses a broader time window  
///     → emphasizes global consistency.
///
/// * `dtmax` –
///   Minimum time span (days) used when evaluating an orbit's RMS over the
///   observation arc. This acts as a floor for the window size computed with `extf`.
///
/// * `dt_min` –
///   Minimum allowed time span (days) between the first and last
///   observation in a triplet.
///
/// * `dt_max_triplet` –
///   Maximum allowed time span (days) for the three observations that
///   form a valid triplet.
///
/// * `optimal_interval_time` –
///   Target time spacing (days) between the three observations within a
///   triplet. Used to favor well-separated observations.
///
/// * `max_obs_for_triplets` –
///   Maximum number of observations to keep when generating triplets.
///   If the dataset contains more observations, it is downsampled
///   uniformly while keeping the first and last.
///
/// * `max_triplets` –
///   Maximum number of triplets to evaluate after filtering and sorting.
///   This prevents combinatorial explosion for large datasets.
///
/// * `gap_max` –
///   Maximum allowed time gap (days) within a batch when applying batch
///   RMS corrections. Helps control error correction behavior.
///
/// # Defaults
///
/// The [`Default`] implementation provides a good starting point:
///
/// ```rust, ignore
/// use outfit::initial_orbit_determination::IODParams;
///
/// let params = IODParams::default();
/// ```
///
/// Default values:
///
/// * `n_noise_realizations`: 20
/// * `noise_scale`: 1.0
/// * `extf`: -1.0 (use broad window)
/// * `dtmax`: 30.0 days
/// * `dt_min`: 0.03 days (~43 minutes)
/// * `dt_max_triplet`: 150.0 days
/// * `optimal_interval_time`: 20.0 days
/// * `max_obs_for_triplets`: 100
/// * `max_triplets`: 10
/// * `gap_max`: 8.0 / 24.0 days (8 hours)
///
/// # See also
///
/// * [`ObservationIOD::estimate_best_orbit`] – Main orbit determination entry point.
#[derive(Debug, Clone)]
pub struct IODParams {
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
}

impl IODParams {
    /// Construct a new [`IODParams`] with sensible default values.
    ///
    /// This is equivalent to calling [`IODParams::default()`].
    pub fn new() -> Self {
        Self::default()
    }

    /// Start building a customized [`IODParams`] using a fluent builder.
    pub fn builder() -> IODParamsBuilder {
        IODParamsBuilder::new()
    }
}

impl Default for IODParams {
    fn default() -> Self {
        IODParams {
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
        }
    }
}

/// Builder for [`IODParams`], with validation.
#[derive(Debug, Clone)]
pub struct IODParamsBuilder {
    params: IODParams,
}

impl IODParamsBuilder {
    /// Create a new builder initialized with default values.
    pub fn new() -> Self {
        Self {
            params: IODParams::default(),
        }
    }

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

    /// Finalize the builder and produce an [`IODParams`] instance.
    ///
    /// This method applies a set of **validation rules** to ensure that the
    /// configured parameters are consistent and safe to use before constructing
    /// the final [`IODParams`] structure.
    ///
    /// # Validation rules
    ///
    /// The following checks are performed:
    ///
    /// * `noise_scale >= 0.0` – the Gaussian noise scaling factor cannot be negative.
    /// * `dt_min >= 0.0`, `dt_max_triplet >= 0.0`, `dtmax >= 0.0` –
    ///   time spans must be non-negative.
    ///
    /// **Special case – `n_noise_realizations = 0`:**
    ///
    /// * Allowed.
    /// * In this case, the noise generation stage does not produce any
    ///   perturbed copies of a triplet.  
    ///   The orbit estimation uses **only the original triplet** (so
    ///   [`generate_noisy_realizations`] returns a vector of size 1).
    ///
    /// **Special case – `max_obs_for_triplets < 3`:**
    ///
    /// * This value is accepted.  
    /// * The downsampling function (`downsample_uniform_with_edges_indices`)
    ///   will still return three indices: first, middle, last.  
    ///   This guarantees that at least one valid triplet can be formed.
    /// * In other words:
    ///   - `max_obs_for_triplets = 1` or `2` behaves the same as `3`.
    ///   - If `max_obs_for_triplets >= n` (number of observations), no
    ///     downsampling occurs.
    ///
    /// # Returns
    ///
    /// * `Ok(IODParams)` if all values are valid and the parameters object
    ///   can safely be used.
    /// * `Err(OutfitError)` if any validation rule fails.
    ///
    /// # Examples
    ///
    /// ```rust, ignore
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
    /// to [`ObservationIOD::estimate_best_orbit`].
    pub fn build(self) -> Result<IODParams, OutfitError> {
        if self.params.noise_scale < 0.0 {
            return Err(OutfitError::InvalidIODParameter(
                "noise_scale must be non-negative".into(),
            ));
        }

        if self.params.dt_min < 0.0 || self.params.dt_max_triplet < 0.0 || self.params.dtmax < 0.0 {
            return Err(OutfitError::InvalidIODParameter(
                "time parameters must be non-negative".into(),
            ));
        }

        Ok(self.params)
    }
}
