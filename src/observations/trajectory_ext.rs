//! Trajectory ingestion and batch Initial Orbit Determination (IOD).
//!
//! Overview
//! -----------------
//! This module provides **high-level utilities** to:
//! * Build a [`TrajectorySet`] from multiple data sources (80-column MPC files, Parquet, ADES,
//!   or in-memory batches).
//! * Append observations to an existing set.
//! * Run a **Gauss-based IOD** over all trajectories and collect **per-object results** in a
//!   single map ([`FullOrbitResult`]).
//!
//! Data model
//! -----------------
//! * A [`TrajectorySet`] is a `HashMap<ObjectNumber, Observations>` storing one time-ordered
//!   list of astrometric observations per object.
//! * [`ObservationBatch`] is a thin, zero-copy view used to create/append observations that
//!   all share the **same observer** and **uniform uncertainties**.
//! * Results of the batch IOD are returned as a [`FullOrbitResult`]:
//!   `HashMap<ObjectNumber, Result<(Option<GaussResult>, f64), OutfitError>, RandomState>`.
//!
//! Ingestion sources
//! -----------------
//! The [`TrajectoryExt`] trait exposes constructors and appenders for:
//! * **80-column MPC files**: [`TrajectoryExt::new_from_80col`], [`TrajectoryExt::add_from_80col`]
//! * **Parquet** files (columns: `"ra"`, `"dec"`, `"jd"`, `"trajectory_id"`):
//!   [`TrajectoryExt::new_from_parquet`], [`TrajectoryExt::add_from_parquet`]
//! * **ADES** (MPC XML/JSON): [`TrajectoryExt::new_from_ades`], [`TrajectoryExt::add_from_ades`]
//! * **In-memory batches** (one observer): [`TrajectoryExt::new_from_vec`], [`TrajectoryExt::add_from_vec`]
//!
//! Batch IOD
//! -----------------
//! Use [`TrajectoryExt::estimate_all_orbits`] to run the full Gauss IOD pipeline on every
//! `(ObjectNumber → Observations)` pair in the set. Each object is processed with the same
//! [`Outfit`] state and [`IODParams`]. The outcome for each object is either:
//! * `Ok(Some(GaussResult))` + RMS – a viable preliminary or corrected orbit,
//! * `Ok(None)` – pipeline executed but no viable solution kept,
//! * `Err(OutfitError)` – a failure specific to that object (does not abort the batch).
//!
//! Assumptions & format expectations
//! -----------------
//! * Right ascension and declination are in **degrees** on input (unless stated otherwise).
//! * Times are **TT** on input where noted and converted internally to **MJD**.
//! * For Parquet ingestion, required columns are `"ra"`, `"dec"`, `"jd"`, `"trajectory_id"`.
//! * For 80-col and ADES, inputs must follow MPC specifications.
//!
//! Errors
//! -----------------
//! Ingestion and IOD may return [`OutfitError`]. Batch execution isolates failures so that
//! one object cannot prevent others from being processed.
//!
//! Example
//! -----------------
//! ```no_run
//! use std::sync::Arc;
//! use camino::Utf8Path;
//! use ahash::RandomState;
//! use rand::SeedableRng;
//! use outfit::outfit::Outfit;
//! use outfit::observations::trajectory_ext::{TrajectoryExt, FullOrbitResult};
//! use outfit::observers::Observer;
//! use outfit::initial_orbit_determination::IODParams;
//! use outfit::constants::TrajectorySet;
//!
//! # fn demo() -> Result<(), outfit::outfit_errors::OutfitError> {
//! let mut state = Outfit::new("horizon:DE440", outfit::error_models::ErrorModel::FCCT14)?;
//! let observer: Arc<Observer> = state.get_observer_from_mpc_code(&"I41".into());
//!
//! // Build a TrajectorySet from Parquet
//! let mut trajs: TrajectorySet = TrajectoryExt::new_from_parquet(
//!     &mut state,
//!     Utf8Path::new("observations.parquet"),
//!     observer.clone(),
//!     0.5, 0.5,
//!     Some(2048),
//! )?;
//!
//! // Run batch IOD
//! let mut rng = rand::rngs::StdRng::from_os_rng();
//! let params = IODParams::builder().max_triplets(32).build()?;
//! let results: FullOrbitResult = trajs.estimate_all_orbits(&state, &mut rng, &params);
//!
//! // Iterate results
//! for (obj, outcome) in results {
//!     match outcome {
//!         Ok((orbit, rms)) => eprintln!("{:?} → {:?}, rms={:.4}", obj, orbit, rms),
//!         Err(err)               => eprintln!("{:?} → error: {}", obj, err),
//!     }
//! }
//! # Ok(()) }
//! ```
//!
//! See also
//! ------------
//! * [`TrajectoryExt`] – Constructors and appenders for [`TrajectorySet`].
//! * [`ObservationIOD::estimate_best_orbit`] – Per-trajectory IOD and scoring.
//! * [`IODParams`] – Tuning for triplet generation and correction.
//! * [`GaussResult`] – Preliminary vs. corrected orbit representations.
//! * [`Outfit`] – Ephemerides, reference frames, and observer registry.
use std::borrow::Cow;
use std::fmt;
use std::{collections::HashMap, sync::Arc};

use crate::constants::{ArcSec, Degree, ObjectNumber, Observations, Radian, TrajectorySet, MJD};
use crate::conversion::arcsec_to_rad;
use crate::initial_orbit_determination::gauss_result::GaussResult;
use crate::initial_orbit_determination::IODParams;
use crate::observations::observation_from_batch;
use crate::observations::observations_ext::ObservationIOD;
use crate::observers::Observer;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use ahash::RandomState;
use camino::Utf8Path;
use rand::Rng;

use super::ades_reader::parse_ades;
use super::extract_80col;
use super::parquet_reader::parquet_to_trajset;

#[cfg(feature = "progress")]
use crate::observations::progress_bar::IterTimer;
#[cfg(feature = "progress")]
use indicatif::{ProgressBar, ProgressStyle};
#[cfg(feature = "progress")]
use std::time::Duration;

/// Full batch orbit determination results.
///
/// Each entry maps an [`ObjectNumber`] to the outcome of a full
/// Initial Orbit Determination (IOD) attempt on its set of observations.
/// The result type is:
///
/// * `Ok((Option<GaussResult>, f64))` – successful execution of the IOD pipeline:
///   * `Option<GaussResult>` – the best preliminary or corrected orbit found
///     for this object (or `None` if no viable orbit was selected),
///   * `f64` – the RMS of normalized astrometric residuals for that solution.
/// * `Err(OutfitError)` – a failure during orbit estimation for this object.
///   Errors are per-object and do not abort the rest of the batch.
///
/// Internally, this is implemented as:
///
/// ```ignore
/// HashMap<ObjectNumber, Result<(Option<GaussResult>, f64), OutfitError>, RandomState>
/// ```
pub type FullOrbitResult =
    HashMap<ObjectNumber, Result<(GaussResult, f64), OutfitError>, RandomState>;

/// Borrow a Gauss solution (if any) and its RMS for a given key.
///
/// Arguments
/// -----------------
/// * `all`: The map of all IOD outcomes.
/// * `key`: The object identifier.
///
/// Return
/// ----------
/// * `Ok(Some((&GaussResult, f64)))` – a solution is present for the key.
/// * `Ok(None)` – key absent OR present but no acceptable solution (`None`).
/// * `Err(&OutfitError)` – the IOD attempt failed for that key.
///
/// See also
/// ------------
/// * [`GaussResult`] – Gauss IOD output structure.
pub fn gauss_result_for<'a>(
    all: &'a FullOrbitResult,
    key: &ObjectNumber,
) -> Result<Option<(&'a GaussResult, f64)>, &'a OutfitError> {
    match all.get(key) {
        None => Ok(None),
        Some(Err(e)) => Err(e),
        Some(Ok((g, rms))) => Ok(Some((g, *rms))),
    }
}

/// Take ownership of the solution for `key`, removing it from the map.
pub fn take_gauss_result(
    all: &mut FullOrbitResult,
    key: &ObjectNumber,
) -> Result<Option<(GaussResult, f64)>, OutfitError> {
    match all.remove(key) {
        None => Ok(None),
        Some(Err(e)) => Err(e),
        Some(Ok((g, rms))) => Ok(Some((g, rms))),
    }
}

/// Summary statistics for per-trajectory observation counts.
///
/// Each [`TrajectorySet`] entry (one object) has an associated
/// [`Observations`] container. This structure stores basic distribution
/// statistics on the **number of observations per trajectory**, as
/// returned by [`obs_count_stats`](crate::observations::trajectory_ext::TrajectoryExt::obs_count_stats).
///
/// Fields
/// -----------------
/// * `min` – smallest number of observations in any trajectory.
/// * `p25` – 25th percentile (first quartile) of observation counts.
/// * `median` – 50th percentile (second quartile).
/// * `p95` – 95th percentile, indicating the upper tail of the distribution.
/// * `max` – largest number of observations in any trajectory.
///
/// Percentiles are computed using the *nearest-rank* method:
/// the index is `round(q × (N-1))` for quantile `q ∈ [0,1]`, clamped to valid range.
/// This convention makes results stable even for small sample sizes.
///
/// Display
/// -----------------
/// * `format!("{}", stats)` – compact single-line summary, e.g.:
///   ```text
///   min=2, p25=4, median=8, p95=15, max=20
///   ```
///
/// * `format!("{:#}", stats)` – pretty multi-line table, e.g.:
///   ```text
///   Observation count per trajectory — summary
///   -----------------------------------------
///   min    : 2
///   p25    : 4
///   median : 8
///   p95    : 15
///   max    : 20
///   ```
///
/// See also
/// ------------
/// * [`obs_count_stats`](crate::observations::trajectory_ext::TrajectoryExt::obs_count_stats) – Computes these statistics from a [`TrajectorySet`].
#[derive(Debug, Clone, Copy)]
pub struct ObsCountStats {
    pub min: usize,
    pub p25: usize,
    pub median: usize,
    pub p95: usize,
    pub max: usize,
}

impl fmt::Display for ObsCountStats {
    /// Compact by default; pretty multi-line when using the alternate flag (`{:#}`).
    ///
    /// See also
    /// ------------
    /// * [`obs_count_stats`](crate::observations::trajectory_ext::TrajectoryExt::obs_count_stats) – Builder of these summary statistics.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if f.alternate() {
            // Pretty, multi-line, aligned output (ASCII-only for portability).
            writeln!(f, "Observation count per trajectory — summary")?;
            writeln!(f, "-----------------------------------------")?;
            writeln!(f, "min    : {}", self.min)?;
            writeln!(f, "p25    : {}", self.p25)?;
            writeln!(f, "median : {}", self.median)?;
            writeln!(f, "p95    : {}", self.p95)?;
            write!(f, "max    : {}", self.max)
        } else {
            // Compact single-line for logs and quick prints.
            write!(
                f,
                "min={}, p25={}, median={}, p95={}, max={}",
                self.min, self.p25, self.median, self.p95, self.max
            )
        }
    }
}

/// Batch of observations from a single observer (angles in **radians**).
///
/// This container groups multiple astrometric measurements sharing the same
/// observer into a single batch, ready to be expanded into
/// [`Observation`](crate::observations::Observation)s and stored in a
/// [`TrajectorySet`].
///
/// Each measurement includes:
/// - A trajectory identifier (`trajectory_id`) so that a single batch can hold
///   observations for multiple objects simultaneously.
/// - Right ascension and declination in **radians**, with uniform 1-σ uncertainties
///   (also in **radians**).
/// - Epochs in **MJD (TT)** (days).
///
/// Fields
/// -----------------
/// * `trajectory_id` — Integer trajectory IDs (object numbers). Length must match `ra`/`dec`/`time`.
/// * `ra` — Right ascension values (**radians**). Length must match `dec`, `time`, and `trajectory_id`.
/// * `error_ra` — 1-σ uncertainty on right ascension (**radians**) applied uniformly to the batch.
/// * `dec` — Declination values (**radians**). Length must match `ra`, `time`, and `trajectory_id`.
/// * `error_dec` — 1-σ uncertainty on declination (**radians**) applied uniformly to the batch.
/// * `time` — Observation epochs as **MJD (TT)** (days). Length must match `ra`/`dec`/`trajectory_id`.
///
/// Invariants
/// -----------------
/// * `trajectory_id.len() == ra.len() == dec.len() == time.len()`
/// * Angles and uncertainties are expressed in **radians**.
/// * Time scale is **TT** (use appropriate conversion if your source data are in UTC/TAI).
///
/// Construction
/// -----------------
/// Prefer the dedicated constructors:
/// * [`ObservationBatch::from_radians_borrowed`] — zero-copy when your inputs are already in radians.
/// * [`ObservationBatch::from_degrees_owned`] — converts degrees/arcseconds once into owned buffers.
///
/// Example
/// -----------------
/// ```rust, no_run
/// # use outfit::observations::trajectory_ext::ObservationBatch;
/// # let (traj_id, ra_deg, dec_deg, mjd) = (vec![0, 0, 1], vec![14.62, 14.63, 15.01], vec![9.98, 10.01, 11.02], vec![43785.35, 43785.36, 43785.40]);
/// // Inputs in degrees / arcseconds (converted once to radians internally):
/// let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &mjd);
///
/// // Or, if you already have radians:
/// // let batch = ObservationBatch::from_radians_borrowed(&ra_rad, &dec_rad, err_ra_rad, err_dec_rad, &mjd);
/// ```
///
/// See also
/// ------------
/// * [`ObservationBatch::from_radians_borrowed`] – Borrow slices already in radians (zero-copy).
/// * [`ObservationBatch::from_degrees_owned`] – Convert degrees/arcseconds → radians once.
/// * [`conversion::arcsec_to_rad`](crate::conversion::arcsec_to_rad) – Arcseconds → radians helper.
#[derive(Debug, Clone)]
pub struct ObservationBatch<'a> {
    pub trajectory_id: Cow<'a, [u32]>,

    /// Right ascension values (**radians**). Must have the same length as `dec` and `time`.
    pub ra: Cow<'a, [Radian]>,

    /// 1-σ uncertainty on right ascension (**radians**), applied uniformly to the batch.
    /// Note: the weighting scheme accounts for the RA geometry (e.g., cos δ factors) downstream.
    pub error_ra: Radian,

    /// Declination values (**radians**). Must have the same length as `ra` and `time`.
    pub dec: Cow<'a, [Radian]>,

    /// 1-σ uncertainty on declination (**radians**), applied uniformly to the batch.
    pub error_dec: Radian,

    /// Observation epochs as **MJD (TT)**, in days. Must have the same length as `ra`/`dec`.
    pub time: Cow<'a, [MJD]>,
}

impl<'a> ObservationBatch<'a> {
    /// Construct a batch by **borrowing** slices that are already in radians.
    ///
    /// The returned batch holds `Cow::Borrowed` views of the provided slices,
    /// performing **no allocation** and **no unit conversion**.
    /// Use this when your upstream pipeline already provides:
    /// - Trajectory identifiers (`trajectory_id`)
    /// - Right ascension / declination in **radians**
    /// - Uncertainties in **radians**
    /// - Epochs in **MJD (TT)**
    ///
    /// Arguments
    /// -----------------
    /// * `trajectory_id` — Integer trajectory IDs; length must match all angle/time slices.
    /// * `ra_rad` — Right ascension values in **radians** (borrowed).
    /// * `dec_rad` — Declination values in **radians** (borrowed).
    /// * `error_ra_rad` — 1-σ uncertainty on RA in **radians**, applied uniformly to the batch.
    /// * `error_dec_rad` — 1-σ uncertainty on DEC in **radians**, applied uniformly to the batch.
    /// * `time_mjd` — Observation epochs as **MJD (TT)** (borrowed).
    ///
    /// Return
    /// ----------
    /// * A batch borrowing the provided slices (**zero-copy**).
    ///
    /// Invariants
    /// ----------
    /// * `trajectory_id.len() == ra_rad.len() == dec_rad.len() == time_mjd.len()`
    ///
    /// Panics
    /// ----------
    /// * Debug builds only: panics if the slice lengths do not match.
    ///
    /// Complexity
    /// ----------
    /// * O(1) — no allocation, no conversion.
    ///
    /// See also
    /// ------------
    /// * [`ObservationBatch::from_degrees_owned`] – Convert degrees/arcseconds → radians and own the buffers.
    /// * [`conversion::arcsec_to_rad`](crate::conversion::arcsec_to_rad) – Arcseconds → radians helper.
    pub fn from_radians_borrowed(
        trajectory_id: &'a [u32],
        ra_rad: &'a [Radian],
        dec_rad: &'a [Radian],
        error_ra_rad: Radian,
        error_dec_rad: Radian,
        time_mjd: &'a [MJD],
    ) -> Self {
        debug_assert_eq!(ra_rad.len(), dec_rad.len(), "RA/DEC length mismatch");
        debug_assert_eq!(ra_rad.len(), time_mjd.len(), "RA/time length mismatch");

        Self {
            trajectory_id: Cow::Borrowed(trajectory_id),
            ra: Cow::Borrowed(ra_rad),
            dec: Cow::Borrowed(dec_rad),
            time: Cow::Borrowed(time_mjd),
            error_ra: error_ra_rad,
            error_dec: error_dec_rad,
        }
    }

    /// Construct a batch from **degrees** (angles) and **arcseconds** (uncertainties),
    /// converting to **radians** and **owning** the resulting buffers.
    ///
    /// Use this when your inputs come from common astrometric formats (e.g., MPC/ADES)
    /// that report RA/DEC in degrees and uncertainties in arcseconds.
    /// Conversion is performed **once** at construction; downstream code operates purely in radians.
    ///
    /// Arguments
    /// -----------------
    /// * `trajectory_id` — Integer trajectory IDs; length must match all angle/time slices.
    /// * `ra_deg` — Right ascension in **degrees** (borrowed); converted to radians.
    /// * `dec_deg` — Declination in **degrees** (borrowed); converted to radians.
    /// * `error_ra_arcsec` — 1-σ uncertainty on RA in **arcseconds**; converted to radians.
    /// * `error_dec_arcsec` — 1-σ uncertainty on DEC in **arcseconds**; converted to radians.
    /// * `time_mjd` — Observation epochs as **MJD (TT)** (borrowed; cloned to owned buffer).
    ///
    /// Return
    /// ----------
    /// * A batch **owning** converted buffers (no dangling slices).
    ///
    /// Invariants
    /// ----------
    /// * `trajectory_id.len() == ra_deg.len() == dec_deg.len() == time_mjd.len()`
    ///
    /// Panics
    /// ----------
    /// * Panics if the slice lengths do not match.
    ///
    /// Complexity
    /// ----------
    /// * O(n) for the degree→radian and arcsec→radian conversions + one `to_vec()` for time.
    ///
    /// See also
    /// ------------
    /// * [`ObservationBatch::from_radians_borrowed`] – Zero-copy constructor when inputs are already in radians.
    /// * [`conversion::arcsec_to_rad`](crate::conversion::arcsec_to_rad) – Arcseconds → radians helper.
    pub fn from_degrees_owned(
        trajectory_id: &'a [u32],
        ra_deg: &[Degree],
        dec_deg: &[Degree],
        error_ra_arcsec: ArcSec,
        error_dec_arcsec: ArcSec,
        time_mjd: &[MJD],
    ) -> Self {
        assert_eq!(ra_deg.len(), dec_deg.len(), "RA/DEC length mismatch");
        assert_eq!(ra_deg.len(), time_mjd.len(), "RA/time length mismatch");

        let ra: Vec<f64> = ra_deg.iter().map(|&d| d.to_radians()).collect();
        let dec: Vec<f64> = dec_deg.iter().map(|&d| d.to_radians()).collect();
        let time: Vec<MJD> = time_mjd.to_vec();

        Self {
            trajectory_id: Cow::Owned(trajectory_id.to_vec()),
            ra: Cow::Owned(ra),
            dec: Cow::Owned(dec),
            time: Cow::Owned(time),
            error_ra: arcsec_to_rad(error_ra_arcsec),
            error_dec: arcsec_to_rad(error_dec_arcsec),
        }
    }
}

/// A trait for the TrajectorySet type def.
/// This trait provides methods to create a TrajectorySet from different sources.
/// It allows to create a TrajectorySet from an 80 column file, a parquet file, or an ADES file.
/// It also allows to add observations to an existing TrajectorySet from these sources.
/// The methods are:
/// * `from_80col`: Create a TrajectorySet from an 80 column file.
/// * `add_80col`: Add observations to a TrajectorySet from an 80 column file.
/// * `new_from_vec`: Create a TrajectorySet from a vector of observations.
/// * `add_from_vec`: Add observations to a TrajectorySet from a vector of observations.
/// * `new_from_parquet`: Create a TrajectorySet from a parquet file.
/// * `add_from_parquet`: Add observations to a TrajectorySet from a parquet file.
/// * `new_from_ades`: Create a TrajectorySet from an ADES file.
/// * `add_from_ades`: Add observations to a TrajectorySet from an ADES file.
///
/// Note
/// ----
/// * Warning: No check is done for duplicated observations for every add method.
///   * The user shoud be careful to not add the same observation or same file twice
pub trait TrajectoryExt {
    /// Create a TrajectorySet from an 80 column file
    /// The trajectory are added in place in the TrajectorySet.
    /// If a trajectory id already exists, the observations are added to the existing trajectory.
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `colfile`: a path to an 80 column file
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations from the 80 column file
    ///
    /// Note
    /// ----
    /// * The 80 column file must respect the MPC format.
    ///   * ref: <https://www.minorplanetcenter.net/iau/info/OpticalObs.html>
    fn new_from_80col(env_state: &mut Outfit, colfile: &Utf8Path) -> Self;

    /// Add a set of trajectories from an 80 column file to a TrajectorySet
    /// The trajectory are added in place in the TrajectorySet.
    /// If a trajectory id already exists, the observations are added to the existing trajectory.
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `colfile`: a path to an 80 column file
    ///
    /// Note
    /// ----
    /// * The 80 column file must respect the MPC format.
    ///   * ref: <https://www.minorplanetcenter.net/iau/info/OpticalObs.html>
    fn add_from_80col(&mut self, env_state: &mut Outfit, colfile: &Utf8Path);

    /// Create a new [`TrajectorySet`] from a batch of observations taken by a single observer.
    ///
    /// This constructor consumes an [`ObservationBatch`] and groups its observations
    /// into trajectories, keyed by their `trajectory_id`.  
    /// Each observation in the batch must have been recorded by the **same observer**,
    /// but may belong to **different objects** (distinguished by `trajectory_id`).
    ///
    /// Arguments
    /// -----------------
    /// * `env_state` — Mutable reference to the global [`Outfit`] state (used for ephemerides, UT1, etc.).
    /// * `batch` — An [`ObservationBatch`] containing RA/DEC/epoch values (radians + MJD/TT) and trajectory IDs.
    /// * `observer` — The observer that recorded all observations in the batch.
    ///
    /// Return
    /// -----------------
    /// * `Ok(Self)` — A new [`TrajectorySet`] containing one or more trajectories populated from the batch.
    /// * `Err(OutfitError)` — If observation construction or position computations fail.
    ///
    /// Invariants
    /// -----------------
    /// * `batch.trajectory_id.len() == batch.ra.len() == batch.dec.len() == batch.time.len()`
    /// * Angles and uncertainties in the batch must already be in **radians**.
    ///
    /// Example
    /// -----------------
    /// ```rust, no_run
    /// # use outfit::observations::trajectory_ext::ObservationBatch;
    /// # use outfit::TrajectorySet;
    /// # use outfit::TrajectoryExt;
    /// # use outfit::{Outfit, ErrorModel};
    /// # use std::sync::Arc;
    /// # let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
    /// # let observer = env.get_observer_from_mpc_code(&"I41".to_string());
    /// # let (traj_id, ra_deg, dec_deg, mjd) = (vec![0, 0, 1], vec![14.62, 14.63, 15.01], vec![9.98, 10.01, 11.02], vec![43785.35, 43785.36, 43785.40]);
    /// let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &mjd);
    ///
    /// // Build a trajectory set directly from the batch:
    /// let ts = TrajectorySet::new_from_vec(&mut env, &batch, observer).unwrap();
    /// ```
    fn new_from_vec(
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<Self, OutfitError>
    where
        Self: Sized;

    /// Add the observations from a batch to an existing [`TrajectorySet`].
    ///
    /// This method inserts all observations from the provided [`ObservationBatch`] into
    /// the current set, grouping them into trajectories by `trajectory_id`.
    /// Each observation in the batch must have been recorded by the **same observer**,
    /// but may belong to multiple distinct objects.
    ///
    /// Arguments
    /// -----------------
    /// * `env_state` — Mutable reference to the global [`Outfit`] state (used for ephemerides, UT1, etc.).
    /// * `batch` — An [`ObservationBatch`] containing RA/DEC/epoch values (radians + MJD/TT) and trajectory IDs.
    /// * `observer` — The observer that recorded all observations in the batch.
    ///
    /// Return
    /// -----------------
    /// * `Ok(())` — If all observations were successfully inserted into the `TrajectorySet`.
    /// * `Err(OutfitError)` — If observation construction or position computations fail.
    ///
    /// Example
    /// -----------------
    /// ```rust, no_run
    /// use outfit::observations::trajectory_ext::ObservationBatch;
    /// use outfit::TrajectorySet;
    /// use outfit::TrajectoryExt;
    /// use outfit::{Outfit, ErrorModel};
    /// use std::sync::Arc;
    /// use ahash::RandomState;
    /// use std::collections::HashMap;
    ///
    /// let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
    /// let observer = env.get_observer_from_mpc_code(&"I41".to_string());
    /// let (traj_id, ra_deg, dec_deg, mjd) = (vec![0, 0, 1], vec![14.62, 14.63, 15.01], vec![9.98, 10.01, 11.02], vec![43785.35, 43785.36, 43785.40]);
    /// let batch = ObservationBatch::from_degrees_owned(&traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &mjd);
    ///
    /// let mut ts = HashMap::with_hasher(RandomState::new());
    /// ts.add_from_vec(&mut env, &batch, observer).unwrap();
    /// ```
    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<(), OutfitError>;

    /// Create a new [`TrajectorySet`] from a Parquet file.
    ///
    /// This function reads a Parquet file containing astrometric observations
    /// and constructs a full [`TrajectorySet`]. Each observation is associated
    /// with the provided `observer` and assigned constant uncertainties in
    /// right ascension and declination.
    ///
    /// Arguments
    /// -----------------
    /// * `env_state` – Global environment providing ephemerides, UT1 provider, and observer mapping.
    /// * `parquet` – Path to the input Parquet file.
    /// * `observer` – Observer metadata (shared reference, resolved once to a compact id).
    /// * `error_ra` – 1-σ uncertainty in right ascension \[arcsec\], applied uniformly.
    /// * `error_dec` – 1-σ uncertainty in declination \[arcsec\], applied uniformly.
    /// * `batch_size` – Record batch size for Parquet reader; defaults to 2048 if `None`.
    ///
    /// Return
    /// ----------
    /// * `Ok(TrajectorySet)` – A new set of trajectories populated from the file.
    /// * `Err(OutfitError)` – If the file cannot be opened, parsed, or contains invalid data.
    ///
    /// Notes
    /// ----------
    /// * The Parquet file must contain the following columns: `"ra"`, `"dec"`, `"jd"`, `"trajectory_id"`.
    /// * The `"jd"` values are assumed to be in TT scale and are converted internally to MJD via [`JDTOMJD`](crate::constants::JDTOMJD).
    /// * The `ra` and `dec` columns have to be in degrees and of type `Float64`.
    /// * The `jd` column has to be in Julian Date (TT) and of type `Float64`.
    /// * The `trajectory_id` column has to be of type `UInt32` and is used to group
    ///   observations by object.
    ///
    /// See also
    /// ------------
    /// * [`add_from_parquet`](crate::observations::trajectory_ext::TrajectoryExt::add_from_parquet) – Adds observations from a Parquet file to an existing set.
    fn new_from_parquet(
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        mpc_code: Arc<Observer>,
        error_ra: ArcSec,
        error_dec: ArcSec,
        batch_size: Option<usize>,
    ) -> Result<Self, OutfitError>
    where
        Self: Sized;

    /// Add observations from a Parquet file to an existing [`TrajectorySet`].
    ///
    /// This function appends new observations (grouped by `trajectory_id`)
    /// to the current set. The same `observer` and astrometric uncertainties
    /// are applied to all ingested rows.
    ///
    /// Arguments
    /// -----------------
    /// * `env_state` – Global environment providing ephemerides, UT1 provider, and observer mapping.
    /// * `parquet` – Path to the input Parquet file.
    /// * `observer` – Observer metadata (shared reference, resolved once to a compact id).
    /// * `error_ra` – 1-σ uncertainty in right ascension \[arcsec\], applied uniformly.
    /// * `error_dec` – 1-σ uncertainty in declination \[arcsec\], applied uniformly.
    /// * `batch_size` – Record batch size for Parquet reader; defaults to 2048 if `None`.
    ///
    /// Return
    /// ----------
    /// * `Ok(())` – On successful ingestion, with the internal set updated in place.
    /// * `Err(OutfitError)` – If the file cannot be opened, parsed, or contains invalid data.
    ///
    /// Notes
    /// ----------
    /// * The Parquet file must contain the following columns: `"ra"`, `"dec"`, `"jd"`, `"trajectory_id"`.
    /// * The `"jd"` values are assumed to be in TT scale and are converted internally to MJD via [`JDTOMJD`](crate::constants::JDTOMJD).
    /// * The `ra` and `dec` columns have to be in degrees and of type `Float64`.
    /// * The `jd` column has to be in Julian Date (TT) and of type `Float64`.
    /// * The `trajectory_id` column has to be of type `UInt32` and is used to group
    ///   observations by object.
    ///
    /// See also
    /// ------------
    /// * [`new_from_parquet`](crate::observations::trajectory_ext::TrajectoryExt::new_from_parquet) – Creates a brand new set from a Parquet file.
    fn add_from_parquet(
        &mut self,
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        error_ra: ArcSec,
        error_dec: ArcSec,
        batch_size: Option<usize>,
    ) -> Result<(), OutfitError>;

    /// Add a set of trajectories to a TrajectorySet from an ADES file
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `ades`: a path to an ADES file
    /// * `error_ra`: the error in right ascension (if some values are given, the error ra is supposed to be the same for all observations)
    /// * `error_dec`: the error in declination (if some values are given, the error dec is supposed to be the same for all observations)
    ///
    /// Note
    /// ----
    /// * The ADES file must respect the MPC format.
    ///   * ref: <https://minorplanetcenter.net/iau/info/ADES.html>
    fn new_from_ades(
        env_state: &mut Outfit,
        ades: &Utf8Path,
        error_ra: Option<ArcSec>,
        error_dec: Option<ArcSec>,
    ) -> Self;

    /// Create a TrajectorySet from an ADES file
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `ades`: a path to an ADES file
    /// * `error_ra`: the error in right ascension (if some values are given, the error ra is supposed to be the same for all observations)
    /// * `error_dec`: the error in declination (if some values are given, the error dec is supposed to be the same for all observations)
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations from the ADES file
    ///
    /// Note
    /// ----
    /// * The ADES file must respect the MPC format.
    ///   * ref: <https://minorplanetcenter.net/iau/info/ADES.html>
    fn add_from_ades(
        &mut self,
        env_state: &mut Outfit,
        ades: &Utf8Path,
        error_ra: Option<ArcSec>,
        error_dec: Option<ArcSec>,
    );

    /// Estimate an initial orbit for **every trajectory** in the set and collect the results.
    ///
    /// This method runs the Gauss-based IOD pipeline on each `(ObjectNumber → Observations)`
    /// pair contained in the trajectory set. All objects are processed with the same
    /// configuration (`state`, `error_model`, `params`) and random number generator.
    /// Results are aggregated into a [`FullOrbitResult`] map.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: The global environment providing ephemerides, constants, and reference frames.
    /// * `rng`: Random number generator used for noisy triplet realizations (e.g., \[`StdRng`\]).
    /// * `params`: Parameters controlling triplet generation, scoring, and correction loops.
    ///
    /// Return
    /// ----------
    /// * A [`FullOrbitResult`] mapping each object to either:
    ///   * `Ok((Option<GaussResult>, f64))` – a valid solution with its RMS,
    ///   * `Err(OutfitError)` – an error diagnostic for that object.
    ///
    /// Notes
    /// ----------
    /// * The method iterates in place on the current set.  
    /// * Observations themselves are not modified, but computation time scales
    ///   with the number of trajectories and candidate triplets.  
    /// * Errors are isolated: one object failing does not prevent others from being processed.
    ///
    /// See also
    /// ------------
    /// * [`ObservationIOD::estimate_best_orbit`] – Per-trajectory IOD with best-orbit selection.
    /// * [`GaussResult`] – Variants for preliminary or corrected orbit solutions.
    /// * [`IODParams`] – Tuning parameters for IOD batch execution.
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult;

    /// Count the total number of [`Observation`](crate::observations::Observation) entries across all trajectories.
    ///
    /// This method iterates once over all values in the [`TrajectorySet`],
    /// summing the length of each [`Observations`] container.
    ///
    /// Return
    /// ----------
    /// * The total number of observations across all objects.
    fn total_observations(&self) -> usize;

    /// Compute distribution statistics for the number of observations per trajectory.
    ///
    /// Each trajectory (one object in the [`TrajectorySet`]) has an associated
    /// [`Observations`] container. This function collects their sizes and computes:
    ///
    /// * `min` – smallest number of observations in any trajectory,
    /// * `p25` – 25th percentile (first quartile),
    /// * `median` – 50th percentile (second quartile),
    /// * `p95` – 95th percentile (upper tail indicator),
    /// * `max` – largest number of observations in any trajectory.
    ///
    /// Percentiles are computed using the *nearest-rank* method:
    /// the index is `round(q × (N-1))` for quantile `q ∈ [0,1]`, clamped to valid range.
    /// This makes results robust even for small datasets.
    ///
    /// Return
    /// ----------
    /// * `None` if the set is empty.
    /// * `Some(ObsCountStats)` containing the summary statistics otherwise.
    ///
    /// See also
    /// ------------
    /// * [`total_observations`](crate::observations::trajectory_ext::TrajectoryExt::total_observations) – Sum of all observations across trajectories.
    fn obs_count_stats(&self) -> Option<ObsCountStats>;
}

impl TrajectoryExt for TrajectorySet {
    fn new_from_vec(
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<Self, OutfitError> {
        let mut traj_set: TrajectorySet = HashMap::default();
        observation_from_batch(&mut traj_set, env_state, batch, observer)?;
        Ok(traj_set)
    }

    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<(), OutfitError> {
        observation_from_batch(self, env_state, batch, observer)?;
        Ok(())
    }

    fn add_from_parquet(
        &mut self,
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        error_ra: ArcSec,
        error_dec: ArcSec,
        batch_size: Option<usize>,
    ) -> Result<(), OutfitError> {
        parquet_to_trajset(
            self, env_state, parquet, observer, error_ra, error_dec, batch_size,
        )
    }

    fn new_from_parquet(
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        error_ra: ArcSec,
        error_dec: ArcSec,
        batch_size: Option<usize>,
    ) -> Result<Self, OutfitError>
    where
        Self: Sized,
    {
        let mut trajs: TrajectorySet = HashMap::default();
        parquet_to_trajset(
            &mut trajs, env_state, parquet, observer, error_ra, error_dec, batch_size,
        )?;
        Ok(trajs)
    }

    fn new_from_80col(env_state: &mut Outfit, colfile: &Utf8Path) -> Self {
        let mut traj_set: TrajectorySet = HashMap::default();
        let (observations, object_number) =
            extract_80col(env_state, colfile).expect("Failed to extract 80col data");
        traj_set.insert(object_number, observations);
        traj_set
    }

    fn add_from_80col(&mut self, env_state: &mut Outfit, colfile: &Utf8Path) {
        let (observations, object_number) =
            extract_80col(env_state, colfile).expect("Failed to extract 80col data");
        self.insert(object_number, observations);
    }

    fn add_from_ades(
        &mut self,
        env_state: &mut Outfit,
        ades: &Utf8Path,
        error_ra: Option<ArcSec>,
        error_dec: Option<ArcSec>,
    ) {
        parse_ades(env_state, ades, self, error_ra, error_dec);
    }

    fn new_from_ades(
        env_state: &mut Outfit,
        ades: &Utf8Path,
        error_ra: Option<ArcSec>,
        error_dec: Option<ArcSec>,
    ) -> Self {
        let mut trajs: TrajectorySet = HashMap::default();
        parse_ades(env_state, ades, &mut trajs, error_ra, error_dec);
        trajs
    }

    #[cfg(feature = "progress")]
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult {
        let total = self.len() as u64;
        let pb = ProgressBar::new(total.max(1));
        pb.set_style(
            ProgressStyle::with_template(
                "{bar:40.cyan/blue} {pos}/{len} ({percent:>3}%) \
             | {per_sec} | ETA {eta_precise} | {msg}",
            )
            .expect("indicatif template"),
        );
        pb.enable_steady_tick(Duration::from_millis(200));

        let mut results: FullOrbitResult = HashMap::default();
        let mut it_timer = IterTimer::new(0.2);

        for (obj, observations) in self.iter_mut() {
            // Time of the previous loop (zero for the first, acceptable).

            use crate::observations::progress_bar::fmt_dur;
            let last = it_timer.tick();
            let avg = it_timer.avg();
            pb.set_message(format!("last: {}, avg: {}", fmt_dur(last), fmt_dur(avg)));

            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);

            pb.inc(1);
        }

        pb.finish_and_clear();
        results
    }

    #[cfg(not(feature = "progress"))]
    fn estimate_all_orbits(
        &mut self,
        state: &Outfit,
        rng: &mut impl Rng,
        params: &IODParams,
    ) -> FullOrbitResult {
        // Output map using the same fast hasher as TrajectorySet.
        let mut results: FullOrbitResult = HashMap::default();

        for (obj, observations) in self.iter_mut() {
            let res = observations.estimate_best_orbit(state, &state.error_model, rng, params);
            results.insert(obj.clone(), res);
        }

        results
    }

    #[inline]
    fn total_observations(&self) -> usize {
        self.values().map(|obs: &Observations| obs.len()).sum()
    }

    fn obs_count_stats(&self) -> Option<ObsCountStats> {
        // Collect sizes (one pass, O(N))
        let mut counts: Vec<usize> = self.values().map(|obs| obs.len()).collect();
        if counts.is_empty() {
            return None;
        }

        // Sort once, O(N log N). `unstable` is fine since we only need order.
        counts.sort_unstable();

        #[inline]
        fn q_index(n: usize, q: f64) -> usize {
            // Nearest-rank on [0, n-1] using linear index; robust for small n.
            let pos = q * (n as f64 - 1.0);
            let idx = pos.round() as isize;
            idx.clamp(0, (n as isize) - 1) as usize
        }

        let n = counts.len();
        let min = counts[0];
        let max = counts[n - 1];
        let p25 = counts[q_index(n, 0.25)];
        let median = counts[q_index(n, 0.50)];
        let p95 = counts[q_index(n, 0.95)];

        Some(ObsCountStats {
            min,
            p25,
            median,
            p95,
            max,
        })
    }
}
