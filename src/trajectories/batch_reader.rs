//! # Single-Observer Astrometric Batch Ingestion
//!
//! This module provides the [`ObservationBatch`] type, which groups multiple
//! astrometric detections from a **single observer** into a compact container.
//! Such a batch can then be expanded into concrete [`Observation`]s and stored
//! in a [`TrajectorySet`].
//!
//! ## Overview
//! -----------------
//! A wide-field survey typically delivers angle-only astrometry (RA/DEC, with
//! per-epoch timestamps). [`ObservationBatch`] wraps such measurements, together
//! with uniform error estimates, into a structured form ready for ingestion into
//! orbit-determination pipelines.
//!
//! To actually turn batches into stored observations, use the trait
//! [`TrajectoryFile`](crate::trajectories::trajectory_file::TrajectoryFile):
//! - [`TrajectoryFile::new_from_vec`](crate::trajectories::trajectory_file::TrajectoryFile::new_from_vec) — build a new [`TrajectorySet`] from a batch.
//! - [`TrajectoryFile::add_from_vec`](crate::trajectories::trajectory_file::TrajectoryFile::add_from_vec) — append a batch into an existing [`TrajectorySet`].
//!
//! Both methods transparently handle the internal expansion of a batch into
//! per-sample [`Observation`]s (site position lookups, heliocentric positions,
//! RA/DEC/error propagation, etc.).
//!
//! ## Units & Conventions
//! -----------------
//! - **Angles:** Right ascension and declination in **radians**.  
//!   If your upstream data are in **degrees/arcseconds**, use
//!   [`ObservationBatch::from_degrees_owned`] to convert once at construction.
//! - **Uncertainties:** 1-σ errors in RA/DEC (radians). For arcsecond inputs,
//!   the degree-based constructor performs the conversion for you.
//! - **Epochs:** Times in **MJD (TT)** (days). Convert UTC/TAI upstream.
//! - **Observer:** All rows in a batch must come from the **same** observer.
//!
//! ## Invariants
//! -----------------
//! - `trajectory_id.len() == ra.len() == dec.len() == time.len()`  
//! - All angles and uncertainties are in **radians**.  
//! - All epochs are in **MJD (TT)**.  
//! - Batch content belongs to a **single observer**.
//!
//! ## Construction Paths
//! -----------------
//! - [`ObservationBatch::from_radians_borrowed`] — zero-copy when your pipeline already
//!   provides radians and MJD (TT).
//! - [`ObservationBatch::from_degrees_owned`] — converts degrees/arcseconds → radians
//!   once and stores owned buffers.
//!
//! ## Example
//! -----------------
//! ```rust,no_run
//! use std::sync::Arc;
//! use outfit::{
//!     Outfit, Observer, TrajectorySet, TrajectoryFile, ErrorModel,
//!     trajectories::batch_reader::ObservationBatch,
//! };
//!
//! // Inputs in degrees / arcseconds (mixed objects: 0 and 1).
//! let traj_id = vec![0_u32, 0, 1];
//! let ra_deg  = vec![210.01, 210.02, 211.00];
//! let dec_deg = vec![-5.00, -4.99, -4.00];
//! let mjd_tt  = vec![60345.12, 60345.13, 60345.20];
//!
//! let batch = ObservationBatch::from_degrees_owned(
//!     &traj_id, &ra_deg, &dec_deg, 0.5, 0.5, &mjd_tt
//! );
//!
//! // Global environment and observer.
//! let mut outfit = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
//! let observer = outfit.get_observer_from_mpc_code(&"I41".to_string());
//!
//! // Build a new TrajectorySet directly from the batch.
//! let traj_set = TrajectorySet::new_from_vec(&mut outfit, &batch, observer.clone())
//!     .expect("ingestion OK");
//!
//! // Or append to an existing set:
//! let mut other = TrajectorySet::default();
//! other.add_from_vec(&mut outfit, &batch, observer).expect("append OK");
//! ```
//!
//! ## Notes
//! -----------------
//! - Internally, ingestion is performed by a crate-private routine
//!   (`observation_from_batch`) which expands the batch into per-sample
//!   [`Observation`]s and caches observer positions per epoch. Users should rely
//!   on the public `TrajectoryFile` methods instead.
//! - For multi-observer datasets, create one [`ObservationBatch`] per observer,
//!   then ingest them separately.
//!
//! ## See also
//! ------------
//! * [`ObservationBatch::from_radians_borrowed`] – Zero-copy construction.  
//! * [`ObservationBatch::from_degrees_owned`] – Convert degrees/arcseconds once.  
//! * [`TrajectoryFile::new_from_vec`](crate::trajectories::trajectory_file::TrajectoryFile::new_from_vec) – Public entry point for batch ingestion.  
//! * [`TrajectoryFile::add_from_vec`](crate::trajectories::trajectory_file::TrajectoryFile::add_from_vec) – Append batch into an existing set.
use std::{borrow::Cow, sync::Arc};

use ahash::RandomState;
use hifitime::Epoch;
use nalgebra::Vector3;
use ordered_float::OrderedFloat;
use smallvec::SmallVec;

use crate::{
    constants::Radian, conversion::arcsec_to_rad, observations::Observation,
    trajectories::parquet_reader::FastHashMap, ArcSec, Degree, ObjectNumber, Observer, Outfit,
    OutfitError, TrajectorySet, MJD,
};

/// Batch of observations from a single observer (angles in **radians**).
///
/// This container groups multiple astrometric measurements sharing the same
/// observer into a single batch, ready to be expanded into
/// [`Observation`]s and stored in a
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
/// # use outfit::trajectories::batch_reader::ObservationBatch;
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
        debug_assert_eq!(ra_deg.len(), dec_deg.len(), "RA/DEC length mismatch");
        debug_assert_eq!(ra_deg.len(), time_mjd.len(), "RA/time length mismatch");

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

/// Expand a single-observer batch into concrete [`Observation`]s and append them into a [`TrajectorySet`].
///
/// This routine ingests an [`ObservationBatch`] whose angles and uncertainties are in **radians**
/// and whose epochs are **MJD (TT)**, then materializes per-sample [`Observation`]s enriched with
/// site geocentric and heliocentric positions. All measurements are assumed to come from the **same
/// observer** (hence a single `observer: Arc<Observer>` argument).
///
/// For performance, it:
/// - resolves the observer into a compact `u16` **once** (hot path),
/// - pre-fetches the UT1 provider **once**,
/// - caches observer positions by epoch: **MJD(TT) → (geo_pos, helio_pos)**,
///   so repeated timestamps incur **no extra** position computation.
///
/// Internally, epoch keys are wrapped in `OrderedFloat` to enable their use in a hash map
/// (total order on `f64` while rejecting `NaN` inputs by construction).
///
/// Arguments
/// -----------------
/// * `trajectories` — Target container to receive observations, bucketed by `trajectory_id`.
/// * `env_state` — Global [`Outfit`] state (ephemerides, EOP/UT1 providers, etc.).
/// * `batch` — Angles and 1-σ uncertainties in **radians**; epochs as **MJD (TT)**; includes `trajectory_id`.
/// * `observer` — The (single) observer for **all** samples in `batch`.
///
/// Return
/// ----------
/// * `Ok(())` if all observations were successfully appended into `trajectories`,
/// * `Err(OutfitError)` if site position or heliocentric position computations fail.
///
/// Panics
/// ----------
/// * **Debug builds only**: length mismatches across `ra/dec/time/trajectory_id` trigger `debug_assert!`.
///
/// Complexity
/// ----------
/// * Time: **O(n)**, with at most **O(u)** geocentric/heliocentric computations where `u` is the number of
///   **unique** epochs in the batch (`u ≤ n`) thanks to the epoch→position cache.
/// * Space: **O(u)** for the epoch→position cache.
///
/// Notes
/// ----------
/// * Input angles (RA/DEC) and uncertainties **must already be in radians**. If your source is degrees/arcsec,
///   build the batch via [`ObservationBatch::from_degrees_owned`] (conversion done once at construction).
/// * Epochs are expected as **TT**. Convert upstream if your pipeline feeds UTC/TAI.
/// * This function mutates `trajectories` and reads from `env_state`. If you need parallelization, consider
///   extracting immutable position providers beforehand, or designing providers that accept shared references.
///
/// See also
/// ------------
/// * [`ObservationBatch::from_radians_borrowed`] – Zero-copy batch when inputs are already radians.
/// * [`ObservationBatch::from_degrees_owned`] – Degree/arcsec → rad conversion once at construction.
/// * [`parquet_to_trajset`] – Parquet ingestion using the same unit/weighting logic.
/// * [`conversion::arcsec_to_rad`] – Arcseconds → radians helper.
pub(crate) fn observation_from_batch(
    trajectories: &mut TrajectorySet,
    env_state: &mut Outfit,
    batch: &ObservationBatch<'_>,
    observer: Arc<Observer>,
) -> Result<(), OutfitError> {
    // --- Fast sanity checks (debug only) ---------------------------------------
    // All slices must be aligned one-to-one. These checks are enforced at construction
    // time for production, but we keep them here in debug builds to catch regressions.
    debug_assert_eq!(batch.ra.len(), batch.dec.len(), "RA/DEC length mismatch");
    debug_assert_eq!(batch.ra.len(), batch.time.len(), "RA/time length mismatch");
    debug_assert_eq!(
        batch.ra.len(),
        batch.trajectory_id.len(),
        "RA/trajectory_id length mismatch"
    );

    // Resolve observer ID once (hot path avoids map lookups later).
    let uint16_obs = env_state.uint16_from_observer(observer.clone());

    // Pre-fetch UT1 provider once.
    let ut1 = env_state.get_ut1_provider();

    // Heuristic: many surveys repeat epochs per exposure → cache a fraction of N.
    let n = batch.ra.len();
    let est_cache_cap = (n / 4).clamp(64, 4096);
    let mut pos_cache: FastHashMap<OrderedFloat<f64>, (Vector3<f64>, Vector3<f64>)> =
        FastHashMap::with_capacity_and_hasher(est_cache_cap, RandomState::default());

    // Safe and fast: iterators avoid per-iteration bounds checks
    let ra_it = batch.ra.iter().copied();
    let dec_it = batch.dec.iter().copied();
    let time_it = batch.time.iter().copied();
    let id_it = batch.trajectory_id.iter().copied();

    for ((ra, dec), (mjd_tt, traj_id)) in ra_it.zip(dec_it).zip(time_it.zip(id_it)) {
        // Use OrderedFloat to permit f64 as a key with a total order.
        let key = OrderedFloat(mjd_tt);

        // Compute (or reuse) positions at this epoch.
        let (geo_pos, helio_pos) = if let Some(&(geo, helio)) = pos_cache.get(&key) {
            (geo, helio)
        } else {
            let epoch = Epoch::from_mjd_in_time_scale(mjd_tt, hifitime::TimeScale::TT);

            // Geocentric position (velocity returned but unused here).
            let (geo, _vel) = observer.pvobs(&epoch, ut1)?;

            // Heliocentric position (uses geocentric position).
            let helio = observer.helio_position(env_state, &epoch, &geo)?;

            pos_cache.insert(key, (geo, helio));
            (geo, helio)
        };

        // Build observation and append to the right trajectory bucket.
        let obs = Observation::with_positions(
            uint16_obs,
            ra,              // radians
            batch.error_ra,  // radians
            dec,             // radians
            batch.error_dec, // radians
            mjd_tt,          // MJD(TT)
            geo_pos,
            helio_pos,
        );

        let obj = ObjectNumber::Int(traj_id);
        trajectories
            .entry(obj)
            .or_insert_with(|| SmallVec::with_capacity(32))
            .push(obs);
    }

    Ok(())
}
