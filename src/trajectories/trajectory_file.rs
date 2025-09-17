//! # Trajectory ingestion and batch Initial Orbit Determination (IOD)
//!
//! High-level utilities to **build and extend** a [`TrajectorySet`] from multiple sources
//! (MPC 80-column, Parquet, ADES, or in-memory batches) and to run a **Gauss-based IOD**
//! over all trajectories.
//!
//! ## Overview
//! -----------------
//! This module exposes the [`TrajectoryFile`] trait implemented for [`TrajectorySet`].
//! It provides:
//! - Constructors that **create** a new set from a given source (`new_from_*`),
//! - Appenders that **extend** an existing set (`add_from_*`),
//! - Convenience methods to ingest **in-memory batches** (single observer) via
//!   [`ObservationBatch`].
//!
//! Internally, ingestion from in-memory batches uses a crate-private routine
//! `observation_from_batch` (unit/scale/caching logic). End users should interact only
//! with the public `new_from_vec` / `add_from_vec` methods.
//!
//! ## Data model
//! -----------------
//! - A [`TrajectorySet`] is a `HashMap<ObjectNumber, Observations>` storing a
//!   time-ordered list of astrometric observations per object.
//! - [`ObservationBatch`] is a thin container (borrowed/owned) for angle-only astrometry
//!   from a **single observer** with uniform uncertainties; it is expanded into concrete
//!   [`Observation`](crate::observations::Observation)s and grouped by `trajectory_id`.
//! - Batch IOD returns a
//!   [`FullOrbitResult`](crate::trajectories::trajectory_fit::FullOrbitResult), i.e. a map
//!   `ObjectNumber → Result<(Option<GaussResult>, f64), OutfitError>`.
//!
//! ## Ingestion sources & signatures
//! -----------------
//! **MPC 80-column**
//! - [`TrajectoryFile::new_from_80col`] → `Self`  
//!   Reads a file, extracts `(Observations, ObjectNumber)` and **inserts** into a new set.
//!   **Panics** on extraction failure (internal `expect`).  
//! - [`TrajectoryFile::add_from_80col`] → `()`  
//!   Reads and **inserts** into an existing set. **Panics** on extraction failure.
//!
//! **Parquet** (`"ra"`, `"dec"`, `"jd"`, `"trajectory_id"`)
//! - [`TrajectoryFile::new_from_parquet`] → `Result<Self, OutfitError>`  
//!   Creates a new set; errors are propagated.
//! - [`TrajectoryFile::add_from_parquet`] → `Result<(), OutfitError>`  
//!   Appends to an existing set; errors are propagated.
//!   *Units on disk:* `ra/dec` in **degrees**, `jd` in **JD (TT)**. Internally converted to
//!   **radians** and **MJD (TT)** (via [`JDTOMJD`](crate::constants::JDTOMJD)). Per-file
//!   uncertainties are passed in **arcseconds**.
//!
//! **ADES (MPC XML/JSON)**
//! - [`TrajectoryFile::new_from_ades`] → `Self`  
//! - [`TrajectoryFile::add_from_ades`] → `()`  
//!   Both delegate to `parse_ades` and **do not return a `Result`** (errors are handled
//!   inside the parser or may panic depending on its policy).
//!
//! **In-memory batches (single observer)**
//! - [`TrajectoryFile::new_from_vec`] → `Result<Self, OutfitError>`  
//! - [`TrajectoryFile::add_from_vec`] → `Result<(), OutfitError>`  
//!   Expand an [`ObservationBatch`] (RA/DEC/σ in **radians**, epochs in **MJD (TT)**)
//!   into per-sample [`Observation`](crate::observations::Observation)s using the shared
//!   [`Outfit`] state and **append/group** by `trajectory_id`.
//!
//! ## Units & time scales
//! -----------------
//! - **Angles**: internal [`Observation`](crate::observations::Observation)s store RA/DEC in **radians**.  
//!   Parquet/80-column/ADES readers perform degree→radian conversions as needed.
//! - **Uncertainties**: expected in **arcseconds** at call-site for Parquet/ADES; for
//!   in-memory batches they must already be in **radians** (uniform per batch).
//! - **Times**: internal epochs are **MJD (TT)**. Parquet `"jd"` values are assumed **TT**
//!   and converted via [`JDTOMJD`](crate::constants::JDTOMJD). 80-col/ADES readers apply their respective conversions.
//!
//! ## Duplicates & ordering
//! -----------------
//! - **No deduplication** is performed by any `add_*` method. Users must avoid re-ingesting
//!   the same file/batch twice if duplicates are undesirable.
//! - Observations are stored **as provided**; ordering by time is not enforced here.
//!
//! ## Error semantics
//! -----------------
//! - Methods returning `Result<_, OutfitError>` propagate I/O/schema/ephemeris errors.
//! - `new_from_80col` / `add_from_80col` use `expect(...)` internally and therefore may **panic**
//!   on parse/read failures (fail-fast behavior).
//! - `new_from_ades` / `add_from_ades` currently **do not** return a `Result`; error handling
//!   is delegated to `parse_ades` (which may log or panic depending on implementation).
//!
//! ## Batch IOD
//! -----------------
//! Use [`crate::trajectories::trajectory_fit::TrajectoryFit::estimate_all_orbits`] to run the
//! full Gauss IOD over each `(ObjectNumber → Observations)` pair. Outcomes per object:
//! - `Ok(Some(GaussResult))` + RMS — a viable preliminary/corrected orbit,
//! - `Ok(None)` — pipeline executed but no acceptable solution kept,
//! - `Err(OutfitError)` — failure **isolated** to that object.
//!
//! ## Example
//! -----------------
//! ```no_run
//! use std::sync::Arc;
//! use camino::Utf8Path;
//! use rand::SeedableRng;
//! use outfit::outfit::Outfit;
//! use outfit::observers::Observer;
//! use outfit::trajectories::trajectory_file::TrajectoryFile;
//! use outfit::TrajectoryFit;
//! use outfit::initial_orbit_determination::IODParams;
//! use outfit::TrajectorySet;
//!
//! # fn demo() -> Result<(), outfit::outfit_errors::OutfitError> {
//! let mut state = Outfit::new("horizon:DE440", outfit::error_models::ErrorModel::FCCT14)?;
//! let observer: Arc<Observer> = state.get_observer_from_mpc_code(&"I41".into());
//!
//! // 1) From Parquet (propagates errors)
//! let mut trajs: TrajectorySet = TrajectorySet::new_from_parquet(
//!     &mut state,
//!     Utf8Path::new("observations.parquet"),
//!     observer.clone(),
//!     0.5, 0.5,
//!     Some(8192),
//! )?;
//!
//! // 2) From MPC 80-column (may panic on parse error)
//! trajs.add_from_80col(&mut state, Utf8Path::new("obs_80col.txt"));
//!
//! // 3) Run batch IOD
//! let mut rng = rand::rngs::StdRng::from_os_rng();
//! let params = IODParams::builder().max_triplets(32).build()?;
//! let results = trajs.estimate_all_orbits(&state, &mut rng, &params);
//! # Ok(()) }
//! ```
//!
//! ## See also
//! ------------
//! * [`TrajectoryFile`] – Public ingestion API surface.
//! * [`ObservationBatch`] – Zero-copy batch container (single observer).
//! * [`crate::trajectories::trajectory_fit::TrajectoryFit::estimate_all_orbits`] – Batch Gauss IOD.
//! * [`Outfit`] – Ephemerides, reference frames, and observer registry.
use std::{collections::HashMap, sync::Arc};

use super::batch_reader::observation_from_batch;
use crate::constants::ArcSec;
use crate::observers::Observer;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use crate::trajectories::batch_reader::ObservationBatch;
use crate::TrajectorySet;
use camino::Utf8Path;

use super::ades_reader::parse_ades;
use super::mpc_80col_reader::extract_80col;
use super::parquet_reader::parquet_to_trajset;

/// A trait for the TrajectorySet type definition.
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
pub trait TrajectoryFile {
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
    /// # use outfit::trajectories::batch_reader::ObservationBatch;
    /// # use outfit::TrajectorySet;
    /// # use outfit::TrajectoryFile;
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
    /// use outfit::trajectories::batch_reader::ObservationBatch;
    /// use outfit::TrajectorySet;
    /// use outfit::TrajectoryFile;
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
    /// * [`add_from_parquet`](crate::trajectories::trajectory_file::TrajectoryFile::add_from_parquet) – Adds observations from a Parquet file to an existing set.
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
    /// * [`new_from_parquet`](crate::trajectories::trajectory_file::TrajectoryFile::new_from_parquet) – Creates a brand new set from a Parquet file.
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
}

impl TrajectoryFile for TrajectorySet {
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
}
