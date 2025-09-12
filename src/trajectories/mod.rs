//! # Trajectories: ingestion, storage, and batch IOD
//!
//! High-level facilities to **ingest**, **store**, and **process** astrometric observations
//! grouped by object. The central type is [`TrajectorySet`], a fast hash map that buckets
//! time-ordered observations per [`ObjectNumber`]. Public helpers let you build a set from
//! multiple formats (MPC 80-column, Parquet, ADES, or in-memory batches) and run a
//! **Gauss-based Initial Orbit Determination (IOD)** over all objects.
//!
//! Modules
//! -----------------
//! * [`batch_reader`](crate::trajectories::batch_reader) – Zero-copy container and routines to expand single-observer batches
//!   into concrete [`Observation`](crate::observations::Observation)s.
//! * [`mpc_80col_reader`](crate::trajectories::mpc_80col_reader) – Minimal MPC **80-column** file reader.
//! * [`parquet_reader`](crate::trajectories::parquet_reader) – Arrow/Parquet-based ingestion (`ra`, `dec`, `jd`, `trajectory_id`).
//! * [`ades_reader`](crate::trajectories::ades_reader) – ADES (MPC XML/JSON) ingestion.
//! * [`trajectory_file`](crate::trajectories::trajectory_file) – **Public** trait exposing `new_from_*` and `add_from_*` helpers
//!   to construct/extend a [`TrajectorySet`] from the above sources.
//! * [`trajectory_fit`](crate::trajectories::trajectory_fit) – Batch Gauss IOD over a set (`TrajectoryFit` trait, results & stats).
//! * *(crate-private)* `progress_bar` – Optional progress UI when the `progress` feature is enabled.
//!
//! Data Model
//! -----------------
//! * **Key:** [`ObjectNumber`] (logical object identifier).
//! * **Value:** `Observations` = `SmallVec<Observation>` time-ordered per object.
//! * **Set:** [`TrajectorySet`] = `HashMap<ObjectNumber, Observations, ahash::RandomState>`
//!   for fast hashing and predictable performance on large catalogs.
//!
//! Ingestion Sources
//! -----------------
//! Use the [`trajectory_file::TrajectoryFile`](crate::trajectories::trajectory_file) trait (implemented for [`TrajectorySet`]):
//! * **80-col MPC** — `new_from_80col`, `add_from_80col` (fail-fast on parse errors).
//! * **Parquet** — `new_from_parquet`, `add_from_parquet` (propagate `OutfitError` on I/O/schema).
//! * **ADES** — `new_from_ades`, `add_from_ades` (error policy handled in the parser).
//! * **In-memory batch** (single observer) — `new_from_vec`, `add_from_vec`
//!   using [`batch_reader::ObservationBatch`](crate::trajectories::batch_reader::ObservationBatch) (angles/σ in **radians**, epochs in **MJD (TT)**).
//!
//! Units & Time Scales
//! -----------------
//! * Internal angles are **radians**; readers convert from **degrees/arcsec** as needed.
//! * Epochs are **MJD (TT)**; Parquet `"jd"` (assumed **TT**) is converted via
//!   [`constants::JDTOMJD`](crate::constants::JDTOMJD).
//! * Single-observer batches carry uniform 1-σ uncertainties (radians) applied per component.
//!
//! Batch IOD
//! -----------------
//! Use [`trajectory_fit::TrajectoryFit::estimate_all_orbits`](crate::trajectories::trajectory_fit::TrajectoryFit::estimate_all_orbits) to run Gauss IOD over the set.
//! Returns a map `ObjectNumber → Result<(GaussResult, rms), OutfitError>`. Errors are **per-object**
//! and do not abort other objects. A cooperative-cancel variant is also available.
//!
//! Performance Notes
//! -----------------
//! * Ingestion paths project only required columns and cache site positions by epoch.
//! * [`TrajectorySet`] uses `ahash` for speed; no deduplication is performed on `add_*` methods.
//! * Ordering is preserved as provided by sources; sorting by time is not enforced here.
//!
//! Feature Flags
//! -----------------
//! * `progress` — Enables a live progress bar and timing during batch IOD. See
//!   [`trajectory_fit`](crate::trajectories::trajectory_fit) for details. The UI is crate-internal and optional.
//!
//! Quick-Start
//! -----------------
//! ```rust,no_run
//! use rand::SeedableRng;
//! use std::sync::Arc;
//! use camino::Utf8Path;
//! use outfit::{Outfit, TrajectorySet};
//! use outfit::observers::Observer;
//! use outfit::trajectories::trajectory_file::TrajectoryFile;
//! use outfit::trajectories::trajectory_fit::TrajectoryFit;
//! use outfit::initial_orbit_determination::IODParams;
//!
//! # fn run() -> Result<(), outfit::outfit_errors::OutfitError> {
//! let mut state = Outfit::new("horizon:DE440", outfit::error_models::ErrorModel::FCCT14)?;
//! let observer: Arc<Observer> = state.get_observer_from_mpc_code(&"I41".into());
//!
//! // Ingest from Parquet, then append 80-column MPC
//! let mut trajs: TrajectorySet = TrajectorySet::new_from_parquet(
//!     &mut state, Utf8Path::new("obs.parquet"), observer.clone(), 0.5, 0.5, Some(8192)
//! )?;
//! trajs.add_from_80col(&mut state, Utf8Path::new("obs_80col.txt"));
//!
//! // Batch IOD
//! let mut rng = rand::rngs::StdRng::from_os_rng();
//! let params = IODParams::builder().max_triplets(32).build()?;
//! let results = trajs.estimate_all_orbits(&state, &mut rng, &params);
//! # Ok(()) }
//! ```
//!
//! See also
//! ------------
//! * [`trajectory_file`](crate::trajectories::trajectory_file) – Public ingestion API.
//! * [`batch_reader`](crate::trajectories::batch_reader) – Batch expansion for single-observer inputs.
//! * [`parquet_reader`](crate::trajectories::parquet_reader), [`mpc_80col_reader`](crate::trajectories::mpc_80col_reader), [`ades_reader`](crate::trajectories::ades_reader) – File readers.
//! * [`trajectory_fit`](crate::trajectories::trajectory_fit) – Batch IOD, results, and statistics.
//! * [`crate::observations::Observation`] – Atomic astrometric sample.
//!
//! ---
use std::collections::HashMap;

use ahash::RandomState;

use crate::{constants::Observations, ObjectNumber};

pub mod ades_reader;
pub mod batch_reader;
pub mod mpc_80col_reader;
pub mod parquet_reader;
pub mod trajectory_file;
pub mod trajectory_fit;

#[cfg(feature = "progress")]
pub(crate) mod progress_bar;

/// A full set of trajectories for multiple objects.
///
/// The key is the [`ObjectNumber`] (identifier of an object).
/// The value is the list of [`Observation`](crate::observations::Observation) associated with this object.
///
/// Uses [`ahash`](https://docs.rs/ahash) for fast hashing.
pub type TrajectorySet = HashMap<ObjectNumber, Observations, RandomState>;
