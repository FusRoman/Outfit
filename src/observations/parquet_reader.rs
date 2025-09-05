use arrow_array::Array;
use hifitime::Epoch;
use nalgebra::Vector3;
use ordered_float::OrderedFloat;
use parquet::errors::ParquetError;
use smallvec::SmallVec;
use std::collections::hash_map::Entry;
use std::io;
use std::sync::Arc;

use crate::constants::{ArcSec, RADEG};
use crate::observers::Observer;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use crate::{
    constants::{ObjectNumber, TrajectorySet, JDTOMJD},
    observations::Observation,
};
use arrow_array::array::{Float64Array, UInt32Array};
use camino::Utf8Path;
use parquet::arrow::{arrow_reader::ParquetRecordBatchReaderBuilder, ProjectionMask};

use ahash::RandomState;
use std::collections::HashMap;

type FastHashMap<K, V> = HashMap<K, V, RandomState>;

/// Load astrometric observations from a Parquet file into an existing [`TrajectorySet`].
///
/// This routine deserializes batches of observations from a Parquet file, converts
/// Julian Dates (JD) to Modified Julian Dates (MJD; TT scale), and constructs [`Observation`]
/// instances with the provided astrometric uncertainties. To avoid redundant and costly
/// ephemeris computations, it caches the observer's geocentric and heliocentric positions
/// per unique `(observer, time)` encountered during the read.
///
/// Arguments
/// -----------------
/// * `trajectories` – The mutable [`TrajectorySet`] to which observations are appended.
/// * `env_state` – Global environment providing ephemerides, Earth orientation data,
///   and observer definitions.
/// * `parquet` – Path to the input Parquet file containing the columns
///   `ra`, `dec`, `jd`, and `trajectory_id`.
///   The `ra` and `dec` columns have to be in degrees and of type `Float64`.
///   The `jd` column has to be in Julian Date (TT) and of type `Float64`.
///   The `trajectory_id` column has to be of type `UInt32` and is used to group
///   observations by object.
/// * `observer` – The [`Observer`] associated with all observations in this file.
/// * `error_ra` – Right ascension astrometric uncertainty (radians).
/// * `error_dec` – Declination astrometric uncertainty (radians).
/// * `batch_size` – Optional Arrow reader batch size (default: 8192 rows).
///
/// Performance notes
/// -----------------
/// * Projects only the required columns and accesses them by **index** to avoid
///   per-batch name lookups.
/// * Uses a per-file cache keyed by MJD (TT) to store `(geo_pos, helio_pos)` and
///   reuse them across batches and rows.
/// * Fast path when all columns are non-null: iterates over raw slices (`&[f64]`, `&[u32]`)
///   for minimal overhead.
/// * Downcasts to concrete Arrow arrays once per batch (not per row).
///
/// Return
/// ----------
/// * No return value. Observations are appended in-place to `trajectories`.
///
/// See also
/// ------------
/// * [`Observation::with_positions`] – Zero-ephemeris constructor used here.
/// * [`Observer::pvobs`] – Geocentric observer position (and velocity).
/// * [`Observer::helio_position`] – Heliocentric observer position.
pub(crate) fn parquet_to_trajset(
    trajectories: &mut TrajectorySet,
    env_state: &mut Outfit,
    parquet: &Utf8Path,
    observer: Arc<Observer>,
    error_ra: ArcSec,
    error_dec: ArcSec,
    batch_size: Option<usize>,
) -> Result<(), OutfitError> {
    // Resolve observer to its compact u16 key once (hot path avoids map lookups later).
    let uint16_obs = env_state.uint16_from_observer(observer);

    // Open file and inspect Parquet metadata (I/O and schema discovery happen here).
    let file = std::fs::File::open(parquet)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;

    let parquet_metadata = builder.metadata();
    let schema_descr = parquet_metadata.file_metadata().schema_descr();

    // Build a stable projection mask to materialize **only** the columns used by Outfit.
    // We rely on the projection order to index columns directly in the hot loop.
    let all_fields = schema_descr.columns();
    let column_names = ["ra", "dec", "jd", "trajectory_id"];
    let projection_indices: Vec<usize> = column_names
        .iter()
        .map(|name| {
            all_fields
                .iter()
                .position(|f| f.name() == *name)
                // If not found, surface a clean error instead of panicking.
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::NotFound,
                        format!("Column '{name}' not found in schema"),
                    )
                })
        })
        .collect::<Result<_, _>>()?;
    let mask = ProjectionMask::leaves(schema_descr, projection_indices);

    // A larger batch often amortizes decompression + Arrow deserialization cost.
    // Tune with benches (8192–65536) depending on your I/O and CPU characteristics.
    let batch_size = batch_size.unwrap_or(8192);
    let reader = builder
        .with_projection(mask)
        .with_batch_size(batch_size)
        .build()?;

    // Pre-fetch shared providers and observer reference to avoid repeated lookups in the hot loop.
    let ut1 = env_state.get_ut1_provider();
    let obs_ref = env_state.get_observer_from_uint16(uint16_obs);

    // Cache MJD(TT) → (geo_pos, helio_pos).
    // Note: we assume here the file contains a **single** observer. If you ingest multiple
    // observers, extend the key to (observer_id, time), e.g. by packing into a u64 or using a tuple.
    let mut pos_cache: FastHashMap<OrderedFloat<f64>, (Vector3<f64>, Vector3<f64>)> =
        FastHashMap::with_capacity_and_hasher(4096, RandomState::default());

    // Iterate over Parquet record batches
    for maybe_batch in reader {
        // I/O boundary; failures here usually indicate corruption or incompatible schema.
        let batch = maybe_batch.map_err(ParquetError::from)?;
        let len = batch.num_rows();

        // Projected columns by index: [0]=ra, [1]=dec, [2]=jd, [3]=trajectory_id
        // We downcast **once** per batch (cheap) and reuse typed views in the row loop.
        let ra_arr = batch
            .column(0)
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "ra must be Float64Array"))?;
        let dec_arr = batch
            .column(1)
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidData, "dec must be Float64Array")
            })?;
        let jd_arr = batch
            .column(2)
            .as_any()
            .downcast_ref::<Float64Array>()
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "jd must be Float64Array"))?;
        let tid_arr = batch
            .column(3)
            .as_any()
            .downcast_ref::<UInt32Array>()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "trajectory_id must be UInt32Array",
                )
            })?;

        // Fast path when all projected columns have no nulls.
        // This unlocks tight, bounds-checked loops on `&[T]` slices without per-row Option unwrapping.
        let no_nulls = ra_arr.nulls().is_none()
            && dec_arr.nulls().is_none()
            && jd_arr.nulls().is_none()
            && tid_arr.nulls().is_none();

        if no_nulls {
            // Raw slice views (zero allocs, no per-element downcast/boxing).
            let ra_vals: &[f64] = ra_arr.values();
            let dec_vals: &[f64] = dec_arr.values();
            let jd_vals: &[f64] = jd_arr.values();
            let tid_vals: &[u32] = tid_arr.values();

            // Hot loop: build `Observation`s with positions sourced from the per-file cache.
            for i in 0..len {
                // Convert degrees → radians (fast path: one multiply each).
                let ra_rad = ra_vals[i] * RADEG;
                let dec_rad = dec_vals[i] * RADEG;

                // Convert JD → MJD(TT). Assumes `jd` is already in TT scale in the file;
                // if not, convert scales here before caching.
                let mjd_time = jd_vals[i] - JDTOMJD;
                // Using OrderedFloat to allow float keys in a hash map (with a total order).
                let key = OrderedFloat(mjd_time);

                // Compute observer positions once per unique epoch.
                // We handle errors with `?` while using the entry API (no closures returning Result).
                let (geo_pos, helio_pos) = match pos_cache.entry(key) {
                    Entry::Occupied(e) => *e.get(),
                    Entry::Vacant(v) => {
                        let epoch =
                            Epoch::from_mjd_in_time_scale(mjd_time, hifitime::TimeScale::TT);

                        // Geocentric position (velocity unused here, but available if needed).
                        let (geo, _vel) = obs_ref.pvobs(&epoch, ut1)?;
                        // Heliocentric position (requires geocentric as input).
                        let helio = obs_ref.helio_position(env_state, &epoch, &geo)?;

                        v.insert((geo, helio));
                        (geo, helio)
                    }
                };

                // Zero-ephemeris constructor: avoids recomputing positions at construction.
                let obs = Observation::with_positions(
                    uint16_obs, ra_rad, error_ra, dec_rad, error_dec, mjd_time, geo_pos, helio_pos,
                );

                // Group observations by `trajectory_id` (ObjectNumber::Int).
                let obj = ObjectNumber::Int(tid_vals[i]);
                trajectories
                    .entry(obj)
                    .or_insert_with(|| SmallVec::with_capacity(32))
                    .push(obs);
            }
        } else {
            // Safety fallback: if any column contains nulls, we check row-by-row and skip incomplete rows.
            // This path is slower, but maintains correctness for sparse/missing data.
            for i in 0..len {
                if ra_arr.is_null(i)
                    || dec_arr.is_null(i)
                    || jd_arr.is_null(i)
                    || tid_arr.is_null(i)
                {
                    continue; // Drop incomplete rows (policy: skip; alternatively, surface an error)
                }

                let ra_rad: f64 = ra_arr.value(i) * RADEG;
                let dec_rad = dec_arr.value(i) * RADEG;
                let mjd_time = jd_arr.value(i) - JDTOMJD;
                let tid = tid_arr.value(i);
                let key = OrderedFloat(mjd_time);

                let (geo_pos, helio_pos) = match pos_cache.entry(key) {
                    Entry::Occupied(e) => *e.get(),
                    Entry::Vacant(v) => {
                        let epoch =
                            Epoch::from_mjd_in_time_scale(mjd_time, hifitime::TimeScale::TT);

                        let (geo, _vel) = obs_ref.pvobs(&epoch, ut1)?;
                        let helio = obs_ref.helio_position(env_state, &epoch, &geo)?;

                        v.insert((geo, helio));
                        (geo, helio)
                    }
                };

                let obs = Observation::with_positions(
                    uint16_obs, ra_rad, error_ra, dec_rad, error_dec, mjd_time, geo_pos, helio_pos,
                );

                trajectories
                    .entry(ObjectNumber::Int(tid))
                    .or_insert_with(|| SmallVec::with_capacity(32))
                    .push(obs);
            }
        }
    }

    Ok(())
}
