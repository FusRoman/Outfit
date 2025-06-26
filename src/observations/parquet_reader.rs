use smallvec::SmallVec;
use std::{fs::File, sync::Arc};

use crate::constants::ArcSec;
use crate::observers::observers::Observer;
use crate::outfit::Outfit;
use crate::{
    constants::{ObjectNumber, TrajectorySet, JDTOMJD},
    observations::observations::Observation,
};
use arrow_array::array::{Float64Array, UInt32Array};
use camino::Utf8Path;
use parquet::arrow::{arrow_reader::ParquetRecordBatchReaderBuilder, ProjectionMask};

/// Reads a Parquet file and converts it to a TrajectorySet.
/// The Parquet file must contain the following columns: "ra", "dec", "jd", and "trajectory_id".
/// The "jd" column is converted to MJD using the JDTOMJD constant.
/// The "ra" and "dec" columns are expected to be in degrees.
///
/// Arguments
/// ---------
/// * `trajectories`: a mutable reference to a TrajectorySet
/// * `env_state`: a mutable reference to an Outfit
/// * `parquet`: a path to a Parquet file
/// * `observer`: an Arc<Observer>
/// * `error_ra`: the error in right ascension (RA) in arcseconds
/// * `error_dec`: the error in declination (DEC) in arcseconds
/// * `batch_size`: an optional batch size to use when reading the Parquet file, If None, the default batch size is 2048
///
/// Return
/// ------
/// * a TrajectorySet containing the observations from the Parquet file
pub(crate) fn parquet_to_trajset(
    trajectories: &mut TrajectorySet,
    env_state: &mut Outfit,
    parquet: &Utf8Path,
    observer: Arc<Observer>,
    error_ra: ArcSec,
    error_dec: ArcSec,
    batch_size: Option<usize>,
) {
    let uint16_obs = env_state.uint16_from_observer(observer);

    let file = File::open(parquet).expect("Unable to open file");

    let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect("Failed to read metadata");

    let parquet_metadata = builder.metadata();
    let schema_descr = parquet_metadata.file_metadata().schema_descr();

    let all_fields = schema_descr.columns();
    let column_names = ["ra", "dec", "jd", "trajectory_id"];
    let projection_indices: Vec<usize> = column_names
        .iter()
        .map(|name| {
            all_fields
                .iter()
                .position(|f| f.name() == *name)
                .unwrap_or_else(|| panic!("Column '{}' not found in schema", name))
        })
        .collect();

    let mask = ProjectionMask::leaves(schema_descr, projection_indices);

    let batch_size = batch_size.unwrap_or(2048);
    let reader = builder
        .with_projection(mask)
        .with_batch_size(batch_size)
        .build()
        .expect("Failed to build reader");

    let mut mjd_time = Vec::with_capacity(batch_size);

    for maybe_batch in reader {
        let batch = maybe_batch.expect("Error reading batch");

        let ra = batch
            .column_by_name("ra")
            .expect("Error getting ra column")
            .as_any()
            .downcast_ref::<Float64Array>()
            .expect("Error downcasting ra column");

        let dec = batch
            .column_by_name("dec")
            .expect("Error getting dec column")
            .as_any()
            .downcast_ref::<Float64Array>()
            .expect("Error downcasting dec column");

        let jd_time = batch
            .column_by_name("jd")
            .expect("Error getting jd column")
            .as_any()
            .downcast_ref::<Float64Array>()
            .expect("Error downcasting jd column");

        mjd_time.extend(
            jd_time
                .into_iter()
                .map(|jd| jd.expect("Expected JD") - JDTOMJD),
        );

        let traj_id = batch
            .column_by_name("trajectory_id")
            .expect("Error getting trajectory_id column")
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("Error downcasting trajectory_id column");

        for ((ra, dec), (mjd_time, traj_id)) in
            ra.into_iter().zip(dec).zip(mjd_time.drain(..).zip(traj_id))
        {
            let obs = Observation::new(
                uint16_obs,
                ra.expect("Expected RA"),
                error_ra,
                dec.expect("Expected DEC"),
                error_dec,
                mjd_time,
            );

            let traj_id = traj_id.expect("Expected TrajID");
            trajectories
                .entry(ObjectNumber::Int(traj_id))
                .or_insert_with(|| SmallVec::with_capacity(10))
                .push(obs);
        }
    }
}
