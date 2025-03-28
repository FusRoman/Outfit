use ahash::{AHasher, RandomState};
use smallvec::SmallVec;
use std::{collections::HashMap, env, fs::File, iter, rc::Rc, sync::Arc};

use crate::observers::observers::Observer;
use crate::outfit::Outfit;
use crate::{
    constants::{Degree, MpcCode, ObjectNumber, Observations, TrajectorySet, JDTOMJD, MJD},
    observations::observations::Observation,
};
use arrow::{
    array::{Float64Array, PrimitiveArray, RecordBatch, StringArray, UInt32Array},
    datatypes::SchemaRef,
};
use camino::Utf8Path;
use hifitime::{efmt::format, Epoch};
use parquet::{
    arrow::{arrow_reader::ParquetRecordBatchReaderBuilder, ProjectionMask},
    file::reader::{FileReader, SerializedFileReader},
    record::{Row, RowAccessor},
};

use super::observations::{extract_80col, observation_from_vec};

/// A trait to extend the TrajectorySet struct
/// It provides methods to create a TrajectorySet from an 80 column file and to add observations from an 80 column file
/// It also provides methods to create a TrajectorySet from vectors of right ascension, declination, time, and observer and to add observations from vectors of right ascension, declination, time, and observer
pub trait TrajectoryExt {
    fn from_80col(env_state: &mut Outfit, colfile: &Utf8Path) -> Self;
    fn add_80col(&mut self, env_state: &mut Outfit, colfile: &Utf8Path);

    fn new_from_vec(
        env_state: &mut Outfit,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: Arc<Observer>,
    ) -> Self;
    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: Arc<Observer>,
    );

    fn new_from_parquet(
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        mpc_code: Arc<Observer>,
        batch_size: Option<usize>,
    ) -> Self;
    fn add_from_parquet(
        &mut self,
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        batch_size: Option<usize>,
    );
}

impl TrajectoryExt for TrajectorySet {
    /// Create a TrajectorySet from an object number, right ascension, declination, time, and one observer.
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `dec`: a vector of declination
    /// * `time`: a vector of time
    /// * `observer`: the observer
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations
    fn new_from_vec(
        env_state: &mut Outfit,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: Arc<Observer>,
    ) -> Self {
        let observations: Observations = observation_from_vec(env_state, ra, dec, time, observer);
        let mut traj_set: TrajectorySet = HashMap::default();
        traj_set.insert(ObjectNumber::String(object_number.into()), observations);
        traj_set
    }

    /// Add the observations of an object number, right ascension, declination, time, and one observer to a TrajectorySet
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `dec`: a vector of declination
    /// * `time`: a vector of time
    /// * `observer`: the observer
    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: Arc<Observer>,
    ) {
        let observations: Observations = observation_from_vec(env_state, ra, dec, time, observer);
        self.insert(
            ObjectNumber::String(object_number.to_string()),
            observations,
        );
    }

    fn add_from_parquet(
        &mut self,
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        batch_size: Option<usize>,
    ) {
        let new_trajs = Self::new_from_parquet(env_state, parquet, observer, batch_size);
        for (traj_id, obs) in new_trajs {
            self.entry(traj_id).or_default().extend(obs);
        }
    }

    fn new_from_parquet(
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        batch_size: Option<usize>,
    ) -> Self {
        let uint16_obs = env_state.uint16_from_observer(observer);

        let file = File::open(parquet).expect("Unable to open file");

        let builder =
            ParquetRecordBatchReaderBuilder::try_new(file).expect("Failed to read metadata");

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

        let mut trajectories: TrajectorySet = HashMap::default();
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
                    dec.expect("Expected DEC"),
                    mjd_time,
                );

                let traj_id = traj_id.expect("Expected TrajID");
                trajectories
                    .entry(ObjectNumber::Int(traj_id))
                    .or_insert_with(|| SmallVec::with_capacity(10))
                    .push(obs);
            }
        }

        trajectories
    }

    /// Create a TrajectorySet from an 80 column file
    ///
    /// Arguments
    /// ---------
    /// * `colfile`: a path to an 80 column file
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations from the 80 column file
    fn from_80col(env_state: &mut Outfit, colfile: &Utf8Path) -> Self {
        let mut traj_set: TrajectorySet = HashMap::default();
        let (observations, object_number) = extract_80col(env_state, colfile);
        traj_set.insert(object_number, observations);
        traj_set
    }

    /// Add the observations of a 80 column file to a TrajectorySet
    ///
    /// Arguments
    /// ---------
    /// * `colfile`: a path to an 80 column file
    fn add_80col(&mut self, env_state: &mut Outfit, colfile: &Utf8Path) {
        let (observations, object_number) = extract_80col(env_state, colfile);
        self.insert(object_number, observations);
    }
}
