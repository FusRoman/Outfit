use std::{collections::HashMap, fs::File};

use crate::constants::{Degree, ObjectNumber, Observations, TrajectorySet, MJD};
use arrow::array::{Float64Array, RecordBatch};
use camino::Utf8Path;
use itertools::Itertools;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

use super::observations::{extract_80col, observation_from_vec};

/// A trait to extend the TrajectorySet struct
/// It provides methods to create a TrajectorySet from an 80 column file and to add observations from an 80 column file
/// It also provides methods to create a TrajectorySet from vectors of right ascension, declination, time, and observer and to add observations from vectors of right ascension, declination, time, and observer
pub trait TrajectoryExt {
    fn from_80col(colfile: &Utf8Path) -> Self;
    fn add_80col(&mut self, colfile: &Utf8Path);

    fn new_from_vec(
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) -> Self;
    fn add_from_vec(
        &mut self,
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    );

    fn new_from_parquet(parquet: &Utf8Path) -> Self;
    fn add_from_parquet(&mut self, parquet: &Utf8Path);
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
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) -> Self {
        let observations: Observations = observation_from_vec(ra, dec, time, observer);
        let mut traj_set: HashMap<ObjectNumber, Observations> = HashMap::new();
        traj_set.insert(object_number.to_string(), observations);
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
        object_number: &str,
        ra: &Vec<Degree>,
        dec: &Vec<Degree>,
        time: &Vec<MJD>,
        observer: &str,
    ) {
        let observations: Observations = observation_from_vec(ra, dec, time, observer);
        self.insert(object_number.to_string(), observations);
    }

    fn add_from_parquet(&mut self, parquet: &Utf8Path) {
        todo!()
    }

    fn new_from_parquet(parquet: &Utf8Path) -> Self {
        let file = File::open(parquet)
            .expect(format!("Could not open file {}", parquet.as_str()).as_str());
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).expect(
            format!(
                "Could not create a parquet reader for file {}",
                parquet.as_str()
            )
            .as_str(),
        );

        // get the fields of the schema
        let schema_field_names = builder
            .schema()
            .fields()
            .into_iter()
            .map(|f| f.name().clone())
            .collect_vec();

        // check if the required fields are present
        // ra, dec, jd, and trajectory_id
        assert!(
            schema_field_names.contains(&"ra".to_string())
                && schema_field_names.contains(&"dec".to_string())
                && schema_field_names.contains(&"jd".to_string())
                && schema_field_names.contains(&"trajectory_id".to_string()),
            "The parquet file does not contain the required columns ra, dec, jd, and trajectory_id, schema fields: {:?}",
            schema_field_names
        );

        let mut ra_vec: Vec<Degree> = Vec::new();
        let dec_vec: Vec<Degree> = Vec::new();
        let time_vec: Vec<MJD> = Vec::new();

        let mut reader = builder.build().unwrap();

        fn batch_col_to_arc_dyn<'a>(
            batch: &'a RecordBatch,
            col_name: &'a str,
            schema_name: &Vec<String>,
        ) -> &'a dyn arrow::array::Array {
            let col_index = batch.schema().index_of(col_name).expect(
                format!(
                    "The schema does not contain the column {}, schema fields: {:?}",
                    col_name, schema_name
                )
                .as_str(),
            );
            batch.column(col_index)
        }

        fn arc_dyn_to_vec(arc: &dyn arrow::array::Array) -> Vec<Degree> {
            arc.as_any()
                .downcast_ref::<Float64Array>()
                .expect("Error downcasting to Float64Array")
                .values()
                .iter()
                .map(|v| *v)
                .collect()
        }

        while let Some(batch) = reader.next() {
            let batch = batch.unwrap();

            let ra_col = batch_col_to_arc_dyn(&batch, &"ra", &schema_field_names);
            let dec_col = batch_col_to_arc_dyn(&batch, &"dec", &schema_field_names);
            let jd_col = batch_col_to_arc_dyn(&batch, &"jd", &schema_field_names);
            let trajectory_id_col =
                batch_col_to_arc_dyn(&batch, &"trajectory_id", &schema_field_names);

            let dec_col = batch.column_by_name("dec").expect(
                format!(
                    "The schema does not contain the column ra, schema fields: {:?}",
                    schema_field_names
                )
                .as_str(),
            );

            let tmp_ra: Vec<Degree> = arc_dyn_to_vec(ra_col);
            ra_vec.extend(tmp_ra);

            println!("{:?}", ra_vec);

            let col_index = batch.schema().index_of("ra").expect(
                format!(
                    "The schema does not contain the column ra, schema fields: {:?}",
                    schema_field_names
                )
                .as_str(),
            );
            let column = batch.column(col_index);
        }

        HashMap::new()
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
    fn from_80col(colfile: &Utf8Path) -> Self {
        let mut traj_set: HashMap<ObjectNumber, Observations> = HashMap::new();
        let (observations, object_number) = extract_80col(colfile);
        traj_set.insert(object_number, observations);
        traj_set
    }

    /// Add the observations of a 80 column file to a TrajectorySet
    ///
    /// Arguments
    /// ---------
    /// * `colfile`: a path to an 80 column file
    fn add_80col(&mut self, colfile: &Utf8Path) {
        let (observations, object_number) = extract_80col(colfile);
        self.insert(object_number, observations);
    }
}