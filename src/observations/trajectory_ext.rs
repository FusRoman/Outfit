use std::{collections::HashMap, fs::File, iter};

use crate::{
    constants::{Degree, MpcCode, ObjectNumber, Observations, TrajectorySet, MJD},
    observations::observations::Observation,
};
use camino::Utf8Path;
use hifitime::Epoch;
use parquet::{
    file::reader::{FileReader, SerializedFileReader},
    record::{Row, RowAccessor},
};

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

    fn new_from_parquet(parquet: &Utf8Path, mpc_code: MpcCode) -> Self;
    fn add_from_parquet(&mut self, parquet: &Utf8Path, observer: MpcCode);
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

    fn add_from_parquet(&mut self, parquet: &Utf8Path, observer: MpcCode) {
        let new_trajs = Self::new_from_parquet(parquet, observer);
        for (traj_id, obs) in new_trajs {
            self.entry(traj_id).or_default().extend(obs);
        }
    }

    fn new_from_parquet(parquet: &Utf8Path, observer: MpcCode) -> Self {
        let file = File::open(parquet)
            .expect(format!("Could not open file {}", parquet.as_str()).as_str());

        let reader = SerializedFileReader::new(file).expect(
            format!(
                "Could not create a SerializedFileReader from file {}",
                parquet.as_str()
            )
            .as_str(),
        );

        let mut trajectories: TrajectorySet = HashMap::new();

        let schema = reader.metadata().file_metadata().schema();
        let fields = schema.get_fields();

        let col_index: HashMap<&str, usize> = fields
            .iter()
            .enumerate()
            .map(|(i, field)| (field.name(), i))
            .collect();

        let iter = reader.get_row_iter(None).expect(
            format!(
                "Could not create a row iterator from file {}",
                parquet.as_str()
            )
            .as_str(),
        );

        let get_column = |row: &Row, column_name: &str| {
            row.get_double(
                *col_index.get(column_name).expect(
                    format!(
                        "The row does not contain the column {}, schema fields: {:?}",
                        column_name, fields
                    )
                    .as_str(),
                ),
            )
            .or_else(|_| {
                row.get_long(
                    *col_index.get(column_name).expect(
                        format!(
                            "The row does not contain the column {}, schema fields: {:?}",
                            column_name, fields
                        )
                        .as_str(),
                    ),
                )
                .map(|x| x as f64)
            })
            .expect(
                format!(
                    "Cannot parse the column value as a double or a long: {}",
                    column_name
                )
                .as_str(),
            )
        };

        for row in iter {
            let row = row.expect("Error reading row");

            let ra = get_column(&row, "ra");
            let dec = get_column(&row, "dec");

            let jd_time = get_column(&row, "jd");
            let mjd_time = Epoch::from_jde_utc(jd_time).to_mjd_utc_days();

            let traj_id = get_column(&row, "trajectory_id");

            let obs = Observation {
                observer: observer.clone(),
                ra: ra,
                dec: dec,
                time: mjd_time,
            };

            trajectories
                .entry(traj_id.to_string())
                .or_insert_with(Vec::new)
                .push(obs);
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
