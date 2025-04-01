use std::{collections::HashMap, sync::Arc};

use crate::constants::{Degree, ObjectNumber, Observations, TrajectorySet, MJD};
use crate::observers::observers::Observer;
use crate::outfit::Outfit;
use camino::Utf8Path;

use super::ades_reader::parse_ades;
use super::observations::{extract_80col, observation_from_vec};
use super::parquet_reader::parquet_to_trajset;

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

    fn new_from_ades(env_state: &mut Outfit, ades: &Utf8Path) -> Self;
    fn add_from_ades(&mut self, env_state: &mut Outfit, ades: &Utf8Path);
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

    /// Add a set of trajectories from a parquet file to a TrajectorySet
    ///
    /// Arguments
    /// ---------
    /// * `parquet`: a path to a parquet file
    /// * `observer`: the observer
    /// * `batch_size`: the batch size to use when reading the parquet file, if None, the default batch size is 2048
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the new observations is added to the existing TrajectorySet
    ///
    /// Note: the parquet file should contain the columns "ra", "dec", "jd", and "trajectory_id"
    /// The "jd" column is converted to MJD using the JDTOMJD constant
    /// The "ra" and "dec" columns are converted to 32 bits for performance
    fn add_from_parquet(
        &mut self,
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        batch_size: Option<usize>,
    ) {
        parquet_to_trajset(self, env_state, parquet, observer, batch_size);
    }

    /// Create a TrajectorySet from a parquet file
    ///
    /// Arguments
    /// ---------
    /// * `parquet`: a path to a parquet file
    /// * `observer`: the observer
    /// * `batch_size`: the batch size to use when reading the parquet file, if None, the default batch size is 2048
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations from the parquet file
    ///
    /// Note: the parquet file should contain the columns "ra", "dec", "jd", and "trajectory_id"
    /// The "jd" column is converted to MJD using the JDTOMJD constant
    /// The "ra" and "dec" columns are converted to 32 bits for performance
    fn new_from_parquet(
        env_state: &mut Outfit,
        parquet: &Utf8Path,
        observer: Arc<Observer>,
        batch_size: Option<usize>,
    ) -> Self {
        let mut trajs: TrajectorySet = HashMap::default();
        parquet_to_trajset(&mut trajs, env_state, parquet, observer, batch_size);
        trajs
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

    fn add_from_ades(&mut self, env_state: &mut Outfit, ades: &Utf8Path) {
        parse_ades(env_state, ades, self);
    }

    fn new_from_ades(env_state: &mut Outfit, ades: &Utf8Path) -> Self {
        let mut trajs: TrajectorySet = HashMap::default();
        parse_ades(env_state, ades, &mut trajs);
        trajs
    }
}
