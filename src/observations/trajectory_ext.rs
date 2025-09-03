use std::{collections::HashMap, sync::Arc};

use crate::constants::{ArcSec, Degree, ObjectNumber, Observations, TrajectorySet, MJD};
use crate::observations::observation_from_batch;
use crate::observers::Observer;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use camino::Utf8Path;

use super::ades_reader::parse_ades;
use super::extract_80col;
use super::parquet_reader::parquet_to_trajset;

pub struct ObservationBatch<'a> {
    pub ra: &'a [Degree],
    pub error_ra: ArcSec,
    pub dec: &'a [Degree],
    pub error_dec: ArcSec,
    pub time: &'a [MJD],
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

    /// Create a TrajectorySet from an object number, right ascension, declination, time, and one observer.
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `error_ra`: the error in right ascension (it is the same for all observations as it is the same observer)
    /// * `dec`: a vector of declination
    /// * `error_dec`: the error in declination (it is the same for all observations as it is the same observer)
    /// * `time`: a vector of time in MJD
    /// * `observer`: the observer
    ///
    /// Return
    /// ------
    /// * a TrajectorySet containing the observations
    fn new_from_vec(
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Self;

    /// Add the observations of an object number, right ascension, declination, time, and one observer to a TrajectorySet
    /// Each observations should have been observed by the same observer.
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: a mutable reference to the Outfit instance
    /// * `object_number`: the object number
    /// * `ra`: a vector of right ascension
    /// * `error_ra`: the error in right ascension (it is the same for all observations as it is the same observer)
    /// * `dec`: a vector of declination
    /// * `error_dec`: the error in declination (it is the same for all observations as it is the same observer)
    /// * `time`: a vector of time in MJD
    /// * `observer`: the observer
    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    );

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
}

impl TrajectoryExt for TrajectorySet {
    fn new_from_vec(
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Self {
        let observations: Observations = observation_from_batch(env_state, batch, observer);
        let mut traj_set: TrajectorySet = HashMap::default();
        traj_set.insert(ObjectNumber::String(object_number.into()), observations);
        traj_set
    }

    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) {
        let observations: Observations = observation_from_batch(env_state, batch, observer);
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
