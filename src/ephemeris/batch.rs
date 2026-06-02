//! Batch ephemeris computation over a [`FullOrbitResult`].
//!
//! This module exposes the [`FullOrbitResultExt`] extension trait, which adds
//! ephemeris generation methods directly to [`FullOrbitResult`].
//!
//! # Sequential vs parallel
//!
//! Two methods are available:
//!
//! | Method | Feature gate | Description |
//! |---|---|---|
//! | [`compute_ephemerides`](FullOrbitResultExt::compute_ephemerides) | *(none)* | Sequential iteration over all orbits |
//! | [`compute_ephemerides_parallel`](FullOrbitResultExt::compute_ephemerides_parallel) | `parallel` | Rayon-parallel iteration over all orbits |
//!
//! # Return type
//!
//! Both methods return a
//! `HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState>`:
//!
//! - `Ok(EphemerisResult<…>)` — the orbit was successfully determined **and**
//!   ephemeris computation was requested. Note that individual
//!   `(epoch, observer)` pairs inside the [`EphemerisResult`] may still carry
//!   per-entry errors.
//! - `Err(OutfitError)` — the orbit determination itself failed for this
//!   trajectory; no ephemeris can be produced.
//!
//! # Example
//!
//! ```rust,ignore
//! use outfit::{FullOrbitResultExt, EphemerisRequest, EphemerisConfig, Combined};
//!
//! // `full_orbit_result` is a FullOrbitResult from fit_full_iod / fit_full_lsq
//! let ephemerides = full_orbit_result.compute_ephemerides(
//!     &EphemerisRequest::<Combined>::new(EphemerisConfig::default())
//!         .add(observer, EphemerisMode::Range { start, end, step }),
//!     &jpl,
//!     &ut1,
//! );
//!
//! for (traj_id, result) in &ephemerides {
//!     match result {
//!         Ok(ephem) => println!("{traj_id:?}: {} entries", ephem.len()),
//!         Err(e)    => eprintln!("{traj_id:?}: orbit error — {e}"),
//!     }
//! }
//! ```

use std::collections::HashMap;

use ahash::RandomState;
use photom::TrajId;

#[cfg(feature = "parallel")]
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    constants::FullOrbitResult,
    ephemeris::{EphemerisOutputKind, EphemerisRequest, EphemerisResult},
    JPLEphem, OutfitError,
};
use hifitime::ut1::Ut1Provider;

// ---------------------------------------------------------------------------
// FullOrbitResultExt
// ---------------------------------------------------------------------------

/// Extension trait that adds batch ephemeris generation to [`FullOrbitResult`].
///
/// Import this trait to call
/// [`compute_ephemerides`](Self::compute_ephemerides) (and, with the
/// `parallel` feature, [`compute_ephemerides_parallel`](Self::compute_ephemerides_parallel))
/// on any [`FullOrbitResult`] value.
pub trait FullOrbitResultExt {
    /// Compute ephemerides for every orbit in the map, sequentially.
    ///
    /// Iterates over all `(traj_id, orbit_result)` pairs:
    ///
    /// - If the orbit determination **succeeded** (`Ok`), the orbital elements
    ///   are used to evaluate the full [`EphemerisRequest`] and the result is
    ///   stored as `Ok(EphemerisResult<…>)`.
    /// - If the orbit determination **failed** (`Err`), the error is forwarded
    ///   as-is (converted to a string and re-wrapped) so the caller always
    ///   gets a complete map keyed by every [`TrajId`] in the input.
    ///
    /// # Arguments
    ///
    /// - `request` — typed ephemeris request (observers, modes, config).
    /// - `jpl`     — JPL planetary ephemeris.
    /// - `ut1`     — UT1 time-scale provider.
    ///
    /// # Returns
    ///
    /// A `HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState>`
    /// with the same key set as `self`.
    fn compute_ephemerides<O: EphemerisOutputKind>(
        &self,
        request: &EphemerisRequest<O>,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
    ) -> HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState>;

    /// Compute ephemerides for every orbit in the map, **in parallel**.
    ///
    /// Identical to [`compute_ephemerides`](Self::compute_ephemerides) but
    /// uses Rayon to process trajectories concurrently.  Each trajectory is
    /// independent, so there are no ordering guarantees on the map entries.
    ///
    /// Enabled only with the `parallel` feature flag.
    ///
    /// # Arguments
    ///
    /// Same as [`compute_ephemerides`](Self::compute_ephemerides).
    ///
    /// # Returns
    ///
    /// Same type as [`compute_ephemerides`](Self::compute_ephemerides).
    #[cfg(feature = "parallel")]
    fn compute_ephemerides_parallel<O>(
        &self,
        request: &EphemerisRequest<O>,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
    ) -> HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState>
    where
        O: EphemerisOutputKind + Send + Sync,
        O::Output: Send;
}

// ---------------------------------------------------------------------------
// impl for FullOrbitResult
// ---------------------------------------------------------------------------

impl FullOrbitResultExt for FullOrbitResult {
    fn compute_ephemerides<O: EphemerisOutputKind>(
        &self,
        request: &EphemerisRequest<O>,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
    ) -> HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState> {
        let mut map = HashMap::with_capacity_and_hasher(self.len(), RandomState::new());

        for (traj_id, orbit_result) in self {
            let entry = match orbit_result {
                Ok(fit) => Ok(fit.orbital_elements().compute(request, jpl, ut1)),
                Err(e) => Err(OutfitError::InvalidConversion(e.to_string())),
            };
            map.insert(traj_id.clone(), entry);
        }

        map
    }

    #[cfg(feature = "parallel")]
    fn compute_ephemerides_parallel<O>(
        &self,
        request: &EphemerisRequest<O>,
        jpl: &JPLEphem,
        ut1: &Ut1Provider,
    ) -> HashMap<TrajId, Result<EphemerisResult<O::Output>, OutfitError>, RandomState>
    where
        O: EphemerisOutputKind + Send + Sync,
        O::Output: Send,
    {
        let new_map = || HashMap::with_hasher(RandomState::new());

        self.par_iter()
            .map(|(traj_id, orbit_result)| {
                let entry = match orbit_result {
                    Ok(fit) => Ok(fit.orbital_elements().compute(request, jpl, ut1)),
                    Err(e) => Err(OutfitError::InvalidConversion(e.to_string())),
                };
                (traj_id.clone(), entry)
            })
            .fold(new_map, |mut map, (k, v)| {
                map.insert(k, v);
                map
            })
            .reduce(new_map, |mut a, b| {
                a.extend(b);
                a
            })
    }
}
