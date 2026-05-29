//! Bulk ephemeris results.
//!
//! The central type is [`EphemerisTable<T>`], which is the return type of all
//! bulk computation methods on [`OrbitalElements`].  Each entry pairs an input
//! [`hifitime::Epoch`] with either a successfully computed value of type `T` or
//! an [`OutfitError`] that describes why the computation failed for that
//! particular instant.
//!
//! # Concrete table types
//!
//! Three type aliases are provided for the most common payload types:
//!
//! | Alias | `T` | Produced by |
//! |---|---|---|
//! | [`ApparentPositionTable`] | [`ApparentPosition`] | `apparent_position_range / _at` |
//! | [`GeometryTable`] | [`BodyGeometry`] | `body_geometry_range / _at` |
//! | [`ApparentPositionAndGeometryTable`] | `(`[`ApparentPosition`]`, `[`BodyGeometry`]`)` | `apparent_position_and_geometry_range / _at` |
//!
//! # Design
//!
//! Errors are collected rather than short-circuited, so a single bad epoch
//! (e.g. one outside the JPL ephemeris coverage window) does not discard
//! results for all other epochs.  The caller can inspect failures via
//! [`EphemerisTable::error_count`] and iterate over successful entries with
//! [`EphemerisTable::successes`].
//!
//! Both `epochs` and `results` are guaranteed to have the same length and to
//! be in the same order as the input time sequence.
//!
//! # Iteration
//!
//! The table implements [`IntoIterator`], yielding
//! `(`[`Epoch`]`, `[`Result`]`<T, `[`OutfitError`]`>)` pairs:
//!
//! ```rust,ignore
//! for (epoch, result) in table {
//!     match result {
//!         Ok(pos) => println!("{epoch}: RA = {:.6}", pos.coord.ra),
//!         Err(e)  => eprintln!("{epoch}: {e}"),
//!     }
//! }
//! ```

use hifitime::Epoch;

use crate::OutfitError;

use super::{ApparentPosition, BodyGeometry};

// ---------------------------------------------------------------------------
// EphemerisTable<T>
// ---------------------------------------------------------------------------

/// A paired sequence of epochs and their computed ephemeris results.
///
/// The type parameter `T` determines what was computed at each epoch:
///
/// - [`ApparentPositionTable`] — sky coordinates + distances.
/// - [`GeometryTable`] — phase angle, elongation, radial velocity, angular rates.
/// - [`ApparentPositionAndGeometryTable`] — both of the above, from a **single**
///   orbit propagation per epoch.
///
/// Both `epochs` and `results` are always the same length and in the same
/// order as the input time sequence.
///
/// # Iteration
///
/// ```rust,ignore
/// for (epoch, result) in table {
///     if let Ok(value) = result {
///         println!("{epoch}: {value:?}");
///     }
/// }
/// ```
#[derive(Debug)]
pub struct EphemerisTable<T> {
    /// Observation epochs in the order they were requested.
    pub epochs: Vec<Epoch>,
    /// Per-epoch computation results, in the same order as `epochs`.
    pub results: Vec<Result<T, OutfitError>>,
}

impl<T> EphemerisTable<T> {
    /// Create an empty [`EphemerisTable`] with pre-allocated capacity.
    pub(crate) fn with_capacity(n: usize) -> Self {
        Self {
            epochs: Vec::with_capacity(n),
            results: Vec::with_capacity(n),
        }
    }

    /// Append a single epoch/result pair.
    pub(crate) fn push(&mut self, epoch: Epoch, result: Result<T, OutfitError>) {
        self.epochs.push(epoch);
        self.results.push(result);
    }

    // -----------------------------------------------------------------------
    // Public API
    // -----------------------------------------------------------------------

    /// Total number of epochs in the table (including failures).
    #[inline]
    pub fn len(&self) -> usize {
        self.epochs.len()
    }

    /// Returns `true` if the table contains no entries.
    #[inline]
    pub fn is_empty(&self) -> bool {
        self.epochs.is_empty()
    }

    /// Number of epochs for which computation succeeded.
    ///
    /// Equals `self.len() - self.error_count()`.
    pub fn success_count(&self) -> usize {
        self.results.iter().filter(|r| r.is_ok()).count()
    }

    /// Number of epochs for which computation failed.
    ///
    /// A value of `0` means every epoch was computed successfully.
    pub fn error_count(&self) -> usize {
        self.results.iter().filter(|r| r.is_err()).count()
    }

    /// Iterate over the epochs that were computed successfully.
    ///
    /// Yields `(&`[`Epoch`]`, &T)` pairs in input order, skipping any epoch
    /// whose computation returned an error.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// for (epoch, value) in table.successes() {
    ///     println!("{epoch}: {value:?}");
    /// }
    /// ```
    pub fn successes(&self) -> impl Iterator<Item = (&Epoch, &T)> {
        self.epochs
            .iter()
            .zip(self.results.iter())
            .filter_map(|(epoch, result)| result.as_ref().ok().map(|v| (epoch, v)))
    }

    /// Iterate over the epochs that failed, together with their errors.
    ///
    /// Yields `(&`[`Epoch`]`, &`[`OutfitError`]`)` pairs in input order.
    /// Returns an empty iterator when all epochs succeeded.
    pub fn errors(&self) -> impl Iterator<Item = (&Epoch, &OutfitError)> {
        self.epochs
            .iter()
            .zip(self.results.iter())
            .filter_map(|(epoch, result)| result.as_ref().err().map(|e| (epoch, e)))
    }

    /// Consume the table and return only the successfully computed entries.
    ///
    /// Drops all epochs that produced an error.  Useful when you need
    /// ownership of the results and are not interested in errors.
    ///
    /// # Example
    ///
    /// ```rust,ignore
    /// let successes: Vec<(Epoch, _)> = table.into_successes();
    /// ```
    pub fn into_successes(self) -> Vec<(Epoch, T)> {
        self.epochs
            .into_iter()
            .zip(self.results)
            .filter_map(|(epoch, result)| result.ok().map(|v| (epoch, v)))
            .collect()
    }
}

// ---------------------------------------------------------------------------
// IntoIterator
// ---------------------------------------------------------------------------

/// Consuming iterator over all epoch/result pairs in the table.
///
/// Yields `(`[`Epoch`]`, `[`Result`]`<T, `[`OutfitError`]`>)` in input order,
/// including both successes and failures.
impl<T> IntoIterator for EphemerisTable<T> {
    type Item = (Epoch, Result<T, OutfitError>);
    type IntoIter =
        std::iter::Zip<std::vec::IntoIter<Epoch>, std::vec::IntoIter<Result<T, OutfitError>>>;

    fn into_iter(self) -> Self::IntoIter {
        self.epochs.into_iter().zip(self.results)
    }
}

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

/// An [`EphemerisTable`] whose results are [`ApparentPosition`] values.
///
/// Produced by [`crate::OrbitalElements::apparent_position_range`] and
/// [`crate::OrbitalElements::apparent_position_at`].
pub type ApparentPositionTable = EphemerisTable<ApparentPosition>;

/// An [`EphemerisTable`] whose results are [`BodyGeometry`] values.
///
/// Produced by [`crate::OrbitalElements::body_geometry_range`] and
/// [`crate::OrbitalElements::body_geometry_at`].
pub type GeometryTable = EphemerisTable<BodyGeometry>;

/// An [`EphemerisTable`] whose results are `(`[`ApparentPosition`]`,
/// `[`BodyGeometry`]`)` tuples.
///
/// Both values are computed from a **single** orbit propagation per epoch,
/// making this more efficient than calling the position and geometry bulk
/// methods separately.
///
/// Produced by [`crate::OrbitalElements::apparent_position_and_geometry_range`] and
/// [`crate::OrbitalElements::apparent_position_and_geometry_at`].
pub type ApparentPositionAndGeometryTable = EphemerisTable<(ApparentPosition, BodyGeometry)>;
