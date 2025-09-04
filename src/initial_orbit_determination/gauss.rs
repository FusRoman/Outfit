//! # Gauss Method for Initial Orbit Determination
//!
//! This module provides an implementation of the **classical Gauss method** for
//! initial orbit determination (IOD) using three optical astrometric observations
//! acquired from the same observing site.
//!
//! ## Core structure: [`GaussObs`]
//!
//! The [`GaussObs`] struct encapsulates all the data required to apply the Gauss algorithm:
//!
//! * Indices of the three selected observations,
//! * Right ascension and declination angles `[rad]`,
//! * Observation epochs (MJD, TT scale),
//! * Observer heliocentric positions at each epoch (in AU, equatorial mean J2000).
//!
//! It serves as the central container for preparing the geometry, solving the
//! Gauss polynomial, and recovering a preliminary orbit.
//!
//! ## Main functionalities
//!
//! - **Construction**
//!   * [`GaussObs::with_observer_position`] – Build from RA/DEC/epochs with explicit observer positions.
//!
//! - **Geometry preparation**
//!   * Assemble the `3×3` matrix of line-of-sight unit vectors,
//!   * Build the time-scaled coefficients and the `a`, `b` vectors required by the Gauss polynomial.
//!
//! - **Orbit determination**
//!   * Solve the 8th-degree Gauss polynomial for the central topocentric distance,
//!   * Reconstruct heliocentric positions,
//!   * Estimate velocity via **Gibbs’ method**,
//!   * Return a [`GaussResult`] containing preliminary [`OrbitalElements`] (not necessarily Keplerian).
//!
//! - **Orbit refinement**
//!   * [`GaussObs::pos_and_vel_correction`] – Iterative update of positions and velocities,
//!   * Apply filters on eccentricity and perihelion distance via [`eccentricity_control`].
//!
//! - **Monte Carlo perturbations**
//!   * [`GaussObs::generate_noisy_realizations`] – Generate perturbed triplets by adding Gaussian noise
//!     to RA/DEC, useful for uncertainty propagation and robustness analysis.
//!
//! ## Algorithm outline
//!
//! 1. Compute unit direction vectors from the three RA/DEC observations.
//! 2. Invert the direction matrix to solve the linear system for topocentric distances.
//! 3. Solve the **8th-degree polynomial** for the central distance `ρ₂` (Aberth’s method).
//! 4. Reconstruct heliocentric positions and estimate velocity via **Gibbs’ method**.
//! 5. Convert position/velocity to **orbital elements** (which may be Keplerian, Equinoctial, or Cometary).
//! 6. Optionally apply iterative corrections.
//!
//! ## Output
//!
//! The primary result is a [`GaussResult`] which stores the computed orbit as
//! a generic [`OrbitalElements`] value.  
//! Depending on refinement, this may represent either:
//! * a **Preliminary orbit** (direct Gauss solution), or
//! * a **Corrected orbit** (after iteration).
//!
//! ## Example
//!
//! ```rust, no_run
//! use outfit::initial_orbit_determination::gauss::GaussObs;
//! use outfit::orbit_type::{OrbitalElements, keplerian_element::KeplerianElements};
//! use nalgebra::Vector3;
//! use outfit::outfit::Outfit;
//! use outfit::error_models::ErrorModel;
//! use outfit::initial_orbit_determination::gauss_result::GaussResult;
//! use outfit::initial_orbit_determination::IODParams;
//!
//! let env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
//!
//! // Build GaussObs (here positions are assumed precomputed)
//! let gauss: GaussObs = unimplemented!("Construct GaussObs from RA/DEC/time and observer state");
//!
//! // Run the preliminary orbit computation
//! let result = gauss.prelim_orbit(&env, &IODParams::default()).unwrap();
//!
//! // Match on the returned orbital-element representation
//! match result {
//!     GaussResult::PrelimOrbit(oe)
//!     | GaussResult::CorrectedOrbit(oe) => {
//!         match oe {
//!             OrbitalElements::Keplerian(kepl) => {
//!                 println!("a [AU] = {}", kepl.semi_major_axis);
//!             }
//!             OrbitalElements::Equinoctial(eq) => {
//!                 println!("lambda [rad] = {}", eq.mean_longitude);
//!             }
//!             OrbitalElements::Cometary(com) => {
//!                 println!("q [AU] = {}", com.perihelion_distance);
//!             }
//!         }
//!     }
//! }
//! ```
//!
//! ## References
//!
//! * Milani & Gronchi (2010) – *Theory of Orbit Determination*
//!
//! ## See also
//!
//! - [`GaussObs`]
//! - [`GaussResult`]
//! - [`OrbitalElements`]
//! - [`KeplerianElements`](crate::orbit_type::keplerian_element::KeplerianElements)
use std::ops::ControlFlow;

use aberth::StopReason;
use nalgebra::Matrix3;
use nalgebra::Vector3;
use rand_distr::StandardNormal;
use smallvec::SmallVec;

use crate::constants::Radian;
use crate::constants::{GAUSS_GRAV, VLIGHT_AU};

use crate::initial_orbit_determination::gauss_result::GaussResult;
use crate::initial_orbit_determination::IODParams;
use crate::kepler::velocity_correction_with_guess;
use crate::orb_elem::eccentricity_control;
use crate::orbit_type::OrbitalElements;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use aberth::aberth;
use rand::Rng;

/// Observation triplet structure for Gauss's initial orbit determination (IOD).
///
/// This struct encapsulates the minimal set of information needed to apply the
/// classical **Gauss method** for orbit determination, based on three optical
/// astrometric observations taken from the same observatory site.
///
/// Fields
/// -----------------
/// * `idx_obs`: Indices of the three observations in the full dataset (useful for traceability).
/// * `ra`: Right ascensions `[rad]` of the three observations, expressed in the ICRS frame.
/// * `dec`: Declinations `[rad]` of the three observations, expressed in the ICRS frame.
/// * `time`: Observation epochs in Modified Julian Date (TT scale).
/// * `observer_helio_position`: `3×3` matrix of heliocentric observer position vectors, with:
///   - **columns** = observation epochs (1 column per observation),
///   - **units** = astronomical units (AU),
///   - **frame** = equatorial mean J2000 (aligned with ICRS).
///
/// Return
/// ----------
/// * A lightweight container [`GaussObs`] holding triplet data and observer positions,
///   ready to be passed to the Gauss IOD solver.
///
/// See also
/// -------------
/// * [`GaussObs::with_observer_position`] – Constructor from RA/DEC/epochs with explicit observer positions.
#[derive(Debug, PartialEq, Clone)]
pub struct GaussObs {
    pub(crate) idx_obs: Vector3<usize>,
    pub(crate) ra: Vector3<Radian>,
    pub(crate) dec: Vector3<Radian>,
    pub(crate) time: Vector3<f64>,
    pub(crate) observer_helio_position: Matrix3<f64>,
}

type GaussPrelimResult = Result<
    (
        f64,
        f64,
        Matrix3<f64>,
        Matrix3<f64>,
        Vector3<f64>,
        Vector3<f64>,
    ),
    OutfitError,
>;

/// Compute Descartes’ sign-variation upper bound for a sparse degree-8 polynomial.
/// Specialized for:
/// `p(x) = c0 + c3·x^3 + c6·x^6 + 1·x^8` (monic leading term).
///
/// This routine applies **Descartes’ rule of signs** to obtain an *upper bound* on the
/// number of **positive real roots** of `p`. For positive roots, the rule counts the
/// sign changes in the coefficient sequence ordered by **decreasing degree**, after
/// removing coefficients that are exactly (or effectively) zero.  Because the polynomial
/// is sparse and monic, the sequence reduces to just four terms:
/// `[ +1 (for x^8), sign(c6), sign(c3), sign(c0) ]`.
///
/// The bound returned is exact **only as an upper bound**: the true number of positive
/// real roots equals this bound minus an **even, non-negative integer**.  In particular,
/// a result of `0` guarantees that **no** positive real root exists and you can safely
/// skip any numerical root-finding for `x>0`.
///
/// Mathematical context
/// --------------------
/// * **Descartes’ rule of signs:** Let `V` be the number of sign variations in the ordered
///   coefficient sequence of `p(x)` (descending degree, zeros removed). Then the number
///   of positive real roots `N₊` (counted with multiplicities) satisfies:  
///   `N₊ ≤ V` and `V − N₊` is an even integer.
/// * The bound is **invariant under positive scaling** of the polynomial (multiplying all
///   coefficients by a positive constant does not change the sign pattern).
/// * Handling near-zero coefficients (`|c| ≤ zero_eps`) improves robustness by avoiding
///   spurious sign flips caused by numerical noise.
///
/// Arguments
/// -----------------
/// * `c0`, `c3`, `c6`: Coefficients of the sparse polynomial `p(x) = c0 + c3·x^3 + c6·x^6 + x^8`.
/// * `zero_eps`: Absolute threshold; coefficients with `|c| ≤ zero_eps` are treated as zero
///   (they are ignored when counting sign variations).
///
/// Return
/// ----------
/// * An upper bound (`u32`) on the number of **positive real roots** of `p(x)`.
///
/// Notes
/// ----------
/// * Time complexity is **O(1)** (constant time): only four signs are inspected.
/// * A return value of `0` implies “no positive real root” (useful as a fast pre-filter
///   before calling a complex root solver such as Aberth).
#[inline]
fn descartes_upper_bound_deg8_sparse(c0: f64, c3: f64, c6: f64, zero_eps: f64) -> u32 {
    #[inline]
    fn s(v: f64, eps: f64) -> i8 {
        if v.abs() <= eps {
            0
        } else if v.is_sign_positive() {
            1
        } else {
            -1
        }
    }
    // Leading coefficient is +1 for x^8.
    let seq = [1_i8, s(c6, zero_eps), s(c3, zero_eps), s(c0, zero_eps)];

    let mut last = 0_i8;
    let mut count = 0_u32;
    for &cur in &seq {
        if cur == 0 {
            continue;
        }
        if last != 0 && cur != last {
            count += 1;
        }
        last = cur;
    }
    count
}

impl GaussObs {
    /// Construct a `GaussObs` object from an observation triplet and corresponding observer positions.
    ///
    /// This constructor builds a [`GaussObs`] instance directly from pre-computed inputs:  
    /// indices of the selected observations, their right ascension/declination, the observation epochs,
    /// and the heliocentric positions of the observer at each epoch.
    ///
    /// Arguments
    /// -----------------
    /// * `idx_obs`: Triplet of indices identifying the three observations used for initial orbit determination.
    /// * `ra`: Right ascension values `[rad]` for the three observations.
    /// * `dec`: Declination values `[rad]` for the three observations.
    /// * `mjd_time`: Observation times in Modified Julian Date (TT scale).
    /// * `observer_helio_position`: Matrix `3×3` of heliocentric observer positions, each column corresponding
    ///   to the observer position (in AU) at the epoch of the matching observation.
    ///
    /// Return
    /// ----------
    /// * A `Result` containing a [`GaussObs`] struct populated with:
    ///   - the input indices, RA/DEC, and epochs,
    ///   - the provided heliocentric observer positions.
    ///
    /// Remarks
    /// ------------
    /// * The reference frame is the equatorial mean J2000 system (AU units).
    /// * Element `i` in `ra`, `dec`, `mjd_time`, and the `i`-th column of `observer_helio_position`
    ///   all correspond to the same astrometric observation.
    ///
    /// See also
    /// ------------
    /// * [`GaussObs`] – Data structure for storing triplet-based IOD input.
    pub fn with_observer_position(
        idx_obs: Vector3<usize>,
        ra: Vector3<f64>,
        dec: Vector3<f64>,
        mjd_time: Vector3<f64>,
        observer_helio_position: Matrix3<f64>,
    ) -> GaussObs {
        GaussObs {
            idx_obs,
            ra,
            dec,
            time: mjd_time,
            observer_helio_position,
        }
    }

    /// Generate a lazy sequence of noisy triplet realizations for orbit determination.
    ///
    /// This iterator always yields the **original triplet first**, followed by
    /// `n_realizations` additional synthetic copies where each RA/DEC coordinate
    /// is perturbed by Gaussian noise. The noise model is:
    ///
    /// * Draw `z ~ N(0, 1)` (standard normal) independently for each coordinate,
    /// * Scale by the 1-σ astrometric uncertainty times `noise_scale`,
    /// * Add to the nominal RA/DEC value.
    ///
    /// This function is designed for **Monte Carlo propagation of measurement
    /// errors** in the initial orbit determination (IOD) stage. The lazy design
    /// avoids building a large intermediate `Vec` and allows consumers to
    /// **short-circuit early** (e.g., once a sufficiently good orbit fit is found).
    ///
    /// Arguments
    /// -----------------
    /// * `errors_ra` – 1-σ RA uncertainties for the three observations (radians).
    /// * `errors_dec` – 1-σ DEC uncertainties for the three observations (radians).
    /// * `n_realizations` – Number of noisy copies to generate (excluding the original).
    /// * `noise_scale` – Scalar multiplier applied to the uncertainties
    ///   (e.g. `1.0` = nominal, `2.0` = twice the quoted error).
    /// * `rng` – Random number generator used to draw standard normals.
    ///
    /// Return
    /// ----------
    /// * An iterator yielding `GaussObs` values on the fly, starting with the
    ///   unperturbed triplet.
    ///
    /// See also
    /// ------------
    /// * [`generate_noisy_realizations`](crate::initial_orbit_determination::gauss::GaussObs::generate_noisy_realizations) – Eager version collecting all realizations into a `Vec`.
    /// * [`GaussObs::prelim_orbit`] – Consumes each realization to compute a Gauss preliminary orbit.
    /// * [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit) – High-level search loop that leverages this iterator with early pruning.
    pub fn realizations_iter<'a, R: Rng + 'a>(
        &'a self,
        errors_ra: &'a Vector3<f64>,
        errors_dec: &'a Vector3<f64>,
        n_realizations: usize,
        noise_scale: f64,
        rng: &'a mut R,
    ) -> impl Iterator<Item = GaussObs> + 'a {
        // Precompute scaled sigmas (σ × scale) to avoid redundant multiplications.
        let s_ra = *errors_ra * noise_scale;
        let s_dec = *errors_dec * noise_scale;

        // Internal counter:
        // i = 0 → original,
        // 1..=n_realizations → noisy copies,
        // > n_realizations → end of iteration.
        let mut i = 0usize;

        std::iter::from_fn(move || {
            if i == 0 {
                i += 1;
                // Yield the exact original triplet first (no noise).
                return Some(self.clone());
            }
            if i > n_realizations {
                // All requested noisy realizations have been produced.
                return None;
            }
            i += 1;

            // --- Draw Gaussian noise ---
            // For each RA/DEC coordinate: z ~ N(0, 1), then scale by σ.
            let (z_ra0, z_ra1, z_ra2): (f64, f64, f64) = (
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
            );
            let (z_dec0, z_dec1, z_dec2): (f64, f64, f64) = (
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
                rng.sample(StandardNormal),
            );

            // Apply perturbations independently to each of the 3 observations.
            let ra = Vector3::new(
                self.ra.x + z_ra0 * s_ra.x,
                self.ra.y + z_ra1 * s_ra.y,
                self.ra.z + z_ra2 * s_ra.z,
            );
            let dec = Vector3::new(
                self.dec.x + z_dec0 * s_dec.x,
                self.dec.y + z_dec1 * s_dec.y,
                self.dec.z + z_dec2 * s_dec.z,
            );

            // Return a synthetic GaussObs realization
            Some(GaussObs {
                idx_obs: self.idx_obs,
                ra,
                dec,
                time: self.time,
                observer_helio_position: self.observer_helio_position,
            })
        })
    }

    /// Generate noisy variants of this `GaussObs` triplet by injecting Gaussian noise
    /// into the right ascension (RA) and declination (DEC) coordinates, and
    /// collect all results into a vector.
    ///
    /// This function is the **eager counterpart** of [`realizations_iter`](crate::initial_orbit_determination::gauss::GaussObs::realizations_iter): it builds
    /// and returns a full `Vec<GaussObs>` in memory, containing:
    ///
    /// - The original, unperturbed triplet as the first element,
    /// - Followed by `n_realizations` noisy copies where each RA/DEC coordinate is
    ///   perturbed independently by Gaussian noise.
    ///
    /// The noise model is:
    /// - Draw `z ~ N(0, 1)` independently for each coordinate,
    /// - Multiply by the corresponding 1-σ uncertainty scaled by `noise_scale`,
    /// - Add to the nominal RA/DEC value.
    ///
    /// Arguments
    /// -----------------
    /// * `errors_ra` – 1-σ RA uncertainties for the three observations (radians).
    /// * `errors_dec` – 1-σ DEC uncertainties for the three observations (radians).
    /// * `n_realizations` – Number of noisy versions to generate (excluding the original).
    /// * `noise_scale` – Factor applied to uncertainties (e.g. `1.0` for nominal,
    ///   >1.0 to simulate larger errors).
    /// * `rng` – Random number generator used to sample standard normals.
    ///
    /// Return
    /// ----------
    /// * A vector of `GaussObs`:
    ///   - First element = original triplet,
    ///   - Remaining = noisy copies.
    ///
    /// See also
    /// ------------
    /// * [`realizations_iter`](crate::initial_orbit_determination::gauss::GaussObs::realizations_iter) – Lazy version that yields realizations on demand and
    ///   supports early-stop pruning.
    /// * [`GaussObs::prelim_orbit`] – Compute a preliminary Gauss solution from each realization.
    /// * [`estimate_best_orbit`](crate::observations::observations_ext::ObservationIOD::estimate_best_orbit) – End-to-end search that consumes realizations.
    pub fn generate_noisy_realizations(
        &self,
        errors_ra: &Vector3<f64>,
        errors_dec: &Vector3<f64>,
        n_realizations: usize,
        noise_scale: f64,
        rng: &mut impl Rng,
    ) -> Vec<GaussObs> {
        // Collect all items from the lazy iterator into a Vec.
        // Less memory-efficient than the iterator version, but convenient when
        // all realizations are needed at once.
        self.realizations_iter(errors_ra, errors_dec, n_realizations, noise_scale, rng)
            .collect()
    }

    /// Compute the direction cosine vector pointing toward the observed celestial body.
    ///
    /// This function converts right ascension (RA) and declination (DEC) coordinates
    /// into a 3D unit vector in the topocentric equatorial frame (ICRS/J2000),
    /// pointing from the observer toward the observed object.
    ///
    /// Arguments
    /// ----------
    /// * `ra` - Right ascension in radians (0 ≤ RA < 2π), measured eastward from the vernal equinox.
    /// * `dec` - Declination in radians (−π/2 ≤ DEC ≤ π/2), measured northward from the celestial equator.
    /// * `cos_dec` - Precomputed cosine of the declination (`cos(dec)`), used for performance optimization.
    ///
    /// Returns
    /// ----------
    /// * `Vector3<f64>` - The unit direction vector `[x, y, z]`, where:
    ///   - `x = cos(ra) * cos(dec)`
    ///   - `y = sin(ra) * cos(dec)`
    ///   - `z = sin(dec)`
    ///
    /// Notes
    /// ------
    /// * The resulting vector has norm ≈ 1 and lies on the celestial sphere.
    /// * This representation is suitable for computing line-of-sight directions or constructing the S matrix in Gauss's method.
    fn unit_vector(ra: f64, dec: f64, cos_dec: f64) -> Vector3<f64> {
        Vector3::new(
            f64::cos(ra) * cos_dec,
            f64::sin(ra) * cos_dec,
            f64::sin(dec),
        )
    }

    /// Construct the matrix of unit direction vectors for the observation triplet.
    ///
    /// This method returns a 3×3 matrix where each column is a unit vector pointing
    /// from the observer toward the target object at one of the three observation epochs.
    /// The unit vectors are computed from the corresponding right ascension (RA) and
    /// declination (DEC) values stored in the `GaussObs` struct.
    ///
    /// Returns
    /// -------
    /// * `Matrix3<f64>` — A 3×3 matrix where:
    ///     - Column 0 corresponds to the direction at the first observation,
    ///     - Column 1 to the second observation,
    ///     - Column 2 to the third observation.
    ///
    /// Each column vector has components:
    ///     - x = cos(RA) * cos(DEC)
    ///     - y = sin(RA) * cos(DEC)
    ///     - z = sin(DEC)
    ///
    /// Notes
    /// ------
    /// * The resulting matrix is expressed in the topocentric equatorial frame (ICRS/J2000).
    /// * This matrix is used in Gauss’s method to solve for the position of the object relative to the observer.
    /// * It is typically inverted to construct the linear system that links line-of-sight geometry to barycentric positions.
    fn unit_matrix(&self) -> Matrix3<f64> {
        Matrix3::from_columns(&[
            GaussObs::unit_vector(self.ra[0], self.dec[0], f64::cos(self.dec[0])),
            GaussObs::unit_vector(self.ra[1], self.dec[1], f64::cos(self.dec[1])),
            GaussObs::unit_vector(self.ra[2], self.dec[2], f64::cos(self.dec[2])),
        ])
    }

    /// Prepare preliminary quantities for Gauss’s method of initial orbit determination.
    ///
    /// This method computes all the geometrical and temporal quantities needed to set up
    /// the 8th-degree polynomial used in Gauss’s method. It includes the normalized
    /// time intervals (`tau1`, `tau3`), the matrix of unit direction vectors (`S`),
    /// its inverse (`S⁻¹`), and the coefficient vectors `a` and `b` used to construct
    /// the polynomial equation for the object's topocentric distance.
    ///
    /// Returns
    /// -------
    /// * `tau1`: Scaled time interval between first and second observations (g * (t1 - t2))
    /// * `tau3`: Scaled time interval between third and second observations (g * (t3 - t2))
    /// * `unit_matrix`: 3×3 matrix of unit direction vectors from observer to object (S)
    /// * `inv_unit_matrix`: Inverse of `unit_matrix`, used in solving for the object’s position
    /// * `vector_a`: Coefficients for linear combination of observer positions
    /// * `vector_b`: Coefficients for quadratic correction based on topocentric distance
    ///
    /// Errors
    /// ------
    /// Returns `OutfitError::SingularDirectionMatrix` if the matrix of direction vectors
    /// is singular (typically when the three directions are nearly coplanar or colinear).
    ///
    /// Notes
    /// ------
    /// * This method implements the initialization described in classic orbit determination literature.
    ///
    /// # See also
    /// * [`GAUSS_GRAV`] – Gaussian gravitational constant used for time normalization.
    pub fn gauss_prelim(&self) -> GaussPrelimResult {
        let tau1 = GAUSS_GRAV * (self.time[0] - self.time[1]);
        let tau3 = GAUSS_GRAV * (self.time[2] - self.time[1]);
        let tau13 = tau3 - tau1;
        let vector_a = Vector3::new(tau3 / tau13, -1.0, -(tau1 / tau13));
        let vector_b = Vector3::new(
            vector_a[0] * (tau13.powi(2) - tau3.powi(2)) / 6.0,
            0.0,
            vector_a[2] * (tau13.powi(2) - tau1.powi(2)) / 6.0,
        );

        let unit_matrix = self.unit_matrix();
        let inv_unit_matrix = unit_matrix
            .try_inverse()
            .ok_or(OutfitError::SingularDirectionMatrix)?;

        Ok((tau1, tau3, unit_matrix, inv_unit_matrix, vector_a, vector_b))
    }

    /// Compute the coefficients of the 8th-degree polynomial in Gauss’s orbit determination method.
    ///
    /// This function calculates the three non-zero coefficients of the 8th-degree scalar polynomial
    /// used to solve for the geocentric distance (r₂) of the object at the central observation time.
    /// The polynomial is of the form:
    ///
    /// ```text
    /// r2^8 + c6 * r2^6 + c3 * r2^3 + c0 = 0
    /// ```
    ///
    /// where the coefficients (c₆, c₃, c₀) encapsulate geometric constraints from the observer’s
    /// positions and directions at the three observation epochs.
    ///
    /// Arguments
    /// ---------
    /// * `unit_matrix` – 3×3 matrix of unit direction vectors (from observer to object), as returned by `unit_matrix()`.
    /// * `unit_invmatrix` – Inverse of `unit_matrix`, used to solve for line-of-sight distances.
    /// * `vector_a` – First-order combination vector from `gauss_prelim()`, for linear interpolation.
    /// * `vector_b` – Second-order combination vector from `gauss_prelim()`, for parabolic correction.
    ///
    /// Returns
    /// --------
    /// * Tuple `(c6, c3, c0)`:
    ///     - `c6`: coefficient of r₂⁶,
    ///     - `c3`: coefficient of r₂³,
    ///     - `c0`: constant term of the polynomial.
    ///
    /// Remarks
    /// --------
    /// * The resulting polynomial is derived from the scalar projection of the interpolated position vector
    ///   onto the unit direction at the central epoch, ensuring consistency between dynamic and geometric models.
    ///
    /// # See also
    /// * [`GaussObs::gauss_prelim`] – for computing `vector_a`, `vector_b`, and the unit direction matrix.
    pub fn coeff_eight_poly(
        &self,
        unit_matrix: &Matrix3<f64>,
        unit_invmatrix: &Matrix3<f64>,
        vector_a: &Vector3<f64>,
        vector_b: &Vector3<f64>,
    ) -> (f64, f64, f64) {
        let observer_position_t = self.observer_helio_position.transpose();
        let ra = self.observer_helio_position * vector_a;
        let rb = self.observer_helio_position * vector_b;

        let second_row_t = unit_invmatrix.row(1).transpose();
        let a2star = second_row_t.dot(&ra);
        let b2star = second_row_t.dot(&rb);

        let r22 = observer_position_t
            .row(1)
            .component_mul(&observer_position_t.row(1))
            .sum();
        let s2r2 = unit_matrix
            .column(1)
            .transpose()
            .dot(&observer_position_t.row(1));

        (
            -(a2star.powi(2)) - r22 - (2.0 * a2star * s2r2),
            -(2.0 * b2star * (a2star + s2r2)),
            -(b2star.powi(2)),
        )
    }

    /// Solve the 8th-degree polynomial for topocentric distance using the Aberth method.
    ///
    /// This function numerically computes the real positive roots of the scalar polynomial
    /// of degree 8 used in Gauss's initial orbit determination. It applies the Aberth method
    /// to find all complex roots and then filters them based on physical criteria.
    ///
    /// Only roots with negligible imaginary part and strictly positive real part
    /// are retained, as required by the physical model (the object must be at a
    /// positive distance from the observer).
    ///
    /// Arguments
    /// ---------
    /// * `polynom` – Array of 9 coefficients `[c₀, ..., c₈]` representing the polynomial:
    ///   c₀ + c₁·r + c₂·r² + ... + c₈·r⁸ = 0
    /// * `max_iterations` – Maximum number of iterations for the Aberth root-finding algorithm.
    /// * `aberth_epsilon` – Convergence threshold for the Aberth iterations.
    /// * `root_acceptance_epsilon` – Threshold on the imaginary part below which a root is considered real.
    ///
    /// Returns
    /// --------
    /// * `Ok(Vec<f64>)` – List of real, positive roots (the possible values of the object's geocentric distance at epoch 2).
    /// * `Err(OutfitError::PolynomialRootFindingFailed)` – If Aberth iteration fails to converge.
    ///
    /// Remarks
    /// --------
    /// * The physical validity of a root is determined by two conditions:
    ///     - its imaginary part is smaller than `root_acceptance_epsilon`,
    ///     - its real part is strictly greater than zero.
    /// * The root-finding is powered by the `aberth` crate.
    ///
    /// # See also
    /// * [`GaussObs::coeff_eight_poly`] – for computing the coefficients of the polynomial.
    pub fn solve_8poly(
        &self,
        polynom: &[f64; 9],
        max_iterations: u32,
        aberth_epsilon: f64,
        root_acceptance_epsilon: f64,
    ) -> Result<Vec<f64>, OutfitError> {
        let roots = aberth(polynom, max_iterations, aberth_epsilon);
        match roots.stop_reason {
            StopReason::Converged(_) | StopReason::MaxIteration(_) => {
                Ok(roots
                    .iter()
                    // root criteria acceptance: root must be real and real part must be positive
                    .filter(|complex| complex.re > 0. && complex.im.abs() < root_acceptance_epsilon)
                    .map(|complex| complex.re)
                    .collect::<Vec<f64>>())
            }
            StopReason::Failed(_) => Err(OutfitError::PolynomialRootFindingFailed),
        }
    }

    /// Compute the asteroid position vectors at each observation time and the corrected reference epoch.
    ///
    /// This function reconstructs the heliocentric position of the observed body at the three observation
    /// epochs using Gauss’s method. It solves for the topocentric distances `ρ₁`, `ρ₂`, `ρ₃` by inverting the
    /// geometry defined by the observer positions and line-of-sight unit vectors.
    ///
    /// Arguments
    /// ---------
    /// * `iod_params` – Parameters controlling the IOD process, including minimum acceptable distance `ρ₂`.
    /// * `unit_matrix` – 3×3 matrix of unit direction vectors from observer to object (one per epoch).
    /// * `unit_matinv` – Inverse of `unit_matrix`, used to solve for the scalar distances ρ₁, ρ₂, ρ₃.
    /// * `vector_c` – Coefficient vector computed from a candidate root of the 8th-degree polynomial;
    ///   used to combine observer positions linearly to match the line-of-sight geometry.
    ///
    /// Returns
    /// --------
    /// * `Ok((position_matrix, reference_epoch))`:
    ///     - `position_matrix`: 3×3 matrix where each column is the heliocentric Cartesian position
    ///       of the asteroid at the corresponding observation time (in AU, J2000 equatorial frame),
    ///     - `reference_epoch`: corrected central epoch (observation time minus light travel time ρ₂ / c).
    ///
    /// Errors
    /// --------
    /// Returns `OutfitError::SpuriousRootDetected` if the computed topocentric distance ρ₂ is too small
    /// (typically < 0.01 AU), indicating a non-physical or numerically unstable root.
    ///
    /// Remarks
    /// --------
    /// * The returned `reference_epoch` accounts for the light-time delay correction (aberration).
    ///
    /// # See also
    /// * [`GaussObs::accept_root`] – filters valid roots using this method.
    /// * [`VLIGHT_AU`] – speed of light in AU/day.
    pub fn position_vector_and_reference_epoch(
        &self,
        iod_params: &IODParams,
        unit_matrix: &Matrix3<f64>,
        unit_matinv: &Matrix3<f64>,
        vector_c: &Vector3<f64>,
    ) -> Result<(Matrix3<f64>, f64), OutfitError> {
        let obs_pos_t = self.observer_helio_position;
        let gcap = obs_pos_t * vector_c;
        let crhom = unit_matinv * gcap;
        let rho: Vector3<f64> = -crhom.component_div(vector_c);
        if rho[1] < iod_params.min_rho2_au {
            return Err(OutfitError::SpuriousRootDetected);
        }
        let rho_unit = Matrix3::from_columns(&[
            rho[0] * unit_matrix.column(0),
            rho[1] * unit_matrix.column(1),
            rho[2] * unit_matrix.column(2),
        ]);
        let reference_epoch = self.time[1] - rho[1] / VLIGHT_AU;

        Ok((obs_pos_t + rho_unit, reference_epoch))
    }

    /// Estimate the asteroid's velocity vector at the central observation time using Gibbs’ method.
    ///
    /// This function computes the velocity of the observed body at the second observation epoch
    /// using a refined form of Gibbs' method. It relies on the positions of the object at three
    /// distinct times and uses normalized time intervals to approximate the derivative of position.
    ///
    /// Arguments
    /// ---------
    /// * `ast_pos_vector` – 3×3 matrix of Cartesian asteroid positions (in AU), where:
    ///     - Column 0: position at first observation time (t1),
    ///     - Column 1: position at second (central) observation time (t2),
    ///     - Column 2: position at third observation time (t3).
    ///
    /// * `tau1` – Scaled time interval between the first and second observations (GAUSS_GRAV × (t1 − t2)).
    /// * `tau3` – Scaled time interval between the third and second observations (GAUSS_GRAV × (t3 − t2)).
    ///
    /// Returns
    /// --------
    /// * `Vector3<f64>` – Estimated heliocentric velocity vector of the asteroid at the time of the second observation (in AU/day).
    ///
    /// Remarks
    /// --------
    /// * The implementation follows the classical formulation used in Gauss’s method with finite-difference
    ///   approximations and gravity-based scaling.
    ///
    /// # See also
    /// * [`GaussObs::position_vector_and_reference_epoch`] – computes `ast_pos_vector` used as input here.
    /// * [`GAUSS_GRAV`] – Gaussian gravitational constant.
    pub fn gibbs_correction(
        &self,
        ast_pos_vector: &Matrix3<f64>,
        tau1: f64,
        tau3: f64,
    ) -> Vector3<f64> {
        let tau13 = tau3 - tau1;

        // Compute inverse cube of distances at each observation time
        let r1m3 = 1.0 / ast_pos_vector.column(0).norm().powi(3); // at t1
        let r2m3 = 1.0 / ast_pos_vector.column(1).norm().powi(3); // at t2
        let r3m3 = 1.0 / ast_pos_vector.column(2).norm().powi(3); // at t3

        // Compute time-dependent scalar weights for finite difference approximation
        let d1 = tau3 * (r1m3 / 12.0 - 1.0 / (tau1 * tau13));
        let d2 = (tau1 + tau3) * (r2m3 / 12.0 - 1.0 / (tau1 * tau3));
        let d3 = -tau1 * (r3m3 / 12.0 + 1.0 / (tau3 * tau13));

        // Combine the three positions using the scalar weights
        let d_vect = Vector3::new(-d1, d2, d3);

        // Compute the velocity vector as a weighted sum of position components
        Vector3::new(
            GAUSS_GRAV * ast_pos_vector.row(0).dot(&d_vect.transpose()), // vx
            GAUSS_GRAV * ast_pos_vector.row(1).dot(&d_vect.transpose()), // vy
            GAUSS_GRAV * ast_pos_vector.row(2).dot(&d_vect.transpose()), // vz
        )
    }

    /// Evaluate whether a root of the 8th-degree polynomial leads to a valid preliminary orbit.
    ///
    /// This method checks if a given root (candidate topocentric distance `r2`) leads to a
    /// physically meaningful orbit according to the Gauss method. It reconstructs the object's
    /// position and velocity vectors, applies the eccentricity and perihelion filters,
    /// and returns the full geometric state if accepted.
    ///
    /// Arguments
    /// ---------
    /// * `iod_params` – Parameters controlling the IOD process, including eccentricity and perihelion limits.
    /// * `root` – A real, positive root of the polynomial equation for r₂ (distance at central epoch).
    /// * `unit_matrix` – 3×3 matrix of unit direction vectors (from observer to object).
    /// * `inv_unit_matrix` – Inverse of the direction matrix.
    /// * `vect_a` – First-order combination vector (from `gauss_prelim()`).
    /// * `vect_b` – Second-order combination vector (from `gauss_prelim()`).
    /// * `tau1` – Scaled time between first and second observation (GAUSS_GRAV × (t1 − t2)).
    /// * `tau3` – Scaled time between third and second observation (GAUSS_GRAV × (t3 − t2)).
    ///
    /// Returns
    /// --------
    /// * `Some((position_matrix, velocity_vector, reference_epoch))` if the root is valid,
    /// * `None` if the resulting solution is not physically acceptable (e.g. eccentricity too high).
    ///
    /// Remarks
    /// --------
    /// * The eccentricity and perihelion distance are checked against internal thresholds
    ///   using the [`eccentricity_control`] function.
    ///
    /// # See also
    /// * [`eccentricity_control`] – checks if an orbit meets physical bounds.
    /// * [`GaussObs::position_vector_and_reference_epoch`] – computes the full 3×3 object position matrix.
    /// * [`GaussObs::gibbs_correction`] – estimates the velocity at central epoch from positions.
    #[allow(clippy::too_many_arguments)]
    pub fn accept_root(
        &self,
        iod_params: &IODParams,
        root: f64,
        unit_matrix: &Matrix3<f64>,
        inv_unit_matrix: &Matrix3<f64>,
        vect_a: &Vector3<f64>,
        vect_b: &Vector3<f64>,
        tau1: f64,
        tau3: f64,
    ) -> Option<(Matrix3<f64>, Vector3<f64>, f64)> {
        // Compute r₂⁻³ once (used in vector c)
        let r2m3 = 1.0 / root.powi(3);

        // Construct the vector c used to linearly combine the observer positions
        let vector_c = Vector3::new(
            vect_a[0] + vect_b[0] * r2m3,
            -1.0,
            vect_a[2] + vect_b[2] * r2m3,
        );

        // Compute asteroid position vectors at each observation time, and the reference epoch
        let Some((ast_pos_all_time, reference_epoch)) = self
            .position_vector_and_reference_epoch(
                iod_params,
                unit_matrix,
                inv_unit_matrix,
                &vector_c,
            )
            .ok()
        else {
            return None; // root leads to invalid geometry (e.g., ρ₂ < 0.01)
        };

        // Extract the position vector at the central epoch (t2)
        let ast_pos_second_time = ast_pos_all_time.column(1).into();

        // Estimate the velocity vector at the second observation via Gibbs correction
        let asteroid_vel = self.gibbs_correction(&ast_pos_all_time, tau1, tau3);

        // Apply eccentricity and perihelion distance control
        let (is_accepted, _, _, _) = eccentricity_control(
            &ast_pos_second_time,
            &asteroid_vel,
            iod_params.max_perihelion_au,
            iod_params.max_ecc,
        )?;

        // Accept only orbits within eccentricity and perihelion limits
        if is_accepted {
            Some((ast_pos_all_time, asteroid_vel, reference_epoch))
        } else {
            None
        }
    }

    /// Convert Cartesian position and velocity vectors to Keplerian orbital elements.
    ///
    /// This function transforms a state vector (position and velocity) from the equatorial
    /// J2000 frame into the mean ecliptic J2000 frame, then computes the corresponding
    /// Keplerian orbital elements using classical two-body dynamics.
    ///
    /// Arguments
    /// ---------
    /// * `state` – The outfit global state containing the rotation matrix
    /// * `asteroid_position` – Cartesian heliocentric position vector of the object (in AU), in equatorial J2000 frame.
    /// * `asteroid_velocity` – Cartesian heliocentric velocity vector of the object (in AU/day), in equatorial J2000 frame.
    /// * `reference_epoch` – Epoch (in MJD TT) corresponding to the state vector, used as the reference time for the elements.
    ///
    /// Returns
    /// --------
    /// * `KeplerianElements` – Struct containing the classical orbital elements of the body, with the following fields:
    ///     - `reference_epoch`: epoch of the elements in Modified Julian Date (MJD TT),
    ///     - `semi_major_axis`: semi-major axis in astronomical units (AU),
    ///     - `eccentricity`: orbital eccentricity (unitless),
    ///     - `inclination`: inclination in radians,
    ///     - `ascending_node_longitude`: longitude of ascending node in radians,
    ///     - `periapsis_argument`: argument of periapsis in radians,
    ///     - `mean_anomaly`: mean anomaly at epoch in radians.
    ///
    /// Remarks
    /// --------
    /// * The transformation matrix from equatorial to ecliptic frame is computed using [`rotpn`].
    /// * Orbital elements are computed using the function [`ccek1`] from the `orb_elem` module.
    ///
    /// See also
    /// --------
    /// * [`rotpn`] – computes the rotation matrix between celestial reference frames.
    /// * [`OrbitalElements::from_orbital_state`] – Computes classical orbital elements from a Cartesian state vector.
    /// * [`KeplerianElements`] – definition of the orbital elements struct.
    fn compute_orbit_from_state(
        &self,
        state: &Outfit,
        &asteroid_position: &Vector3<f64>,
        &asteroid_velocity: &Vector3<f64>,
        reference_epoch: f64,
    ) -> Result<OrbitalElements, OutfitError> {
        // get the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let roteqec = state.get_rot_equmj2000_to_eclmj2000();

        // Apply the transformation to position and velocity vectors
        let matrix_elc_transform = roteqec.transpose();
        let ecl_pos = matrix_elc_transform * asteroid_position;
        let ecl_vel = matrix_elc_transform * asteroid_velocity;

        // Compute the classical orbital elements from the state vector and return orbital elements in a structured form
        Ok(OrbitalElements::from_orbital_state(
            &ecl_pos,
            &ecl_vel,
            reference_epoch,
        ))
    }

    /// Visit all real-positive roots of the fixed 8th-degree Gauss polynomial and call a visitor for each candidate.
    ///
    /// This helper runs the Aberth–Ehrlich solver on a degree-8 polynomial (coefficients `[a, b, ..., x^8]`),
    /// then iterates **once** over the computed complex roots without allocating. For each root whose
    /// imaginary part is below `root_imag_eps` in magnitude and whose real part is strictly positive,
    /// the visitor `on_root` is invoked with the real value `r`.
    ///
    /// The visitor controls iteration:
    /// return `ControlFlow::Continue(())` to keep scanning, or `ControlFlow::Break(())` to short-circuit
    /// as soon as a satisfactory root has been processed.
    ///
    /// Arguments
    /// -----------------
    /// * `poly`: The 9 coefficients of the 8th-degree polynomial `[c0, c1, …, c8]`.
    /// * `max_iterations`: Upper bound on Aberth iterations.
    /// * `aberth_epsilon`: Convergence threshold for Aberth updates (delta on roots).
    /// * `root_imag_eps`: Tolerance on `|Im(z)|` used to deem a complex root “real”.
    /// * `on_root`: Visitor called for each real-positive root (in solver order).  
    ///   Return `ControlFlow::Break(())` to stop early, or `ControlFlow::Continue(())` to keep iterating.
    ///
    /// Return
    /// ----------
    /// * `Ok(())` if the solver converged or reached `max_iterations` (all admissible roots were visited).
    /// * `Err(OutfitError::PolynomialRootFindingFailed)` if Aberth reported failure (NaN/Inf, etc.).
    ///
    /// Notes
    /// ----------
    /// * Only **real** (within `root_imag_eps`) and **positive** roots are passed to the visitor.
    /// * The order of visited roots is the solver’s order and is not guaranteed to be sorted by magnitude.
    /// * This function performs no heap allocation and visits each root at most once.
    ///
    /// See also
    /// ------------
    /// * [`aberth::aberth`] – Const-generic Aberth–Ehrlich root finder.
    /// * [`GaussObs::coeff_eight_poly`] – Builds the (sparse) coefficients of the degree-8 polynomial.
    /// * [`GaussObs::prelim_orbit_all`] – Scans all roots and collects up to three Gauss solutions.
    /// * [`GaussObs::prelim_orbit`] – Picks a single “best” solution (corrected preferred).
    /// * [`GaussObs::accept_root`] – Physical admissibility test and preliminary state construction.
    #[inline]
    fn visit_real_positive_roots(
        poly: &[f64; 9],
        max_iterations: u32,
        aberth_epsilon: f64,
        root_imag_eps: f64,
        mut on_root: impl FnMut(f64) -> ControlFlow<(), ()>,
    ) -> Result<(), OutfitError> {
        let roots = aberth(poly, max_iterations, aberth_epsilon);

        match roots.stop_reason {
            aberth::StopReason::Converged(_) | aberth::StopReason::MaxIteration(_) => {
                for z in roots.iter() {
                    if z.re > 0.0 && z.im.abs() < root_imag_eps {
                        if let ControlFlow::Break(()) = on_root(z.re) {
                            break;
                        }
                    }
                }
                Ok(())
            }
            aberth::StopReason::Failed(_) => Err(OutfitError::PolynomialRootFindingFailed),
        }
    }

    /// Try to accept a single Gauss root and build a preliminary dynamical state `(r(t_i), v(t2), epoch)`.
    ///
    /// Given a candidate heliocentric distance at the central epoch `t2` (`r2`), this helper applies the
    /// same physical/admissibility checks as [`accept_root`] and, if successful, reconstructs:
    /// - the heliocentric positions of the body at the three observation epochs (columns correspond to `t1,t2,t3`),
    /// - the heliocentric velocity at `t2`,
    /// - and the reference epoch used to interpret the state.
    ///
    /// No heap allocation is performed; inputs are passed by reference and the function is `#[inline]`.
    ///
    /// Arguments
    /// ---------
    /// * `iod_params` – Parameters controlling the IOD process, including minimum acceptable distance `ρ₂`.
    /// * `r2`: Heliocentric distance at the second epoch `t2` (AU).
    /// * `unit_m`: Direction (line-of-sight) matrix `S`, columns are the unit vectors for each observation.
    /// * `inv_unit_m`: Precomputed inverse of `S`.  
    ///   _Note_: consider passing a linear-system solver (e.g. LU/QR) instead of an explicit inverse for better numerical stability.
    /// * `a`, `b`: Gauss interpolation vectors (depend on `tau1`, `tau3`).
    /// * `tau1`, `tau3`: Time offsets (`gk * (t1 - t2)` and `gk * (t3 - t2)`) in Gauss units.
    ///
    /// Return
    /// ----------
    /// * `Some((positions_all_times, velocity_t2, epoch_ref))` if the root passes physical checks and reconstruction succeeds.  
    ///   - `positions_all_times`: `3×3` matrix of heliocentric positions (AU), columns = `[t1 | t2 | t3]`.  
    ///   - `velocity_t2`: heliocentric velocity at `t2` (AU/day).  
    ///   - `epoch_ref`: reference epoch used for the state (usually MJD(TT)).
    /// * `None` if the root is rejected (e.g. spurious geometry, `rho2` too small, eccentricity/perihelion limits).
    ///
    /// Notes
    /// ----------
    /// * The acceptance criteria are delegated to [`accept_root`] (e.g. positivity of ranges, minimum `rho2`,
    ///   eccentricity and perihelion bounds, and geometric consistency).
    /// * The middle column `positions_all_times.column(1)` corresponds to the state at `t2`.
    ///
    /// See also
    /// ------------
    /// * [`accept_root`] – Physical admissibility check and preliminary state construction.
    /// * [`coeff_eight_poly`] – Builds the sparse coefficients of the 8th-degree Gauss polynomial.
    /// * [`prelim_orbit_all`] – Scans all admissible roots and collects up to three solutions.
    /// * [`velocity_correction`] – Lagrange-based velocity update from the three-epoch geometry.
    #[inline]
    #[allow(clippy::too_many_arguments)]
    fn try_accept_root(
        &self,
        iod_params: &IODParams,
        r2: f64,
        unit_m: &Matrix3<f64>,
        inv_unit_m: &Matrix3<f64>,
        a: &Vector3<f64>,
        b: &Vector3<f64>,
        tau1: f64,
        tau3: f64,
    ) -> Option<(nalgebra::Matrix3<f64>, Vector3<f64>, f64)> {
        self.accept_root(iod_params, r2, unit_m, inv_unit_m, a, b, tau1, tau3)
    }

    /// Build a `GaussResult` from a position/velocity state at `t2`, applying aberration timing in
    /// [`compute_orbit_from_state`]. Returns `None` if orbit construction fails.
    ///
    /// Arguments
    /// -----------------
    /// * `pos_all_time`: `3×3` heliocentric positions at `t1|t2|t3` (AU).
    /// * `vel_t2`: heliocentric velocity at `t2` (AU/day).
    /// * `epoch`: reference epoch for the state (MJD TT).
    /// * `corrected`: set `true` if the state is post-correction, `false` if preliminary.
    ///
    /// Return
    /// ----------
    /// * `Some(GaussResult::CorrectedOrbit|PrelimOrbit)` on success, `None` otherwise.
    ///
    /// See also
    /// ------------
    /// * [`compute_orbit_from_state`] – Converts `(r(t2), v(t2), epoch)` to orbital elements with aberration.
    #[inline]
    fn build_result(
        &self,
        state: &Outfit,
        pos_all_time: &Matrix3<f64>,
        vel_t2: &Vector3<f64>,
        epoch: f64,
        corrected: bool,
    ) -> Option<GaussResult> {
        let r_t2: Vector3<f64> = pos_all_time.column(1).into();
        match self.compute_orbit_from_state(state, &r_t2, vel_t2, epoch) {
            Ok(orbit) if corrected => Some(GaussResult::CorrectedOrbit(orbit)),
            Ok(orbit) => Some(GaussResult::PrelimOrbit(orbit)),
            Err(_) => None,
        }
    }

    /// Scan all real-positive roots of the Gauss polynomial and collect up to three acceptable orbits.
    ///
    /// This routine enumerates admissible roots of the 8th-degree distance polynomial, reconstructs a
    /// preliminary state for each candidate root, attempts a velocity correction, then records either the
    /// corrected orbit (on success) or the preliminary one (if correction fails). The collection is capped
    /// at three solutions to keep the output concise and representative.
    ///
    /// Algorithm
    /// -----------------
    /// 1. Compute core quantities for Gauss’ method: direction matrix `S` (unit line-of-sight vectors),
    ///    its inverse, time offsets `tau₁, tau₃`, and the interpolation vectors `A/B`.
    /// 2. Build the sparse 8th-degree polynomial coefficients `(c₆, c₃, c₀)`.
    /// 3. Visit real-positive roots using the Aberth–Ehrlich solver; optionally reject via a cheap `r₂`
    ///    plausibility window before expensive checks.
    /// 4. For each root:
    ///    - Apply geometric/physical tests and reconstruct positions `(t₁ | t₂ | t₃)`, `v(t₂)`, and the epoch.
    ///    - Attempt an iterative velocity correction; on success, store a **CorrectedOrbit**, otherwise a **PrelimOrbit**.
    /// 5. Stop when three solutions have been collected or when all roots are exhausted.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: Global context (ephemerides, constants, settings).
    /// * `iod_params`: Parameters controlling the IOD process, including root filtering and correction settings.
    ///
    /// Return
    /// ----------
    /// * `Ok(Vec<GaussResult>)` containing **1..=3** solutions in discovery order.
    /// * `Err(OutfitError::GaussNoRootsFound)` if no admissible root yields a solution.
    ///
    /// Notes
    /// ----------
    /// * Only roots with `Im(z)` below `ROOT_IMAG_EPS` and `Re(z) > 0` are considered; a light plausibility filter
    ///   `[R2_MIN, R2_MAX]` may further prune candidates.
    /// * Numerical stability downstream can improve by replacing `inv(S)` with a linear solve (LU/QR).
    /// * Root scanning does not allocate; solution accumulation uses a `SmallVec<[..; 3]>`.
    ///
    /// See also
    /// ------------
    /// * [`prelim_orbit`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit) – Picks a single “best” solution (corrected preferred).
    /// * [`pos_and_vel_correction`](crate::initial_orbit_determination::gauss::GaussObs::pos_and_vel_correction) – Iterative velocity update (Lagrange f/g).
    /// * [`IODParams`] – Configuration parameters for the IOD process.
    pub fn prelim_orbit_all(
        &self,
        state: &Outfit,
        iod_params: &IODParams,
    ) -> Result<Vec<GaussResult>, OutfitError> {
        // 1) Core Gauss quantities
        let (tau1, tau3, unit_m, inv_unit_m, vec_a, vec_b) = self.gauss_prelim()?;

        // 2) Degree-8 polynomial (only 3 non-zero coefficients)
        let (c6, c3, c0) = self.coeff_eight_poly(&unit_m, &inv_unit_m, &vec_a, &vec_b);
        let poly = [c0, 0.0, 0.0, c3, 0.0, 0.0, c6, 0.0, 1.0];

        // --- Cheap prefilter: Descartes' rule of signs on (0, +∞).
        // If zero, we know there is NO positive real root; skip Aberth entirely.
        let ub = descartes_upper_bound_deg8_sparse(c0, c3, c6, 0.0);
        if ub == 0 {
            return Err(OutfitError::GaussNoRootsFound);
        }

        // Accumulate up to three solutions (stack-first, then into_vec())
        let mut solutions: SmallVec<[GaussResult; 4]> = SmallVec::new();

        // 3–5) Enumerate admissible roots → accept → correct (if possible) → push
        Self::visit_real_positive_roots(
            &poly,
            iod_params.aberth_max_iter,
            iod_params.aberth_eps,
            iod_params.root_imag_eps,
            |r2| {
                // Cheap plausibility filter on r2 (optional)
                if !(iod_params.r2_min_au..=iod_params.r2_max_au).contains(&r2) {
                    return ControlFlow::Continue(());
                }

                match self.try_accept_root(
                    iod_params,
                    r2,
                    &unit_m,
                    &inv_unit_m,
                    &vec_a,
                    &vec_b,
                    tau1,
                    tau3,
                ) {
                    None => ControlFlow::Continue(()),

                    Some((pos_all, v_pre, epoch_ref)) => {
                        // Attempt velocity correction first; prefer corrected orbits
                        if let Some((pos_cor, v_cor, epoch_cor)) = self.pos_and_vel_correction(
                            iod_params,
                            &pos_all,
                            &v_pre,
                            &unit_m,
                            &inv_unit_m,
                            iod_params.max_perihelion_au,
                            iod_params.max_ecc,
                            iod_params.newton_eps,
                            iod_params.newton_max_it,
                        ) {
                            if let Some(res) =
                                self.build_result(state, &pos_cor, &v_cor, epoch_cor, true)
                            {
                                if solutions.len() < iod_params.max_tested_solutions {
                                    solutions.push(res);
                                }
                            }
                        } else if let Some(res) =
                            self.build_result(state, &pos_all, &v_pre, epoch_ref, false)
                        {
                            if solutions.len() < iod_params.max_tested_solutions {
                                solutions.push(res);
                            }
                        }

                        if solutions.len() >= iod_params.max_tested_solutions {
                            ControlFlow::Break(())
                        } else {
                            ControlFlow::Continue(())
                        }
                    }
                }
            },
        )?;

        if solutions.is_empty() {
            Err(OutfitError::GaussNoRootsFound)
        } else {
            Ok(solutions.into_vec())
        }
    }

    /// Estimate a single initial orbit from three astrometric observations using Gauss’ method.
    ///
    /// This wrapper delegates to [`prelim_orbit_all`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit_all) to enumerate and evaluate all admissible
    /// roots of the 8th-degree Gauss polynomial, then applies a simple selection policy to
    /// return **one** orbit:
    /// 1) prefer any **corrected** solution (velocity correction converged),
    /// 2) otherwise return the first **preliminary** solution.
    ///
    /// Arguments
    /// -----------------
    /// * `state`: The global context (ephemerides, constants, settings).
    /// * `iod_params`: Parameters controlling the IOD process, including root filtering and correction settings.
    ///
    /// Return
    /// ----------
    /// * `Ok(GaussResult::CorrectedOrbit)` if at least one corrected solution exists.
    /// * `Ok(GaussResult::PrelimOrbit)` if no corrected solution exists but a preliminary one does.
    /// * `Err(OutfitError::GaussNoRootsFound)` if no admissible root yields a solution.
    ///
    /// Notes
    /// ----------
    /// * The ordering of candidates is the discovery order used in [`prelim_orbit_all`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit_all).
    /// * If you need access to **all** acceptable solutions (for scoring, diagnostics, or tie-breaking),
    ///   call [`prelim_orbit_all`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit_all) directly and implement your own selection.
    ///
    /// See also
    /// ------------
    /// * [`prelim_orbit_all`](crate::initial_orbit_determination::gauss::GaussObs::prelim_orbit_all) – Enumerates and returns up to three acceptable solutions.
    /// * [`pos_and_vel_correction`](crate::initial_orbit_determination::gauss::GaussObs::pos_and_vel_correction) – Iterative velocity update (Lagrange f/g).
    /// * [`IODParams`] – Configuration parameters for the IOD process.
    pub fn prelim_orbit(
        &self,
        state: &Outfit,
        iod_params: &IODParams,
    ) -> Result<GaussResult, OutfitError> {
        let all = self.prelim_orbit_all(state, iod_params)?;
        if let Some(best_corr) = all
            .iter()
            .find(|s| matches!(s, GaussResult::CorrectedOrbit(_)))
        {
            return Ok(best_corr.clone());
        }
        all.into_iter().next().ok_or(OutfitError::GaussNoRootsFound)
    }

    /// Iteratively refine the asteroid's velocity and positions using two-sided Lagrange f–g updates.
    ///
    /// This routine performs Gauss’s correction loop centered on the middle observation (index 1).
    /// At each iteration it:
    /// 1) Updates the **central velocity** by applying two f–g velocity corrections built from
    ///    (`r1`,`r2`, `dt_01`) and (`r3`,`r2`, `dt_21`), optionally warm-started on the universal anomaly `χ`.
    /// 2) Recomputes the **Gauss C-vector** and derives the **three positions** and the **reference epoch**
    ///    via [`GaussObs::position_vector_and_reference_epoch`].
    /// 3) Applies a **dynamic acceptability** test on the central state and checks **convergence** on the
    ///    relative Frobenius norm of the position update.
    ///
    /// Arguments
    /// ---------
    /// * `iod_params` – Parameters controlling the IOD process, including minimum acceptable distance `ρ₂`.
    /// * `asteroid_position` – Initial 3×3 positions at epochs (t₁, t₂, t₃), columns are the position vectors (AU).
    /// * `asteroid_velocity` – Initial velocity at the central epoch t₂ (AU/day).
    /// * `unit_matrix` – 3×3 matrix of line-of-sight unit vectors (observer → object).
    /// * `inv_unit_matrix` – Inverse of `unit_matrix`.
    /// * `peri_max` – Maximum allowed perihelion distance (AU) used by dynamic acceptability.
    /// * `ecc_max` – Maximum allowed eccentricity used by dynamic acceptability.
    /// * `err_max` – Relative convergence threshold on positions (Frobenius norm).
    /// * `itmax` – Maximum number of correction iterations.
    ///
    /// Return
    /// ------
    /// * `Some((corrected_positions, corrected_velocity, corrected_epoch))` on success.
    /// * `None` if the correction diverges (ill-conditioned geometry, `g` issues, etc.)
    ///   or if the orbit is dynamically rejected.
    ///
    /// See also
    /// --------
    /// * [`velocity_correction_with_guess`] – Lagrange-based velocity update with χ warm-start.
    /// * [`eccentricity_control`] – Filters solutions based on dynamic acceptability.
    /// * [`GaussObs::position_vector_and_reference_epoch`] – Recomputes positions at each epoch.
    #[allow(clippy::too_many_arguments)]
    pub fn pos_and_vel_correction(
        &self,
        iod_params: &IODParams,
        asteroid_position: &Matrix3<f64>,
        asteroid_velocity: &Vector3<f64>,
        unit_matrix: &Matrix3<f64>,
        inv_unit_matrix: &Matrix3<f64>,
        peri_max: f64,
        ecc_max: f64,
        err_max: f64,
        itmax: usize,
    ) -> Option<(Matrix3<f64>, Vector3<f64>, f64)> {
        // --- Loop state (positions, velocity, reference epoch)
        let mut pos = *asteroid_position;
        let mut vel = *asteroid_velocity;
        let mut epoch = 0.0_f64;

        // Optional warm-starts for the universal anomaly χ on both sides of the middle epoch.
        // If velocity_correction_with_guess returns χ, we cache it here to speed up the next iteration.
        let mut chi_guess_01: Option<f64> = None; // for dt_01 = t1 - t2
        let mut chi_guess_21: Option<f64> = None; // for dt_21 = t3 - t2

        // Time offsets relative to the central observation (2 ≡ index 1).
        let dt_01 = self.time[0] - self.time[1];
        let dt_21 = self.time[2] - self.time[1];

        // Early exit if one of the time gaps is (near) zero → would drive g → 0 in f–g update.
        if dt_01.abs() <= f64::EPSILON || dt_21.abs() <= f64::EPSILON {
            return None;
        }

        // Iterate until convergence on the relative Frobenius norm or until the iteration budget is exhausted.
        for _iter in 0..itmax {
            // Snapshot the three position columns (stack-owned vectors).
            let r1 = pos.column(0).into_owned();
            let r2 = pos.column(1).into_owned();
            let r3 = pos.column(2).into_owned();

            // --- Two velocity corrections around the central position r2.

            // If velocity_correction_with_guess returns χ, we keep it and warm-start the next iteration.
            let res_left = velocity_correction_with_guess(
                &r1,
                &r2,
                &vel,
                dt_01,
                peri_max,
                ecc_max,
                chi_guess_01,
                iod_params.kepler_eps,
            );
            let res_right = velocity_correction_with_guess(
                &r3,
                &r2,
                &vel,
                dt_21,
                peri_max,
                ecc_max,
                chi_guess_21,
                iod_params.kepler_eps,
            );

            // Both sides must succeed; otherwise, try another outer iteration (or eventually give up).
            let ((v1, f1, g1, chi1), (v2, f2, g2, chi2)) = match (res_left, res_right) {
                (Ok(a), Ok(b)) => (a, b),
                _ => continue, // e.g., ill-conditioned local geometry this iteration
            };

            // Cache χ for warm-starting the universal Kepler solver on the next pass.
            chi_guess_01 = Some(chi1);
            chi_guess_21 = Some(chi2);

            // Defensive check: g should already be validated inside velocity_correction_with_guess.
            if !g1.is_finite() || !g2.is_finite() {
                continue;
            }

            // Average the two velocity estimates; this is robust against mild asymmetries.
            // Optional alternative (uncomment to enable): |g|-weighted average instead of plain mean.
            // let w1 = g1.abs();
            // let w2 = g2.abs();
            // let new_vel = (v1 * w1 + v2 * w2) / (w1 + w2);
            let new_vel = (v1 + v2) * 0.5;

            // --- Build the Gauss C-vector from (f, g).
            // One reciprocal to avoid two divisions.
            let fl = f1 * g2 - f2 * g1;
            if !fl.is_finite() || fl.abs() < f64::EPSILON {
                // Degenerate combination of f–g coefficients; skip this iteration.
                continue;
            }
            let inv_f = 1.0 / fl;
            let c_vec = Vector3::new(g2 * inv_f, -1.0, -g1 * inv_f);

            // --- Recompute asteroid positions and light-time-corrected epoch.
            let (new_pos, new_epoch) = match self.position_vector_and_reference_epoch(
                iod_params,
                unit_matrix,
                inv_unit_matrix,
                &c_vec,
            ) {
                Ok(v) => v,
                _ => continue, // e.g., degenerate ρ₂ estimate
            };

            // --- Dynamic acceptability at the middle epoch.
            let r_mid = new_pos.column(1).into_owned();
            let accepted = match eccentricity_control(&r_mid, &new_vel, peri_max, ecc_max) {
                Some((ok, _, _, _)) => ok,
                None => return None, // pathological case (e.g., zero angular momentum)
            };
            if !accepted {
                return None; // reject dynamically unacceptable solution
            }

            // --- Convergence check on the relative Frobenius norm of the position update.
            let denom = new_pos.norm();
            if !denom.is_finite() || denom <= f64::EPSILON {
                // Degenerate position matrix; skip this iteration and try again.
                continue;
            }
            let rel_err = (new_pos - pos).norm() / denom;

            // Commit the iteration state.
            pos = new_pos;
            vel = new_vel;
            epoch = new_epoch;

            if rel_err <= err_max {
                break;
            }
        }

        Some((pos, vel, epoch))
    }
}

#[cfg(test)]
#[cfg(feature = "jpl-download")]
pub(crate) mod gauss_test {

    use super::*;
    use crate::{
        orbit_type::{keplerian_element::KeplerianElements, orbit_type_test::approx_equal},
        unit_test_global::OUTFIT_HORIZON_TEST,
    };

    #[test]
    fn test_gauss_prelim() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(
                1.082_468_037_385_525,
                0.943_580_504_794_621_6,
                0.827_376_240_789_998_6,
            ),
            time: Vector3::new(
                57028.479297592596,
                57_049.245_147_592_59,
                57_063.977_117_592_59,
            ),
            observer_helio_position: Matrix3::zeros(),
        };

        let (tau1, tau3, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        assert_eq!(tau1, -0.35721620648079105);
        assert_eq!(tau3, 0.25342080566844405);

        let unit_mat_array: [f64; 9] = unit_matrix
            .as_slice()
            .try_into()
            .expect("Conversion failed");
        assert_eq!(
            unit_mat_array,
            [
                -0.05549934652247514,
                0.46585594034226024,
                0.8831183756345503,
                -0.06972979004485365,
                0.5827357012279389,
                0.8096646582966821,
                -0.12245931009139571,
                0.6656387438390606,
                0.7361581216507068
            ]
        );

        let inv_unit_mat_array: [f64; 9] = inv_unit_matrix
            .as_slice()
            .try_into()
            .expect("Conversion failed");
        assert_eq!(
            inv_unit_mat_array,
            [
                -18.774792915974594,
                41.814279122702025,
                -23.466669573973437,
                -8.16479071034311,
                11.489343729350427,
                -2.8418335594428186,
                4.259482782736117,
                -3.432964304649723,
                0.024345794753282718
            ]
        );

        let vect_a_array: [f64; 3] = vector_a.as_slice().try_into().expect("Conversion failed");
        assert_eq!(
            vect_a_array,
            [0.41501055557783634, -1.0, 0.5849894444221637]
        );

        let vect_b_array: [f64; 3] = vector_b.as_slice().try_into().expect("Conversion failed");
        assert_eq!(
            vect_b_array,
            [0.021349212036493866, 0.0, 0.023913797385599792]
        );
    }

    #[test]
    fn test_coeff_8poly() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(
                1.082_468_037_385_525,
                0.943_580_504_794_621_6,
                0.827_376_240_789_998_6,
            ),
            time: Vector3::new(
                57028.479297592596,
                57_049.245_147_592_59,
                57_063.977_117_592_59,
            ),
            observer_helio_position: Matrix3::new(
                -0.26456661713915464,
                0.868_935_164_369_495,
                0.376_699_621_109_192_2,
                -0.589_163_185_217_412_7,
                0.723_887_251_679_477_7,
                0.313_818_651_652_458_5,
                -0.774_387_443_796_959_6,
                0.561_288_470_926_116_4,
                0.24334971075289916,
            )
            .transpose(),
        };

        let (_, _, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        let (coeff_6, coeff_3, coeff_0) =
            gauss.coeff_eight_poly(&unit_matrix, &inv_unit_matrix, &vector_a, &vector_b);

        assert_eq!(coeff_6, -2.615803718759013);
        assert_eq!(coeff_3, 2.0305173353541064);
        assert_eq!(coeff_0, -0.4771346939201045);
    }

    #[test]
    fn test_solving_polynom() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::zeros(),
            dec: Vector3::zeros(),
            time: Vector3::zeros(),
            observer_helio_position: Matrix3::zeros(),
        };

        let polynomial = [
            -0.477_134_693_920_104_8,
            0.,
            0.,
            2.0305173353541064,
            0.,
            0.,
            -2.615_803_718_759_011,
            0.,
            1.,
        ];

        let roots = gauss.solve_8poly(&polynomial, 50, 1e-6, 1e-6).unwrap();

        assert_eq!(
            roots,
            vec![1.3856312487504954, 0.7328107254669438, 0.9540135094917113]
        );
    }

    #[test]
    fn test_asteroid_position() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(
                1.082_468_037_385_525,
                0.943_580_504_794_621_6,
                0.827_376_240_789_998_6,
            ),
            time: Vector3::new(
                57028.479297592596,
                57_049.245_147_592_59,
                57_063.977_117_592_59,
            ),
            observer_helio_position: Matrix3::new(
                -0.26456661713915464,
                0.868_935_164_369_495,
                0.376_699_621_109_192_2,
                -0.589_163_185_217_412_7,
                0.723_887_251_679_477_7,
                0.313_818_651_652_458_5,
                -0.774_387_443_796_959_6,
                0.561_288_470_926_116_4,
                0.24334971075289916,
            )
            .transpose(),
        };

        let (_, _, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        let first_root: f64 = 0.732_810_725_466_943_7;
        let r2m3 = 1. / first_root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vector_a[0] + vector_b[0] * r2m3,
            -1.,
            vector_a[2] + vector_b[2] * r2m3,
        );
        let ast_pos = gauss.position_vector_and_reference_epoch(
            &IODParams::default(),
            &unit_matrix,
            &inv_unit_matrix,
            &vector_c,
        );
        assert!(ast_pos.is_err());

        let second_root: f64 = 1.3856312487504951;
        let r2m3 = 1. / second_root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vector_a[0] + vector_b[0] * r2m3,
            -1.,
            vector_a[2] + vector_b[2] * r2m3,
        );
        let (ast_pos, ref_epoch) = gauss
            .position_vector_and_reference_epoch(
                &IODParams::default(),
                &unit_matrix,
                &inv_unit_matrix,
                &vector_c,
            )
            .unwrap();

        let ast_pos_slice: [f64; 9] = ast_pos
            .as_slice()
            .try_into()
            .expect("test_asteroid_position result matrix into slice failed");

        assert_eq!(
            ast_pos_slice,
            [
                -0.28811969067349597,
                1.06663729794052,
                0.7514815481797275,
                -0.6235500510031637,
                1.0112601855976917,
                0.713100363506241,
                -0.8445850475187664,
                0.9428539454255418,
                0.6653391541170498
            ]
        );

        assert_eq!(ref_epoch, 57049.24229942721);
    }

    #[test]
    fn test_gibbs_correction() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(
                1.082_468_037_385_525,
                0.943_580_504_794_621_6,
                0.827_376_240_789_998_6,
            ),
            time: Vector3::new(
                57028.479297592596,
                57_049.245_147_592_59,
                57_063.977_117_592_59,
            ),
            observer_helio_position: Matrix3::zeros(),
        };
        let (tau1, tau3, _, _, _, _) = gauss.gauss_prelim().unwrap();

        let ast_pos = Matrix3::new(
            -0.28811969067349597,
            1.06663729794052,
            0.7514815481797275,
            -0.6235500510031637,
            1.0112601855976917,
            0.713100363506241,
            -0.8445850475187664,
            0.9428539454255418,
            0.6653391541170498,
        )
        .transpose();
        let asteroid_velocity = gauss.gibbs_correction(&ast_pos, tau1, tau3);

        let ast_vel_slice = asteroid_velocity.as_slice();
        assert_eq!(
            ast_vel_slice,
            [
                -0.015549845137774663,
                -0.003876936109837664,
                -0.0027014074002979886
            ]
        )
    }

    #[test]
    fn test_solve_orbit() {
        let env = &OUTFIT_HORIZON_TEST.0;
        let tol = 1e-13;

        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6894680985108945, 1.6898614520910629, 1.7526450904422723),
            dec: Vector3::new(
                1.0825984522657437,
                0.943_679_018_934_623_1,
                0.827_517_321_571_201_4,
            ),
            time: Vector3::new(
                57_028.454_047_592_59,
                57_049.231_857_592_59,
                57_063.959_487_592_59,
            ),
            observer_helio_position: Matrix3::new(
                -0.264_135_633_607_079,
                -0.588_973_552_650_573_5,
                -0.774_192_148_350_372,
                0.869_046_620_910_086,
                0.724_011_718_791_646,
                0.561_510_219_548_918_2,
                0.376_746_685_666_572_5,
                0.313_873_420_677_094,
                0.243_444_791_401_658_5,
            ),
        };

        let binding = gauss.prelim_orbit(env, &IODParams::default()).unwrap();
        let prelim_orbit = binding.get_orbit();

        // This is the expected orbit based on the Orbfit software
        // The values are obtained from the Orbfit output for the same observations
        // The values are very close to the ones obtained from the Rust implementation
        // The floating point differences is very close to one ulp
        // (unit in the last place) of the floating point representation
        let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 57_049.229_045_244_22,
            semi_major_axis: 1.8014943988486352,
            eccentricity: 0.283_514_142_249_080_7,
            inclination: 0.20264170920820326,
            ascending_node_longitude: 8.118_562_444_269_591E-3,
            periapsis_argument: 1.244_795_311_814_302,
            mean_anomaly: 0.44065425435816186,
        });

        // Compare the prelim_orbit with the expected orbit
        assert!(approx_equal(prelim_orbit, &expected_orbit, tol));

        let a = GaussObs {
            idx_obs: [[23, 24, 33]].into(),
            ra: [[1.6893715963476699, 1.689861452091063, 1.7527345385664372]].into(),
            dec: [[1.082468037385525, 0.9436790189346231, 0.8273762407899986]].into(),
            time: [[57028.479297592596, 57049.2318575926, 57063.97711759259]].into(),
            observer_helio_position: [
                [-0.2645666171486676, 0.8689351643673471, 0.3766996211112465],
                [-0.5889735526502539, 0.7240117187952059, 0.3138734206791042],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        let binding = a.prelim_orbit(env, &IODParams::default()).unwrap();
        let prelim_orbit_a = binding.get_orbit();

        let expected_orbit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 57049.22904525282,
            semi_major_axis: 1.801490008178814,
            eccentricity: 0.28350961635625993,
            inclination: 0.20264261257939395,
            ascending_node_longitude: 0.008105552171682476,
            periapsis_argument: 1.244832121745955,
            mean_anomaly: 0.4406444535028061,
        });

        // Compare the prelim_orbit_a with the expected orbit
        assert!(approx_equal(prelim_orbit_a, &expected_orbit, tol));

        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(
                1.6894680552416277,
                1.689_861_821_442_152,
                1.7526488678231147,
            ),
            dec: Vector3::new(
                1.0825994437405373,
                0.943_679_863_334_145,
                0.827_517_360_507_228_6,
            ),
            time: Vector3::new(
                57_028.454_047_592_59,
                57_049.231_857_592_59,
                57_063.959_487_592_59,
            ),
            observer_helio_position: Matrix3::new(
                -0.264_135_633_607_079,
                -0.588_973_552_650_573_5,
                -0.774_192_148_350_372,
                0.869_046_620_910_086,
                0.724_011_718_791_646,
                0.561_510_219_548_918_2,
                0.376_746_685_666_572_5,
                0.313_873_420_677_094,
                0.243_444_791_401_658_5,
            ),
        };

        let binding = gauss.prelim_orbit(env, &IODParams::default()).unwrap();
        let prelim_orbit_b = binding.get_orbit();

        // This is the expected orbit based on the Orbfit software
        // The values are obtained from the Orbfit output for the same observations
        let expected_orbfit = OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch: 57_049.229_045_608_86,
            semi_major_axis: 1.8013098187420686,
            eccentricity: 0.28347096712267805,
            inclination: 0.202_617_665_872_441_2,
            ascending_node_longitude: 8.194_805_420_465_082E-3,
            periapsis_argument: 1.2446747244785052,
            mean_anomaly: 0.44073731381184733,
        });

        // Compare the prelim_orbit_b with the expected orbit
        // The epsilon value is set to 1e-13 for a very close match between the OrbFit and the Rust implementation
        assert!(approx_equal(prelim_orbit_b, &expected_orbfit, tol));
    }

    #[test]
    fn test_orbit_correction() {
        let gauss = GaussObs {
            idx_obs: Vector3::new(0, 1, 2),
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(
                1.082_468_037_385_525,
                0.943_580_504_794_621_6,
                0.827_376_240_789_998_6,
            ),
            time: Vector3::new(
                57028.479297592596,
                57_049.245_147_592_59,
                57_063.977_117_592_59,
            ),
            observer_helio_position: Matrix3::new(
                -0.26456661713915464,
                0.868_935_164_369_495,
                0.376_699_621_109_192_2,
                -0.589_163_185_217_412_7,
                0.723_887_251_679_477_7,
                0.313_818_651_652_458_5,
                -0.774_387_443_796_959_6,
                0.561_288_470_926_116_4,
                0.24334971075289916,
            )
            .transpose(),
        };

        let ast_pos = Matrix3::new(
            -0.288_119_690_673_496,
            1.0666372979405205,
            0.751_481_548_179_728_3,
            -0.623_550_051_003_163_9,
            1.011_260_185_597_692,
            0.713_100_363_506_241_4,
            -0.844_585_047_518_766_5,
            0.942_853_945_425_542_1,
            0.665_339_154_117_050_1,
        )
        .transpose();

        let ast_vel = Vector3::new(
            -1.554_984_513_777_466_3E-2,
            -3.876_936_109_837_657_7E-3,
            -2.701_407_400_297_996_4E-3,
        );

        let unit_matrix = Matrix3::new(
            -5.549_934_652_247_514E-2,
            0.46585594034226024,
            0.883_118_375_634_550_3,
            -6.972_979_004_485_365E-2,
            0.582_735_701_227_938_9,
            0.809_664_658_296_682_1,
            -0.12245931009139571,
            0.665_638_743_839_060_6,
            0.736_158_121_650_706_8,
        )
        .transpose();

        let inv_unit_matrix = Matrix3::new(
            -18.774792915974604,
            41.814_279_122_702_03,
            -23.466669573973437,
            -8.164_790_710_343_109,
            11.489343729350425,
            -2.8418335594428177,
            4.259_482_782_736_117,
            -3.432_964_304_649_721,
            2.434_579_475_328_179E-2,
        )
        .transpose();

        let (new_ast_pos, new_ast_vel, corrected_epoch) = gauss
            .pos_and_vel_correction(
                &IODParams::default(),
                &ast_pos,
                &ast_vel,
                &unit_matrix,
                &inv_unit_matrix,
                1e3,
                1.,
                1e-10,
                50,
            )
            .unwrap();

        assert_eq!(
            new_ast_pos.as_slice(),
            [
                -0.2878540141559046,
                1.06440723593647,
                0.7472540422835181,
                -0.6231216182863288,
                1.0076797497536536,
                0.7081256342111117,
                -0.8435611164802848,
                0.9372882749205874,
                0.6591838430228918
            ]
        );

        assert_eq!(
            new_ast_vel.as_slice(),
            [
                -0.015524309979972159,
                -0.003984105628190921,
                -0.0027640157742952693
            ]
        );

        assert_eq!(corrected_epoch, 57049.24233491307);
    }

    mod test_generate_noisy_realizations {
        use super::*;

        use nalgebra::vector;
        use rand::{rngs::StdRng, SeedableRng};

        #[test]
        fn test_generate_noisy_realizations_count_and_identity() {
            let gauss = GaussObs {
                idx_obs: Vector3::new(0, 1, 2),
                ra: vector![1.0, 2.0, 3.0],
                dec: vector![0.1, 0.2, 0.3],
                time: vector![59000.0, 59001.0, 59002.0],
                observer_helio_position: Matrix3::identity(),
            };

            let errors_ra = vector![1e-6, 1e-6, 1e-6];
            let errors_dec = vector![2e-6, 2e-6, 2e-6];

            let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

            let realizations =
                gauss.generate_noisy_realizations(&errors_ra, &errors_dec, 5, 1.0, &mut rng);

            assert_eq!(realizations.len(), 6); // 1 original + 5 noisy

            // Check the original is preserved
            assert_eq!(realizations[0].ra, gauss.ra);
            assert_eq!(realizations[0].dec, gauss.dec);

            // Check that at least one noisy version differs
            let some_different = realizations[1..]
                .iter()
                .any(|g| g.ra != gauss.ra || g.dec != gauss.dec);
            assert!(
                some_different,
                "Expected at least one noisy realization to differ"
            );
        }

        #[test]
        fn test_generate_noisy_realizations_with_zero_noise() {
            let gauss = GaussObs {
                idx_obs: Vector3::new(0, 1, 2),
                ra: vector![1.5, 1.6, 1.7],
                dec: vector![0.9, 1.0, 1.1],
                time: vector![59000.0, 59001.0, 59002.0],
                observer_helio_position: Matrix3::identity(),
            };

            let errors_ra = vector![0.0, 0.0, 0.0];
            let errors_dec = vector![0.0, 0.0, 0.0];

            let mut rng = StdRng::seed_from_u64(123);

            let realizations =
                gauss.generate_noisy_realizations(&errors_ra, &errors_dec, 3, 1.0, &mut rng);

            for g in realizations {
                assert_eq!(g.ra, gauss.ra);
                assert_eq!(g.dec, gauss.dec);
            }
        }

        #[test]
        fn test_no_noise() {
            let gauss = GaussObs {
                idx_obs: Vector3::new(0, 1, 2),
                ra: vector![1.5, 1.6, 1.7],
                dec: vector![0.9, 1.0, 1.1],
                time: vector![59000.0, 59001.0, 59002.0],
                observer_helio_position: Matrix3::identity(),
            };

            let errors_ra = vector![0.0, 0.0, 0.0];
            let errors_dec = vector![0.0, 0.0, 0.0];

            let mut rng = StdRng::seed_from_u64(123);

            let realizations =
                gauss.generate_noisy_realizations(&errors_ra, &errors_dec, 0, 1.0, &mut rng);

            assert_eq!(realizations.len(), 1); // Only the original observation
            assert_eq!(realizations[0], gauss);
        }

        #[test]
        fn test_generate_noisy_realizations_with_custom_scale() {
            let gauss = GaussObs {
                idx_obs: Vector3::new(0, 1, 2),
                ra: vector![1.0, 1.1, 1.2],
                dec: vector![0.5, 0.6, 0.7],
                time: vector![59000.0, 59001.0, 59002.0],
                observer_helio_position: Matrix3::identity(),
            };

            let errors_ra = vector![1e-6, 1e-6, 1e-6];
            let errors_dec = vector![1e-6, 1e-6, 1e-6];

            let mut rng_low = StdRng::seed_from_u64(42);
            let mut rng_high = StdRng::seed_from_u64(42); // same seed

            let low_noise =
                gauss.generate_noisy_realizations(&errors_ra, &errors_dec, 1, 0.1, &mut rng_low);
            let high_noise =
                gauss.generate_noisy_realizations(&errors_ra, &errors_dec, 1, 10.0, &mut rng_high);

            let diff_low = (low_noise[1].ra - gauss.ra).norm();
            let diff_high = (high_noise[1].ra - gauss.ra).norm();

            assert!(
                diff_high > diff_low,
                "High noise should cause greater deviation than low noise"
            );
        }
    }
}
