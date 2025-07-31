use aberth::StopReason;
use nalgebra::Matrix3;
use nalgebra::Vector3;

use crate::constants::Radian;
use crate::constants::{GAUSS_GRAV, VLIGHT_AU};

use crate::initial_orbit_determination::gauss_result::GaussResult;
use crate::kepler::velocity_correction;
use crate::keplerian_element::KeplerianElements;
use crate::observers::helio_obs_pos;
use crate::observers::Observer;
use crate::orb_elem::ccek1;
use crate::orb_elem::eccentricity_control;
use crate::outfit::Outfit;
use crate::outfit_errors::OutfitError;
use crate::ref_system::rotpn;
use crate::ref_system::RefEpoch;
use crate::ref_system::RefSystem;
use aberth::aberth;
use rand::Rng;
use rand_distr::{Distribution, Normal};

/// Observation triplet structure for Gauss's initial orbit determination.
///
/// This struct holds the minimal set of data required to apply the classical Gauss method
/// for orbit determination using three optical astrometric observations taken from the same site.
///
/// Fields
/// -------
/// * `idx_obs`: indices of the three observations in the full dataset (typically used for traceability).
/// * `ra`: right ascensions \[radians\] of the three observations, in the ICRS frame.
/// * `dec`: declinations \[radians\] of the three observations, in the ICRS frame.
/// * `time`: observation times \[MJD TT\], i.e., Modified Julian Dates in Terrestrial Time scale.
/// * `observer_position`: 3×3 matrix of the observer’s heliocentric position vectors at each observation time,
///   with:
///   - columns = observation epochs,
///   - units = astronomical units (AU),
///   - frame = equatorial mean J2000 (same as ICRS).
///
/// Remarks
/// -------
/// * All three observations are assumed to come from the same observing site.
/// * The observer position matrix is computed via [`helio_obs_pos`] during construction.
///
/// # See also
/// * [`GaussObs::with_observer_position`] – constructor that computes this struct from RA/DEC/times
/// * [`helio_obs_pos`] – computes heliocentric observer positions
#[derive(Debug, PartialEq, Clone)]
pub struct GaussObs {
    pub(crate) idx_obs: Vector3<usize>,
    pub(crate) ra: Vector3<Radian>,
    pub(crate) dec: Vector3<Radian>,
    pub(crate) time: Vector3<f64>,
    pub(crate) observer_position: Matrix3<f64>,
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

impl GaussObs {
    /// Construct a `GaussObs` object from observation data and observer position.
    ///
    /// This constructor initializes a `GaussObs` instance using raw observational inputs (RA/DEC/time)
    /// and computes the heliocentric position of the observer at each observation epoch.
    ///
    /// Arguments
    /// ---------
    /// * `state`: reference to the [`Outfit`] context, providing access to JPL ephemerides and Earth orientation.
    /// * `idx_obs`: index triplet identifying the observation triplet used for initial orbit determination.
    /// * `ra`: vector of right ascension values \[radians\] corresponding to the three observations.
    /// * `dec`: vector of declination values \[radians\] corresponding to the three observations.
    /// * `mjd_time`: observation times in Modified Julian Date (TT scale).
    /// * `observer`: reference to the [`Observer`] corresponding to the observing site.
    ///
    /// Returns
    /// --------
    /// * A `GaussObs` struct populated with:
    ///     - the input observation indices, RA/DEC, and times,
    ///     - the observer's heliocentric position matrix computed via [`helio_obs_pos`] (3×3).
    ///
    /// Remarks
    /// -------
    /// * The observer positions are expressed in the equatorial mean J2000 frame (in AU),
    ///   consistent with the requirements of Gauss’s method for orbit determination.
    /// * The order of `ra`, `dec`, `time`, and `observer_position` is consistent: element `i` in each
    ///   vector corresponds to a single astrometric observation.
    ///
    /// # See also
    /// * [`GaussObs`] – data structure for storing triplet-based IOD input
    /// * [`helio_obs_pos`] – used internally to compute heliocentric observer positions
    pub fn with_observer_position(
        state: &Outfit,
        idx_obs: Vector3<usize>,
        ra: Vector3<f64>,
        dec: Vector3<f64>,
        mjd_time: Vector3<f64>,
        observer: [&Observer; 3],
    ) -> Result<GaussObs, OutfitError> {
        let observer_position = helio_obs_pos(observer, &mjd_time, state)?;
        Ok(GaussObs {
            idx_obs,
            ra,
            dec,
            time: mjd_time,
            observer_position,
        })
    }

    /// Generate noisy variants of this `GaussObs` triplet by injecting Gaussian noise
    /// into the right ascension (RA) and declination (DEC) coordinates.
    ///
    /// The function returns a vector of `GaussObs` instances consisting of:
    /// - The original, unperturbed triplet as the first element,
    /// - Followed by `n_realizations` perturbed copies with Gaussian noise applied
    ///   independently to each RA and DEC coordinate, using their corresponding
    ///   standard deviations and a global scaling factor.
    ///
    /// This method is functionally equivalent to OrbFit’s IOD behavior governed by the
    /// `iodnit` and `iodksi` parameters, and is useful for testing the robustness
    /// of initial orbit determination with respect to observational uncertainties.
    ///
    /// # Arguments
    /// ----------
    /// * `errors_ra` - A `Vector3<f64>` containing 1σ errors (in radians) on the RA
    ///   of the three observations.
    /// * `errors_dec` - A `Vector3<f64>` containing 1σ errors (in radians) on the DEC
    ///   of the three observations.
    /// * `n_realizations` - The number of noisy versions to generate (excluding the original).
    /// * `noise_scale` - A multiplicative factor applied to each σ (e.g., 1.0 for nominal noise,
    ///   >1.0 to test instability under exaggerated error).
    /// * `rng` - A mutable reference to a random number generator implementing `rand::Rng`.
    ///
    /// # Returns
    /// ----------
    /// A `Result` containing a vector of `GaussObs` instances:
    /// - `Ok(vec)` if all noise samples were generated successfully,
    /// - `Err(OutfitError::NoiseInjectionError)` if any distribution parameters were invalid.
    ///
    /// # Errors
    /// ----------
    /// Returns an error if any call to `Normal::new` fails due to non-finite or negative variance.
    pub fn generate_noisy_realizations(
        &self,
        errors_ra: &Vector3<f64>,
        errors_dec: &Vector3<f64>,
        n_realizations: usize,
        noise_scale: f64,
        rng: &mut impl Rng,
    ) -> Result<Vec<GaussObs>, OutfitError> {
        let noisy_iter = (0..n_realizations).map(|_| {
            let ra_vec: Vec<f64> = self
                .ra
                .iter()
                .zip(errors_ra.iter())
                .map(|(&ra, &sigma)| {
                    let normal = Normal::new(0.0, sigma * noise_scale)?;
                    Ok(ra + normal.sample(rng))
                })
                .collect::<Result<Vec<f64>, OutfitError>>()?;

            let dec_vec: Vec<f64> = self
                .dec
                .iter()
                .zip(errors_dec.iter())
                .map(|(&dec, &sigma)| {
                    let normal = Normal::new(0.0, sigma * noise_scale)?;
                    Ok(dec + normal.sample(rng))
                })
                .collect::<Result<Vec<f64>, OutfitError>>()?;

            Ok(GaussObs {
                idx_obs: self.idx_obs,
                ra: Vector3::from_column_slice(&ra_vec),
                dec: Vector3::from_column_slice(&dec_vec),
                time: self.time,
                observer_position: self.observer_position,
            })
        });

        std::iter::once(Ok(self.clone()))
            .chain(noisy_iter)
            .collect::<Result<Vec<_>, _>>()
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
        let observer_position_t = self.observer_position.transpose();
        let ra = self.observer_position * vector_a;
        let rb = self.observer_position * vector_b;

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
        unit_matrix: &Matrix3<f64>,
        unit_matinv: &Matrix3<f64>,
        vector_c: &Vector3<f64>,
    ) -> Result<(Matrix3<f64>, f64), OutfitError> {
        let obs_pos_t = self.observer_position;
        let gcap = obs_pos_t * vector_c;
        let crhom = unit_matinv * gcap;
        let rho: Vector3<f64> = -crhom.component_div(vector_c);
        if rho[1] < 0.01 {
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
            .position_vector_and_reference_epoch(unit_matrix, inv_unit_matrix, &vector_c)
            .ok()
        else {
            return None; // root leads to invalid geometry (e.g., ρ₂ < 0.01)
        };

        // Extract the position vector at the central epoch (t2)
        let ast_pos_second_time = ast_pos_all_time.column(1).into();

        // Estimate the velocity vector at the second observation via Gibbs correction
        let asteroid_vel = self.gibbs_correction(&ast_pos_all_time, tau1, tau3);

        // Apply eccentricity and perihelion distance control
        let (is_accepted, _, _, _) =
            eccentricity_control(&ast_pos_second_time, &asteroid_vel, 1e3, 5.0)?;

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
    /// # See also
    /// * [`rotpn`] – computes the rotation matrix between celestial reference frames.
    /// * [`ccek1`] – converts position and velocity vectors to orbital elements.
    /// * [`KeplerianElements`] – definition of the orbital elements struct.
    fn compute_orbit_from_state(
        &self,
        &asteroid_position: &Vector3<f64>,
        &asteroid_velocity: &Vector3<f64>,
        reference_epoch: f64,
    ) -> Result<KeplerianElements, OutfitError> {
        // Compute the rotation matrix from equatorial mean J2000 to ecliptic mean J2000
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2)?;

        // Apply the transformation to position and velocity vectors
        let matrix_elc_transform = Matrix3::from(roteqec).transpose();
        let ecl_pos = matrix_elc_transform * asteroid_position;
        let ecl_vel = matrix_elc_transform * asteroid_velocity;

        // Merge position and velocity into a single array expected by `ccek1`
        let ast_pos_vel: [f64; 6] = [
            ecl_pos.x, ecl_pos.y, ecl_pos.z, ecl_vel.x, ecl_vel.y, ecl_vel.z,
        ];

        // Prepare buffers for orbital element output
        let mut elem = [0.0; 6];
        let mut type_ = String::new();

        // Compute the classical orbital elements from the state vector
        ccek1(&mut elem, &mut type_, &ast_pos_vel);

        // Return orbital elements in a structured form
        Ok(KeplerianElements {
            reference_epoch,
            semi_major_axis: elem[0],
            eccentricity: elem[1],
            inclination: elem[2],
            ascending_node_longitude: elem[3],
            periapsis_argument: elem[4],
            mean_anomaly: elem[5],
        })
    }

    /// Estimate an initial orbit using Gauss’s method from three astrometric observations.
    ///
    /// This method performs a full Gauss orbit determination on a triplet of observations (`GaussObs`),
    /// returning either a preliminary or a corrected orbital solution. The process includes:
    ///
    /// - Computation of direction vectors and time intervals,
    /// - Construction of the 8th-degree polynomial for the topocentric distance,
    /// - Root solving using the Aberth method,
    /// - Evaluation and filtering of roots based on physical criteria (e.g. eccentricity),
    /// - Optional iterative refinement of the orbit via velocity correction.
    ///
    /// Returns
    /// --------
    /// * `Ok(GaussResult::PrelimOrbit)` if a valid solution is found but the correction step fails or is skipped.
    /// * `Ok(GaussResult::CorrectedOrbit)` if the solution converges after correction.
    /// * `Err(OutfitError)` if no valid root is found or the direction matrix is singular.
    ///
    /// # See also
    /// * [`GaussObs`] – Struct holding the observation triplet used in this computation.
    /// * [`GaussResult`] – Enum indicating whether the orbit is preliminary or corrected.
    /// * [`GaussObs::pos_and_vel_correction`] – Performs the refinement step on the orbital state.
    /// * [`GaussObs::solve_8poly`] – Finds valid roots of the 8th-degree distance polynomial.
    pub fn prelim_orbit(&self) -> Result<GaussResult, OutfitError> {
        // Compute core quantities: time intervals, direction matrix and its inverse, and coefficients
        let (tau1, tau3, unit_matrix, inv_unit_matrix, vect_a, vect_b) = self.gauss_prelim()?;

        // Build the coefficients of the 8th-degree polynomial in r2
        let (coeff_6, coeff_3, coeff_0) =
            self.coeff_eight_poly(&unit_matrix, &inv_unit_matrix, &vect_a, &vect_b);
        let polynomial = [coeff_0, 0.0, 0.0, coeff_3, 0.0, 0.0, coeff_6, 0.0, 1.0];

        // Solve for real, positive roots using the Aberth method
        let roots = self.solve_8poly(&polynomial, 50, 1e-6, 1e-6)?;

        // Try to find the first physically valid root (rho2 > 0.01, eccentricity bounds, etc.)
        let Some((asteroid_pos_all_time, asteroid_vel, reference_epoch)) =
            roots.into_iter().find_map(|root| {
                self.accept_root(
                    root,
                    &unit_matrix,
                    &inv_unit_matrix,
                    &vect_a,
                    &vect_b,
                    tau1,
                    tau3,
                )
            })
        else {
            return Err(OutfitError::GaussNoRootsFound);
        };

        // Try to refine the orbital state through iterative velocity correction
        let Some((corrected_pos, corrected_vel, corrected_epoch)) = self.pos_and_vel_correction(
            &asteroid_pos_all_time,
            &asteroid_vel,
            &unit_matrix,
            &inv_unit_matrix,
            1e3,   // perihelion max
            5.0,   // eccentricity max
            1e-10, // convergence threshold
            50,    // max iterations
        ) else {
            // If correction failed or diverged, return preliminary orbit
            return Ok(GaussResult::PrelimOrbit(self.compute_orbit_from_state(
                &asteroid_pos_all_time.column(1).into(),
                &asteroid_vel,
                reference_epoch,
            )?));
        };

        // Correction succeeded; return refined orbit
        Ok(GaussResult::CorrectedOrbit(self.compute_orbit_from_state(
            &corrected_pos.column(1).into(),
            &corrected_vel,
            corrected_epoch,
        )?))
    }

    /// Iteratively refine the asteroid's velocity and position at the central observation time.
    ///
    /// This function performs the correction loop of Gauss’s method by adjusting the velocity
    /// estimate using Lagrange coefficients and recomputing the corresponding positions.
    /// The iterations continue until the position correction converges or until the iteration
    /// limit is reached.
    ///
    /// Arguments
    /// ---------
    /// * `asteroid_position` – Initial 3×3 matrix of asteroid positions at each observation time (in AU).
    /// * `asteroid_velocity` – Initial velocity vector at the central time (in AU/day).
    /// * `unit_matrix` – 3×3 matrix of line-of-sight unit vectors from observer to object.
    /// * `inv_unit_matrix` – Inverse of the unit direction matrix.
    /// * `peri_max` – Maximum allowed perihelion distance (used in eccentricity filtering).
    /// * `ecc_max` – Maximum allowed orbital eccentricity.
    /// * `err_max` – Relative convergence threshold on position vectors.
    /// * `itmax` – Maximum number of iterations allowed.
    ///
    /// Returns
    /// --------
    /// * `Some((corrected_positions, corrected_velocity, corrected_epoch))` if convergence succeeded
    /// * `None` if the solution failed due to divergence or non-physical orbits (e.g., hyperbolic).
    ///
    /// Remarks
    /// --------
    /// * Each iteration updates the velocity using Lagrange coefficient matching (`velocity_correction`)
    ///   and recomputes the associated topocentric distances and positions.
    /// * The updated position matrix is compared to the previous one using relative RMS error.
    /// * This routine corresponds to the velocity refinement loop in OrbFit’s `gaussn.f90`.
    ///
    /// # See also
    /// * [`velocity_correction`] – Lagrange-based velocity update.
    /// * [`eccentricity_control`] – Filters solutions based on dynamic acceptability.
    /// * [`GaussObs::position_vector_and_reference_epoch`] – Recomputes full asteroid positions at each epoch.
    #[allow(clippy::too_many_arguments)]
    pub fn pos_and_vel_correction(
        &self,
        asteroid_position: &Matrix3<f64>,
        asteroid_velocity: &Vector3<f64>,
        unit_matrix: &Matrix3<f64>,
        inv_unit_matrix: &Matrix3<f64>,
        peri_max: f64,
        ecc_max: f64,
        err_max: f64,
        itmax: i32,
    ) -> Option<(Matrix3<f64>, Vector3<f64>, f64)> {
        // Initialize state with uncorrected values
        let mut previous_ast_pos = *asteroid_position;
        let mut previous_ast_vel = *asteroid_velocity;
        let mut corrected_epoch = 0.0;

        // Main correction loop
        for _ in 0..itmax {
            // Compute velocity correction using first observation
            let Ok((ast_vel1, f1, g1)) = velocity_correction(
                &previous_ast_pos.column(0).into(),
                &previous_ast_pos.column(1).into(),
                &previous_ast_vel,
                self.time[0] - self.time[1],
                peri_max,
                ecc_max,
            ) else {
                continue; // skip iteration if invalid geometry
            };

            // Compute velocity correction using third observation
            let Ok((ast_vel2, f2, g2)) = velocity_correction(
                &previous_ast_pos.column(2).into(),
                &previous_ast_pos.column(1).into(),
                &previous_ast_vel,
                self.time[2] - self.time[1],
                peri_max,
                ecc_max,
            ) else {
                continue; // skip iteration if invalid geometry
            };

            // Average the two velocity estimates
            let corrected_velocity = (ast_vel1 + ast_vel2) / 2.0;

            // Compute updated vector C (Gauss coefficients)
            let f_lagrange = f1 * g2 - f2 * g1;
            let c_vector = Vector3::new(g2 / f_lagrange, -1.0, -(g1 / f_lagrange));

            // Recompute asteroid positions and light-time-corrected epoch
            let Some((new_ast_pos, it_epoch)) = self
                .position_vector_and_reference_epoch(unit_matrix, inv_unit_matrix, &c_vector)
                .ok()
            else {
                continue; // skip if positions are invalid (e.g., rho2 too small)
            };

            // Compute relative error between new and previous position matrices
            let diff_pos = new_ast_pos - previous_ast_pos;
            let diff_pos_squared = diff_pos.component_mul(&diff_pos).sum();
            let sca = new_ast_pos.component_mul(&new_ast_pos).sum();
            let ast_pos_err = (diff_pos_squared / sca).sqrt();

            // Check orbital acceptability (eccentricity, energy, perihelion)
            let Some((is_accepted, _, _, _)) =
                eccentricity_control(&new_ast_pos.column(1).into(), &corrected_velocity, 1e3, 5.0)
            else {
                return None; // orbit is dynamically invalid
            };

            if !is_accepted {
                return None;
            }

            // Update state for next iteration or final output
            previous_ast_pos = new_ast_pos;
            previous_ast_vel = corrected_velocity;
            corrected_epoch = it_epoch;

            // Check convergence
            if ast_pos_err <= err_max {
                break;
            }
        }

        // Return the refined state if successful
        Some((previous_ast_pos, previous_ast_vel, corrected_epoch))
    }
}

#[cfg(test)]
pub(crate) mod gauss_test {

    use super::*;
    use crate::keplerian_element::test_keplerian_element::assert_orbit_close;

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
            observer_position: Matrix3::zeros(),
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
            observer_position: Matrix3::new(
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
            observer_position: Matrix3::zeros(),
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
            observer_position: Matrix3::new(
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
        let ast_pos =
            gauss.position_vector_and_reference_epoch(&unit_matrix, &inv_unit_matrix, &vector_c);
        assert!(ast_pos.is_err());

        let second_root: f64 = 1.3856312487504951;
        let r2m3 = 1. / second_root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vector_a[0] + vector_b[0] * r2m3,
            -1.,
            vector_a[2] + vector_b[2] * r2m3,
        );
        let (ast_pos, ref_epoch) = gauss
            .position_vector_and_reference_epoch(&unit_matrix, &inv_unit_matrix, &vector_c)
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
            observer_position: Matrix3::zeros(),
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
            observer_position: Matrix3::new(
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

        let prelim_orbit = gauss.prelim_orbit().unwrap();

        // This is the expected orbit based on the Orbfit software
        // The values are obtained from the Orbfit output for the same observations
        // The values are very close to the ones obtained from the Rust implementation
        // The floating point differences is very close to one ulp
        // (unit in the last place) of the floating point representation
        let expected_orbit = KeplerianElements {
            reference_epoch: 57_049.229_045_244_22,
            semi_major_axis: 1.8014943988486352,
            eccentricity: 0.283_514_142_249_080_7,
            inclination: 0.20264170920820326,
            ascending_node_longitude: 8.118_562_444_269_591E-3,
            periapsis_argument: 1.244_795_311_814_302,
            mean_anomaly: 0.44065425435816186,
        };

        // Compare the prelim_orbit with the expected orbit
        assert_orbit_close(&prelim_orbit, &expected_orbit, 1e-14);

        let a = GaussObs {
            idx_obs: [[23, 24, 33]].into(),
            ra: [[1.6893715963476699, 1.689861452091063, 1.7527345385664372]].into(),
            dec: [[1.082468037385525, 0.9436790189346231, 0.8273762407899986]].into(),
            time: [[57028.479297592596, 57049.2318575926, 57063.97711759259]].into(),
            observer_position: [
                [-0.2645666171486676, 0.8689351643673471, 0.3766996211112465],
                [-0.5889735526502539, 0.7240117187952059, 0.3138734206791042],
                [-0.7743874438017259, 0.5612884709246775, 0.2433497107566823],
            ]
            .into(),
        };

        let prelim_orbit_a = a.prelim_orbit().unwrap();

        let expected_orbit = KeplerianElements {
            reference_epoch: 57049.22904525282,
            semi_major_axis: 1.801490008178814,
            eccentricity: 0.28350961635625993,
            inclination: 0.20264261257939395,
            ascending_node_longitude: 0.008105552171682476,
            periapsis_argument: 1.244832121745955,
            mean_anomaly: 0.4406444535028061,
        };

        // Compare the prelim_orbit_a with the expected orbit
        assert_orbit_close(&prelim_orbit_a, &expected_orbit, 1e-13);

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
            observer_position: Matrix3::new(
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

        let prelim_orbit_b = gauss.prelim_orbit().unwrap();

        // This is the expected orbit based on the Orbfit software
        // The values are obtained from the Orbfit output for the same observations
        let expected_orbfit = KeplerianElements {
            reference_epoch: 57_049.229_045_608_86,
            semi_major_axis: 1.8013098187420686,
            eccentricity: 0.28347096712267805,
            inclination: 0.202_617_665_872_441_2,
            ascending_node_longitude: 8.194_805_420_465_082E-3,
            periapsis_argument: 1.2446747244785052,
            mean_anomaly: 0.44073731381184733,
        };

        // Compare the prelim_orbit_b with the expected orbit
        // The epsilon value is set to 1e-10 for a very close match between the OrbFit and the Rust implementation
        assert_orbit_close(&prelim_orbit_b, &expected_orbfit, 1e-13);
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
            observer_position: Matrix3::new(
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
                observer_position: Matrix3::identity(),
            };

            let errors_ra = vector![1e-6, 1e-6, 1e-6];
            let errors_dec = vector![2e-6, 2e-6, 2e-6];

            let mut rng = StdRng::seed_from_u64(42_u64); // seed for reproducibility

            let realizations = gauss
                .generate_noisy_realizations(&errors_ra, &errors_dec, 5, 1.0, &mut rng)
                .unwrap();

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
                observer_position: Matrix3::identity(),
            };

            let errors_ra = vector![0.0, 0.0, 0.0];
            let errors_dec = vector![0.0, 0.0, 0.0];

            let mut rng = StdRng::seed_from_u64(123);

            let realizations = gauss
                .generate_noisy_realizations(&errors_ra, &errors_dec, 3, 1.0, &mut rng)
                .unwrap();

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
                observer_position: Matrix3::identity(),
            };

            let errors_ra = vector![0.0, 0.0, 0.0];
            let errors_dec = vector![0.0, 0.0, 0.0];

            let mut rng = StdRng::seed_from_u64(123);

            let realizations = gauss
                .generate_noisy_realizations(&errors_ra, &errors_dec, 0, 1.0, &mut rng)
                .unwrap();

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
                observer_position: Matrix3::identity(),
            };

            let errors_ra = vector![1e-6, 1e-6, 1e-6];
            let errors_dec = vector![1e-6, 1e-6, 1e-6];

            let mut rng_low = StdRng::seed_from_u64(42);
            let mut rng_high = StdRng::seed_from_u64(42); // same seed

            let low_noise = gauss
                .generate_noisy_realizations(&errors_ra, &errors_dec, 1, 0.1, &mut rng_low)
                .unwrap();
            let high_noise = gauss
                .generate_noisy_realizations(&errors_ra, &errors_dec, 1, 10.0, &mut rng_high)
                .unwrap();

            let diff_low = (low_noise[1].ra - gauss.ra).norm();
            let diff_high = (high_noise[1].ra - gauss.ra).norm();

            assert!(
                diff_high > diff_low,
                "High noise should cause greater deviation than low noise"
            );
        }
    }
}
