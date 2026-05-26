use nalgebra::Matrix2;
use photom::observation_dataset::observation::Observation;

pub mod least_square;
pub mod obs_dataset_api;

/// Returns the 2x2 weight matrix for a given observation, based on its reported RA and Dec errors.
/// If the observation is rejected, returns a zero matrix.
///
/// The weight matrix is diagonal, with entries equal to the inverse of the squared errors in RA and Dec.
///
/// # Arguments
///
/// - `obs`: The observation for which to compute the weight matrix.
/// - `rejected`: A boolean indicating whether the observation is rejected. If true, the weight matrix will be zero, effectively excluding the observation from orbit determination computations.
///
/// # Returns
///
/// A 2x2 diagonal matrix where the diagonal entries are the inverse of the squared RA and Dec errors, or a zero matrix if the observation is rejected.
pub fn observation_weight(obs: &Observation, rejected: bool) -> Matrix2<f64> {
    if rejected {
        return Matrix2::zeros();
    }

    let sa2 = obs.equ_coord().ra_error.powf(2.0);
    let sd2 = obs.equ_coord().dec_error.powf(2.0);
    Matrix2::new(1.0 / sa2, 0.0, 0.0, 1.0 / sd2)
}
