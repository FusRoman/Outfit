use approx::abs_diff_eq;
use outfit::{orbit_type::uncertainty::OrbitalCovariance, OrbitalElements};

#[allow(dead_code)]
pub fn approx_equal(current: &OrbitalElements, other: &OrbitalElements, tol: f64) -> bool {
    match (current, other) {
        (
            OrbitalElements::Keplerian {
                elements: ke1,
                uncertainty: unc1,
                covariance: cov1,
            },
            OrbitalElements::Keplerian {
                elements: ke2,
                uncertainty: unc2,
                covariance: cov2,
            },
        ) => {
            let elements_eq = abs_diff_eq!(ke1.semi_major_axis, ke2.semi_major_axis, epsilon = tol)
                && abs_diff_eq!(ke1.eccentricity, ke2.eccentricity, epsilon = tol)
                && abs_diff_eq!(ke1.inclination, ke2.inclination, epsilon = tol)
                && abs_diff_eq!(
                    ke1.ascending_node_longitude,
                    ke2.ascending_node_longitude,
                    epsilon = tol
                )
                && abs_diff_eq!(
                    ke1.periapsis_argument,
                    ke2.periapsis_argument,
                    epsilon = tol
                )
                && abs_diff_eq!(ke1.mean_anomaly, ke2.mean_anomaly, epsilon = tol);

            let uncertainty_eq = match (unc1, unc2) {
                (None, None) => true,
                (Some(u1), Some(u2)) => {
                    abs_diff_eq!(u1.semi_major_axis, u2.semi_major_axis, epsilon = tol)
                        && abs_diff_eq!(u1.eccentricity, u2.eccentricity, epsilon = tol)
                        && abs_diff_eq!(u1.inclination, u2.inclination, epsilon = tol)
                        && abs_diff_eq!(
                            u1.ascending_node_longitude,
                            u2.ascending_node_longitude,
                            epsilon = tol
                        )
                        && abs_diff_eq!(u1.periapsis_argument, u2.periapsis_argument, epsilon = tol)
                        && abs_diff_eq!(u1.mean_anomaly, u2.mean_anomaly, epsilon = tol)
                }
                _ => false,
            };

            let covariance_eq = approx_equal_covariance(cov1, cov2, tol);

            elements_eq && uncertainty_eq && covariance_eq
        }

        (
            OrbitalElements::Equinoctial {
                elements: ee1,
                uncertainty: unc1,
                covariance: cov1,
            },
            OrbitalElements::Equinoctial {
                elements: ee2,
                uncertainty: unc2,
                covariance: cov2,
            },
        ) => {
            let elements_eq = abs_diff_eq!(ee1.semi_major_axis, ee2.semi_major_axis, epsilon = tol)
                && abs_diff_eq!(
                    ee1.eccentricity_sin_lon,
                    ee2.eccentricity_sin_lon,
                    epsilon = tol
                )
                && abs_diff_eq!(
                    ee1.eccentricity_cos_lon,
                    ee2.eccentricity_cos_lon,
                    epsilon = tol
                )
                && abs_diff_eq!(
                    ee1.tan_half_incl_sin_node,
                    ee2.tan_half_incl_sin_node,
                    epsilon = tol
                )
                && abs_diff_eq!(
                    ee1.tan_half_incl_cos_node,
                    ee2.tan_half_incl_cos_node,
                    epsilon = tol
                )
                && abs_diff_eq!(ee1.mean_longitude, ee2.mean_longitude, epsilon = tol);

            let uncertainty_eq = match (unc1, unc2) {
                (None, None) => true,
                (Some(u1), Some(u2)) => {
                    abs_diff_eq!(u1.semi_major_axis, u2.semi_major_axis, epsilon = tol)
                        && abs_diff_eq!(
                            u1.eccentricity_sin_lon,
                            u2.eccentricity_sin_lon,
                            epsilon = tol
                        )
                        && abs_diff_eq!(
                            u1.eccentricity_cos_lon,
                            u2.eccentricity_cos_lon,
                            epsilon = tol
                        )
                        && abs_diff_eq!(
                            u1.tan_half_incl_sin_node,
                            u2.tan_half_incl_sin_node,
                            epsilon = tol
                        )
                        && abs_diff_eq!(
                            u1.tan_half_incl_cos_node,
                            u2.tan_half_incl_cos_node,
                            epsilon = tol
                        )
                        && abs_diff_eq!(u1.mean_longitude, u2.mean_longitude, epsilon = tol)
                }
                _ => false,
            };

            let covariance_eq = approx_equal_covariance(cov1, cov2, tol);

            elements_eq && uncertainty_eq && covariance_eq
        }

        (
            OrbitalElements::Cometary {
                elements: ce1,
                uncertainty: unc1,
                covariance: cov1,
            },
            OrbitalElements::Cometary {
                elements: ce2,
                uncertainty: unc2,
                covariance: cov2,
            },
        ) => {
            let elements_eq = abs_diff_eq!(
                ce1.perihelion_distance,
                ce2.perihelion_distance,
                epsilon = tol
            ) && abs_diff_eq!(ce1.eccentricity, ce2.eccentricity, epsilon = tol)
                && abs_diff_eq!(ce1.inclination, ce2.inclination, epsilon = tol)
                && abs_diff_eq!(
                    ce1.ascending_node_longitude,
                    ce2.ascending_node_longitude,
                    epsilon = tol
                )
                && abs_diff_eq!(
                    ce1.periapsis_argument,
                    ce2.periapsis_argument,
                    epsilon = tol
                )
                && abs_diff_eq!(ce1.true_anomaly, ce2.true_anomaly, epsilon = tol);

            let uncertainty_eq = match (unc1, unc2) {
                (None, None) => true,
                (Some(u1), Some(u2)) => {
                    abs_diff_eq!(
                        u1.perihelion_distance,
                        u2.perihelion_distance,
                        epsilon = tol
                    ) && abs_diff_eq!(u1.eccentricity, u2.eccentricity, epsilon = tol)
                        && abs_diff_eq!(u1.inclination, u2.inclination, epsilon = tol)
                        && abs_diff_eq!(
                            u1.ascending_node_longitude,
                            u2.ascending_node_longitude,
                            epsilon = tol
                        )
                        && abs_diff_eq!(u1.periapsis_argument, u2.periapsis_argument, epsilon = tol)
                        && abs_diff_eq!(u1.true_anomaly, u2.true_anomaly, epsilon = tol)
                }
                _ => false,
            };

            let covariance_eq = approx_equal_covariance(cov1, cov2, tol);

            elements_eq && uncertainty_eq && covariance_eq
        }

        _ => false,
    }
}

/// Compare two optional [`OrbitalCovariance`] matrices entry-wise within `tol`.
///
/// * `(None, None)` → `true`
/// * `(Some, Some)` → all 36 entries within `tol`
/// * mixed         → `false`
#[allow(dead_code)]
fn approx_equal_covariance(
    cov1: &Option<OrbitalCovariance>,
    cov2: &Option<OrbitalCovariance>,
    tol: f64,
) -> bool {
    match (cov1, cov2) {
        (None, None) => true,
        (Some(c1), Some(c2)) => c1
            .matrix
            .iter()
            .zip(c2.matrix.iter())
            .all(|(a, b)| abs_diff_eq!(a, b, epsilon = tol)),
        _ => false,
    }
}
