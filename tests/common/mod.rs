use approx::abs_diff_eq;
use outfit::orbit_type::OrbitalElements;

pub fn approx_equal(current: &OrbitalElements, other: &OrbitalElements, tol: f64) -> bool {
    match (current, other) {
        (OrbitalElements::Keplerian(ke1), OrbitalElements::Keplerian(ke2)) => {
            abs_diff_eq!(ke1.semi_major_axis, ke2.semi_major_axis, epsilon = tol)
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
                && abs_diff_eq!(ke1.mean_anomaly, ke2.mean_anomaly, epsilon = tol)
        }
        (OrbitalElements::Equinoctial(ee1), OrbitalElements::Equinoctial(ee2)) => {
            abs_diff_eq!(ee1.semi_major_axis, 0.0, epsilon = tol)
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
                && abs_diff_eq!(ee1.mean_longitude, ee2.mean_longitude, epsilon = tol)
        }
        (OrbitalElements::Cometary(ce1), OrbitalElements::Cometary(ce2)) => {
            abs_diff_eq!(
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
                && abs_diff_eq!(ce1.true_anomaly, ce2.true_anomaly, epsilon = tol)
        }
        _ => false, // Different types cannot be equal
    }
}
