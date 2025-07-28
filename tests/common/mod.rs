use approx::assert_relative_eq;
use outfit::keplerian_element::KeplerianElements;

pub fn assert_orbit_close(actual: &KeplerianElements, expected: &KeplerianElements, epsilon: f64) {
    assert_relative_eq!(
        actual.reference_epoch,
        expected.reference_epoch,
        epsilon = epsilon
    );
    assert_relative_eq!(
        actual.semi_major_axis,
        expected.semi_major_axis,
        epsilon = epsilon
    );
    assert_relative_eq!(
        actual.eccentricity,
        expected.eccentricity,
        epsilon = epsilon
    );
    assert_relative_eq!(actual.inclination, expected.inclination, epsilon = epsilon);
    assert_relative_eq!(
        actual.ascending_node_longitude,
        expected.ascending_node_longitude,
        epsilon = epsilon
    );
    assert_relative_eq!(
        actual.periapsis_argument,
        expected.periapsis_argument,
        epsilon = epsilon
    );
    assert_relative_eq!(
        actual.mean_anomaly,
        expected.mean_anomaly,
        epsilon = epsilon
    );
}
