#[derive(Debug, PartialEq)]
pub struct KeplerianOrbit {
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub ascending_node_longitude: f64,
    pub periapsis_argument: f64,
    pub mean_anomaly: f64
}
