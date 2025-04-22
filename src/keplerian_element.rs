
/// Keplerian orbital elements
/// Units:
/// * `reference_epoch`: MJD (Modified Julian Date)
/// * `semi_major_axis`: AU (Astronomical Units)
/// * `eccentricity`: unitless
/// * `inclination`: degrees
/// * `ascending_node_longitude`: degrees
/// * `periapsis_argument`: degrees
/// * `mean_anomaly`: degrees
#[derive(Debug, PartialEq)]
pub struct KeplerianElements {
    pub reference_epoch: f64,
    pub semi_major_axis: f64,
    pub eccentricity: f64,
    pub inclination: f64,
    pub ascending_node_longitude: f64,
    pub periapsis_argument: f64,
    pub mean_anomaly: f64
}
