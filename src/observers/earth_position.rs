use crate::{
    jpl_ephem::jpl_request::{deserialize_vector, request_vector, PosRecord},
    outfit::Outfit,
};

/// Request the JPL Horizon API to get the position vector of Earth
/// with respect to the Sun at different time.
///
/// Argument
/// --------
/// * mjd_list: a vector of date in modified julian date format (MJD)
///
/// Return
/// ------
/// * a vector of PosRecord, the position vector component are in astronomical units
pub(in crate::observers) fn earth_position_from_jpl_horizon(
    mjd_list: &Vec<f64>,
    env_state: &Outfit,
) -> Vec<PosRecord> {
    let response_data = request_vector(mjd_list, env_state);
    deserialize_vector(&response_data)
}

#[cfg(test)]
mod earth_pos_tests {
    use crate::{error_models::ErrorModel, time::date_to_mjd};

    use super::*;
    use crate::outfit::Outfit;

    #[test]
    #[ignore]
    fn test_get_earth_pos() {
        let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let jd_list = date_to_mjd(&date_list);
        let earth_vector = earth_position_from_jpl_horizon(&jd_list, &state);
        assert_eq!(
            earth_vector,
            vec![
                PosRecord {
                    jd: 2459400.032916666,
                    date: "A.D.2021-Jul-0412:47:24.0000".into(),
                    x: 0.2195672929244244,
                    y: -0.9108330730147444,
                    z: -0.3948423288985838
                },
                PosRecord {
                    jd: 2460672.574629629,
                    date: "A.D.2024-Dec-2801:47:28.0000".into(),
                    x: -0.1107728032684787,
                    y: 0.8965650072539966,
                    z: 0.388651757715346
                }
            ]
        );
    }

    #[test]
    #[ignore]
    fn test_earth_pos_with_mjd() {
        let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let test_mjd = vec![57028.479297592596, 57_049.245_147_592_59, 57_063.977_117_592_59];
        let earth_vector = earth_position_from_jpl_horizon(&test_mjd, &state);
        assert_eq!(
            earth_vector,
            vec![
                PosRecord {
                    jd: 2457028.979297593,
                    date: "A.D.2015-Jan-0611:30:11.3120".into(),
                    x: -0.2645457550344549,
                    y: 0.8689011473732227,
                    z: 0.3766846006816814
                },
                PosRecord {
                    jd: 2457049.745147592,
                    date: "A.D.2015-Jan-2705:53:00.7520".into(),
                    x: -0.5891845254273143,
                    y: 0.7238535049291953,
                    z: 0.3138036934112774
                },
                PosRecord {
                    jd: 2457064.477117592,
                    date: "A.D.2015-Feb-1023:27:02.9600".into(),
                    x: -0.7743644470250609,
                    y: 0.5612696626811766,
                    z: 0.2433192212445796
                }
            ]
        )
    }
}
