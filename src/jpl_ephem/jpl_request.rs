use itertools::Itertools;
use nalgebra::Vector3;
use regex::Regex;

use crate::outfit::Outfit;

/// Contains the informations from the JPL Horizons vector state query
/// x,y,z are the components of the position vector at the time contained in the jd and date field
#[derive(Debug, serde::Deserialize, PartialEq)]
pub(crate) struct PosRecord {
    #[serde(rename = "JDTDB")]
    pub(crate) jd: f64,
    #[serde(rename = "CalendarDate(TDB)")]
    pub(crate) date: String,
    #[serde(rename = "X")]
    pub(crate) x: f64, // km
    #[serde(rename = "Y")]
    pub(crate) y: f64, // km
    #[serde(rename = "Z")]
    pub(crate) z: f64, // km
}

impl PosRecord {
    pub fn pos_vector(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }
}

fn jd_tlist(mjd_list: &Vec<f64>) -> String {
    mjd_list.iter().join(",")
}

/// Request the JPL HORIZON API to get the Earth position vector with respect to the Sun.
///
/// Argument
/// --------
/// * mjd_list: a list of date in modified julian date format (MJD)
///
/// Return
/// ------
/// * The JPL API raw response
pub(crate) fn request_vector(mjd_list: &Vec<f64>, env_state: &Outfit) -> String {
    let requested_params = format!(
        "
!$$SOF
COMMAND='399'
OBJ_DATA='NO'
MAKE_EPHEM='YES'
TABLE_TYPE='VECTORS'
CENTER='500@10'
TLIST_TYPE=MJD
TLIST={}
CSV_FORMAT=YES
REF_SYSTEM=ICRF
OUT_UNITS=AU-D
REF_PLANE=FRAME
VEC_TABLE=1
",
        jd_tlist(mjd_list)
    );
    env_state.post_url(
        "https://ssd.jpl.nasa.gov/api/horizons_file.api",
        [("format", "text"), ("input", &requested_params)],
    )
}

/// Parse the JPL raw response and return a vector of PosRecord
/// containg the vector component of the Earth position vector
///
/// Argument
/// --------
/// * jpl_response: the raw JPL response from the API
///
/// Return
/// ------
/// * a vector of PosRecord
pub(crate) fn deserialize_vector(jpl_response: &String) -> Vec<PosRecord> {
    // regex to match the data part of the jpl horizon response
    let data_regex = Regex::new(r"\$\$SOE\n([^]]*),\n\$\$EOE").unwrap();
    // regex to match the header part of the jpl horizon response
    let header_regex = Regex::new(r"ICRF\n\*{79}([^]]*)\*{122}\n\$\$SOE").unwrap();

    let error_msg_data = format!(
        "JPL deserializer: Error when matching the data regex\njpl_response:\n{}",
        jpl_response
    );

    let error_msg_header = format!(
        "JPL deserializer: Error when matching the header regex\njpl_response:\n{}",
        &jpl_response
    );

    let match_data = data_regex
        .captures(&jpl_response)
        .expect(&error_msg_data)
        .get(1)
        .expect("JPL deserializer: Error when trying to get the match at index 1 for data regex")
        .as_str()
        .replace(" ", "")
        .replace(",\n", "\n");
    let mut match_header = header_regex
        .captures(&jpl_response)
        .expect(&error_msg_header)
        .get(1)
        .expect("JPL deserializer: Error when trying to get the match at index 1 for header regex")
        .as_str()
        .replace(" ", "");

    // remove the last two characters of the headers corresponding to '\n' and ','
    // The last ',' makes the CSV deserializer think there is an empty column at the end.
    match_header.pop();
    match_header.pop();

    // reconstruct the csv string by adding the header above the data.
    let data = format!("{}\n{}", match_header, match_data);

    // Use serde to deserialize the string into a vector of PosRecord.
    let mut csv_reader = csv::Reader::from_reader(data.as_bytes());
    let record_it = csv_reader.deserialize::<PosRecord>();
    record_it
        .map(|x| {
            x.expect(
                "JPL deserializer: Error during the csv iterator mapping into a Vec<PosRecord>",
            )
        })
        .collect()
}

#[cfg(test)]
mod jpl_request_test {
    use crate::{error_models::ErrorModel, time::{date_to_mjd, mjd_to_jd}};

    use super::*;

    #[test]
    fn test_jd_list() {
        let jd_list = vec![0.0, 1.5, 2.6];
        assert_eq!(jd_tlist(&jd_list), "0,1.5,2.6")
    }

    #[test]
    fn test_date_to_mjd() {
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let mjd_list = date_to_mjd(&date_list);
        assert_eq!(mjd_list, vec![59399.53291666666, 60672.07462962963]);
        let jd_list = mjd_to_jd(&mjd_list);
        assert_eq!(jd_list, vec![2459400.0329166665, 2460672.5746296295])
    }

    #[test]
    #[ignore]
    fn test_jplvector_request() {
        let state = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let mjd_list = date_to_mjd(&date_list);
        let response_data = request_vector(&mjd_list, &state);
        assert!(response_data.contains(
            "$$SOE
2459400.032916666, A.D. 2021-Jul-04 12:47:24.0000,  2.195672929244244E-01, -9.108330730147444E-01, -3.948423288985838E-01,
2460672.574629629, A.D. 2024-Dec-28 01:47:28.0000, -1.107728032684787E-01,  8.965650072539966E-01,  3.886517577153460E-01,
$$EOE"
        ));
    }

    #[test]
    fn test_deserialize_vector() {
        let fake_jpl = "
Reference frame : ICRF
*******************************************************************************
            JDTDB,            Calendar Date (TDB),                      X,                      Y,                      Z,
**************************************************************************************************************************
$$SOE
2459400.032916666, A.D. 2021-Jul-04 12:47:24.0000,  2.195672929244244E-01, -9.108330730147444E-01, -3.948423288985838E-01,
2460672.574629629, A.D. 2024-Dec-28 01:47:28.0000, -1.107728032684787E-01,  8.965650072539966E-01,  3.886517577153460E-01,
$$EOE
**************************************************************************************************************************
";
        let vec_earth = deserialize_vector(&fake_jpl.into());
        assert_eq!(
            vec_earth,
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
}
