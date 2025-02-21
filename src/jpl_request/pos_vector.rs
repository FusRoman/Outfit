use chrono::NaiveDateTime;
use itertools::Itertools;
use julian_day_converter::JulianDay;
use nalgebra::Vector3;
use regex::Regex;
use reqwest::Client;
use std::str::FromStr;

fn jd_tlist(jd_list: &Vec<f64>) -> String {
    jd_list.iter().join(",")
}

async fn request_vector(jd_list: &Vec<f64>) -> String {
    let requested_params = format!(
        "
!$$SOF
COMMAND='399'
OBJ_DATA='NO'
MAKE_EPHEM='YES'
TABLE_TYPE='VECTORS'
CENTER='500@10'
TLIST_TYPE=JD
TLIST={}
CSV_FORMAT=YES
REF_SYSTEM=ICRF
OUT_UNITS=AU-D
REF_PLANE=FRAME
VEC_TABLE=1
",
        jd_tlist(jd_list)
    );
    let client = Client::new();
    client
        .post("https://ssd.jpl.nasa.gov/api/horizons_file.api")
        .form(&[("format", "text"), ("input", &requested_params)])
        .send()
        .await
        .expect("")
        .text()
        .await
        .expect("")
}

/// Contains the informations from the JPL Horizons vector state query
/// x,y,z are the components of the position vector at the time contained in the jd and date field
#[derive(Debug, serde::Deserialize, PartialEq)]
pub struct PosRecord {
    #[serde(rename = "JDTDB")]
    jd: f64,
    #[serde(rename = "CalendarDate(TDB)")]
    date: String,
    #[serde(rename = "X")]
    x: f64, // km
    #[serde(rename = "Y")]
    y: f64, // km
    #[serde(rename = "Z")]
    z: f64, // km
}

impl PosRecord {
    pub fn pos_vector(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }
}

fn deserialize_vector(jpl_response: &String) -> Vec<PosRecord> {
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

/// Request the JPL Horizon API to get the position vector of Earth
/// with respect to the Sun at different time.
///
/// Argument: a vector of date in Julian Date format
/// Return: a vector of PosRecord
pub async fn get_earth_position(jd_list: &Vec<f64>) -> Vec<PosRecord> {
    let response_data = request_vector(jd_list).await;
    deserialize_vector(&response_data)
}

pub fn date_to_jd(date: &Vec<&str>) -> Vec<f64> {
    date.iter()
        .map(|x| {
            NaiveDateTime::from_str(&x)
                .expect("Error while parsing date to jd")
                .to_jd()
        })
        .collect::<Vec<f64>>()
}

/// Transformation from modified julian date (MJD) in julian date (JD)
///
/// Argument: a vector of MJD
/// Return: a vector of jd
pub fn mjd_to_jd(mjd: &Vec<f64>) -> Vec<f64> {
    mjd.iter().map(|mjd| mjd + 2_400_000.5).collect()
}

#[cfg(test)]
mod jpl_tests {
    use super::*;

    #[test]
    fn test_jd_list() {
        let jd_list = vec![0.0, 1.5, 2.6];
        assert_eq!(jd_tlist(&jd_list), "0,1.5,2.6")
    }

    #[tokio::test]
    async fn test_jplvector_request() {
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let jd_list = date_to_jd(&date_list);
        let response_data = request_vector(&jd_list).await;
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

    #[tokio::test]
    async fn test_get_earth_pos() {
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let jd_list = date_to_jd(&date_list);
        println!("{jd_list:?}");
        let earth_vector = get_earth_position(&jd_list).await;
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

    #[tokio::test]
    async fn test_earth_pos_with_mjd() {
        let test_date = [57028.479297592596, 57049.245147592592, 57063.977117592593];
        let test_jd: Vec<f64> = mjd_to_jd(&test_date.into());
        let earth_vector = get_earth_position(&test_jd).await;
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
                    x: -0.5891845254244624,
                    y: 0.7238535049311241,
                    z: 0.3138036934121133
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
