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
/// x,y,z are the components of the sun position vector at the time contained in the jd and date field
/// vx,vy,vz are the components of the sun velocity vector
/// light_time is the one-way newtonian light-time to the target
/// range is the distance to the target
/// range_rate is the radial velocity with respect to the coordinate center
#[derive(Debug, serde::Deserialize, PartialEq)]
pub struct HelioRecord {
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
    #[serde(rename = "VX")]
    vx: f64, // km/sec
    #[serde(rename = "VY")]
    vy: f64, // km/sec
    #[serde(rename = "VZ")]
    vz: f64, //km/sec
    #[serde(rename = "LT")]
    light_time: f64, //sec
    #[serde(rename = "RG")]
    range: f64, //km
    #[serde(rename = "RR")]
    range_rate: f64, //km/sec
}

impl HelioRecord {
    pub fn pos_vector(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    pub fn vel_vector(&self) -> Vector3<f64> {
        Vector3::new(self.vx, self.vy, self.vz)
    }
}

fn deserialize_vector(jpl_response: &String) -> Vec<HelioRecord> {
    // regex to match the data part of the jpl horizon response
    let data_regex = Regex::new(r"\$\$SOE\n([^]]*),\n\$\$EOE").unwrap();
    // regex to match the header part of the jpl horizon response
    let header_regex = Regex::new(r"J2000\.0\n\*{79}([^]]*)\*{266}\n\$\$SOE").unwrap();

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

    // Use serde to deserialize the string into a vector of HelioRecord.
    let mut csv_reader = csv::Reader::from_reader(data.as_bytes());
    let record_it = csv_reader.deserialize::<HelioRecord>();
    record_it
        .map(|x| {
            x.expect(
                "JPL deserializer: Error during the csv iterator mapping into a Vec<HelioRecord>",
            )
        })
        .collect()
}

/// Request the JPL Horizon API to get vector state of the sun
/// at different time.
///
/// Argument: a vector of date in Julian Date format
/// Return: a vector of HelioRecord
pub async fn get_helio_pos(jd_list: &Vec<f64>) -> Vec<HelioRecord> {
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
2459400.032916666, A.D. 2021-Jul-04 12:47:24.0000,  3.284679949685707E+07, -1.485106329685027E+08,  7.156044432416558E+03,  2.860834584138197E+01,  6.313673173977221E+00, -1.324421108332530E-03,  5.073500062298830E+02,  1.520997054339719E+08,  1.344506803319922E-02,
2460672.574629629, A.D. 2024-Dec-28 01:47:28.0000, -1.657137550043441E+07,  1.461839130215131E+08, -7.788453816257417E+03, -3.009513750045504E+01, -3.462277804314375E+00,  8.254731385908265E-04,  4.907400928461101E+02,  1.471201786734835E+08, -5.037717899031938E-02,
$$EOE"
        ));
    }

    #[test]
    fn test_deserialize_vector() {
        let fake_jpl = "
Reference frame : Ecliptic of J2000.0
*******************************************************************************
            JDTDB,            Calendar Date (TDB),                      X,                      Y,                      Z,                     VX,                     VY,                     VZ,                     LT,                     RG,                     RR,
**************************************************************************************************************************************************************************************************************************************************************************
$$SOE
2459400.032916666, A.D. 2021-Jul-04 12:47:24.0000, -3.285213437722069E+07,  1.485094907917824E+08, -1.045064456634223E+04, -2.862749646371449E+01, -6.670100246343959E+00,  1.559012932249821E-01,  5.073501301870917E+02,  1.520997425954083E+08, -3.293921397855185E-01,
2460672.574629629, A.D. 2024-Dec-28 01:47:28.0000,  1.656606640874423E+07, -1.461858433966990E+08,  4.837673705182970E+03,  3.013868677111592E+01,  3.107607789780265E+00,  1.528420680554570E-01,  4.907444960886756E+02,  1.471214987323954E+08,  3.058149938328106E-01,
$$EOE
**************************************************************************************************************************************************************************************************************************************************************************
";
        let vec_helio = deserialize_vector(&fake_jpl.into());
        assert_eq!(
            vec_helio,
            vec![
                HelioRecord {
                    jd: 2459400.032916666,
                    date: "A.D.2021-Jul-0412:47:24.0000".into(),
                    x: -32852134.37722069,
                    y: 148509490.7917824,
                    z: -10450.64456634223,
                    vx: -28.62749646371449,
                    vy: -6.670100246343959,
                    vz: 0.1559012932249821,
                    light_time: 507.3501301870917,
                    range: 152099742.5954083,
                    range_rate: -0.3293921397855185
                },
                HelioRecord {
                    jd: 2460672.574629629,
                    date: "A.D.2024-Dec-2801:47:28.0000".into(),
                    x: 16566066.40874423,
                    y: -146185843.396699,
                    z: 4837.67370518297,
                    vx: 30.13868677111592,
                    vy: 3.107607789780265,
                    vz: 0.152842068055457,
                    light_time: 490.7444960886756,
                    range: 147121498.7323954,
                    range_rate: 0.3058149938328106
                }
            ]
        );
    }

    #[tokio::test]
    async fn test_get_helio_pos() {
        let date_list = vec!["2021-07-04T12:47:24", "2024-12-28T01:47:28"];
        let jd_list = date_to_jd(&date_list);
        let helio_vector = get_helio_pos(&jd_list).await;
        assert_eq!(
            helio_vector,
            vec![
                HelioRecord {
                    jd: 2459400.032916666,
                    date: "A.D.2021-Jul-0412:47:24.0000".into(),
                    x: 32846799.49685707,
                    y: -148510632.9685027,
                    z: 7156.044432416558,
                    vx: 28.60834584138197,
                    vy: 6.313673173977221,
                    vz: -0.00132442110833253,
                    light_time: 507.350006229883,
                    range: 152099705.4339719,
                    range_rate: 0.01344506803319922
                },
                HelioRecord {
                    jd: 2460672.574629629,
                    date: "A.D.2024-Dec-2801:47:28.0000".into(),
                    x: -16571375.50043441,
                    y: 146183913.0215131,
                    z: -7788.453816257417,
                    vx: -30.09513750045504,
                    vy: -3.462277804314375,
                    vz: 0.0008254731385908265,
                    light_time: 490.7400928461101,
                    range: 147120178.6734835,
                    range_rate: -0.05037717899031938
                }
            ]
        );
    }
}
