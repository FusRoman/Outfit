use nom::{
    bytes::complete::take,
    multi::count,
    number::complete::{le_f64, le_i32, le_u32},
    IResult, Parser,
};

use crate::{constants::MJD, jpl_ephem::download_jpl_file::EphemFilePath};

use super::{
    horizon_ids::HorizonID, horizon_records::HorizonRecord, horizon_version::JPLHorizonVersion,
    interpolation_result::InterpResult,
};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read, Seek},
};

/// Type for the JPL ephemeris file
/// The file is a binary file containing the ephemeris data
/// The file is divided into blocks, each block contains the data for a specific time interval
/// Each block is in the following vec type and can be accessed by the index
/// The HashMap contains the body number (0-14) as the key
/// and a vector of HorizonRecord as the value
type HorizonRecords = Vec<HashMap<u8, Vec<HorizonRecord>>>;

/// Header containing informations for each solar system body
/// The first index is the body number (0-14)
/// The second index is a three element array
/// The first element is the offset in the block corresponding to the body
/// The second element is the number of Tchebyshev coefficients
/// The third element is the number of subintervals
pub type IPT = [[u32; 3]; 15];

#[derive(Debug, PartialEq, Clone)]
pub struct HorizonHeader {
    jpl_version: String,
    ipt: IPT,
    start_period: f64,
    end_period: f64,
    period_lenght: f64,
    recsize: usize,
    earth_moon_mass_ratio: f64,
}

/// Structure to hold the data from the JPL ephemeris file
/// The header contains the JPL version and the IPT table
/// The records are stored in a vector of hashmaps
/// where the key is the body number (0-14) and the value is a vector of HorizonRecord
/// Each HorizonRecord contains the start and end JD, and the coefficients for the Tchebyshev polynomial
#[derive(Debug, Clone)]
pub struct HorizonData {
    header: HorizonHeader,
    records: HorizonRecords,
}

/// Dimension for each record body
///
/// Arguments
/// ---------
/// * `index` : Index of the body
///
/// Returns
/// -------
/// * Dimension of the body
fn dimension(index: usize) -> usize {
    match index {
        0..=10 => 3, // planètes + Soleil
        11 => 2,     // Earth-Moon barycenter
        12 => 3,     // Lunar librations
        13 => 3,     // Lunar Euler angle rates
        14 => 1,     // TT-TDB offsets
        _ => 0,      // par sécurité
    }
}

/// Compute the size of the record in bytes
///
/// Arguments
/// ---------
/// * `ipt` : IPT table
///
/// Returns
/// -------
/// * Size of the record in bytes
fn compute_recsize(ipt: IPT) -> usize {
    let kernel_size: usize = 4
        + (0..15)
            .map(|i| {
                let n_subintervals = ipt[i][1] as usize;
                let n_coeffs = ipt[i][2] as usize;
                let dim = dimension(i);
                2 * n_subintervals * n_coeffs * dim
            })
            .sum::<usize>();

    kernel_size * 4
}

/// Parse the IPT table from the binary file
///
/// Arguments
/// ---------
/// * `input` : Input buffer containing the IPT table in bytes to decode
///
/// Returns
/// -------
/// * A tuple containing the remaining input and the IPT table
fn parse_ipt(input: &[u8]) -> IResult<&[u8], (IPT, u32)> {
    let mut ipt: IPT = [[0; 3]; 15];

    let mut input_remaining = input;
    for i in 0..36 {
        let (rest, val) = le_u32(input_remaining)?;

        input_remaining = rest;
        let row = i / 3;
        let col = i % 3;
        ipt[row][col] = val;
    }

    let (remaining_input, numde) = le_u32(input_remaining)?;

    let (input_remaining, lpt) = count(le_u32, 3).parse(remaining_input)?;

    ipt[12] = lpt
        .try_into()
        .expect("Expected lpt to have exactly 3 elements");

    Ok((input_remaining, (ipt, numde)))
}

/// Parse the remaining part of the IPT table
///
/// Arguments
/// ---------
/// * `input` : Input buffer containing the remaining IPT table in bytes to decode
///
/// Returns
/// -------
/// * A tuple containing the remaining input and the IPT 13 and 14 elements
fn parse_ipt_13_14(input: &[u8]) -> IResult<&[u8], [[u32; 3]; 2]> {
    fn parse_three(mut input: &[u8]) -> IResult<&[u8], [u32; 3]> {
        let mut result = [0u32; 3];
        for slot in &mut result {
            let (rest, val) = le_u32(input)?;
            *slot = val;
            input = rest;
        }
        Ok((input, result))
    }

    let (input, ipt_13) = parse_three(input)?;
    let (input, ipt_14) = parse_three(input)?;

    Ok((input, [ipt_13, ipt_14]))
}

/// Read the IPT\[13\] and IPT\[14\] elements from the file
///
/// Arguments
/// ---------
/// * `file` : File to read from
/// * `ncon` : Number of constants in the jpl file
/// * `de_version` : JPL Horizon version
///
/// Returns
/// -------
/// * IPT\[13\] and IPT\[14\] elements
fn read_ipt_13_14(
    file: &mut BufReader<File>,
    ncon: u32,
    de_version: &JPLHorizonVersion,
) -> std::io::Result<Option<[[u32; 3]; 2]>> {
    if *de_version < JPLHorizonVersion::DE440 || ncon <= 400 {
        return Ok(None); // IPT[13] and IPT[14] are not used in version earlier than DE440 or for ncon <= 400
    }

    const START_400TH_CONSTANT_NAME: u64 = 2856;
    let offset = START_400TH_CONSTANT_NAME + (ncon as u64 - 400) * 6;

    file.seek(std::io::SeekFrom::Start(offset))?;
    let mut buffer = [0u8; 24];
    file.read_exact(&mut buffer)?;

    let (_, ipt_extra) = parse_ipt_13_14(&buffer).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Failed to parse ipt[13][14]",
        )
    })?;

    Ok(Some(ipt_extra))
}

/// Parse a block of data from the binary file
///
/// Arguments
/// ---------
/// * `input` : Input buffer containing the block in bytes to decode
///
/// Returns
/// -------
/// * A tuple containing the remaining input and the block data
fn parse_block(input: &[u8], ncoeff: usize) -> IResult<&[u8], Vec<f64>> {
    count(le_f64, ncoeff).parse(input)
}

/// Extract the body records from a block of data
///
/// Arguments
/// ---------
/// * `block` : Block of data containing the coefficients
/// * `ipt` : IPT table for the body
///
/// Returns
/// -------
/// * A vector of HorizonRecord containing the coefficients for the body
fn extract_body_records(block: &[f64], ipt: [u32; 3]) -> Vec<HorizonRecord> {
    let offset = ipt[0] as usize;
    let n_coeffs = ipt[1] as usize;
    let n_subs = ipt[2] as usize;

    let jd_start = block[0];
    let jd_end = block[1];

    let coeffs = &block[2..]; // skip JD start/end
    let mut records = Vec::with_capacity(n_subs);

    for subinterval_number in 0..n_subs {
        let record = HorizonRecord::new(
            jd_start,
            jd_end,
            coeffs,
            offset,
            subinterval_number,
            n_coeffs,
        );

        records.push(record);
    }

    records
}

/// Parse all blocks from the binary file
///
/// Arguments
/// ---------
/// * `file` : File to read from
/// * `recsize` : Size of each record in bytes
/// * `ncoeff` : Number of coefficients in each record
/// * `ipt` : IPT table
///
/// Returns
/// -------
/// * A vector of hashmaps containing the records for each body
fn parse_all_blocks(
    file: &mut BufReader<File>,
    recsize: usize,
    ncoeff: usize,
    ipt: IPT,
) -> std::io::Result<HorizonRecords> {
    let offset = 2 * recsize;
    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    let n_blocks = buffer.len() / recsize;

    (0..n_blocks)
        .map(|i| {
            let start = i * recsize;
            let end = start + recsize;
            let slice = &buffer[start..end];

            let (_, block) = parse_block(slice, ncoeff).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse block {i}"),
                )
            })?;

            let records = (0..12)
                .map(|body| {
                    let ipt_body = ipt[body];
                    let body_rec = extract_body_records(&block, ipt_body);
                    (body as u8, body_rec)
                })
                .collect::<HashMap<u8, _>>();

            Ok(records)
        })
        .collect()
}

impl HorizonData {
    /// Read the JPL ephemeris file and parse the data
    ///
    /// Arguments
    /// ---------
    /// * `ephem_path` : Path to the JPL ephemeris file
    ///
    /// Returns
    /// -------
    /// * A HorizonData struct containing the header and records
    pub fn read_horizon_file(ephem_path: &EphemFilePath) -> Self {
        let version = match ephem_path {
            EphemFilePath::JPLHorizon(_, version) => version,
            _ => panic!("Invalid ephemeris file path"),
        };

        const OLD_MAX: usize = 400;

        /// Parse une chaîne Fortran CHAR*6 (6 bytes)
        fn parse_char6(input: &[u8]) -> IResult<&[u8], String> {
            let (rest, raw) = take(6usize)(input)?;
            Ok((rest, String::from_utf8_lossy(raw).trim_end().to_string()))
        }

        /// Parse les en-têtes TTL (14×3 CHAR*6)
        fn parse_ttl(input: &[u8]) -> IResult<&[u8], Vec<String>> {
            count(parse_char6, 14 * 3).parse(input)
        }

        /// Parse le champ CNAM (400 chaînes CHAR*6)
        fn parse_cnam(input: &[u8]) -> IResult<&[u8], Vec<String>> {
            count(parse_char6, OLD_MAX).parse(input)
        }

        let mut file = BufReader::new(
            File::open(ephem_path.path())
                .unwrap_or_else(|_| panic!("Failed to open the JPL ephemeris file: {ephem_path}")),
        );

        let mut file_data = vec![0u8; 1 << 12];
        file.read_exact(&mut file_data)
            .expect("Failed to read the JPL ephemeris file.");

        let (input, _) = parse_ttl(&file_data).expect("failed to parse TTL");

        let (input, _) = parse_cnam(input).expect("failed to parse CNAM");

        let (input, ss) = count(le_f64::<_, nom::error::Error<_>>, 3)
            .parse(input)
            .expect("failed to parse SS");

        let (input, ncon) = le_i32::<_, nom::error::Error<_>>(input).expect("failed to parse NCON");

        let (input, _) =
            take::<_, _, nom::error::Error<_>>(8usize)(input).expect("failed to parse 4 bytes");
        let (_, earth_moon_mass_ratio) = le_f64::<_, nom::error::Error<_>>(input)
            .expect("failed to parse Earth-Moon mass ratio");

        file.seek(std::io::SeekFrom::Start(2696))
            .expect("Failed to seek to the JPL ephemeris file IPT.");
        let mut buffer = [0u8; 204]; // JPL_HEADER_SIZE
        file.read_exact(&mut buffer)
            .expect("Failed to read the JPL ephemeris file IPT.");

        let (_, (mut ipt, numde)) = parse_ipt(&buffer).expect("failed to parse IPT header");

        let ipt_extra = read_ipt_13_14(&mut file, ncon as u32, version)
            .expect("Failed to read IPT[13] and IPT[14]");

        ipt[13] = ipt_extra
            .as_ref()
            .map(|ipt_extra| ipt_extra[0])
            .unwrap_or([0; 3]);

        ipt[14] = ipt_extra
            .as_ref()
            .map(|ipt_extra| ipt_extra[1])
            .unwrap_or([0; 3]);

        let recsize = compute_recsize(ipt);

        let ncoeff = recsize / 8;
        let blocks = parse_all_blocks(&mut file, recsize, ncoeff, ipt).unwrap();

        HorizonData {
            header: HorizonHeader {
                jpl_version: format!("DE{numde}"),
                ipt,
                start_period: ss[0],
                end_period: ss[1],
                period_lenght: ss[2],
                recsize,
                earth_moon_mass_ratio,
            },
            records: blocks,
        }
    }

    /// Get the index of the record corresponding to the given ephemeris time
    /// The index nr is the number of the block in the file corresponding the the requested time.
    /// The tau value is the time fraction in the block.
    ///
    /// Arguments
    /// ---------
    /// * `et` : Ephemeris time (MJD)
    ///
    /// Returns
    /// -------
    /// * A tuple containing the index of the record and the tau value
    fn get_record_index(&self, et: MJD) -> (usize, f64) {
        // ephem_start and ephem_end are in JD
        let (ephem_start, ephem_end, ephem_step) = (
            self.header.start_period,
            self.header.end_period,
            self.header.period_lenght,
        );

        let jd_base = 2_400_000.5;
        let et_jd = jd_base + et.trunc();

        if et_jd < ephem_start || et_jd > ephem_end {
            panic!("Time outside ephemeris range");
        }

        let mut nr = ((et_jd - ephem_start) / ephem_step).floor() as usize;

        if (et_jd - ephem_end).abs() < 1e-10 {
            nr -= 1;
        }

        let interval_start = (nr as f64) * ephem_step + ephem_start;
        let tau = ((et_jd - interval_start) + et.fract()) / ephem_step;
        (nr, tau)
    }

    /// Get the record corresponding to the given ephemeris time and body
    ///
    /// Arguments
    /// ---------
    /// * `body` : Body number (0-14)
    /// * `et` : Ephemeris time (MJD)
    ///
    /// Returns
    /// -------
    /// * A tuple containing the record and the tau value
    fn get_record_horizon(&self, body: u8, et: MJD) -> Option<(&HorizonRecord, f64)> {
        let (nr, tau) = self.get_record_index(et);

        let records = &self.records[nr];
        let record_body = records
            .get(&body)
            .unwrap_or_else(|| panic!("Failed to get record for body {body} in block {nr}"));

        let ipt_body = self.header.ipt[body as usize];
        let n_subs = ipt_body[2] as usize;

        let sub_index = (tau * n_subs as f64).floor().min(n_subs as f64 - 1.0) as usize;

        Some((&record_body[sub_index], tau))
    }

    /// Get the ephemeris data for a given target and center body at a given ephemeris time
    /// The target and center bodies are specified by their IDs (0-14).
    /// The ephemeris data is computed by interpolating the records for the target and center bodies
    /// at the given ephemeris time.
    /// The result is the difference between the target and center ephemeris data.
    ///
    /// Arguments
    /// ---------
    /// * `target` : Target body number (0-14)
    /// * `center` : Center body number (0-14)
    /// * `et` : Ephemeris time (MJD)
    /// * `compute_velocity` : Flag to compute velocity
    /// * `compute_acceleration` : Flag to compute acceleration
    ///
    /// Returns
    /// -------
    /// * An InterpResult struct containing the position, velocity, and acceleration
    ///   of the target body relative to the center body at the given ephemeris time.
    ///   The position is in Km, velocity in Km/day, and acceleration in Km/day^2.
    /// * The velocity and acceleration are optional, depending on the user request.
    pub fn ephemeris(
        &self,
        target: HorizonID,
        center: HorizonID,
        et: MJD,
        compute_velocity: bool,
        compute_acceleration: bool,
    ) -> InterpResult {
        let ipt_target = self.header.ipt[target as usize];
        let ipt_center = self.header.ipt[center as usize];
        let (record, tau) = self.get_record_horizon(target.into(), et).unwrap();
        let mut interp_target = record.interpolate(
            tau,
            compute_velocity,
            compute_acceleration,
            ipt_target[2] as usize,
        );

        if target == HorizonID::Earth {
            let ipt_moon = self.header.ipt[HorizonID::Moon as usize];
            let (record, tau) = self.get_record_horizon(HorizonID::Moon.into(), et).unwrap();
            let interp_moon = record.interpolate(
                tau,
                compute_velocity,
                compute_acceleration,
                ipt_moon[2] as usize,
            );
            interp_target = interp_target - interp_moon / (1. + self.header.earth_moon_mass_ratio);
        }

        let (record, tau) = self.get_record_horizon(center.into(), et).unwrap();
        let interp_center = record.interpolate(
            tau,
            compute_velocity,
            compute_acceleration,
            ipt_center[2] as usize,
        );

        interp_target - interp_center
    }
}

#[cfg(test)]
mod test_horizon_reader {
    #[cfg(feature = "jpl-download")]
    use super::*;

    #[cfg(feature = "jpl-download")]
    use crate::unit_test_global::JPL_EPHEM_HORIZON;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_reader_from_horizon() {
        assert_eq!(
            JPL_EPHEM_HORIZON.header,
            HorizonHeader {
                jpl_version: "DE440".to_string(),
                ipt: [
                    [3, 14, 4,],
                    [171, 10, 2,],
                    [231, 13, 2,],
                    [309, 11, 1,],
                    [342, 8, 1,],
                    [366, 7, 1,],
                    [387, 6, 1,],
                    [405, 6, 1,],
                    [423, 6, 1,],
                    [441, 13, 8,],
                    [753, 11, 2,],
                    [819, 10, 4,],
                    [899, 10, 4,],
                    [1019, 0, 0,],
                    [1019, 0, 0,],
                ],
                start_period: 2287184.5,
                end_period: 2688976.5,
                period_lenght: 32.0,
                recsize: 8144,
                earth_moon_mass_ratio: 81.30056822149722,
            }
        );

        assert_eq!(JPL_EPHEM_HORIZON.records.len(), 12556);
        assert_eq!(
            JPL_EPHEM_HORIZON
                .records
                .iter()
                .fold(0, |acc, x| acc + x.len()),
            150672
        );

        assert_eq!(
            &JPL_EPHEM_HORIZON.records[0].get(&0).unwrap()[0],
            &HorizonRecord {
                start_jd: 2287184.5,
                end_jd: 2287216.5,
                x: vec![
                    -45337704.29199142,
                    -11420952.2182182,
                    1231640.71525489,
                    13474.74253284046,
                    -4516.491868368539,
                    255.9277522073883,
                    -0.6852443171165454,
                    -1.186291661304585,
                    0.1035908763014567,
                    -0.002445022374643364,
                    -0.0003798510942101873,
                    4.866105543177588e-5,
                    -2.092814231764934e-6,
                    -1.118766077984058e-7
                ],
                y: vec![
                    19369902.32537584,
                    -12637978.7311934,
                    -560491.8087798061,
                    75269.29929452346,
                    -1778.157215889909,
                    -139.1111450981129,
                    18.58942667472901,
                    -0.8470743363574242,
                    -0.02513657259988719,
                    0.006832859866311067,
                    -0.0004590542406706136,
                    5.317330201500909e-7,
                    2.764140697666422e-6,
                    -2.579709214481493e-7
                ],
                z: vec![
                    15085546.93080799,
                    -5541484.32081567,
                    -428281.1733889182,
                    38734.17578050081,
                    -474.2794185515511,
                    -101.0728095721337,
                    9.98761508428064,
                    -0.3272813126206938,
                    -0.02428417941794788,
                    0.003901383011452982,
                    -0.000204980181894197,
                    -4.825494149048871e-6,
                    1.694428982754255e-6,
                    -1.260377702474755e-7
                ]
            }
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_get_record_from_horizon() {
        let (record, tau) = JPL_EPHEM_HORIZON
            .get_record_horizon(4, 57028.479297592596)
            .unwrap();

        assert_eq!(tau, 0.6399780497686152);

        assert_eq!(
            record,
            &HorizonRecord {
                start_jd: 2457008.5,
                end_jd: 2457040.5,
                x: vec![
                    -558258575.175456,
                    -13081137.99303625,
                    70236.79190208673,
                    240.5982017483768,
                    -0.8633795271054904,
                    -0.0002571950240619499,
                    4.544021170810181e-6,
                    -1.839536840631949e-8
                ],
                y: vec![
                    515866454.9405554,
                    -10984765.90508825,
                    -64874.14676298688,
                    261.5078409285777,
                    0.4491808227783525,
                    -0.001923759413667307,
                    2.250631406524948e-6,
                    -6.23548505589978e-8
                ],
                z: vec![
                    234693317.6472925,
                    -4389890.166800962,
                    -29516.72247187927,
                    106.2324198779676,
                    0.2135457077378355,
                    -0.000818018169817293,
                    9.216197575953068e-7,
                    -2.484668250013335e-8
                ]
            },
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_get_record_index() {
        let (index, tau) = JPL_EPHEM_HORIZON.get_record_index(57028.479297592596);

        assert_eq!(index, 5307);
        assert_eq!(tau, 0.6399780497686152);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_interpolation_from_horizon() {
        let (record, tau) = JPL_EPHEM_HORIZON
            .get_record_horizon(10, 57028.479297592596)
            .unwrap();

        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [[428149.04652929713, -105270.2354841389, -68083.3417807072]].into(),
                velocity: Some([[589.5451313546943, 729.3492107653134, 300.3651374864013]].into()),
                acceleration: Some(
                    [[-0.919215789187584, 0.8829808730528121, 0.38984144066801635]].into()
                )
            }
        );

        let (record, tau) = JPL_EPHEM_HORIZON
            .get_record_horizon(10, 57_049.231_857_592_59)
            .unwrap();

        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [[440183.15997455275, -89933.41046658876, -61760.61450611751]].into(),
                velocity: Some([[569.7900066838628, 749.1773000798205, 309.2398409126585]].into()),
                acceleration: Some(
                    [[-1.038549415842939, 0.9990469757389817, 0.4553928281168678]].into()
                )
            }
        );

        let (record, tau) = JPL_EPHEM_HORIZON
            .get_record_horizon(10, 60781.51949044435)
            .unwrap();
        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [[-742814.3000341875, -727671.3536844663, -288321.53733077314]].into(),
                velocity: Some(
                    [[1085.6256324761375, -327.2648113240611, -162.13166222548898]].into()
                ),
                acceleration: Some(
                    [[-0.06525161089395584, 1.261197328255844, 0.5376907387973575]].into()
                )
            }
        )
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_target_center_interpolation() {
        let interp = JPL_EPHEM_HORIZON.ephemeris(
            HorizonID::Earth,
            HorizonID::Sun,
            60781.51949044435,
            true,
            true,
        );

        assert_eq!(
            interp.to_au(),
            InterpResult {
                position: [[
                    -0.8988034555610618,
                    -0.4096429669081428,
                    -0.17756458835186803
                ]]
                .into(),
                velocity: Some(
                    [[
                        0.007368756100393145,
                        -0.014193450423040196,
                        -0.006152419296313352
                    ]]
                    .into()
                ),
                acceleration: Some(
                    [[
                        0.0002624924486094171,
                        0.00011873474463671509,
                        5.13298993717986e-5
                    ]]
                    .into()
                )
            }
        );

        let interp = JPL_EPHEM_HORIZON.ephemeris(
            HorizonID::Mars,
            HorizonID::Sun,
            60781.51949044435,
            true,
            true,
        );

        assert_eq!(
            interp.to_au(),
            InterpResult {
                position: [[-1.5211490583088887, 0.601249019351933, 0.31680918257233887]].into(),
                velocity: Some(
                    [[
                        -0.005170287251815057,
                        -0.010584792921766313,
                        -0.0047155422859661905
                    ]]
                    .into()
                ),
                acceleration: Some(
                    [[
                        9.733944178212551e-5,
                        -3.84704427625744e-5,
                        -2.0271140253173187e-5
                    ]]
                    .into()
                )
            }
        );

        let interp = JPL_EPHEM_HORIZON.ephemeris(
            HorizonID::Earth,
            HorizonID::Sun,
            52550.18467592593,
            true,
            true,
        );

        assert_eq!(
            interp.to_au(),
            InterpResult {
                position: [[0.9861809593158523, 0.1557130752386633, 0.06750574448131498]].into(),
                velocity: Some(
                    [[
                        -0.003193635816506268,
                        0.015502851910183767,
                        0.006721021432383671
                    ]]
                    .into()
                ),
                acceleration: Some(
                    [[
                        -0.0002926932324959577,
                        -4.506599368500431e-5,
                        -1.9370434372826522e-5
                    ]]
                    .into()
                )
            }
        );
    }
}
