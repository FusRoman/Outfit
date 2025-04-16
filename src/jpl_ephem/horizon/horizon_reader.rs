use camino::Utf8Path;
use nom::{
    bytes::complete::take,
    multi::count,
    number::complete::{le_f64, le_i32, le_u32},
    IResult, Parser,
};

use crate::jpl_ephem::download_jpl_file::{EphemFilePath, JPLHorizonVersion};

use super::horizon_records::HorizonRecord;
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
type HorizonRecords = Vec<HashMap<u32, Vec<HorizonRecord>>>;

/// Header containing informations for each solar system body
/// The first index is the body number (0-14)
/// The second index is a three element array
/// The first element is the offset in the block corresponding to the body
/// The second element is the number of Tchebyshev coefficients
/// The third element is the number of subintervals
pub type IPT = [[u32; 3]; 15];

#[derive(Debug, PartialEq)]
pub struct HorizonHeader {
    jpl_version: String,
    ipt: IPT,
    start_period: f64,
    end_period: f64,
    period_lenght: f64,
    recsize: usize,
}

/// Structure to hold the data from the JPL ephemeris file
/// The header contains the JPL version and the IPT table
/// The records are stored in a vector of hashmaps
/// where the key is the body number (0-14) and the value is a vector of HorizonRecord
/// Each HorizonRecord contains the start and end JD, and the coefficients for the Tchebyshev polynomial
#[derive(Debug)]
pub struct HorizonData {
    header: HorizonHeader,
    records: HorizonRecords,
}

/// Dimensions pour chaque index IPT (indice du corps)
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

/// Calcule la taille en octets d'un enregistrement (recsize)
pub fn compute_recsize(ipt: IPT) -> usize {
    let mut kernel_size = 4; // mots de 32 bits

    for i in 0..15 {
        let n_subintervals = ipt[i][1] as usize;
        let n_coeffs = ipt[i][2] as usize;
        let dim = dimension(i);

        kernel_size += 2 * n_subintervals * n_coeffs * dim;
    }

    // 1 mot = 4 octets → recsize en octets
    kernel_size * 4
}

/// Parse les 204 octets du header (à partir de offset 2652) pour récupérer ipt[0..13]
pub fn parse_ipt(input: &[u8]) -> IResult<&[u8], (IPT, u32)> {
    let mut ipt: IPT = [[0; 3]; 15];

    // On saute les 5 premiers doubles (5 * 8 = 40 octets)
    // let (input, _) = take(40usize)(input)?;

    // On lit ensuite 40 entiers u32 (40 * 4 = 160 octets)
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

pub fn parse_ipt_13_14(input: &[u8]) -> IResult<&[u8], [[u32; 3]; 2]> {
    let mut ipt_13 = [0u32; 3];
    let mut ipt_14 = [0u32; 3];
    let mut remaining = input;

    for i in 0..3 {
        let (rest, val) = le_u32(remaining)?;
        ipt_13[i] = val;
        remaining = rest;
    }

    for i in 0..3 {
        let (rest, val) = le_u32(remaining)?;
        ipt_14[i] = val;
        remaining = rest;
    }

    Ok((remaining, [ipt_13, ipt_14]))
}

fn read_ipt_13_14(
    file: &mut BufReader<File>,
    ncon: u32,
    de_version: &JPLHorizonVersion,
) -> std::io::Result<Option<[[u32; 3]; 2]>> {
    if *de_version < JPLHorizonVersion::DE440 || ncon <= 400 {
        return Ok(None); // pas présent
    }

    const START_400TH_CONSTANT_NAME: u64 = 2856;
    let offset = START_400TH_CONSTANT_NAME + (ncon as u64 - 400) * 6;

    file.seek(std::io::SeekFrom::Start(offset))?;
    let mut buffer = [0u8; 24]; // 6 * 4
    file.read_exact(&mut buffer)?;

    let (_, ipt_extra) = parse_ipt_13_14(&buffer).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "Failed to parse ipt[13][14]",
        )
    })?;

    Ok(Some(ipt_extra))
}

/// Parse un bloc complet : ncoeff * f64
fn parse_block(input: &[u8], ncoeff: usize) -> IResult<&[u8], Vec<f64>> {
    count(le_f64, ncoeff).parse(input)
}

pub fn extract_body_records(block: &[f64], ipt: [u32; 3]) -> Vec<HorizonRecord> {
    let offset = ipt[0] as usize;
    let n_coeffs = ipt[1] as usize;
    let n_subs = ipt[2] as usize;

    let jd_start = block[0];
    let jd_end = block[1];

    let coeffs = &block[2..]; // skip JD start/end
    let mut records = Vec::with_capacity(n_subs);

    for i in 0..n_subs {
        let mut x = Vec::with_capacity(n_coeffs);
        let mut y = Vec::with_capacity(n_coeffs);
        let mut z = Vec::with_capacity(n_coeffs);

        let base = offset - 3 + i * n_coeffs * 3;

        for j in 0..n_coeffs {
            let idx = base + j;

            let x_val = coeffs.get(idx).copied().expect("x coeff out of bounds");
            let y_val = coeffs
                .get(idx + n_coeffs)
                .copied()
                .expect("y coeff out of bounds");
            let z_val = coeffs
                .get(idx + n_coeffs * 2)
                .copied()
                .expect("z coeff out of bounds");

            x.push(x_val);
            y.push(y_val);
            z.push(z_val);
        }

        records.push(HorizonRecord {
            start_jd: jd_start,
            end_jd: jd_end,
            x,
            y,
            z,
        });
    }

    records
}

pub fn parse_all_blocks(
    file: &mut BufReader<File>,
    recsize: usize,
    ncoeff: usize,
    ipt: IPT,
) -> std::io::Result<HorizonRecords> {
    // Skip les 2 premiers blocs (header texte + header binaire)
    let offset = 2 * recsize;
    file.seek(std::io::SeekFrom::Start(offset as u64))?;

    // Lis tout le reste du fichier
    let mut buffer = Vec::new();
    file.read_to_end(&mut buffer)?;

    // Vérifie que la taille est un multiple du recsize
    let n_blocks = buffer.len() / recsize;
    let mut blocks: HorizonRecords = Vec::with_capacity(n_blocks);

    for i in 0..n_blocks {
        let start = i * recsize;
        let end = start + recsize;
        let slice = &buffer[start..end];

        let (_, block) = parse_block(slice, ncoeff).map_err(|_| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Failed to parse block {i}"),
            )
        })?;

        let mut records = HashMap::new();
        for body in 0..12 {
            let ipt_body = ipt[body];

            let body_rec = extract_body_records(&block, ipt_body);
            records.insert(body as u32, body_rec);
        }
        blocks.push(records);
    }

    Ok(blocks)
}

impl HorizonData {
    fn read_horizon_file(ephem_path: &EphemFilePath) -> Self {
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
                .expect(format!("Failed to open the JPL ephemeris file: {}", ephem_path).as_str()),
        );

        let mut file_data = vec![0u8; 1 << 12];
        file.read_exact(&mut file_data)
            .expect("Failed to read the JPL ephemeris file.");

        let (input, _) = parse_ttl(&file_data).expect("failed to parse TTL");

        let (input, _) = parse_cnam(input).expect("failed to parse CNAM");

        let (input, ss) = count(le_f64::<_, nom::error::Error<_>>, 3)
            .parse(input)
            .expect("failed to parse SS");

        let (_, ncon) = le_i32::<_, nom::error::Error<_>>(input).expect("failed to parse NCON");

        file.seek(std::io::SeekFrom::Start(2696))
            .expect("Failed to seek to the JPL ephemeris file IPT.");
        let mut buffer = [0u8; 204]; // JPL_HEADER_SIZE
        file.read_exact(&mut buffer)
            .expect("Failed to read the JPL ephemeris file IPT.");

        let (_, (mut ipt, numde)) = parse_ipt(&buffer).expect("failed to parse IPT header");

        let ipt_extra = read_ipt_13_14(&mut file, ncon as u32, &version)
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
                jpl_version: format!("DE{}", numde),
                ipt,
                start_period: ss[0],
                end_period: ss[1],
                period_lenght: ss[2],
                recsize: recsize,
            },
            records: blocks,
        }
    }

    fn get_record_index(&self, et: f64) -> (usize, f64) {
        let (ephem_start, ephem_end, ephem_step) = (
            self.header.start_period,
            self.header.end_period,
            self.header.period_lenght,
        );
        if et < ephem_start || et > ephem_end {
            panic!("Time outside ephemeris range");
        }

        let mut nr = ((et - ephem_start) / ephem_step).floor() as usize;

        if (et - ephem_end).abs() < 1e-10 {
            nr -= 1;
        }

        let interval_start = (nr as f64) * ephem_step + ephem_start;
        let tau = (et - interval_start) / ephem_step;

        (nr, tau)
    }

    pub fn get_record_horizon(&self, body: u32, et: f64) -> Option<(&HorizonRecord, f64)> {
        let (nr, tau) = self.get_record_index(et);

        let records = &self.records[nr];
        let record_body = records
            .get(&body)
            .expect(format!("Failed to get record for body {} in block {}", body, nr).as_str());

        let ipt_body = self.header.ipt[body as usize];
        let n_subs = ipt_body[2] as usize;

        let sub_index = (tau * n_subs as f64).floor().min(n_subs as f64 - 1.0) as usize;

        let local_tau = tau * n_subs as f64 - sub_index as f64;
        Some((&record_body[sub_index], local_tau))
    }
}

#[cfg(test)]
mod test_horizon_reader {
    use super::*;
    use crate::jpl_ephem::download_jpl_file::{EphemFilePath, EphemFileSource};
    use crate::jpl_ephem::horizon::horizon_records::InterpResult;
    use hifitime::Epoch;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_reader_from_horizon() {
        let file_source: EphemFileSource = Some("horizon:DE440".try_into().unwrap()).unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = HorizonData::read_horizon_file(&file_path);

        assert_eq!(
            jpl_ephem.header,
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
            }
        );

        assert_eq!(jpl_ephem.records.len(), 12556);
        assert_eq!(
            jpl_ephem.records.iter().fold(0, |acc, x| acc + x.len()),
            150672
        );

        assert_eq!(
            &jpl_ephem.records[0].get(&0).unwrap()[0],
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
        let file_source: EphemFileSource = Some("horizon:DE440".try_into().unwrap()).unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = HorizonData::read_horizon_file(&file_path);

        let epoch1 = Epoch::from_mjd_in_time_scale(57028.479297592596, hifitime::TimeScale::TT);
        let (record, tau) = jpl_ephem
            .get_record_horizon(4, epoch1.to_jde_et_days())
            .unwrap();

        assert_eq!(tau, 0.639978049788624);

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
        let file_source: EphemFileSource = Some("horizon:DE440".try_into().unwrap()).unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = HorizonData::read_horizon_file(&file_path);

        let epoch1 = Epoch::from_mjd_in_time_scale(57028.479297592596, hifitime::TimeScale::TT);
        let (index, tau) = jpl_ephem.get_record_index(epoch1.to_jde_et_days());

        assert_eq!(index, 5307);
        assert_eq!(tau, 0.639978049788624);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_interpolation_from_horizon() {
        let file_source: EphemFileSource = Some("horizon:DE440".try_into().unwrap()).unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = HorizonData::read_horizon_file(&file_path);

        let epoch1 = Epoch::from_mjd_in_time_scale(57028.479297592596, hifitime::TimeScale::TT);
        let (record, tau) = jpl_ephem
            .get_record_horizon(10, epoch1.to_jde_et_days())
            .unwrap();

        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [428149.04652967455, -105270.23548367192, -68083.3417805149,].into(),
                velocity: Some([589.5451313541057, 729.3492107658788, 300.3651374866509,].into(),),
                acceleration: Some(
                    [-0.9192157891864692, 0.8829808730566571, 0.3898414406697089,].into(),
                ),
            }
        );

        let epoch1 = Epoch::from_mjd_in_time_scale(57049.231857592589, hifitime::TimeScale::TT);
        let (record, tau) = jpl_ephem
            .get_record_horizon(10, epoch1.to_jde_et_days())
            .unwrap();

        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [440183.15997859894, -89933.41046126859, -61760.6145039215].into(),
                velocity: Some([569.7900066764879, 749.1773000869151, 309.2398409158924].into()),
                acceleration: Some(
                    [-1.038549415912712, 0.9990469757280209, 0.45539282811508386].into()
                )
            }
        );

        let epoch1 = Epoch::from_mjd_in_time_scale(60781.51949044435, hifitime::TimeScale::TT);
        let (record, tau) = jpl_ephem
            .get_record_horizon(10, epoch1.to_jde_et_days())
            .unwrap();
        let res = record.interpolate(tau, true, true, 2);

        assert_eq!(
            res,
            InterpResult {
                position: [-742814.3000137291, -727671.3536906336, -288321.5373338285].into(),
                velocity: Some([1085.625632474908, -327.2648113002942, -162.1316622153563].into()),
                acceleration: Some(
                    [-0.06525161081610019, 1.2611973281426831, 0.5376907387399421].into()
                )
            }
        )
    }
}
