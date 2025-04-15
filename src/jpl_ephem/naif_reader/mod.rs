mod daf_header;
mod directory;
mod ephemeris_record;
mod jpl_ephem_header;
mod summary_record;

use camino::Utf8Path;
use core::num;
use nalgebra::Vector3;
use nom::error::ErrorKind;
use nom::multi::count;
use nom::number::complete::{be_f64, be_i32, le_f64, le_i32, le_u32};
use nom::{IResult, Input, Parser};
use parquet::{file, record};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek};

use nom::bytes::complete::take;

use crate::jpl_ephem::{
    naif_ids::NaifIds,
    naif_reader::{
        daf_header::DAFHeader, directory::DirectoryData, ephemeris_record::EphemerisRecord,
        jpl_ephem_header::JPLEphemHeader, summary_record::Summary,
    },
};

use super::download_jpl_file::{EphemFilePath, JPLHorizonVersion};

type daf_records = HashMap<(i32, i32), (Summary, Vec<EphemerisRecord>, DirectoryData)>;

type horizon_records = Vec<HashMap<u32, Vec<EphemerisRecord>>>;

/// Structure représentant les 15 triplets IPT
pub type IPT = [[u32; 3]; 15];

#[derive(Debug)]
pub struct NaifData {
    daf_header: DAFHeader,
    header: JPLEphemHeader,
    jpl_data: daf_records,
}

#[derive(Debug, PartialEq)]
pub struct HorizonHeader {
    jpl_version: String,
    ipt: [[u32; 3]; 15],
    start_period: f64,
    end_period: f64,
    period_lenght: f64,
    ksize: usize,
}

#[derive(Debug)]
pub struct HorizonData {
    header: HorizonHeader,
    records: horizon_records,
}

#[derive(Debug)]
pub enum JPLEphem {
    NAIF(NaifData),
    JPLHorizon(HorizonData),
}

pub fn print_hex_dump(data: &[u8]) {
    const BYTES_PER_LINE: usize = 40;

    for (i, chunk) in data.chunks(BYTES_PER_LINE).enumerate() {
        // Offset
        print!("{:08x}: ", i * BYTES_PER_LINE);

        // Hex representation
        for byte in chunk {
            print!("{:02x} ", byte);
        }

        // Padding if chunk is not full
        if chunk.len() < BYTES_PER_LINE {
            for _ in 0..(BYTES_PER_LINE - chunk.len()) {
                print!("   ");
            }
        }

        // ASCII representation
        print!("|");
        for byte in chunk {
            let ascii = if byte.is_ascii_graphic() || *byte == b' ' {
                *byte as char
            } else {
                '.'
            };
            print!("{}", ascii);
        }
        println!("|");
    }
}

impl JPLEphem {
    fn read_naif_file(file_path: &Utf8Path) -> Self {
        let mut file = BufReader::new(
            File::open(file_path)
                .expect(format!("Failed to open the JPL ephemeris file: {}", file_path).as_str()),
        );

        let mut buffer = [0u8; 1 << 10];
        file.read_exact(&mut buffer)
            .expect("Failed to read the DAF header. (first 1024 bytes)");

        let (_, daf_header) =
            DAFHeader::parse(&buffer).expect("Failed to parse the DAF header with nom !");

        let end_ascii_comment = (daf_header.fward as usize - 1) * 1024;
        let mut ascii_buffer = vec![0u8; end_ascii_comment];
        file.read_exact(&mut ascii_buffer)
            .expect("Failed to read the ASCII comment area.");

        let binding = String::from_utf8_lossy(&ascii_buffer).replace('\0', "\n");
        let ascii_comment = binding.as_str();

        let (_, jpl_header) = JPLEphemHeader::parse(&ascii_comment)
            .expect("Failed to parse the JPL header with nom !");

        let offset_bytes = (daf_header.fward as usize - 1) * 1024;
        file.seek(std::io::SeekFrom::Start(offset_bytes as u64))
            .unwrap();

        let mut buffer = [0u8; 1 << 10];
        file.read_exact(&mut buffer).unwrap();

        let (input, _) =
            take::<_, _, nom::error::Error<&[u8]>>(16usize)(buffer.as_slice()).unwrap();

        // nsum is the number of summary records
        let (_, nsum) = le_f64::<_, nom::error::Error<_>>(input).unwrap();

        // ss is the size of the summary record
        let ss = daf_header.nd as usize + ((daf_header.ni as usize + 1) / 2);

        let mut jpl_data: daf_records = HashMap::new();

        // Read the summary records and the element records
        for i in 0..(nsum as usize) {
            let start = 24 + i * ss * 8;
            let end = start + ss * 8;

            let summary_bytes = &buffer[start..end];

            let (_, summary) =
                Summary::parse(summary_bytes).expect("Failed to parse the summary with nom !");

            let dir_data = DirectoryData::parse(&mut file, summary.final_addr as usize);

            let records = EphemerisRecord::parse(
                &mut file,
                summary.initial_addr as usize,
                dir_data.rsize,
                dir_data.n_records,
            );
            jpl_data.insert(
                (summary.target, summary.center),
                (summary, records, dir_data),
            );
        }
        JPLEphem::NAIF(NaifData {
            daf_header: daf_header,
            header: jpl_header,
            jpl_data: jpl_data,
        })
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
    pub fn compute_recsize(ipt: [[u32; 3]; 15]) -> usize {
        let mut kernel_size = 4; // mots de 32 bits

        for i in 0..15 {
            let n_subintervals = ipt[i][1] as usize;
            let n_coeffs = ipt[i][2] as usize;
            let dim = Self::dimension(i);

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

    pub fn parse_ipt_13_14(input: &[u8]) -> IResult<&[u8], ([[u32; 3]; 2])> {
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

        dbg!(offset);

        file.seek(std::io::SeekFrom::Start(offset))?;
        let mut buffer = [0u8; 24]; // 6 * 4
        file.read_exact(&mut buffer)?;

        let (_, ipt_extra) = Self::parse_ipt_13_14(&buffer).map_err(|_| {
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

    pub fn extract_body_records(block: &[f64], ipt: [u32; 3]) -> Vec<EphemerisRecord> {
        let offset = ipt[0] as usize;
        let n_coeffs = ipt[1] as usize;
        let n_subs = ipt[2] as usize;

        let jd_start = block[0];
        let jd_end = block[1];
        let sub_span = (jd_end - jd_start) / n_subs as f64;

        let coeffs = &block[2..]; // skip JD start/end

        let mut records = Vec::with_capacity(n_subs);

        for i in 0..n_subs {
            let mut x = Vec::with_capacity(n_coeffs);
            let mut y = Vec::with_capacity(n_coeffs);
            let mut z = Vec::with_capacity(n_coeffs);

            for j in 0..n_coeffs {
                let base = offset - 1 + (i * 3 * n_coeffs);
                let x_val = coeffs.get(base + j).copied().expect(
                    format!(
                        "Failed to get x coefficient at index {} for body {}",
                        base + j,
                        i
                    )
                    .as_str(),
                );
                let y_val = coeffs.get(base + j + n_coeffs).copied().expect(
                    format!(
                        "Failed to get y coefficient at index {} for body {}",
                        base + j + n_coeffs,
                        i
                    )
                    .as_str(),
                );
                let z_val = coeffs.get(base + j + 2 * n_coeffs).copied().expect(
                    format!(
                        "Failed to get z coefficient at index {} for body {}",
                        base + j + 2 * n_coeffs,
                        i
                    )
                    .as_str(),
                );

                x.push(x_val);
                y.push(y_val);
                z.push(z_val);
            }

            let start = jd_start + i as f64 * sub_span;
            let end = start + sub_span;
            let mid = (start + end) / 2.0;
            let radius = sub_span / 2.0;

            records.push(EphemerisRecord {
                mid,
                radius,
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
    ) -> std::io::Result<horizon_records> {
        // Skip les 2 premiers blocs (header texte + header binaire)
        let offset = 2 * recsize;
        file.seek(std::io::SeekFrom::Start(offset as u64))?;

        // Lis tout le reste du fichier
        let mut buffer = Vec::new();
        file.read_to_end(&mut buffer)?;

        // Vérifie que la taille est un multiple du recsize
        let n_blocks = buffer.len() / recsize;
        let mut blocks: horizon_records = Vec::with_capacity(n_blocks);

        for i in 0..n_blocks {
            let start = i * recsize;
            let end = start + recsize;
            let slice = &buffer[start..end];

            let (_, block) = Self::parse_block(slice, ncoeff).map_err(|_| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse block {i}"),
                )
            })?;

            let mut records = HashMap::new();
            for body in 0..12 {
                let ipt_body = ipt[body];

                let body_rec = Self::extract_body_records(&block, ipt_body);
                records.insert(body as u32, body_rec);
            }
            blocks.push(records);
        }

        Ok(blocks)
    }

    fn read_horizon_file(file_path: &Utf8Path, version: &JPLHorizonVersion) -> Self {
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
            File::open(file_path)
                .expect(format!("Failed to open the JPL ephemeris file: {}", file_path).as_str()),
        );

        let mut file_data = vec![0u8; 1 << 12];
        file.read_exact(&mut file_data)
            .expect("Failed to read the JPL ephemeris file.");

        let (input, ttl) = parse_ttl(&file_data).expect("failed to parse TTL");

        let (input, cnam) = parse_cnam(input).expect("failed to parse CNAM");

        let (input, ss) = count(le_f64::<_, nom::error::Error<_>>, 3)
            .parse(input)
            .expect("failed to parse SS");

        let (_, ncon) = le_i32::<_, nom::error::Error<_>>(input).expect("failed to parse NCON");

        file.seek(std::io::SeekFrom::Start(2696))
            .expect("Failed to seek to the JPL ephemeris file IPT.");
        let mut buffer = [0u8; 204]; // JPL_HEADER_SIZE
        file.read_exact(&mut buffer)
            .expect("Failed to read the JPL ephemeris file IPT.");

        let (_, (mut ipt, numde)) = Self::parse_ipt(&buffer).expect("failed to parse IPT header");

        let ipt_extra = Self::read_ipt_13_14(&mut file, ncon as u32, version)
            .expect("Failed to read IPT[13] and IPT[14]");

        ipt[13] = ipt_extra
            .as_ref()
            .map(|ipt_extra| ipt_extra[0])
            .unwrap_or([0; 3]);

        ipt[14] = ipt_extra
            .as_ref()
            .map(|ipt_extra| ipt_extra[1])
            .unwrap_or([0; 3]);

        let recsize = JPLEphem::compute_recsize(ipt);

        let ncoeff = recsize / 8;
        let blocks = Self::parse_all_blocks(&mut file, recsize, ncoeff, ipt).unwrap();

        JPLEphem::JPLHorizon(HorizonData {
            header: HorizonHeader {
                jpl_version: format!("DE{}", numde),
                ipt: ipt,
                start_period: ss[0],
                end_period: ss[1],
                period_lenght: ss[2],
                ksize: 0,
            },
            records: blocks,
        })
    }

    pub fn read(ephem_file: &EphemFilePath) -> Self {
        match ephem_file {
            EphemFilePath::NAIF(file_path, version) => JPLEphem::read_naif_file(file_path),
            EphemFilePath::JPLHorizon(file_path, version) => {
                JPLEphem::read_horizon_file(file_path, version)
            }
        }
    }

    fn get_records(
        &self,
        target: NaifIds,
        center: NaifIds,
    ) -> Option<(&Summary, &Vec<EphemerisRecord>, &DirectoryData)> {
        match self {
            JPLEphem::NAIF(naif_data) => naif_data
                .jpl_data
                .get(&(target.to_id(), center.to_id()))
                .map(|(summary, records, dir_data)| (summary, records, dir_data)),
            JPLEphem::JPLHorizon(_) => None,
        }
    }

    fn get_record(
        &self,
        target: NaifIds,
        center: NaifIds,
        et_seconds: f64,
    ) -> Option<&EphemerisRecord> {
        let (summary, records, directory) = self.get_records(target, center)?;
        if et_seconds < summary.start_epoch || et_seconds > summary.end_epoch {
            return None;
        }

        let idx = ((et_seconds - directory.init) / directory.intlen as f64).floor() as usize;

        records.get(idx)
    }

    /// Interpolates the ephemeris record for the given target and center at the specified epoch.
    ///
    /// Arguments
    /// ---------
    /// * `target`: The target NAIF ID.
    /// * `center`: The center NAIF ID.
    /// * `et_seconds`: The epoch in seconds since the J2000 epoch.
    ///
    /// Returns
    /// --------
    /// * A tuple containing the position and velocity vectors in the J2000 frame.
    /// The position and velocity vectors are in kilometers and kilometers per second, respectively.
    pub fn ephemeris_prediction(
        &self,
        target: NaifIds,
        center: NaifIds,
        et_seconds: f64,
    ) -> (Vector3<f64>, Vector3<f64>) {
        let record = self.get_record(target, center, et_seconds).expect(
            format!(
                "Failed to get ephemeris record for target: {:?}, center: {:?} at epoch: {}",
                target, center, et_seconds
            )
            .as_str(),
        );

        record.interpolate(et_seconds)
    }

    /// Prints information about the available targets, centers, and epoch ranges in the ephemeris file.
    pub fn info(&self) {
        if let JPLEphem::NAIF(naif_data) = self {
            println!("+{:-^78}+", " Ephemeris File Information ");
            println!("{}", naif_data.daf_header); // Use `Display` for `DAFHeader`
            println!("{}", naif_data.header); // Use `Display` for `JPLEphemHeader`
            println!("+{:-^78}+", " Available Records ");

            for ((target, center), (summary, _, directory)) in &naif_data.jpl_data {
                println!(
                    "+{:-^78}+",
                    format!(" Target: {}, Center: {} ", target, center)
                );
                println!("{}", summary); // Use `Display` for `Summary`
                println!(
                    "| {:<76} |",
                    format!(
                        "Directory Data: Init = {:.2}, Interval Length = {:.2}, Records = {}",
                        directory.init, directory.intlen, directory.n_records
                    )
                );
                println!("+{:-^78}+", "");
            }
        }
    }
}

#[cfg(test)]
mod jpl_reader_test {
    use super::*;
    use std::io::BufReader;

    use hifitime::Epoch;

    use crate::jpl_ephem::{
        download_jpl_file::{EphemFilePath, EphemFileSource},
        naif_ids::{planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds},
        naif_reader::JPLEphem,
    };

    use crate::constants::AU;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_daf_header() {
        let file_source: EphemFileSource = "naif:DE440".try_into().unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();

        let mut file = BufReader::new(File::open(file_path.path()).unwrap());
        let mut buffer = [0u8; 1024];
        file.read_exact(&mut buffer).unwrap();
        let (_, daf_header) = DAFHeader::parse(&buffer).unwrap();

        assert_eq!(daf_header.idword, "DAF/SPK");
        assert_eq!(daf_header.internal_filename, "NIO2SPK");
        assert_eq!(daf_header.fward, 62);
        assert_eq!(daf_header.bward, 62);
        assert_eq!(daf_header.free, 14974889);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_reader_from_naif() {
        let file_source: EphemFileSource = "naif:DE440".try_into().unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_data = JPLEphem::read(&file_path);

        let JPLEphem::NAIF(jpl_ephem) = &jpl_data else {
            panic!("Expected JPL ephemeris data");
        };

        assert_eq!(
            jpl_ephem.daf_header,
            DAFHeader {
                idword: "DAF/SPK".to_string(),
                internal_filename: "NIO2SPK".to_string(),
                nd: 2,
                ni: 6,
                fward: 62,
                bward: 62,
                free: 14974889,
                locfmt: "LTL-IEEE".to_string(),
                fptstr: "FTPSTR:\r:\n:\r\n:\r\0:�:\u{10}�:ENDFTP".to_string()
            }
        );

        assert_eq!(
            jpl_ephem.header,
            JPLEphemHeader {
                version: "DE440".to_string(),
                creation_date: "25 June 2020".to_string(),
                start_ephem: "31-DEC-1549 00:00".to_string(),
                end_ephem: "25-JAN-2650 00:00".to_string(),
                start_jd: 2287184.5,
                end_jd: 2688976.5
            }
        );

        let record_earth_sun = jpl_data
            .get_records(
                NaifIds::PB(PlanetaryBary::EarthMoon),
                NaifIds::SSB(SolarSystemBary::SSB),
            )
            .unwrap();

        assert_eq!(
            record_earth_sun.0,
            &Summary {
                start_epoch: -14200747200.0,
                end_epoch: 20514081600.0,
                target: 3,
                center: 0,
                frame_id: 1,
                data_type: 2,
                initial_addr: 3021513,
                final_addr: 4051108
            }
        );

        assert_eq!(
            record_earth_sun.1.len(),
            25112,
            "Expected 25112 records for Earth-Moon barycenter"
        );

        assert_eq!(
            record_earth_sun.2,
            &DirectoryData {
                init: -14200747200.0,
                intlen: 1382400,
                rsize: 41,
                n_records: 25112
            }
        );

        let first_tchebychev_coeff = record_earth_sun.1[0].clone();
        assert_eq!(
            first_tchebychev_coeff,
            EphemerisRecord {
                mid: -14200056000.0,
                radius: 691200.0,
                x: vec![
                    -59117487.054044664,
                    -19163216.532728795,
                    291991.27938009636,
                    15847.329699283478,
                    -133.03948110729542,
                    -4.459284869049275,
                    0.03379900481247174,
                    0.0011716375873243507,
                    1.4852185006919311e-5,
                    -1.1096435596643423e-5,
                    -3.3277738887706986e-6,
                    -1.7115381406088932e-7,
                    1.8759402940064767e-7
                ],
                y: vec![
                    122675629.90130842,
                    -7590835.1393347755,
                    -614598.6274020968,
                    6451.054990886844,
                    265.301332213569,
                    -1.9960672202990268,
                    -0.05594323453310299,
                    0.00034547006847617906,
                    -7.545272774250726e-5,
                    -4.915547872121608e-6,
                    2.9685818368805314e-6,
                    9.978975536291768e-7,
                    4.606523495199457e-8
                ],
                z: vec![
                    53352759.40834735,
                    -3301893.184660649,
                    -267186.62682516687,
                    2806.0081419486182,
                    115.33401330537684,
                    -0.8678118085367849,
                    -0.02423486566672119,
                    0.00012129101423135178,
                    -4.0581575470237784e-5,
                    -1.4813269496375633e-6,
                    1.7525667147619721e-6,
                    4.996616729595978e-7,
                    8.008233021055763e-9
                ]
            },
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_get_record() {
        let file_source: EphemFileSource = "naif:DE440".try_into().unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = JPLEphem::read(&file_path);

        let date_str = "2024-04-10T12:30:45";
        let epoch = Epoch::from_gregorian_str(date_str).unwrap();

        let record = jpl_ephem
            .get_record(
                NaifIds::PB(PlanetaryBary::EarthMoon),
                NaifIds::SSB(SolarSystemBary::SSB),
                epoch.to_et_seconds(),
            )
            .unwrap();

        assert_eq!(
            record.clone(),
            EphemerisRecord {
                mid: 765806400.0,
                radius: 691200.0,
                x: vec![
                    -142642883.69463348,
                    6126046.973280299,
                    669712.3503474531,
                    -5579.567590155643,
                    -253.10830382052777,
                    2.0009704562042607,
                    0.03311672381410127,
                    -0.00019248035934527232,
                    -0.00013487467910304937,
                    -2.0877032284606548e-5,
                    9.494283770641653e-6,
                    1.4492164476747281e-6,
                    -6.919414478143085e-7
                ],
                y: vec![
                    -43487170.26426978,
                    -17970535.021323554,
                    203314.02126174027,
                    13906.974989740016,
                    -103.50456842568748,
                    -2.9828833120574494,
                    0.031229191202504284,
                    -4.8443077899008055e-5,
                    -7.758663306889506e-5,
                    1.4377880407035774e-5,
                    4.749883816396439e-6,
                    -1.0111886939662983e-6,
                    -3.469900956265658e-7
                ],
                z: vec![
                    -18816839.385616597,
                    -7790117.51618049,
                    88129.26429267193,
                    6028.524885779582,
                    -44.86708885740644,
                    -1.2926124793403566,
                    0.013597938686248697,
                    -6.981467955749496e-5,
                    -3.8067570619539636e-5,
                    9.97884013115029e-6,
                    2.3463112012575143e-6,
                    -7.071057213710581e-7,
                    -1.7106237886381444e-7
                ]
            }
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_ephemeris() {
        let epoch1 = Epoch::from_mjd_in_time_scale(57028.479297592596, hifitime::TimeScale::TT);
        let file_source: EphemFileSource = "naif:DE440".try_into().unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_ephem = JPLEphem::read(&file_path);

        let (position, velocity) = jpl_ephem.ephemeris_prediction(
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::SSB(SolarSystemBary::SSB),
            epoch1.to_et_seconds(),
        );

        assert_eq!(
            position / AU,
            Vector3::new(
                -0.26169997917112875,
                0.8682243128558709,
                0.37623815974362734
            )
        );

        assert_eq!(
            velocity / AU,
            Vector3::new(
                -3.8995165075699485e-7,
                -9.957615232661774e-8,
                -4.316895883931796e-8
            )
        );

        let epoch2 = Epoch::from_mjd_in_time_scale(57049.231857592589, hifitime::TimeScale::TT);
        let (position, velocity) = jpl_ephem.ephemeris_prediction(
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::SSB(SolarSystemBary::SSB),
            epoch2.to_et_seconds(),
        );
        assert_eq!(
            position / AU,
            Vector3::new(-0.5860307419898751, 0.7233961430776997, 0.31345193147254585)
        );

        assert_eq!(
            velocity / AU,
            Vector3::new(
                -3.2554490509465264e-7,
                -2.1982148078907505e-7,
                -9.529706060142567e-8
            )
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_reader_from_horizon() {
        let file_source: EphemFileSource = Some("horizon:DE440".try_into().unwrap()).unwrap();

        let file_path = EphemFilePath::get_ephemeris_file(file_source).unwrap();
        let jpl_data = JPLEphem::read(&file_path);

        let JPLEphem::JPLHorizon(jpl_ephem) = &jpl_data else {
            panic!("Expected JPL ephemeris data");
        };

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
                ksize: 0,
            }
        );

        assert_eq!(jpl_ephem.records.len(), 12556);
        assert_eq!(
            jpl_ephem.records.iter().fold(0, |acc, x| acc + x.len()),
            150672
        );

        assert_eq!(
            &jpl_ephem.records[0].get(&0).unwrap()[0],
            &EphemerisRecord {
                mid: 2287188.5,
                radius: 4.0,
                x: vec![
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
                    -1.118766077984058e-7,
                    19369902.32537584,
                    -12637978.7311934,
                ],
                y: vec![
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
                    -2.579709214481493e-7,
                    15085546.93080799,
                    -5541484.32081567,
                ],
                z: vec![
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
                    -1.260377702474755e-7,
                    -58449689.86096714,
                    -1800104.423508977,
                ],
            }
        );
    }
}
