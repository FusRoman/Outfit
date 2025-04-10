mod daf_header;
mod directory;
mod ephemeris_record;
mod jpl_ephem_header;
mod summary_record;

use camino::Utf8Path;
use nom::number::complete::le_f64;
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

type daf_records = HashMap<(i32, i32), (Summary, Vec<EphemerisRecord>, DirectoryData)>

#[derive(Debug)]
pub struct JPLEphem {
    daf_header: DAFHeader,
    header: JPLEphemHeader,
    jpl_data: daf_records,
}

impl JPLEphem {
    pub fn read(jpl_path: &Utf8Path) -> Self {
        let mut file = BufReader::new(File::open(jpl_path).unwrap());

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
        JPLEphem {
            daf_header: daf_header,
            header: jpl_header,
            jpl_data: jpl_data,
        }
    }

    fn get_records(
        &self,
        target: NaifIds,
        center: NaifIds,
    ) -> Option<(&Summary, &Vec<EphemerisRecord>, &DirectoryData)> {
        self.jpl_data
            .get(&(target.to_id(), center.to_id()))
            .map(|(summary, records, dir_data)| (summary, records, dir_data))
    }
}

#[cfg(test)]
mod jpl_reader_test {
    use super::*;
    use crate::jpl_ephem::download_jpl_file::get_ephemeris_file;
    use std::io::BufReader;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_daf_header() {
        let version = Some("de440");
        let user_path = None;

        let file_path = get_ephemeris_file(version, user_path).unwrap();

        let mut file = BufReader::new(File::open(file_path).unwrap());
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
    fn test_jpl_ephem() {
        use crate::jpl_ephem::naif_ids::{
            planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary,
        };

        let version = Some("de440");
        let user_path = None;

        let file_path = get_ephemeris_file(version, user_path).unwrap();
        let jpl_ephem = JPLEphem::read(&file_path);

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

        let record_earth_sun = jpl_ephem
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
        )
    }
}
