use nom::{bytes::complete::take, number::complete::le_f64};
use std::{
    collections::HashMap,
    fs::File,
    io::{BufReader, Read, Seek},
};

use crate::jpl_ephem::{
    download_jpl_file::EphemFilePath, horizon::interpolation_result::InterpResult,
};

use super::{
    daf_header::DAFHeader, directory::DirectoryData, ephemeris_record::EphemerisRecord,
    jpl_ephem_header::JPLEphemHeader, naif_ids::NaifIds, summary_record::Summary,
};

type DafRecords = HashMap<(i32, i32), (Summary, Vec<EphemerisRecord>, DirectoryData)>;

#[derive(Debug, Clone)]
pub struct NaifData {
    daf_header: DAFHeader,
    header: JPLEphemHeader,
    jpl_data: DafRecords,
}

impl NaifData {
    /// Reads the JPL ephemeris file and parses the DAF header, JPL header, and ephemeris records.
    ///
    /// Arguments
    /// ---------
    /// * `file_path`: The path to the JPL ephemeris file.
    ///
    /// Returns
    /// -------
    /// * A `NaifData` instance containing the parsed data.
    pub fn read_naif_file(file_path: &EphemFilePath) -> Self {
        let mut file =
            BufReader::new(File::open(file_path.path()).unwrap_or_else(|_| {
                panic!("Failed to open the JPL ephemeris file: {}", file_path)
            }));

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

        let (_, jpl_header) = JPLEphemHeader::parse(ascii_comment)
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

        let mut jpl_data: DafRecords = HashMap::new();

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
        NaifData {
            daf_header,
            header: jpl_header,
            jpl_data,
        }
    }

    /// Retrieves the ephemeris records for the given target and center NAIF IDs.
    ///
    /// Arguments
    /// ---------
    /// * `target`: The target NAIF ID.
    /// * `center`: The center NAIF ID.
    ///
    /// Returns
    /// -------
    /// * An `Option` containing a tuple of the summary, ephemeris records, and directory data.
    /// If the target and center combination is not found, it returns `None`.
    fn get_records(
        &self,
        target: NaifIds,
        center: NaifIds,
    ) -> Option<(&Summary, &Vec<EphemerisRecord>, &DirectoryData)> {
        self.jpl_data
            .get(&(target.to_id(), center.to_id()))
            .map(|(summary, records, dir_data)| (summary, records, dir_data))
    }

    /// Retrieves the ephemeris record for the given target and center at the specified epoch.
    ///
    /// Arguments
    /// ---------
    /// * `target`: The target NAIF ID.
    /// * `center`: The center NAIF ID.
    /// * `et_seconds`: The epoch in seconds since the J2000 epoch.
    ///
    /// Returns
    /// -------
    /// * An `Option` containing a reference to the ephemeris record.
    /// If the target and center combination is not found, or if the epoch is out of range,
    /// it returns `None`.
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
    pub fn ephemeris(&self, target: NaifIds, center: NaifIds, et_seconds: f64) -> InterpResult {
        let record = self
            .get_record(target, center, et_seconds)
            .unwrap_or_else(|| {
                panic!(
                    "Failed to get ephemeris record for target: {:?}, center: {:?} at epoch: {}",
                    target, center, et_seconds
                )
            });

        let (position, velocity) = record.interpolate(et_seconds);
        InterpResult {
            position,
            velocity: Some(velocity),
            acceleration: None,
        }
    }

    /// Prints information about the available targets, centers, and epoch ranges in the ephemeris file.
    pub fn info(&self) {
        println!("+{:-^78}+", " Ephemeris File Information ");
        println!("{}", self.daf_header); // Use `Display` for `DAFHeader`
        println!("{}", self.header); // Use `Display` for `JPLEphemHeader`
        println!("+{:-^78}+", " Available Records ");

        for ((target, center), (summary, _, directory)) in &self.jpl_data {
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

#[cfg(test)]
mod test_naif_file {
    use super::*;

    #[cfg(feature = "jpl-download")]
    use crate::jpl_ephem::naif::naif_ids::{
        planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary,
    };
    #[cfg(feature = "jpl-download")]
    use crate::unit_test_global::JPL_EPHEM_NAIF;
    #[cfg(feature = "jpl-download")]
    use hifitime::Epoch;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_jpl_reader_from_naif() {
        assert_eq!(
            JPL_EPHEM_NAIF.daf_header,
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
            JPL_EPHEM_NAIF.header,
            JPLEphemHeader {
                version: "DE440".to_string(),
                creation_date: "25 June 2020".to_string(),
                start_ephem: "31-DEC-1549 00:00".to_string(),
                end_ephem: "25-JAN-2650 00:00".to_string(),
                start_jd: 2287184.5,
                end_jd: 2688976.5
            }
        );

        let record_earth_sun = JPL_EPHEM_NAIF
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
        let date_str = "2024-04-10T12:30:45";
        let epoch = Epoch::from_gregorian_str(date_str).unwrap();

        let record = JPL_EPHEM_NAIF
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

        let interp = JPL_EPHEM_NAIF.ephemeris(
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::SSB(SolarSystemBary::SSB),
            epoch1.to_et_seconds(),
        );

        assert_eq!(
            interp.to_au(),
            InterpResult {
                position: [[
                    -0.26169997917112875,
                    0.8682243128558709,
                    0.37623815974362734
                ]]
                .into(),
                velocity: Some(
                    [[
                        -3.8995165075699485e-7,
                        -9.957615232661774e-8,
                        -4.316895883931796e-8
                    ]]
                    .into()
                ),
                acceleration: None
            }
        );

        let epoch2 = Epoch::from_mjd_in_time_scale(57049.231857592589, hifitime::TimeScale::TT);
        let interp = JPL_EPHEM_NAIF.ephemeris(
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::SSB(SolarSystemBary::SSB),
            epoch2.to_et_seconds(),
        );
        assert_eq!(
            interp.to_au(),
            InterpResult {
                position: [[-0.5860307419898751, 0.7233961430776997, 0.31345193147254585]].into(),
                velocity: Some(
                    [[
                        -3.2554490509465264e-7,
                        -2.1982148078907505e-7,
                        -9.529706060142567e-8
                    ]]
                    .into()
                ),
                acceleration: None
            }
        );
    }
}
