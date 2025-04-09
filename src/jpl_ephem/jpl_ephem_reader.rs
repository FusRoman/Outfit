use camino::Utf8Path;
use chrono::Duration;
use chrono::NaiveDateTime;
use nom::number::complete::le_f64;
use nom::number::complete::le_i32;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;
use std::io::Seek;

use nom::{
    bytes::complete::{tag, take, take_until},
    character::complete::{line_ending, not_line_ending},
    number::complete::double,
    IResult, Parser,
};

#[derive(Debug, PartialEq)]
struct DAFHeader {
    idword: String,
    internal_filename: String,
    nd: i32,
    ni: i32,
    fward: i32,
    bward: i32,
    free: i32,
    locfmt: String,
    fptstr: String,
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

fn parse_daf_header(input: &[u8]) -> IResult<&[u8], DAFHeader> {
    let (input, id_word) = take(8usize)(input)?; // "DAF/SPK "
    let (input, nd_bytes) = le_i32(input)?; // ND
    let (input, ni_bytes) = le_i32(input)?; // NI
    let (input, ifname) = take(60usize)(input)?; // internal file name
    let (input, fwd) = le_i32(input)?; // forward ptr
    let (input, bwd) = le_i32(input)?; // backward ptr
    let (input, free) = le_i32(input)?; // first free address
    let (input, locfmt) = take(8usize)(input)?; // location format
    let (input, _) = take(603usize)(input)?; // reserved
    let (input, ftpstr) = take(28usize)(input)?; // ftp string
    Ok((
        input,
        DAFHeader {
            idword: String::from_utf8_lossy(id_word).trim().to_string(),
            internal_filename: String::from_utf8_lossy(ifname).trim().to_string(),
            nd: nd_bytes,
            ni: ni_bytes,
            fward: fwd,
            bward: bwd,
            free: free,
            locfmt: String::from_utf8_lossy(locfmt).trim().to_string(),
            fptstr: String::from_utf8_lossy(ftpstr).trim().to_string(),
        },
    ))
}

#[derive(Debug, PartialEq)]
struct JPLEphemHeader {
    version: String,
    creation_date: String,
    start_ephem: String,
    end_ephem: String,
    start_jd: f64,
    end_jd: f64,
}

fn parse_version(input: &str) -> IResult<&str, &str> {
    let (input, _) = take_until("JPL planetary and lunar ephemeris")(input)?;
    let (input, _) = tag("JPL planetary and lunar ephemeris ")(input)?;
    let (input, version) = not_line_ending(input)?;
    let (input, _) = line_ending(input)?;

    Ok((input, version.trim()))
}

fn parse_creation_date(input: &str) -> IResult<&str, &str> {
    let (input, _) = take_until("Integrated ")(input)?;
    let (input, _) = tag("Integrated ")(input)?;
    let (input, creation_date) = not_line_ending(input)?;
    let (input, _) = line_ending(input)?;

    Ok((input, creation_date.trim()))
}

fn parse_date_range(input: &str) -> IResult<&str, (&str, &str)> {
    let (input, _) = take_until("Time span covered by ephemeris:")(input)?;
    let (input, _) = tag("Time span covered by ephemeris:")(input)?;
    let (input, _) = (line_ending, line_ending).parse(input)?;
    let (input, (start_ephem, _, end_ephem)) =
        (take_until("to"), tag("to   "), not_line_ending).parse(input)?;

    Ok((input, (start_ephem.trim(), end_ephem.trim())))
}

fn parse_jd_range(input: &str) -> IResult<&str, (f64, f64)> {
    let (input, (_, _, start_jd, _, end_jd)) = (
        line_ending,
        tag("JD   "),
        |s| double(s),
        tag("   to   JD   "),
        |s| double(s),
    )
        .parse(input)?;

    Ok((input, (start_jd, end_jd)))
}

fn parse_ephem_header(input: &str) -> IResult<&str, JPLEphemHeader> {
    let (input, version) = parse_version(input)?;
    let (input, creation_date) = parse_creation_date(input)?;

    let (input, (start_ephem, end_ephem)) = parse_date_range(input)?;
    let (input, (start_jd, end_jd)) = parse_jd_range(input)?;
    Ok((
        input,
        JPLEphemHeader {
            version: version.to_string(),
            creation_date: creation_date.trim().to_string(),
            start_ephem: start_ephem.trim().to_string(),
            end_ephem: end_ephem.trim().to_string(),
            start_jd,
            end_jd,
        },
    ))
}

#[derive(Debug)]
struct Summary {
    start_epoch: f64,
    end_epoch: f64,
    target: i32,
    center: i32,
    frame_id: i32,
    data_type: i32,
    initial_addr: i32,
    final_addr: i32,
}
fn parse_summary(input: &[u8]) -> IResult<&[u8], Summary> {
    let (input, start_epoch) = le_f64(input)?;
    let (input, end_epoch) = le_f64(input)?;

    let (input, target) = le_i32(input)?;
    let (input, center) = le_i32(input)?;
    let (input, frame_id) = le_i32(input)?;
    let (input, data_type) = le_i32(input)?;
    let (input, initial_addr) = le_i32(input)?;
    let (input, final_addr) = le_i32(input)?;
    Ok((
        input,
        Summary {
            start_epoch,
            end_epoch,
            target: target,
            center: center,
            frame_id: frame_id,
            data_type: data_type,
            initial_addr: initial_addr,
            final_addr: final_addr,
        },
    ))
}

#[derive(Debug)]
struct EphemerisRecord {
    mid: f64,
    radius: f64,
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
}

fn parse_ephemeris_record(
    input: &[u8],
    ncoeff: usize,
) -> Result<EphemerisRecord, Box<dyn std::error::Error>> {
    let mut offset = 0;

    // Helper closure pour lire un f64 et avancer dans le buffer
    let mut read_f64 = || -> Result<f64, Box<dyn std::error::Error>> {
        let val = f64::from_le_bytes(input[offset..offset + 8].try_into().unwrap());
        offset += 8;
        Ok(val)
    };

    let mid = read_f64()?;
    let radius = read_f64()?;

    let mut x = Vec::with_capacity(ncoeff);
    let mut y = Vec::with_capacity(ncoeff);
    let mut z = Vec::with_capacity(ncoeff);

    for _ in 0..ncoeff {
        x.push(read_f64()?);
    }
    for _ in 0..ncoeff {
        y.push(read_f64()?);
    }
    for _ in 0..ncoeff {
        z.push(read_f64()?);
    }

    Ok(EphemerisRecord {
        mid,
        radius,
        x,
        y,
        z,
    })
}

fn parse_segment_records(
    file: &mut BufReader<File>,
    segment_start_addr: usize,
    rsize: usize,
    n_records: usize,
) -> Vec<EphemerisRecord> {
    let mut records = Vec::with_capacity(n_records);

    let record_byte_size = rsize * 8;
    let mut buf = vec![0u8; record_byte_size];

    // L'offset du premier record
    let start_byte_offset = (segment_start_addr - 1) * 8;

    let n_coeffs = (rsize - 2) / 3;

    for i in 0..n_records {
        let byte_offset = start_byte_offset + i * record_byte_size;

        file.seek(std::io::SeekFrom::Start(byte_offset as u64))
            .unwrap();
        file.read_exact(&mut buf).unwrap();

        let input = buf.as_slice();

        // MID, RADIUS
        let record = parse_ephemeris_record(input, n_coeffs).unwrap();

        records.push(record);
    }

    records
}

fn format_record_bounds(mid: f64, radius: f64) -> (String, String, String) {
    let j2000_epoch = NaiveDateTime::new(
        chrono::NaiveDate::from_ymd_opt(2000, 1, 1).unwrap(),
        chrono::NaiveTime::from_hms_opt(12, 0, 0).unwrap(),
    );

    let format_time = |seconds: f64| {
        let secs = seconds.trunc() as i64;
        let nanos = ((seconds.fract()) * 1_000_000_000.0).round() as i64;
        (j2000_epoch + Duration::seconds(secs) + Duration::nanoseconds(nanos))
            .format("%Y-%m-%d %H:%M:%S")
            .to_string()
    };

    let start = mid - radius;
    let end = mid + radius;

    (format_time(start), format_time(mid), format_time(end))
}

#[derive(Debug)]
struct DirectoryData {
    init: f64,
    intlen: f64,
    rsize: f64,
    n_records: f64,
}

fn parse_directory_record(file: &mut BufReader<File>, end_addr: usize) -> DirectoryData {
    let directory_offset_bytes = (end_addr - 4) * 8;
    let mut dir_buf = [0u8; 32]; // 4 f64 = 32 octets
    file.seek(std::io::SeekFrom::Start(directory_offset_bytes as u64))
        .unwrap();
    file.read_exact(&mut dir_buf).unwrap();
    let (_, init) = le_f64::<_, nom::error::Error<_>>(&dir_buf[0..8]).unwrap();
    let (_, intlen) = le_f64::<_, nom::error::Error<_>>(&dir_buf[8..16]).unwrap();
    let (_, rsize) = le_f64::<_, nom::error::Error<_>>(&dir_buf[16..24]).unwrap();
    let (_, n_records) = le_f64::<_, nom::error::Error<_>>(&dir_buf[24..32]).unwrap();
    DirectoryData {
        init,
        intlen,
        rsize,
        n_records,
    }
}

#[derive(Debug)]
struct JPLEphem {
    daf_header: DAFHeader,
    header: JPLEphemHeader,
    jpl_data: HashMap<(i32, i32), (Summary, Vec<EphemerisRecord>, DirectoryData)>,
}

impl JPLEphem {
    fn new(jpl_path: &Utf8Path) -> Self {
        let mut file = BufReader::new(File::open(jpl_path).unwrap());

        let mut buffer = [0u8; 1 << 10];
        file.read_exact(&mut buffer)
            .expect("Failed to read the DAF header. (first 1024 bytes)");

        let (_, daf_header) =
            parse_daf_header(&buffer).expect("Failed to parse the DAF header with nom !");

        let end_ascii_comment = (daf_header.fward as usize - 1) * 1024;
        let mut ascii_buffer = vec![0u8; end_ascii_comment];
        file.read_exact(&mut ascii_buffer)
            .expect("Failed to read the ASCII comment area.");

        let binding = String::from_utf8_lossy(&ascii_buffer).replace('\0', "\n");
        let ascii_comment = binding.as_str();

        let (_, jpl_header) =
            parse_ephem_header(&ascii_comment).expect("Failed to parse the JPL header with nom !");

        let offset_bytes = (daf_header.fward as usize - 1) * 1024;
        file.seek(std::io::SeekFrom::Start(offset_bytes as u64))
            .unwrap();

        let mut buffer = [0u8; 1 << 10];
        file.read_exact(&mut buffer).unwrap();

        let (input, _) =
            take::<_, _, nom::error::Error<&[u8]>>(16usize)(buffer.as_slice()).unwrap();
        let (_, nsum) = le_f64::<_, nom::error::Error<_>>(input).unwrap();

        // Taille d'un summary
        let ss = daf_header.nd as usize + ((daf_header.ni as usize + 1) / 2);

        let mut jpl_data: HashMap<(i32, i32), (Summary, Vec<EphemerisRecord>, DirectoryData)> =
            HashMap::new();

        // Lire les summaries
        for i in 0..(nsum as usize) {
            let start = 24 + i * ss * 8;
            let end = start + ss * 8;

            let summary_bytes = &buffer[start..end];

            let (_, summary) =
                parse_summary(summary_bytes).expect("Failed to parse the summary with nom !");

            let dir_data = parse_directory_record(&mut file, summary.final_addr as usize);

            let records = parse_segment_records(
                &mut file,
                summary.initial_addr as usize,
                dir_data.rsize as usize,
                dir_data.n_records as usize,
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
        let (_, daf_header) = parse_daf_header(&buffer).unwrap();

        assert_eq!(daf_header.idword, "DAF/SPK");
        assert_eq!(daf_header.internal_filename, "NIO2SPK");
        assert_eq!(daf_header.fward, 62);
        assert_eq!(daf_header.bward, 62);
        assert_eq!(daf_header.free, 14974889);
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_new_jpl_ephem() {
        let version = Some("de440");
        let user_path = None;

        let file_path = get_ephemeris_file(version, user_path).unwrap();
        let jpl_ephem = JPLEphem::new(&file_path);

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
    }
}
