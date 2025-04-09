use camino::Utf8Path;
use nom::number::complete::le_i32;
use std::fs::File;
use std::io::BufReader;
use std::io::Read;

use nom::{
    bytes::complete::{tag, take, take_until},
    character::complete::{line_ending, not_line_ending},
    number::complete::double,
    IResult, Parser,
};

#[derive(Debug)]
struct DAFHeader {
    idword: String,
    internal_filename: String,
    fward: i32,
    bward: i32,
    free: i32,
}

fn parse_daf_header(input: &[u8]) -> IResult<&[u8], DAFHeader> {
    let (input, id_word) = take(8usize)(input)?; // "DAF/SPK "
    let (input, nd_bytes) = le_i32(input)?; // ND
    let (input, ni_bytes) = le_i32(input)?; // NI
    let (input, ifname) = take(60usize)(input)?; // internal file name
    let (input, fwd) = le_i32(input)?; // forward ptr
    let (input, bwd) = le_i32(input)?; // backward ptr
    let (input, free) = le_i32(input)?; // first free address
    Ok((
        input,
        DAFHeader {
            idword: String::from_utf8_lossy(id_word).trim().to_string(),
            internal_filename: String::from_utf8_lossy(ifname).trim().to_string(),
            fward: fwd,
            bward: bwd,
            free: free,
        },
    ))
}

#[derive(Debug)]
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
struct JPLEphem {
    daf_header: DAFHeader,
    header: JPLEphemHeader,
}

impl JPLEphem {
    fn new(jpl_path: &Utf8Path) -> Self {
        let mut file = BufReader::new(File::open(jpl_path).unwrap());

        let mut buffer: Vec<u8> = Vec::new();
        let mut temp = [0u8; 1 << 10];
        file.read_exact(&mut temp)
            .expect("Failed to read the DAF header. (first 1024 bytes)");
        buffer.extend_from_slice(&temp[..]);
        let (_, daf_header) =
            parse_daf_header(&buffer).expect("Failed to parse the DAF header with nom !");

        let end_ascii_comment = ((daf_header.fward as usize) - 1) * 1024;
        let mut ascii_buffer = vec![0u8; end_ascii_comment];
        file.read_exact(&mut ascii_buffer)
            .expect("Failed to read the ASCII comment area.");

        buffer.extend_from_slice(&ascii_buffer[..]);

        let binding = String::from_utf8_lossy(&buffer).replace('\0', "\n");
        let ascii_comment = binding.as_str();

        let (_, jpl_header) =
            parse_ephem_header(&ascii_comment).expect("Failed to parse the JPL header with nom !");
        JPLEphem {
            daf_header: daf_header,
            header: jpl_header,
        }
    }
}

#[cfg(test)]
mod jpl_reader_test {
    use std::io::BufReader;
    use super::*;
    use crate::jpl_ephem::download_jpl_file::get_ephemeris_file;

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

        assert_eq!(jpl_ephem.daf_header.idword, "DAF/SPK");
        assert_eq!(jpl_ephem.header.version, "DE440");
    }
}
