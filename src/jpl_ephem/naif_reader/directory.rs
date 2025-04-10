use std::{
    fs::File,
    io::{BufReader, Read, Seek},
};

use hifitime::{Duration, Epoch};
use nom::number::complete::le_f64;

/// Struct to hold the directory data of a JPL ephemeris file.
/// The directory data is located at the end of an ephemeris record and contains
/// information about the ephemeris records.
/// The data is stored in a binary format and is read from the file
/// using the `nom` library.
/// The data is stored in the following format:
/// - `init`: Initial epoch in seconds since J2000
/// - `intlen`: Interval length in seconds
/// - `rsize`: Size of the ephemeris record in bytes
/// - `n_records`: Number of ephemeris records
/// The data is read from the file in little-endian format.
#[derive(Debug, PartialEq)]
pub struct DirectoryData {
    pub init: f64,
    pub intlen: usize,
    pub rsize: usize,
    pub n_records: usize,
}

impl DirectoryData {
    /// Parses the directory data from the file.
    /// The data is read from the file in little-endian format.
    ///
    /// Arguments
    /// ---------
    /// `file`: A mutable reference to a `BufReader<File>` object that represents the file to read from.
    /// `end_addr`: The end address of the ephemeris record.
    ///
    /// Returns
    /// -------
    /// A `DirectoryData` object containing the parsed data.
    pub fn parse(file: &mut BufReader<File>, end_addr: usize) -> Self {
        let directory_offset_bytes = (end_addr - 4) * 8;
        let mut dir_buf = [0u8; 32]; // 4 f64 = 32 octets
        file.seek(std::io::SeekFrom::Start(directory_offset_bytes as u64))
            .unwrap();
        file.read_exact(&mut dir_buf).unwrap();
        let (input, init) = le_f64::<_, nom::error::Error<_>>(dir_buf.as_slice()).unwrap();
        let (input, intlen) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        let (input, rsize) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        let (_, n_records) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        DirectoryData {
            init,
            intlen: intlen as usize,
            rsize: rsize as usize,
            n_records: n_records as usize,
        }
    }
}

impl std::fmt::Display for DirectoryData {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let epoch = Epoch::from_et_seconds(self.init);
        let record_length = Duration::from_seconds(self.intlen as f64);

        writeln!(f, "+----------------+----------------------------+")?;
        writeln!(f, "| {:<14} | {:<26} |", "Field", "Value")?;
        writeln!(f, "+----------------+----------------------------+")?;
        writeln!(
            f,
            "| {:<14} | {:<26} |",
            "init (epoch)",
            format!("{}", epoch)
        )?;
        writeln!(
            f,
            "| {:<14} | {:<26} |",
            "intlen",
            format!("{}", record_length)
        )?;
        writeln!(f, "| {:<14} | {:<26} |", "rsize", format!("{}", self.rsize))?;
        writeln!(
            f,
            "| {:<14} | {:<26} |",
            "n_records",
            format!("{}", self.n_records)
        )?;
        writeln!(f, "+----------------+----------------------------+")?;

        Ok(())
    }
}

#[cfg(test)]
mod test_directory {
    use super::*;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_directory_display() {
        use crate::jpl_ephem::download_jpl_file::get_ephemeris_file;

        let version = Some("de440");
        let user_path = None;

        let file_path = get_ephemeris_file(version, user_path).unwrap();
        let mut file = BufReader::new(File::open(file_path).unwrap());
        let end_addr = 2217924;

        let dir_data = DirectoryData::parse(&mut file, end_addr);
        let output = format!("{}", dir_data);
        let expected_output = r#"+----------------+----------------------------+
| Field          | Value                      |
+----------------+----------------------------+
| init (epoch)   | 1549-12-31T00:00:00 ET     |
| intlen         | 8 days                     |
| rsize          | 44                         |
| n_records      | 50224                      |
+----------------+----------------------------+
"#;

        assert_eq!(output, expected_output);
    }
}
