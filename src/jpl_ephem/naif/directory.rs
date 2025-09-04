//! Directory section reader for NAIF/JPL ephemeris records (DAF/SPK).
//!
//! This module extracts the **directory data** stored at the tail of a DAF
//! record (typically an SPK ephemeris record). The directory is a compact
//! four‐value footer that encodes:
//!
//! * `init` — initial epoch of the first record (ET seconds from J2000 TDB),
//! * `intlen` — length of each record’s time span (seconds),
//! * `rsize` — **record size in double-precision words**, not bytes,
//! * `n_records` — total number of records.
//!
//! # Endianness & layout
//! Values are read as little‑endian `f64` (IEEE‑754). A DAF “address” counts
//! **double‑precision words (8‑byte units)** starting at 1. To locate the
//! directory, this module seeks to `(end_addr - 4) * 8` bytes and reads the
//! last four `f64` values of the record.
//!
//! # Example
//! ```rust, ignore
//! use std::fs::File;
//! use std::io::{BufReader};
//! use outfit::jpl_ephem::naif::directory::DirectoryData;
//!
//! let f = File::open("de440.bsp")?;
//! let mut reader = BufReader::new(f);
//!
//! // Suppose `end_addr` is known from the DAF summary (in DP-words, 1-based).
//! let end_addr: usize = /* ... */ 0;
//! let dir = DirectoryData::parse(&mut reader, end_addr);
//! println!("{dir}");
//! assert!(dir.n_records > 0);
//! ```
//!
//! # See also
//! ------------
//! * `daf_header` module — parses the leading DAF header (idword, ND/NI, pointers).
//! * NAIF SPK Required Reading — details on directory/footer semantics.

use std::{
    fs::File,
    io::{BufReader, Read, Seek},
};

use hifitime::{Duration, Epoch};
use nom::number::complete::le_f64;

/// Directory footer of a DAF/SPK ephemeris record.
///
/// The four values mirror the NAIF directory layout and are exposed with
/// convenient Rust types. Note that `rsize` is expressed in **double-precision
/// words** (8‑byte units), not bytes.
///
/// * `init` is an epoch in **ET seconds from J2000 TDB**.
/// * `intlen` is the duration per record, in **seconds**.
/// * `rsize` is the record size (number of `f64` words per record).
/// * `n_records` is the count of records in the file.
///
/// See also
/// ------------
/// * [`DirectoryData::parse`] – Reads the four directory values from a file.
/// * NAIF DAF/SPK docs – Authoritative description of addresses and layout.
#[derive(Debug, PartialEq, Clone)]
pub struct DirectoryData {
    pub init: f64,
    pub intlen: usize,
    pub rsize: usize,
    pub n_records: usize,
}

impl DirectoryData {
    /// Parse the 4‑value directory footer from a DAF/SPK record.
    ///
    /// Arguments
    /// -----------------
    /// * `file`: A mutable [`BufReader<File>`] positioned anywhere; this function will
    ///   seek to the directory location based on `end_addr`.
    /// * `end_addr`: End address of the record in **double‑precision words** (1‑based).
    ///
    /// Return
    /// ----------
    /// * A [`DirectoryData`] populated from the last four `f64` values of the record:
    ///   `(init, intlen, rsize, n_records)`. The floating‑point values are cast to
    ///   `usize` for `intlen`, `rsize`, and `n_records`.
    ///
    /// See also
    /// ------------
    /// * `daf_header` module – To obtain `ND/NI` and summary pointers required to
    ///   compute addresses.
    /// * NAIF SPK Required Reading – Directory/footer interpretation and units.
    pub fn parse(file: &mut BufReader<File>, end_addr: usize) -> Self {
        // Last four f64 values live at addresses end_addr-3 .. end_addr (DP-words).
        let directory_offset_bytes = (end_addr - 4) * 8;
        let mut dir_buf = [0u8; 32]; // 4 f64 = 32 bytes
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
    /// Pretty‑print the directory as a small, fixed‑width table.
    ///
    /// Arguments
    /// -----------------
    /// * `f`: The standard Rust formatter sink.
    ///
    /// Return
    /// ----------
    /// * A [`std::fmt::Result`] indicating success or failure when writing.
    ///
    /// See also
    /// ------------
    /// * [`DirectoryData::parse`] – Produces the values displayed here.
    /// * `hifitime::Epoch`/`Duration` – Formatting of ET epochs and durations.
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
    fn test_directory_display() {
        let dir_data = DirectoryData {
            init: -14200747200.0,
            intlen: 1382400,
            rsize: 41,
            n_records: 25112,
        };

        let output = format!("{dir_data}");
        let expected_output = r#"+----------------+----------------------------+
| Field          | Value                      |
+----------------+----------------------------+
| init (epoch)   | 1549-12-31T00:00:00 ET     |
| intlen         | 16 days                    |
| rsize          | 41                         |
| n_records      | 25112                      |
+----------------+----------------------------+
"#;

        assert_eq!(output, expected_output);
    }
}
