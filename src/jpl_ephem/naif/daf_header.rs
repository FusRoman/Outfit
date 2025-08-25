//! DAF (Double Precision Array File) header parsing utilities.
//!
//! This module provides a lightweight parser and a pretty-printer for the
//! fixed-size DAF/SPK file header used by NAIF/SPICE ephemerides.
//! It reads the first record of a DAF container (typically an SPK kernel),
//! extracts structural metadata (summary layout, directory pointers, address
//! of the first free word, etc.), and exposes them through the [`DAFHeader`]
//! struct.
//!
//! # What the header contains
//!
//! * **`idword`**: Format identifier (e.g. `"DAF/SPK "`), eight ASCII bytes.
//! * **`nd`** / **`ni`**: Number of double-precision / integer components
//!   in each array summary. For SPK this is commonly `nd = 2`, `ni = 6`,
//!   but other values exist for other DAF-based kernels.
//! * **`fward`** / **`bward`**: Record numbers (1-based) of the first and last
//!   summary/descriptor record. They define the doubly-linked directory of
//!   arrays stored in the file.
//! * **`free`**: Address (1-based, in double-precision words) of the first free
//!   location in the file heap where new data could be appended.
//! * **`internal_filename`**: Human-readable kernel name (60 bytes, padded).
//! * **`locfmt`**: Binary platform tag (e.g. `"BIG-IEEE"`, `"LTL-IEEE"`) telling
//!   how numeric data are encoded *inside the file*.
//! * **`fptstr`**: NAIF FTP sentinel string used for transfer integrity checks.
//!
//! # Endianness & safety notes
//!
//! This parser currently reads header integers using **little-endian** (`le_i32`)
//! because most modern SPK you encounter on common platforms ship as `"LTL-IEEE"`.
//! If you need to support `"BIG-IEEE"` files on a little-endian host, you should
//! detect `locfmt` first and dispatch to the appropriate integer/float readers.
//!
//! The header itself is fixed-size and includes reserved bytes. This module
//! skips them as opaque padding.
//!
//! # Example
//!
//! ```rust, no_run
//! use std::fs::File;
//! use std::io::Read;
//! use outfit::jpl_ephem::naif::daf::DAFHeader; // adjust the path to your crate
//!
//! let mut buf = vec![0u8; 1024]; // the first DAF record is 1024 bytes
//! File::open("de440.bsp")?.read_exact(&mut buf)?;
//! let (_rest, header) = DAFHeader::parse(&buf).expect("valid DAF header");
//! println!("{header}");
//! assert_eq!(header.idword, "DAF/SPK");
//! ```
//!
//! # See also
//! ------------
//! * [`DAFHeader::parse`] – Binary decoder for the first DAF record.
//! * [`core::fmt::Display`] for [`DAFHeader`] – Nicely formatted, fixed-width summary.
//! * NAIF/SPICE DAF and SPK required reading (file layout, summaries, addresses).

use nom::{bytes::complete::take, number::complete::le_i32, IResult};

/// In-memory representation of the DAF/SPK header (first 1024‑byte record).
///
/// The fields mirror the canonical NAIF layout and are already trimmed of
/// trailing padding where applicable (e.g., `internal_filename`, `idword`,
/// `locfmt`, `fptstr`).
#[derive(Debug, PartialEq, Clone)]
pub struct DAFHeader {
    /// 8‑byte identifier, typically `"DAF/SPK"`.
    pub idword: String,
    /// 60‑byte, padded internal kernel name.
    pub internal_filename: String,
    /// Number of double-precision components in each summary (ND).
    pub nd: i32,
    /// Number of integer components in each summary (NI).
    pub ni: i32,
    /// Record index of the first summary record (forward pointer).
    pub fward: i32,
    /// Record index of the last summary record (backward pointer).
    pub bward: i32,
    /// First free address (in double-precision words, 1‑based).
    pub free: i32,
    /// Platform tag describing numeric representation (e.g. `"LTL-IEEE"`).
    pub locfmt: String,
    /// NAIF FTP sentinel string.
    pub fptstr: String,
}

impl DAFHeader {
    /// Parse the first 1024‑byte DAF record into a [`DAFHeader`].
    ///
    /// Arguments
    /// -----------------
    /// * `input`: A byte slice starting at the beginning of the file, at least 1024 bytes long.
    ///
    /// Return
    /// ----------
    /// * An [`IResult`] whose value is a tuple `(remaining, header)`. On success,
    ///   `remaining` points to the rest of the input after the header; `header`
    ///   contains all extracted fields with trailing spaces removed.
    ///
    /// See also
    /// ------------
    /// * [`Self::parse`] (this function) – Binary decoder for the DAF header.
    /// * [`core::fmt::Display`] for [`DAFHeader`] – Human‑readable rendering.
    pub fn parse(input: &[u8]) -> IResult<&[u8], Self> {
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
                free,
                locfmt: String::from_utf8_lossy(locfmt).trim().to_string(),
                fptstr: String::from_utf8_lossy(ftpstr).trim().to_string(),
            },
        ))
    }
}

use std::fmt;

impl fmt::Display for DAFHeader {
    /// Render a fixed-width table summarizing the DAF header fields.
    ///
    /// Arguments
    /// -----------------
    /// * `f`: Standard formatter provided by the Rust formatting machinery.
    ///
    /// Return
    /// ----------
    /// * A [`fmt::Result`] indicating whether writing to the formatter succeeded.
    ///
    /// See also
    /// ------------
    /// * [`DAFHeader::parse`] – Produces the data displayed here.
    /// * [`core::fmt::Display`] – Rust’s formatting traits and conventions.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const LABEL_WIDTH: usize = 18;
        const VALUE_WIDTH: usize = 50;

        let border = format!(
            "+{:-<label$}+{:-<value$}+",
            "",
            "",
            label = LABEL_WIDTH + 1,
            value = VALUE_WIDTH + 1
        );

        writeln!(f, "{border}")?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "DAF File Header",
            "",
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(f, "{border}")?;

        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "ID Word",
            format!("{} (Format ID)", self.idword),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Internal Name",
            format!("{}", self.internal_filename),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "ND (doubles)",
            format!("{} double precision summary components", self.nd),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "NI (integers)",
            format!("{} integer summary components", self.ni),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Forward Ptr",
            format!("Record # of first summary: {}", self.fward),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Backward Ptr",
            format!("Record # of last summary: {}", self.bward),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Free Addr",
            format!("Next free address: {}", self.free),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Binary Format",
            format!("{} (e.g., BIG-IEEE or LTL-IEEE)", self.locfmt),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;

        writeln!(f, "{border}")
    }
}

#[cfg(test)]
mod test_daf_header {
    use super::*;

    #[test]
    fn test_display_daf_header() {
        let header = DAFHeader {
            idword: "DAF/SPK".to_string(),
            internal_filename: "NIO2SPK".to_string(),
            nd: 2,
            ni: 6,
            fward: 62,
            bward: 62,
            free: 14974889,
            locfmt: "LTL-IEEE".to_string(),
            fptstr: "FTPSTR:\r:\n:\r\n:\r\0:�:\u{10}�:ENDFTP".to_string(),
        };

        let expected = r#"+-------------------+---------------------------------------------------+
| DAF File Header   |                                                   |
+-------------------+---------------------------------------------------+
| ID Word           | DAF/SPK (Format ID)                               |
| Internal Name     | NIO2SPK                                           |
| ND (doubles)      | 2 double precision summary components             |
| NI (integers)     | 6 integer summary components                      |
| Forward Ptr       | Record # of first summary: 62                     |
| Backward Ptr      | Record # of last summary: 62                      |
| Free Addr         | Next free address: 14974889                       |
| Binary Format     | LTL-IEEE (e.g., BIG-IEEE or LTL-IEEE)             |
+-------------------+---------------------------------------------------+
"#;
        let output = format!("{header}");
        assert_eq!(output, expected);
    }
}
