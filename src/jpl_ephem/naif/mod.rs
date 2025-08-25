//! NAIF/JPL SPK ephemeris support (binary DAF container + high-level access).
//!
//! This module is the entry point for reading NAIF/SPK kernels and exposing a
//! clean API to query ephemerides. Internally it follows the canonical DAF/SPK
//! layout: a binary DAF header, an ASCII human header, segment summaries,
//! per‑segment directory/footer, and per‑record Chebyshev coefficients.
//!
//! ### Units & time scales
//! * Epochs are in **ET/TDB seconds from J2000** (SPK convention),
//! * Positions are **kilometers** and velocities are **km/s**,
//! * DAF “addresses” are in **double‑precision words** (8‑byte units, 1‑based).
//!
//! ### Endianness
//! Most public NAIF kernels ship as `"LTL-IEEE"`; integers are read little‑endian.
//! If you need `"BIG-IEEE"`, dispatch the integer/float readers accordingly.
//!
//! # File layout (NAIF DAF/SPK)
//!
//! A SPK file is a DAF container made of **fixed 1024-byte records**.
//! *DAF addresses* are in **double-precision words** (DP-words = 8 bytes), 1-based.
//! One 1024-byte record contains 128 DP-words. Address 1 corresponds to the
//! **first word** of **Record #1** (the File Record).
//!
//! ```text
//! +-------------------------------------------------------------------------------------------------+
//! |                                       DAF / SPK KERNEL                                          |
//! +-------------------------------------------------------------------------------------------------+
//! |                                          RECORD #1                                              |
//! |                                        (FILE RECORD)                                            |
//! +-------------------------------------------------------------------------------------------------+
//! | IDWORD[8]                 | "DAF/SPK "                                                          |
//! | ND(i32), NI(i32)          | #doubles / #integers per summary                                    |
//! | IFNAME[60]                | Internal file name                                                  |
//! | FWARD(i32), BWARD(i32)    | Record numbers (1-based) of first/last Summary Records              |
//! | FREE(i32)                 | First free address (in DP-words)                                    |
//! | LOCFMT[8]                 | "LTL-IEEE" or "BIG-IEEE"                                            |
//! | RESERVED[603]             | Padding                                                             |
//! | FTPSTR[28]                | NAIF FTP sentinel                                                   |
//! +-------------------------------------------------------------------------------------------------+
//! |                         RECORDS #2 .. (FWARD-1) : COMMENT RECORDS                               |
//! |  ASCII area (includes the JPL “ephemeris header”: version, dates, coverage calendar + JD, etc.) |
//! +-------------------------------------------------------------------------------------------------+
//! |                                     RECORD #(FWARD) : SUMMARY RECORD                            |
//! |     (and possibly subsequent ones up to BWARD if the table does not fit in a single record)     |
//! +-------------------------------------------------------------------------------------------------+
//! | Control area (at the beginning of the Summary Record)                                           |
//! |   NEXT(i32), PREV(i32), NSUM(i32 or f64 depending on impl.) *                                   |
//! |   * Some parsers read NSUM as f64 for DP-word alignment                                         |
//! |                                                                                                 |
//! | Array of SUMMARIES (fixed size per entry = SS DP-words)                                         |
//! |   SS = ND + ceil(NI/2)  (integers are packed two per DP-word)                                   |
//! |                                                                                                 |
//! |  For SPK (ND=2, NI=6), one summary encodes:                                                     |
//! |    start_epoch     (f64, ET seconds from J2000)                                                 |
//! |    end_epoch       (f64, ET seconds from J2000)                                                 |
//! |    target          (i32, NAIF ID)                                                               |
//! |    center          (i32, NAIF ID)                                                               |
//! |    frame_id        (i32, e.g. J2000 = 1)                                                        |
//! |    data_type       (i32, e.g. Type 2 = Chebyshev position only)                                 |
//! |    initial_addr    (i32, start address of segment, DP-words, 1-based)                           |
//! |    final_addr      (i32, end   address of segment, DP-words, 1-based)                           |
//! +-------------------------------------------------------------------------------------------------+
//! |                                     NAME RECORD(S) (optional)                                   |
//! |  Each summary may have an associated segment name.                                              |
//! +-------------------------------------------------------------------------------------------------+
//! |                                   ELEMENT (DATA) RECORDS                                        |
//! |  Each **SPK segment** (pointed by a summary) occupies an address range [initial..final].        |
//! |  In a **Type 2** segment (Chebyshev position only), the data are a sequence of                  |
//! |  **Ephemeris Records** of fixed size `rsize` (in DP-words).                                     |
//! |                                                                                                 |
//! |  EPHEMERIS RECORD (Type 2)                                                                      |
//! |    mid     (f64)  : midpoint ET (s)                                                             |
//! |    radius  (f64)  : half interval length (s)                                                    |
//! |    X coeffs (f64 × ncoeff)  -> position X (km)                                                  |
//! |    Y coeffs (f64 × ncoeff)  -> position Y (km)                                                  |
//! |    Z coeffs (f64 × ncoeff)  -> position Z (km)                                                  |
//! |                                                                                                 |
//! |  DIRECTORY FOOTER (4 f64) – **always stored at the very end of the segment**                    |
//! |    init      (f64) : epoch of first record (ET s)                                               |
//! |    intlen    (f64) : interval length per record (s)                                             |
//! |    rsize     (f64) : record size in DP-words (8 B * rsize = bytes per record)                   |
//! |    n_records (f64) : number of records in the segment                                           |
//! +-------------------------------------------------------------------------------------------------+
//!
//! Addressing and record geometry
//! -----------------------------
//! • 1 DP-word = 8 bytes ; 1 record = 1024 bytes = 128 DP-words.  
//! • Record #k (1-based) covers addresses ((k-1)*128 + 1) .. (k*128).  
//! • Byte offset for DAF address `A`:  `offset_bytes = (A-1) * 8`.  
//!
//! Interpolation (Type 2 segments)
//! -------------------------------
//! • Normalized time  :  `t = clamp( (et - mid) / radius , -1, 1 )`  
//! • Position [km]    :  Σ c_n · T_n(t) for each axis (X,Y,Z)  
//! • Velocity [km/s]  :  Σ c_n · T'_n(t) · (2 / radius)  
//! ```
//!
//! ## Notes
//! * Most public kernels are `"LTL-IEEE"` (little-endian).  
//! * Integers in summaries are packed into DP-words → `SS = ND + ceil(NI/2)`.  
//! * Segment *Name Records* exist but are not required for interpolation.
//!
//! The high‑level API you are expected to use is provided by:
//! * [`naif_data`]  – load a BSP and query ephemerides,
//! * [`naif_ids`]   – typed NAIF identifiers (targets/centers, SPK types),
//! * [`naif_version`] – enum of known JPL kernel versions.

mod daf_header; // DAF header reader (binary, private)
mod directory; // Directory/footer per segment (private)
mod ephemeris_record; // Chebyshev records + interpolation (private)
mod jpl_ephem_header; // ASCII JPL header parser (private)
mod summary_record; // Segment summary record (private)

pub mod naif_data; // High-level loader (public API)
pub mod naif_ids; // NAIF identifiers & types (public API)
pub mod naif_version; // Known JPL ephemeris versions (public API)

/// Print a hex + ASCII dump of a byte slice (debug/inspection).
///
/// Arguments
/// -----------------
/// * `data`: Byte slice to dump.
///
/// Return
/// ----------
/// * `()` – the function writes to `stdout`.
///
/// See also
/// ------------
/// * `naif_data` – for reading raw blocks before decoding (useful when debugging).
/// * `daf_header`/`summary_record` – to correlate offsets with decoded structures.
pub fn print_hex_dump(data: &[u8]) {
    const BYTES_PER_LINE: usize = 40;

    for (i, chunk) in data.chunks(BYTES_PER_LINE).enumerate() {
        // Offset (left column)
        print!("{:08x}: ", i * BYTES_PER_LINE);

        // Hex representation
        for byte in chunk {
            print!("{byte:02x} ");
        }

        // Padding if line is shorter than BYTES_PER_LINE
        if chunk.len() < BYTES_PER_LINE {
            for _ in 0..(BYTES_PER_LINE - chunk.len()) {
                print!("   ");
            }
        }

        // ASCII representation (right column)
        print!("|");
        for byte in chunk {
            let ascii = if byte.is_ascii_graphic() || *byte == b' ' {
                *byte as char
            } else {
                '.'
            };
            print!("{ascii}");
        }
        println!("|");
    }
}
