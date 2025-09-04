//! Legacy JPL DE binary ephemerides (“Horizons”) reader.
//!
//! This module provides full support for parsing and interpolating
//! the **legacy binary DE files** (distributed by JPL Horizons under
//! `.../eph/planets/Linux/`).
//!
//! The architecture is split into several submodules, each handling a
//! specific aspect of the format:
//!
//! - [`horizon_data`] — High-level entry point. Loads a DE binary,
//!   parses headers and records, and provides the public API
//!   [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris)
//!   to query state vectors (pos/vel/acc).
//!
//! - [`horizon_ids`] — Identifiers for bodies and centers used in the DE files.
//!   Defines the `HorizonID` enum mapping to official JPL numbering
//!   (Earth, Sun, Moon, EMB, etc.).
//!
//! - [`horizon_records`] — Low-level Chebyshev records. Contains
//!   the coefficient blocks extracted from the file, with methods
//!   to evaluate the polynomial at arbitrary epochs.
//!
//! - [`horizon_version`] — Handles mapping from DE version (e.g., DE440)
//!   to official filenames and metadata. Provides
//!   [`JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion).
//!
//! - [`interpolation_result`] — Unified struct for the output of an
//!   interpolation: position, velocity, acceleration, with conversion
//!   helpers to astronomical units.
//!
//! # Typical workflow
//! 1. Select a
//!    [`JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion)
//!    and resolve the file path (via
//!    \[`crate::jpl_ephem::download_jpl_file::download_big_file`\], gated by the `jpl-download` feature).
//! 2. Construct a
//!    [`HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData)
//!    by reading the binary file.
//! 3. Call
//!    [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris)
//!    with a target, center, epoch, and interpolation options.
//! 4. Use
//!    [`interpolation_result::InterpResult::to_au`](crate::jpl_ephem::horizon::interpolation_result::InterpResult::to_au)
//!    to convert the output into AU / AU/day.
//!
//! # File layout (legacy)
//!
//! The legacy Horizons/DE binary is organized as a header followed by a
//! sequence of **DATA RECORDS** (fixed time-span blocks) containing Chebyshev
//! coefficients for each supported body.
//!
//! ```text
//! +-----------------------------------------------------------------------------------+
//! |                                HORIZONS LEGACY FILE                               |
//! +-----------------------------------------------------------------------------------+
//! | Preamble / constants (optional, format-dependent)                                 |
//! |   (may include ASCII banner, constants table, etc. — often skipped by parser)     |
//! +-----------------------------------------------------------------------------------+
//! |                                HEADER RECORD                                      |
//! +-----------------------------------------------------------------------------------+
//! | NUMDE (u32)                      |  Ephemeris identifier (e.g., 440 for DE440)    |
//! | IPT[15][3] (u32)                 |  Index Pointer Table per body:                 |
//! |   for b = 0..14:                 |    IPT[b] = [ offset, n_coeff, n_sub ]         |
//! |                                  |    - offset  : start position (in 4-byte words |
//! |                                  |                 from beginning of DATA RECORD) |
//! |                                  |    - n_coeff : #Chebyshev coeff. per component |
//! |                                  |    - n_sub   : #sub-intervals per block        |
//! | LPT patch (u32 × 3)              |  Special replacement for IPT[12] (libration)   |
//! | EMRAT (f64)                      |  Earth–Moon mass ratio                         |
//! | START_PERIOD (f64)               |  File coverage start (JD TDB)                  |
//! | END_PERIOD (f64)                 |  File coverage end (JD TDB)                    |
//! | PERIOD_LENGTH (f64)              |  Duration of one DATA RECORD (days)            |
//! | RECSIZE (u32 or computed)        |  Size (bytes) of one DATA RECORD               |
//! +-----------------------------------------------------------------------------------+
//!
//!   Legend for bodies (b = 0..14):
//!     0..10 : Planets + Sun (3 components: x,y,z)
//!     11    : Earth–Moon barycenter (2 components)
//!     12    : Lunar librations (3 components)
//!     13    : Lunar Euler angle rates (3 components)
//!     14    : TT–TDB offset (1 component)
//!
//! +-----------------------------------------------------------------------------------+
//! |                               DATA RECORD  #k                                     |
//! |                         Covers time [T_k, T_{k+1}]                                |
//! +-----------------------------------------------------------------------------------+
//! | T_START (f64)                   | Block start epoch (JD TDB)                      |
//! | T_END   (f64)                   | Block end   epoch (JD TDB)                      |
//! |                                                                                   |
//! |  Body 0 payload (at offset IPT[0][0] * 4 bytes from DATA RECORD start)            |
//! |    For sub = 0 .. IPT[0][2]-1:                                                    |
//! |      For comp = 0 .. dim(0)-1:                                                    |
//! |        Chebyshev coeffs[ IPT[0][1] ]  (f64 × n_coeff)                             |
//! |                                                                                   |
//! |  Body 1 payload (at offset IPT[1][0] * 4)                                         |
//! |    For sub = 0 .. IPT[1][2]-1:                                                    |
//! |      For comp = 0 .. dim(1)-1:                                                    |
//! |        Chebyshev coeffs[ IPT[1][1] ]  (f64 × n_coeff)                             |
//! |                                                                                   |
//! |  ...                                                                              |
//! |                                                                                   |
//! |  Body 14 payload (at offset IPT[14][0] * 4)                                       |
//! |    For sub = 0 .. IPT[14][2]-1:                                                   |
//! |      For comp = 0 .. dim(14)-1:                                                   |
//! |        Chebyshev coeffs[ IPT[14][1] ]  (f64 × n_coeff)                            |
//! +-----------------------------------------------------------------------------------+
//!
//! +-----------------------------------------------------------------------------------+
//! |                               DATA RECORD  #k+1                                   |
//! |                                   (same layout)                                   |
//! +-----------------------------------------------------------------------------------+
//! |                                     ...                                           |
//! +-----------------------------------------------------------------------------------+
//! ```
//!
//! **Notes**
//! 1) Endianness: integers and floats are little-endian in this parser (`le_u32`, `le_f64`).  
//! 2) Offsets in the IPT are expressed in 4‑byte words. Multiply by 4 to get bytes.  
//! 3) Total record size (RECSIZE) can be recomputed as:
//!    ```text
//!    kernel_words = 4
//!                  + Σ_{b=0..14} [ 2 * IPT[b][2] * IPT[b][1] * dim(b) ]
//!    RECSIZE = kernel_words * 4
//!    ```
//!    where `dim(b)` is the component count for body `b` (see legend), and the `2×`
//!    factor accounts for position & velocity sets in the legacy split.
//! 4) Special rows: `b = 12` (librations) uses the LPT patch row after `NUMDE`;
//!    `b = 14` (TT–TDB) is 1D time-offset series.
//! 5) Interpolation: within a DATA RECORD, select the sub-interval by time, then evaluate
//!    the Chebyshev polynomial using the `n_coeff` coefficients for each component.
//!
//! # See also
//! * [`crate::jpl_ephem::naif`] — Alternative backend for NAIF SPK/DAF kernels.
//! * JPL Horizons FTP: <https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/>

/// High-level API: load DE binaries, manage records, query state vectors.
pub mod horizon_data;
/// Body/center identifiers (HorizonID) mapped to official JPL numbering.
pub mod horizon_ids;
/// Low-level Chebyshev coefficient records extracted from the binary file.
pub mod horizon_records;
/// Version management: DE labels to official filenames and metadata.
pub mod horizon_version;
/// Interpolation outputs (position, velocity, acceleration) with unit conversions.
pub mod interpolation_result;
