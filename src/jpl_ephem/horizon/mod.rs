//! Legacy JPL DE binary ephemerides (“Horizons”) reader.
//!
//! This module provides the full support for parsing and interpolating
//! the **legacy binary DE files** (distributed by JPL Horizons under
//! `.../eph/planets/Linux/`).  
//!
//! The architecture is split into several submodules, each handling a
//! specific aspect of the format:
//!
//! - [`horizon_data`] — High-level entry point. Loads a DE binary,
//!   parses headers and records, and provides the public API
//!   [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris) to query state vectors (pos/vel/acc).
//!
//! - [`horizon_ids`] — Identifiers for bodies and centers used in the DE files.
//!   Defines the `HorizonID` enum mapping to official JPL numbering
//!   (Earth, Sun, Moon, EMB, etc.).
//!
//! - [`horizon_records`] — Low-level Chebyshev records. Contains
//!   the coefficient blocks extracted from the file, with methods
//!   to evaluate the polynomial at arbitrary epochs.
//!
//! - [`horizon_version`] — Handles mapping from DE version (e.g. DE440)
//!   to official filenames and metadata. Provides [`JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion).
//!
//! - [`interpolation_result`] — Unified struct for the output of an
//!   interpolation: position, velocity, acceleration, with conversion
//!   helpers to astronomical units.
//!
//! # Typical workflow
//! 1. Select a [`JPLHorizonVersion`](crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion) and resolve the file path (via
//!    [`crate::jpl::download_jpl_file`](crate::jpl_ephem::download_jpl_file::download_big_file)).
//! 2. Construct a [`HorizonData`](crate::jpl_ephem::horizon::horizon_data::HorizonData) by reading the binary file.
//! 3. Call [`HorizonData::ephemeris`](crate::jpl_ephem::horizon::horizon_data::HorizonData::ephemeris) with a target, center, epoch,
//!    and interpolation options.
//! 4. Use [`interpolation_result::InterpResult::to_au`](crate::jpl_ephem::horizon::interpolation_result::InterpResult::to_au) to convert
//!    the output into AU / AU/day.
//!
//! # See also
//! * [`crate::jpl::naif`](crate::jpl_ephem::naif) — Alternative backend for NAIF SPK/DAF kernels.
//! * JPL Horizons FTP: <https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/>

pub mod horizon_data;
pub mod horizon_ids;
pub mod horizon_records;
pub mod horizon_version;
pub mod interpolation_result;
