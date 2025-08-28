//! # Outfit error types
//!
//! This module defines the [`OutfitError`](crate::outfit_errors::OutfitError) enum, the unified error type
//! for the **Outfit** crate.  
//!
//! # Outfit error taxonomy (thematic groups)
//!
//! To make diagnostics and handling more systematic, we group variants into themes below.
//!
//! ## 1) Ephemerides & SPK handling
//!
//! **When**: Selecting the JPL source/version, opening/reading files, decoding SPK content.
//!
//! - [`OutfitError::InvalidJPLStringFormat`](crate::outfit_errors::OutfitError::InvalidJPLStringFormat)
//! - [`OutfitError::InvalidJPLEphemFileSource`](crate::outfit_errors::OutfitError::InvalidJPLEphemFileSource)
//! - [`OutfitError::InvalidJPLEphemFileVersion`](crate::outfit_errors::OutfitError::InvalidJPLEphemFileVersion)
//! - [`OutfitError::JPLFileNotFound`](crate::outfit_errors::OutfitError::JPLFileNotFound)
//! - [`OutfitError::InvalidSpkDataType`](crate::outfit_errors::OutfitError::InvalidSpkDataType)
//!
//! **Typical causes**: malformed `"horizon:DE###"` strings; unsupported ephemeris revision;
//! missing local file; unexpected SPK segment type.
//!
//! **Remediation**: validate the source string early; log available ephemeris revisions;
//! verify file paths; assert known SPK types during parsing.
//!
//! ## 2) Network & filesystem I/O
//!
//! **When**: Downloading MPC/JPL resources, touching files/directories, handling paths.
//!
//! - [`OutfitError::InvalidUrl`](crate::outfit_errors::OutfitError::InvalidUrl)
//! - [`OutfitError::UreqHttpError`](crate::outfit_errors::OutfitError::UreqHttpError)
//! - [`OutfitError::ReqwestError`](crate::outfit_errors::OutfitError::ReqwestError) *(feature: `jpl-download`)*
//! - [`OutfitError::IoError`](crate::outfit_errors::OutfitError::IoError)
//! - [`OutfitError::UnableToCreateBaseDir`](crate::outfit_errors::OutfitError::UnableToCreateBaseDir)
//! - [`OutfitError::Utf8PathError`](crate::outfit_errors::OutfitError::Utf8PathError)
//!
//! **Typical causes**: malformed URL; HTTP errors (4xx/5xx, timeouts); missing permissions;
//! non-UTF-8 paths on some platforms.
//!
//! **Remediation**: sanitize/format URLs; add retries/timeouts; create directories recursively;
//! normalize/encode paths.
//!
//! ## 3) Parsing & data ingestion
//!
//! **When**: Parsing observations, fixed-width legacy records, or internal byte streams.
//!
//! - [`OutfitError::NomParsingError`](crate::outfit_errors::OutfitError::NomParsingError)
//! - [`OutfitError::Parsing80ColumnFileError`](crate::outfit_errors::OutfitError::Parsing80ColumnFileError)
//!
//! **Typical causes**: schema drift; corrupted inputs; unexpected whitespace/locale separators.
//!
//! **Remediation**: strengthen parsers with tolerant combinators; surface line/column context;
//! round-trip tests with golden files.
//!
//! ## 4) Numerical methods & stochastic routines
//!
//! **When**: Root-finding, complex polynomial solvers, Gaussian noise injection.
//!
//! - [`OutfitError::RootFindingError`](crate::outfit_errors::OutfitError::RootFindingError)
//! - [`OutfitError::PolynomialRootFindingFailed`](crate::outfit_errors::OutfitError::PolynomialRootFindingFailed)
//! - [`OutfitError::NoiseInjectionError`](crate::outfit_errors::OutfitError::NoiseInjectionError)
//!
//! **Typical causes**: poor initial guesses; ill-conditioned polynomials; invalid σ parameters.
//!
//! **Remediation**: guard preconditions and parameter ranges; implement fallbacks
//! or multiple initial seeds; log residuals and iteration counts.
//!
//! ## 5) Reference frames & orbit determination
//!
//! **When**: Building frame rotations, IOD preconditions, sanity checks on states.
//!
//! - [`OutfitError::InvalidIODParameter`](crate::outfit_errors::OutfitError::InvalidIODParameter)
//! - [`OutfitError::InvalidRefSystem`](crate::outfit_errors::OutfitError::InvalidRefSystem)
//! - [`OutfitError::VelocityCorrectionError`](crate::outfit_errors::OutfitError::VelocityCorrectionError)
//! - [`OutfitError::InvalidOrbit`](crate::outfit_errors::OutfitError::InvalidOrbit)
//! - [`OutfitError::SingularDirectionMatrix`](crate::outfit_errors::OutfitError::SingularDirectionMatrix)
//! - [`OutfitError::SpuriousRootDetected`](crate::outfit_errors::OutfitError::SpuriousRootDetected)
//! - [`OutfitError::GaussNoRootsFound`](crate::outfit_errors::OutfitError::GaussNoRootsFound)
//!
//! **Typical causes**: coplanar/degenerate geometry; negative/near-zero ranges;
//! inconsistent elements; unsupported frame transitions.
//!
//! **Remediation**: pre-filter observations (geometry checks); enforce bounds;
//! provide multiple IOD seeds; degrade gracefully to alternative solvers.
//!
//! ## 6) Observation catalog & indexing
//!
//! **When**: Looking up observations by index or key in internal stores.
//!
//! - [`OutfitError::ObservationNotFound`](crate::outfit_errors::OutfitError::ObservationNotFound)
//!
//! **Typical causes**: stale indices; filtered/compacted buffers.
//!
//! **Remediation**: validate indices; prefer keyed lookups; add debug metadata in errors.
//!
//! ---
//!
//! ## Testing & equality
//!
//! - [`OutfitError`](crate::outfit_errors::OutfitError) implements [`PartialEq`] for *deterministic* variants; variants that
//!   wrap opaque external errors (e.g., I/O/HTTP) compare equal by *kind* only.
//! - Use pattern matches in tests to assert the right **class** of failure without relying
//!   on OS-/locale-dependent strings.
//!
//! ```rust,ignore
//! match result {
//!     Err(OutfitError::InvalidRefSystem(msg)) => assert!(msg.contains("Equm")),
//!     Err(e) => panic!("unexpected error: {e}"),
//!     Ok(_) => {}
//! }
//! ```
//!
//! `OutfitError` implements [`thiserror::Error`] for user-friendly messages,
//! as well as [`PartialEq`] for convenient testing and assertions.
//!
//! ## See also
//! -------------
//! * [`Outfit`](crate::outfit::Outfit) – main context using this error type.
//! * [`ParseObsError`](crate::observations::ParseObsError) – specific parser error re-used here.
//! * [`roots::SearchError`] – root-finding errors wrapped into this type.
//! * [`rand_distr::NormalError`] – Gaussian noise generation errors.
//!
//! ## Example
//! ```rust,ignore
//! use outfit::outfit_errors::OutfitError;
//!
//! fn do_something() -> Result<(), OutfitError> {
//!     Err(OutfitError::InvalidUrl("http:/bad_url".into()))
//! }
//!
//! match do_something() {
//!     Err(OutfitError::InvalidUrl(msg)) => eprintln!("Bad URL: {msg}"),
//!     Err(e) => eprintln!("General error: {e}"),
//!     Ok(_) => println!("Success"),
//! }
//! ```

use thiserror::Error;

use crate::observations::ParseObsError;

/// All possible errors raised within the Outfit crate.
///
/// Each variant corresponds to a specific failure mode, either external
/// (I/O, HTTP, parsing) or internal (root-finding, invalid parameters).
#[derive(Error, Debug)]
pub enum OutfitError {
    /// Invalid JPL string format (e.g., `"horizonDE440"` instead of `"horizon:DE440"`).
    #[error("Invalid JPL string format: {0}")]
    InvalidJPLStringFormat(String),

    /// Invalid JPL ephemeris file source (bad scheme or unknown provider).
    #[error("Invalid JPL ephemeris file source: {0}")]
    InvalidJPLEphemFileSource(String),

    /// Unsupported or unrecognized JPL ephemeris file version.
    #[error("Invalid JPL ephemeris file version: {0}")]
    InvalidJPLEphemFileVersion(String),

    /// Invalid URL string (parse failure).
    #[error("Invalid URL: {0}")]
    InvalidUrl(String),

    /// HTTP error raised by [`ureq`] client.
    #[error("HTTP ureq error: {0}")]
    UreqHttpError(#[from] ureq::Error),

    /// Filesystem I/O error (open/read/write).
    #[error("Unable to perform file operation: {0}")]
    IoError(#[from] std::io::Error),

    /// HTTP error raised by [`reqwest`] client (only with `jpl-download` feature).
    #[cfg(feature = "jpl-download")]
    #[error("HTTP reqwest error: {0}")]
    ReqwestError(#[from] reqwest::Error),

    /// Unable to create the base directory for JPL ephemeris file storage.
    #[error("Base dir creation error for JPL ephemeris file: {0}")]
    UnableToCreateBaseDir(String),

    /// UTF-8 error when decoding a filesystem path.
    #[error("UTF-8 Path error: {0}")]
    Utf8PathError(String),

    /// JPL ephemeris file not found at the expected location.
    #[error("JPL File not found at: {0}")]
    JPLFileNotFound(String),

    /// Failure during polynomial or numerical root finding (from `roots` crate).
    #[error("ROOTS finding error: {0}")]
    RootFindingError(#[from] roots::SearchError),

    /// Observation not found (index out of range or missing data).
    #[error("Observation not found: {0}")]
    ObservationNotFound(usize),

    /// Invalid error model specification (e.g., unknown identifier).
    #[error("Invalid Error model: {0}")]
    InvalidErrorModel(String),

    /// Invalid file path provided for an error model file.
    #[error("Invalid File Path for Error model: {0}")]
    InvalidErrorModelFilePath(String),

    /// Error during parsing with `nom` parser combinators.
    #[error("Error during the nom parsing: {0}")]
    NomParsingError(String),

    /// Error while parsing legacy 80-column formatted observation files.
    #[error("Error during the 80 column file parsing: {0}")]
    Parsing80ColumnFileError(ParseObsError),

    /// Gaussian noise generation failed (invalid parameters for normal distribution).
    #[error("Gaussian noise generation failed: {0:?}")]
    NoiseInjectionError(rand_distr::NormalError),

    /// The unit direction matrix is singular (cannot invert; likely coplanar observations).
    #[error(
        "Unit direction matrix is singular (cannot be inverted); observations may be coplanar"
    )]
    SingularDirectionMatrix,

    /// Aberth–Ehrlich method failed to find acceptable complex polynomial roots.
    #[error("Aberth–Ehrlich method failed to find acceptable complex roots")]
    PolynomialRootFindingFailed,

    /// Spurious root detected (e.g., negative or near-zero geocentric distance).
    #[error("Spurious root detected (e.g., negative or near-zero geocentric distance)")]
    SpuriousRootDetected,

    /// Gauss method for initial orbit determination failed to find roots.
    #[error("Gauss method failed to find roots")]
    GaussNoRootsFound,

    /// Invalid SPK segment data type in JPL ephemeris file.
    #[error("Invalid SPK data type: {0}")]
    InvalidSpkDataType(i32),

    /// Invalid input parameter for initial orbit determination.
    #[error("Invalid IOD parameter: {0}")]
    InvalidIODParameter(String),

    /// Invalid reference system or unsupported transformation request.
    #[error("Invalid reference system: {0}")]
    InvalidRefSystem(String),

    /// Failure during velocity correction procedure.
    #[error("Velocity correction error: {0}")]
    VelocityCorrectionError(String),

    /// Invalid or inconsistent orbital elements/state vector.
    #[error("Invalid orbit: {0}")]
    InvalidOrbit(String),

    /// Generic invalid input error with a descriptive message.
    #[error("Invalid input: {0}")]
    InvalidConversion(String),

    /// Error indicating a floating-point value is NaN (Not a Number).
    #[error("Invalid floating-point value: {0}")]
    InvalidFloatValue(ordered_float::FloatIsNan),
}

impl From<rand_distr::NormalError> for OutfitError {
    fn from(err: rand_distr::NormalError) -> Self {
        OutfitError::NoiseInjectionError(err)
    }
}

impl From<ordered_float::FloatIsNan> for OutfitError {
    fn from(err: ordered_float::FloatIsNan) -> Self {
        OutfitError::InvalidFloatValue(err)
    }
}

impl PartialEq for OutfitError {
    fn eq(&self, other: &Self) -> bool {
        use OutfitError::*;
        match (self, other) {
            (InvalidJPLStringFormat(a), InvalidJPLStringFormat(b)) => a == b,
            (InvalidJPLEphemFileSource(a), InvalidJPLEphemFileSource(b)) => a == b,
            (InvalidJPLEphemFileVersion(a), InvalidJPLEphemFileVersion(b)) => a == b,
            (InvalidUrl(a), InvalidUrl(b)) => a == b,

            // Error variants that cannot be compared directly
            (UreqHttpError(_), UreqHttpError(_)) => true,
            (IoError(_), IoError(_)) => true,
            #[cfg(feature = "jpl-download")]
            (ReqwestError(_), ReqwestError(_)) => true,

            (UnableToCreateBaseDir(a), UnableToCreateBaseDir(b)) => a == b,
            (Utf8PathError(a), Utf8PathError(b)) => a == b,
            (JPLFileNotFound(a), JPLFileNotFound(b)) => a == b,
            (RootFindingError(a), RootFindingError(b)) => a == b,
            (ObservationNotFound(a), ObservationNotFound(b)) => a == b,
            (InvalidErrorModel(a), InvalidErrorModel(b)) => a == b,
            (InvalidErrorModelFilePath(a), InvalidErrorModelFilePath(b)) => a == b,
            (NomParsingError(a), NomParsingError(b)) => a == b,
            (Parsing80ColumnFileError(a), Parsing80ColumnFileError(b)) => a == b,
            (NoiseInjectionError(a), NoiseInjectionError(b)) => a == b,
            (InvalidSpkDataType(a), InvalidSpkDataType(b)) => a == b,
            (InvalidIODParameter(a), InvalidIODParameter(b)) => a == b,
            (InvalidRefSystem(a), InvalidRefSystem(b)) => a == b,
            (VelocityCorrectionError(a), VelocityCorrectionError(b)) => a == b,
            (InvalidOrbit(a), InvalidOrbit(b)) => a == b,
            (InvalidConversion(a), InvalidConversion(b)) => a == b,

            // Variantes unitaires
            (SingularDirectionMatrix, SingularDirectionMatrix) => true,
            (PolynomialRootFindingFailed, PolynomialRootFindingFailed) => true,
            (SpuriousRootDetected, SpuriousRootDetected) => true,
            (GaussNoRootsFound, GaussNoRootsFound) => true,

            _ => false,
        }
    }
}
