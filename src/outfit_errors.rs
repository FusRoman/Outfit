//! # Outfit error types
//!
//! This module defines the [`OutfitError`](crate::outfit_errors::OutfitError) enum,
//! the **unified error type** for the Outfit crate.  
//! Each variant represents a distinct failure mode, either from external dependencies
//! (I/O, HTTP, parsing) or internal algorithms (orbit determination, root finding).
//!
//! # Outfit error taxonomy (thematic groups)
//!
//! To help diagnostics and error handling, we classify variants into thematic groups.
//!
//! ## 1) Ephemerides & SPK handling
//!
//! **When**: Selecting JPL ephemeris sources, reading SPK files, decoding binary records.
//!
//! - [`InvalidJPLStringFormat`](crate::outfit_errors::OutfitError::InvalidJPLStringFormat)
//! - [`InvalidJPLEphemFileSource`](crate::outfit_errors::OutfitError::InvalidJPLEphemFileSource)
//! - [`InvalidJPLEphemFileVersion`](crate::outfit_errors::OutfitError::InvalidJPLEphemFileVersion)
//! - [`JPLFileNotFound`](crate::outfit_errors::OutfitError::JPLFileNotFound)
//! - [`InvalidSpkDataType`](crate::outfit_errors::OutfitError::InvalidSpkDataType)
//!
//! **Typical causes**: malformed `"horizon:DE###"` strings; unsupported ephemeris revisions;
//! missing local SPK files; unexpected segment types.
//!
//! **Remediation**: validate early; log available revisions; check paths; assert known SPK types.
//!
//! ## 2) Network & filesystem I/O
//!
//! **When**: Downloading MPC/JPL data, touching files or directories.
//!
//! - [`InvalidUrl`](crate::outfit_errors::OutfitError::InvalidUrl)
//! - [`UreqHttpError`](crate::outfit_errors::OutfitError::UreqHttpError)
//! - [`ReqwestError`](crate::outfit_errors::OutfitError::ReqwestError) *(feature: `jpl-download`)*
//! - [`IoError`](crate::outfit_errors::OutfitError::IoError)
//! - [`UnableToCreateBaseDir`](crate::outfit_errors::OutfitError::UnableToCreateBaseDir)
//! - [`Utf8PathError`](crate::outfit_errors::OutfitError::Utf8PathError)
//!
//! **Typical causes**: malformed URL; HTTP 4xx/5xx; permission denied; non-UTF-8 paths.
//!
//! **Remediation**: sanitize URLs; retry with backoff; create directories recursively;
//! normalize/encode paths.
//!
//! ## 3) Parsing & data ingestion
//!
//! **When**: Parsing observations, fixed-width records, or internal streams.
//!
//! - [`NomParsingError`](crate::outfit_errors::OutfitError::NomParsingError)
//! - [`Parsing80ColumnFileError`](crate::outfit_errors::OutfitError::Parsing80ColumnFileError)
//!
//! **Typical causes**: schema drift; corrupted inputs; locale-specific formats.
//!
//! **Remediation**: strengthen parsers; surface line/column context; round-trip with golden files.
//!
//! ## 4) Numerical methods & stochastic routines
//!
//! **When**: Root finding, polynomial solving, Gaussian noise injection.
//!
//! - [`RootFindingError`](crate::outfit_errors::OutfitError::RootFindingError)
//! - [`PolynomialRootFindingFailed`](crate::outfit_errors::OutfitError::PolynomialRootFindingFailed)
//! - [`NoiseInjectionError`](crate::outfit_errors::OutfitError::NoiseInjectionError)
//! - [`InvalidFloatValue`](crate::outfit_errors::OutfitError::InvalidFloatValue)
//!
//! **Typical causes**: poor initial guesses; ill-conditioned polynomials; invalid σ; NaN propagation.
//!
//! **Remediation**: guard parameter ranges; multiple seeds; log residuals; check NaN early.
//!
//! ## 5) Reference frames & orbit determination
//!
//! **When**: Building rotation matrices, running IOD, checking orbital states.
//!
//! - [`InvalidIODParameter`](crate::outfit_errors::OutfitError::InvalidIODParameter)
//! - [`InvalidRefSystem`](crate::outfit_errors::OutfitError::InvalidRefSystem)
//! - [`VelocityCorrectionError`](crate::outfit_errors::OutfitError::VelocityCorrectionError)
//! - [`InvalidOrbit`](crate::outfit_errors::OutfitError::InvalidOrbit)
//! - [`SingularDirectionMatrix`](crate::outfit_errors::OutfitError::SingularDirectionMatrix)
//! - [`SpuriousRootDetected`](crate::outfit_errors::OutfitError::SpuriousRootDetected)
//! - [`GaussNoRootsFound`](crate::outfit_errors::OutfitError::GaussNoRootsFound)
//! - [`InvalidConversion`](crate::outfit_errors::OutfitError::InvalidConversion)
//!
//! **Typical causes**: coplanar geometry; invalid orbital elements; unsupported frame conversions.
//!
//! **Remediation**: pre-filter observations; enforce numeric bounds; fall back to alternative solvers.
//!
//! ## 6) Observation catalog & indexing
//!
//! **When**: Looking up observations in internal stores.
//!
//! - [`ObservationNotFound`](crate::outfit_errors::OutfitError::ObservationNotFound)
//!
//! **Typical causes**: stale indices; filtered/compacted buffers.
//!
//! **Remediation**: validate indices; prefer keyed lookups; embed metadata in errors.
//!
//! ## 7) Error model configuration
//!
//! **When**: Loading or selecting observational error models.
//!
//! - [`InvalidErrorModel`](crate::outfit_errors::OutfitError::InvalidErrorModel)
//! - [`InvalidErrorModelFilePath`](crate::outfit_errors::OutfitError::InvalidErrorModelFilePath)
//!
//! **Typical causes**: wrong model identifier; bad configuration file path.
//!
//! **Remediation**: validate identifiers against registry; check file paths upfront.
//!
//! ---
//!
//! ## Testing & equality
//!
//! - [`OutfitError`](crate::outfit_errors::OutfitError) implements [`PartialEq`] for deterministic variants.
//! - Variants wrapping opaque errors (e.g., I/O, HTTP) compare equal by kind only.
//!
//! ```rust,ignore
//! match result {
//!     Err(OutfitError::InvalidRefSystem(msg)) => assert!(msg.contains("Equm")),
//!     Err(e) => panic!("unexpected error: {e}"),
//!     Ok(_) => {}
//! }
//! ```
//!
//! ## See also
//! -------------
//! * [`Outfit`](crate::outfit::Outfit) – main crate context using this error type.
//! * [`ParseObsError`](crate::observations::ParseObsError) – observation parsing errors.
//! * [`roots::SearchError`] – wrapped in [`RootFindingError`](crate::outfit_errors::OutfitError::RootFindingError).
//! * [`rand_distr::NormalError`] – wrapped in [`NoiseInjectionError`](crate::outfit_errors::OutfitError::NoiseInjectionError).
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
    #[error("Invalid URL string: {0}")]
    InvalidUrl(String),

    /// HTTP error raised by [`ureq`] client.
    #[error("HTTP request failed (ureq): {0}")]
    UreqHttpError(#[from] ureq::Error),

    /// Filesystem I/O error (open/read/write).
    #[error("Filesystem I/O error: {0}")]
    IoError(#[from] std::io::Error),

    /// HTTP error raised by [`reqwest`] client (only with `jpl-download` feature).
    #[cfg(feature = "jpl-download")]
    #[error("HTTP request failed (reqwest): {0}")]
    ReqwestError(#[from] reqwest::Error),

    /// Unable to create the base directory for JPL ephemeris file storage.
    #[error("Failed to create base directory for JPL ephemeris file: {0}")]
    UnableToCreateBaseDir(String),

    /// UTF-8 error when decoding a filesystem path.
    #[error("Filesystem path is not valid UTF-8: {0}")]
    Utf8PathError(String),

    /// JPL ephemeris file not found at the expected location.
    #[error("JPL ephemeris file not found: {0}")]
    JPLFileNotFound(String),

    /// Failure during polynomial or numerical root finding (from `roots` crate).
    #[error("Numerical root finding failed: {0}")]
    RootFindingError(#[from] roots::SearchError),

    /// Observation not found (index out of range or missing data).
    #[error("Observation not found at index: {0}")]
    ObservationNotFound(usize),

    /// Invalid error model specification (e.g., unknown identifier).
    #[error("Invalid error model identifier: {0}")]
    InvalidErrorModel(String),

    /// Invalid file path provided for an error model file.
    #[error("Invalid error model file path: {0}")]
    InvalidErrorModelFilePath(String),

    /// Error during parsing with `nom` parser combinators.
    #[error("Parsing error (nom): {0}")]
    NomParsingError(String),

    /// Error while parsing legacy 80-column formatted observation files.
    #[error("Parsing error in 80-column observation file: {0}")]
    Parsing80ColumnFileError(ParseObsError),

    /// Gaussian noise generation failed (invalid parameters for normal distribution).
    #[error("Gaussian noise generation failed: {0:?}")]
    NoiseInjectionError(rand_distr::NormalError),

    /// The unit direction matrix is singular (cannot invert; likely coplanar observations).
    #[error("Singular direction matrix (cannot invert); observations may be coplanar")]
    SingularDirectionMatrix,

    /// Aberth–Ehrlich method failed to find acceptable complex polynomial roots.
    #[error("Polynomial root finding failed (Aberth–Ehrlich method did not converge)")]
    PolynomialRootFindingFailed,

    /// Spurious root detected (e.g., negative or near-zero geocentric distance).
    #[error("Spurious root detected (negative or near-zero geocentric distance)")]
    SpuriousRootDetected,

    /// Gauss method for initial orbit determination failed to find roots.
    #[error("Initial orbit determination (Gauss method) failed to find valid roots")]
    GaussNoRootsFound,

    /// Invalid SPK segment data type in JPL ephemeris file.
    #[error("Invalid SPK segment data type: {0}")]
    InvalidSpkDataType(i32),

    /// Invalid input parameter for initial orbit determination.
    #[error("Invalid parameter for initial orbit determination: {0}")]
    InvalidIODParameter(String),

    /// Invalid reference system or unsupported transformation request.
    #[error("Invalid reference system: {0}")]
    InvalidRefSystem(String),

    /// Failure during velocity correction procedure.
    #[error("Velocity correction procedure failed: {0}")]
    VelocityCorrectionError(String),

    /// Invalid or inconsistent orbital elements/state vector.
    #[error("Invalid orbital state or inconsistent elements: {0}")]
    InvalidOrbit(String),

    /// Generic invalid input error with a descriptive message.
    #[error("Invalid input conversion: {0}")]
    InvalidConversion(String),

    /// Error indicating a floating-point value is NaN (Not a Number).
    #[error("Invalid floating-point value (NaN encountered): {0}")]
    InvalidFloatValue(ordered_float::FloatIsNan),

    /// RMS computation failed (e.g., no valid observations).
    #[error("RMS computation failed: {0}")]
    RmsComputationFailed(String),
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
