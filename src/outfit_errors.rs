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
//! **When**: Parsing observations, fixed-width records, ADES, or internal streams.
//!
//! - [`NomParsingError`](crate::outfit_errors::OutfitError::NomParsingError)
//! - [`Parsing80ColumnFileError`](crate::outfit_errors::OutfitError::Parsing80ColumnFileError)
//! - [`Parquet`](crate::outfit_errors::OutfitError::Parquet)
//!
//! **Typical causes**: schema drift; corrupted inputs; locale-specific formats;
//! columnar format mismatches.
//!
//! **Remediation**: strengthen parsers; surface line/column context; round-trip with golden files;
//! validate Parquet schema and logical types.
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
//! - [`RmsComputationFailed`](crate::outfit_errors::OutfitError::RmsComputationFailed)
//! - [`GaussPrelimOrbitFailed`](crate::outfit_errors::OutfitError::GaussPrelimOrbitFailed)
//! - [`NoViableOrbit`](crate::outfit_errors::OutfitError::NoViableOrbit)
//! - [`NoFeasibleTriplets`](crate::outfit_errors::OutfitError::NoFeasibleTriplets)
//!
//! **Typical causes**: coplanar geometry; invalid orbital elements; unsupported frame conversions;
//! all candidate triplets/realizations failing to produce a viable orbit.
//!
//! **Remediation**: pre-filter observations; enforce numeric bounds; fall back to alternative solvers;
//! surface last/aggregated failure details when no orbit can be determined.
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
//! - Variants wrapping opaque errors (e.g., I/O, HTTP, Parquet) compare equal by kind only.
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

use crate::observations::ParseObsError;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum OutfitError {
    #[error("Invalid JPL string format: {0}")]
    InvalidJPLStringFormat(String),

    #[error("Invalid JPL ephemeris file source: {0}")]
    InvalidJPLEphemFileSource(String),

    #[error("Invalid JPL ephemeris file version: {0}")]
    InvalidJPLEphemFileVersion(String),

    #[error("Invalid URL string: {0}")]
    InvalidUrl(String),

    #[error("HTTP request failed (ureq): {0}")]
    UreqHttpError(#[from] ureq::Error),

    #[error("Filesystem I/O error: {0}")]
    IoError(#[from] std::io::Error),

    #[cfg(feature = "jpl-download")]
    #[error("HTTP request failed (reqwest): {0}")]
    ReqwestError(#[from] reqwest::Error),

    #[error("Failed to create base directory for JPL ephemeris file: {0}")]
    UnableToCreateBaseDir(String),

    #[error("Filesystem path is not valid UTF-8: {0}")]
    Utf8PathError(String),

    #[error("JPL ephemeris file not found: {0}")]
    JPLFileNotFound(String),

    #[error("Numerical root finding failed: {0}")]
    RootFindingError(#[from] roots::SearchError),

    #[error("Observation not found at index: {0}")]
    ObservationNotFound(usize),

    #[error("Invalid error model identifier: {0}")]
    InvalidErrorModel(String),

    #[error("Invalid error model file path: {0}")]
    InvalidErrorModelFilePath(String),

    #[error("Parsing error (nom): {0}")]
    NomParsingError(String),

    #[error("Parsing error in 80-column observation file: {0}")]
    Parsing80ColumnFileError(ParseObsError),

    #[error("Gaussian noise generation failed: {0:?}")]
    NoiseInjectionError(rand_distr::NormalError),

    #[error("Singular direction matrix (cannot invert); observations may be coplanar")]
    SingularDirectionMatrix,

    #[error("Polynomial root finding failed (Aberth–Ehrlich method did not converge)")]
    PolynomialRootFindingFailed,

    #[error("Spurious root detected (negative or near-zero geocentric distance)")]
    SpuriousRootDetected,

    #[error("Initial orbit determination (Gauss method) failed to find valid roots")]
    GaussNoRootsFound,

    #[error("Invalid SPK segment data type: {0}")]
    InvalidSpkDataType(i32),

    #[error("Invalid parameter for initial orbit determination: {0}")]
    InvalidIODParameter(String),

    #[error("Invalid reference system: {0}")]
    InvalidRefSystem(String),

    #[error("Velocity correction procedure failed: {0}")]
    VelocityCorrectionError(String),

    #[error("Invalid orbital state or inconsistent elements: {0}")]
    InvalidOrbit(String),

    #[error("Invalid input conversion: {0}")]
    InvalidConversion(String),

    #[error("Invalid floating-point value (NaN encountered): {0}")]
    InvalidFloatValue(ordered_float::FloatIsNan),

    #[error("RMS computation failed: {0}")]
    RmsComputationFailed(String),

    #[error("Gauss preliminary orbit determination failed: {0}")]
    GaussPrelimOrbitFailed(String),

    #[error(transparent)]
    Parquet(#[from] parquet::errors::ParquetError),

    #[error("No viable orbit could be determined after {attempts} attempts: {cause}")]
    NoViableOrbit {
        cause: Box<OutfitError>,
        attempts: usize,
    },

    #[error(
        "No feasible triplets (span={span:.6} d, n_obs={n_obs}, dt_min={dt_min}, dt_max={dt_max})"
    )]
    NoFeasibleTriplets {
        span: f64,
        n_obs: usize,
        dt_min: f64,
        dt_max: f64,
    },
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

            // Opaque external error kinds compare equal by variant only
            (UreqHttpError(_), UreqHttpError(_)) => true,
            (IoError(_), IoError(_)) => true,
            #[cfg(feature = "jpl-download")]
            (ReqwestError(_), ReqwestError(_)) => true,
            (Parquet(_), Parquet(_)) => true,

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
            (InvalidFloatValue(a), InvalidFloatValue(b)) => a == b,
            (RmsComputationFailed(a), RmsComputationFailed(b)) => a == b,
            (GaussPrelimOrbitFailed(a), GaussPrelimOrbitFailed(b)) => a == b,
            (
                NoViableOrbit {
                    cause: a,
                    attempts: na,
                },
                NoViableOrbit {
                    cause: b,
                    attempts: nb,
                },
            ) => a == b && na == nb,
            (
                NoFeasibleTriplets {
                    span: a,
                    n_obs: na,
                    dt_min: da,
                    dt_max: ma,
                },
                NoFeasibleTriplets {
                    span: b,
                    n_obs: nb,
                    dt_min: db,
                    dt_max: mb,
                },
            ) => a == b && na == nb && da == db && ma == mb,

            // Unit-like variants
            (SingularDirectionMatrix, SingularDirectionMatrix) => true,
            (PolynomialRootFindingFailed, PolynomialRootFindingFailed) => true,
            (SpuriousRootDetected, SpuriousRootDetected) => true,
            (GaussNoRootsFound, GaussNoRootsFound) => true,

            _ => false,
        }
    }
}
