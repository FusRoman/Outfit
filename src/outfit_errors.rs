use thiserror::Error;

use crate::observations::observations::ParseObsError;

#[derive(Error, Debug)]
pub enum OutfitError {
    #[error("Invalid JPL string format: {0}")]
    InvalidJPLStringFormat(String),

    #[error("Invalid JPL ephemeris file source: {0}")]
    InvalidJPLEphemFileSource(String),

    #[error("Invalid JPL ephemeris file version: {0}")]
    InvalidJPLEphemFileVersion(String),

    #[error("Invalid URL: {0}")]
    InvalidUrl(String),

    #[error("HTTP ureq error: {0}")]
    UreqHttpError(#[from] ureq::Error),

    #[error("Unable to perform file operation: {0}")]
    IoError(#[from] std::io::Error),

    #[cfg(feature = "jpl-download")]
    #[error("HTTP reqwest error: {0}")]
    ReqwestError(#[from] reqwest::Error),

    #[error("Base dir creation error for JPL ephemeris file: {0}")]
    UnableToCreateBaseDir(String),

    #[error("UTF-8 Path error: {0}")]
    Utf8PathError(String),

    #[error("JPL File not found at: {0}")]
    JPLFileNotFound(String),

    #[error("ROOTS finding error: {0}")]
    RootFindingError(#[from] roots::SearchError),

    #[error("Observation not found: {0}")]
    ObservationNotFound(usize),

    #[error("Invalid Error model: {0}")]
    InvalidErrorModel(String),

    #[error("Invalid File Path for Error model: {0}")]
    InvalidErrorModelFilePath(String),

    #[error("Error during the nom parsing: {0}")]
    NomParsingError(String),

    #[error("Error during the 80 column file parsing: {0}")]
    Parsing80ColumnFileError(ParseObsError),

    #[error("Gaussian noise generation failed: {0:?}")]
    NoiseInjectionError(rand_distr::NormalError),

    #[error(
        "Unit direction matrix is singular (cannot be inverted); observations may be coplanar"
    )]
    SingularDirectionMatrix,

    #[error("Aberth–Ehrlich method failed to find acceptable complex roots")]
    PolynomialRootFindingFailed,

    #[error("Spurious root detected (e.g., negative or near-zero geocentric distance)")]
    SpuriousRootDetected,

    #[error("Gauss method failed to find roots")]
    GaussNoRootsFound,

    #[error("Invalid SPK data type: {0}")]
    InvalidSpkDataType(i32),
}

impl From<rand_distr::NormalError> for OutfitError {
    fn from(err: rand_distr::NormalError) -> Self {
        OutfitError::NoiseInjectionError(err)
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

            // Ces erreurs ne sont pas comparables : égalité si même variant
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

            // Variantes unitaires
            (SingularDirectionMatrix, SingularDirectionMatrix) => true,
            (PolynomialRootFindingFailed, PolynomialRootFindingFailed) => true,
            (SpuriousRootDetected, SpuriousRootDetected) => true,
            (GaussNoRootsFound, GaussNoRootsFound) => true,

            _ => false,
        }
    }
}
