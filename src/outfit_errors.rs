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
}

impl From<rand_distr::NormalError> for OutfitError {
    fn from(err: rand_distr::NormalError) -> Self {
        OutfitError::NoiseInjectionError(err)
    }
}
