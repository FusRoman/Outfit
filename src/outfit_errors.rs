use thiserror::Error;

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
}
