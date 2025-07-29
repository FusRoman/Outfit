use camino::{Utf8Path, Utf8PathBuf};
use directories::BaseDirs;
use std::{fs, str::FromStr};

#[cfg(feature = "jpl-download")]
use tokio::{fs::File, io::AsyncWriteExt};
#[cfg(feature = "jpl-download")]
use tokio_stream::StreamExt;

use super::{horizon::horizon_version::JPLHorizonVersion, naif::naif_version::NAIFVersion};
use crate::outfit_errors::OutfitError;

#[derive(Debug, Clone)]
pub enum EphemFileSource {
    JPLHorizon(JPLHorizonVersion),
    NAIF(NAIFVersion),
}

impl TryFrom<&str> for EphemFileSource {
    type Error = OutfitError;

    fn try_from(value: &str) -> Result<Self, Self::Error> {
        let parts: Vec<&str> = value.split(':').collect();
        if parts.len() != 2 {
            return Err(OutfitError::InvalidJPLStringFormat(
                "Expected format: {source}:{version}, example: 'naif:DE440' or 'horizon:DE440'"
                    .to_string(),
            ));
        }

        match parts[0].to_lowercase().as_str() {
            "horizon" => {
                if let Ok(version) = JPLHorizonVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::JPLHorizon(version))
                } else {
                    Err(OutfitError::InvalidJPLEphemFileVersion(format!(
                        "Invalid JPL Horizon version: {}",
                        parts[1]
                    )))
                }
            }
            "naif" => {
                if let Ok(version) = NAIFVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::NAIF(version))
                } else {
                    Err(OutfitError::InvalidJPLEphemFileVersion(format!(
                        "Invalid NAIF version: {}",
                        parts[1]
                    )))
                }
            }
            _ => Err(OutfitError::InvalidJPLEphemFileSource(format!(
                "Unknown ephemeris file source: {}",
                parts[0]
            ))),
        }
    }
}

impl EphemFileSource {
    fn get_baseurl(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(_) => "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/",
            EphemFileSource::NAIF(_) => {
                "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
            }
        }
    }

    fn get_version_url(&self) -> String {
        let base_url = self.get_baseurl();
        match self {
            EphemFileSource::JPLHorizon(jpl_version) => {
                format!("{}{}", base_url, jpl_version.get_filename())
            }
            EphemFileSource::NAIF(version) => format!("{}{}", base_url, version.get_filename()),
        }
    }

    fn get_cache_dir(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(_) => "jpl_horizon",
            EphemFileSource::NAIF(_) => "naif",
        }
    }

    fn filename(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(version) => version.to_filename(),
            EphemFileSource::NAIF(version) => version.get_filename(),
        }
    }
}

/// Download a large file from a URL
/// Uses reqwest to download the file in chunks
/// and saves it to the specified path using tokio's async file I/O
/// and stream processing.
///
/// Arguments
/// ---------
/// * `url`: the URL of the file to download
/// * `path`: the path to save the downloaded file
///
/// Return
/// ------
/// * An error if the download fails
/// * Ok(()) if the download is successful
#[cfg(feature = "jpl-download")]
async fn download_big_file(url: &str, path: &Utf8Path) -> Result<(), OutfitError> {
    let mut file = File::create(path).await?;
    println!("Downloading {}...", url);

    let mut stream = reqwest::get(url).await?.bytes_stream();

    while let Some(chunk_result) = stream.next().await {
        let chunk = chunk_result?;
        file.write_all(&chunk).await?;
    }

    file.flush().await?;

    println!("Downloaded {}", url);
    Ok(())
}

pub enum EphemFilePath {
    JPLHorizon(Utf8PathBuf, JPLHorizonVersion),
    NAIF(Utf8PathBuf, NAIFVersion),
}

impl std::fmt::Display for EphemFilePath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EphemFilePath::JPLHorizon(path, version) => {
                write!(f, "JPL Horizon: {} ({})", path, version.get_filename())
            }
            EphemFilePath::NAIF(path, version) => {
                write!(f, "NAIF: {} ({})", path, version.get_filename())
            }
        }
    }
}

impl EphemFilePath {
    /// Get the path of the ephemeris file.
    /// The file should be located in the OS cache directory.
    /// On linux, the cache directory is ~/.cache/outfit_cache/jpl_ephem
    /// On windows, the cache directory is C:\Users\<username>\AppData\Local\outfit_cache\jpl_ephem
    ///
    /// The full path after jpl_ephem should be:
    /// * jpl_horizon/<filename> for JPL Horizon files
    /// * naif/<filename> for NAIF files
    ///
    /// If the file does not exist and the feature 'jpl-download' is activated,
    /// the file will be downloaded from the JPL website.
    /// Otherwise, an error will be returned.
    ///
    /// Arguments
    /// ---------
    /// * `file_source`: the source of the ephemeris file
    ///
    /// Return
    /// ------
    /// * The path to the ephemeris file
    /// * An error if the file does not exist and the feature 'jpl-download' is not activated
    pub fn get_ephemeris_file(file_source: &EphemFileSource) -> Result<EphemFilePath, OutfitError> {
        let local_file = EphemFilePath::try_from(file_source.clone())?;

        if local_file.exists() {
            Ok(local_file)
        } else {
            #[cfg(feature = "jpl-download")]
            {
                let url = file_source.get_version_url();

                let rt = tokio::runtime::Runtime::new().expect("Failed to create runtime");
                rt.block_on(async { download_big_file(&url, &local_file.path()).await })?;

                return Ok(local_file);
            }
            #[cfg(not(feature = "jpl-download"))]
            {
                Err(OutfitError::JPLFileNotFound(local_file.path().to_string()))
            }
        }
    }

    pub fn exists(&self) -> bool {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.exists(),
            EphemFilePath::NAIF(path, _) => path.exists(),
        }
    }

    pub fn path(&self) -> &Utf8Path {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path,
            EphemFilePath::NAIF(path, _) => path,
        }
    }

    pub fn file_name(&self) -> Option<&str> {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.file_name(),
            EphemFilePath::NAIF(path, _) => path.file_name(),
        }
    }

    pub fn extension(&self) -> Option<&str> {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.extension(),
            EphemFilePath::NAIF(path, _) => path.extension(),
        }
    }
}

impl TryFrom<EphemFileSource> for EphemFilePath {
    type Error = OutfitError;

    fn try_from(value: EphemFileSource) -> Result<Self, Self::Error> {
        let base_dir = BaseDirs::new().ok_or_else(|| {
            OutfitError::UnableToCreateBaseDir("Failed to retrieve base directory".to_string())
        })?;

        let cache_path = Utf8Path::from_path(base_dir.cache_dir()).ok_or_else(|| {
            OutfitError::Utf8PathError(
                "Failed to convert base directory path to Utf8Path".to_string(),
            )
        })?;
        let cache_path = cache_path.join("outfit_cache").join("jpl_ephem");
        let cache_dir = value.get_cache_dir();
        let cache_path = cache_path.join(cache_dir);
        fs::create_dir_all(&cache_path)?;

        let local_file = cache_path.join(value.filename());

        match value {
            EphemFileSource::JPLHorizon(version) => {
                Ok(EphemFilePath::JPLHorizon(local_file, version))
            }
            EphemFileSource::NAIF(version) => Ok(EphemFilePath::NAIF(local_file, version)),
        }
    }
}

#[cfg(test)]
mod jpl_reader_test {

    #[test]
    #[cfg(not(feature = "jpl-download"))]
    fn test_no_feature_download_jpl_ephem() {
        use super::*;
        let file_source = "naif:DE442".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(&file_source);
        assert!(
            result.is_err(),
            "feature jpl-download is enabled, weird ..."
        );
    }
}
