use camino::{Utf8Path, Utf8PathBuf};
use directories::BaseDirs;
use std::{
    fs,
    io::{self},
    str::FromStr,
};

#[cfg(feature = "jpl-download")]
use tokio::{fs::File, io::AsyncWriteExt};
#[cfg(feature = "jpl-download")]
use tokio_stream::StreamExt;

use super::{horizon::horizon_version::JPLHorizonVersion, naif::naif_version::NAIFVersion};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(Debug, Clone)]
pub enum EphemFileSource {
    JPLHorizon(JPLHorizonVersion),
    NAIF(NAIFVersion),
}

impl TryFrom<&str> for EphemFileSource {
    type Error = String;

    fn try_from(value: &str) -> std::result::Result<EphemFileSource, std::string::String> {
        let parts: Vec<&str> = value.split(':').collect();
        if parts.len() != 2 {
            return Err("Expected format: {source}:{version}, example: 'naif:DE440' or 'horizon:DE440'".to_string());
        }

        match parts[0].to_lowercase().as_str() {
            "horizon" => {
                if let Ok(version) = JPLHorizonVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::JPLHorizon(version))
                } else {
                    Err(format!("Invalid JPL Horizon version: {}", parts[1]))
                }
            }
            "naif" => {
                if let Ok(version) = NAIFVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::NAIF(version))
                } else {
                    Err(format!("Invalid NAIF version: {}", parts[1]))
                }
            }
            _ => Err(format!("Unknown ephemeris file source: {}", parts[0])),
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
async fn download_big_file(url: &str, path: &Utf8Path) -> Result<()> {
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
    /// Get the ephemeris file from the JPL
    ///
    /// Arguments
    /// ---------
    /// * `env_state`: the Outfit environment state containing the http_client
    /// * `version`: the version of the ephemeris file to download, example: "de442" or "de442s"
    /// * `user_path`: an optional user-provided path to the ephemeris file
    ///
    /// Return
    /// ------
    /// * The path to the ephemeris file
    /// * An error if the file cannot be found or downloaded
    pub fn get_ephemeris_file(file_source: &EphemFileSource) -> io::Result<Self> {
        let local_file = EphemFilePath::try_from(file_source.clone())
            .map_err(|e| io::Error::new(io::ErrorKind::Other, e))?;

        if local_file.exists() {
            return Ok(local_file);
        } else {
            #[cfg(feature = "jpl-download")]
            {
                let url = file_source.get_version_url();

                let rt = tokio::runtime::Runtime::new().expect("Failed to create runtime");
                rt.block_on(async { download_big_file(&url, &local_file.path()).await })
                    .map_err(|_| {
                        io::Error::new(
                            io::ErrorKind::Other,
                            "Failed to download JPL ephemeris file",
                        )
                    })?;

                return Ok(local_file);
            }
            #[cfg(not(feature = "jpl-download"))]
            {
                return Err(io::Error::new(
                    io::ErrorKind::NotFound,
                    "JPL ephemeris file not found",
                ));
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
    type Error = String;

    fn try_from(value: EphemFileSource) -> std::result::Result<Self, Self::Error> {
        let base_dir =
            BaseDirs::new().ok_or_else(|| "Failed to retrieve base directory".to_string())?;
        let cache_path = Utf8Path::from_path(base_dir.cache_dir())
            .ok_or_else(|| "Failed to convert base directory path to Utf8Path".to_string())?;
        let cache_path = cache_path.join("outfit_cache").join("jpl_ephem");
        let cache_dir = value.get_cache_dir();
        let cache_path = cache_path.join(cache_dir);
        fs::create_dir_all(&cache_path).map_err(|e| e.to_string())?;

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
    use super::*;

    #[test]
    #[cfg(not(feature = "jpl-download"))]
    fn test_no_feature_download_jpl_ephem() {
        let file_source = "naif:DE442".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(&file_source);
        assert!(
            result.is_err(),
            "feature jpl-download is enabled, weird ..."
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_feature_download_jpl_ephem_from_naif() {
        let file_source = "naif:DE442".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(&file_source);
        assert!(result.is_ok(), "Failed to download JPL ephemeris file");
        let path = result.unwrap();
        assert!(path.exists(), "JPL ephemeris file not found");
        assert!(
            path.extension().unwrap() == "bsp",
            "JPL ephemeris file is not a .bsp file"
        );
        assert!(
            path.file_name().unwrap() == "de442.bsp",
            "JPL ephemeris file is not de442.bsp"
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_feature_download_jpl_ephem_from_horizon() {
        let file_source = "horizon:DE440".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(&file_source);
        assert!(result.is_ok(), "Failed to download JPL ephemeris file");
        let path = result.unwrap();
        assert!(path.exists(), "JPL ephemeris file not found");
        assert!(
            path.extension().unwrap() == "bsp",
            "JPL ephemeris file is not a .bsp file"
        );
        assert!(
            path.file_name().unwrap() == "DE440.bsp",
            "JPL ephemeris file is not DE440.bsp"
        );
    }
}
