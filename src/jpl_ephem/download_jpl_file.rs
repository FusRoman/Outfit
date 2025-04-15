use camino::{Utf8Path, Utf8PathBuf};
use directories::BaseDirs;
use std::{
    fs,
    io::{self},
    str::FromStr,
};

use crate::outfit::Outfit;
#[cfg(feature = "jpl-download")]
use tokio::{fs::File, io::AsyncWriteExt};
#[cfg(feature = "jpl-download")]
use tokio_stream::StreamExt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum JPLHorizonVersion {
    DE102,
    DE200,
    DE202,
    DE403,
    DE405,
    DE406,
    DE410,
    DE413,
    DE414,
    DE418,
    DE421,
    DE422,
    DE423,
    DE430,
    DE430t,
    DE431,
    DE440,
    DE440t,
    DE441,
}

impl JPLHorizonVersion {
    fn get_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "de102/lnxm1410p3002.102",
            JPLHorizonVersion::DE200 => "de200/lnxm1600p2170.200",
            JPLHorizonVersion::DE202 => "de202/lnxp1900p2050.202",
            JPLHorizonVersion::DE403 => "de403/lnxp1600p2200.403",
            JPLHorizonVersion::DE405 => "de405/lnxp1600p2200.405",
            JPLHorizonVersion::DE406 => "de406/lnxm3000p3000.406",
            JPLHorizonVersion::DE410 => "de410/lnxp1960p2020.410",
            JPLHorizonVersion::DE413 => "de413/lnxp1900p2050.413",
            JPLHorizonVersion::DE414 => "de414/lnxp1600p2200.414",
            JPLHorizonVersion::DE418 => "de418/lnxp1900p2050.418",
            JPLHorizonVersion::DE421 => "de421/lnxp1900p2053.421",
            JPLHorizonVersion::DE422 => "de422/lnxm3000p3000.422",
            JPLHorizonVersion::DE423 => "de423/lnxp1800p2200.423",
            JPLHorizonVersion::DE430 => "de430/linux_p1550p2650.430",
            JPLHorizonVersion::DE430t => "de430t/linux_p1550p2650.430t",
            JPLHorizonVersion::DE431 => "de431/lnxm13000p17000.431",
            JPLHorizonVersion::DE440 => "de440/linux_p1550p2650.440",
            JPLHorizonVersion::DE440t => "de440t/linux_p1550p2650.440t",
            JPLHorizonVersion::DE441 => "de441/linux_m13000p17000.441",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE102" => Some(JPLHorizonVersion::DE102),
            "DE200" => Some(JPLHorizonVersion::DE200),
            "DE202" => Some(JPLHorizonVersion::DE202),
            "DE403" => Some(JPLHorizonVersion::DE403),
            "DE405" => Some(JPLHorizonVersion::DE405),
            "DE406" => Some(JPLHorizonVersion::DE406),
            "DE410" => Some(JPLHorizonVersion::DE410),
            "DE413" => Some(JPLHorizonVersion::DE413),
            "DE414" => Some(JPLHorizonVersion::DE414),
            "DE418" => Some(JPLHorizonVersion::DE418),
            "DE421" => Some(JPLHorizonVersion::DE421),
            "DE422" => Some(JPLHorizonVersion::DE422),
            "DE423" => Some(JPLHorizonVersion::DE423),
            "DE430" => Some(JPLHorizonVersion::DE430),
            "DE430t" => Some(JPLHorizonVersion::DE430t),
            "DE431" => Some(JPLHorizonVersion::DE431),
            "DE440" => Some(JPLHorizonVersion::DE440),
            "DE440t" => Some(JPLHorizonVersion::DE440t),
            "DE441" => Some(JPLHorizonVersion::DE441),
            _ => None,
        }
    }

    fn to_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "DE102.bsp",
            JPLHorizonVersion::DE200 => "DE200.bsp",
            JPLHorizonVersion::DE202 => "DE202.bsp",
            JPLHorizonVersion::DE403 => "DE403.bsp",
            JPLHorizonVersion::DE405 => "DE405.bsp",
            JPLHorizonVersion::DE406 => "DE406.bsp",
            JPLHorizonVersion::DE410 => "DE410.bsp",
            JPLHorizonVersion::DE413 => "DE413.bsp",
            JPLHorizonVersion::DE414 => "DE414.bsp",
            JPLHorizonVersion::DE418 => "DE418.bsp",
            JPLHorizonVersion::DE421 => "DE421.bsp",
            JPLHorizonVersion::DE422 => "DE422.bsp",
            JPLHorizonVersion::DE423 => "DE423.bsp",
            JPLHorizonVersion::DE430 => "DE430.bsp",
            JPLHorizonVersion::DE430t => "DE430t.bsp",
            JPLHorizonVersion::DE431 => "DE431.bsp",
            JPLHorizonVersion::DE440 => "DE440.bsp",
            JPLHorizonVersion::DE440t => "DE440t.bsp",
            JPLHorizonVersion::DE441 => "DE441.bsp",
        }
    }
}

impl FromStr for JPLHorizonVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        JPLHorizonVersion::from_str(s).ok_or_else(|| format!("Invalid JPL Horizon version: {}", s))
    }
}

#[derive(Debug, Clone)]
pub enum NAIFVersion {
    DE430,
    DE431p1,
    DE431p2,
    DE432,
    DE435,
    DE438,
    DE440,
    DE440s,
    DE441p1,
    DE441p2,
    DE442,
}

impl NAIFVersion {
    fn get_filename(&self) -> &str {
        match self {
            NAIFVersion::DE430 => "de430.bsp",
            NAIFVersion::DE431p1 => "de431_part-1.bsp",
            NAIFVersion::DE431p2 => "de431_part-2.bsp",
            NAIFVersion::DE432 => "de432.bsp",
            NAIFVersion::DE435 => "de435.bsp",
            NAIFVersion::DE438 => "de438.bsp",
            NAIFVersion::DE440 => "de440.bsp",
            NAIFVersion::DE440s => "de440s.bsp",
            NAIFVersion::DE441p1 => "de441_part-1.bsp",
            NAIFVersion::DE441p2 => "de441_part-2.bsp",
            NAIFVersion::DE442 => "de442.bsp",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE430" => Some(NAIFVersion::DE430),
            "DE431_part-1" => Some(NAIFVersion::DE431p1),
            "DE431_part-2" => Some(NAIFVersion::DE431p2),
            "DE432" => Some(NAIFVersion::DE432),
            "DE435" => Some(NAIFVersion::DE435),
            "DE438" => Some(NAIFVersion::DE438),
            "DE440" => Some(NAIFVersion::DE440),
            "DE440s" => Some(NAIFVersion::DE440s),
            "DE441_part-1" => Some(NAIFVersion::DE441p1),
            "DE441_part-2" => Some(NAIFVersion::DE441p2),
            "DE442" => Some(NAIFVersion::DE442),
            _ => None,
        }
    }
}

impl FromStr for NAIFVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        NAIFVersion::from_str(s).ok_or_else(|| format!("Invalid NAIF version: {}", s))
    }
}

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
            return Err("Expected format: 'naif:DE440' or 'jpl:DE440t'".to_string());
        }

        match parts[0].to_lowercase().as_str() {
            "horizon" => {
                if let Some(version) = JPLHorizonVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::JPLHorizon(version))
                } else {
                    Err(format!("Invalid JPL Horizon version: {}", parts[1]))
                }
            }
            "naif" => {
                if let Some(version) = NAIFVersion::from_str(parts[1]) {
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
    #[cfg(feature = "jpl-download")]
    pub fn get_ephemeris_file(file_source: EphemFileSource) -> io::Result<Self> {
        let local_file = EphemFilePath::try_from(file_source.clone())
            .map_err(|_| io::Error::new(io::ErrorKind::NotFound, "JPL ephemeris file not found"))?;

        if local_file.exists() {
            return Ok(local_file);
        }

        let url = file_source.get_version_url();

        let rt = tokio::runtime::Runtime::new().expect("Failed to create runtime");
        rt.block_on(async { download_big_file(&url, &local_file.path()).await })
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to download JPL ephemeris file",
                )
            })?;

        Ok(local_file)
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
        let base_dir = BaseDirs::new().expect("Cannot find the base directory");
        let cache_path = Utf8Path::from_path(base_dir.cache_dir()).expect("Invalid cache path");
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
        let version = Some("de442");
        let user_path = None;

        let result = get_ephemeris_file(version, user_path);
        assert!(
            result.is_err(),
            "feature jpl-download is enabled, weird ..."
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_feature_download_jpl_ephem_from_naif() {
        let file_source = "naif:DE442".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(file_source);
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
        fs::remove_file(path.path()).expect("Failed to remove JPL ephemeris file after test");
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_feature_download_jpl_ephem_from_horizon() {
        let file_source = "horizon:DE440".try_into().unwrap();

        let result = EphemFilePath::get_ephemeris_file(file_source);
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
        fs::remove_file(path.path()).expect("Failed to remove JPL ephemeris file after test");
    }
}
