use camino::{Utf8Path, Utf8PathBuf};
use directories::BaseDirs;
use std::{fs, io::{self}};

use crate::outfit::Outfit;
use tokio::{fs::File, io::AsyncWriteExt};

use tokio_stream::StreamExt;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

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
pub fn get_ephemeris_file(
    env_state: &Outfit,
    version: Option<&str>,
    user_path: Option<&str>,
) -> io::Result<Utf8PathBuf> {
    // If user_path is provided, check if the file exists
    if let Some(path_str) = user_path {
        let path = Utf8Path::new(path_str);
        if path.exists() {
            return Ok(path.to_path_buf());
        } else {
            return Err(io::Error::new(
                io::ErrorKind::NotFound,
                "Fichier JPL non trouvé à ce chemin.",
            ));
        }
    }

    #[cfg(feature = "jpl-download")]
    {
        use super::download_jpl_file::download_big_file;

        let base_dir = BaseDirs::new().expect("Cannot find the base directory");

        // otherwise use the default cache directory and download the file
        let cache_path = Utf8Path::from_path(base_dir.cache_dir()).expect("Invalid cache path");
        let cache_path = cache_path.join("outfit_cache").join("jpl_ephem");
        fs::create_dir_all(&cache_path)?;

        let Some(version) = version else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                "Version de l'éphéméride non spécifiée.",
            ));
        };

        let filename = format!("{}.bsp", version);
        let local_file = cache_path.join(&filename);

        if local_file.exists() {
            return Ok(local_file);
        }

        let url = format!(
            "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/{}.bsp",
            version
        );

        let rt = tokio::runtime::Runtime::new().expect("Failed to create runtime");
        rt.block_on(async { download_big_file(&url, &local_file).await })
            .map_err(|_| {
                io::Error::new(
                    io::ErrorKind::Other,
                    "Failed to download JPL ephemeris file",
                )
            })?;

        Ok(local_file)
    }

    #[cfg(not(feature = "jpl-download"))]
    {
        Err(io::Error::new(
            io::ErrorKind::NotFound,
            "JPL ephemeris file not found and download feature is disabled.",
        ))
    }
}

#[cfg(test)]
mod jpl_reader_test {
    use super::*;

    #[test]
    #[cfg(not(feature = "jpl-download"))]
    fn test_no_feature_download_jpl_ephem() {
        let env_state = Outfit::new();
        let version = Some("de442");
        let user_path = None;

        let result = get_ephemeris_file(&env_state, version, user_path);
        assert!(
            result.is_err(),
            "feature jpl-download is enabled, weird ..."
        );
    }

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_feature_download_jpl_ephem() {
        use crate::outfit::Outfit;

        let env_state = Outfit::new();
        let version = Some("de442");
        let user_path = None;

        let result = get_ephemeris_file(&env_state, version, user_path);
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
        fs::remove_file(path).expect("Failed to remove JPL ephemeris file after test");
    }
}
