//! Ephemeris file resolution, caching, and (optional) download.
//!
//! This module determines where to find the JPL ephemeris file required by the
//! `jpl` layer, placing a unified cache in the user's OS cache directory and,
//! when the `jpl-download` feature is enabled, downloading the missing file
//! from the official JPL locations (Horizons legacy DE binaries and NAIF SPK/DAF).
//!
//! # What this module does
//! - Parse a high-level ephemeris source specification (backend + version).
//! - Resolve a **cache path** under the OS cache directory (via `directories`).
//! - Materialize a typed handle [`EphemFilePath`] used by readers (`horizon`, `naif`).
//! - Optionally **download** the file if it is not present (feature `jpl-download`).
//!
//! # Cache layout
//! Cache root: `<os-cache>/outfit_cache/jpl_ephem`
//!
//! Subdirectories:
//! - `jpl_horizon/` — legacy JPL DE binaries (TTL/CNAM/IPT layout).
//! - `naif/`       — NAIF SPK/DAF kernels.
//!
//! File names are provided by the version enums
//! ([`crate::jpl_ephem::horizon::horizon_version::JPLHorizonVersion::to_filename`],
//!  [`crate::jpl_ephem::naif::naif_version::NaifVersion::get_filename`]) and match official names.
//!
//! # String parsing
//! You can create a source from a string with `TryFrom<&str>`:
//!
//! - `horizon:DE440`  → legacy DE binary
//! - `naif:DE440`     → NAIF SPK/DAF
//!
//! # Platform paths (examples)
//! - **Linux:** `~/.cache/outfit_cache/jpl_ephem/...`
//! - **Windows:** `C:\Users\<user>\AppData\Local\outfit_cache\jpl_ephem\...`
//! - **macOS:** `~/Library/Caches/outfit_cache/jpl_ephem/...`
//!
//! # See also
//! * [`crate::jpl_ephem::horizon`] — Legacy DE binary reader.
//! * [`crate::jpl_ephem::naif`] — NAIF SPK/DAF reader.
//! * [`EphemFileSource`] — User-facing source specification.
//! * [`EphemFilePath`] — Resolved on-disk path and version.

use camino::{Utf8Path, Utf8PathBuf};
use directories::BaseDirs;
use std::{fs, str::FromStr};

#[cfg(feature = "jpl-download")]
use tokio::{fs::File, io::AsyncWriteExt};
#[cfg(feature = "jpl-download")]
use tokio_stream::StreamExt;

use super::{horizon::horizon_version::JPLHorizonVersion, naif::naif_version::NaifVersion};
use crate::outfit_errors::OutfitError;

/// Ephemeris file source selector (backend + version).
///
/// This enum encodes **what** file we want, independent from **where** it lives
/// on disk. Use [`EphemFilePath::get_ephemeris_file`] to resolve it.
///
/// Variants
/// --------
/// - [`EphemFileSource::JPLHorizon`] — Legacy JPL DE binaries (e.g. `DE440`).
/// - [`EphemFileSource::Naif`] — NAIF SPK/DAF kernels (e.g. `DE440`).
///
/// See also
/// --------
/// * [`JPLHorizonVersion`] — Maps DE labels to official legacy filenames.
/// * [`NaifVersion`] — Maps DE labels to official SPK filenames.
#[derive(Debug, Clone)]
pub enum EphemFileSource {
    JPLHorizon(JPLHorizonVersion),
    Naif(NaifVersion),
}

/// Parse a source string like `"horizon:DE440"` or `"naif:DE442"`.
///
/// Expected format: `"{source}:{version}"`.
///
/// Errors
/// ------
/// Returns [`OutfitError`] when the string does not match the expected format,
/// or when the version token is unknown for the specified source.
///
/// Examples
/// --------
/// ```
/// use outfit::jpl::download_jpl_file::EphemFileSource;
/// let src: EphemFileSource = "naif:DE440".try_into().unwrap();
/// ```
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
                if let Ok(version) = NaifVersion::from_str(parts[1]) {
                    Ok(EphemFileSource::Naif(version))
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
    /// Return the official base URL for the given backend (only with `jpl-download`).
    ///
    /// Horizons legacy binaries live under:
    /// `https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/`
    ///
    /// NAIF SPK kernels live under:
    /// `https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/`
    ///
    /// See also
    /// --------
    /// * [`EphemFileSource::get_version_url`]
    #[cfg(feature = "jpl-download")]
    fn get_baseurl(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(_) => "https://ssd.jpl.nasa.gov/ftp/eph/planets/Linux/",
            EphemFileSource::Naif(_) => {
                "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/"
            }
        }
    }

    /// Compose the full URL for the concrete version file (only with `jpl-download`).
    ///
    /// This uses the backend‑specific filename returned by the version enums.
    ///
    /// See also
    /// --------
    /// * [`JPLHorizonVersion::get_filename`]
    /// * [`NaifVersion::get_filename`]
    #[cfg(feature = "jpl-download")]
    pub fn get_version_url(&self) -> String {
        let base_url = self.get_baseurl();
        match self {
            EphemFileSource::JPLHorizon(jpl_version) => {
                format!("{}{}", base_url, jpl_version.get_filename())
            }
            EphemFileSource::Naif(version) => format!("{}{}", base_url, version.get_filename()),
        }
    }

    /// Name of the subdirectory under the cache root for this backend.
    ///
    /// - `jpl_horizon/` for legacy binaries
    /// - `naif/` for SPK/DAF files
    fn get_cache_dir(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(_) => "jpl_horizon",
            EphemFileSource::Naif(_) => "naif",
        }
    }

    /// Canonical filename for the requested version (backend‑specific).
    ///
    /// See also
    /// --------
    /// * [`JPLHorizonVersion::to_filename`]
    /// * [`NaifVersion::get_filename`]
    fn filename(&self) -> &str {
        match self {
            EphemFileSource::JPLHorizon(version) => version.to_filename(),
            EphemFileSource::Naif(version) => version.get_filename(),
        }
    }
}

/// Download a (potentially large) file to `path` (feature `jpl-download`).
///
/// Uses `reqwest` to stream the HTTP body in chunks and writes it asynchronously
/// with Tokio's `File` implementation.
///
/// Parameters
/// ----------
/// * `url` — Remote file URL (HTTPS).
/// * `path` — Destination file path.
///
/// Errors
/// ------
/// Returns [`OutfitError`] if the HTTP request, stream, or file I/O fails.
///
/// See also
/// --------
/// * [`EphemFileSource::get_version_url`] — Compose the URL for a versioned file.
#[cfg(feature = "jpl-download")]
pub async fn download_big_file(url: &str, path: &Utf8Path) -> Result<(), OutfitError> {
    let mut file = File::create(path).await?;
    println!("Downloading {url}...");

    let mut stream = reqwest::get(url).await?.bytes_stream();

    while let Some(chunk_result) = stream.next().await {
        let chunk = chunk_result?;
        file.write_all(&chunk).await?;
    }

    file.flush().await?;

    println!("Downloaded {url}");
    Ok(())
}

/// A typed path to a concrete ephemeris file on disk (plus its version).
///
/// This is what the low‑level readers (`horizon` / `naif`) expect to open.
///
/// Variants
/// --------
/// - [`EphemFilePath::JPLHorizon`] — legacy DE binary file.
/// - [`EphemFilePath::Naif`] — NAIF SPK/DAF file.
///
/// See also
/// --------
/// * [`EphemFilePath::get_ephemeris_file`] — Resolve and (optionally) download.
/// * [`EphemFileSource`] — User‑facing selection.
pub enum EphemFilePath {
    JPLHorizon(Utf8PathBuf, JPLHorizonVersion),
    Naif(Utf8PathBuf, NaifVersion),
}

impl std::fmt::Display for EphemFilePath {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            EphemFilePath::JPLHorizon(path, version) => {
                write!(f, "JPL Horizon: {} ({})", path, version.get_filename())
            }
            EphemFilePath::Naif(path, version) => {
                write!(f, "NAIF: {} ({})", path, version.get_filename())
            }
        }
    }
}

impl EphemFilePath {
    /// Resolve the on‑disk path for a given [`EphemFileSource`].
    ///
    /// Behavior:
    /// 1. Build the cache directory under the OS cache root if needed.
    /// 2. Compose the full local filename for the requested version.
    /// 3. If the file **exists**, return its typed path.
    /// 4. If the file is **missing**:
    ///    - with feature `jpl-download`: download it to the cache and return the path,
    ///    - otherwise: return an error.
    ///
    /// Errors
    /// ------
    /// Returns [`OutfitError`] if the base directory cannot be created/resolved,
    /// if the path cannot be represented as UTF‑8, or if the file is missing and
    /// downloads are disabled (or the download fails).
    ///
    /// Examples
    /// --------
    /// ```no_run
    /// use outfit::jpl::download_jpl_file::{EphemFilePath, EphemFileSource};
    /// let source: EphemFileSource = "horizon:DE440".try_into().unwrap();
    /// let handle = EphemFilePath::get_ephemeris_file(&source)?;
    /// println!("Using ephemeris: {}", handle);
    /// # Ok::<(), outfit::outfit_errors::OutfitError>(())
    /// ```
    ///
    /// See also
    /// --------
    /// * [`EphemFileSource`] — Backend + version selector.
    /// * [`download_big_file`] — Async downloader (gated by `jpl-download`).
    pub fn get_ephemeris_file(file_source: &EphemFileSource) -> Result<EphemFilePath, OutfitError> {
        let local_file = EphemFilePath::try_from(file_source.clone())?;

        if local_file.exists() {
            Ok(local_file)
        } else {
            #[cfg(feature = "jpl-download")]
            {
                let url = file_source.get_version_url();

                // Note: since this is a sync API, we spin a Tokio runtime here.
                // If you already run Tokio, consider moving downloads outside
                // or providing an async variant to avoid nested runtimes.
                let rt = tokio::runtime::Runtime::new().expect("Failed to create runtime");
                rt.block_on(async { download_big_file(&url, local_file.path()).await })?;

                Ok(local_file)
            }
            #[cfg(not(feature = "jpl-download"))]
            {
                Err(OutfitError::JPLFileNotFound(local_file.path().to_string()))
            }
        }
    }

    /// Whether the concrete file exists on disk.
    pub fn exists(&self) -> bool {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.exists(),
            EphemFilePath::Naif(path, _) => path.exists(),
        }
    }

    /// Borrow the OS path of the file.
    pub fn path(&self) -> &Utf8Path {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path,
            EphemFilePath::Naif(path, _) => path,
        }
    }

    /// The last path component (file name), if any.
    pub fn file_name(&self) -> Option<&str> {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.file_name(),
            EphemFilePath::Naif(path, _) => path.file_name(),
        }
    }

    /// The file extension (e.g., `bsp`, `bsp.gz`, …), if any.
    pub fn extension(&self) -> Option<&str> {
        match self {
            EphemFilePath::JPLHorizon(path, _) => path.extension(),
            EphemFilePath::Naif(path, _) => path.extension(),
        }
    }
}

/// Build a typed local file path for the requested source, creating the cache dir.
///
/// This `TryFrom` only sets up the directory structure and returns the
/// **intended** file path. It does **not** check for existence beyond
/// creating the parent directories. Use [`EphemFilePath::get_ephemeris_file`]
/// if you need existence + optional download.
///
/// # Errors
/// Returns [`OutfitError`] if the OS cache directory cannot be resolved or
/// converted to a UTF‑8 path, or if creating directories fails.
///
/// # See also
/// * [`directories::BaseDirs`] — Cross‑platform cache dir discovery.
/// * [`EphemFileSource`] — Which backend and version to materialize.
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
            EphemFileSource::Naif(version) => Ok(EphemFilePath::Naif(local_file, version)),
        }
    }
}

#[cfg(test)]
mod jpl_reader_test {
    /// If `jpl-download` is **not** enabled, requesting a missing file must error.
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
