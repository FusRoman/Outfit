//! [![github]](https://github.com/FusRoman/Outfit)&ensp;[![crates-io]](https://crates.io/crates/outfit)&ensp;[![docs-rs]](crate)
//!
//! [github]: https://img.shields.io/badge/github-8da0cb?style=for-the-badge&labelColor=555555&logo=github
//! [crates-io]: https://img.shields.io/badge/crates.io-fc8d62?style=for-the-badge&labelColor=555555&logo=rust
//! [docs-rs]: https://img.shields.io/badge/docs.rs-66c2a5?style=for-the-badge&labelColor=555555&logo=docs.rs
//!
//! <br>
//!
//! # Outfit
//!
//! **Outfit** is a modern Rust library for orbit determination and propagation
//! of small Solar System bodies.  
//! It is inspired by and validated against the [OrbFit](http://adams.dm.unipi.it/orbfit/)
//! Fortran package, but rewritten from scratch with a modern, safe architecture.
//!
//! ## Overview
//!
//! Outfit provides algorithms and data structures to transform raw astrometric
//! observations into preliminary orbits, refine those orbits, and propagate them
//! through time.  
//!
//! Key features currently include:
//!
//! - **Initial Orbit Determination (IOD):**
//!   - Classical **Gauss method** for 3 observations
//! - **Keplerian orbital elements** output
//! - **Reference frame transformations:**
//!   - Precession, nutation, aberration, light-time correction
//!   - Ecliptic/equatorial conversions
//! - **Ephemerides handling:**
//!   - Built-in support for JPL DE440 data, both NAIF kernels and Horizons files
//! - **Astrometric preprocessing:**
//!   - Parsing RA/DEC, time conversions, observatory geodetic location
//! - **Extensible framework** for error models and RMS residual filtering
//!
//! ### Planned features
//! - **Vaisala method** for short observational arcs
//! - Advanced RMS residual filtering and orbit refinement
//! - Full hyperbolic orbit support
//!
//! ## Typical Workflow
//!
//! 1. **Load your observational data** (RA, DEC, time, observatory position).
//! 2. **Initialize** an [`Outfit`](crate::outfit::Outfit) environment with a JPL ephemeris.
//! 3. **Generate triplets** and use the Gauss method for preliminary orbit determination.
//! 4. **Optionally refine** the orbit using residual filtering and RMS computation.
//! 5. **Propagate** the resulting orbit forward or backward in time.
//!
//! ## Example: Orbit Determination from MPC 80-column files
//!
//! ```rust, ignore
//! use camino::Utf8Path;
//! use rand::SeedableRng;
//! use rand::rngs::StdRng;
//! use outfit::outfit::Outfit;
//! use outfit::constants::{TrajectorySet, ObjectNumber};
//! use outfit::error_models::ErrorModel;
//! use outfit::initial_orbit_determination::{IODParams, gauss_result::GaussResult};
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Initialize environment with JPL Horizons ephemerides (DE440).
//!     let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//!
//!     // Load astrometric observations from MPC 80-column formatted files.
//!     let mut traj_set = TrajectorySet::new_from_80col(
//!         &mut env,
//!         Utf8Path::new("tests/data/2015AB.obs"),
//!     );
//!     traj_set.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
//!
//!     // Select the object of interest by its MPC designation.
//!     let object_id = ObjectNumber::String("K09R05F".into());
//!
//!     // Configure IOD parameters (triplets, noise model, etc.).
//!     let obs = traj_set.get_mut(&object_id).unwrap();
//!     let params = IODParams::builder()
//!         .n_noise_realizations(10)
//!         .noise_scale(1.1)
//!         .max_obs_for_triplets(obs.len())
//!         .max_triplets(30)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!
//!     // Estimate the best preliminary orbit using the Gauss method.
//!     let (best_orbit, best_rms) = obs.estimate_best_orbit(&mut env, &ErrorModel::FCCT14, &mut rng, &params)?;
//!
//!     println!("Best orbit: {:?}", best_orbit);
//!     println!("RMS residuals: {}", best_rms);
//!
//!     Ok(())
//! }
//! ```
//!
//! This example demonstrates how to:
//! - Set up an environment (`Outfit`) with a JPL ephemeris.
//! - Load observations from MPC 80-column files.
//! - Run the Gauss IOD algorithm to obtain a preliminary orbit and RMS residuals.
//!
//! ## Optional features
//!
//! The following Cargo feature can be enabled to extend Outfit's functionality:
//!
//! - **`jpl-download`**: Enables automatic downloading of JPL ephemeris files (NAIF/Horizons) when needed, including for unit and integration tests. If this feature is not enabled, you must manually place the required ephemeris file in the cache directory on your computer.
//!
//!   - **On Linux:** `~/.cache/outfit_cache/jpl_ephem/`
//!   - **On Windows:** The standard user cache directory (e.g., `C:\Users\<username>\AppData\Local\outfit_cache\jpl_ephem\`)
//!
//!   The cache directory contains two subfolders:
//!   - `jpl_horizon/` for Horizon-format files (e.g., `DE440.bsp`)
//!   - `naif/` for NAIF/SPICE files (e.g., `de440.bsp`)
//!
//!   Example structure:
//!   ```text
//!   ~/.cache/outfit_cache/jpl_ephem/
//!   ├── jpl_horizon/
//!   │   └── DE440.bsp
//!   └── naif/
//!       └── de440.bsp
//!   ```
//!
//!   For more details about the cache structure and ephemeris file handling, see [`src/jpl_ephem/download_jpl_file.rs`](src/jpl_ephem/download_jpl_file.rs).
//!
//!   This feature pulls in the `reqwest`, `tokio`, and `tokio-stream` dependencies and is required for some integration tests and for the `unit_test_global` module.
//!
//! Enable the feature in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! outfit = { version = "...", features = ["jpl-download"] }
//! ```
//!
//! ## Scientific Background
//!
//! The algorithms implemented are based on:
//!
//! - Milani, A. & Gronchi, G. F. (2010), *Theory of Orbit Determination*,
//!   Cambridge University Press.
//! - The open-source [OrbFit](http://adams.dm.unipi.it/orbfit/) software suite.
//!
//! ## License
//!
//! Outfit is distributed under the MIT license. See the [LICENSE](https://github.com/FusRoman/Outfit)
//! file for details.
//!
//! ## Authors
//!
//! Roman and contributors
//!
//! ## Links
//!
//! - [Documentation](https://docs.rs/outfit)
//! - [Source code](https://github.com/FusRoman/Outfit)
//!
/// Constants and astronomical unit conversions.
pub mod constants;

/// Coordinate conversions (RA/DEC, spherical ↔ cartesian, etc.).
pub mod conversion;

/// Earth orientation parameters and related corrections (nutation, precession).
pub mod earth_orientation;

/// Environment state: ephemerides, dynamical models, configuration.
pub mod env_state;

/// Equinoctial orbital elements and related conversions.
pub mod equinoctial_element;

/// Error models used for weighting astrometric residuals.
pub mod error_models;

/// Initial Orbit Determination algorithms (Gauss method).
pub mod initial_orbit_determination;

/// JPL ephemerides management (Horizon/NAIF kernels).
pub mod jpl_ephem;

/// Keplerian solver for propagation.
pub mod kepler;

/// Classical Keplerian elements structure and utilities.
pub mod keplerian_element;

/// Observation handling (RA, DEC, times).
pub mod observations;

/// Observers and observatory positions.
pub mod observers;

/// Orbital elements utilities (conversion, normalization).
pub mod orb_elem;

/// Main Outfit struct: central orchestrator for orbit determination.
pub mod outfit;

/// Errors returned by Outfit operations.
pub mod outfit_errors;

/// Reference frame transformations.
pub mod ref_system;

/// Time management and conversions (UTC, TDB, TT).
pub mod time;

#[cfg(all(test, feature = "jpl-download"))]
pub(crate) mod unit_test_global {
    use std::sync::LazyLock;

    use camino::Utf8Path;

    use crate::{
        constants::TrajectorySet,
        error_models::ErrorModel,
        jpl_ephem::{horizon::horizon_data::HorizonData, naif::naif_data::NaifData},
        observations::trajectory_ext::TrajectoryExt,
        outfit::Outfit,
    };

    pub(crate) static OUTFIT_NAIF_TEST: LazyLock<Outfit> =
        LazyLock::new(|| Outfit::new("naif:DE440", ErrorModel::FCCT14).unwrap());

    pub(crate) static OUTFIT_HORIZON_TEST: LazyLock<(Outfit, TrajectorySet)> =
        LazyLock::new(|| {
            let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14).unwrap();

            let path_file = Utf8Path::new("tests/data/2015AB.obs");
            let traj_set = TrajectorySet::new_from_80col(&mut env, path_file);
            (env, traj_set)
        });

    pub(crate) static JPL_EPHEM_HORIZON: LazyLock<&HorizonData> = LazyLock::new(|| {
        let jpl_ephem = OUTFIT_HORIZON_TEST.0.get_jpl_ephem().unwrap();
        match jpl_ephem {
            crate::jpl_ephem::JPLEphem::HorizonFile(horizon_data) => horizon_data,
            _ => panic!("JPL ephemeris is not a Horizon file"),
        }
    });

    pub(crate) static JPL_EPHEM_NAIF: LazyLock<&NaifData> = LazyLock::new(|| {
        let jpl_ephem = OUTFIT_NAIF_TEST.get_jpl_ephem().unwrap();
        match jpl_ephem {
            crate::jpl_ephem::JPLEphem::NaifFile(naif_data) => naif_data,
            _ => panic!("JPL ephemeris is not a Naif file"),
        }
    });
}
