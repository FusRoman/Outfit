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
//! **Outfit** is a modern **Rust library for orbit determination and propagation**
//! of small Solar System bodies.  
//!
//! It provides a **memory-safe, modular, and extensible** implementation of
//! algorithms that transform raw astrometric observations into orbital solutions,
//! and propagate them reliably through time.  
//!
//! ## Why Outfit?
//!
//! Outfit is designed to provide a **modern, reliable, and extensible toolkit**
//! for orbit determination and propagation.  
//!
//! Its core principles are:
//!
//! - A **safe Rust implementation**, leveraging strong typing and memory safety,
//! - A **modular design**, allowing interchangeable ephemerides, error models, and IOD algorithms,
//! - A **scientifically validated foundation**, consistent with established astronomical standards,
//! - Built-in support for **scalable and reproducible workflows**, suitable for survey-scale processing
//!   and long-term research pipelines.
//!
//! ## Features
//!
//! Outfit currently provides:
//!
//! - **Initial Orbit Determination (IOD):**
//!   - Classical **Gauss method** for 3 observations.
//! - **Orbital elements:**
//!   - Classical Keplerian elements,
//!   - Non-singular equinoctial elements (planned).
//! - **Reference frame transformations:**
//!   - Precession, nutation, aberration, and light-time correction,
//!   - Conversions between ecliptic and equatorial frames.
//! - **Ephemerides handling:**
//!   - Built-in support for **JPL DE440** (both NAIF kernels and Horizons-format files).
//! - **Astrometric preprocessing:**
//!   - RA/DEC parsing, time system conversions, and observatory geodetic coordinates.
//! - **Residual filtering framework:**
//!   - Extensible error models, preliminary RMS computation.
//!
//! ### Planned extensions
//!
//! - **Vaisala method** for short arcs,  
//! - Advanced RMS residual filtering and iterative orbit refinement,  
//! - Full support for **hyperbolic trajectories** (e.g., ‘Oumuamua),  
//! - Batch and parallelized pipelines for large surveys.  
//!
//! ## Workflow at a Glance
//!
//! 1. **Load observational data** (RA, DEC, times, observatory position).
//! 2. **Initialize** an [`Outfit`](crate::outfit::Outfit) environment with JPL ephemerides.
//! 3. **Form triplets** of observations and run the Gauss IOD solver.
//! 4. **Evaluate residuals** and optionally refine the orbit.
//! 5. **Propagate** forward/backward in time using Keplerian dynamics.
//!
//! ## Example: Orbit Determination from MPC 80-column Files
//!
//! ```rust, ignore
//! use camino::Utf8Path;
//! use rand::{SeedableRng, rngs::StdRng};
//! use outfit::outfit::Outfit;
//! use outfit::constants::{TrajectorySet, ObjectNumber};
//! use outfit::error_models::ErrorModel;
//! use outfit::initial_orbit_determination::IODParams;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Initialize environment with DE440 ephemerides.
//!     let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//!
//!     // Load astrometric observations.
//!     let mut traj_set = TrajectorySet::new_from_80col(
//!         &mut env,
//!         Utf8Path::new("tests/data/2015AB.obs"),
//!     );
//!     traj_set.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
//!
//!     // Select an object by its MPC designation.
//!     let object_id = ObjectNumber::String("K09R05F".into());
//!     let obs = traj_set.get_mut(&object_id).unwrap();
//!
//!     // Configure IOD parameters.
//!     let params = IODParams::builder()
//!         .n_noise_realizations(10)
//!         .noise_scale(1.1)
//!         .max_obs_for_triplets(obs.len())
//!         .max_triplets(30)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!
//!     // Estimate the preliminary orbit.
//!     let (best_orbit, best_rms) = obs.estimate_best_orbit(
//!         &mut env,
//!         &ErrorModel::FCCT14,
//!         &mut rng,
//!         &params,
//!     )?;
//!
//!     println!("Best orbit: {:?}", best_orbit);
//!     println!("RMS residuals: {}", best_rms);
//!
//!     Ok(())
//! }
//! ```
//!
//! ## Cargo Features
//!
//! Outfit supports optional extensions via Cargo features:
//!
//! - **`jpl-download`**  
//!   Automatically download JPL ephemerides (NAIF/Horizons) into a local cache.  
//!   If not enabled, files must be placed manually.  
//!
//!   Cache directory structure (Linux example):
//!   ```text
//!   ~/.cache/outfit_cache/jpl_ephem/
//!   ├── jpl_horizon/
//!   │   └── DE440.bsp
//!   └── naif/
//!       └── de440.bsp
//!   ```
//!
//!   This feature enables some integration tests and pulls in `reqwest`, `tokio`, and `tokio-stream`.
//!
//! ```toml
//! [dependencies]
//! outfit = { version = "...", features = ["jpl-download"] }
//! ```
//!
//! ## Scientific Background
//!
//! Outfit’s algorithms follow the formalism of:
//!
//! - Milani, A. & Gronchi, G. F. (2010), *Theory of Orbit Determination*,
//!   Cambridge University Press.
//! - The [OrbFit](http://adams.dm.unipi.it/orbfit/) software suite (open-source).
//!
//! ## License
//!
//! Distributed under the **CeCILL-C** license. See [LICENSE](https://github.com/FusRoman/Outfit/blob/main/LICENSE).
//!
//! ## Authors
//!
//! Developed by **FusRoman** and contributors.
//!
//! ## Links
//!
//! - [Documentation](https://docs.rs/outfit)
//! - [Source code](https://github.com/FusRoman/Outfit)

/// Constants and astronomical unit conversions.
pub mod constants;

/// Coordinate conversions (RA/DEC, spherical ↔ cartesian, etc.).
pub mod conversion;

/// Earth orientation parameters and related corrections (nutation, precession).
pub mod earth_orientation;

/// Environment state: ephemerides, dynamical models, configuration.
pub mod env_state;

/// Error models used for weighting astrometric residuals.
pub mod error_models;

/// Initial Orbit Determination algorithms (Gauss method).
pub mod initial_orbit_determination;

/// JPL ephemerides management (Horizon/NAIF kernels).
pub mod jpl_ephem;

/// Keplerian solver for propagation.
pub mod kepler;

/// Orbital types and conversions between them.
pub mod orbit_type;

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
