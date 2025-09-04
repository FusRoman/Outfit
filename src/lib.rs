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
//! and propagate them through time.
//!
//! ## Why Outfit?
//!
//! - **Safe Rust**: strong typing, memory safety, and explicit error handling.
//! - **Modular**: interchangeable ephemerides, error models, and IOD algorithms.
//! - **Scientific grounding**: consistent with established astronomical standards.
//! - **Scalable workflows**: designed for survey-scale batch processing.
//!
//! ## Features
//!
//! - **Initial Orbit Determination (IOD)**:
//!   - Classical **Gauss method** for three observations.
//! - **Orbital elements**:
//!   - Classical **Keplerian elements**,
//!   - **Equinoctial elements** with conversions and two-body solver.
//! - **Reference frames & preprocessing**:
//!   - Precession, nutation (IAU 1980), aberration, and light-time correction,
//!   - Ecliptic ↔ equatorial conversions, RA/DEC parsing, time systems.
//! - **Ephemerides**:
//!   - Built-in support for **JPL DE440** (NAIF/SPICE kernels and Horizons format).
//! - **Observation I/O**:
//!   - **MPC 80-column**, **ADES XML**, and **Parquet** (batch) readers,
//!   - Optimized batched loading, automatic per-observer error assignments.
//! - **Observer management**:
//!   - Build from **MPC observatory code** or custom geodetic coordinates.
//! - **Residuals & quality metrics**:
//!   - RMS computation of normalized astrometric residuals, filtering utilities.
//! - **Examples & benches**:
//!   - End-to-end examples in `examples/`, Criterion benchmarks for IOD.
//!
//! ### Planned extensions
//!
//! - **Vaisalä method** for short arcs,
//! - Full support for **hyperbolic trajectories**,
//!
//! ## Workflow at a Glance
//!
//! 1. **Load observations** (MPC/ADES/Parquet).
//! 2. **Initialize** an [`Outfit`](crate::outfit::Outfit) environment with JPL ephemerides.
//! 3. **Form triplets** and run the Gauss IOD solver.
//! 4. **Evaluate residuals** (RMS) and select the best candidate.
//! 5. **Propagate** using Keplerian dynamics or convert to equinoctial elements as needed.
//!
//! ## Example (MPC 80-column)
//!
//! ```rust
//! use camino::Utf8Path;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::outfit::Outfit;
//! use outfit::constants::{TrajectorySet, ObjectNumber};
//! use outfit::error_models::ErrorModel;
//! use outfit::initial_orbit_determination::IODParams;
//! use outfit::observations::trajectory_ext::TrajectoryExt;
//! use outfit::observations::observations_ext::ObservationIOD;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//!
//!     // Load observations (80-col MPC).
//!     let mut trajs = TrajectorySet::new_from_80col(
//!         &mut env,
//!         Utf8Path::new("tests/data/2015AB.obs"),
//!     );
//!     trajs.add_from_80col(&mut env, Utf8Path::new("tests/data/8467.obs"));
//!
//!     // Select an object.
//!     let obj = ObjectNumber::String("K09R05F".into());
//!     let obs = trajs.get_mut(&obj).expect("object not found");
//!
//!     // Configure IOD.
//!     let params = IODParams::builder()
//!         .n_noise_realizations(10)
//!         .noise_scale(1.1)
//!         .max_obs_for_triplets(obs.len())
//!         .max_triplets(30)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!
////!     // Run Gauss IOD.
//!     let (best_orbit, best_rms) = obs.estimate_best_orbit(
//!         &mut env,
//!         &ErrorModel::FCCT14,
//!         &mut rng,
//!         &params,
//!     )?;
//!
//!     println!("Best orbit: {:?}", best_orbit);
//!     println!("RMS: {}", best_rms);
//!     Ok(())
//! }
//! ```
//!
//! For more end-to-end flows, see the `examples/` folder (e.g. `parquet_to_orbit.rs`).
//!
//! ## Data Formats
//!
//! - **MPC 80-column** — standard fixed-width astrometry  
//!   <https://minorplanetcenter.net/iau/info/OpticalObs.html>
//! - **ADES XML** — IAU’s Astrometric Data Exchange Standard  
//!   <https://www.iau.org/static/science/scientific_bodies/commissions/f1/ADES-Specification.html>
//! - **Parquet** — columnar format for large batch processing (typical columns: `ra`, `dec`, `jd|mjd`, `trajectory_id`)  
//!   <https://parquet.apache.org/docs/>
//!
//! ## Cargo Features
//!
//! - **`jpl-download`**  
//!   Automatically download JPL ephemerides (NAIF/Horizons) into a local cache.
//!
//!   Cache layout (Linux):
//!   ```text
//!   ~/.cache/outfit_cache/jpl_ephem/
//!   ├── jpl_horizon/
//!   │   └── DE440.bsp
//!   └── naif/
//!       └── de440.bsp
//!   ```
//!   This feature enables some integration tests and pulls `reqwest`, `tokio`, `tokio-stream`.
//!
//! - **`progress`** *(optional)*  
//!   Enables lightweight progress bars (via `indicatif`) and loop timing utilities for long-running jobs.
//!
//! ```toml
//! [dependencies]
//! outfit = { version = "...", features = ["jpl-download"] }
//! # or, with progress indicators
//! outfit = { version = "...", features = ["jpl-download", "progress"] }
//! ```
//!
//! ## Error Handling
//!
//! Outfit uses [`thiserror`] for descriptive error types and idiomatic `Result<T, E>` returns,
//! with transparent conversions via `#[from]` where appropriate.
//!
//! ## MSRV & Platforms
//!
//! - Tested on **stable Rust** (see `Cargo.toml` and CI for the current MSRV policy).
//! - Linux/macOS are primary targets; Windows is expected to work but may have different I/O performance.
//!
//! ## Scientific Background
//!
//! - Milani, A. & Gronchi, G. F. (2010), *Theory of Orbit Determination*, Cambridge University Press.
//! - The open-source [OrbFit](http://adams.dm.unipi.it/orbfit/) software suite.
//!
//! ## License
//!
//! Distributed under the **CeCILL-C** license. See
//! [`LICENSE`](https://github.com/FusRoman/Outfit/blob/main/LICENSE).
//!
//! ## Authors
//!
//! Developed by **FusRoman** and contributors.
//!
//! ## Useful Modules (See also)
//! ------------
//! - [`initial_orbit_determination`] — Gauss IOD and helpers.
//! - [`observations`] — Readers (MPC/ADES/Parquet) and trajectory utilities.
//! - [`observers`] — Observatory registry and geodetic utilities.
//! - [`jpl_ephem`] — Ephemerides backends (Horizons/NAIF).
//! - [`ref_system`] — Reference frame conversions and rotations.

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
