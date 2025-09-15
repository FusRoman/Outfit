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
//! - **Parallel IOD (feature `parallel`)**:
//!   - Batched multi-core execution via **Rayon**,
//!     optional global progress bar when combined with the `progress` feature.
//!
//! ### Planned extensions
//!
//! - **Vaisalä method** for short arcs,
//! - Full support for **hyperbolic trajectories**,
//!
//! ## Workflow at a Glance
//!
//! 1. **Load observations** (MPC/ADES/Parquet).
//! 2. **Initialize** an [`Outfit`](crate::Outfit) environment with JPL ephemerides.
//! 3. **Form triplets** and run the Gauss IOD solver.
//! 4. **Evaluate residuals** (RMS) and select the best candidate.
//! 5. **Propagate** using Keplerian dynamics or convert to equinoctial elements as needed.
//!
//! ## Example (MPC 80-column)
//!
//! ```rust
//! use camino::Utf8Path;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::{Outfit, ErrorModel, IODParams};
//! use outfit::constants::ObjectNumber;
//! use outfit::TrajectorySet;
//! use outfit::prelude::*; // TrajectoryExt, ObservationIOD
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
//!     // Run Gauss IOD.
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
//! For more end-to-end flows, see the [`examples/`](https://github.com/FusRoman/Outfit/tree/main/examples) folder (e.g. `parquet_to_orbit.rs`).
//!
//! ## Example (Parallel batched IOD)
//!
//! This example requires the `parallel` feature (and optionally `progress` for a global progress bar).
//!
//! ```rust
//! use camino::Utf8Path;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::{Outfit, ErrorModel, IODParams};
//! use outfit::constants::ObjectNumber;
//! use outfit::TrajectorySet;
//! use outfit::TrajectoryFit;
//! use outfit::prelude::*; // TrajectoryExt, ObservationIOD
//!
//! # #[cfg(all(feature = "parallel", feature = "jpl-download"))]
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let mut env = Outfit::new("horizon:DE440", ErrorModel::FCCT14)?;
//!     
//!     let test_data = "tests/data/test_from_fink.parquet";
//!     let path_file = Utf8Path::new(test_data);
//!     
//!     let ztf_observer = env.get_observer_from_mpc_code(&"I41".into());
//!     let mut traj_set = TrajectorySet::new_from_parquet(
//!         &mut env,
//!         path_file,
//!         ztf_observer,
//!         0.5,
//!         0.5,
//!         None
//!     )?;
//!
//!     let params = IODParams::builder()
//!         .max_obs_for_triplets(12)
//!         .n_noise_realizations(10)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!     let batch_size = 256; // tune for locality & memory
//!
//!     let results = traj_set.estimate_all_orbits_in_batches_parallel(&env, &mut rng, &params, batch_size);
//!
//!     for (obj, res) in results {
//!         match res {
//!             Ok((gauss, rms)) => {
//!                 println!("{} → RMS = {rms:.4}", obj);
//!                 println!("{gauss}");
//!             }
//!             Err(e) => eprintln!("{} → error: {e}", obj),
//!         }
//!     }
//!
//!     Ok(())
//! }
//!
//! # #[cfg(not(feature = "parallel"))]
//! # fn main() {}
//! ```
//!
//! **Notes**
//! - Batches are processed **in parallel**; each batch is handled **sequentially** to preserve cache locality.
//! - Per-object RNG seeds are **deterministically derived** from a single base seed.
//! - Set `RAYON_NUM_THREADS=N` to cap threads if needed.
//!
//! ## Data Formats
//!
//! - **MPC 80-column** — standard fixed-width astrometry  
//!   <https://minorplanetcenter.net/iau/info/OpticalObs.html>
//! - **ADES XML** — IAU’s Astrometric Data Exchange Standard  
//!   <https://www.iau.org/static/science/scientific_bodies/commissions/f1/ADES-Specification.html>
//! - **Parquet** — columnar format for large batch processing (typical columns: `ra`, `dec`, `jd`, `trajectory_id`)  
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
//! - **`parallel`** *(optional)*  
//!   Enables **multi-core** IOD via **Rayon**, exposing
//!   [`TrajectorySet::estimate_all_orbits_in_batches_parallel`](crate::trajectories::trajectory_fit::TrajectoryFit::estimate_all_orbits_in_batches_parallel).
//!   Combine with `progress` for a thread-safe global progress bar.
//!
//! ```toml
//! [dependencies]
//! outfit = { version = "...", features = ["jpl-download"] }
//! # with progress indicators
//! outfit = { version = "...", features = ["jpl-download", "progress"] }
//! # with multi-core IOD
//! outfit = { version = "...", features = ["parallel"] }
//! # combine as needed
//! outfit = { version = "...", features = ["jpl-download", "progress", "parallel"] }
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

// === Modules (internals). Keep public modules as they are; the facade is built via `pub use` below.

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

/// Trajectory management and file I/O.
pub mod trajectories;

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

// === Public API FACADE =====================================================
// Re-export carefully curated symbols for a simple, stable top-level API.
// Users can import from `outfit::...` without diving into deep module paths.

// Core orchestrator
pub use crate::outfit::Outfit;

// Core data types & units
pub use crate::constants::{ArcSec, Degree, ObjectNumber, MJD};
pub use crate::observers::Observer;
pub use crate::trajectories::TrajectorySet;

// Orbital element representations
pub use crate::orbit_type::{
    cometary_element::CometaryElements, equinoctial_element::EquinoctialElements,
    keplerian_element::KeplerianElements, OrbitalElements,
};

// Error handling and models
pub use crate::error_models::ErrorModel;
pub use crate::outfit_errors::OutfitError;

// IOD (Gauss) key types
pub use crate::initial_orbit_determination::gauss_result::GaussResult;
pub use crate::initial_orbit_determination::IODParams;

// Frequently-used extension traits (ergonomic entry points) and key types
pub use crate::observations::observations_ext::ObservationIOD;
pub use crate::trajectories::trajectory_file::TrajectoryFile;
pub use crate::trajectories::trajectory_fit::FullOrbitResult;
pub use crate::trajectories::trajectory_fit::TrajectoryFit;

// Selected constants that are widely useful
pub use crate::constants::{
    AU, GAUSS_GRAV, RADEG, RADH, RADSEC, SECONDS_PER_DAY, T2000, VLIGHT_AU,
};

// JPL ephemeris enum for runtime inspection (optional but convenient)
pub use crate::jpl_ephem::JPLEphem;

// A convenient crate-wide Result alias.
pub type Result<T> = core::result::Result<T, OutfitError>;

/// Prelude with common imports for quick-start users.
///
/// Bring the most frequently used types and traits into scope in one line:
/// ```rust
/// use outfit::prelude::*;
/// ```
pub mod prelude {
    pub use crate::{
        ArcSec, Degree, ErrorModel, FullOrbitResult, GaussResult, IODParams, JPLEphem,
        ObjectNumber, ObservationIOD, Observer, Outfit, OutfitError, TrajectoryFile, TrajectorySet,
        MJD,
    };
    // Optionally include widely-used constants:
    pub use crate::{AU, GAUSS_GRAV, RADEG, RADH, RADSEC, SECONDS_PER_DAY, T2000, VLIGHT_AU};
}

// === Tests support ==========================================================
#[cfg(all(test, feature = "jpl-download"))]
pub(crate) mod unit_test_global {
    use std::sync::LazyLock;

    use camino::Utf8Path;

    use crate::{
        error_models::ErrorModel,
        jpl_ephem::{horizon::horizon_data::HorizonData, naif::naif_data::NaifData},
        outfit::Outfit,
        trajectories::trajectory_file::TrajectoryFile,
        trajectories::TrajectorySet,
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
