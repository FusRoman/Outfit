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
//! - **Observer management**:
//!   - Build from **MPC observatory code** or custom geodetic coordinates.
//! - **Residuals & quality metrics**:
//!   - RMS computation of normalized astrometric residuals, filtering utilities.
//! - **Examples**:
//!   - End-to-end examples in `examples/`.
//! - **Parallel IOD (feature `parallel`)**:
//!   - Batched multi-core execution via **Rayon**,
//!     optional global progress bar when combined with the `progress` feature.
//!
//! ### Planned extensions
//!
//! - **Vaisalä method** for short arcs,
//! - Full support for **hyperbolic trajectories**.
//!
//! ## Workflow at a Glance
//!
//! 1. **Load observations** into an [`ObsDataset`](photom::observation_dataset::ObsDataset)
//!    (MPC 80-column, ADES XML, or Parquet).
//! 2. **Load JPL ephemerides** via [`JPLEphem`] and prepare a UT1 provider.
//! 3. **Run IOD** by calling [`FitIOD::fit_iod`](obs_dataset::FitIOD::fit_iod) on the dataset.
//! 4. **Evaluate residuals** (RMS) using the returned [`IODRMS`](obs_dataset::IODRMS).
//! 5. **Propagate** using Keplerian dynamics or convert to equinoctial elements as needed.
//!
//! ## Example (single-trajectory IOD from MPC 80-column)
//!
//! ```rust,no_run
//! use photom::observation_dataset::ObsDataset;
//! use photom::observer::error_model::ObsErrorModel;
//! use hifitime::ut1::Ut1Provider;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::obs_dataset::FitIOD;
//! use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem};
//! use outfit::IODParams;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     // Load observations from an MPC 80-column file.
//!     let dataset = ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")?;
//!
//!     // Load JPL DE440 ephemeris (Horizons format).
//!     let jpl_source: EphemFileSource = "horizon:DE440".try_into()?;
//!     let jpl = JPLEphem::new(&jpl_source)?;
//!
//!     // Obtain a UT1 provider for Earth orientation corrections.
//!     let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long")?;
//!
//!     // Configure IOD.
//!     let params = IODParams::builder()
//!         .n_noise_realizations(10)
//!         .max_obs_for_triplets(50)
//!         .max_triplets(30)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!
//!     // Run Gauss IOD for a single trajectory.
//!     let (best_orbit, best_rms) = dataset.fit_iod(
//!         "K09R05F",
//!         &jpl,
//!         &ut1,
//!         &params,
//!         ObsErrorModel::FCCT14,
//!         &mut rng,
//!     )?;
//!
//!     println!("Best orbit: {best_orbit}");
//!     println!("RMS: {best_rms:.6}");
//!     Ok(())
//! }
//! ```
//!
//! For more end-to-end flows, see the [`examples/`](https://github.com/FusRoman/Outfit/tree/main/examples) folder.
//!
//! ## Example (batch IOD — all trajectories at once)
//!
//! ```rust,no_run
//! use photom::observation_dataset::ObsDataset;
//! use photom::observer::error_model::ObsErrorModel;
//! use hifitime::ut1::Ut1Provider;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::obs_dataset::{FitIOD, FullOrbitResult};
//! use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem};
//! use outfit::IODParams;
//!
//! fn main() -> Result<(), Box<dyn std::error::Error>> {
//!     let dataset = ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")?;
//!
//!     let jpl_source: EphemFileSource = "horizon:DE440".try_into()?;
//!     let jpl = JPLEphem::new(&jpl_source)?;
//!     let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long")?;
//!
//!     let params = IODParams::builder()
//!         .n_noise_realizations(10)
//!         .max_obs_for_triplets(50)
//!         .build()?;
//!
//!     let mut rng = StdRng::seed_from_u64(42);
//!
//!     // Run Gauss IOD for every trajectory in the dataset.
//!     let results: FullOrbitResult = dataset.fit_full_iod(
//!         &jpl,
//!         &ut1,
//!         &params,
//!         ObsErrorModel::FCCT14,
//!         &mut rng,
//!     )?;
//!
//!     for (traj_id, res) in &results {
//!         match res {
//!             Ok((gauss, rms)) => println!("{traj_id} → RMS = {rms:.4}\n{gauss}"),
//!             Err(e) => eprintln!("{traj_id} → error: {e}"),
//!         }
//!     }
//!     Ok(())
//! }
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
//! ## See also
//!
//! - [`initial_orbit_determination`] — Gauss IOD algorithm, triplet generation, and IOD parameters.
//! - [`obs_dataset`] — [`FitIOD`](obs_dataset::FitIOD) trait and batch/single-trajectory IOD entry points.
//! - [`jpl_ephem`] — Ephemerides backends (Horizons/NAIF DE440).
//! - [`orbit_type`] — Orbital element representations (Keplerian, Equinoctial, Cometary).
//! - [`ref_system`] — Reference frame transformations.
//! - [`constants`] — Physical constants and unit conversions.
//! - [`conversion`] — RA/DEC parsing and coordinate utilities.
//! - [`earth_orientation`] — Precession, nutation, and obliquity models.
//! - [`kepler`] — Universal Kepler propagator and Lagrange f–g solver.
//! - [`orb_elem`] — State-vector to orbital-elements conversion.
//! - [`outfit_errors`] — Unified error enum for the whole crate.
//! - [`time`] — Time scale conversions and sidereal time.
//! - [`cache`] — Precomputed observer position cache.
//! - [`observation_ephemeris`] — Apparent position computation and ephemeris residuals.
//! - [`observer_extension`] — Geocentric and heliocentric observer position routines.
//! - [`trajectory`](crate::obs_dataset) — Core IOD pipeline trait over observation slices.

// === Modules (internals). Keep public modules as they are; the facade is built via `pub use` below.

/// Constants and astronomical unit conversions.
pub mod constants;

/// Coordinate conversions (RA/DEC, spherical ↔ cartesian, etc.).
pub mod conversion;

/// Earth orientation parameters and related corrections (nutation, precession).
pub mod earth_orientation;

/// Initial Orbit Determination algorithms (Gauss method).
pub mod initial_orbit_determination;

/// JPL ephemerides management (Horizon/NAIF kernels).
pub mod jpl_ephem;

/// Keplerian solver for propagation.
pub mod kepler;

/// Orbital types and conversions between them.
pub mod orbit_type;

/// Orbital elements utilities (conversion, normalization).
pub mod orb_elem;

/// Errors returned by Outfit operations.
pub mod outfit_errors;

/// Reference frame transformations.
pub mod ref_system;

/// Time management and conversions (UTC, TDB, TT).
pub mod time;

/// Precomputed observer position cache used throughout the fitting pipeline.
pub mod cache;

/// [`FitIOD`] trait and batch/single-trajectory IOD entry points on [`photom::observation_dataset::ObsDataset`].
pub mod obs_dataset;

/// Apparent position computation and astrometric residuals for individual observations.
pub mod observation_ephemeris;

/// Ground-observer geometry: body-fixed and heliocentric position routines.
pub mod observer_extension;

/// Core IOD pipeline trait over sorted observation slices.
pub(crate) mod trajectory;

// === Public API FACADE =====================================================
// Re-export carefully curated symbols for a simple, stable top-level API.
// Users can import from `outfit::...` without diving into deep module paths.

// Orbital element representations
pub use crate::orbit_type::{
    cometary_element::CometaryElements, equinoctial_element::EquinoctialElements,
    keplerian_element::KeplerianElements, OrbitalElements,
};

pub use crate::outfit_errors::OutfitError;

// IOD (Gauss) key types
pub use crate::initial_orbit_determination::gauss_result::GaussResult;
pub use crate::initial_orbit_determination::IODParams;

// Selected constants that are widely useful
pub use crate::constants::{
    AU, GAUSS_GRAV, RADEG, RADH, RADSEC, SECONDS_PER_DAY, T2000, VLIGHT_AU,
};

// JPL ephemeris enum for runtime inspection (optional but convenient)
pub use crate::jpl_ephem::JPLEphem;

// IOD entry points and result types
pub use crate::obs_dataset::{FitIOD, FullOrbitResult, IODRMS};

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
        FitIOD, FullOrbitResult, GaussResult, IODParams, JPLEphem, OutfitError, IODRMS,
    };
    // Optionally include widely-used constants:
    pub use crate::{AU, GAUSS_GRAV, RADEG, RADH, RADSEC, SECONDS_PER_DAY, T2000, VLIGHT_AU};
}

// === Tests support ==========================================================
#[cfg(test)]
pub(crate) mod test_fixture {
    use std::sync::LazyLock;

    use hifitime::ut1::Ut1Provider;
    use photom::observation_dataset::ObsDataset;

    use crate::{jpl_ephem::download_jpl_file::EphemFileSource, JPLEphem};

    pub(crate) static UT1_PROVIDER: LazyLock<Ut1Provider> = LazyLock::new(|| {
        Ut1Provider::download_from_jpl("latest_eop2.long")
            .expect("Download of the JPL short time scale UT1 data failed")
    });

    pub(crate) static JPL_EPHEM_HORIZON: LazyLock<JPLEphem> = LazyLock::new(|| {
        let jpl_file: EphemFileSource = "horizon:DE440"
            .try_into()
            .expect("Failed to parse JPL ephemeris source");
        JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Horizon")
    });

    pub(crate) static JPL_EPHEM_NAIF: LazyLock<JPLEphem> = LazyLock::new(|| {
        let jpl_file: EphemFileSource = "naif:DE440"
            .try_into()
            .expect("Failed to parse JPL ephemeris source");
        JPLEphem::new(&jpl_file).expect("Failed to load JPL ephemeris from Naif")
    });

    pub(crate) static DATASET_2015AB: LazyLock<ObsDataset> = LazyLock::new(|| {
        ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")
            .expect("Failed to load test dataset 2015AB")
    });
}
