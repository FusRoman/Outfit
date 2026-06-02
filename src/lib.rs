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
//!   - **Equinoctial elements** with conversions and two-body solver,
//!   - **Cometary elements** for parabolic and hyperbolic trajectories.
//! - **Uncertainty propagation**:
//!   - Full **covariance matrix propagation** between orbital element representations,
//!   - Standard deviations and correlations tracked via **Jacobian transformations**,
//!   - Rigorous linear uncertainty propagation for orbit conversions.
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
//! - **Parallel IOD (feature `parallel`)**
//!
//! ### Planned extensions
//!
//! - **Vaisalä method** for short arcs,
//! - Full support for **hyperbolic trajectories**.
//!
//! ## Why Track Uncertainties?
//!
//! Orbital element uncertainties quantify the confidence in derived orbits and are essential for:
//!
//! - **Orbit quality assessment** — distinguishing reliable solutions from poorly-constrained fits.
//! - **Prediction uncertainty** — forecasting position errors for future ephemerides.
//! - **Data association** — deciding whether observations belong to the same object.
//! - **Differential correction** — weighting observations and detecting outliers.
//! - **Mission planning** — evaluating pointing and timing requirements for follow-up observations.
//!
//! Outfit propagates uncertainties rigorously through coordinate transformations using **Jacobian
//! matrices** and **linear covariance propagation**, ensuring that correlations between orbital
//! parameters are preserved when converting between Keplerian, equinoctial, and cometary representations.
//!
//! For mathematical details and usage examples, see the [`orbit_type`] module documentation.
//!
//! ## Workflow at a Glance
//!
//! 1. **Load observations** into an [`ObsDataset`](photom::observation_dataset::ObsDataset)
//!    (MPC 80-column, ADES XML, or Parquet).
//! 2. **Load JPL ephemerides** via [`JPLEphem`] and prepare a UT1 provider.
//! 3. **Run IOD** by calling [`FitIOD::fit_iod`](crate::FitIOD::fit_iod) on the dataset.
//! 4. **Evaluate residuals** (RMS) using the returned [`IODRMS`](crate::IODRMS).
//! 5. **Propagate** using Keplerian dynamics or convert to equinoctial elements as needed.
//!
//! ## Example (single-trajectory IOD from MPC 80-column)
//!
//! ```rust,no_run
//! use photom::observation_dataset::ObsDataset;
//! use photom::observer::error_model::ObsErrorModel;
//! use hifitime::ut1::Ut1Provider;
//! use rand::{rngs::StdRng, SeedableRng};
//! use outfit::FitIOD;
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
//!     let fit_result = dataset.fit_iod(
//!         "K09R05F",
//!         &jpl,
//!         &ut1,
//!         &params,
//!         ObsErrorModel::FCCT14,
//!         &mut rng,
//!     )?;
//!
//!     println!("Best orbit: {}", fit_result.orbital_elements());
//!     println!("Quality: {:.6}", fit_result.orbit_quality());
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
//! use outfit::{FitIOD, FullOrbitResult};
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
//!             Ok(fit) => println!("{traj_id} → quality = {:.4}\n{}", fit.orbit_quality(), fit.orbital_elements()),
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
//! - [`initial_orbit_determination`] — Gauss IOD algorithm, triplet generation, IOD parameters, public api to perform iod on an [`ObsDataset`](photom::observation_dataset::ObsDataset).
//! - [`jpl_ephem`] — Ephemerides backends (Horizons/NAIF DE440).
//! - [`orbit_type`] — Orbital element representations (Keplerian, Equinoctial, Cometary) with full covariance propagation and uncertainty tracking for conversions between representations.
//! - [`ref_system`] — Reference frame transformations.
//! - [`constants`] — Physical constants and unit conversions.
//! - [`conversion`] — RA/DEC parsing and coordinate utilities.
//! - [`earth_orientation`] — Precession, nutation, and obliquity models.
//! - [`kepler`] — Universal Kepler propagator and Lagrange f–g solver.
//! - [`orb_elem`] — State-vector to orbital-elements conversion.
//! - [`outfit_errors`] — Unified error enum for the whole crate.
//! - [`time`] — Time scale conversions and sidereal time.
//! - [`cache`] — Precomputed observer position cache.
//! - [`ephemeris`] — Apparent position computation and ephemeris residuals.
//! - [`observer_extension`] — Geocentric and heliocentric observer position routines.

// === Modules (internals). Keep public modules as they are; the facade is built via `pub use` below.

/// Constants and astronomical unit conversions.
pub mod constants;

/// Coordinate conversions (RA/DEC, spherical ↔ cartesian, etc.).
pub mod conversion;

/// Earth orientation parameters and related corrections (nutation, precession).
pub mod earth_orientation;

/// Initial Orbit Determination algorithms (Gauss method).
pub mod initial_orbit_determination;

/// Differential orbit correction utilities.
pub mod differential_orbit_correction;

/// JPL ephemerides management (Horizon/NAIF kernels).
pub mod jpl_ephem;

/// Ephemeris interpolation and access.
pub mod ephemeris;

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

/// Ground-observer geometry: body-fixed and heliocentric position routines.
pub mod observer_extension;

/// Orbit propagation strategies (TwoBody, N-body DOP853).
pub mod propagator;

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
    FullOrbitResult, AU, GAUSS_GRAV, IODRMS, RADEG, RADH, RADSEC, SECONDS_PER_DAY, T2000, VLIGHT_AU,
};

// JPL ephemeris enum for runtime inspection (optional but convenient)
pub use crate::jpl_ephem::JPLEphem;

// Ephemeris façade
pub use crate::ephemeris::{
    AberrationOrder, ApparentPosition, BodyGeometry, Combined, EphemerisConfig, EphemerisEntry,
    EphemerisMode, EphemerisOutputKind, EphemerisRequest, EphemerisResult, FullOrbitResultExt,
    Geometry, ObserverRequest, Position,
};

// IOD entry points and result types
pub use crate::initial_orbit_determination::obs_dataset_api::FitIOD;

// Differential correction entry points
pub use crate::differential_orbit_correction::obs_dataset_api::FitLSQ;
pub use crate::differential_orbit_correction::{
    DifferentialCorrectionConfig, DifferentialCorrectionOutput,
};

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
