<p align="center">
  <img src="https://raw.githubusercontent.com/FusRoman/Outfit/main/outfit_logo.jpeg" alt="Outfit Logo" width="160"/>
</p>

<div align="center">
<h1>Outfit</h1>

A fast, safe, and modular Rust library for **initial orbit determination (IOD)** of Solar System small bodies. Outfit focuses exclusively on **orbital dynamics and computation**: Gauss IOD, Keplerian propagation, reference frame transformations, and JPL ephemeris support. Observation I/O and data structures are provided by the companion crate [**photom**](https://crates.io/crates/photom).
</div>

<p align="center">
  <a href="https://crates.io/crates/outfit">
    <img src="https://img.shields.io/crates/v/outfit.svg" alt="Crates.io"/>
  </a>
  <a href="https://docs.rs/outfit">
    <img src="https://docs.rs/outfit/badge.svg" alt="Docs.rs"/>
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-CeCILL--C-blue.svg" alt="License: CeCILL-C"/>
  </a>
  <a href="https://github.com/FusRoman/Outfit/actions">
    <img src="https://github.com/FusRoman/Outfit/actions/workflows/ci.yml/badge.svg" alt="CI"/>
  </a>
  <a href="https://codecov.io/gh/FusRoman/Outfit">
    <img src="https://codecov.io/gh/FusRoman/Outfit/branch/main/graph/badge.svg" alt="codecov"/>
  </a>
  <a href="https://www.rust-lang.org/">
    <img src="https://img.shields.io/badge/MSRV-1.94%2B-orange" alt="MSRV"/>
  </a>
</p>

> **Ecosystem split**  
> Outfit v3 is a focused **orbital dynamics engine**. Everything related to observation loading (MPC 80-column, ADES XML, Parquet), observer management, and astrometric data structures now lives in the separate **[photom](https://crates.io/crates/photom)** crate. Outfit depends on `photom` and exposes the `FitIOD` trait that bridges photom's `ObsDataset` into the IOD pipeline.

---

## Table of Contents

- [Features](#features)
- [Ecosystem: Outfit + photom](#ecosystem-outfit--photom)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Initial Orbit Determination](#initial-orbit-determination)
- [Cargo Feature Flags](#cargo-feature-flags)
- [Performance & Reproducibility](#performance--reproducibility)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Features

- **Initial Orbit Determination (IOD)**
  - **Gauss method** on observation triplets
  - Iterative velocity correction with Lagrange coefficients
  - Dynamic acceptability filters (perihelion, eccentricity, geometry)
  - RMS evaluation on extended arcs
- **Orbital elements**
  - Classical **Keplerian elements**
  - **Equinoctial elements** with conversions and two-body solver
  - **Cometary elements**
- **Ephemerides**
  - Interface with **JPL DE440** (Horizons and NAIF/SPICE formats)
- **Reference frames & preprocessing**
  - Precession, nutation (IAU 1980), aberration, light-time correction
  - Ecliptic ↔ equatorial conversions, RA/DEC parsing, time systems
- **Observer geometry**
  - Geocentric and heliocentric observer positions in AU, J2000
- **Batch IOD**
  - Single-trajectory and full-dataset IOD via the `FitIOD` trait
  - Optional parallel batch execution with Rayon (`parallel` feature)

---

## Ecosystem: Outfit + photom

Outfit is one half of a two-crate pipeline:

| Crate | Responsibility |
|-------|---------------|
| [**photom**](https://crates.io/crates/photom) | Observation I/O (MPC 80-col, ADES XML, Parquet), data structures (`ObsDataset`, `Observer`, error models), trajectory grouping |
| **outfit** | Pure orbital computation: Gauss IOD, Keplerian/equinoctial elements, JPL ephemerides, reference frames, residuals |

A typical workflow:

1. Use **photom** to load observations into an `ObsDataset`.
2. Call **outfit**'s `FitIOD::fit_iod` (single trajectory) or `FitIOD::fit_full_iod` (all trajectories) on the dataset.
3. Inspect the returned `GaussResult` and RMS values.

---

## Installation

Add Outfit and photom to your `Cargo.toml`:

~~~toml
[dependencies]
outfit = "3.0.0"
photom = { version = "0.1.0", features = ["mpc_80_col"] }
~~~

Enable optional features as needed:

~~~toml
[dependencies]
outfit  = "3.0.0"
photom  = { version = "0.1.0", features = ["mpc_80_col", "ades", "polars"] }
~~~

---

## Quick Start

### Single-trajectory IOD from an MPC 80-column file

```rust,no_run
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use hifitime::ut1::Ut1Provider;
use rand::{rngs::StdRng, SeedableRng};
use outfit::obs_dataset::FitIOD;
use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem};
use outfit::IODParams;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Load observations (photom handles I/O).
    let dataset = ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")?;

    // Load JPL DE440 ephemeris.
    let jpl_source: EphemFileSource = "horizon:DE440".try_into()?;
    let jpl = JPLEphem::new(&jpl_source)?;

    // Earth orientation corrections.
    let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long")?;

    // Configure IOD parameters.
    let params = IODParams::builder()
        .n_noise_realizations(10)
        .max_obs_for_triplets(50)
        .max_triplets(30)
        .build()?;

    let mut rng = StdRng::seed_from_u64(42);

    // Run Gauss IOD — outfit does the orbital computation.
    let (best_orbit, best_rms) = dataset.fit_iod(
        "K09R05F",
        &jpl,
        &ut1,
        &params,
        ObsErrorModel::FCCT14,
        &mut rng,
    )?;

    println!("Best orbit: {best_orbit}");
    println!("RMS: {best_rms:.6}");
    Ok(())
}
```

### Batch IOD — all trajectories at once

```rust,no_run
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use hifitime::ut1::Ut1Provider;
use rand::{rngs::StdRng, SeedableRng};
use outfit::obs_dataset::{FitIOD, FullOrbitResult};
use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem};
use outfit::IODParams;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dataset = ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")?;

    let jpl_source: EphemFileSource = "horizon:DE440".try_into()?;
    let jpl = JPLEphem::new(&jpl_source)?;
    let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long")?;

    let params = IODParams::builder()
        .n_noise_realizations(10)
        .max_obs_for_triplets(50)
        .build()?;

    let mut rng = StdRng::seed_from_u64(42);

    // Run Gauss IOD for every trajectory in the dataset.
    let results: FullOrbitResult = dataset.fit_full_iod(
        &jpl,
        &ut1,
        &params,
        ObsErrorModel::FCCT14,
        &mut rng,
    )?;

    for (traj_id, res) in &results {
        match res {
            Ok((gauss, rms)) => println!("{traj_id} → RMS = {rms:.4}\n{gauss}"),
            Err(e)           => eprintln!("{traj_id} → error: {e}"),
        }
    }
    Ok(())
}
```

---

## Initial Orbit Determination

Outfit implements the **Gauss method** on observation triplets:

1. Select candidate triplets `(i, j, k)` with spacing/quality heuristics.
2. Solve for **topocentric distances** and the heliocentric state at the central epoch.
3. Apply **dynamic acceptability** filters (perihelion, eccentricity, geometry).
4. **Refine** velocity with Lagrange coefficients (iterative correction loop).
5. Evaluate **RMS of normalized astrometric residuals** on an extended arc to rank solutions.

Designed for **robustness and speed** in survey-scale use (LSST, ZTF, ...).

---

## Performance & Reproducibility

- **Deterministic runs** by default (set RNG seeds explicitly when noise is used).
- **Precompute observer positions** to avoid ephemeris I/O in hot paths.
- Compile with `--release` for production.
- Keep ephemerides cached locally to avoid I/O stalls on repeated runs.
- For large trajectory sets, tune `IODParams.batch_size` (with `parallel`) so each batch fits comfortably in cache (start with 4–8 and benchmark).
- Handle errors explicitly: a non-finite RMS surfaces as `OutfitError::NonFiniteScore`.

---

## Roadmap

- **Full least-squares orbit fitting** across full arcs
- **Hyperbolic & parabolic orbits** (e ≥ 1) for interstellar candidates
- **Vaisälä method** for short arcs
- **Additional ephemerides backends** (e.g., ANISE/SPICE integration)

See the issue tracker for details and discussion.

---

## Contributing

Contributions are welcome!  
Please open an issue to discuss substantial changes. For PRs, please:

- Add tests for new functionality
- Keep public APIs documented
- Ensure `cargo fmt` and `cargo clippy` pass (`-D warnings` recommended)
- Include benchmarks when touching hot paths

A typical dev loop:

~~~bash
cargo fmt --all
cargo clippy --all-targets --all-features
cargo test --all-features
~~~

---

## License

This project is licensed under the **CeCILL-C** Free Software License Agreement. See the `LICENSE` file.

---

## Acknowledgements

- The **OrbFit** project and its authors for foundational algorithms and references.
- The maintainers of **JPL** ephemerides and the broader open-source ecosystem (linear algebra, I/O, benchmarking crates) that Outfit builds upon.
