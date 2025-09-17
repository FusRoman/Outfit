<p align="center">
  <img src="https://raw.githubusercontent.com/FusRoman/Outfit/main/outfit_logo.jpeg" alt="Outfit Logo" width="160"/>
</p>

<div align="center">
<h1>Outfit</h1>

A fast, safe, and extensible Rust library for **managing astrometric observations** and **determining preliminary orbits** of small bodies. Outfit reads common observation formats (MPC 80-column, ADES XML, Parquet), performs **initial orbit determination (IOD)** with the **Gauss method**, manages observers (topocentric geometry), and interfaces with **JPL ephemerides** for accurate state propagation.
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
    <img src="https://img.shields.io/badge/MSRV-1.82%2B-orange" alt="MSRV"/>
  </a>
</p>

> **Why Outfit?**  
> Modern asteroid pipelines need a library that is **fast (Rust)**, **reproducible**, and **easy to integrate** in data-intensive workflows (batch files, Parquet, CI/benchmarks). Outfit re-implements classic OrbFit IOD logic with a **memory-safe**, **modular** design and production-grade ergonomics (features, docs, tests, benches). It is built to:
> - ingest **large datasets** efficiently (columnar Parquet, batch APIs);
> - run **deterministic IOD** with controlled noise and repeatable seeds;
> - interface **cleanly** with JPL ephemerides (e.g., DE440);
> - provide a **clean API** that composes well across projects.

---

## Table of Contents

- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Data Formats](#data-formats)
- [Initial Orbit Determination](#initial-orbit-determination)
- [Observers & Reference Frames](#observers--reference-frames)
- [Cargo Feature Flags](#cargo-feature-flags)
- [Performance & Reproducibility](#performance--reproducibility)
- [Roadmap](#roadmap)
- [Contributing](#contributing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Features

- **Observation I/O**
  - MPC **80-column** files
  - **ADES XML** files
  - **Parquet** batches for high-throughput pipelines
- **Observer management**
  - Lookup by **MPC code**
  - Topocentric geometry (geocentric & heliocentric positions, AU, J2000)
- **Initial Orbit Determination**
  - **Gauss method** on observation triplets
  - Iterative velocity correction with Lagrange coefficients
  - Dynamic acceptability filters (perihelion, eccentricity, geometry)
  - RMS evaluation on extended arcs
- **Ephemerides**
  - Interface with **JPL** ephemerides (e.g., **DE440**)
- **Error models**
  - Built-in astrometric uncertainty models (e.g., FCCT14)
- **Batch processing & benchmarking**
  - Stream triplets, evaluate, and rank candidates
  - Criterion-based micro/macro benchmarks
  - Optional parallel batch IOD using Rayon (feature: `parallel`)

---

## Installation

Add Outfit to your `Cargo.toml`:

~~~toml
[dependencies]
outfit = "2.0.0"
~~~

Enable automatic ephemeris download (JPL DE440) with the `jpl-download` feature:

~~~toml
[dependencies]
outfit = { version = "2.0.0", features = ["jpl-download"] }
~~~

Enable a CLI-style progress bar for long loops:

~~~toml
[dependencies]
outfit = { version = "2.0.0", features = ["progress"] }
~~~

Combine features as needed (example):

~~~toml
[dependencies]
outfit = { version = "2.0.0", features = ["jpl-download", "progress", "parallel"] }
~~~

---

## Quick Start

The crate ships with several **ready-to-run examples** in the [`examples/`](examples) directory.  
They demonstrate end-to-end workflows such as:

- Reading observations (MPC 80-column, ADES XML, Parquet)
- Building a `TrajectorySet` (now via `TrajectoryFile` ingestion helpers)
- Running Gauss initial orbit determination (single triplet or full batch)
- Inspecting and printing orbital elements with RMS statistics

Run an example directly with Cargo:

```bash
cargo run --release --example parquet_to_orbit --features jpl-download
```

This will:

- Load observations from a Parquet file,

- Identify the observer by MPC code,

- Run the Gauss IOD pipeline,

- Print the best-fit orbit with its RMS value.

---

## Data Formats

- **MPC 80-column** – [Minor Planet Center fixed-width astrometry format](https://minorplanetcenter.net/iau/info/OpticalObs.html)  
- **ADES XML** – [IAU Astrometric Data Exchange Standard (ADES)](https://minorplanetcenter.net/iau/info/ADES.html)  
- **Parquet** – [Apache Parquet](https://parquet.apache.org/) columnar format for large batch processing  
  Typical columns: `ra`, `dec`, `jd` or `mjd`, `trajectory_id` (configurable)


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

### Parallel batches

Compile with `--features parallel` and prefer the adaptive batching API for very large sets:

```bash
cargo run --release --example parquet_to_orbit --features "jpl-download,parallel"
```

Then call `estimate_all_orbits_in_batches_parallel` (same signature + extra `batch_size` argument) for improved throughput when trajectory sets exceed CPU cache friendliness.

### Result helpers

Access individual solutions ergonomically:

```rust
use outfit::trajectories::trajectory_fit::{gauss_result_for, take_gauss_result};
// Borrow without moving:
if let Some(Ok((g, rms))) = gauss_result_for(&results, &some_object) {
   println!("Semi-major axis: {} AU (rms={rms:.3})", g.keplerian.a);
}
```

Errors are explicit (no `Option` inside `Result`); pathological numeric scores (`NaN`, `inf`) surface as `OutfitError::NonFiniteScore`.

---

## Observers & Reference Frames

- Observer lookup by **MPC code**.
- Geocentric and heliocentric positions in **AU**, **J2000** (equatorial).
- Earth orientation (nutation, precession) and **aberration** corrections via internal reference-system utilities.

---

## Cargo Feature Flags

The crate keeps the feature surface intentionally small and orthogonal. All core parsing (MPC 80-col, ADES XML, Parquet) and Gauss IOD logic are always compiled; features only toggle optional runtime dependencies.

| Feature        | Adds                                                    | Notes |
|----------------|---------------------------------------------------------|-------|
| `jpl-download` | `reqwest`, `tokio` (multi-thread RT), on-demand fetch   | Downloads and caches JPL ephemerides (e.g. DE440) in the user data dir on first access. Offline re-use afterward. |
| `progress`     | `indicatif` progress bars + moving average timing       | Enabled in long batch IOD loops; zero cost when not used. |
| `parallel`     | `rayon`                                                 | Enables `TrajectoryFit::estimate_all_orbits_in_batches_parallel`; set `IODParams.batch_size` to tune batch granularity. |

Planned (not yet implemented): least-squares refinement and alternative IOD methods will likely get their own feature gates.

---

## Performance & Reproducibility

- **Deterministic runs** by default (set RNG seeds explicitly when noise is used).
- **Batch-friendly APIs** (Parquet; streaming triplets).
- Avoid ephemeris I/O in hot paths by **precomputing observer positions**.
- Benchmarks via **criterion** (see below).

**Tips**
- Compile with `--release` for production.
- Keep ephemerides cached locally (with `jpl-download` enabled) to avoid I/O stalls.
- Use the `progress` feature to instrument long loops without cluttering core logic.
- `ObservationBatch` conversions (deg/arcsec → rad) happen once; avoid re-normalizing angles inside tight loops.
- Group observations by unique epoch to reuse cached observer positions (already done by built-in ingestion pipelines).
- Handle errors rather than filtering silently: a non-finite RMS is raised early as `OutfitError::NonFiniteScore` to prevent propagating invalid states.
- For very large trajectory sets, tune `IODParams.batch_size` (with `parallel`) so each batch fits comfortably in cache (start with 4–8 and benchmark).

---

## Roadmap

- **Full least-squares orbit fitting** across full arcs
- **Hyperbolic & parabolic orbits** (e ≥ 1) for interstellar candidates
- **Alternative IOD methods** (e.g., **Vaisala**)
- **Ephemerides backends** (e.g., ANISE/SPICE integration)

See the issue tracker for the latest details and discussion.

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
cargo bench --features jpl-download
~~~

---

## License

This project is licensed under the **CeCILL-C** Free Software License Agreement. See the `LICENSE` file.

---

## Acknowledgements

- The **OrbFit** project and its authors for foundational algorithms and references.
- The maintainers of **JPL** ephemerides and the broader open-source ecosystem (linear algebra, I/O, benchmarking crates) that Outfit builds upon.
