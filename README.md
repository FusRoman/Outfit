# Outfit

A fast, safe, and extensible Rust library for **managing astrometric observations** and **determining preliminary orbits** of small bodies. Outfit reads common observation formats (MPC 80-column, ADES XML, Parquet), performs **initial orbit determination (IOD)** with the **Gauss method**, manages observers (topocentric geometry), and interfaces with **JPL ephemerides** for accurate state propagation.

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

---

## Installation

Add Outfit to your `Cargo.toml`:

~~~toml
[dependencies]
outfit = "1.0.0"
~~~

Enable automatic ephemeris download (JPL DE440) with the `jpl-download` feature:

~~~toml
[dependencies]
outfit = { version = "1.0.0", features = ["jpl-download"] }
~~~

Enable a CLI-style progress bar for long loops:

~~~toml
[dependencies]
outfit = { version = "1.0.0", features = ["progress"] }
~~~

Combine features as needed (example):

~~~toml
[dependencies]
outfit = { version = "1.0.0", features = ["jpl-download", "progress", "parquet", "ades", "mpc80"] }
~~~

---

## Quick Start

The crate ships with several **ready-to-run examples** in the [`examples/`](examples) directory.  
They demonstrate end-to-end workflows such as:

- Reading observations (MPC 80-column, Parquet)
- Building a `TrajectorySet`
- Running Gauss initial orbit determination
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

## Observers & Reference Frames

- Observer lookup by **MPC code**.
- Geocentric and heliocentric positions in **AU**, **J2000** (equatorial).
- Earth orientation (nutation, precession) and **aberration** corrections via internal reference-system utilities.

---

## Cargo Feature Flags

| Feature         | Description                                                                 |
|-----------------|-----------------------------------------------------------------------------|
| `jpl-download`  | Automatically fetch JPL ephemerides on first use (e.g., DE440).            |
| `progress`      | Lightweight progress indicators (e.g., via `indicatif`) for long pipelines.|

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
