---
title: "Outfit: a fast, safe Rust library for astrometric observations and preliminary orbit determination"
tags:
  - Rust
  - astronomy
  - orbit determination
  - astrometry
  - ephemerides
authors:
  - name: Roman Le Montagner
    orcid: 0000-0002-6099-8939
    affiliation: "1"
affiliations:
  - index: 1
    name: IJCLab (CNRS/IN2P3), Université Paris-Saclay, France
date: 5 December 2025
bibliography: paper.bib
---

# Summary

Upcoming wide-field surveys—most notably the Vera C. Rubin Observatory’s Legacy Survey of Space and Time (LSST)—are expected to increase the number of known small Solar System bodies by an order of magnitude. This scale demands efficient, robust algorithms that can derive preliminary orbits quickly and reproducibly from large, heterogeneous datasets. Outfit is a Rust library for reading astrometric observations and computing preliminary orbits of small bodies using an implementation of the Gauss method. It provides high-throughput ingestion of common formats (Minor Planet Center 80-column, ADES XML, and Parquet), explicit observer management (MPC codes, topocentric geometry), and accurate state propagation via JPL planetary ephemerides (e.g., DE440). The library emphasizes speed, memory safety, reproducibility, and modular design: observation parsing and trajectory batching are decoupled from initial orbit determination (IOD) logic, and all public APIs are documented and tested. Outfit re-implements classic OrbFit IOD behavior in a modern, production-grade codebase, suitable for survey-scale pipelines (e.g., LSST/ZTF-like workloads).

# Statement of need

Astrometric pipelines for near-Earth and small-body surveys must ingest heterogeneous observation formats, manage observers, and determine preliminary orbits efficiently and reproducibly. Existing tools (e.g., OrbFit [@OrbFit]) are widely used and trusted, but often expose monolithic interfaces and legacy implementations that are harder to embed in data-intensive workflows. Outfit fills this gap by offering:

- A memory-safe, high-performance core implemented in Rust.
- Unified I/O across MPC 80-column [@MPC80col] and ADES XML [@ADES] with a columnar option (Parquet/Arrow [@Arrow; @Parquet]) for large batches.
- Explicit observer handling with MPC codes and topocentric geometry.
- Deterministic Gauss IOD with controlled numerical behavior and explicit error types.
- Clean integration with JPL ephemerides (DE series) [@DE440; @NAIF] for accurate state propagation
- Accurate time-scale handling via `hifitime` [@hifitime].
- Optional parallel computation for Gauss IOD (`--features parallel`) to scale across large batches.

This combination enables reproducible, scalable preliminary orbit determination and downstream evaluation (RMS of normalized residuals) within modern data processing stacks. Outfit is intended for researchers building survey pipelines, teaching orbit determination, or benchmarking algorithmic variants on large datasets.

# Software description

## Functionality

Outfit provides:

- Observation parsers: MPC 80-column, ADES XML, and Parquet.
- Trajectory batching utilities and adaptive batch APIs for large sets.
- Gauss IOD on triplets with acceptability filters and Lagrange coefficient refinement.
- Observer management (MPC codes, geodetic/geocentric geometry) and pretty-printing utilities.
- Ephemeris integration via JPL (DE440 and similar) with optional automatic download.
- Display helpers for compact/wide/ISO tables of observations.

## Architecture

The crate is organized into modules reflecting responsibilities: `observations` (parsers, types), `trajectories` (batch readers, fit helpers), `initial_orbit_determination` (Gauss method), `observers` (MPC observatories, geometry), `jpl_ephem` (horizon/NAIF interfaces), and `orbit_type` (keplerian/equinoctial/cometary elements). Public APIs prefer explicit types and error handling (`OutfitError`). Performance-critical routines use `nalgebra` for vector math [@nalgebra] and optional parallelism via `rayon` [@rayon]. Columnar ingestion relies on Arrow/Parquet [@Arrow; @Parquet]. Time handling uses `hifitime` [@hifitime] for TT/UTC conversions and Earth orientation utilities.

## Quality, documentation, and testing

Outfit includes unit tests and integration tests covering parsing (MPC/ADES/Parquet), observer handling, Gauss IOD numerical behavior, and trajectory batching (see `tests/`). In addition, we use property‑based testing (via `proptest`) to stress numerical invariants and edge cases over wide parameter ranges, improving robustness beyond hand‑crafted cases. For accuracy, we run regression tests on real observation datasets and compare against OrbFit outputs [@OrbFit]: notably the objects 2015AB, 8467, and 33803 are validated with MPC 80‑column inputs (see `tests/data/` and `orbfit_tests/`).

Software quality is enforced through continuous integration (CI): we run automatic tests on each change, track code coverage, and perform semantic‑versioning checks (e.g., `semver` compatibility) to ensure API stability across releases. Version history is maintained via git with a clear changelog and tagged releases; this structured release process is not available in OrbFit’s upstream distribution. Examples in `examples/` demonstrate end‑to‑end workflows (e.g., reading Parquet and running Gauss IOD). The repository follows semantic versioning and maintains a comprehensive changelog [@SemVer; @KeepAChangelog]. The crate documentation is published on docs.rs and the README provides installation, quick‑start, feature flags (including `parallel`), and performance notes.

## Performance and reproducibility

 Rust enables predictable performance and memory safety. For survey-scale workloads, Outfit offers Parquet ingestion and optional parallel Gauss IOD (`--features parallel`) with a configurable batch size set via `IODParams.batch_size`. Numeric stability is enforced via explicit acceptability checks and non-finite score detection. Determinism is promoted through clear random seeding paths and ordered iteration for reporting.

# Demonstration

The crate ships runnable examples demonstrating observation ingestion and Gauss IOD. For instance, running the Parquet example with automatic ephemeris download:

```bash
cargo run --release --example parquet_to_orbit --features "jpl-download"
```

Illustrative outputs include orbital elements and RMS of normalized residuals, and optional progress reporting (`--features progress`). Examples are documented in the README and serve as minimal, reproducible entry points for reviewers to validate functionality locally as per JOSS guidelines.

# Acknowledgements

We acknowledge the maintainers of OrbFit for foundational work in preliminary orbit determination, the JPL Horizons/NAIF teams for ephemerides and kernels, and contributors and users who provided feedback on API ergonomics and tests. Development benefited from the Rust ecosystem and libraries cited below.

# References

