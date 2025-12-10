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
    name: Université Paris-Saclay, CNRS/IN2P3, IJCLab, 91405 Orsay, France
date: 5 December 2025
bibliography: paper.bib
---

# Summary

The Vera C. Rubin Observatory will conduct the Legacy Survey of Space and Time (LSST), which will produce an unprecedented optical alert stream exceeding ten million alerts per night. A substantial fraction of these alerts will correspond to moving Solar System objects. LSST is expected to increase the number of known asteroids by at least an order of magnitude—from roughly 1.3 million today to several million new discoveries over the survey lifetime. Handling this volume requires fast and reliable tools to read astrometric observations and estimate preliminary orbits automatically.

Alert brokers such as Fink receive the LSST alert stream ahead of all other users and classify alerts in real time. Their position at the earliest stage of data dissemination places them on the front line for recognising and flagging previously unknown asteroids. Early identification of such objects improves alert classification, enables rapid computation of ephemerides for follow-up, and supports linking future detections to newly derived orbits. It also allows brokers to filter out moving objects when searching for fast optical transients, an important capability for multi-messenger astronomy.

# Statement of need

Preliminary orbit determination is a foundational step in Solar System science. Researchers, survey teams, alert brokers, and even advanced amateur observers rely on software that can ingest astrometric measurements and rapidly estimate orbits for newly detected objects. Within the astronomical community, tools such as OrbFit have long played an essential role by providing scientifically trusted implementations of classical orbit-determination algorithms. They remain widely used in observatories and research groups around the world.

However, the software landscape in which these tools operate is changing. Modern surveys, real-time alert brokers, and large-scale data-processing systems require libraries that can be integrated directly into automated workflows, support programmatic access, and operate reliably at scale. Developers building next-generation pipelines need open-source components that are composable, reproducible, and easy to embed within larger scientific software stacks—whether for nightly survey operations, real-time alert classification, or educational and citizen-science applications.

**Outfit** was created to meet this need. It provides a modern, open-source library for preliminary orbit determination that can be called directly from within automated pipelines, allowing users to process observations in memory and estimate orbits without relying on external executables. By offering a clear and well-tested programming interface, Outfit enables the community to integrate classical orbit-determination techniques into contemporary survey pipelines and alert-broker infrastructures, ensuring that these methods remain accessible and usable as astronomical data volumes continue to grow.

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

The crate is organised into focused modules with clear responsibilities:

![Outfit crate architecture overview.\label{fig:architecture}](outfit_architecture.png){ width=80% }

- observations — parsing MPC/ADES/Parquet, data validation, time-scale normalisation. It provides robust parsers for standard formats, normalises times to consistent scales, and validates fields and units to ensure inputs are well-formed before IOD.
- observers — MPC observatory codes, geodetic conversion, topocentric geometry. This module resolves MPC codes to observer locations, handles Earth models and coordinate conversions, and computes observer-centric vectors used in line-of-sight calculations.
- initial_orbit_determination — Gauss method, Lagrange coefficient iteration, acceptability filters, error reporting. It implements the Gauss triplet solver, refines f–g (Lagrange) coefficients iteratively, and applies numerical and geometric filters; all failures surface explicit, typed errors.
- trajectories — batching, grouping, and utilities for large candidate sets. It assembles observation triplets into candidate trajectories, supports adaptive batching for large datasets, and exposes helpers to evaluate and sort solutions.
- jpl_ephem — DE ephemeris access (local kernels or automatic download). Provides planetary states from JPL DE kernels with caching and optional automated download, ensuring accurate Sun/planet positions for IOD and propagation.
- orbit_type — Keplerian, equinoctial, and cometary elements with conversions. Provides element representations and conversions between standard orbital parameterisations, along with formatting utilities for reporting and downstream analysis.

The public API favours explicit types, immutable inputs where possible, and strongly typed time and angle units. Performance-critical routines rely on nalgebra, while parallelism is opt-in via a crate feature. Columnar ingestion uses Arrow/Parquet to support HPC and distributed workflows.

## Quality, documentation, and testing

The Outfit package includes unit tests and integration tests that cover parsing (MPC/ADES/Parquet), observer handling, Gauss IOD numerical behaviour, and trajectory batching (see `tests/`). In addition, we use property-based testing (via `proptest`) to stress numerical invariants and edge cases across wide parameter ranges, thereby improving robustness beyond hand-crafted cases. For accuracy, we run regression tests on real observation datasets and compare the results against OrbFit outputs [@OrbFit]. Notably, the objects 2015AB, 8467, and 33803 are validated using MPC 80-column inputs (see `tests/data/` and `orbfit_tests/`).

Software quality is enforced through continuous integration (CI): we run automated tests on each change, track code coverage, and perform semantic versioning checks (e.g., `semver` compatibility) to ensure API stability across releases. Version history is maintained via git with a clear changelog and tagged releases; this structured release process is not available in OrbFit's upstream distribution. Examples in `examples/` demonstrate end‑to‑end workflows (e.g., reading Parquet and running Gauss IOD). The repository follows semantic versioning and maintains a comprehensive changelog [@SemVer; @KeepAChangelog]. The crate documentation is published on docs.rs, and the README provides installation instructions, a quick start guide, feature flags (including `parallel`), and performance notes. Outfit is available on crates.io for easy installation and inclusion in other Rust projects via Cargo's dependency system.

## Performance and reproducibility

 Rust enables predictable performance and memory safety. For survey-scale workloads, Outfit offers Parquet ingestion and optional parallel Gauss IOD (`--features parallel`) with a configurable batch size set via `IODParams.batch_size`. Numeric stability is enforced via explicit acceptability checks and non-finite score detection. Determinism is promoted through clear random seeding paths and ordered iteration for reporting.

## Python bindings (pyOutfit)

Outfit is also exposed to Python via pyOutfit, a binding built with PyO3 and distributed with maturin: https://github.com/FusRoman/pyOutfit. It provides a thin, idiomatic Python interface to the Rust core (parsing MPC/ADES/Parquet, observer handling, and Gauss IOD), enabling integration in notebooks and Python pipelines. The binding is available on PyPI and ships precompiled Rust manylinux wheels for Python 3.10, 3.11, and 3.12.

Project documentation is published on GitHub Pages using MkDocs Material, and the public API is type‑annotated via .pyi stub files to improve IDE completion and static type checking. Examples mirror the Rust examples to ease cross‑language adoption. The package includes a test suite, and the documentation examples are sourced from Python files and executed automatically in continuous integration to keep the documentation accurate and up to date.

# Demonstration

The crate ships runnable examples demonstrating observation ingestion and Gauss IOD. For instance, running the Parquet example with automatic ephemeris download:

```bash
cargo run --release --example parquet_to_orbit --features "jpl-download"
```

Illustrative outputs include orbital elements and RMS of normalised residuals, and optional progress reporting (`--features progress`). Examples are documented in the README and serve as minimal, reproducible entry points to validate functionality locally.

# Limitations and future work

Current limitations:

- Initial orbit determination is limited to the Gauss method on optical astrometry; no radar/radio observables are supported yet.
- No full least-squares refinement (differential corrections), hence no formal covariance or robust uncertainty quantification.
- Hyperbolic trajectories are not yet supported in fitting, limiting interstellar object handling.
- Efficient ephemerides generation from estimated orbits (initial or refined) is not yet exposed as a high-level API.

Planned work:

- Add full orbital estimation via a least-squares method (as in OrbFit), including covariance output and convergence diagnostics.
- Add an alternative initial orbit method (Vaisala/Väisälä) to complement Gauss.
- Add support for hyperbolic orbits to fit interstellar candidates (e.g., 1I/'Oumuamua, 2I/Borisov, 3I/ATLAS).
- Add efficient ephemerides computation from estimated orbits (initial or full) for prediction and residual analysis.
- Add support for radio/radar observations (e.g., delay and Doppler), alongside optical astrometry.

# Acknowledgements

We acknowledge the maintainers of OrbFit for foundational work in preliminary orbit determination, the JPL Horizons/NAIF teams for ephemerides and kernels, and contributors and users who provided feedback on API ergonomics and tests. Development benefited from the Rust ecosystem and libraries cited below.

# References

