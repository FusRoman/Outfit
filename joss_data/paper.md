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

The NSF-DOE Vera C. Rubin Observatory [@LSST_from_science_to_design] will conduct the Legacy Survey of Space and Time (LSST), which will produce an unprecedented optical alert stream exceeding ten million alerts per night. A substantial fraction of these alerts will correspond to moving Solar System objects. LSST is expected to increase the number of known asteroids by at least an order of magnitude—from roughly 1.3 million today to several million new discoveries over the survey lifetime [@rubinobservatorylsstsolarsystemsciencecollaboration2020scientificimpactverac]. Handling this volume requires fast and reliable tools to read astrometric observations and estimate preliminary orbits automatically.

Alert brokers such as Fink [@Fink_2020] receive the LSST alert stream ahead of all other users and classify alerts in real time. Their position at the earliest stage of data dissemination places them on the front line for recognising and flagging previously unknown asteroids. Early identification of such objects improves alert classification, enables rapid computation of ephemerides for follow-up, and supports linking future detections to newly derived orbits. It also allows brokers to filter out moving objects when searching for fast optical transients, an important capability for multi-messenger astronomy.

# Statement of need

Preliminary orbit determination is a foundational step in Solar System science. Researchers, survey teams, alert brokers, and even advanced amateur observers rely on software that can ingest astrometric measurements and rapidly estimate orbits for newly detected objects. Within the astronomical community, tools such as OrbFit [@OrbFit] or Find_Orb [@FindOrb] have long played an essential role by providing scientifically trusted implementations of classical orbit-determination algorithms. They remain widely used in observatories and research groups around the world.

However, the software landscape in which these tools operate is changing. Modern surveys, real-time alert brokers, and large-scale data-processing systems require libraries that can be integrated directly into automated workflows, support programmatic access, and operate reliably at scale. Developers building next-generation pipelines need open-source components that are composable, reproducible, and easy to embed within larger scientific software stacks—whether for nightly survey operations, real-time alert classification, or educational and citizen-science applications.

**Outfit[^outfit]** was created to meet this need. It provides a modern, open-source library for preliminary orbit determination that can be called directly from within automated pipelines, allowing users to process observations in memory and estimate orbits without relying on external executables. By offering a clear and well-tested programming interface, Outfit enables the community to integrate classical orbit-determination techniques into contemporary survey pipelines and alert-broker infrastructures, ensuring that these methods remain accessible and usable as astronomical data volumes continue to grow.

[^outfit]: <https://crates.io/crates/outfit>


# Software description

## Functionality

Outfit transforms raw astrometric observations into preliminary orbital solutions for small Solar System bodies and propagates these solutions to other epochs. It combines robust input handling, explicit modelling of observatories and reference frames, and classical initial-orbit-determination algorithms in a form suited to automated, high-throughput pipelines.

Designed primarily as an embeddable library, Outfit accepts observations directly as in-memory batches generated upstream, enabling programmatic orbit fitting without external executables or temporary files. For interoperability with existing workflows, it also ingests standard formats—including MPC 80-column [@MPC80col], ADES XML [@ADES], and columnar layouts such as Parquet/Arrow [@Arrow; @Parquet]—with validation and time-scale normalisation performed on load.

Initial orbit determination follows the classical Gauss method on triplets of optical observations. The solver performs iterative refinement, numerical and geometric acceptability checks, and RMS residual evaluation. The eighth-degree polynomial arising in Gauss’s formulation is solved using an Aberth–Ehrlich method (via the aberth crate [@Aberth]), ensuring stable and efficient root finding. Observatory geometry is handled explicitly via MPC codes or geodetic coordinates, which are converted into consistent geocentric and topocentric frames.

Core numerical operations rely on nalgebra [@nalgebra] for efficient linear-algebra routines. Outfit supports multiple orbital element representations (Keplerian, equinoctial, and cometary), two-body propagation, and planetary ephemerides from JPL DE kernels (e.g., DE440) [@DE440; @NAIF]. Time handling uses hifitime [@hifitime] for consistent UTC/TT/TDB conversions.

For survey-scale workloads, Outfit provides batched—and optionally parallel—processing of trajectories. Parallelism is implemented with Rayon [@rayon], allowing multi-core evaluation of independent orbit fits while preserving deterministic behaviour when required. This design lets users tune memory use and throughput effectively across a wide range of data volumes.

## Quality, documentation, and testing

Outfit includes unit, integration, and property-based tests covering parsing, observer handling, Gauss IOD, and trajectory batching, as well as regression tests on real astrometric datasets (e.g., 2015AB, 8467, 33803) to ensure numerical stability. Continuous integration enforces automated testing, coverage tracking, and semantic-versioning checks[@SemVer]. The repository includes additional guides and a versioned changelog [@KeepAChangelog].

Documentation is available on docs.rs and supported by runnable examples demonstrating end-to-end usage, such as reading Parquet observations and running Gauss IOD:

```bash
cargo run --release --example parquet_to_orbit --features "jpl-download"
```

These examples provide small, reproducible entry points for validating behaviour.

## Python bindings (pyOutfit)

Outfit is also exposed to Python via **pyOutfit[^pyoutfit]**, a binding built with PyO3 and distributed with maturin. It provides a thin, idiomatic Python interface to the Rust core, enabling integration in notebooks and Python pipelines. The binding is available on PyPI and ships precompiled Rust manylinux wheels for Python 3.10, 3.11, and 3.12.

[^pyoutfit]: <https://github.com/FusRoman/pyOutfit>

Project documentation is published on GitHub Pages using MkDocs Material, and the public API is type‑annotated via .pyi stub files to improve IDE completion and static type checking. Examples mirror the Rust examples to ease cross‑language adoption. The package includes a test suite, and the documentation examples are sourced from Python files and executed automatically in continuous integration to keep the documentation accurate and up to date.

# Limitations and future work

Outfit currently focuses on optical astrometry and provides initial orbit determination through the classical Gauss method. It does not yet support radar or radio observables, hyperbolic orbit fitting, or full least-squares refinement, and therefore cannot produce formal covariances or robust uncertainty estimates. In addition, while planetary ephemerides are fully integrated, the library does not yet expose a high-level API for generating predicted ephemerides from estimated orbits.

Future development will extend the IOD capabilities beyond Gauss, notably through the addition of a Väisälä solver and a full differential-corrections pipeline with covariance output. Support for hyperbolic trajectories will enable the analysis of interstellar candidates, and a unified API for orbit-based ephemeris prediction will facilitate downstream residual analysis. Longer-term plans include incorporating radar and Doppler measurements to complement optical data and broaden Outfit’s applicability across Solar System studies.

# Acknowledgements

We acknowledge the maintainers of OrbFit for foundational work in preliminary orbit determination, the JPL Horizons/NAIF teams for ephemerides and kernels, and contributors and users who provided feedback on API ergonomics and tests. Development benefited from the Rust ecosystem and libraries cited.

# References

