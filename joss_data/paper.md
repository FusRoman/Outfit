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
  - name: Hadrien Grasland
    orcid: 0000-0003-3379-7100
    affiliation: "1"
  - name: Julien Peloton
    orcid: 0000-0002-8560-4449
    affiliation: "1"
affiliations:
  - index: 1
    name: Université Paris-Saclay, CNRS/IN2P3, IJCLab, 91405 Orsay, France
date: 5 December 2025
bibliography: paper.bib
---

# Summary

The NSF-DOE Vera C. Rubin Observatory [@LSST_from_science_to_design] conducts the Legacy Survey of Space and Time (LSST), which produces an unprecedented optical alert stream exceeding ten million alerts per night. A substantial fraction of these alerts corresponds to moving Solar System objects. LSST is expected to increase the number of known asteroids by at least an order of magnitude—from roughly 1.3 million today to several million new discoveries over the survey lifetime [@rubinobservatorylsstsolarsystemsciencecollaboration2020scientificimpactverac]. Handling this volume requires fast and reliable tools to read astrometric observations and estimate preliminary orbits automatically.

Alert brokers such as Fink [@Fink_2020] receive the LSST alert stream ahead of all other users and classify alerts in real time. Their position at the earliest stage of data dissemination places them on the front line for recognising and flagging previously unknown asteroids. Early identification of such objects improves alert classification, enables rapid computation of ephemerides for follow-up, and supports linking future detections to newly derived orbits. It also allows brokers to filter out moving objects when searching for fast optical transients, an important capability for multi-messenger astronomy.

**Outfit[^outfit]** is a software library that performs these calculations automatically. It reads standardised observation records, applies classical orbit-determination algorithms, and returns preliminary orbital solutions that can be used for object identification, follow-up planning, and alert classification. Outfit is designed to be embedded directly into automated data-processing pipelines, enabling real-time analysis of large-scale survey data.

[^outfit]: <https://crates.io/crates/outfit>

# Statement of need

Preliminary orbit determination is a foundational step in Solar System science. Researchers, survey teams, alert brokers, and even advanced amateur observers rely on software that can ingest astrometric measurements and rapidly estimate orbits for newly detected objects.

Modern surveys operate at unprecedented scales, and their data-processing infrastructures require software components that are composable, reliable, and easy to integrate programmatically. Developers building next-generation pipelines need open-source components that are composable, reproducible, and easy to embed within larger scientific software stacks—whether for nightly survey operations, real-time alert classification, or educational and citizen-science applications.

Outfit addresses this need by providing a modern library interface for preliminary orbit determination. It allows users to process observations in memory and estimate orbits programmatically, making classical orbit-determination methods accessible and usable within contemporary survey pipelines and alert-broker infrastructures. This design is particularly important for real-time systems such as Fink [@Fink_2020], which receive the LSST alert stream ahead of all other users and must classify alerts within seconds. Early identification of moving objects improves alert classification, enables rapid computation of ephemerides for follow-up, supports linking future detections to newly derived orbits, and allows brokers to filter out asteroids when searching for fast optical transients—an important capability for multi-messenger astronomy.

# State of the field

Within the astronomical community, tools such as OrbFit [@OrbFit] and Find_Orb [@FindOrb] have long played an essential role by providing scientifically trusted implementations of classical orbit-determination algorithms. OrbFit, developed in Fortran, is a mature and rigorously validated tool that has been widely adopted in observatories and research institutions. Find_Orb, implemented in C, is similarly well established and provides comparable functionality. Both tools are distributed as standalone executables that read input files and write output files, a design pattern that has served the community well for decades. Recent validation efforts such as Orb_It [@orb_it] provide end-to-end testing frameworks for these tools, confirming their internal consistency and accuracy.

However, these tools were not designed for programmatic integration into automated pipelines. Invoking them from within a larger application requires launching external processes, writing temporary input files, parsing textual output, and managing inter-process communication—all sources of friction that become significant bottlenecks at the scales required by modern surveys. Furthermore, neither tool exposes a library API that would allow direct in-memory processing.

While OrbFit and Find_Orb remain the gold standard for scientific accuracy and validation, their executable-based design makes them ill-suited for the use cases Outfit targets. Rather than contribute additional features to these existing tools, we chose to build a new library from the ground up with a different architectural philosophy: an embeddable, memory-safe library with a clear API designed for automated workflows. This design choice represents a meaningful scholarly contribution, as it enables classical orbit-determination methods to be integrated into modern data-processing systems that would otherwise struggle to adopt executable-based tools.

# Software design

The architecture of Outfit reflects a deliberate set of trade-offs made to balance scientific accuracy, computational performance, and ease of integration. The most fundamental design decision was to implement Outfit as a library rather than an executable. This choice prioritises programmatic access over interactive use, enabling developers to embed orbit-determination capabilities directly into their applications without inter-process communication overhead.

Outfit is implemented in Rust [@rust_book], a programming language that provides memory safety without garbage collection and zero-cost abstractions for performance-critical code. This choice was motivated by the need for both speed and reliability in a library intended for automated, high-throughput pipelines. Rust's ownership system eliminates entire classes of runtime errors (null pointer dereferences, use-after-free, data races) without imposing runtime overhead, making it well-suited for scientific computing applications that must be both fast and correct.

The library accepts observations directly as in-memory data structures, avoiding file I/O in the critical path. For interoperability with existing workflows, Outfit also supports standard formats including MPC 80-column [@MPC80col], ADES XML [@ADES], and columnar layouts such as Parquet/Arrow [@Arrow; @Parquet]. These parsers perform validation and time-scale normalisation on load, ensuring that downstream orbit-fitting code operates on consistently formatted data.

Initial orbit determination follows the classical Gauss method on triplets of optical observations. The solver performs iterative refinement with numerical and geometric acceptability checks, evaluating RMS residuals to assess fit quality. The eighth-degree polynomial arising in Gauss's formulation is solved using an Aberth–Ehrlich root-finding method (via the aberth crate [@Aberth]), which provides numerical stability and efficiency. Observatory geometry is handled explicitly via MPC codes or geodetic coordinates, converted into consistent geocentric and topocentric reference frames. Outfit supports multiple orbital element representations (Keplerian, equinoctial, and cometary), two-body propagation, and planetary ephemerides from JPL DE kernels (e.g., DE440) [@DE440; @NAIF]. Time handling uses the hifitime library [@hifitime] for consistent UTC/TT/TDB conversions.

For survey-scale workloads, Outfit provides batched and optionally parallel processing of trajectories. Parallelism is implemented with Rayon [@rayon], allowing multi-core evaluation of independent orbit fits while preserving deterministic behaviour when required. This design lets users tune memory use and throughput effectively across a wide range of data volumes. Core numerical operations rely on nalgebra [@nalgebra] for efficient linear-algebra routines, leveraging a mature and well-tested foundation for vector and matrix operations.

To maximise accessibility, Outfit is also exposed to Python via **pyOutfit[^pyoutfit]**, a binding built with PyO3 [@pyo3] and distributed with maturin [@maturin]. This provides a thin, idiomatic Python interface to the Rust core, enabling integration in notebooks and Python-based pipelines. The binding is available on PyPI and ships precompiled manylinux wheels for Python 3.10, 3.11, and 3.12 at the time of writing. Project documentation is published on GitHub Pages using MkDocs Material, and the public API is type-annotated via .pyi stub files to improve IDE completion and static type checking.

[^pyoutfit]: <https://github.com/FusRoman/pyOutfit>

These design choices matter for research applications because they enable orbit determination to be embedded directly into real-time alert classification systems (such as Fink), distributed computing frameworks (such as Apache Spark), and other automated pipelines where launching external executables would introduce unacceptable latency or complexity. By providing a library interface with low integration friction, Outfit makes classical orbit-determination methods accessible to a broader range of research workflows.

# Research impact statement

Outfit has demonstrated realized impact through integration into operational research infrastructure and scholarly publication. The library is actively used within the Fink broker [@Fink_2020], which processes the LSST alert stream in real time and serves as a critical component of the survey's early data-dissemination pipeline. This operational deployment validates Outfit's design for high-throughput, low-latency applications.

The research need that motivated Outfit's development is documented in [@Le_Montagner_2023], which used OrbFit to fit tens of thousands of trajectory hypotheses by invoking the Fortran executable through Apache Spark across multiple compute nodes—a workflow that required dozens of CPU cores and took several minutes to complete. With Outfit's library interface, the same analysis runs in seconds on a single multi-core machine, demonstrating the performance gains enabled by in-memory processing and eliminating the overhead of inter-process communication and file I/O.

To ensure numerical reliability, Outfit includes regression tests on real astrometric datasets (2015AB, 8467, 33803) that validate orbit-fitting accuracy against known solutions returned by OrbFit. These tests run automatically in continuous integration and serve as a reproducible benchmark for the library's correctness. The library is distributed on crates.io (the Rust package registry) and PyPI (the Python Package Index), with precompiled binaries for multiple platforms, lowering the barrier to adoption. Documentation is available on docs.rs and GitHub Pages, and the repository includes runnable examples demonstrating end-to-end usage:

```bash
cargo run --release --example parquet_to_orbit --features "jpl-download"
```

# AI usage disclosure

GitHub Copilot's auto-completion feature was used in VS Code during implementation, and the GitHub Copilot coding agent assisted with writing testing infrastructure, documentation, and refactoring. AI-assisted language editing was used to improve wording and clarity in parts of the manuscript.  All code and text generated or modified by AI tools was carefully reviewed and proofread by the authors, and all AI-generated code was extensively tested.

# Limitations and future work

Outfit currently focuses on optical astrometry and provides initial orbit determination through the classical Gauss method. It does not yet support radar or radio observables, hyperbolic orbit fitting, or full least-squares refinement, and therefore cannot produce formal covariances or robust uncertainty estimates. In addition, while planetary ephemerides are fully integrated, the library does not yet expose a high-level API for generating predicted ephemerides from estimated orbits.

Future development will extend the IOD capabilities beyond Gauss, notably through the addition of a Väisälä solver [@Williams1997] and a full differential-corrections pipeline with covariance output. Support for hyperbolic trajectories will enable the analysis of interstellar candidates, and a unified API for orbit-based ephemeris prediction will facilitate downstream residual analysis.

# Acknowledgements

We acknowledge the maintainers of OrbFit for foundational work in preliminary orbit determination, the JPL Horizons/NAIF teams for ephemerides and kernels, and contributors. Development benefited from the Rust ecosystem and libraries cited.

# References
