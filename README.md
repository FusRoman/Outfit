<p align="center">
  <img src="https://raw.githubusercontent.com/FusRoman/Outfit/main/outfit_logo.jpeg" alt="Outfit Logo" width="160"/>
</p>

<div align="center">
<h1>Outfit</h1>

A fast, safe, and modular Rust library for **orbit determination** of Solar System small bodies. Outfit focuses exclusively on **orbital dynamics and computation**: Gauss IOD, differential orbit correction, ephemeris generation, Keplerian propagation, reference frame transformations, and JPL ephemeris support. Observation I/O and data structures are provided by the companion crate [**photom**](https://crates.io/crates/photom).
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
> Outfit v3 is a focused **orbital dynamics engine**. Everything related to observation loading (MPC 80-column, ADES XML, Parquet), observer management, and astrometric data structures now lives in the separate **[photom](https://crates.io/crates/photom)** crate. Outfit depends on `photom` and exposes the `FitIOD` trait that bridges photom's `ObsDataset` into the orbit fitting pipeline.

---

## Table of Contents

- [Features](#features)
- [Ecosystem: Outfit + photom](#ecosystem-outfit--photom)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Initial Orbit Determination](#initial-orbit-determination)
- [Differential Orbit Correction](#differential-orbit-correction)
- [Uncertainty Propagation](#uncertainty-propagation)
- [Ephemeris Generation](#ephemeris-generation)
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
- **Differential Orbit Correction**
  - Iterative **Newton–Raphson least-squares** refinement of equinoctial elements
  - Projection-based **outlier rejection** loop (chi-squared per observation)
  - Covariance matrix estimation with posterior uncertainty rescaling
  - Configurable free/fixed element mask (useful for short arcs)
  - Stagnation and divergence detection with typed error variants
  - `FitLSQ` trait: per-trajectory pipeline (IOD seed → differential correction) on any `ObsDataset`
- **Ephemeris generation**
  - Predict apparent sky positions `(RA, Dec)` and distances from any `OrbitalElements`
  - Compute geometric quantities: **phase angle**, **solar elongation**, **radial velocity**, apparent angular rates
  - Combined mode computes position and geometry in a single propagation
  - Three generation modes per observer: `Single`, `Range` (uniform grid), `At` (arbitrary epoch list)
  - Multiple observers in one typed `EphemerisRequest`; per-epoch errors collected without aborting the batch
  - Choice of **two-body (Keplerian)** or **N-body (DOP853)** propagator; first- or second-order aberration correction
- **Orbital elements**
  - Classical **Keplerian elements**
  - **Equinoctial elements** with conversions and two-body solver
  - **Cometary elements**
- **JPL ephemerides**
  - Interface with **JPL DE440** (Horizons and NAIF/SPICE formats)
- **Reference frames & preprocessing**
  - Precession, nutation (IAU 1980), aberration, light-time correction
  - Ecliptic ↔ equatorial conversions, RA/DEC parsing, time systems
- **Observer geometry**
  - Geocentric and heliocentric observer positions in AU, J2000
- **Batch processing**
  - Single-trajectory and full-dataset IOD via the `FitIOD` trait
  - Full least-squares fitting for all trajectories via the `FitLSQ` trait
  - Optional parallel batch execution with Rayon (`parallel` feature)

---

## Ecosystem: Outfit + photom

Outfit is one half of a two-crate pipeline:

| Crate | Responsibility |
|-------|---------------|
| [**photom**](https://crates.io/crates/photom) | Observation I/O (MPC 80-col, ADES XML, Parquet), data structures (`ObsDataset`, `Observer`, error models), trajectory grouping |
| **outfit** | Pure orbital computation: Gauss IOD, differential correction, ephemeris generation, Keplerian/equinoctial elements, JPL ephemerides, reference frames, residuals |

A typical workflow:

1. Use **photom** to load observations into an `ObsDataset`.
2. Call **outfit**'s `FitIOD::fit_iod` / `FitIOD::fit_full_iod` for initial orbit determination, or `FitLSQ::fit_lsq` for a full least-squares fit.
3. Inspect the returned `GaussResult` / `DifferentialCorrectionOutput` and RMS values.
4. Optionally, call `OrbitalElements::compute` with an `EphemerisRequest` to generate predicted positions.

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

## Differential Orbit Correction

Starting from an initial orbit (e.g. from the Gauss IOD step), outfit can refine the solution via **weighted least-squares differential correction**:

1. Compute predicted `(RA, Dec)` for each observation using the current elements.
2. Form the **design matrix** of partial derivatives of the predicted coordinates with respect to the six equinoctial elements.
3. Solve the **normal equations** to obtain the element correction vector δx.
4. Apply the correction, check convergence (`‖δx‖ < threshold`), and iterate (Newton–Raphson).
5. After convergence, apply **projection-based outlier rejection**: observations whose chi-squared contribution exceeds the configured threshold are flagged and the inner loop is re-run on the reduced set. Iterate until the selection is stable.
6. Rescale the **covariance matrix** by the posterior normalised RMS.

The `FitLSQ` trait exposes this pipeline at the dataset level:

```rust,no_run
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use hifitime::ut1::Ut1Provider;
use rand::{rngs::StdRng, SeedableRng};
use outfit::differential_orbit_correction::{DifferentialCorrectionConfig, FitLSQ};
use outfit::jpl_ephem::{download_jpl_file::EphemFileSource, JPLEphem};
use outfit::IODParams;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let dataset = ObsDataset::from_mpc_80_col("tests/data/2015AB.obs")?;

    let jpl_source: EphemFileSource = "horizon:DE440".try_into()?;
    let jpl = JPLEphem::new(&jpl_source)?;
    let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long")?;

    let iod_params = IODParams::builder().build()?;
    let dc_config = DifferentialCorrectionConfig::default();

    let mut rng = StdRng::seed_from_u64(42);

    // Run IOD + differential correction for every trajectory.
    let results = dataset.fit_lsq(
        &jpl,
        &ut1,
        ObsErrorModel::FCCT14,
        &iod_params,
        &dc_config,
        None,   // no pre-computed initial orbits
        &mut rng,
    )?;

    for (traj_id, res) in &results {
        match res {
            Ok(fit) => println!("{traj_id} → normalised RMS = {:.4}", fit.orbit_quality()),
            Err(e)  => eprintln!("{traj_id} → error: {e}"),
        }
    }
    Ok(())
}
```

---

## Uncertainty Propagation

Outfit tracks orbital uncertainties throughout the full pipeline — from the least-squares fit through to element representation conversions.

### Generation from differential correction

The covariance matrix is produced directly by the `FitLSQ` pipeline. At the end of each Newton–Raphson convergence, `solve_weighted_least_squares` builds the **normal matrix** G⊤WG and inverts it (Cholesky, with QR fallback) to obtain the **6×6 covariance matrix** Γ = (G⊤WG)⁻¹ in equinoctial element space:

```
Γ[j,k] = (G⊤WG)⁻¹[j,k]    (6×6, equinoctial basis)
```

This raw covariance is then **rescaled** by a posterior inflation factor μ that accounts for the degrees of freedom and the quality of the fit:

- if normalised RMS ≤ 1 (good fit): μ = √(n_measurements / (n_measurements − n_free))
- if normalised RMS > 1 (poor fit): μ = normalised_rms × √(n_measurements / (n_measurements − n_free))

The final covariance is stored in `DifferentialCorrectionOutput::uncertainty` and embedded in the returned `OrbitalElements::Equinoctial` variant.

### Propagation between element representations

Each `OrbitalElements` variant carries an optional `covariance: Option<OrbitalCovariance>` and an optional `uncertainty: Option<*Uncertainty>` (per-element 1-σ standard deviations). When converting between representations, the covariance is propagated via **first-order linear (Jacobian) propagation**:

```
Σ_y = J · Σ_x · Jᵀ
```

where J = ∂y/∂x is the 6×6 Jacobian of the transformation evaluated at the nominal elements. The following Jacobians are computed analytically:

| Conversion | Jacobian |
|---|---|
| Keplerian → Equinoctial | ∂(a,h,k,p,q,λ) / ∂(a,e,i,Ω,ω,M) |
| Equinoctial → Keplerian | ∂(a,e,i,Ω,ω,M) / ∂(a,h,k,p,q,λ) |
| Cometary → Keplerian    | ∂(a,e,i,Ω,ω,M) / ∂(q,e,i,Ω,ω,ν) |
| Cometary → Equinoctial  | chain rule via Keplerian             |

Conversions near singularities (e → 0 for Keplerian, i → 0 for equatorial) set undefined partial derivatives to zero; **equinoctial elements** are preferred for nearly circular or equatorial orbits to avoid these singular regions.

```rust,no_run
use outfit::orbit_type::{OrbitalElements, keplerian_element::KeplerianElements};
use outfit::orbit_type::uncertainty::{KeplerianUncertainty, OrbitalCovariance};
use nalgebra::Matrix6;

// Keplerian elements with a diagonal covariance (simplified example).
let cov = OrbitalCovariance { matrix: Matrix6::from_diagonal_element(1e-6) };
let oe = OrbitalElements::Keplerian {
    elements: kep,
    uncertainty: Some(KeplerianUncertainty::from_covariance(&cov)),
    covariance: Some(cov),
};

// Convert to equinoctial — covariance is automatically propagated via Jacobian.
let oe_eq = oe.to_equinoctial().unwrap();

if let OrbitalElements::Equinoctial { uncertainty, .. } = oe_eq {
    if let Some(unc) = uncertainty {
        println!("σ(h) = {:.2e}", unc.eccentricity_sin_lon);
    }
}
```

---

## Ephemeris Generation

Given any `OrbitalElements`, outfit can predict the apparent sky position and geometric state of the body as seen from one or more observers:

```rust,no_run
use outfit::{
    OrbitalElements,
    ephemeris::{EphemerisConfig, EphemerisMode, EphemerisRequest, request::Combined},
};
use hifitime::{Epoch, Duration};

// `elements` obtained from IOD or differential correction.
let result = elements.compute(
    &EphemerisRequest::<Combined>::new(EphemerisConfig::default())
        .add(observer, EphemerisMode::Range {
            start: Epoch::from_mjd_tt(60310.0),
            end:   Epoch::from_mjd_tt(60340.0),
            step:  Duration::from_days(1.0),
        }),
    &jpl,
    &ut1,
);

for entry in result.successes() {
    let (pos, geo) = entry.result.as_ref().unwrap();
    println!(
        "{}: RA={:.4} Dec={:.4}  phase={:.2}°  elongation={:.2}°",
        entry.epoch, pos.coord.ra, pos.coord.dec,
        geo.phase_angle.to_degrees(), geo.solar_elongation.to_degrees(),
    );
}
```

The pipeline converts elements to **equinoctial form**, propagates with the selected propagator (two-body or N-body), rotates from ecliptic to equatorial J2000, and applies aberration correction. Errors at individual epochs are stored in the result rather than aborting the computation.

---

## Performance & Reproducibility

- **Deterministic runs** by default (set RNG seeds explicitly when noise is used).
- **Precompute observer positions** to avoid ephemeris I/O in hot paths.
- Compile with `--release` for production.
- Keep ephemerides cached locally to avoid I/O stalls on repeated runs.
- For large trajectory sets, tune `IODParams.batch_size` (with `parallel`) so each batch fits comfortably in cache (start with 4–8 and benchmark).
- Handle errors explicitly: a non-finite RMS surfaces as `OutfitError::NonFiniteScore`.

---


## Documentation

To compile the documentation locally, run the following command in the terminal:
```bash
RUSTDOCFLAGS="--html-in-header $(pwd)/katex-header.html" cargo doc --no-deps --all-features
```

---

## Roadmap

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
