# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.1.0] - 2025-09-08

### Added
- **BiMap iteration API**:
  - `iter()` → `(&K, &V)` over the forward map
  - `iter_rev()` → `(&V, &K)` over the reverse map
  - `keys()`, `values()` convenience views
  - `IntoIterator` impls for `&BiMap`, `&mut BiMap`, and `BiMap` (owned)
- **Observer inversion**:
  - `Observer::geodetic_lat_height_wgs84()` — recover **geodetic latitude (deg)** and **ellipsoidal height (m)** from parallax constants `(ρ·cosφ, ρ·sinφ)` using Bowring’s method
- **Observatories pretty-printing**:
  - `impl Display for Observatories` — prints *User-defined observers* first, then *MPC observers* (if available); longitude/latitude in **degrees**, elevation in **kilometers**
  - `Outfit::show_observatories()` — zero-allocation display adaptor for pretty-printing
  - `Outfit::show_observatories_string()` — convenience method returning a formatted `String`
  - `Outfit::write_observatories<W: io::Write>()` — write formatted listing to any writer
- **Observer pretty-printing**:
  - `impl Display for Observer` with **compact** one-line default (`{}`) and **verbose** multi-line **alternate** formatting (`{:#}`) including parallax constants, geocentric latitude `φ_geo`, geocentric distance `ρ` (RE), and optional RA/DEC 1-σ accuracies (arcsec)

### Changed
- Documentation expanded and clarified:
  - Docstrings for the inversion routine (units, ellipsoid, numerical notes)
  - Docstrings for `ObservatoriesView` and the Outfit display helpers, including usage examples and “See also” sections
- Display output now explicitly states units and section ordering; internal indices remain hidden from users

### Fixed
- **Critical bug – `parquet_to_trajset` RA/DEC uncertainties**:
  - Observational uncertainties for RA/DEC provided in **arcseconds** were being treated as **radians**.  
    The routine now **converts arcseconds → radians** (`ArcSec::to_radians()`) once per file before constructing `Observation`s, ensuring correct weighting in the fit and realistic residuals.
- **Documentation**: corrected `Observer` docs — elevation is in **meters** (previously incorrectly stated as kilometers). APIs and internal computations were already in meters; this is a documentation fix only (no breaking change).
- Added unit tests for:
  - BiMap iteration symmetry (`iter` vs `iter_rev`)
  - Observatories `Display` output (presence of headers, user/MPC sections)
  - Outfit display helpers (`show_observatories` vs `show_observatories_string` parity)
  - Geodetic inversion round-trips on real sites (Haleakalā, Mauna Kea, Paranal, Cerro Pachón, La Silla, Kitt Peak, Roque de los Muchachos), with tolerant assertions via `approx`

### Notes
- Iteration order for `BiMap`/`Observatories` is hash-based and not guaranteed to be stable; switch to `indexmap` if deterministic insertion order is required.

---

## [1.0.0] - 2025-09-04

### Added
- Extensive new tests across the codebase, including:
  - Property-based testing (proptest) for critical functions: ephemeris computation, velocity correction, universal Kepler solver, and orbit determination
  - Unit and integration tests for observation parsing (80-column, ADES, Parquet), JPL ephemeris, and error models
  - Test data files for real and synthetic observations
- Support for RMS (root mean square) error estimation for 80-column files, with automatic error parsing from MPC format
- JPL ephemeris file reading:
  - Support for both HORIZON and NAIF/SPICE binary formats
  - Interpolation and extraction of solar system body positions
- Parquet file reading for batch observations, with optimized batch loading and error assignment
- ADES (XML) file reading and parsing, supporting both structured and flat ADES formats
- Observer management:
  - Build observers from MPC code or custom coordinates
  - New observer struct and observatory registry
- CI checks for clippy and fmt, and pre-commit hooks for code quality
- Batch processing and benchmarking infrastructure using Criterion
- Equinoctial elements and conversion utilities:
  - Conversion between Keplerian and equinoctial elements
  - Two-body solver and Jacobian computation for equinoctial elements
- Progress indicators (feature-gated): optional real-time progress bars with `indicatif` and loop performance measurement via `IterTimer`
- Extensive documentation:
  - Unified docstrings in English with *Arguments*, *Return*, and *See also* sections
  - Examples integrated directly in Rustdoc (e.g. with `cfg(feature = "jpl-download")`)
- Examples folder: runnable end-to-end examples (`examples/parquet_to_orbit.rs`, `examples/iod.rs`) to demonstrate usage

### Changed
- Major refactor of `kepler.rs` (Keplerian motion, universal variable solver, anomaly conversion), `orb_elem.rs`, `equinoctial_element.rs`, and `initial_orbit_determination/gauss.rs` for clarity, maintainability, and performance
- Refactor and modularization of the orbit determination pipeline: moved and split code into `initial_orbit_determination/`, `keplerian_element.rs`, and `equinoctial_element.rs`
- Refactor and modularization of the observation handling: split `observations.rs` into `observations/mod.rs`, `observations/ades_reader.rs`, `observations/parquet_reader.rs`, `observations/trajectory_ext.rs`, and `observations/observations_ext.rs`
- Refactor and modularization of observer management: new `observers/` folder with `mod.rs`, `bimap.rs`, and `observatories.rs`
- Refactor and improvement of reference frame and rotation utilities in `ref_system.rs`
- Refactor and improvement of error models in `error_models/` and `outfit_errors.rs`
- Refactor and improvement of JPL ephemeris handling: new `jpl_ephem/` folder with submodules for HORIZON and NAIF file support
- Refactor and improve documentation throughout the codebase
- Refactor and optimize triplet selection and Gauss IOD routines
- Refactor and optimize use of nalgebra for vector/matrix operations
- Refactor for better performance and code clarity
- Update README and license (now CeCILL-C)
- Make almost everything public for documentation purposes
- Benchmarks: improved Criterion benches, especially for `prelim_orbit`
- Error handling: switched fully to [`thiserror`] with `#[from]` conversions and improved error messages
- Profiles: added optimized `[profile.bench]` with `required-features = ["jpl-download"]`

### Fixed
- Fix and improve documentation links
- Fix unwraps in production code
- Improve error reporting and handling
- Test infrastructure: fixed JPL-dependent tests by clarifying ephemeris file requirements and adding fallback
- Cargo manifest: audited and removed unused dependencies (via `cargo-udeps`)

### Removed
- Remove dead and unused code

### Security
- No security fixes in this release

---

For a full commit history, see the repository log between `main` and `release-1.0.0`.
