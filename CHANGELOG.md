# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.0.0] - 2025-09-10

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
- **Angle conversions**:
  - `conversion::arcsec_to_rad(f64) -> f64` — utility to convert **arcseconds → radians**, implemented via `.to_radians()` for numerical consistency
- **ObjectNumber ordering traits**:
  - `ObjectNumber` now derives `PartialOrd` and `Ord` (in addition to `Debug`, `Clone`, `PartialEq`, `Eq`, `Hash`):
    ```rust
    #[derive(Debug, Clone, PartialEq, Eq, Hash, PartialOrd, Ord)]
    pub enum ObjectNumber {
        /// Integer-based MPC designation (e.g. 1, 433…)
        Int(u32),
        /// String-based designation (provisional, comet, etc.)
        String(String),
    }
    ```
  - This allows deterministic sorting of objects (useful for stable iteration in tests & reports).

### Changed
- Documentation expanded and clarified:
  - Docstrings for the inversion routine (units, ellipsoid, numerical notes)
  - Docstrings for `ObservatoriesView` and the Outfit display helpers, including usage examples and “See also” sections
- Display output now explicitly states units and section ordering; internal indices remain hidden from users
- **ObservationBatch now stores angles in radians and supports zero-copy or owned buffers**:
  - Representation switched to **`Cow<'a, [f64]>`** for `ra`, `dec`, and `time`
  - New constructors:
    - `ObservationBatch::from_radians_borrowed(&[f64], &[f64], err_ra_rad, err_dec_rad, &[MJD])`
    - `ObservationBatch::from_degrees_owned(&[f64], &[f64], err_ra_arcsec, err_dec_arcsec, &[MJD])`
  - Internally, RA/DEC and their uncertainties are **radians**; degree/arcsec inputs are converted **once** at construction
  - **Migration note** (source-compatible pattern):
    - **Before**
      ```rust
      let batch = ObservationBatch {
          ra: &ra_deg,
          error_ra: 0.5,           // arcsec
          dec: &dec_deg,
          error_dec: 0.5,          // arcsec
          time: &time_mjd,
      };
      ```
    - **After** (choose one)
      ```rust
      // If your inputs are in degrees / arcseconds:
      let batch = ObservationBatch::from_degrees_owned(&ra_deg, &dec_deg, 0.5, 0.5, &time_mjd);

      // If your inputs are already radians:
      let batch = ObservationBatch::from_radians_borrowed(&ra_rad, &dec_rad, err_ra_rad, err_dec_rad, &time_mjd);
      ```

- **`observation_from_batch` now fills the entire `TrajectorySet`**:
  - Previously returned a single `Observations` vector for one trajectory.
  - Now takes a `&mut TrajectorySet` and **inserts observations for all trajectory IDs found in the batch** (grouping by `trajectory_id`).
  - See **Breaking Changes** for details and migration snippet.
- **`observation_from_batch` performance/clarity improvements** (parity with `parquet_to_trajset`):
  - Resolve observer → compact `u16` **once** before the loop
  - Pre-fetch UT1 provider once
  - Introduce an **epoch-position cache**: unique MJD(TT) → `(geo_pos, helio_pos)` so positions are computed once per timestamp
  - Hot loop no longer performs angle unit conversions; it assumes batch angles/uncertainties are already **radians** (as guaranteed by `ObservationBatch` constructors)

- **`estimate_best_orbit` signature changed**:
  - **Before**
    ```rust
    fn estimate_best_orbit(
        &mut self,
        state: &Outfit,
        error_model: &ErrorModel,
        rng: &mut impl rand::Rng,
        params: &IODParams,
    ) -> Result<(Option<GaussResult>, f64), OutfitError>
    ```
  - **After**
    ```rust
    fn estimate_best_orbit(
        &mut self,
        state: &Outfit,
        error_model: &ErrorModel,
        rng: &mut impl rand::Rng,
        params: &IODParams,
    ) -> Result<(GaussResult, f64), OutfitError>
    ```
  - **Migration note**:
    - Previously, you had to handle `Ok((None, rms))` as a “no solution” case.
    - Now, a successful `Ok` **always** contains a valid `GaussResult`.  
      Absence of a solution is reported as an `Err(OutfitError)`.
    - Typical update:
      ```rust
      // Before
      match obs.estimate_best_orbit(&env, &error_model, &mut rng, &params) {
          Ok((Some(gauss), rms)) => { /* use orbit */ }
          Ok((None, rms)) => { /* no solution */ }
          Err(e) => { /* error */ }
      }

      // After
      match obs.estimate_best_orbit(&env, &error_model, &mut rng, &params) {
          Ok((gauss, rms)) => { /* always valid orbit */ }
          Err(e) => { /* error (includes no-solution case) */ }
      }
      ```


### Fixed
- **Critical bug – `parquet_to_trajset` RA/DEC uncertainties**:
  - Observational uncertainties for RA/DEC provided in **arcseconds** were being treated as **radians**.  
    The routine now **converts arcseconds → radians** (`ArcSec::to_radians()`) once per file before constructing `Observation`s, ensuring correct weighting in the fit and realistic residuals.
- **Documentation**: corrected `Observer` docs — elevation is in **meters** (previously incorrectly stated as kilometers). APIs and internal computations were already in meters; this is a documentation fix only (no breaking change).
- Added unit/integration tests for:
  - BiMap iteration symmetry (`iter` vs `iter_rev`)
  - Observatories `Display` output (presence of headers, user/MPC sections)
  - Outfit display helpers (`show_observatories` vs `show_observatories_string` parity)
  - Geodetic inversion round-trips on real sites (Haleakalā, Mauna Kea, Paranal, Cerro Pachón, La Silla, Kitt Peak, Roque de los Muchachos), with tolerant assertions via `approx`
  - **Integration test `tests/vec_to_iod.rs`** — end-to-end pipeline:
    `Vec<deg/arcsec> → ObservationBatch::from_degrees_owned → TrajectorySet → estimate_all_orbits (Gauss IOD) → gauss_result_for → KeplerianElements`.  
    Uses `approx` with angle wrap-around, seeded RNG (`StdRng::seed_from_u64(42)`) for reproducibility, and `jd_to_mjd` for epoch conversion.

### Breaking Changes
- **`TrajectorySet::new_from_vec` and `add_from_vec` signatures changed**:
  - **Before**
    ```rust
    fn new_from_vec(
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Self;

    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        object_number: &str,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    );
    ```
  - **After**
    ```rust
    fn new_from_vec(
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<Self, OutfitError>
    where
        Self: Sized;

    fn add_from_vec(
        &mut self,
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<(), OutfitError>;
    ```
  - **Migration note**:
    - The `object_number` parameter was removed; the object identity is now taken from `batch.trajectory_id`.
    - Error handling switched to `Result<..., OutfitError>`.
    - Typical update:
      ```rust
      // Before
      let mut ts = TrajectorySet::new_from_vec(&mut env, "33803", &batch, observer);

      // After
      let mut ts = TrajectorySet::new_from_vec(&mut env, &batch, observer)?;
      ```

- **`observation_from_batch` now populates all trajectories**:
  - **Before**: returned `Observations` for a single trajectory.
    ```rust
    pub(crate) fn observation_from_batch(
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Observations
    ```
  - **After**: returns `Result<(), OutfitError>` and appends observations directly into a `TrajectorySet`, grouping by `batch.trajectory_id`.
    ```rust
    pub(crate) fn observation_from_batch(
        trajectories: &mut TrajectorySet,
        env_state: &mut Outfit,
        batch: &ObservationBatch<'_>,
        observer: Arc<Observer>,
    ) -> Result<(), OutfitError>
    ```
  - **Migration note**:
    - If you previously did:
      ```rust
      let obs = observation_from_batch(&mut env, &batch, observer);
      ts.insert(object_number, obs);
      ```
    - You should now do:
      ```rust
      observation_from_batch(&mut ts, &mut env, &batch, observer)?;
      ```

### Notes
- Iteration order for `BiMap`/`Observatories` is hash-based and not guaranteed to be stable; switch to `indexmap` if deterministic insertion order is required.
- With `ObjectNumber: Ord`, you can now sort keys for stable iteration in tests and reports:
  ```rust
  let mut keys: Vec<_> = traj_set.keys().cloned().collect();
  keys.sort();

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
