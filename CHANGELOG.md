# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [2.1.0] - 2025-09-18

### Added
- **Observations tabular display (`observations::display`)**
  - New **display adaptor** `ObservationsDisplay` with ergonomic builders:
    - `observations.show()` – compact, fixed-width table,
    - `observations.table_wide()` – diagnostic table,
    - `observations.table_iso()` – timestamp-centric table.
  - **Sorting:** `.sorted()` prints rows **sorted by epoch** (MJD TT asc.); the `#` column always shows the **original index**.
  - **Precision knobs:** `.with_seconds_precision(p)` (sexagesimal + ISO seconds) and `.with_distance_precision(p)` (AU columns in wide mode).
  - **Site labels:** `.with_env(&Outfit)` resolves observer **names**; renders `"Name (#id)"` when available, or falls back to the numeric site id.
  - **Wide/ISO layouts** now use **`comfy-table`** for neat borders and alignment; compact layout keeps a lightweight fixed-width formatter.
  - Convenience: `observations.show_string()` returns an owned string in compact mode.

- **Observation pretty-printing**
  - `impl Display for Observation`:
    - **Compact** one-liner (`{}`):  
      `Obs(site=…, MJD=… TT, RA=HHhMMmSS.sss ± σ", DEC=±DD°MM’SS.sss" ± σ", r_geo=[…] AU, r_hel=[…] AU)`
    - **Pretty** multi-line **alternate** (`{:#}`): adds **ISO TT** and **ISO UTC** timestamps (via `hifitime`), aligned sections, explicit units.

- **Sexagesimal & formatting helpers**
  - `time::fmt_ss(seconds, prec) -> String` — force `"SS.sss"` with a **two-digit** integer part.
  - `conversion::ra_hms_prec(rad, prec) -> (u32, u32, f64)` — RA (rad) → `(HH, MM, SS.sss)` with rounding & carry.
  - `conversion::dec_sdms_prec(rad, prec) -> (char, u32, u32, f64)` — DEC (rad) → `(sign, DD, MM, SS.sss)` with rounding & carry, clamped at **±90°**.
  - `conversion::fmt_vec3_au(&Vector3, prec) -> String` — `[ x, y, z ] AU` in fixed-point with `prec` decimals.

- **`hifitime` ISO renderers**
  - `time::iso_tt_from_epoch(epoch_tt, sec_prec) -> String` — `YYYY-MM-DDThh:mm:SS.sss TT`.
  - `time::iso_utc_from_epoch(epoch_tt, sec_prec) -> String` — `YYYY-MM-DDThh:mm:SS.sssZ`.

- **ObjectNumber conversions (ergonomics)**
  - `From<u8/u16/u32>`, `From<&u32>`, `From<String>`, `From<&String>`, `From<&str>`, `From<Cow<'_, str>>`,
  - `TryFrom<i64/i32/i16/i8>` fallible integer conversions.
  - `Display` renders either the integer or the stored string form.

- **Tests**
  - Added unit tests for **Observations** display:
    - default headers & fixed-width formatting,
    - ISO mode headers and `Z` suffix,
    - wide mode headers, radians & AU distances,
    - `sorted()` ordering while preserving original index in `#`,
    - seconds & distance precision knobs,
    - negative DEC sign preservation.
  - Extended tests for `Observation` display (compact & pretty), RA wrap to 24h, DEC sign, AU vector precision, arcsec uncertainties.

### Changed
- **Wide/ISO rendering** now uses **`comfy-table`** (keeps compact/default mode as fixed-width). This removes manual spacing bugs and improves readability.
- **Sorting builder API:** `ObservationsDisplay::sorted(true)` ➜ **`ObservationsDisplay::sorted()`**.  
  *Migration:* replace `.sorted(true)` with `.sorted()`; remove `.sorted(false)` (default is unsorted).
- **Seconds formatting** is strictly zero-padded: `SS.sss` (e.g., `00.000`).
- **Carry-safe rounding** for RA/DEC seconds (e.g., `59.9995s` → `60.000s` then carried to min/hr or min/deg).
- **Docs** across modules expanded with **Units**, **See also**, **Notes**, and runnable **Examples**.

### Fixed
- **Table alignment**: wide tables had header/data misalignment; `comfy-table` fixes this and ensures right-aligned numerics.
- **Consistent sexagesimal output**: negative DEC strings keep the leading sign and zero-padded seconds.
- Edge cases near RA=24h and DEC=±90° now produce canonical outputs (`00h00m00.000s`, `±90°00'00.000"`).

### Documentation
- **Top-level docs updated**:
  - `observations::display`: overview, layouts (Default/Wide/ISO), sorting, precision, `with_env`, performance notes, and multiple examples.
  - `conversion.rs`: parsing overview, accuracy estimation from decimals, sexagesimal utilities, vector formatting.
  - `time.rs`: calendar/JD/MJD conversions, ISO TT/UTC rendering via `hifitime`, batch scale changes (UTC→TT), and GMST.
- Function-level docstrings enriched with explicit **Arguments/Returns** and **See also** links.

### Dependencies
- New runtime dependency: **`comfy-table`** (actively maintained; used for Wide/ISO layouts).

### Notes
- No behavior-breaking changes in numeric results; **string formatting** is stricter (zero-padding & carry), which may affect golden-file tests and logs.

---

## [2.0.0] - 2025-09-17

### Added
- **Trajectories module split (new folder `trajectories/`)**  
  Public submodules and their roles:
  - `trajectories::batch_reader` — `ObservationBatch` (zero-copy or owned), `observation_from_batch` (crate-private).
  - `trajectories::trajectory_file` — public ingestion trait with `new_from_*` / `add_from_*` entry points.
  - `trajectories::trajectory_fit` — batch Gauss IOD (`TrajectoryFit`), results map, stats, helpers.
  - `trajectories::parquet_reader` — high-throughput Arrow/Parquet ingestion with projection & epoch cache.
  - `trajectories::mpc_80col_reader` — MPC 80-column reader with precise field slicing and CCD filtering.
- **Batch IOD APIs on `TrajectorySet`** (`TrajectoryFit`):
  - `estimate_all_orbits(&mut self, &Outfit, &mut impl Rng, &IODParams) -> FullOrbitResult`
  - `estimate_all_orbits_with_cancel(&mut self, &Outfit, &mut impl Rng, &IODParams, should_cancel: impl FnMut() -> bool) -> FullOrbitResult`
  - `total_observations(&self) -> usize` — total samples across all trajectories.
  - `obs_count_stats(&self) -> Option<ObsCountStats>` — min/p25/median/p95/max of per-trajectory counts.
- **Parallel batch IOD** *(new feature flag `parallel`)*  
  - New dependency: [`rayon`](https://crates.io/crates/rayon) (optional).  
  - `estimate_all_orbits_in_batches_parallel(&mut self, &Outfit, &mut impl Rng, &IODParams, batch_size: usize) -> FullOrbitResult`  
    Runs Initial Orbit Determination (IOD) in **parallel across batches**, each processed sequentially for locality.  
    - Seeds are deterministically derived per object from a global random base seed.
    - Observations are moved by batch, avoiding clones.
    - Supports global, thread-safe progress bar under `progress` feature.
- **Parallel IOD parameter** in [`IODParams`]:
  - New field `batch_size` (active only with feature `parallel`).
  - Controls the number of trajectory batches to process concurrently
    when multi-threading is enabled via Rayon (`batch_size > 1`).
  - Default value: **4**.
- **Result helpers & types**
  - `type FullOrbitResult = HashMap<ObjectNumber, Result<(GaussResult, f64), OutfitError>, RandomState>`
  - `gauss_result_for(&FullOrbitResult, &ObjectNumber)` — borrow solution & RMS if present.
  - `take_gauss_result(&mut FullOrbitResult, &ObjectNumber)` — move solution & RMS out of the map.
  - `ObsCountStats` + `Display` (compact `{}` / pretty `{:#}`).
- **Progress feature**
  - With `--features progress`: live progress bar via `indicatif` + moving-average timing, and a timer-based cooperative cancel in `estimate_all_orbits_with_cancel`.
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
- **TrajectorySet convenience**:
  - `number_of_trajectories(&self) -> usize` — returns the total number of trajectories
    stored in the set (simple wrapper around `len()`).
- **Error handling**:
  - New `OutfitError::NonFiniteScore(f64)` variant — raised when a non-finite RMS
    value (`NaN`, `±∞`) is produced during orbit scoring in `estimate_best_orbit`.
    Provides the offending floating-point value for diagnostics.

### Changed
- **Folder & paths**  
  Trajectory logic moved out of `observations` into dedicated `trajectories/`. Public docs and examples now reference:
  - `crate::trajectories::batch_reader::ObservationBatch`
  - `crate::trajectories::trajectory_file::TrajectoryFile`
  - `crate::trajectories::trajectory_fit::{TrajectoryFit, FullOrbitResult, ObsCountStats}`
- **Documentation (top-level & modules)**  
  Added/overhauled top-level docs for:
  - `trajectories/mod.rs` — overview, data model, ingestion sources, units/time scales, feature flags, quick-start.
  - `trajectories/trajectory_file.rs` — public ingestion surface (80-col, Parquet, ADES, in-memory batches), error semantics.
  - `trajectories/parquet_reader.rs` — expected schema, null policy, projection & cache, performance notes.
  - `trajectories/mpc_80col_reader.rs` — field map, CCD filter, RA/Dec/time parsing, uncertainty handling.
  - `trajectories/tra
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
- **`estimate_best_orbit` robustness**:
  - Now explicitly rejects **non-finite RMS scores** (`NaN`, `±∞`) during orbit scoring.
  - Such cases are recorded as `OutfitError::NonFiniteScore(value)` and skipped,
    instead of being considered as valid candidates.
  - This prevents spurious “best orbit” selections when numerical instabilities
    produce invalid residuals.
- **MPC 80-column parsing docs & weighting clarity**  
  Docstrings now explicitly document the `cos δ` factor on RA uncertainties, TT/MJD handling, and the CCD-only policy; errors surface as `OutfitError::Parsing80ColumnFileError`.



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
- **`estimate_best_orbit` now returns a *definite* solution on `Ok`**
  This cascades to batch outcomes: `FullOrbitResult` holds `Result<(GaussResult, f64), OutfitError>` (no `Option`).
  - **Migration**:
    ```rust
    // Before:
    match results.get(&obj) {
        Some(Ok((Some(g), rms))) => { /* use */ }
        Some(Ok((None, _)))      => { /* no solution */ }
        Some(Err(e))            => { /* error */ }
        _ => {}
    }

    // After:
    match results.get(&obj) {
        Some(Ok((g, rms))) => { /* always a valid orbit */ }
        Some(Err(e))       => { /* error (includes no-solution cases) */ }
        _ => {}
    }
    ```
- **Ingestion from in-memory batches**  
  `TrajectoryFile::{new_from_vec, add_from_vec}` switched to `Result<…>` and derive object IDs from `ObservationBatch::trajectory_id`.  
  `observation_from_batch` now appends **all** trajectories found in the batch into a provided `&mut TrajectorySet`.


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
