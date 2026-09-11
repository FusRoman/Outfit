# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- **Alternative planetary ephemeris backend via ANISE (`ephem-anise` feature)**
  - The planetary ephemeris backend is now chosen at compile time by Cargo
    feature. `ephem-builtin` (in the `default` set) keeps the in-house JPL DE /
    NAIF SPK reader; `ephem-anise` swaps in the [`anise`](https://crates.io/crates/anise)
    toolkit, a Rust reimplementation of the NAIF SPICE toolkit validated against
    SPICE to machine precision. This is motivated by the still-open velocity-chain
    bug in the in-house NAIF reader (see the ignored regression test
    `naif_body_ephemeris_velocity_is_the_position_derivative`): the ANISE backend
    computes correct velocities, and a new regression test
    (`body_velocity_matches_central_difference`) locks that in.
  - New constructors on `JPLEphem`: `from_builtin` and `from_anise`. `JPLEphem::new`
    is unchanged and selects the backend from the enabled features — the ANISE
    backend when both are on; use `from_builtin` to force the in-house reader.
    Building with neither `ephem-builtin` nor `ephem-anise` is a compile error.
  - The ANISE backend reads NAIF SPK kernels only; a `horizon:` (legacy DE binary)
    source returns `OutfitError::InvalidJPLEphemFileSource`. It reuses Outfit's
    existing ephemeris downloader and on-disk cache (no `anise/metaload`).
  - New error variant `OutfitError::AniseEphemerisError(String)` (feature-gated).
  - `earth_ephemeris` under ANISE returns the Earth **geocenter** relative to the
    Sun, matching the built-in Horizon reader (the built-in NAIF reader returns
    the Earth–Moon barycenter relative to the Solar System Barycenter). Parity
    tests check ANISE against the Horizon reader to `1e-8` AU on planet positions
    and `1e-9` AU / `1e-10` AU·day⁻¹ on Earth.
  - `hifitime` requirement bumped to `4.3` (shared, unchanged epoch type). The
    `nalgebra` requirement is unchanged; the ANISE state is converted at the
    backend boundary.

- **Main-belt asteroids as N-body perturbers (ANISE backend)**
  - `NaifIds` gained an `AST(AsteroidNumber)` variant identifying a numbered
    main-belt asteroid by its official minor-planet number (e.g. `1` for
    Ceres); `AsteroidNumber::CERES` / `PALLAS` / `VESTA` are provided for
    convenience.
  - New constructor `JPLEphem::from_anise_with_main_belt_asteroids` loads DE440
    plus a supplementary NAIF SPK kernel
    (`codes_300ast_20100725.bsp`, downloaded and cached like the primary
    kernel) covering 300 numbered asteroids. Once loaded, any of them can be
    added to `NBodyConfig::perturbing_bodies` like any other body.
  - `propagator::planet_gm::known_main_belt_asteroids()` returns all 300 as
    ready-to-use `NaifIds::AST(_)` perturbers, each with a gravitational
    parameter from the kernel's own published mass table — pick as many or as
    few as needed.
  - This is ANISE-specific: the in-house reader has no code path for
    supplementary small-body kernels, so `NaifIds::AST(_)` only resolves under
    the `ephem-anise` backend.

### Changed

- **Planetary ephemeris readers are behind `feature = "ephem-builtin"`** (breaking)
  - The `jpl_ephem::horizon` and `jpl_ephem::naif` reader submodules,
    `JPLEphem::HorizonFile` / `JPLEphem::NaifFile`, and
    `JPLEphem::try_into_horizon` / `try_into_naif` now require the
    `ephem-builtin` feature, which is enabled by the `default` feature set. The
    lightweight `horizon_version`, `naif_version` and `naif_ids` modules stay
    available for every feature combination.

- **`JPLEphem::earth_ephemeris` / `body_ephemeris` — explicit output frame** (breaking)
  - Both methods now take an `EphemerisFrame` argument (`Equatorial` or
    `Ecliptic`). JPL DE and NAIF SPK data are stored in equatorial mean J2000
    (ICRF); `EphemerisFrame::Equatorial` returns it unchanged,
    `EphemerisFrame::Ecliptic` rotates the state to ecliptic mean J2000 with
    `ROT_EQUMJ2000_TO_ECLMJ2000`. `EphemerisFrame` is re-exported at the crate
    root and in the prelude.
  - Fixes an out-of-plane perturbation error: the N-body propagator integrates in
    ecliptic mean J2000 but `perturber_ephemeris` sampled `body_ephemeris`
    without rotating, so every perturber sat up to ~2 AU off the ecliptic. It now
    samples in `EphemerisFrame::Ecliptic`. The stray force coupled the fitted
    inclination and ascending node; the ephemeris residuals for the N-body
    reference objects drop accordingly (e.g. K09R05F per-site median from ~2.7″
    to ~0.2″).
  - The four `earth_ephemeris` call sites already treated the output as
    equatorial and now pass `EphemerisFrame::Equatorial` — no numerical change.
    The mislabeled "ecliptic J2000" docs on `body_ephemeris`,
    `HorizonData::ephemeris` and `NaifData::ephemeris` were corrected.
  - Reference values in `test_diff_cor_nbody_nonregression` were regenerated; the
    `tests/test_ephemeris.rs` `*_nbody` thresholds are now 0.3″ / 0.4″ / 0.5″.
    Added frame-consistency unit and property tests in `jpl_ephem` and a
    `perturber_positions_are_ecliptic` guard in `perturber_ephemeris`. The
    Sun-only propagation path is unchanged (rotation-invariant).

### Fixed

- **`NaifData::ephemeris` (built-in NAIF SPK backend) — panicked instead of
  returning an error**
  - Querying a `(target, center)` pair absent from the loaded kernel, or an
    epoch outside a segment's coverage, made `naif_data::NaifData::ephemeris`
    panic. `NaifIds::AST(_)` (main-belt asteroids) made this easy to trigger,
    since `de440.bsp`-style kernels never carry small bodies, but the bug is
    general — any unsupported body/epoch combination on this backend hit it.
    `ephemeris` now returns `Result<InterpResult, OutfitError>`
    (`OutfitError::EphemerisBodyNotSupported` on failure); the two
    `JPLEphem::NaifFile` call sites in `jpl_ephem::mod` propagate it through
    the existing `body_ephemeris` `Result`, matching how the built-in Horizon
    reader and the ANISE backend already report the same situation.
    `earth_ephemeris` keeps its infallible contract (Earth is always present
    in a real kernel) by panicking with a clear message on failure instead,
    consistent with the ANISE backend's documented behaviour for Earth. This
    is a breaking change to `NaifData::ephemeris`'s public signature.
  - Documented, on `NaifIds`, exactly which variant each backend can resolve
    and what happens when it can't (always an error, never a panic, except
    for the documented `earth_ephemeris` case above).

- **N-body propagator — inverted sign on the indirect perturbation term**
  - `indirect_acceleration` returned `+GM·r_p/|r_p|³` instead of
    `−GM·r_p/|r_p|³`. In the heliocentric equation of motion the indirect term
    subtracts the acceleration imparted by each perturber to the (non-inertial)
    Sun; with the wrong sign the perturbation was added twice in the in-plane
    direction, so enabling planets *degraded* the fit and the residual grew with
    arc length. Fixed the sign; module and function docs updated to match.
  - The direct term and the gravity gradient `∂a/∂r` are unaffected (the
    indirect term does not depend on the small-body position), so the state
    transition matrix and covariance are unchanged. The Sun-only path is
    unchanged bit-for-bit (the Sun's heliocentric position is identically zero,
    so its indirect term is skipped).
  - Added a `#[cfg(test)] mod tests` in `propagator::nbody` with unit and
    property-based tests: the net perturbation (direct + indirect) vanishes as
    the small body approaches the Sun, the indirect term is antiparallel to the
    perturber position, and the gravity gradient is untouched.
  - Reference values in `test_diff_cor_nbody_nonregression` were regenerated and
    the `*_nbody` thresholds in `tests/test_ephemeris.rs` were tightened to the
    residuals now reached (see also the frame fix below).
  - Also realigned the `direct_acceleration` argument doc with the code
    (`r_asteroid − r_perturber`).

- **`JPLEphem::body_ephemeris` — spurious ×86400 on Horizon velocity**
  - The Horizon branch multiplied the returned velocity by `86400.0` (comment
    "AU/s → AU/day"), but the Horizon backend already yields AU/day after
    `InterpResult::to_au()` — there is no per-second quantity anywhere in that
    path. The velocity was therefore 86400× too large. Removed the scaling so
    the Horizon branch matches `earth_ephemeris` and the documented AU/day
    contract; added regression tests (`body_ephemeris(EarthMoon)` vs
    `earth_ephemeris`, and velocity vs a central finite difference of the
    position). No production result changed: the only callers of
    `body_ephemeris` read the position and discard the velocity.
  - Also documented (ignored regression test
    `naif_body_ephemeris_velocity_is_the_position_derivative`): the **NAIF**
    velocity chain is independently wrong — `EphemerisRecord::interpolate`
    scales by `2.0 / radius` instead of `1.0 / radius`, and the NAIF branches of
    `body_ephemeris` / `earth_ephemeris` divide by 86400 where they must
    multiply. This is only reachable with the NAIF backend selected and is left
    for a dedicated fix (it needs the `test_record_interpolation` /
    `test_jpl_ephemeris` oracles regenerated).

- **N-body propagator — time-independent perturber positions**
  - `NBodyOde::diff` ignored its time argument: every perturbing body's
    heliocentric position was sampled once at `t0` (`build_perturber_snapshots`)
    and held fixed for the whole integration arc. The integrated force field was
    therefore static, which is only acceptable for very short arcs (the previous
    doc comment claimed ≲ 30 days). On the observation arcs actually used by the
    fitting pipeline (median ~90 days, up to ~500 days for the backward
    ground-truth propagation) this *degraded* the fit when planetary
    perturbations were enabled instead of improving it.
  - Perturber positions are now evaluated at the current integration epoch from
    a per-arc **piecewise Chebyshev interpolation**
    (`propagator::perturber_ephemeris`), built once per trajectory in
    `run_differential_correction` and reused by every observation propagation.
    The ephemeris-request path (`EphemerisConfig` / `apparent_position`) builds an
    equivalent table on the fly per call.
  - The Sun keeps an exact special case (identically-zero heliocentric position),
    so the Sun-only path is unchanged bit-for-bit; a new
    `tests/test_nbody_propagator.rs` pins it against the analytic two-body
    solution. Reference values in `test_diff_cor_nbody_nonregression` were
    regenerated (Jupiter is now time-dependent).
  - Backward propagation (negative spans) is fully supported: the interpolation
    window covers `[min(t0, t1), max(t0, t1)]` plus a small margin.
  - Fixed an inverted doc comment on `direct_acceleration` (argument is
    `r_asteroid − r_perturber`).

### Added

- **`NBodyConfig::perturber_interp_degree` / `perturber_panel_days`**
  - Tuning knobs for the perturber-position interpolation: Chebyshev degree per
    panel (default `12`) and maximum panel length in days (default `16.0`, which
    keeps the worst body — Mercury — below `1e-10` AU).
- **`propagator::perturber_ephemeris`** module
  - `PerturberEphemerisSet` / `PerturberEphemeris`: per-arc Chebyshev
    interpolation of heliocentric perturber positions, with pure, individually
    tested helpers (Chebyshev–Gauss nodes, series fit, Clenshaw evaluation, panel
    layout/selection) plus property-based tests and an ephemeris oracle test.

## [4.1.0] - 2026-07-20

### Fixed

- **Universal Kepler solver — `alpha`/`mu` scaling bug**
  - `UniversalKeplerParams::alpha` was computed as the raw vis-viva energy
    `2E = v² − 2μ/r`, but the Stumpff engine, the Kepler residual/derivative
    equations, the Lagrange f/g coefficients, and the `prelim_elliptic` /
    `prelim_hyperbolic` / `prelim_parabolic` initial-guess heuristics all
    require the reciprocal-semi-major-axis convention `alpha = -1/a`
    (`= 2E / μ`). The missing division by `μ` made `propagate_universal`
    fail to converge — or converge to a spurious multi-revolution root — on
    real Kalman-filter state vectors (fink-fat), especially for gaps of
    several weeks or more.
  - Fixed across `propagation.rs`, `newton_solver.rs`,
    `brent_dekker_solver.rs`, and `prelim_kepler/{prelim_elliptic,
    prelim_hyperbolic, prelim_parabolic}.rs`.
  - Verified against an independent two-body ODE integration
    (`test_propag`..`test_propag4` in `src/kepler/propagation.rs`).
- **`velocity_correction` — missing `sqrt(mu)` factor on `sig0`**
  - `src/kepler/velocity.rs` computed the radial-velocity proxy as
    `x2.dot(v2)` instead of `x2.dot(v2) / sqrt(mu)`, a second, independent
    bug in the same Lagrange f/g machinery used by the Gauss IOD velocity
    refinement loop (`initial_orbit_determination/gauss.rs`).
  - Fixing it together with the `alpha` bug above restores agreement with
    the Orbfit reference orbit in `gauss_test::test_solve_orbit` to the
    original 1e-13 tolerance (it silently diverged while only one of the
    two bugs was fixed).

### Added

- **`propagate_universal` test coverage** (`src/kepler/propagation.rs`)
  - `tests_propagation_edge_cases` — deterministic tests for quasi-circular,
    high-eccentricity, near-parabolic (both signs of `alpha`), and
    hyperbolic orbits; `dt ≈ 0`; negative `dt`; survey-cadence gaps from 35
    to 400+ days (including multi-revolution gaps); continuity across a gap
    boundary; degenerate (zero-norm) position input; warm-start vs
    cold-start consistency; and agreement between `NewtonRaphson`,
    `BrentDecker`, and `Auto` solver kinds.
  - `tests_propagation_invariants` — property-based tests (`proptest`)
    checking invariants that must hold for any valid two-body state: no
    panics, the Lagrange identity `f·ġ − g·ḟ = 1`, conservation of specific
    orbital energy and angular momentum, forward/backward round-trip
    consistency, solver-kind agreement, chained small steps vs. a single
    large step, and continuity under a tiny input perturbation. Orbits and
    time-of-flight gaps are sampled from distributions matching real
    fink-fat/Kalman-filter usage (periapsis/eccentricity up to `e = 5`,
    survey-cadence gaps biased toward the >30-day inter-lunation and
    inter-season regimes actually seen with ZTF/LSST).

## [4.0.0] - 2026-07-16

### Added

- **`kepler` module split and rewritten** (`src/kepler.rs` → `src/kepler/`)
  - The former monolithic `kepler.rs` (~1800 lines) is now a module tree:
    `angles.rs`, `orbit_type.rs`, `params.rs`, `stumpff.rs`,
    `newton_solver.rs`, `brent_dekker_solver.rs`,
    `universal_kepler_solution.rs`, `velocity.rs`, `propagation.rs`, and
    `prelim_kepler/` (elliptic, hyperbolic, parabolic).
  - New **BrentDekker solver** (`solve_kepuni_brent_dekker`) for the
    universal anomaly ψ, selectable alongside the existing Newton–Raphson
    solver via `SolverKind` / `SolverType`.
  - New **parabolic case** handling in the universal anomaly preliminary
    guess (`prelim_kepler/prelim_parabolic.rs`), previously unsupported.
  - `SolverParams`, `SolverKind`, `SolverType` — new configuration types
    letting callers pick and tune the Kepler solver used during
    propagation/fitting.
  - `UniversalKeplerSolution` — new result container for the universal
    Kepler solve.
  - New `OutfitError` variants: `DegenerateState`,
    `NewtonRaphsonKeplerConvergence`, `BrentDekkerKeplerConvergence`, giving
    solver-specific convergence failures instead of a generic error.

### Changed

- **`serde` feature** — `SolverParams`, `SolverKind`, and `SolverType` now
  derive `Serialize`/`Deserialize` under the new optional `serde` feature,
  so Kepler solver configuration can be read from / written to top-level
  config files.
  - New optional dependency: [`serde`](https://crates.io/crates/serde)
    (`derive` feature).
- **Brent–Dekker solver** now accepts an initial ψ guess when one is
  available, instead of always bracketing from scratch.
- `rayon` is now declared as `dep:rayon` (matching `serde`'s
  `dep:serde`) rather than an implicit optional dependency.

### Breaking

- Major version bump to reflect the `kepler` module restructuring
  (`src/kepler.rs` is gone; all prior `kepler::*` imports still resolve via
  the new `src/kepler/mod.rs` re-exports, but internal submodule paths have
  changed).

## [3.1.0] - 2026-06-15

### Added

- **`differential_correction` public entry point**
  - New standalone `differential_correction` function exposed in the
    `differential_orbit_correction` module.
  - Enables calling the differential correction pipeline directly on a single
    trajectory, outside of the `FitLSQ` trait machinery.
  - Useful for custom pipeline construction or one-off corrections without
    going through `ObsDataset`.

- **`OrbitalElements::ref_epoch`**
  - New method returning the reference epoch of any `OrbitalElements` variant
    (Keplerian, equinoctial, etc.).
  - Provides a uniform accessor regardless of the underlying orbit
    representation.

- **`FitOrbitResult::epoch`**
  - New method returning the epoch associated with a `FitOrbitResult`.
  - Mirrors `OrbitalElements::ref_epoch` at the result level for convenience.

### Changed

- **`TrajectoryFit` trait is now public**
  - Previously a crate-internal trait; now exposed as `pub` for downstream
    users who need finer control over the fitting pipeline.
  - **Note**: the exposed methods are low-level and intended for advanced use
    only. Prefer `FitLSQ` for standard orbit fitting workflows.


## [3.0.0] - 2026-06-03

### Added

- **N-body propagator (`propagator/`)**
  - New `propagator` module exposing two propagation strategies selectable at runtime via `PropagatorKind`:
    - `PropagatorKind::TwoBody` — analytic Keplerian propagation (default, no extra dependencies).
    - `PropagatorKind::NBody(NBodyConfig)` — numerical N-body integration via an 8th-order Dormand–Prince (DOP853) integrator from the `differential-equations` crate.
  - `NBodyConfig` — configuration struct selecting which planetary perturbers to include.
  - `NBodyResult` — propagation output holding the final heliocentric state and the 6×6 State Transition Matrix (STM) Φ.
  - The STM is propagated alongside the trajectory (variational equations `dΦ/dt = A(t)·Φ`), enabling N-body-quality partial derivatives for the differential corrector.
  - `propagator::planet_gm` — planetary GM constants table.
  - New dependency: [`differential-equations`](https://crates.io/crates/differential-equations) (DOP853 integrator).
  - Integration tests for N-body propagation added in `tests/test_diff_cor.rs`.

- **Ephemeris generation (`ephemeris/`)**
  - Entirely new `ephemeris` module providing a complete, type-safe pipeline to predict apparent sky positions and geometric quantities from an orbit:
    - **`EphemerisRequest<O>`** — typed builder pairing observers with generation modes; `O` is one of three zero-cost marker types:
      - `Position` — computes apparent (RA, Dec), geocentric/heliocentric distances → `ApparentPosition`.
      - `Geometry` — computes phase angle, solar elongation, radial velocity, angular rates → `BodyGeometry`.
      - `Combined` — both of the above from a single propagation → `(ApparentPosition, BodyGeometry)`.
    - **`EphemerisMode`** — epoch generation strategy:
      - `Single` — single epoch.
      - `Range { start, end, step }` — uniform grid.
      - `At(Vec<Epoch>)` — arbitrary epoch list.
    - **`EphemerisResult<O>`** — result container; per-entry errors are recorded without short-circuiting the full computation.
      - `successes()`, `errors()`, `by_observer()` accessors.
    - **`EphemerisConfig`** — global configuration (aberration corrections, etc.).
    - **`OrbitalElements::compute`** — entry point: `elements.compute(&request, &jpl, &ut1)`.
    - **`ApparentPosition`** — apparent RA/Dec, geocentric distance, heliocentric distance, light-travel-time corrected.
    - **`BodyGeometry`** — phase angle, solar elongation, radial velocity, apparent angular rates.
    - **Aberration correction** (`aberration.rs`) — stellar aberration applied to apparent positions.
    - **Batch ephemeris (`batch.rs`)** — `FullOrbitResultExt` extension trait on `FullOrbitResult`:
      - `compute_ephemerides(&request, &jpl, &ut1)` — sequential computation over all orbits in a `FullOrbitResult`.
      - `compute_ephemerides_parallel(&request, &jpl, &ut1)` *(feature `parallel`)* — Rayon-parallel variant.
    - Integration tests in `tests/test_ephemeris.rs` (≈780 lines).

- **Orbital uncertainty & Jacobian propagation (`orbit_type/uncertainty.rs`)**
  - New `uncertainty` submodule providing first-order linear covariance propagation across all three orbital element representations:
    - `KeplerianUncertainty` — 1σ standard deviations for (a, e, i, Ω, ω, M).
    - `EquinoctialUncertainty` — 1σ standard deviations for (a, h, k, p, q, λ).
    - `CometaryUncertainty` — 1σ standard deviations for (q, e, i, Ω, ω, ν).
    - Each built via `from_covariance(&OrbitalCovariance)` extracting diagonal variances.
    - `OrbitalCovariance` — full 6×6 symmetric covariance matrix Σ with `propagate(jacobian)` → `Σ' = J Σ Jᵀ`.
  - **Jacobian methods** on orbital element structs:
    - `EquinoctialElements::jacobian_to_keplerian()` — ∂(a,e,i,Ω,ω,M)/∂(a,h,k,p,q,λ).
    - `KeplerianElements::jacobian_to_equinoctial()` and `jacobian_to_cometary()`.
    - `CometaryElements::jacobian_to_keplerian()` and `jacobian_to_equinoctial()`.
  - `OrbitalElements` (the top-level enum) gains:
    - `uncertainty() -> Option<OrbitalCovariance>` — returns the attached covariance if present.
    - Covariance is propagated and converted automatically when calling `into_keplerian()`, `into_equinoctial()`, `into_cometary()`, preserving uncertainty across representations.
  - `DifferentialCorrectionOutput` updated to expose `OrbitalCovariance` through `OrbitalElements` directly.

- **Differential orbit correction (`differential_orbit_correction`)**
  - New module implementing a full **least-squares differential correction** pipeline over optical astrometry, using 2-body (Keplerian) propagation:
    - `DifferentialCorrectionConfig` — configuration struct (max iterations, convergence threshold, RMS divergence ratio, outlier rejection settings).
    - `FitLSQ` trait — top-level entry point on `ObsDataset`: `obs_dataset.fit_lsq(...)` returns a `FullOrbitResult` with converged equinoctial elements and per-observation fit data.
    - `DifferentialCorrectionOutput` — result type holding final `EquinoctialElements`, per-observation `ObsFitData` (residuals, σ, χ, selection status), `OrbitalUncertainty` (normal matrix, covariance, inversion flag), normalised RMS, iteration count, and measurement count.
    - `ObsFitData` — per-observation fit metadata: `residual_ra`, `residual_dec`, `sigma_ra`, `sigma_dec`, `bias_ra`, `bias_dec`, `chi`, and `ObsSelection` (`Active` / `Rejected`).
    - `OrbitalUncertainty` — 6×6 normal matrix and covariance matrix with `inversion_succeeded` flag.
    - **Outlier rejection** (`outlier_rejection.rs`) — iterative χ²-based rejection with configurable threshold.
    - **Single Newton iteration** (`single_iteration.rs`) — computes the design matrix G (∂ρ/∂x), normal matrix GᵀWG, right-hand side, applies the Δx correction to equinoctial elements, and returns `correction_norm`.
    - **Least-squares accumulator** (`least_square.rs`) — builds the weighted normal equations from active observations, computes normalised RMS.
    - `OutfitError::DifferentialCorrectionDiverged` and `OutfitError::DifferentialCorrectionNotConverged` — new error variants.

- **Observation ephemeris (`observation_ephemeris.rs`)**
  - Massively extended module (≈ 800 lines added):
    - Partial derivatives of apparent (topocentric) RA/DEC with respect to equinoctial orbital elements: `∂ρ/∂a`, `∂ρ/∂h`, `∂ρ/∂k`, `∂ρ/∂p`, `∂ρ/∂q` — used to construct the design matrix G for differential correction.
    - Topocentric coordinate computation from equinoctial elements and observer position.
    - Observation weighting via error-model sigmas (`ObsWeight`).

- **Cache extensions (`cache/`)**
  - `observer_centric_cache.rs` extended with observer heliocentric velocity (`helio_velocity`) computation and storage — required by the design matrix G for velocity-dependent partial derivatives.
  - `ObserverCentricCache` now derives `Debug`.

- **`EquinoctialElements` additions**
  - `is_bizarre() -> bool` — returns `true` if the equinoctial orbit has unphysical or degenerate parameters (used during differential correction to reject diverging solutions).
  - Additional conversion utilities and `Display` improvements.

- **`JPLEphem` ergonomics**
  - `impl TryFrom<&str> for JPLEphem` and `impl TryFrom<String> for JPLEphem` — construct directly from a source string (e.g. `"horizon:DE440".try_into()`), without explicitly building an `EphemFileSource` first.
  - `JPLEphem::new` signature generalised to `impl Into<EphemFileSource>` — accepts both `EphemFileSource` values and references.
  - `impl From<&EphemFileSource> for EphemFileSource` added in `download_jpl_file.rs` (clone-based).

- **Constants (`constants.rs`)**
  - Extended with additional physical and astronomical constants required by the differential correction (e.g. Gaussian gravitational constant, AU/day conversions).

- **Integration tests**
  - `tests/test_diff_cor.rs` — non-regression test for the full differential correction pipeline on three real asteroids (2015 AB / K09R05F, 33803, 8467) with seed-42 oracles; uses `approx_equal` from `tests/common/mod.rs`.
  - `tests/test_gauss_iod.rs` and `tests/test_iod_from_polars.rs` — migrated to use `approx_equal` for orbit comparisons.

- **CI improvements**
  - `fmt`, `clippy`, and `semver-pr` jobs now run **only on pull requests** (skipped on push to `main` after merge).
  - `coverage` job runs on both PR and `main`, but only after `fmt`, `clippy`, and `semver-pr` are green (or skipped).
  - Removed the redundant `test` job — coverage via `cargo llvm-cov` already executes all tests.
  - Added a **feature matrix** (`default`, `parallel`) to `clippy` and `coverage` jobs.
  - Codecov uploads use per-feature `flags` (`default` / `parallel`) for granular coverage reporting.

### Changed

- **Orbit display simplified** — `Display` implementations for `KeplerianElements`, `EquinoctialElements`, and `CometaryElements` simplified; `GaussResult` display updated accordingly. Verbose per-field multi-line output replaced by a more compact single-block format.

- **`PropagatorKind` integrated into differential corrector** — `DifferentialCorrectionConfig` now accepts a `PropagatorKind` field; the corrector dispatches to two-body or N-body propagation transparently. Existing configs default to `PropagatorKind::TwoBody` (no migration required).

- **`observation_ephemeris.rs` moved** — `src/observation_ephemeris.rs` relocated to `src/ephemeris/observation_ephemeris.rs`; public re-exports updated in `src/lib.rs`.

- **CI: test thread count capped** — `RUST_TEST_THREADS` set in CI workflow to limit concurrency during `cargo test`, reducing peak memory usage from heavy doctests.

- **`photom` crate extracted** — observation parsing, error models, observer management, and trajectory ingestion have been moved into a dedicated `photom` crate (published on crates.io). Outfit now depends on `photom` as an external dependency.
  - Removed from Outfit: `src/observations/`, `src/observers/`, `src/trajectories/`, `src/error_models/`, `src/outfit.rs`, `src/env_state.rs` and associated modules.
  - Public re-exports updated in `src/lib.rs`; user-facing types (`ObsDataset`, `ObsErrorModel`, `TrajId`, `Observer`, `Observatories`, …) are now accessed via `photom`.

- **`FullOrbitResult` type alias updated**
  - Now defined as `HashMap<TrajId, Option<DifferentialCorrectionOutput>>` — keyed by `TrajId` (from `photom`), value is `None` if differential correction did not converge.

- **`GaussResult` / IOD result types**
  - `GaussResult` now stores `EquinoctialElements` (instead of `KeplerianElements`) as the primary output of Gauss IOD, feeding directly into the differential corrector.

- **`JPLEphem::new` signature** — changed from `fn new(file_source: &EphemFileSource)` to `fn new(source: impl Into<EphemFileSource>)` (see **Added** above).

- **Examples** — `run_full_iod.rs` and `run_full_iod_parallel.rs` updated to use the new `photom`-based API and the `fit_lsq` / `fit_full_iod` entry points.

- **`tests/common/mod.rs`** — fixed a bug in `approx_equal` for `OrbitalElements::Equinoctial`: `ee1.semi_major_axis` was incorrectly compared against `0.0` instead of `ee2.semi_major_axis`.

### Removed

- Removed benchmarks (`benches/`) — `gauss_prelim_orbit`, `load_parquet`, `outfit_gauss_iod`, `solve_kepler_equation` benches removed as part of the photom extraction refactor.
- Removed old examples (`examples/gauss_iod_once.rs`, `examples/parquet_to_orbit.rs`).
- Removed old test files (`tests/outfit_struct_test.rs`, `tests/reader_80col_test.rs`, `tests/test_read_ades.rs`, `tests/trajectories_from_parquet.rs`, `tests/trajectories_from_vec.rs`, `tests/vec_to_iod.rs`).

### Breaking Changes

- **`photom` dependency required** — types previously in `outfit::observations`, `outfit::observers`, `outfit::trajectories`, `outfit::error_models` are now in the `photom` crate. Update imports:
  ```rust
  // Before
  use outfit::observations::ObsDataset;
  use outfit::error_models::ObsErrorModel;

  // After
  use photom::observation_dataset::ObsDataset;
  use photom::observer::error_model::ObsErrorModel;
  ```

- **`FullOrbitResult` keyed by `TrajId` instead of `ObjectNumber`** — update map lookups accordingly.

- **`JPLEphem::new` accepts `impl Into<EphemFileSource>`** — existing call sites passing `&EphemFileSource` continue to work (via `From<&EphemFileSource>`); call sites passing an owned `EphemFileSource` no longer need the `&`.

### Notes

- The differential corrector now supports both **2-body (Keplerian)** and **N-body** propagation via `PropagatorKind`. The default remains `TwoBody`. For well-observed main-belt asteroids, select `NBody` with the appropriate planetary perturbers for higher accuracy. The default `rms_divergence_ratio` (1.5) in `DifferentialCorrectionConfig` can be raised if needed during convergence.
- The `dev` profile no longer emits debug symbols (`debug = false`) to reduce RAM consumption when running doctests.
- Ephemeris computation is designed to be fault-tolerant: individual `(epoch, observer)` errors are captured per-entry inside `EphemerisResult` without aborting the whole batch.

---

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
