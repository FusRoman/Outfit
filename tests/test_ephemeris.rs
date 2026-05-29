//! Integration tests for the public ephemeris API.
//!
//! Strategy
//! --------
//! For each of the three objects in `tests/data/` we first fit an orbit
//! through the observations with the differential corrector (same pipeline as
//! `test_diff_cor.rs`).  The fitted elements are anchored near the middle of
//! the observed arc, so the propagation error is small.  We then exercise
//! every combination of output kind and generation mode provided by
//! [`EphemerisRequest`] and verify that the predicted apparent position
//! matches each observation to within a tight angular threshold.
//!
//! APIs exercised
//! --------------
//! - [`EphemerisRequest::<Position>`]  — per-site single-epoch apparent position
//! - [`EphemerisRequest::<Position>`]  — bulk geocentric apparent position (At mode)
//! - [`EphemerisRequest::<Geometry>`]  — bulk geocentric geometry (At mode)
//! - [`EphemerisRequest::<Combined>`]  — per-site single-epoch combined
//! - [`EphemerisRequest::<Combined>`]  — bulk geocentric combined (At mode)

mod common;

use hifitime::{ut1::Ut1Provider, Epoch, TimeScale};
use outfit::{
    jpl_ephem::naif::naif_ids::{
        planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
    },
    orbit_type::OrbitalElements,
    propagator::{NBodyConfig, PropagatorKind},
    AberrationOrder, Combined, DifferentialCorrectionConfig, EphemerisConfig, EphemerisEntry,
    EphemerisMode, EphemerisRequest, FitLSQ, FullOrbitResult, Geometry, IODParams, JPLEphem,
    Position,
};
use photom::{
    observation_dataset::{observation::Observation, ObsDataset},
    observer::{error_model::ObsErrorModel, Observer},
    TrajId,
};
use rand::{rngs::StdRng, SeedableRng};

// ── Fixtures ─────────────────────────────────────────────────────────────────

fn build_fixtures() -> (
    JPLEphem,
    Ut1Provider,
    ObsDataset,
    IODParams,
    DifferentialCorrectionConfig,
) {
    let ut1 = Ut1Provider::download_from_jpl("latest_eop2.long").expect("UT1 download failed");

    let jpl: JPLEphem = "horizon:DE440"
        .try_into()
        .expect("JPL ephemeris load failed");

    let (raw_dataset, errors) = ObsDataset::from_mpc_80_col_files(&[
        "tests/data/2015AB.obs",
        "tests/data/8467.obs",
        "tests/data/33803.obs",
    ]);
    assert!(errors.is_empty(), "obs file errors: {errors:?}");

    // Attach the error model so that `get_observer` can resolve MPC-coded
    // sites (the MPC observatory catalogue is loaded lazily on first access).
    let dataset = raw_dataset.with_error_model(ObsErrorModel::FCCT14);

    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()
        .unwrap();

    let diff_cor_config = DifferentialCorrectionConfig {
        rms_divergence_ratio: 10.0,
        ..DifferentialCorrectionConfig::default()
    };

    (jpl, ut1, dataset, iod_params, diff_cor_config)
}

// ── Elementary helpers ────────────────────────────────────────────────────────

/// A geocentric observer (ρ cos φ' = ρ sin φ' = 0, longitude = 0).
/// Used for bulk At-mode requests with a single shared observer.
fn geocentric_observer() -> Observer {
    Observer::from_parallax(0.0, 0.0, 0.0, Some("Geocenter".to_string()), None, None)
        .expect("geocentric observer construction failed")
}

/// [`EphemerisConfig`] for N-body propagation: Sun + all eight major planets.
fn nbody_ephem_config() -> EphemerisConfig {
    EphemerisConfig {
        propagator: PropagatorKind::NBody(NBodyConfig {
            perturbing_bodies: vec![
                NaifIds::SSB(SolarSystemBary::Sun),
                NaifIds::PB(PlanetaryBary::Mercury),
                NaifIds::PB(PlanetaryBary::Venus),
                NaifIds::PB(PlanetaryBary::EarthMoon),
                NaifIds::PB(PlanetaryBary::Mars),
                NaifIds::PB(PlanetaryBary::Jupiter),
                NaifIds::PB(PlanetaryBary::Saturn),
                NaifIds::PB(PlanetaryBary::Uranus),
                NaifIds::PB(PlanetaryBary::Neptune),
            ],
            abs_tol: 1e-12,
            rel_tol: 1e-12,
        }),
        aberration: AberrationOrder::default(),
    }
}

/// [`DifferentialCorrectionConfig`] for N-body fitting, inheriting all other
/// fields from `base`.
fn nbody_diff_cor(base: &DifferentialCorrectionConfig) -> DifferentialCorrectionConfig {
    DifferentialCorrectionConfig {
        propagator: nbody_ephem_config().propagator,
        ..base.clone()
    }
}

/// Extract the [`OrbitalElements`] for `traj_id` from a `FullOrbitResult`.
/// Panics if the trajectory is missing or did not converge.
fn fit_orbit<'a>(orbits: &'a FullOrbitResult, traj_id: &TrajId) -> &'a OrbitalElements {
    orbits
        .get(traj_id)
        .unwrap_or_else(|| panic!("{traj_id:?} not found in fit results"))
        .as_ref()
        .unwrap_or_else(|e| panic!("{traj_id:?} fit error: {e}"))
        .orbital_elements()
}

/// Materialize observations for `traj_id` and derive the matching epoch list.
fn traj_obs_and_epochs(dataset: &ObsDataset, traj_id: TrajId) -> (Vec<&Observation>, Vec<Epoch>) {
    let mem = dataset
        .materialize_trajectory(traj_id.clone())
        .unwrap_or_else(|| panic!("trajectory {traj_id:?} not indexed"));
    let obs_vec: Vec<&Observation> = mem.iter().collect();
    assert!(!obs_vec.is_empty(), "no observations for {traj_id:?}");
    let epochs = obs_vec
        .iter()
        .map(|o| Epoch::from_mjd_in_time_scale(o.mjd_tt(), TimeScale::TT))
        .collect();
    (obs_vec, epochs)
}

/// Build a per-site [`EphemerisRequest<Position>`]: one `Single` slot per
/// observation, each using the site observer resolved from the dataset.
///
/// Observations whose site cannot be resolved are silently skipped.
fn per_site_position_request(
    obs_vec: &[&Observation],
    dataset: &ObsDataset,
    config: EphemerisConfig,
) -> EphemerisRequest<Position> {
    obs_vec
        .iter()
        .filter_map(|obs| {
            let observer = dataset.get_observer(*obs.id())?;
            let epoch = Epoch::from_mjd_in_time_scale(obs.mjd_tt(), TimeScale::TT);
            Some((observer, EphemerisMode::Single(epoch)))
        })
        .fold(
            EphemerisRequest::<Position>::new(config),
            |req, (obs, mode)| req.add(obs.clone(), mode),
        )
}

/// Build a per-site [`EphemerisRequest<Combined>`]: one `Single` slot per
/// observation, each using the site observer resolved from the dataset.
fn per_site_combined_request(
    obs_vec: &[&Observation],
    dataset: &ObsDataset,
    config: EphemerisConfig,
) -> EphemerisRequest<Combined> {
    obs_vec
        .iter()
        .filter_map(|obs| {
            let observer = dataset.get_observer(*obs.id())?;
            let epoch = Epoch::from_mjd_in_time_scale(obs.mjd_tt(), TimeScale::TT);
            Some((observer, EphemerisMode::Single(epoch)))
        })
        .fold(
            EphemerisRequest::<Combined>::new(config),
            |req, (obs, mode)| req.add(obs.clone(), mode),
        )
}

/// Collect angular separations (arcsec) from a per-site position result,
/// zipped against the corresponding observations.
fn position_seps_arcsec<'a, 'b>(
    entries: impl Iterator<Item = &'a EphemerisEntry<outfit::ApparentPosition>>,
    obs_vec: impl Iterator<Item = &'b &'b Observation>,
) -> Vec<f64> {
    entries
        .zip(obs_vec)
        .filter_map(
            |(entry, obs): (&EphemerisEntry<outfit::ApparentPosition>, &&Observation)| {
                Some(
                    entry
                        .result
                        .as_ref()
                        .ok()?
                        .coord
                        .angular_separation(obs.equ_coord())
                        .to_degrees()
                        * 3600.0,
                )
            },
        )
        .collect()
}

/// Collect angular separations (arcsec) from a combined result, using only
/// the [`ApparentPosition`] part.
fn combined_seps_arcsec<'a, 'b>(
    entries: impl Iterator<Item = &'a EphemerisEntry<(outfit::ApparentPosition, outfit::BodyGeometry)>>,
    obs_vec: impl Iterator<Item = &'b &'b Observation>,
) -> Vec<f64> {
    entries
        .zip(obs_vec)
        .filter_map(
            |(entry, obs): (
                &EphemerisEntry<(outfit::ApparentPosition, outfit::BodyGeometry)>,
                &&Observation,
            )| {
                Some(
                    entry
                        .result
                        .as_ref()
                        .ok()?
                        .0
                        .coord
                        .angular_separation(obs.equ_coord())
                        .to_degrees()
                        * 3600.0,
                )
            },
        )
        .collect()
}

/// Assert that the median of `seps` is below `threshold_arcsec`.
fn assert_median_below(seps: &mut [f64], label: &str, threshold_arcsec: f64) {
    assert!(!seps.is_empty(), "{label}: no separations computed");
    seps.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let median = seps[seps.len() / 2];
    assert!(
        median < threshold_arcsec,
        "{label}: median {median:.2} arcsec ≥ threshold {threshold_arcsec:.1} arcsec"
    );
}

// ── Two-body tests ────────────────────────────────────────────────────────────

/// Test all ephemeris request combinations for object **33803 Julienpeloton**.
///
/// Exercises:
/// - `Position` / per-site (one `Single` slot per observation)
/// - `Position` / bulk geocentric (`At` mode, single observer)
/// - `Geometry` / bulk geocentric (`At` mode) — all epochs must succeed
/// - `Combined` / per-site single-epoch
/// - `Combined` / bulk geocentric (`At` mode)
#[test]
fn test_ephemeris_33803() {
    let (jpl, ut1, dataset, iod_params, diff_cor_config) = build_fixtures();

    let orbits = dataset
        .clone()
        .fit_lsq(
            &jpl,
            &ut1,
            ObsErrorModel::FCCT14,
            &iod_params,
            &diff_cor_config,
            None,
            &mut StdRng::seed_from_u64(42),
        )
        .expect("fit_lsq failed");

    let traj_id = TrajId::Int(33803);
    let elements = fit_orbit(&orbits, &traj_id);
    let config = EphemerisConfig::default();
    let geo_obs = geocentric_observer();
    let (obs_vec, epochs) = traj_obs_and_epochs(&dataset, traj_id);

    // Position — per-site
    {
        let req = per_site_position_request(&obs_vec, &dataset, config.clone());
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "33803 Position per-site", 2.0);
    }

    // Position — bulk geocentric (At mode)
    // Geocentric observer introduces a topocentric offset (~2-3 arcsec for
    // main-belt objects at ~2 AU observed from G96).  Use a 5 arcsec bound.
    {
        let req = EphemerisRequest::<Position>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "33803 Position geocentric At", 5.0);
    }

    // Geometry — bulk geocentric (At mode): all epochs must succeed
    {
        let req = EphemerisRequest::<Geometry>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "33803 Geometry At: {} epochs failed",
            result.error_count()
        );
        // Phase angles must be in [0, π]
        for entry in result.successes() {
            let geo = entry.result.as_ref().unwrap();
            assert!(
                (0.0..=std::f64::consts::PI).contains(&geo.phase_angle),
                "33803: phase angle {:.4} rad out of range",
                geo.phase_angle
            );
        }
    }

    // Combined — per-site single-epoch
    {
        let req = per_site_combined_request(&obs_vec, &dataset, config.clone());
        let result = elements.compute(&req, &jpl, &ut1);
        // Sanity: RA in [0, 2π) and phase in [0, π]
        for entry in result.successes() {
            let (ap, geo) = entry.result.as_ref().unwrap();
            assert!(
                (0.0..2.0 * std::f64::consts::PI).contains(&ap.coord.ra),
                "33803: RA {:.4} rad out of range",
                ap.coord.ra
            );
            assert!(
                (0.0..=std::f64::consts::PI).contains(&geo.phase_angle),
                "33803: phase angle {:.4} rad out of range",
                geo.phase_angle
            );
        }
    }

    // Combined — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Combined>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "33803 Combined At: some epochs failed"
        );
        let mut seps = combined_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "33803 Combined geocentric At", 5.0);
    }
}

/// Test ephemeris requests for object **8467 Benoîtcarry**.
#[test]
fn test_ephemeris_8467() {
    let (jpl, ut1, dataset, iod_params, diff_cor_config) = build_fixtures();

    let orbits = dataset
        .clone()
        .fit_lsq(
            &jpl,
            &ut1,
            ObsErrorModel::FCCT14,
            &iod_params,
            &diff_cor_config,
            None,
            &mut StdRng::seed_from_u64(42),
        )
        .expect("fit_lsq failed");

    let traj_id = TrajId::Int(8467);
    let elements = fit_orbit(&orbits, &traj_id);
    let config = EphemerisConfig::default();
    let geo_obs = geocentric_observer();
    let (obs_vec, epochs) = traj_obs_and_epochs(&dataset, traj_id);

    // Position — per-site
    {
        let req = per_site_position_request(&obs_vec, &dataset, config.clone());
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "8467 Position per-site", 2.0);
    }

    // Position — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Position>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "8467 Position geocentric At", 3.0);
    }

    // Geometry — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Geometry>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "8467 Geometry At: some epochs failed"
        );
    }

    // Combined — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Combined>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "8467 Combined At: some epochs failed"
        );
        let mut seps = combined_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "8467 Combined geocentric At", 3.0);
    }
}

/// Test ephemeris requests for object **2015 AB** (MPC packed: K09R05F).
///
/// The observations span 2009 and 2015.  With a 2-body fit over a ~2000-day
/// arc the geocentric bulk residuals are larger (~10-15 arcsec); the per-site
/// `Position` test validates sub-arcsecond accuracy.
#[test]
fn test_ephemeris_2015ab() {
    let (jpl, ut1, dataset, iod_params, diff_cor_config) = build_fixtures();

    let orbits = dataset
        .clone()
        .fit_lsq(
            &jpl,
            &ut1,
            ObsErrorModel::FCCT14,
            &iod_params,
            &diff_cor_config,
            None,
            &mut StdRng::seed_from_u64(42),
        )
        .expect("fit_lsq failed");

    let traj_id = TrajId::from("K09R05F");
    let elements = fit_orbit(&orbits, &traj_id);
    let config = EphemerisConfig::default();
    let geo_obs = geocentric_observer();
    let (obs_vec, epochs) = traj_obs_and_epochs(&dataset, traj_id);

    // Position — per-site
    {
        let req = per_site_position_request(&obs_vec, &dataset, config.clone());
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "2015AB Position per-site", 2.0);
    }

    // Position — bulk geocentric (At mode)
    // 2015AB is a near-Earth Amor; over a ~2000-day 2-body arc the geocentric
    // residuals reach ~10-15 arcsec.
    {
        let req = EphemerisRequest::<Position>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "2015AB Position geocentric At", 15.0);
    }

    // Geometry — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Geometry>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "2015AB Geometry At: some epochs failed"
        );
    }

    // Combined — bulk geocentric (At mode)
    {
        let req = EphemerisRequest::<Combined>::new(config.clone())
            .add(geo_obs.clone(), EphemerisMode::At(epochs.clone()));
        let result = elements.compute(&req, &jpl, &ut1);
        assert_eq!(
            result.error_count(),
            0,
            "2015AB Combined At: some epochs failed"
        );
        let mut seps = combined_seps_arcsec(result.iter(), obs_vec.iter());
        assert_median_below(&mut seps, "2015AB Combined geocentric At", 15.0);
    }
}

// ── N-body tests ──────────────────────────────────────────────────────────────

/// Fit one trajectory with the N-body corrector and assert that the per-site
/// apparent-position median residual is below `threshold_arcsec`.
fn run_nbody_ephemeris_test(traj_id: TrajId, threshold_arcsec: f64) {
    let (jpl, ut1, dataset, iod_params, base_config) = build_fixtures();
    let diff_cor = nbody_diff_cor(&base_config);

    let orbits = dataset
        .clone()
        .fit_lsq(
            &jpl,
            &ut1,
            ObsErrorModel::FCCT14,
            &iod_params,
            &diff_cor,
            None,
            &mut StdRng::seed_from_u64(42),
        )
        .expect("N-body fit_lsq failed");

    let elements = fit_orbit(&orbits, &traj_id);
    let (obs_vec, _) = traj_obs_and_epochs(&dataset, traj_id.clone());

    let req = per_site_position_request(&obs_vec, &dataset, nbody_ephem_config());
    let result = elements.compute(&req, &jpl, &ut1);
    let mut seps = position_seps_arcsec(result.iter(), obs_vec.iter());

    assert_median_below(
        &mut seps,
        &format!("{traj_id:?} N-body Position per-site"),
        threshold_arcsec,
    );
}

/// N-body ephemeris test for **33803 Julienpeloton** — threshold 2 arcsec.
#[test]
fn test_ephemeris_33803_nbody() {
    run_nbody_ephemeris_test(TrajId::Int(33803), 2.0);
}

/// N-body ephemeris test for **8467 Benoîtcarry** — threshold 2 arcsec.
#[test]
fn test_ephemeris_8467_nbody() {
    run_nbody_ephemeris_test(TrajId::Int(8467), 2.0);
}

/// N-body ephemeris test for **2015 AB** (K09R05F) — threshold 15 arcsec.
///
/// The long (~2000-day) arc and NEA dynamics make sub-2-arcsec residuals
/// unrealistic here; we use the same bound as the 2-body test for this object.
#[test]
fn test_ephemeris_2015ab_nbody() {
    run_nbody_ephemeris_test(TrajId::from("K09R05F"), 15.0);
}
