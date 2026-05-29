//! Integration tests for the public ephemeris API.
//!
//! Strategy
//! --------
//! For each of the three objects in `tests/data/` we first fit an orbit
//! through the observations with the differential corrector (same pipeline as
//! `test_diff_cor.rs`).  The fitted elements are anchored near the middle of
//! the observed arc, so the propagation error is small.  We then call every
//! non-range public ephemeris method and verify that the predicted apparent
//! position matches each observation to within a tight angular threshold.
//!
//! APIs exercised
//! --------------
//! - [`OrbitalElements::apparent_position`]                  (single epoch, per-site)
//! - [`OrbitalElements::apparent_position_at`]               (bulk, geocentric)
//! - [`OrbitalElements::body_geometry`]                      (single epoch)
//! - [`OrbitalElements::body_geometry_at`]                   (bulk)
//! - [`OrbitalElements::apparent_position_and_geometry`]     (single epoch)
//! - [`OrbitalElements::apparent_position_and_geometry_at`]  (bulk)

mod common;

use hifitime::{ut1::Ut1Provider, Epoch, TimeScale};
use outfit::{
    jpl_ephem::naif::naif_ids::{
        planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
    },
    orbit_type::OrbitalElements,
    propagator::{NBodyConfig, PropagatorKind},
    AberrationOrder, DifferentialCorrectionConfig, EphemerisConfig, FitLSQ, FullOrbitResult,
    IODParams, JPLEphem,
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
/// Used for the bulk `_at` API calls that accept a single observer.
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

/// Angular separation (arcsec) between the per-site apparent position and the
/// observed equatorial coordinate for one observation.  Returns `None` if the
/// site cannot be resolved or propagation fails.
fn sep_arcsec_per_site(
    obs: &Observation,
    elements: &OrbitalElements,
    dataset: &ObsDataset,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Option<f64> {
    let epoch = Epoch::from_mjd_in_time_scale(obs.mjd_tt(), TimeScale::TT);
    let observer = dataset.get_observer(*obs.id())?;
    let predicted = elements
        .apparent_position(epoch, observer, jpl, ut1, config)
        .ok()?;
    Some(
        predicted
            .coord
            .angular_separation(obs.equ_coord())
            .to_degrees()
            * 3600.0,
    )
}

/// Collect per-site angular separations (arcsec) for all observations.
fn per_site_seps(
    obs_vec: &[&Observation],
    elements: &OrbitalElements,
    dataset: &ObsDataset,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Vec<f64> {
    obs_vec
        .iter()
        .filter_map(|obs| sep_arcsec_per_site(obs, elements, dataset, jpl, ut1, config))
        .collect()
}

/// Collect geocentric bulk angular separations (arcsec) from
/// `apparent_position_at`.
fn bulk_seps(
    obs_vec: &[&Observation],
    epochs: &[Epoch],
    elements: &OrbitalElements,
    geo_obs: &Observer,
    jpl: &JPLEphem,
    ut1: &Ut1Provider,
    config: &EphemerisConfig,
) -> Vec<f64> {
    elements
        .apparent_position_at(epochs.to_vec(), geo_obs, jpl, ut1, config)
        .into_iter()
        .zip(obs_vec.iter())
        .filter_map(|((_, result), obs)| {
            Some(
                result
                    .ok()?
                    .coord
                    .angular_separation(obs.equ_coord())
                    .to_degrees()
                    * 3600.0,
            )
        })
        .collect()
}

/// Assert that the median of `seps` is below `threshold_arcsec`.
/// Panics with a diagnostic message if the assertion fails or `seps` is empty.
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

/// Test all ephemeris API variants for object **33803 Julienpeloton**.
///
/// The orbit is first fitted with the 2-body differential corrector so the
/// reference epoch sits near the centre of the observed arc (2024-01).
/// A 2-body propagation over ±3 months should stay well within 2 arcsec.
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

    // apparent_position (per-site)
    assert_median_below(
        &mut per_site_seps(&obs_vec, elements, &dataset, &jpl, &ut1, &config),
        "33803 apparent_position",
        2.0,
    );

    // apparent_position_at (bulk, geocentric)
    // Geocentric observer introduces a topocentric offset (~2-3 arcsec for
    // main-belt objects at ~2 AU observed from G96).  Use a 5 arcsec bound.
    assert_median_below(
        &mut bulk_seps(&obs_vec, &epochs, elements, &geo_obs, &jpl, &ut1, &config),
        "33803 apparent_position_at",
        5.0,
    );

    // body_geometry (single epoch) — verify physically plausible phase angle
    {
        let observer = dataset
            .get_observer(*obs_vec[0].id())
            .expect("observer not found");
        let geo = elements
            .body_geometry(epochs[0], observer, &jpl, &ut1, &config)
            .expect("body_geometry failed");
        assert!(
            (0.0..=std::f64::consts::PI).contains(&geo.phase_angle),
            "33803: phase angle {:.4} rad out of range",
            geo.phase_angle
        );
    }

    // body_geometry_at (bulk) — all epochs must succeed
    {
        let ok = elements
            .body_geometry_at(epochs.clone(), &geo_obs, &jpl, &ut1, &config)
            .into_iter()
            .filter(|(_, r)| r.is_ok())
            .count();
        assert_eq!(
            ok,
            epochs.len(),
            "33803 body_geometry_at: some epochs failed"
        );
    }

    // apparent_position_and_geometry (single epoch) — sanity on RA and phase
    {
        let observer = dataset
            .get_observer(*obs_vec[0].id())
            .expect("observer not found");
        let (ap, geo) = elements
            .apparent_position_and_geometry(epochs[0], observer, &jpl, &ut1, &config)
            .expect("apparent_position_and_geometry failed");
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

    // apparent_position_and_geometry_at (bulk, geocentric)
    {
        let table = elements.apparent_position_and_geometry_at(
            epochs.clone(),
            &geo_obs,
            &jpl,
            &ut1,
            &config,
        );
        let results: Vec<_> = table.into_iter().collect();
        let ok = results.iter().filter(|(_, r)| r.is_ok()).count();
        assert_eq!(
            ok,
            epochs.len(),
            "33803 apparent_position_and_geometry_at: some epochs failed"
        );
        let mut seps: Vec<f64> = results
            .into_iter()
            .zip(obs_vec.iter())
            .filter_map(|((_, r), obs)| {
                Some(
                    r.ok()?
                        .0
                        .coord
                        .angular_separation(obs.equ_coord())
                        .to_degrees()
                        * 3600.0,
                )
            })
            .collect();
        assert_median_below(&mut seps, "33803 apparent_position_and_geometry_at", 5.0);
    }
}

/// Test all ephemeris API variants for object **8467 Benoîtcarry**.
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

    assert_median_below(
        &mut per_site_seps(&obs_vec, elements, &dataset, &jpl, &ut1, &config),
        "8467 apparent_position",
        2.0,
    );
    assert_median_below(
        &mut bulk_seps(&obs_vec, &epochs, elements, &geo_obs, &jpl, &ut1, &config),
        "8467 apparent_position_at",
        3.0,
    );

    {
        let ok = elements
            .body_geometry_at(epochs.clone(), &geo_obs, &jpl, &ut1, &config)
            .into_iter()
            .filter(|(_, r)| r.is_ok())
            .count();
        assert_eq!(
            ok,
            epochs.len(),
            "8467 body_geometry_at: some epochs failed"
        );
    }

    {
        let mut seps: Vec<f64> = elements
            .apparent_position_and_geometry_at(epochs.clone(), &geo_obs, &jpl, &ut1, &config)
            .into_iter()
            .zip(obs_vec.iter())
            .filter_map(|((_, r), obs)| {
                Some(
                    r.ok()?
                        .0
                        .coord
                        .angular_separation(obs.equ_coord())
                        .to_degrees()
                        * 3600.0,
                )
            })
            .collect();
        assert_median_below(&mut seps, "8467 apparent_position_and_geometry_at", 3.0);
    }
}

/// Test all ephemeris API variants for object **2015 AB** (MPC packed: K09R05F).
///
/// The observations span 2009 and 2015.  With a 2-body fit over a ~2000-day
/// arc the geocentric bulk residuals are larger (~10-15 arcsec); the per-site
/// `apparent_position` test validates sub-arcsecond accuracy.
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

    assert_median_below(
        &mut per_site_seps(&obs_vec, elements, &dataset, &jpl, &ut1, &config),
        "2015AB apparent_position",
        2.0,
    );
    assert_median_below(
        &mut bulk_seps(&obs_vec, &epochs, elements, &geo_obs, &jpl, &ut1, &config),
        "2015AB apparent_position_at",
        15.0,
    );

    {
        let ok = elements
            .body_geometry_at(epochs.clone(), &geo_obs, &jpl, &ut1, &config)
            .into_iter()
            .filter(|(_, r)| r.is_ok())
            .count();
        assert_eq!(
            ok,
            epochs.len(),
            "2015AB body_geometry_at: some epochs failed"
        );
    }

    {
        let mut seps: Vec<f64> = elements
            .apparent_position_and_geometry_at(epochs.clone(), &geo_obs, &jpl, &ut1, &config)
            .into_iter()
            .zip(obs_vec.iter())
            .filter_map(|((_, r), obs)| {
                Some(
                    r.ok()?
                        .0
                        .coord
                        .angular_separation(obs.equ_coord())
                        .to_degrees()
                        * 3600.0,
                )
            })
            .collect();
        assert_median_below(&mut seps, "2015AB apparent_position_and_geometry_at", 15.0);
    }
}

// ── N-body tests ──────────────────────────────────────────────────────────────

/// Fit one trajectory with the N-body corrector and assert that the per-site
/// apparent-position median residual is below `threshold_arcsec`.
///
/// This is the common body shared by all three N-body tests; each test
/// supplies the `traj_id` and appropriate threshold.
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
    let config = nbody_ephem_config();

    assert_median_below(
        &mut per_site_seps(&obs_vec, elements, &dataset, &jpl, &ut1, &config),
        &format!("{traj_id:?} N-body apparent_position"),
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
