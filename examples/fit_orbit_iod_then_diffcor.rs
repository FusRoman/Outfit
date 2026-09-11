//! Minimal walkthrough of Outfit's two-stage orbit-fitting pipeline:
//! **Initial Orbit Determination (IOD)**, then **differential correction**.
//!
//! IOD (the classical Gauss method) turns three well-spaced observations
//! into a rough preliminary orbit — fast, but only as accurate as those
//! three points. Differential correction then refines that seed with a
//! weighted least-squares fit using *every* observation, which is what you
//! want for a science-grade orbit. This example runs both stages explicitly
//! on a single object, then shows how to opt into N-body propagation for the
//! refinement step.
//!
//! Run it with:
//! ```text
//! cargo run --example fit_orbit_iod_then_diffcor
//! ```

use hifitime::ut1::Ut1Provider;
use outfit::jpl_ephem::naif::naif_ids::{
    planet_bary::PlanetaryBary, solar_system_bary::SolarSystemBary, NaifIds,
};
use outfit::propagator::{NBodyConfig, PropagatorKind};
use outfit::{DifferentialCorrectionConfig, FitIOD, FitLSQ, IODParams, JPLEphem, OutfitError};
use photom::observation_dataset::ObsDataset;
use photom::observer::error_model::ObsErrorModel;
use photom::TrajId;
use rand::{rngs::StdRng, SeedableRng};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Load observations for one object from an MPC 80-column file.
    //
    // A real orbit-fitting workflow usually starts here: astrometric
    // measurements (RA, Dec, epoch, observing site) for a single object,
    // spanning one or more nights.
    let (obs_dataset, errors) = ObsDataset::from_mpc_80_col_files(&["tests/data/33803.obs"]);
    assert!(
        errors.is_empty(),
        "failed to parse observations: {errors:?}"
    );
    let traj_id = TrajId::Int(33803);

    // Step 2 — Load the DE440 planetary ephemeris and Earth-orientation data.
    //
    // Both are needed to convert an observer's line of sight into a
    // heliocentric geometry: the ephemeris gives the Earth's (and any
    // perturbing planet's) position, `Ut1Provider` gives the Earth's
    // rotation angle at each observation epoch.
    let jpl_ephem: JPLEphem = "naif:DE440".try_into()?;
    let ut1_provider = Ut1Provider::download_from_jpl("latest_eop2.long")
        .expect("failed to download UT1 Earth-orientation data");

    // Step 3 — Initial Orbit Determination (Gauss method).
    //
    // `fit_full_iod` tries triplets of observations, solves Gauss's
    // classical problem for each, and keeps the candidate with the lowest
    // RMS. Each triplet is tried with several Monte-Carlo noise realizations
    // to average out astrometric noise; `max_obs_for_triplets` caps how many
    // of this trajectory's observations are considered when forming triplets.
    let iod_params = IODParams::builder()
        .n_noise_realizations(10)
        .noise_scale(1.1)
        .max_obs_for_triplets(130)
        .max_triplets(30)
        .build()?;
    let iod_result = obs_dataset.clone().fit_full_iod(
        &jpl_ephem,
        &ut1_provider,
        &iod_params,
        ObsErrorModel::FCCT14,
        &mut StdRng::seed_from_u64(42),
    )?;

    let seed_orbit = iod_result
        .get(&traj_id)
        .expect("trajectory 33803 was in the dataset")
        .as_ref()
        .expect("IOD should find a preliminary orbit for 33803");
    println!("=== After Gauss IOD (preliminary orbit) ===");
    println!("{:#?}", seed_orbit.orbital_elements());
    println!(
        "orbit quality (RMS of normalized residuals): {:.4}\n",
        seed_orbit.orbit_quality()
    );

    // Step 4 — Differential correction: refine the IOD seed with a
    // Newton-Raphson weighted least-squares fit using every observation.
    //
    // Passing `iod_result` as `initial_orbits` reuses the seed found above
    // instead of re-running IOD internally. `DifferentialCorrectionConfig`
    // defaults to two-body (Keplerian) propagation — fast, and accurate
    // enough for arcs where planetary perturbations are negligible. Note the
    // normalized RMS dropping toward 1.0: that means the residuals are now
    // consistent with the reported astrometric uncertainties.
    let diff_cor_config = DifferentialCorrectionConfig::default();
    let refined = obs_dataset.clone().fit_lsq(
        &jpl_ephem,
        &ut1_provider,
        ObsErrorModel::FCCT14,
        &iod_params,
        &diff_cor_config,
        Some(&iod_result),
        &mut StdRng::seed_from_u64(42),
    )?;
    let refined_orbit = refined
        .get(&traj_id)
        .unwrap()
        .as_ref()
        .expect("differential correction should converge for 33803");
    println!("=== After differential correction (two-body) ===");
    println!("{:#?}", refined_orbit.orbital_elements());
    println!(
        "normalized RMS (post-fit residuals / uncertainties): {:.4}\n",
        refined_orbit.orbit_quality()
    );

    // Step 5 — Optional: switch to N-body propagation for the refinement.
    //
    // Objects whose true orbit is measurably perturbed by a nearby massive
    // planet benefit from modeling that perturbation explicitly instead of
    // assuming an unperturbed two-body orbit. Swapping `propagator` to
    // `PropagatorKind::NBody` is the only change needed — everything else
    // (the seed orbit, the observations, the error model) stays the same.
    let nbody_config = NBodyConfig {
        perturbing_bodies: vec![
            NaifIds::SSB(SolarSystemBary::Sun),
            NaifIds::PB(PlanetaryBary::Jupiter),
        ],
        ..NBodyConfig::default()
    };
    let nbody_diff_cor_config = DifferentialCorrectionConfig {
        propagator: PropagatorKind::NBody(nbody_config),
        ..DifferentialCorrectionConfig::default()
    };
    let refined_nbody = obs_dataset.fit_lsq(
        &jpl_ephem,
        &ut1_provider,
        ObsErrorModel::FCCT14,
        &iod_params,
        &nbody_diff_cor_config,
        Some(&iod_result),
        &mut StdRng::seed_from_u64(42),
    )?;
    let refined_nbody_orbit = refined_nbody
        .get(&traj_id)
        .unwrap()
        .as_ref()
        .expect("N-body differential correction should converge for 33803");
    println!("=== After differential correction (N-body: Sun + Jupiter) ===");
    println!("{:#?}", refined_nbody_orbit.orbital_elements());
    println!(
        "normalized RMS (post-fit residuals / uncertainties): {:.4}",
        refined_nbody_orbit.orbit_quality()
    );

    Ok(())
}
