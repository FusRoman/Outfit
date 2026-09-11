//! Minimal walkthrough of the **ANISE** ephemeris backend.
//!
//! ANISE is a Rust reimplementation of the NAIF SPICE toolkit (`ephem-anise`
//! Cargo feature, not enabled by default). It reads the same NAIF SPK/DAF
//! kernels as the `ephem_naif.rs` example, so its results agree with the
//! built-in reader's to high precision — but it is validated against SPICE
//! itself and, unlike the built-in reader, can also load *supplementary*
//! kernels, such as the one covering 300 numbered main-belt asteroids.
//!
//! Run it with:
//! ```text
//! cargo run --example ephem_anise --no-default-features --features ephem-anise
//! ```

use hifitime::{Epoch, TimeScale};
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::jpl_ephem::naif::naif_ids::main_belt::AsteroidNumber;
use outfit::jpl_ephem::naif::naif_ids::NaifIds;
use outfit::{EphemerisFrame, JPLEphem, OutfitError};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Pick the ephemeris source.
    //
    // ANISE reads NAIF SPK kernels only, so the source must use the "naif:"
    // token; a "horizon:" (legacy DE binary) source would be rejected.
    let source: EphemFileSource = "naif:DE440".try_into()?;

    // Step 2 — Build the ephemeris handle.
    //
    // `JPLEphem::from_anise` forces the ANISE backend even when
    // `ephem-builtin` is also enabled (unlike `JPLEphem::new`, which prefers
    // ANISE automatically whenever it is compiled in). The kernel is
    // downloaded into the same local cache directory as the other backends
    // on first use.
    let jpl_ephem = JPLEphem::from_anise(source)?;

    // Step 3 — Pick an evaluation epoch and query Earth, exactly as with the
    // other two backends: the public `earth_ephemeris` / `body_ephemeris` API
    // is identical across every backend, only construction differs.
    let epoch = Epoch::from_mjd_in_time_scale(60_000.0, TimeScale::TT);
    let (earth_pos, earth_vel) =
        jpl_ephem.earth_ephemeris(&epoch, EphemerisFrame::Equatorial, true);
    let earth_vel = earth_vel.expect("velocity was requested");
    println!(
        "Earth heliocentric position (AU): x={:.9} y={:.9} z={:.9}",
        earth_pos.x, earth_pos.y, earth_pos.z
    );
    println!(
        "Earth heliocentric velocity (AU/day): x={:.9} y={:.9} z={:.9}",
        earth_vel.x, earth_vel.y, earth_vel.z
    );

    // Step 4 — Try a main-belt asteroid on the plain DE440 handle.
    //
    // `NaifIds::AST(_)` identifies a numbered main-belt asteroid; DE440 does
    // not carry any of them, so this fails cleanly with
    // `EphemerisBodyNotSupported` rather than panicking.
    let ceres = NaifIds::AST(AsteroidNumber::CERES);
    match jpl_ephem.body_ephemeris(ceres, &epoch, EphemerisFrame::Equatorial) {
        Ok(_) => unreachable!("DE440 alone does not carry asteroid ephemerides"),
        Err(err) => println!("\nQuerying {ceres} on plain DE440 failed as expected: {err}"),
    }

    // Step 5 — Load the main-belt asteroid supplementary kernel and retry.
    //
    // `from_anise_with_main_belt_asteroids` loads DE440 plus a second SPK
    // kernel covering 300 numbered asteroids (downloaded once, then cached).
    // Once loaded, `NaifIds::AST(_)` bodies resolve like any other body.
    let source: EphemFileSource = "naif:DE440".try_into()?;
    let jpl_ephem_with_asteroids = JPLEphem::from_anise_with_main_belt_asteroids(source)?;
    let (ceres_pos, ceres_vel) =
        jpl_ephem_with_asteroids.body_ephemeris(ceres, &epoch, EphemerisFrame::Equatorial)?;
    println!(
        "{ceres} heliocentric position (AU): x={:.9} y={:.9} z={:.9}",
        ceres_pos.x, ceres_pos.y, ceres_pos.z
    );
    println!(
        "{ceres} heliocentric velocity (AU/day): x={:.9} y={:.9} z={:.9}",
        ceres_vel.x, ceres_vel.y, ceres_vel.z
    );

    Ok(())
}
