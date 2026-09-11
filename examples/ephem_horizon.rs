//! Minimal walkthrough of the **legacy JPL Horizon** ephemeris backend.
//!
//! This is the oldest of Outfit's ephemeris backends: it reads the binary
//! `.bsp`-named files JPL Horizons distributes for its DE planetary
//! ephemerides (TTL/CNAM/IPT header layout, Chebyshev polynomial segments),
//! not the NAIF SPK/DAF format used by the other two backends. It is part of
//! the in-house reader, enabled by the `ephem-builtin` Cargo feature
//! (the default).
//!
//! Run it with:
//! ```text
//! cargo run --example ephem_horizon
//! ```

use hifitime::{Epoch, TimeScale};
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::jpl_ephem::naif::naif_ids::{planet_bary::PlanetaryBary, NaifIds};
use outfit::{EphemerisFrame, JPLEphem, OutfitError};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Pick the ephemeris source.
    //
    // The source string has the form "<backend>:<DE version>". The Horizon
    // backend is selected with the "horizon:" token; "DE440" is the JPL
    // planetary ephemeris version (the current standard, covering 1550-2650).
    let source: EphemFileSource = "horizon:DE440".try_into()?;

    // Step 2 — Build the ephemeris handle.
    //
    // `JPLEphem::from_builtin` forces the in-house reader even when the
    // `ephem-anise` feature is also enabled (unlike `JPLEphem::new`, which
    // would prefer ANISE in that case). The first call resolves the file:
    // it is downloaded into a local cache directory if not already present,
    // then parsed into per-body Chebyshev segments.
    let jpl_ephem = JPLEphem::from_builtin(source)?;

    // Step 3 — Pick an evaluation epoch.
    //
    // Outfit's ephemeris entry points take a `hifitime::Epoch`. The Horizon
    // backend internally works in MJD (Terrestrial Time), so building the
    // epoch directly in that time scale avoids any UTC/leap-second handling.
    // 60000.0 is MJD TT for 2023-02-25.
    let epoch = Epoch::from_mjd_in_time_scale(60_000.0, TimeScale::TT);

    // Step 4 — Query Earth's heliocentric state.
    //
    // `earth_ephemeris` always resolves Earth (no body identifier needed).
    // It returns AU / AU/day, equatorial mean J2000, heliocentric (Earth
    // geocenter relative to the Sun for this backend). `true` also asks for
    // the velocity; passing `false` would return `None` for it.
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

    // Step 5 — Query an arbitrary body: Mars' barycenter.
    //
    // `body_ephemeris` takes a `NaifIds` identifier so the same call works
    // unchanged across every backend. It is fallible: a body absent from
    // the loaded kernel — or unsupported by this backend — surfaces as
    // `OutfitError::EphemerisBodyNotSupported` rather than panicking.
    let mars = NaifIds::PB(PlanetaryBary::Mars);
    let (mars_pos, mars_vel) =
        jpl_ephem.body_ephemeris(mars, &epoch, EphemerisFrame::Equatorial)?;
    println!(
        "\n{mars} heliocentric position (AU): x={:.9} y={:.9} z={:.9}",
        mars_pos.x, mars_pos.y, mars_pos.z
    );
    println!(
        "{mars} heliocentric velocity (AU/day): x={:.9} y={:.9} z={:.9}",
        mars_vel.x, mars_vel.y, mars_vel.z
    );

    Ok(())
}
