//! Minimal walkthrough of the **built-in NAIF SPK/DAF** ephemeris backend.
//!
//! This is the second half of Outfit's in-house reader (`ephem-builtin`
//! Cargo feature, the default): a pure-Rust parser for NAIF's SPK/DAF binary
//! kernel format (the modern format SPICE/ANISE use), decoded directly from
//! the DAF header and summary/directory records into Chebyshev segments —
//! no dependency on the `anise` crate.
//!
//! Run it with:
//! ```text
//! cargo run --example ephem_naif
//! ```

use hifitime::{Epoch, TimeScale};
use outfit::jpl_ephem::download_jpl_file::EphemFileSource;
use outfit::jpl_ephem::naif::naif_ids::{solar_system_bary::SolarSystemBary, NaifIds};
use outfit::{EphemerisFrame, JPLEphem, OutfitError};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Pick the ephemeris source.
    //
    // The "naif:" token selects a NAIF SPK/DAF kernel rather than a legacy
    // Horizon binary. "DE440" is the JPL planetary ephemeris version; the
    // same token also works unchanged with the ANISE backend (see
    // `ephem_anise.rs`) — only the reader parsing the bytes differs.
    let source: EphemFileSource = "naif:DE440".try_into()?;

    // Step 2 — Build the ephemeris handle.
    //
    // `JPLEphem::from_builtin` picks the reader from the file kind resolved
    // from `source`: a "naif:" source is parsed as an SPK/DAF kernel here.
    // The kernel is downloaded into a local cache directory on first use.
    let jpl_ephem = JPLEphem::from_builtin(source)?;

    // Step 3 — Pick an evaluation epoch.
    //
    // 60000.0 is MJD TT for 2023-02-25. Both built-in readers accept the
    // same `hifitime::Epoch`; the NAIF reader converts it to ET seconds
    // internally.
    let epoch = Epoch::from_mjd_in_time_scale(60_000.0, TimeScale::TT);

    // Step 4 — Query Earth's heliocentric state.
    //
    // Unlike the Horizon backend, the built-in NAIF reader's DE440 kernel
    // stores the Earth–Moon barycenter relative to the Solar System
    // Barycenter; `earth_ephemeris` derives the Earth geocenter from it
    // internally (using the Earth–Moon mass ratio), so the returned state
    // is heliocentric like every other backend.
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

    // Step 5 — Query the Sun directly.
    //
    // `NaifIds::SSB(SolarSystemBary::Sun)` is the Sun's NAIF id (10). Every
    // backend can resolve it, since it is present in every DE440 file.
    let sun = NaifIds::SSB(SolarSystemBary::Sun);
    let (sun_pos, sun_vel) = jpl_ephem.body_ephemeris(sun, &epoch, EphemerisFrame::Equatorial)?;
    println!(
        "\n{sun} heliocentric position (AU): x={:.9} y={:.9} z={:.9}",
        sun_pos.x, sun_pos.y, sun_pos.z
    );
    println!(
        "{sun} heliocentric velocity (AU/day): x={:.9} y={:.9} z={:.9}",
        sun_vel.x, sun_vel.y, sun_vel.z
    );
    // The Sun is the reference body of every state this API returns, so its
    // own "heliocentric" position and velocity are numerically ~0.

    Ok(())
}
