//! Minimal walkthrough of the **recommended, portable entry point**:
//! [`JPLEphem::new`] / `TryFrom<&str>`, which picks the best available
//! backend at compile time without the caller needing to know which one.
//!
//! Compare this with `ephem_horizon.rs`, `ephem_naif.rs`, `ephem_anise.rs`,
//! which force one specific backend — useful to study or compare them, but
//! not the entry point most applications should use.
//!
//! Run it with:
//! ```text
//! cargo run --example ephem_unified
//! ```

use hifitime::{Epoch, TimeScale};
use outfit::jpl_ephem::naif::naif_ids::{planet_bary::PlanetaryBary, NaifIds};
use outfit::{EphemerisFrame, JPLEphem, OutfitError};

fn main() -> Result<(), OutfitError> {
    // Step 1 — Build the ephemeris handle with `JPLEphem::new`.
    //
    // `JPLEphem::new` resolves the source and builds whichever backend was
    // selected at compile time: the ANISE backend when `ephem-anise` is
    // enabled (even if `ephem-builtin` is also on), otherwise the in-house
    // reader. The "naif:" source token works with every backend, so code
    // written against `JPLEphem::new` and a "naif:" source is portable
    // across every feature combination without any `#[cfg(...)]`.
    //
    // `TryFrom<&str>` (used here through `.try_into()?`) is shorthand for the
    // same thing: parse the source string, then call `JPLEphem::new`.
    let jpl_ephem: JPLEphem = "naif:DE440".try_into()?;

    // Step 2 — Pick an evaluation epoch.
    //
    // Every entry point takes a `hifitime::Epoch`; building it directly in
    // Terrestrial Time (TT) avoids any UTC/leap-second handling.
    let epoch = Epoch::from_mjd_in_time_scale(60_000.0, TimeScale::TT);

    // Step 3 — Query Earth and a planet, in both supported frames.
    //
    // `EphemerisFrame::Equatorial` returns the state as stored in the kernel
    // (equatorial mean J2000 / ICRF); `EphemerisFrame::Ecliptic` rotates it
    // to ecliptic mean J2000 by the mean obliquity of the ecliptic at J2000.
    for frame in [EphemerisFrame::Equatorial, EphemerisFrame::Ecliptic] {
        let (earth_pos, _) = jpl_ephem.earth_ephemeris(&epoch, frame, false);
        println!(
            "Earth heliocentric position, {frame:?} (AU): x={:.9} y={:.9} z={:.9}",
            earth_pos.x, earth_pos.y, earth_pos.z
        );
    }

    // Step 4 — Query an arbitrary body through the same, backend-agnostic
    // `body_ephemeris` call.
    let jupiter = NaifIds::PB(PlanetaryBary::Jupiter);
    let (jup_pos, jup_vel) =
        jpl_ephem.body_ephemeris(jupiter, &epoch, EphemerisFrame::Equatorial)?;
    println!(
        "\n{jupiter} heliocentric position (AU): x={:.9} y={:.9} z={:.9}",
        jup_pos.x, jup_pos.y, jup_pos.z
    );
    println!(
        "{jupiter} heliocentric velocity (AU/day): x={:.9} y={:.9} z={:.9}",
        jup_vel.x, jup_vel.y, jup_vel.z
    );

    Ok(())
}
