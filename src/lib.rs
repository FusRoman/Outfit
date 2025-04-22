pub mod constants;
mod conversion;
pub mod env_state;
pub mod initial_orbit_determination;
pub mod jpl_ephem;
mod kepler;
mod keplerian_element;
mod equinoctial_element;
pub mod observations;
pub mod observers;
mod orb_elem;
pub mod outfit;
pub mod outfit_errors;
mod ref_system;
pub mod time;

#[cfg(all(test, feature = "jpl-download"))]
pub(crate) mod unit_test_global {
    use std::sync::LazyLock;

    use crate::{
        jpl_ephem::{horizon::horizon_data::HorizonData, naif::naif_data::NaifData},
        outfit::Outfit,
    };


    pub(crate) static OUTFIT_NAIF_TEST: LazyLock<Outfit> =
        LazyLock::new(|| Outfit::new("naif:DE440"));

    pub(crate) static OUTFIT_HORIZON_TEST: LazyLock<Outfit> =
        LazyLock::new(|| Outfit::new("horizon:DE440"));

    pub(crate) static JPL_EPHEM_HORIZON: LazyLock<&HorizonData> = LazyLock::new(|| {
        let jpl_ephem = OUTFIT_HORIZON_TEST.get_jpl_ephem().unwrap();
        match jpl_ephem {
            crate::jpl_ephem::JPLEphem::HorizonFile(horizon_data) => horizon_data,
            _ => panic!("JPL ephemeris is not a Horizon file"),
        }
    });

    pub(crate) static JPL_EPHEM_NAIF: LazyLock<&NaifData> = LazyLock::new(|| {
        let jpl_ephem = OUTFIT_NAIF_TEST.get_jpl_ephem().unwrap();
        match jpl_ephem {
            crate::jpl_ephem::JPLEphem::NaifFile(naif_data) => naif_data,
            _ => panic!("JPL ephemeris is not a Naif file"),
        }
    });
}
