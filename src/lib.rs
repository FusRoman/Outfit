pub mod constants;
mod conversion;
pub mod env_state;
mod equinoctial_element;
pub mod initial_orbit_determination;
pub mod jpl_ephem;
mod kepler;
mod keplerian_element;
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

    use camino::Utf8Path;

    use crate::{
        constants::TrajectorySet,
        jpl_ephem::{horizon::horizon_data::HorizonData, naif::naif_data::NaifData},
        observations::trajectory_ext::TrajectoryExt,
        outfit::Outfit,
    };

    pub(crate) static OUTFIT_NAIF_TEST: LazyLock<Outfit> =
        LazyLock::new(|| Outfit::new("naif:DE440"));

    pub(crate) static OUTFIT_HORIZON_TEST: LazyLock<(Outfit, TrajectorySet)> =
        LazyLock::new(|| {
            let mut env = Outfit::new("horizon:DE440");

            let path_file = Utf8Path::new("tests/data/2015AB.obs");
            let traj_set = TrajectorySet::new_from_80col(&mut env, &path_file);
            dbg!(traj_set.keys());
            (env, traj_set)
        });

    pub(crate) static JPL_EPHEM_HORIZON: LazyLock<&HorizonData> = LazyLock::new(|| {
        let jpl_ephem = OUTFIT_HORIZON_TEST.0.get_jpl_ephem().unwrap();
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
