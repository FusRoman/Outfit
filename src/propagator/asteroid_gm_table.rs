//! Gravitational parameters for the 300 main-belt asteroids covered by the
//! `codes_300ast_20100725.bsp` supplementary ephemeris kernel.
//!
//! # Source
//!
//! The table below reproduces, verbatim, the `mu.txt` mass file published in
//! that kernel's comment file
//! (`https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/asteroids/codes_300ast_20100725.cmt`,
//! retrieved 2026-09-11): 300 lines of `<asteroid number> <GM ratio>`, one per
//! body covered by the kernel. Each ratio is `GM_body / GM_sun` (dimensionless,
//! the convention used by the CODES integration that produced the kernel).
//! Cross-checked against the kernel's own segment list
//! (`aa_summaries.txt` in the same directory): the 300 numbers here are
//! exactly the 300 targets stored in the kernel, no more, no fewer.

use crate::propagator::planet_gm::GM_SUN;

/// `(asteroid number, GM_body / GM_sun)`, sorted by asteroid number.
///
/// Values are dimensionless mass ratios as published in `mu.txt` (see the
/// module documentation for provenance); [`asteroid_gm_au3_day2`] converts a
/// ratio to Outfit's AU³/day² convention using [`GM_SUN`], so the whole crate
/// shares one solar GM.
#[rustfmt::skip]
pub(crate) const ASTEROID_GM_RATIOS: &[(u32, f64)] = &[
    (1, 4.760e-10), (2, 1.060e-10), (3, 1.500e-11), (4, 1.340e-10), (5, 1.260e-12),
    (6, 6.300e-12), (7, 7.300e-12), (8, 3.400e-12), (9, 7.400e-12), (10, 4.400e-11),
    (11, 2.290e-12), (12, 1.020e-12), (13, 7.500e-12), (14, 3.500e-12), (15, 1.590e-11),
    (16, 1.140e-11), (17, 7.000e-13), (18, 1.510e-12), (19, 4.190e-12), (20, 1.900e-12),
    (21, 1.300e-12), (22, 4.070e-12), (23, 8.840e-13), (24, 4.000e-12), (25, 3.010e-13),
    (26, 6.110e-13), (27, 1.600e-12), (28, 1.260e-12), (29, 7.800e-12), (30, 9.510e-13),
    (31, 2.700e-11), (32, 3.740e-13), (34, 5.160e-13), (35, 3.870e-13), (36, 4.150e-13),
    (37, 9.040e-13), (38, 5.500e-13), (39, 2.800e-12), (40, 8.860e-13), (41, 1.860e-12),
    (42, 7.150e-13), (43, 2.030e-13), (44, 2.510e-13), (45, 2.860e-12), (46, 6.740e-13),
    (47, 1.090e-12), (48, 1.300e-11), (49, 1.540e-12), (50, 3.510e-13), (51, 1.140e-12),
    (52, 1.320e-11), (53, 5.420e-13), (54, 1.610e-12), (56, 5.120e-13), (57, 1.010e-12),
    (58, 2.880e-13), (59, 1.580e-12), (62, 3.060e-13), (63, 7.800e-13), (65, 6.900e-12),
    (68, 1.310e-12), (69, 3.680e-12), (70, 6.430e-13), (71, 4.130e-13), (72, 2.250e-13),
    (74, 5.900e-13), (75, 2.410e-13), (76, 2.190e-12), (77, 4.630e-13), (78, 6.190e-13),
    (80, 3.420e-13), (81, 5.960e-13), (83, 7.520e-13), (84, 1.750e-13), (85, 1.310e-12),
    (86, 6.180e-13), (87, 7.430e-12), (88, 9.100e-12), (89, 2.470e-12), (90, 4.170e-13),
    (91, 4.670e-13), (92, 2.820e-12), (93, 9.880e-13), (94, 3.030e-12), (95, 8.880e-13),
    (96, 1.730e-12), (97, 7.930e-13), (98, 4.020e-13), (99, 1.310e-13), (102, 2.020e-13),
    (103, 5.390e-13), (104, 6.670e-13), (105, 5.960e-13), (106, 1.110e-12), (107, 5.630e-12),
    (109, 2.520e-13), (110, 8.900e-13), (111, 3.060e-12), (112, 1.330e-13), (114, 3.490e-13),
    (115, 3.620e-13), (117, 1.160e-12), (120, 1.860e-12), (121, 2.710e-12), (124, 3.160e-13),
    (127, 6.620e-13), (128, 2.350e-12), (129, 2.720e-12), (130, 3.320e-12), (134, 6.610e-13),
    (135, 6.940e-13), (137, 1.080e-12), (139, 1.350e-12), (140, 4.670e-13), (141, 7.940e-13),
    (143, 2.570e-13), (144, 1.000e-12), (145, 1.220e-12), (146, 8.150e-13), (147, 8.290e-13),
    (148, 3.290e-13), (150, 1.220e-12), (154, 2.230e-12), (156, 6.250e-13), (159, 6.880e-13),
    (160, 1.890e-13), (162, 6.920e-13), (163, 1.350e-13), (164, 4.070e-13), (165, 1.320e-12),
    (168, 1.150e-12), (171, 5.600e-13), (173, 1.290e-12), (175, 3.630e-13), (176, 6.260e-13),
    (181, 8.460e-13), (185, 1.380e-12), (187, 7.980e-13), (191, 3.640e-13), (192, 7.830e-13),
    (194, 1.690e-12), (195, 2.220e-13), (196, 1.800e-12), (200, 7.460e-13), (201, 4.460e-13),
    (203, 5.540e-13), (205, 1.870e-13), (206, 5.090e-13), (209, 1.440e-12), (210, 2.290e-13),
    (211, 1.040e-12), (212, 8.900e-13), (213, 2.020e-13), (216, 3.440e-12), (221, 7.960e-13),
    (224, 3.300e-13), (225, 6.170e-13), (230, 9.200e-13), (233, 3.830e-13), (236, 4.550e-13),
    (238, 1.150e-12), (240, 3.960e-13), (241, 1.700e-12), (247, 8.570e-13), (250, 7.080e-13),
    (259, 2.010e-12), (266, 4.580e-13), (268, 9.660e-13), (275, 5.360e-13), (276, 6.340e-13),
    (283, 6.940e-13), (287, 2.200e-13), (303, 3.450e-13), (304, 1.100e-13), (308, 9.820e-13),
    (313, 3.150e-13), (322, 4.960e-13), (324, 5.100e-12), (326, 2.840e-13), (328, 1.320e-12),
    (329, 1.660e-13), (334, 1.330e-12), (335, 2.490e-13), (336, 1.170e-13), (337, 2.880e-13),
    (338, 3.510e-13), (344, 8.160e-13), (345, 2.940e-13), (346, 8.590e-13), (347, 1.880e-13),
    (349, 1.940e-12), (350, 5.850e-13), (354, 2.660e-12), (356, 7.990e-13), (357, 4.210e-13),
    (358, 2.520e-13), (360, 5.470e-13), (362, 3.320e-13), (363, 1.340e-13), (365, 4.190e-13),
    (366, 2.910e-13), (369, 3.010e-13), (372, 2.370e-12), (373, 3.100e-13), (375, 3.550e-12),
    (377, 2.660e-13), (381, 6.180e-13), (385, 5.450e-13), (386, 1.580e-12), (387, 7.220e-13),
    (388, 5.250e-13), (389, 3.500e-13), (393, 3.210e-13), (404, 3.290e-13), (405, 6.870e-13),
    (407, 3.030e-13), (409, 1.490e-12), (410, 6.650e-13), (412, 2.650e-13), (416, 4.440e-13),
    (419, 7.570e-13), (420, 9.940e-13), (423, 3.210e-12), (424, 2.340e-13), (426, 7.240e-13),
    (431, 3.030e-13), (433, 3.360e-15), (442, 1.000e-13), (444, 3.600e-12), (449, 2.210e-13),
    (451, 4.020e-12), (454, 1.910e-13), (455, 2.120e-13), (466, 5.440e-13), (469, 6.980e-13),
    (471, 1.720e-12), (476, 5.610e-13), (481, 5.120e-13), (488, 1.190e-12), (489, 9.560e-13),
    (490, 1.330e-12), (491, 3.250e-13), (498, 7.910e-13), (505, 5.360e-13), (506, 4.190e-13),
    (508, 1.020e-12), (511, 1.930e-11), (514, 4.220e-13), (521, 5.460e-13), (532, 8.680e-12),
    (535, 1.460e-13), (536, 1.220e-12), (545, 4.860e-13), (554, 3.110e-13), (566, 1.680e-12),
    (568, 2.320e-13), (595, 4.580e-13), (596, 5.140e-13), (602, 6.840e-13), (618, 6.140e-13),
    (626, 3.610e-13), (635, 3.340e-13), (654, 7.290e-13), (663, 3.620e-13), (674, 6.560e-13),
    (683, 1.940e-13), (690, 8.690e-13), (691, 2.380e-13), (694, 2.640e-13), (702, 3.040e-12),
    (704, 1.960e-11), (705, 8.530e-13), (709, 3.180e-13), (712, 7.320e-13), (713, 4.140e-13),
    (739, 4.370e-13), (740, 2.650e-13), (747, 1.790e-12), (751, 4.760e-13), (762, 7.040e-13),
    (769, 4.250e-13), (772, 5.750e-13), (773, 3.110e-13), (776, 1.220e-12), (780, 2.970e-13),
    (788, 3.930e-13), (790, 1.740e-12), (791, 3.910e-13), (804, 2.060e-12), (814, 4.640e-13),
    (849, 3.300e-13), (895, 1.010e-12), (909, 5.570e-13), (914, 1.590e-13), (980, 4.550e-13),
    (1015, 3.210e-13), (1021, 3.460e-13), (1036, 2.260e-14), (1093, 5.610e-13), (1467, 4.960e-13),
];

/// GM, in AU³/day², of a numbered main-belt asteroid.
///
/// # Arguments
///
/// * `number` – official minor-planet number (e.g. `1` for Ceres).
///
/// # Returns
///
/// `Some(gm)` if `number` is one of the 300 bodies covered by
/// `codes_300ast_20100725.bsp`, `None` otherwise.
pub(crate) fn asteroid_gm_au3_day2(number: u32) -> Option<f64> {
    ASTEROID_GM_RATIOS
        .binary_search_by_key(&number, |&(n, _)| n)
        .ok()
        .map(|i| ASTEROID_GM_RATIOS[i].1 * GM_SUN)
}

#[cfg(test)]
mod asteroid_gm_table_tests {
    use super::*;

    #[test]
    fn table_has_exactly_300_sorted_unique_entries() {
        assert_eq!(ASTEROID_GM_RATIOS.len(), 300);
        assert!(ASTEROID_GM_RATIOS.windows(2).all(|w| w[0].0 < w[1].0));
    }

    #[test]
    fn known_bodies_match_published_ratios() {
        // Ratios are exactly the published mu.txt values (see module docs).
        let ceres_ratio = ASTEROID_GM_RATIOS.iter().find(|&&(n, _)| n == 1).unwrap().1;
        assert!((ceres_ratio - 4.76e-10).abs() < 1e-15);

        // Converted GM is within ~2% of the independently-sourced JPL Horizons
        // value for Ceres (62.6284 km^3/s^2), a sanity check that the ratio
        // convention (GM_body / GM_sun) and GM_SUN scale are both correct.
        const AU_KM: f64 = 1.495_978_707e8;
        const KM3_S2_TO_AU3_DAY2: f64 = (86400.0 * 86400.0) / (AU_KM * AU_KM * AU_KM);
        let ceres_km3_s2 = asteroid_gm_au3_day2(1).unwrap() / KM3_S2_TO_AU3_DAY2;
        let horizons_ceres_km3_s2 = 62.6284;
        let rel = (ceres_km3_s2 - horizons_ceres_km3_s2).abs() / horizons_ceres_km3_s2;
        assert!(rel < 0.02, "Ceres GM relative diff from Horizons: {rel:.3}");
    }

    #[test]
    fn unknown_number_returns_none() {
        assert_eq!(asteroid_gm_au3_day2(999_999), None);
        assert_eq!(asteroid_gm_au3_day2(3_000_000), None);
    }

    #[test]
    fn every_entry_yields_a_finite_positive_gm() {
        for &(number, _) in ASTEROID_GM_RATIOS {
            let gm = asteroid_gm_au3_day2(number).unwrap();
            assert!(gm.is_finite() && gm > 0.0, "asteroid {number}: gm={gm}");
        }
    }
}
