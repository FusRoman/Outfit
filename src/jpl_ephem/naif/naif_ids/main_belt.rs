use std::fmt;

use super::ErrorId;

/// NAIF encodes a numbered small body as `2_000_000 + <minor-planet number>`.
const NAIF_ID_OFFSET: i32 = 2_000_000;

/// NAIF small-body identifiers occupy `2_000_000..=2_999_999`.
const NAIF_ID_RANGE_END: i32 = 3_000_000;

/// `(asteroid number, name)`, sorted by asteroid number, for the 300 bodies
/// covered by the `codes_300ast_20100725.bsp` supplementary kernel.
///
/// Source: the `NAIF_BODY_NAME` / `NAIF_BODY_CODE` pairs published in that
/// kernel's companion frame kernel
/// (`codes_300ast_20100725.tf`, JPL/NAIF, retrieved 2026-09-11), with each
/// name title-cased from its all-caps original (e.g. `'1 CERES'` → `Ceres`).
/// Covers exactly the same 300 numbers as
/// [`crate::propagator::asteroid_gm_table`]'s mass table — checked by
/// `propagator::planet_gm`'s tests.
#[rustfmt::skip]
const ASTEROID_NAMES: &[(u32, &str)] = &[
    (1, "Ceres"), (2, "Pallas"), (3, "Juno"), (4, "Vesta"),
    (5, "Astraea"), (6, "Hebe"), (7, "Iris"), (8, "Flora"),
    (9, "Metis"), (10, "Hygiea"), (11, "Parthenope"), (12, "Victoria"),
    (13, "Egeria"), (14, "Irene"), (15, "Eunomia"), (16, "Psyche"),
    (17, "Thetis"), (18, "Melpomene"), (19, "Fortuna"), (20, "Massalia"),
    (21, "Lutetia"), (22, "Kalliope"), (23, "Thalia"), (24, "Themis"),
    (25, "Phocaea"), (26, "Proserpina"), (27, "Euterpe"), (28, "Bellona"),
    (29, "Amphitrite"), (30, "Urania"), (31, "Euphrosyne"), (32, "Pomona"),
    (34, "Circe"), (35, "Leukothea"), (36, "Atalante"), (37, "Fides"),
    (38, "Leda"), (39, "Laetitia"), (40, "Harmonia"), (41, "Daphne"),
    (42, "Isis"), (43, "Ariadne"), (44, "Nysa"), (45, "Eugenia"),
    (46, "Hestia"), (47, "Aglaja"), (48, "Doris"), (49, "Pales"),
    (50, "Virginia"), (51, "Nemausa"), (52, "Europa"), (53, "Kalypso"),
    (54, "Alexandra"), (56, "Melete"), (57, "Mnemosyne"), (58, "Concordia"),
    (59, "Elpis"), (62, "Erato"), (63, "Ausonia"), (65, "Cybele"),
    (68, "Leto"), (69, "Hesperia"), (70, "Panopaea"), (71, "Niobe"),
    (72, "Feronia"), (74, "Galatea"), (75, "Eurydike"), (76, "Freia"),
    (77, "Frigga"), (78, "Diana"), (80, "Sappho"), (81, "Terpsichore"),
    (83, "Beatrix"), (84, "Klio"), (85, "Io"), (86, "Semele"),
    (87, "Sylvia"), (88, "Thisbe"), (89, "Julia"), (90, "Antiope"),
    (91, "Aegina"), (92, "Undina"), (93, "Minerva"), (94, "Aurora"),
    (95, "Arethusa"), (96, "Aegle"), (97, "Klotho"), (98, "Ianthe"),
    (99, "Dike"), (102, "Miriam"), (103, "Hera"), (104, "Klymene"),
    (105, "Artemis"), (106, "Dione"), (107, "Camilla"), (109, "Felicitas"),
    (110, "Lydia"), (111, "Ate"), (112, "Iphigenia"), (114, "Kassandra"),
    (115, "Thyra"), (117, "Lomia"), (120, "Lachesis"), (121, "Hermione"),
    (124, "Alkeste"), (127, "Johanna"), (128, "Nemesis"), (129, "Antigone"),
    (130, "Elektra"), (134, "Sophrosyne"), (135, "Hertha"), (137, "Meliboea"),
    (139, "Juewa"), (140, "Siwa"), (141, "Lumen"), (143, "Adria"),
    (144, "Vibilia"), (145, "Adeona"), (146, "Lucina"), (147, "Protogeneia"),
    (148, "Gallia"), (150, "Nuwa"), (154, "Bertha"), (156, "Xanthippe"),
    (159, "Aemilia"), (160, "Una"), (162, "Laurentia"), (163, "Erigone"),
    (164, "Eva"), (165, "Loreley"), (168, "Sibylla"), (171, "Ophelia"),
    (173, "Ino"), (175, "Andromache"), (176, "Iduna"), (181, "Eucharis"),
    (185, "Eunike"), (187, "Lamberta"), (191, "Kolga"), (192, "Nausikaa"),
    (194, "Prokne"), (195, "Eurykleia"), (196, "Philomela"), (200, "Dynamene"),
    (201, "Penelope"), (203, "Pompeja"), (205, "Martha"), (206, "Hersilia"),
    (209, "Dido"), (210, "Isabella"), (211, "Isolda"), (212, "Medea"),
    (213, "Lilaea"), (216, "Kleopatra"), (221, "Eos"), (224, "Oceana"),
    (225, "Henrietta"), (230, "Athamantis"), (233, "Asterope"), (236, "Honoria"),
    (238, "Hypatia"), (240, "Vanadis"), (241, "Germania"), (247, "Eukrate"),
    (250, "Bettina"), (259, "Aletheia"), (266, "Aline"), (268, "Adorea"),
    (275, "Sapientia"), (276, "Adelheid"), (283, "Emma"), (287, "Nephthys"),
    (303, "Josephina"), (304, "Olga"), (308, "Polyxo"), (313, "Chaldaea"),
    (322, "Phaeo"), (324, "Bamberga"), (326, "Tamara"), (328, "Gudrun"),
    (329, "Svea"), (334, "Chicago"), (335, "Roberta"), (336, "Lacadiera"),
    (337, "Devosa"), (338, "Budrosa"), (344, "Desiderata"), (345, "Tercidina"),
    (346, "Hermentaria"), (347, "Pariana"), (349, "Dembowska"), (350, "Ornamenta"),
    (354, "Eleonora"), (356, "Liguria"), (357, "Ninina"), (358, "Apollonia"),
    (360, "Carlova"), (362, "Havnia"), (363, "Padua"), (365, "Corduba"),
    (366, "Vincentina"), (369, "Aeria"), (372, "Palma"), (373, "Melusina"),
    (375, "Ursula"), (377, "Campania"), (381, "Myrrha"), (385, "Ilmatar"),
    (386, "Siegena"), (387, "Aquitania"), (388, "Charybdis"), (389, "Industria"),
    (393, "Lampetia"), (404, "Arsinoe"), (405, "Thia"), (407, "Arachne"),
    (409, "Aspasia"), (410, "Chloris"), (412, "Elisabetha"), (416, "Vaticana"),
    (419, "Aurelia"), (420, "Bertholda"), (423, "Diotima"), (424, "Gratia"),
    (426, "Hippo"), (431, "Nephele"), (433, "Eros"), (442, "Eichsfeldia"),
    (444, "Gyptis"), (449, "Hamburga"), (451, "Patientia"), (454, "Mathesis"),
    (455, "Bruchsalia"), (466, "Tisiphone"), (469, "Argentina"), (471, "Papagena"),
    (476, "Hedwig"), (481, "Emita"), (488, "Kreusa"), (489, "Comacina"),
    (490, "Veritas"), (491, "Carina"), (498, "Tokio"), (505, "Cava"),
    (506, "Marion"), (508, "Princetonia"), (511, "Davida"), (514, "Armida"),
    (521, "Brixia"), (532, "Herculina"), (535, "Montague"), (536, "Merapi"),
    (545, "Messalina"), (554, "Peraga"), (566, "Stereoskopia"), (568, "Cheruskia"),
    (595, "Polyxena"), (596, "Scheila"), (602, "Marianna"), (618, "Elfriede"),
    (626, "Notburga"), (635, "Vundtia"), (654, "Zelinda"), (663, "Gerlinde"),
    (674, "Rachele"), (683, "Lanzia"), (690, "Wratislavia"), (691, "Lehigh"),
    (694, "Ekard"), (702, "Alauda"), (704, "Interamnia"), (705, "Erminia"),
    (709, "Fringilla"), (712, "Boliviana"), (713, "Luscinia"), (739, "Mandeville"),
    (740, "Cantabia"), (747, "Winchester"), (751, "Faina"), (762, "Pulcova"),
    (769, "Tatjana"), (772, "Tanete"), (773, "Irmintraud"), (776, "Berbericia"),
    (780, "Armenia"), (788, "Hohensteina"), (790, "Pretoria"), (791, "Ani"),
    (804, "Hispania"), (814, "Tauris"), (849, "Ara"), (895, "Helio"),
    (909, "Ulla"), (914, "Palisana"), (980, "Anacostia"), (1015, "Christa"),
    (1021, "Flammario"), (1036, "Ganymed"), (1093, "Freda"), (1467, "Mashona"),
];

/// A numbered main-belt asteroid, identified by its official minor-planet
/// number (e.g. `1` for Ceres, `4` for Vesta).
///
/// This wraps the raw number rather than enumerating every known asteroid: the
/// supplementary ephemeris kernel Outfit loads for this body family
/// (`codes_300ast_20100725.bsp`, see
/// [`crate::propagator::planet_gm::known_main_belt_asteroids`]) covers 300
/// bodies, far more than is practical to name one by one.
///
/// Wrapped in [`crate::jpl_ephem::naif::naif_ids::NaifIds::AST`], this id is
/// resolvable only by the ANISE backend with the supplementary kernel loaded
/// — see that enum's "Ephemeris backend support" section.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AsteroidNumber(pub u32);

impl AsteroidNumber {
    /// 1 Ceres.
    pub const CERES: Self = Self(1);
    /// 2 Pallas.
    pub const PALLAS: Self = Self(2);
    /// 4 Vesta.
    pub const VESTA: Self = Self(4);

    /// Build an [`AsteroidNumber`] from a raw NAIF small-body integer code.
    ///
    /// Arguments
    /// -----------------
    /// * `id`: the raw NAIF integer identifier (`2_000_000 + number`).
    ///
    /// Return
    /// ----------
    /// * `Ok(AsteroidNumber)` if `id` falls in the NAIF small-body range
    ///   `2_000_000..=2_999_999`.
    /// * `Err(ErrorId)` otherwise.
    pub fn from_id(id: i32) -> Result<Self, ErrorId> {
        if (NAIF_ID_OFFSET..NAIF_ID_RANGE_END).contains(&id) {
            Ok(Self((id - NAIF_ID_OFFSET) as u32))
        } else {
            Err(ErrorId::InvalidAsteroidId(id))
        }
    }

    /// Convert back to the raw NAIF integer code (`2_000_000 + number`).
    pub fn to_id(self) -> i32 {
        NAIF_ID_OFFSET + self.0 as i32
    }

    /// The asteroid's common name (e.g. `"Ceres"`), if it is one of the 300
    /// bodies covered by `codes_300ast_20100725.bsp`.
    ///
    /// # Returns
    ///
    /// `Some(name)` for a known asteroid, `None` otherwise — a number outside
    /// the embedded table (for example a hypothetical future kernel with
    /// wider coverage) is not itself an error, it simply has no printable
    /// name here.
    pub fn name(self) -> Option<&'static str> {
        ASTEROID_NAMES
            .binary_search_by_key(&self.0, |&(n, _)| n)
            .ok()
            .map(|i| ASTEROID_NAMES[i].1)
    }
}

impl From<AsteroidNumber> for i32 {
    fn from(asteroid: AsteroidNumber) -> Self {
        asteroid.to_id()
    }
}

impl TryFrom<i32> for AsteroidNumber {
    type Error = ErrorId;

    fn try_from(id: i32) -> Result<Self, Self::Error> {
        AsteroidNumber::from_id(id)
    }
}

impl fmt::Display for AsteroidNumber {
    /// Formats as `"<Name> (<number>)"` (e.g. `"Ceres (1)"`) when the name is
    /// known, or `"Asteroid <number>"` otherwise.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.name() {
            Some(name) => write!(f, "{name} ({})", self.0),
            None => write!(f, "Asteroid {}", self.0),
        }
    }
}

#[cfg(test)]
mod test_main_belt {
    use super::*;

    #[test]
    fn from_id_round_trips_to_id() {
        for number in [1_u32, 2, 4, 1467] {
            let asteroid = AsteroidNumber(number);
            assert_eq!(AsteroidNumber::from_id(asteroid.to_id()).unwrap(), asteroid);
        }
    }

    #[test]
    fn known_bodies_map_to_expected_naif_ids() {
        assert_eq!(AsteroidNumber::CERES.to_id(), 2_000_001);
        assert_eq!(AsteroidNumber::PALLAS.to_id(), 2_000_002);
        assert_eq!(AsteroidNumber::VESTA.to_id(), 2_000_004);
    }

    #[test]
    fn from_id_rejects_ids_outside_the_small_body_range() {
        assert!(AsteroidNumber::from_id(0).is_err());
        assert!(AsteroidNumber::from_id(399).is_err());
        assert!(AsteroidNumber::from_id(1_999_999).is_err());
        assert!(AsteroidNumber::from_id(3_000_000).is_err());
    }

    #[test]
    fn display_shows_the_name_when_known() {
        assert_eq!(AsteroidNumber::CERES.to_string(), "Ceres (1)");
        assert_eq!(AsteroidNumber(1467).to_string(), "Mashona (1467)");
    }

    #[test]
    fn display_falls_back_to_the_generic_label_when_unknown() {
        assert_eq!(AsteroidNumber(999_999).to_string(), "Asteroid 999999");
    }

    #[test]
    fn name_resolves_known_bodies_and_rejects_unknown_ones() {
        assert_eq!(AsteroidNumber::CERES.name(), Some("Ceres"));
        assert_eq!(AsteroidNumber::PALLAS.name(), Some("Pallas"));
        assert_eq!(AsteroidNumber::VESTA.name(), Some("Vesta"));
        assert_eq!(AsteroidNumber(999_999).name(), None);
    }

    #[test]
    fn name_table_has_exactly_300_sorted_unique_entries() {
        assert_eq!(ASTEROID_NAMES.len(), 300);
        assert!(ASTEROID_NAMES.windows(2).all(|w| w[0].0 < w[1].0));
    }

    #[test]
    fn try_from_matches_from_id() {
        let via_from_id = AsteroidNumber::from_id(2_000_004);
        let via_try_from: Result<AsteroidNumber, ErrorId> = 2_000_004.try_into();
        assert_eq!(via_from_id.unwrap(), via_try_from.unwrap());
    }
}
