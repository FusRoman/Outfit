use std::str::FromStr;

#[derive(Debug, Clone, PartialEq, PartialOrd)]
pub enum JPLHorizonVersion {
    DE102,
    DE200,
    DE202,
    DE403,
    DE405,
    DE406,
    DE410,
    DE413,
    DE414,
    DE418,
    DE421,
    DE422,
    DE423,
    DE430,
    DE430t,
    DE431,
    DE440,
    DE440t,
    DE441,
}

impl JPLHorizonVersion {
    pub fn get_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "de102/lnxm1410p3002.102",
            JPLHorizonVersion::DE200 => "de200/lnxm1600p2170.200",
            JPLHorizonVersion::DE202 => "de202/lnxp1900p2050.202",
            JPLHorizonVersion::DE403 => "de403/lnxp1600p2200.403",
            JPLHorizonVersion::DE405 => "de405/lnxp1600p2200.405",
            JPLHorizonVersion::DE406 => "de406/lnxm3000p3000.406",
            JPLHorizonVersion::DE410 => "de410/lnxp1960p2020.410",
            JPLHorizonVersion::DE413 => "de413/lnxp1900p2050.413",
            JPLHorizonVersion::DE414 => "de414/lnxp1600p2200.414",
            JPLHorizonVersion::DE418 => "de418/lnxp1900p2050.418",
            JPLHorizonVersion::DE421 => "de421/lnxp1900p2053.421",
            JPLHorizonVersion::DE422 => "de422/lnxm3000p3000.422",
            JPLHorizonVersion::DE423 => "de423/lnxp1800p2200.423",
            JPLHorizonVersion::DE430 => "de430/linux_p1550p2650.430",
            JPLHorizonVersion::DE430t => "de430t/linux_p1550p2650.430t",
            JPLHorizonVersion::DE431 => "de431/lnxm13000p17000.431",
            JPLHorizonVersion::DE440 => "de440/linux_p1550p2650.440",
            JPLHorizonVersion::DE440t => "de440t/linux_p1550p2650.440t",
            JPLHorizonVersion::DE441 => "de441/linux_m13000p17000.441",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE102" => Some(JPLHorizonVersion::DE102),
            "DE200" => Some(JPLHorizonVersion::DE200),
            "DE202" => Some(JPLHorizonVersion::DE202),
            "DE403" => Some(JPLHorizonVersion::DE403),
            "DE405" => Some(JPLHorizonVersion::DE405),
            "DE406" => Some(JPLHorizonVersion::DE406),
            "DE410" => Some(JPLHorizonVersion::DE410),
            "DE413" => Some(JPLHorizonVersion::DE413),
            "DE414" => Some(JPLHorizonVersion::DE414),
            "DE418" => Some(JPLHorizonVersion::DE418),
            "DE421" => Some(JPLHorizonVersion::DE421),
            "DE422" => Some(JPLHorizonVersion::DE422),
            "DE423" => Some(JPLHorizonVersion::DE423),
            "DE430" => Some(JPLHorizonVersion::DE430),
            "DE430t" => Some(JPLHorizonVersion::DE430t),
            "DE431" => Some(JPLHorizonVersion::DE431),
            "DE440" => Some(JPLHorizonVersion::DE440),
            "DE440t" => Some(JPLHorizonVersion::DE440t),
            "DE441" => Some(JPLHorizonVersion::DE441),
            _ => None,
        }
    }

    pub fn to_filename(&self) -> &str {
        match self {
            JPLHorizonVersion::DE102 => "DE102.bsp",
            JPLHorizonVersion::DE200 => "DE200.bsp",
            JPLHorizonVersion::DE202 => "DE202.bsp",
            JPLHorizonVersion::DE403 => "DE403.bsp",
            JPLHorizonVersion::DE405 => "DE405.bsp",
            JPLHorizonVersion::DE406 => "DE406.bsp",
            JPLHorizonVersion::DE410 => "DE410.bsp",
            JPLHorizonVersion::DE413 => "DE413.bsp",
            JPLHorizonVersion::DE414 => "DE414.bsp",
            JPLHorizonVersion::DE418 => "DE418.bsp",
            JPLHorizonVersion::DE421 => "DE421.bsp",
            JPLHorizonVersion::DE422 => "DE422.bsp",
            JPLHorizonVersion::DE423 => "DE423.bsp",
            JPLHorizonVersion::DE430 => "DE430.bsp",
            JPLHorizonVersion::DE430t => "DE430t.bsp",
            JPLHorizonVersion::DE431 => "DE431.bsp",
            JPLHorizonVersion::DE440 => "DE440.bsp",
            JPLHorizonVersion::DE440t => "DE440t.bsp",
            JPLHorizonVersion::DE441 => "DE441.bsp",
        }
    }

    pub fn from_filename(filename: &str) -> Option<Self> {
        match filename {
            "DE102.bsp" => Some(JPLHorizonVersion::DE102),
            "DE200.bsp" => Some(JPLHorizonVersion::DE200),
            "DE202.bsp" => Some(JPLHorizonVersion::DE202),
            "DE403.bsp" => Some(JPLHorizonVersion::DE403),
            "DE405.bsp" => Some(JPLHorizonVersion::DE405),
            "DE406.bsp" => Some(JPLHorizonVersion::DE406),
            "DE410.bsp" => Some(JPLHorizonVersion::DE410),
            "DE413.bsp" => Some(JPLHorizonVersion::DE413),
            "DE414.bsp" => Some(JPLHorizonVersion::DE414),
            "DE418.bsp" => Some(JPLHorizonVersion::DE418),
            "DE421.bsp" => Some(JPLHorizonVersion::DE421),
            "DE422.bsp" => Some(JPLHorizonVersion::DE422),
            "DE423.bsp" => Some(JPLHorizonVersion::DE423),
            "DE430.bsp" => Some(JPLHorizonVersion::DE430),
            "DE430t.bsp" => Some(JPLHorizonVersion::DE430t),
            "DE431.bsp" => Some(JPLHorizonVersion::DE431),
            "DE440.bsp" => Some(JPLHorizonVersion::DE440),
            "DE440t.bsp" => Some(JPLHorizonVersion::DE440t),
            "DE441.bsp" => Some(JPLHorizonVersion::DE441),
            _ => None,
        }
    }
}

impl FromStr for JPLHorizonVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        JPLHorizonVersion::from_str(s).ok_or_else(|| format!("Invalid JPL Horizon version: {s}"))
    }
}
