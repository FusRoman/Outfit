use std::str::FromStr;

#[derive(Debug, Clone)]
pub enum NaifVersion {
    DE430,
    DE431p1,
    DE431p2,
    DE432,
    DE435,
    DE438,
    DE440,
    DE440s,
    DE441p1,
    DE441p2,
    DE442,
}

impl NaifVersion {
    pub fn get_filename(&self) -> &str {
        match self {
            NaifVersion::DE430 => "de430.bsp",
            NaifVersion::DE431p1 => "de431_part-1.bsp",
            NaifVersion::DE431p2 => "de431_part-2.bsp",
            NaifVersion::DE432 => "de432.bsp",
            NaifVersion::DE435 => "de435.bsp",
            NaifVersion::DE438 => "de438.bsp",
            NaifVersion::DE440 => "de440.bsp",
            NaifVersion::DE440s => "de440s.bsp",
            NaifVersion::DE441p1 => "de441_part-1.bsp",
            NaifVersion::DE441p2 => "de441_part-2.bsp",
            NaifVersion::DE442 => "de442.bsp",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE430" => Some(NaifVersion::DE430),
            "DE431_part-1" => Some(NaifVersion::DE431p1),
            "DE431_part-2" => Some(NaifVersion::DE431p2),
            "DE432" => Some(NaifVersion::DE432),
            "DE435" => Some(NaifVersion::DE435),
            "DE438" => Some(NaifVersion::DE438),
            "DE440" => Some(NaifVersion::DE440),
            "DE440s" => Some(NaifVersion::DE440s),
            "DE441_part-1" => Some(NaifVersion::DE441p1),
            "DE441_part-2" => Some(NaifVersion::DE441p2),
            "DE442" => Some(NaifVersion::DE442),
            _ => None,
        }
    }
}

impl FromStr for NaifVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        NaifVersion::from_str(s).ok_or_else(|| format!("Invalid NAIF version: {}", s))
    }
}
