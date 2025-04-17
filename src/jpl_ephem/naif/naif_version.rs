use std::str::FromStr;

#[derive(Debug, Clone)]
pub enum NAIFVersion {
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

impl NAIFVersion {
    pub fn get_filename(&self) -> &str {
        match self {
            NAIFVersion::DE430 => "de430.bsp",
            NAIFVersion::DE431p1 => "de431_part-1.bsp",
            NAIFVersion::DE431p2 => "de431_part-2.bsp",
            NAIFVersion::DE432 => "de432.bsp",
            NAIFVersion::DE435 => "de435.bsp",
            NAIFVersion::DE438 => "de438.bsp",
            NAIFVersion::DE440 => "de440.bsp",
            NAIFVersion::DE440s => "de440s.bsp",
            NAIFVersion::DE441p1 => "de441_part-1.bsp",
            NAIFVersion::DE441p2 => "de441_part-2.bsp",
            NAIFVersion::DE442 => "de442.bsp",
        }
    }

    fn from_str(s: &str) -> Option<Self> {
        match s {
            "DE430" => Some(NAIFVersion::DE430),
            "DE431_part-1" => Some(NAIFVersion::DE431p1),
            "DE431_part-2" => Some(NAIFVersion::DE431p2),
            "DE432" => Some(NAIFVersion::DE432),
            "DE435" => Some(NAIFVersion::DE435),
            "DE438" => Some(NAIFVersion::DE438),
            "DE440" => Some(NAIFVersion::DE440),
            "DE440s" => Some(NAIFVersion::DE440s),
            "DE441_part-1" => Some(NAIFVersion::DE441p1),
            "DE441_part-2" => Some(NAIFVersion::DE441p2),
            "DE442" => Some(NAIFVersion::DE442),
            _ => None,
        }
    }
}

impl FromStr for NAIFVersion {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        NAIFVersion::from_str(s).ok_or_else(|| format!("Invalid NAIF version: {}", s))
    }
}