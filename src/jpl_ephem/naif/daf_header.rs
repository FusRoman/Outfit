use nom::{bytes::complete::take, number::complete::le_i32, IResult};

#[derive(Debug, PartialEq)]
pub struct DAFHeader {
    pub idword: String,
    pub internal_filename: String,
    pub nd: i32,
    pub ni: i32,
    pub fward: i32,
    pub bward: i32,
    pub free: i32,
    pub locfmt: String,
    pub fptstr: String,
}

impl DAFHeader {
    pub(super) fn new() -> Self {
        DAFHeader {
            idword: String::new(),
            internal_filename: String::new(),
            nd: 0,
            ni: 0,
            fward: 0,
            bward: 0,
            free: 0,
            locfmt: String::new(),
            fptstr: String::new(),
        }
    }

    pub fn parse(input: &[u8]) -> IResult<&[u8], Self> {
        let (input, id_word) = take(8usize)(input)?; // "DAF/SPK "
        let (input, nd_bytes) = le_i32(input)?; // ND
        let (input, ni_bytes) = le_i32(input)?; // NI
        let (input, ifname) = take(60usize)(input)?; // internal file name
        let (input, fwd) = le_i32(input)?; // forward ptr
        let (input, bwd) = le_i32(input)?; // backward ptr
        let (input, free) = le_i32(input)?; // first free address
        let (input, locfmt) = take(8usize)(input)?; // location format
        let (input, _) = take(603usize)(input)?; // reserved
        let (input, ftpstr) = take(28usize)(input)?; // ftp string
        Ok((
            input,
            DAFHeader {
                idword: String::from_utf8_lossy(id_word).trim().to_string(),
                internal_filename: String::from_utf8_lossy(ifname).trim().to_string(),
                nd: nd_bytes,
                ni: ni_bytes,
                fward: fwd,
                bward: bwd,
                free: free,
                locfmt: String::from_utf8_lossy(locfmt).trim().to_string(),
                fptstr: String::from_utf8_lossy(ftpstr).trim().to_string(),
            },
        ))
    }
}

use std::fmt;

impl fmt::Display for DAFHeader {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        const LABEL_WIDTH: usize = 18;
        const VALUE_WIDTH: usize = 50;

        let border = format!(
            "+{:-<label$}+{:-<value$}+",
            "",
            "",
            label = LABEL_WIDTH + 1,
            value = VALUE_WIDTH + 1
        );

        writeln!(f, "{}", border)?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "DAF File Header",
            "",
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(f, "{}", border)?;

        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "ID Word",
            format!("{} (Format ID)", self.idword),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Internal Name",
            format!("{}", self.internal_filename),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "ND (doubles)",
            format!("{} double precision summary components", self.nd),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "NI (integers)",
            format!("{} integer summary components", self.ni),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Forward Ptr",
            format!("Record # of first summary: {}", self.fward),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Backward Ptr",
            format!("Record # of last summary: {}", self.bward),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Free Addr",
            format!("Next free address: {}", self.free),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;
        writeln!(
            f,
            "| {:<label$}| {:<value$}|",
            "Binary Format",
            format!("{} (e.g., BIG-IEEE or LTL-IEEE)", self.locfmt),
            label = LABEL_WIDTH,
            value = VALUE_WIDTH
        )?;

        writeln!(f, "{}", border)
    }
}

#[cfg(test)]
mod test_daf_header {
    use super::*;

    #[test]
    fn test_display_daf_header() {
        let header = DAFHeader {
            idword: "DAF/SPK".to_string(),
            internal_filename: "NIO2SPK".to_string(),
            nd: 2,
            ni: 6,
            fward: 62,
            bward: 62,
            free: 14974889,
            locfmt: "LTL-IEEE".to_string(),
            fptstr: "FTPSTR:\r:\n:\r\n:\r\0:�:\u{10}�:ENDFTP".to_string(),
        };

        let expected = r#"+-------------------+---------------------------------------------------+
| DAF File Header   |                                                   |
+-------------------+---------------------------------------------------+
| ID Word           | DAF/SPK (Format ID)                               |
| Internal Name     | NIO2SPK                                           |
| ND (doubles)      | 2 double precision summary components             |
| NI (integers)     | 6 integer summary components                      |
| Forward Ptr       | Record # of first summary: 62                     |
| Backward Ptr      | Record # of last summary: 62                      |
| Free Addr         | Next free address: 14974889                       |
| Binary Format     | LTL-IEEE (e.g., BIG-IEEE or LTL-IEEE)             |
+-------------------+---------------------------------------------------+
"#;
        let output = format!("{}", header);
        assert_eq!(output, expected);
    }
}
