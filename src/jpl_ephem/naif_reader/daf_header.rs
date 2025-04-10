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
