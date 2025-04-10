use std::{
    fs::File,
    io::{BufReader, Read, Seek},
};

use nom::number::complete::le_f64;

#[derive(Debug, PartialEq)]
pub struct DirectoryData {
    pub init: f64,
    pub intlen: usize,
    pub rsize: usize,
    pub n_records: usize,
}

impl DirectoryData {
    pub fn parse(file: &mut BufReader<File>, end_addr: usize) -> Self {
        let directory_offset_bytes = (end_addr - 4) * 8;
        let mut dir_buf = [0u8; 32]; // 4 f64 = 32 octets
        file.seek(std::io::SeekFrom::Start(directory_offset_bytes as u64))
            .unwrap();
        file.read_exact(&mut dir_buf).unwrap();
        let (input, init) = le_f64::<_, nom::error::Error<_>>(dir_buf.as_slice()).unwrap();
        let (input, intlen) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        let (input, rsize) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        let (_, n_records) = le_f64::<_, nom::error::Error<_>>(input).unwrap();
        DirectoryData {
            init,
            intlen: intlen as usize,
            rsize: rsize as usize,
            n_records: n_records as usize,
        }
    }
}
