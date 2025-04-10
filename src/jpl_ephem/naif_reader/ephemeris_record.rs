use std::{
    fs::File,
    io::{BufReader, Read, Seek},
};

#[derive(Debug, Clone, PartialEq)]
pub struct EphemerisRecord {
    pub mid: f64,
    pub radius: f64,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

impl EphemerisRecord {
    fn parse_ephemeris_record(
        input: &[u8],
        ncoeff: usize,
    ) -> Result<EphemerisRecord, Box<dyn std::error::Error>> {
        let mut offset = 0;

        // Helper closure pour lire un f64 et avancer dans le buffer
        let mut read_f64 = || -> Result<f64, Box<dyn std::error::Error>> {
            let val = f64::from_le_bytes(input[offset..offset + 8].try_into().unwrap());
            offset += 8;
            Ok(val)
        };

        let mid = read_f64()?;
        let radius = read_f64()?;

        let mut x = Vec::with_capacity(ncoeff);
        let mut y = Vec::with_capacity(ncoeff);
        let mut z = Vec::with_capacity(ncoeff);

        for _ in 0..ncoeff {
            x.push(read_f64()?);
        }
        for _ in 0..ncoeff {
            y.push(read_f64()?);
        }
        for _ in 0..ncoeff {
            z.push(read_f64()?);
        }

        Ok(EphemerisRecord {
            mid,
            radius,
            x,
            y,
            z,
        })
    }

    pub fn parse(
        file: &mut BufReader<File>,
        segment_start_addr: usize,
        rsize: usize,
        n_records: usize,
    ) -> Vec<Self> {
        let mut records = Vec::with_capacity(n_records);

        let record_byte_size = rsize * 8;
        let mut buf = vec![0u8; record_byte_size];

        // L'offset du premier record
        let start_byte_offset = (segment_start_addr - 1) * 8;

        let n_coeffs = (rsize - 2) / 3;

        for i in 0..n_records {
            let byte_offset = start_byte_offset + i * record_byte_size;

            file.seek(std::io::SeekFrom::Start(byte_offset as u64))
                .unwrap();
            file.read_exact(&mut buf).unwrap();

            let input = buf.as_slice();

            // MID, RADIUS
            let record = Self::parse_ephemeris_record(input, n_coeffs).unwrap();

            records.push(record);
        }

        records
    }
}
