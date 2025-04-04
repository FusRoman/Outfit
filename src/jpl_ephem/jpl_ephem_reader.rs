use std::collections::HashMap;
use std::fs::File as StdFile;
use std::io::{self, Read, Seek, SeekFrom};

use camino::Utf8Path;

const OLDMAX: usize = 400;
const NPLANETS: usize = 12;

#[derive(Debug)]
struct EphemerisHeader {
    titles: Vec<[u8; 6]>,
    constants: Vec<String>,
    start_jd: f64,
    end_jd: f64,
    step_jd: f64,
    num_constants: i32,
    au: f64,
    emrat: f64,
    ipt: [[i32; 3]; NPLANETS],
    numde: i32,
    lpt: [i32; 3],
}

fn read_f64_le<R: Read>(reader: &mut R) -> io::Result<f64> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(f64::from_le_bytes(buf))
}

fn read_i32_le<R: Read>(reader: &mut R) -> io::Result<i32> {
    let mut buf = [0u8; 4];
    reader.read_exact(&mut buf)?;
    Ok(i32::from_le_bytes(buf))
}

fn read_jpl_header(file_path: &Utf8Path) -> io::Result<EphemerisHeader> {
    let mut file = StdFile::open(file_path)?;

    let mut titles = Vec::new();
    for _ in 0..(14 * 3) {
        let mut buf = [0u8; 6];
        file.read_exact(&mut buf)?;
        titles.push(buf);
    }

    let mut constants = Vec::new();
    for _ in 0..OLDMAX {
        let mut buf = [0u8; 6];
        file.read_exact(&mut buf)?;
        let s = String::from_utf8_lossy(&buf).trim().to_string();
        constants.push(s);
    }

    let start_jd = read_f64_le(&mut file)?;
    let end_jd = read_f64_le(&mut file)?;
    let step_jd = read_f64_le(&mut file)?;

    let num_constants = read_i32_le(&mut file)?;
    let au = read_f64_le(&mut file)?;
    let emrat = read_f64_le(&mut file)?;

    let mut ipt = [[0i32; 3]; NPLANETS];
    for i in 0..NPLANETS {
        for j in 0..3 {
            ipt[i][j] = read_i32_le(&mut file)?;
        }
    }

    let numde = read_i32_le(&mut file)?;

    let mut lpt = [0i32; 3];
    for i in 0..3 {
        lpt[i] = read_i32_le(&mut file)?;
    }

    Ok(EphemerisHeader {
        titles,
        constants,
        start_jd,
        end_jd,
        step_jd,
        num_constants,
        au,
        emrat,
        ipt,
        numde,
        lpt,
    })
}

#[derive(Debug)]
struct ChebyshevInterval {
    start_jd: f64,
    end_jd: f64,
    coeffs_x: Vec<f64>,
    coeffs_y: Vec<f64>,
    coeffs_z: Vec<f64>,
}

fn read_chebyshev_data(
    file_path: &str,
    header: &EphemerisHeader,
    record_len: usize,
) -> io::Result<HashMap<u8, Vec<ChebyshevInterval>>> {
    let mut file = StdFile::open(file_path)?;

    let header_bytes = record_len * 8;
    file.seek(SeekFrom::Start(header_bytes as u64))?;

    let mut planet_data: HashMap<u8, Vec<ChebyshevInterval>> = HashMap::new();
    let total_days = header.end_jd - header.start_jd;
    let block_span = header.step_jd;
    let n_blocks = (total_days / block_span).ceil() as usize;

    for block_idx in 0..n_blocks {
        let block_start = header.start_jd + block_idx as f64 * block_span;

        for (planet_id, ipt_data) in header.ipt.iter().enumerate() {
            let n_coeff = ipt_data[0] as usize;
            let n_subintervals = ipt_data[1] as usize;
            let subinterval_len = block_span / n_subintervals as f64;

            if n_coeff == 0 || n_subintervals == 0 {
                continue;
            }

            for i in 0..n_subintervals {
                let start_jd = block_start + i as f64 * subinterval_len;
                let end_jd = start_jd + subinterval_len;

                let mut coeffs_x = Vec::with_capacity(n_coeff);
                let mut coeffs_y = Vec::with_capacity(n_coeff);
                let mut coeffs_z = Vec::with_capacity(n_coeff);

                for _ in 0..n_coeff {
                    coeffs_x.push(read_f64_le(&mut file)?);
                }
                for _ in 0..n_coeff {
                    coeffs_y.push(read_f64_le(&mut file)?);
                }
                for _ in 0..n_coeff {
                    coeffs_z.push(read_f64_le(&mut file)?);
                }

                planet_data
                    .entry(planet_id as u8 + 1)
                    .or_default()
                    .push(ChebyshevInterval {
                        start_jd,
                        end_jd,
                        coeffs_x,
                        coeffs_y,
                        coeffs_z,
                    });
            }
        }
    }

    Ok(planet_data)
}

#[cfg(test)]
mod jpl_reader_test {
    use super::*;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_read_jpl_header() {
        use crate::{jpl_ephem::download_jpl_file::get_ephemeris_file, outfit::Outfit};

        let env_state = Outfit::new();
        let version = Some("de442");
        let user_path = None;

        let file_path = get_ephemeris_file(&env_state, version, user_path).unwrap(); // Replace with your actual file path
        let header = super::read_jpl_header(&file_path).unwrap();
        println!("{:?}", header);
        assert_eq!(header.start_jd, 2451545.0);
        assert_eq!(header.end_jd, 2469800.0);
        assert_eq!(header.step_jd, 0.1);
    }
}
