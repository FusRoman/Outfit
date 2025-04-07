use camino::Utf8Path;
use regex::Regex;
use std::fs::File;
use std::io::{self, Read, Seek, SeekFrom};

fn open_and_read_record(
    file_path: &str,
    record_number: u64,
    record_length: usize,
) -> io::Result<Vec<u8>> {
    // 1. Ouvrir le fichier (équivalent de status='old')
    let mut file = File::open(file_path)?;

    // 2. Calculer l’offset du record (record_number commence à 1 comme en Fortran)
    let offset = (record_number - 1) * record_length as u64;
    file.seek(SeekFrom::Start(offset))?;

    // 3. Lire exactement record_length octets (comme un record)
    let mut buffer = vec![0u8; record_length];
    file.read_exact(&mut buffer)?;

    Ok(buffer)
}

fn skip_ftp_header(buffer: &[u8]) -> &[u8] {
    let marker = b"ENDFTP";
    if let Some(pos) = buffer.windows(marker.len()).position(|w| w == marker) {
        // Avance après le marqueur ENDFTP
        let mut start = pos + marker.len();

        // Aligne sur le prochain multiple de 1024 (ou autre taille de record si connue)
        let alignment = 1024; // à ajuster si besoin
        start = ((start + alignment - 1) / alignment) * alignment;

        if start < buffer.len() {
            &buffer[start..]
        } else {
            &[]
        }
    } else {
        // Pas de header détecté, on retourne tout
        buffer
    }
}

#[derive(Debug)]
struct JPLEphemHeader {
    version: String,
    creation_date: String,
    start_ephem: String,
    end_ephem: String,
    start_jd: f64,
    end_jd: f64,
}

fn read_jpl_header(jpl_file_path: &Utf8Path) -> JPLEphemHeader {
    let record = open_and_read_record(jpl_file_path.as_str(), 1, 4096).unwrap();

    let without_ftp_header = skip_ftp_header(&record);

    // Essayons d'afficher uniquement les caractères imprimables
    let ascii_text = String::from_utf8_lossy(&without_ftp_header).to_string();

    let fake_lines = ascii_text
        .split(|c| c == '\r' || c == '\n' || c == 0 as char) // 0 = NULL, souvent dans binaires
        .map(|s| s.trim())
        .filter(|s| !s.is_empty());

    let mut version = String::new();
    let mut creation_date = String::new();
    let mut start_ephem = String::new();
    let mut end_ephem = String::new();
    let mut start_jd = 0.0;
    let mut end_jd = 0.0;

    let ephem_date_span = Regex::new(
        r"(?-u)(\d{2}-[A-Z]{3}-\d{4} \d{2}:\d{2})\s+to\s+(\d{2}-[A-Z]{3}-\d{4} \d{2}:\d{2})",
    )
    .unwrap();

    let ephem_jd_span =
        Regex::new(r"(?-u)JD\s+([0-9]+\.[0-9]+)\s+to\s+JD\s+([0-9]+\.[0-9]+)").unwrap();
    for line in fake_lines {
        if version.is_empty() && line.to_lowercase().contains("ephemeris") {
            version = line
                .trim()
                .to_string()
                .split("ephemeris")
                .nth(1)
                .expect("Failed to split string to get ephemeris version")
                .trim()
                .to_string();
        }

        if creation_date.is_empty() && line.to_lowercase().contains("integrated") {
            creation_date = line
                .trim()
                .to_string()
                .split("ed")
                .nth(1)
                .expect("Failed to split string to get creation date")
                .trim()
                .to_string();
        }

        if start_ephem.is_empty() && end_ephem.is_empty() {
            if let Some(caps) = ephem_date_span.captures(line) {
                start_ephem = caps[1].to_string();
                end_ephem = caps[2].to_string();
            }
        }

        if start_jd == 0. && end_jd == 0. {
            if let Some(caps) = ephem_jd_span.captures(line) {
                start_jd = caps[1].parse().unwrap();
                end_jd = caps[2].parse().unwrap();
            }
        }
    }
    JPLEphemHeader {
        version,
        creation_date,
        start_ephem,
        end_ephem,
        start_jd,
        end_jd,
    }
}

#[cfg(test)]
mod jpl_reader_test {
    use super::*;
    use crate::jpl_ephem::download_jpl_file::get_ephemeris_file;
    use crate::outfit::Outfit;

    #[test]
    #[cfg(feature = "jpl-download")]
    fn test_read_jpl_header() {
        let env_state = Outfit::new();
        let version = Some("de440");
        let user_path = None;

        let file_path = get_ephemeris_file(&env_state, version, user_path).unwrap();
        let jpl_header = read_jpl_header(&file_path);
        assert_eq!(jpl_header.version, "DE440");
        assert_eq!(jpl_header.creation_date, "25 June 2020");
        assert_eq!(jpl_header.start_ephem, "31-DEC-1549 00:00");
        assert_eq!(jpl_header.end_ephem, "25-JAN-2650 00:00");
        assert_eq!(jpl_header.start_jd, 2287184.5);
        assert_eq!(jpl_header.end_jd, 2688976.5);
    }
}
