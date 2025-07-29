mod daf_header;
mod directory;
mod ephemeris_record;
mod jpl_ephem_header;
pub mod naif_data;
pub mod naif_ids;
pub mod naif_version;
mod summary_record;

pub fn print_hex_dump(data: &[u8]) {
    const BYTES_PER_LINE: usize = 40;

    for (i, chunk) in data.chunks(BYTES_PER_LINE).enumerate() {
        // Offset
        print!("{:08x}: ", i * BYTES_PER_LINE);

        // Hex representation
        for byte in chunk {
            print!("{byte:02x} ");
        }

        // Padding if chunk is not full
        if chunk.len() < BYTES_PER_LINE {
            for _ in 0..(BYTES_PER_LINE - chunk.len()) {
                print!("   ");
            }
        }

        // ASCII representation
        print!("|");
        for byte in chunk {
            let ascii = if byte.is_ascii_graphic() || *byte == b' ' {
                *byte as char
            } else {
                '.'
            };
            print!("{ascii}");
        }
        println!("|");
    }
}
