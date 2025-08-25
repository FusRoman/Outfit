//! Ephemeris record decoding and Chebyshev interpolation.
//!
//! This module defines [`EphemerisRecord`], a compact container for a single
//! SPK/DAF *ephemeris record* composed of the record midpoint `mid` (ET seconds),
//! the half‑interval `radius` (seconds), and three arrays of Chebyshev
//! coefficients (`x`, `y`, `z`) for position (km). It also provides:
//!
//! * a binary parser that reads contiguous records from a DAF segment, and
//! * an interpolation routine that evaluates position and velocity at an ET epoch.
//!
//! ## Record layout (per axis)
//! Each record contains, in little‑endian `f64`:
//! 1. `mid` (ET seconds from J2000 TDB),
//! 2. `radius` (seconds),
//! 3. `ncoeff` Chebyshev coefficients for X (km),
//! 4. `ncoeff` for Y (km),
//! 5. `ncoeff` for Z (km).
//!
//! The normalized time is `t = (et - mid) / radius`, clamped to `[-1, 1]`.
//! Position uses `T_n(t)`; velocity uses `T'_n(t)` scaled by `2/radius`.
//!
//! ## Units & time scales
//! * `mid` and `radius` are in **seconds** of ET/TDB (consistent with SPK).
//! * Interpolated **position** is in **kilometers**; **velocity** in **km/s**.
//!
//! ## See also
//! ------------
//! * `daf_header` – header and ND/NI info needed to size records.
//! * `directory` – directory/footer providing `rsize` and `n_records`.
//! * NAIF SPK Required Reading – authoritative spec for record structure.
use std::{
    fmt,
    fs::File,
    io::{BufReader, Read, Seek},
};

use hifitime::{Duration, Epoch};
use nalgebra::Vector3;

/// One SPK ephemeris record (midpoint, half‑width, and Chebyshev coefficients).
///
/// The `x`, `y`, and `z` arrays hold Chebyshev coefficients for the position
/// (in kilometers). `mid` is the center ET epoch of the record; `radius` is
/// the half‑interval length in seconds.
///
/// See also
/// ------------
/// * [`EphemerisRecord::parse`] – Bulk reader for contiguous records in a segment.
/// * [`EphemerisRecord::interpolate`] – Evaluate position and velocity at a given ET.
#[derive(Debug, Clone, PartialEq)]
pub struct EphemerisRecord {
    /// Midpoint of the record time span (ET seconds from J2000 TDB).
    pub mid: f64,
    /// Half‑width of the record interval (seconds).
    pub radius: f64,
    /// Chebyshev coefficients for X (km).
    pub x: Vec<f64>,
    /// Chebyshev coefficients for Y (km).
    pub y: Vec<f64>,
    /// Chebyshev coefficients for Z (km).
    pub z: Vec<f64>,
}

impl EphemerisRecord {
    /// Decode one ephemeris record from a raw byte slice.
    ///
    /// The slice **must** contain exactly the serialized layout:
    /// `mid`, `radius`, then `ncoeff` `f64` for each axis X, Y, Z (little‑endian).
    ///
    /// Arguments
    /// -----------------
    /// * `input`: Byte slice positioned at the start of the record.
    /// * `ncoeff`: Number of Chebyshev coefficients per axis.
    ///
    /// Return
    /// ----------
    /// * `Result<EphemerisRecord, Box<dyn std::error::Error>>` with the decoded record.
    ///
    /// See also
    /// ------------
    /// * [`EphemerisRecord::parse`] – Reads many records directly from a file segment.
    /// * NAIF SPK docs – Record/segment sizes and addressing rules.
    fn parse_ephemeris_record(
        input: &[u8],
        ncoeff: usize,
    ) -> Result<EphemerisRecord, Box<dyn std::error::Error>> {
        let mut offset = 0;

        // Helper to read one LE f64 and advance the cursor
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

    /// Read `n_records` contiguous records from a DAF segment.
    ///
    /// Records are assumed to be tightly packed, starting at
    /// `segment_start_addr` (DAF address in double‑precision words, 1‑based),
    /// with a fixed record size `rsize` (in double‑precision words).
    ///
    /// Arguments
    /// -----------------
    /// * `file`: Buffered file reader to seek and read from.
    /// * `segment_start_addr`: Start address of the segment in **DP-words** (1‑based).
    /// * `rsize`: Record size in **DP-words** (each word = 8 bytes).
    /// * `n_records`: Number of records to read.
    ///
    /// Return
    /// ----------
    /// * `Vec<EphemerisRecord>` containing the decoded records in order.
    ///
    /// See also
    /// ------------
    /// * `directory` module – Provides `rsize` and `n_records` from the footer.
    /// * `daf_header` – For ND/NI and navigation across summary records.
    pub fn parse(
        file: &mut BufReader<File>,
        segment_start_addr: usize,
        rsize: usize,
        n_records: usize,
    ) -> Vec<Self> {
        let mut records = Vec::with_capacity(n_records);

        let record_byte_size = rsize * 8;
        let mut buf = vec![0u8; record_byte_size];

        // Byte offset of first record: (DAF address is 1‑based DP-words)
        let start_byte_offset = (segment_start_addr - 1) * 8;

        // Two header scalars (mid, radius), then 3 * ncoeff coefficients.
        let n_coeffs = (rsize - 2) / 3;

        for i in 0..n_records {
            let byte_offset = start_byte_offset + i * record_byte_size;

            file.seek(std::io::SeekFrom::Start(byte_offset as u64))
                .unwrap();
            file.read_exact(&mut buf).unwrap();

            let input = buf.as_slice();

            let record = Self::parse_ephemeris_record(input, n_coeffs).unwrap();
            records.push(record);
        }

        records
    }

    /// Interpolate Cartesian **position** and **velocity** at an ET epoch.
    ///
    /// The time is first normalized to `t = (et - mid) / radius` and clamped to
    /// `[-1, 1]`. Position uses Chebyshev polynomials `T_n(t)`. Velocity uses
    /// the derivatives `T'_n(t)` and is scaled by `2 / radius`.
    ///
    /// Arguments
    /// -----------------
    /// * `ephem_time`: Target epoch in ET seconds (TDB) at which to interpolate.
    ///
    /// Return
    /// ----------
    /// * `(Vector3<f64>, Vector3<f64>)` where the first vector is **position \[km\]**
    ///   and the second is **velocity [km/s]**.
    ///
    /// See also
    /// ------------
    /// * [`Self::parse`] – To load records before interpolation.
    /// * NAIF SPK docs – Chebyshev series for states and the `2/radius` factor.
    pub fn interpolate(&self, ephem_time: f64) -> (Vector3<f64>, Vector3<f64>) {
        let normalized_time = (ephem_time - self.mid) / self.radius;
        let clamped_time = normalized_time.clamp(-1.0, 1.0);

        let num_coefficients = self.x.len();
        let mut chebyshev_polynomials = vec![0.0; num_coefficients];
        chebyshev_polynomials[0] = 1.0;

        // T_0(t) = 1, T_1(t) = t, T_n(t) = 2 t T_{n-1}(t) - T_{n-2}(t)
        if num_coefficients > 1 {
            chebyshev_polynomials[1] = clamped_time;
            for degree in 2..num_coefficients {
                chebyshev_polynomials[degree] =
                    2.0 * clamped_time * chebyshev_polynomials[degree - 1]
                        - chebyshev_polynomials[degree - 2];
            }
        }

        // Position as linear combination of T_n(t)
        let position = Vector3::new(
            self.x
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(c, p)| c * p)
                .sum(),
            self.y
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(c, p)| c * p)
                .sum(),
            self.z
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(c, p)| c * p)
                .sum(),
        );

        // Velocity using T'_n(t) and scale 2/radius
        let mut velocity = Vector3::zeros();
        if num_coefficients > 1 {
            let mut chebyshev_derivatives = vec![0.0; num_coefficients];
            chebyshev_derivatives[1] = 1.0;

            if num_coefficients > 2 {
                chebyshev_derivatives[2] = 4.0 * clamped_time;
                for degree in 3..num_coefficients {
                    // Recurrence for derivatives:
                    // T'_n = 2 t T'_{n-1} + 2 T_{n-1} - T'_{n-2}
                    chebyshev_derivatives[degree] =
                        2.0 * clamped_time * chebyshev_derivatives[degree - 1]
                            + 2.0 * chebyshev_polynomials[degree - 1]
                            - chebyshev_derivatives[degree - 2];
                }
            }

            let velocity_scaling_factor = 2.0 / self.radius;

            velocity[0] = self
                .x
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(c, d)| c * d)
                .sum::<f64>()
                * velocity_scaling_factor;
            velocity[1] = self
                .y
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(c, d)| c * d)
                .sum::<f64>()
                * velocity_scaling_factor;
            velocity[2] = self
                .z
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(c, d)| c * d)
                .sum::<f64>()
                * velocity_scaling_factor;
        }

        (position, velocity)
    }
}

impl fmt::Display for EphemerisRecord {
    /// Pretty‑print a record with midpoint/radius and formatted coefficients.
    ///
    /// Arguments
    /// -----------------
    /// * `f`: Output formatter sink.
    ///
    /// Return
    /// ----------
    /// * A [`fmt::Result`] indicating success or failure.
    ///
    /// See also
    /// ------------
    /// * [`EphemerisRecord::parse`] – Produces the values displayed here.
    /// * `hifitime::Epoch`/`Duration` – Formatting of ET epochs and durations.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mid = Epoch::from_et_seconds(self.mid);
        let radius = Duration::from_seconds(self.radius);

        // Pre-format for consistent column widths
        let mid_str = format!("{mid}");
        let radius_str = format!("{radius}");

        let label_width = 16;
        let value_width = mid_str.len().max(radius_str.len()).max(55);

        let border_header = format!(
            "+{:-<label$}+{:-<value$}+",
            "",
            "",
            label = label_width + 2,
            value = value_width + 2
        );

        writeln!(
            f,
            "+{:^label$}+{:^value$}+",
            "Ephemeris Record",
            "",
            label = label_width + 2,
            value = value_width + 2
        )?;
        writeln!(f, "{border_header}")?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Midpoint",
            mid_str,
            label = label_width,
            value = value_width
        )?;
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Radius",
            radius_str,
            label = label_width,
            value = value_width
        )?;

        writeln!(f, "{border_header}")?;

        // Coefficients section
        writeln!(
            f,
            "| {:<label$} | {:<value$} |",
            "Axis",
            "Chebyshev Coefficients",
            label = label_width,
            value = value_width
        )?;
        writeln!(f, "{border_header}")?;

        for (axis, coeffs) in &[("X", &self.x), ("Y", &self.y), ("Z", &self.z)] {
            writeln!(
                f,
                "| {:<label$} | {:value$} |",
                *axis,
                "",
                label = label_width,
                value = value_width
            )?;
            for chunk in coeffs.chunks(4) {
                let line = chunk
                    .iter()
                    .map(|c| format!("{c:>12.4e}"))
                    .collect::<Vec<_>>()
                    .join(" ");
                writeln!(
                    f,
                    "| {:<label$} | {:<value$} |",
                    "",
                    line,
                    label = label_width,
                    value = value_width
                )?;
            }
            writeln!(f, "{border_header}")?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod test_ephemeris_record {
    use super::*;

    #[test]
    fn test_ephem_record_display() {
        let record = EphemerisRecord {
            mid: -14200056000.0,
            radius: 691200.0,
            x: vec![
                -59117487.054044664,
                -19163216.532728795,
                291991.27938009636,
                15847.329699283478,
                -133.03948110729542,
                -4.459284869049275,
                0.03379900481247174,
                0.0011716375873243507,
                1.4852185006919311e-5,
                -1.1096435596643423e-5,
                -3.3277738887706986e-6,
                -1.7115381406088932e-7,
                1.8759402940064767e-7,
            ],
            y: vec![
                122675629.90130842,
                -7590835.1393347755,
                -614598.6274020968,
                6451.054990886844,
                265.301332213569,
                -1.9960672202990268,
                -0.05594323453310299,
                0.00034547006847617906,
                -7.545272774250726e-5,
                -4.915547872121608e-6,
                2.9685818368805314e-6,
                9.978975536291768e-7,
                4.606523495199457e-8,
            ],
            z: vec![
                53352759.40834735,
                -3301893.184660649,
                -267186.62682516687,
                2806.0081419486182,
                115.33401330537684,
                -0.8678118085367849,
                -0.02423486566672119,
                0.00012129101423135178,
                -4.0581575470237784e-5,
                -1.4813269496375633e-6,
                1.7525667147619721e-6,
                4.996616729595978e-7,
                8.008233021055763e-9,
            ],
        };

        let expected = r#"+ Ephemeris Record +                                                         +
+------------------+---------------------------------------------------------+
| Midpoint         | 1550-01-08T00:00:00 ET                                  |
| Radius           | 8 days                                                  |
+------------------+---------------------------------------------------------+
| Axis             | Chebyshev Coefficients                                  |
+------------------+---------------------------------------------------------+
| X                |                                                         |
|                  |    -5.9117e7    -1.9163e7     2.9199e5     1.5847e4     |
|                  |    -1.3304e2    -4.4593e0    3.3799e-2    1.1716e-3     |
|                  |    1.4852e-5   -1.1096e-5   -3.3278e-6   -1.7115e-7     |
|                  |    1.8759e-7                                            |
+------------------+---------------------------------------------------------+
| Y                |                                                         |
|                  |     1.2268e8    -7.5908e6    -6.1460e5     6.4511e3     |
|                  |     2.6530e2    -1.9961e0   -5.5943e-2    3.4547e-4     |
|                  |   -7.5453e-5   -4.9155e-6    2.9686e-6    9.9790e-7     |
|                  |    4.6065e-8                                            |
+------------------+---------------------------------------------------------+
| Z                |                                                         |
|                  |     5.3353e7    -3.3019e6    -2.6719e5     2.8060e3     |
|                  |     1.1533e2   -8.6781e-1   -2.4235e-2    1.2129e-4     |
|                  |   -4.0582e-5   -1.4813e-6    1.7526e-6    4.9966e-7     |
|                  |    8.0082e-9                                            |
+------------------+---------------------------------------------------------+
"#;
        let output = format!("{record}");
        assert_eq!(output, expected);
    }

    #[test]
    fn test_record_interpolation() {
        let record = EphemerisRecord {
            mid: -14200056000.0,
            radius: 691200.0,
            x: vec![
                -59117487.054044664,
                -19163216.532728795,
                291991.27938009636,
                15847.329699283478,
                -133.03948110729542,
                -4.459284869049275,
                0.03379900481247174,
                0.0011716375873243507,
                1.4852185006919311e-5,
                -1.1096435596643423e-5,
                -3.3277738887706986e-6,
                -1.7115381406088932e-7,
                1.8759402940064767e-7,
            ],
            y: vec![
                122675629.90130842,
                -7590835.1393347755,
                -614598.6274020968,
                6451.054990886844,
                265.301332213569,
                -1.9960672202990268,
                -0.05594323453310299,
                0.00034547006847617906,
                -7.545272774250726e-5,
                -4.915547872121608e-6,
                2.9685818368805314e-6,
                9.978975536291768e-7,
                4.606523495199457e-8,
            ],
            z: vec![
                53352759.40834735,
                -3301893.184660649,
                -267186.62682516687,
                2806.0081419486182,
                115.33401330537684,
                -0.8678118085367849,
                -0.02423486566672119,
                0.00012129101423135178,
                -4.0581575470237784e-5,
                -1.4813269496375633e-6,
                1.7525667147619721e-6,
                4.996616729595978e-7,
                8.008233021055763e-9,
            ],
        };

        let date_str = "2024-04-10T12:30:45 UTC";
        let epoch = Epoch::from_gregorian_str(date_str)
            .unwrap()
            .to_time_scale(hifitime::TimeScale::UTC);

        let (position, velocity) = record.interpolate(epoch.to_et_seconds());

        assert_eq!(
            position,
            Vector3::new(-77973002.44148895, 114476910.4391533, 49786600.04705159)
        );

        assert_eq!(
            velocity,
            Vector3::new(
                -51.663380462484454,
                -28.897503577030022,
                -12.568179510296154
            )
        )
    }
}
