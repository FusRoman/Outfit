use std::{
    fmt,
    fs::File,
    io::{BufReader, Read, Seek},
};

use hifitime::{Duration, Epoch};
use nalgebra::Vector3;

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
            // optimisation avec chunk_exacts(8)
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

    /// Interpolate the position and velocity at a given ephemeris time.
    ///
    /// Arguments
    /// ---------
    /// * `ephem_time` - Ephemeris time at which compute the position and the velocity. The reference time frame is BDT (Barycentric Dynamical Time).
    ///
    /// Returns
    /// -------
    /// * A tuple containing the position and velocity vectors.
    ///
    /// The position and velocity are computed using Chebyshev polynomials.
    /// The Chebyshev coefficients are stored in the `x`, `y`, and `z` fields of the struct.
    /// The interpolation is done in the range [-1, 1].
    /// where `T_n(t)` is the Chebyshev polynomial of degree `n` and `t` is the normalized time.
    /// The velocity is computed using the derivative of the Chebyshev polynomial
    /// The position is in kilometers and the velocity is in kilometers per second.
    /// The Chebyshev coefficients are assumed to be in the order of decreasing degree.
    /// The Chebyshev polynomials are evaluated using Horner's method for efficiency.
    pub fn interpolate(&self, ephem_time: f64) -> (Vector3<f64>, Vector3<f64>) {
        let normalized_time = (ephem_time - self.mid) / self.radius;
        let clamped_time = normalized_time.clamp(-1.0, 1.0);

        let num_coefficients = self.x.len();
        let mut chebyshev_polynomials = vec![0.0; num_coefficients];
        chebyshev_polynomials[0] = 1.0;

        // Chebyshev polynomials are defined for n >= 0
        // T_0(t) = 1
        if num_coefficients > 1 {
            chebyshev_polynomials[1] = clamped_time;
            for degree in 2..num_coefficients {
                // T_n(t) = 2 * t * T_(n-1)(t) - T_(n-2)(t)
                // Horner's method for Chebyshev polynomials
                chebyshev_polynomials[degree] =
                    2.0 * clamped_time * chebyshev_polynomials[degree - 1]
                        - chebyshev_polynomials[degree - 2];
            }
        }

        // Compute the position using the Chebyshev coefficients
        // The position is a linear combination of the Chebyshev polynomials
        // The coefficients are the Chebyshev coefficients for each axis
        let position = Vector3::new(
            self.x
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(coefficient, polynomial)| coefficient * polynomial)
                .sum(),
            self.y
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(coefficient, polynomial)| coefficient * polynomial)
                .sum(),
            self.z
                .iter()
                .zip(&chebyshev_polynomials)
                .map(|(coefficient, polynomial)| coefficient * polynomial)
                .sum(),
        );

        // Compute the velocity using the Chebyshev derivatives
        // The velocity is a linear combination of the Chebyshev derivatives
        let mut velocity = Vector3::zeros();
        if num_coefficients > 1 {
            let mut chebyshev_derivatives = vec![0.0; num_coefficients];
            chebyshev_derivatives[1] = 1.0;

            if num_coefficients > 2 {
                chebyshev_derivatives[2] = 4.0 * clamped_time;
                for degree in 3..num_coefficients {
                    // T'_n(t) = 2 * T'_(n-1)(t) + 2 * t * T_(n-1)(t) - T'_(n-2)(t)
                    // Horner's method for Chebyshev derivatives
                    chebyshev_derivatives[degree] =
                        2.0 * clamped_time * chebyshev_derivatives[degree - 1]
                            + 2.0 * chebyshev_polynomials[degree - 1]
                            - chebyshev_derivatives[degree - 2];
                }
            }

            // The velocity is scaled by a factor of 2 / (b - a) = 2 / radius
            // where b and a are the bounds of the interval [-1, 1]
            // The radius is the length of the interval
            let velocity_scaling_factor = 2.0 / self.radius;

            velocity[0] = self
                .x
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(coefficient, derivative)| coefficient * derivative)
                .sum::<f64>()
                * velocity_scaling_factor;
            velocity[1] = self
                .y
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(coefficient, derivative)| coefficient * derivative)
                .sum::<f64>()
                * velocity_scaling_factor;
            velocity[2] = self
                .z
                .iter()
                .zip(&chebyshev_derivatives)
                .skip(1)
                .map(|(coefficient, derivative)| coefficient * derivative)
                .sum::<f64>()
                * velocity_scaling_factor;
        }

        (position, velocity)
    }
}

impl fmt::Display for EphemerisRecord {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mid = Epoch::from_et_seconds(self.mid);
        let radius = Duration::from_seconds(self.radius);

        // Conversion des champs complexes en String pour mesurer leur largeur
        let mid_str = format!("{mid}");
        let radius_str = format!("{radius}");

        // Détermination de la largeur max pour les valeurs à gauche et droite
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

        // Affichage des coefficients Chebyshev
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
