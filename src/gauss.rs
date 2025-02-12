use core::fmt;
use std::vec;

use nalgebra::Matrix3;
use nalgebra::Vector3;

use super::jpl_request::pos_vector::{date_to_jd, get_helio_pos};

/// Gaussian gravitational constant
const GaussGrav: f64 = 0.0172020989484;

struct OrbitGauss {
    ra: Vector3<f64>,
    dec: Vector3<f64>,
    time: Vector3<f64>,
    sun_pos: Matrix3<f64>,
}

#[derive(Debug, Clone)]
struct GaussSingMatrix;

impl fmt::Display for GaussSingMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The unit matrix cannot be inverse, orbot could be coplanar"
        )
    }
}

impl OrbitGauss {
    /// Initialise the struct used for the Gauss method.
    /// Use only three observations to estimate an initial orbit
    pub fn new(
        ra: Vector3<f64>,
        dec: Vector3<f64>,
        time: Vector3<f64>,
        sun_pos: Matrix3<f64>,
    ) -> OrbitGauss {
        OrbitGauss {
            ra: ra,
            dec: dec,
            time: time,
            sun_pos: sun_pos,
        }
    }

    /// Compute the orbiting body direction cosine vector
    /// (from observations point to orbiting body in Topocentric equatorial coordinate system)
    fn unit_vector(ra: f64, dec: f64, cos_dec: f64) -> Vector3<f64> {
        Vector3::new(
            f64::cos(ra) * cos_dec,
            f64::sin(ra) * cos_dec,
            f64::sin(dec),
        )
    }

    /// Compute the matrix containing the unit vector
    fn unit_matrix(&self) -> Matrix3<f64> {
        Matrix3::from_columns(&[
            OrbitGauss::unit_vector(self.ra[0], self.dec[0], f64::cos(self.dec[0])),
            OrbitGauss::unit_vector(self.ra[1], self.dec[1], f64::cos(self.dec[1])),
            OrbitGauss::unit_vector(self.ra[2], self.dec[2], f64::cos(self.dec[2])),
        ])
    }

    fn compute_A_B_vector(&self) -> (Vector3<f64>, Vector3<f64>) {
        let tau1 = GaussGrav * (self.time[0] - self.time[1]);
        let tau3 = GaussGrav * (self.time[2] - self.time[1]);
        let tau13 = tau3 - tau1;
        let vector_a = Vector3::new(tau3 / tau13, -1.0, tau1 / tau13);
        let vector_b = Vector3::new(
            vector_a[0] * (tau13.powi(2) - tau3.powi(2) / 6.0),
            0.0,
            vector_a[2] * (tau13.powi(2) - tau1.powi(2) / 6.0),
        );
        (vector_a, vector_b)
    }

    fn coeff_eight_poly(&self) -> Result<(f64, f64, f64), GaussSingMatrix> {
        let unit_matrix = self.unit_matrix();
        if let Some(inv_unit_matrix) = unit_matrix.try_inverse() {
            let (a_vector, b_vector) = self.compute_A_B_vector();

            let ra = self.sun_pos * a_vector;
            let rb = self.sun_pos * b_vector;

            let a2star = inv_unit_matrix.column(1).dot(&ra);
            let b2star = inv_unit_matrix.column(1).dot(&rb);

            let r22 = self.sun_pos.row(1).norm();
            let s2r2 = unit_matrix.row(1).dot(&self.sun_pos.row(1));
            Ok((
                -(a2star.powi(2) - r22 - (2.0 * a2star * s2r2)),
                -(2.0 * b2star * (a2star + s2r2)),
                -(b2star.powi(2)),
            ))
        } else {
            Err(GaussSingMatrix)
        }
    }
}

#[cfg(test)]
mod gauss_test {
    use super::*;

    #[tokio::test]
    async fn test_unit_matrix() {
        let date_list = vec![
            "2021-07-04T12:47:24",
            "2021-07-04T13:47:24",
            "2021-07-04T14:47:24",
        ];
        let jd_time = date_to_jd(&date_list);
        let pos_vector = get_helio_pos(&jd_time).await;

        let sun_pos_matrix = Matrix3::from_columns(&[
            pos_vector.get(0).unwrap().pos_vector(),
            pos_vector.get(1).unwrap().pos_vector(),
            pos_vector.get(2).unwrap().pos_vector(),
        ]);

        let gauss = OrbitGauss {
            ra: Vector3::new(0.0, 1.0, 2.0),
            dec: Vector3::new(0.0, 1.0, 2.0),
            time: Vector3::from_vec(jd_time),
            sun_pos: sun_pos_matrix,
        };

        let (coeff_8, coeff_6, coeff_3) = gauss.coeff_eight_poly().expect("");
        println!("{coeff_8}");
        println!("{coeff_6}");
        println!("{coeff_3}");
    }
}
