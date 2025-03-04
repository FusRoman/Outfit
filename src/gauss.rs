use aberth::StopReason;
use core::fmt;
use nalgebra::Matrix3;
use nalgebra::Vector3;

use super::constants::GaussGrav;

use aberth::aberth;

use super::env_state::OutfitState;
use super::jpl_request::observer_pos::helio_obs_pos;
use super::kepler::velocity_correction;
use super::keplerian_orbit::KeplerianOrbit;
use super::orb_elem::ccek1;
use super::orb_elem::eccentricity_control;
use super::ref_system::rotpn;

/// Gauss struct data
/// ra is right ascension in degree
/// dec is declination in degree
/// time is the observation time in modified julian date (MJD
/// observer_position is the position of the observer from where the observation have been taken.
/// The observer position is in equatorial mean J2000 reference frame and units is astronomical unit
struct GaussObs {
    ra: Vector3<f64>,
    dec: Vector3<f64>,
    time: Vector3<f64>,
    observer_position: Matrix3<f64>,
}

/// Define the errors that the Gauss method could return during the execution
///
/// GaussSingMatrix means that the matrix made of unit_vector towards the orbiting body cannot be inverse.
#[derive(Debug, Clone)]
struct GaussSingMatrix;

/// Solve8PolyFailed is used in case the Aberth–Ehrlich method failed to return roots.
#[derive(Debug, Clone)]
struct Solve8PolyFailed;

/// Spurious root, root not accepted for orbital estimation
#[derive(Debug, Clone)]
struct SpuriousRoot;

impl fmt::Display for GaussSingMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The unit matrix cannot be inverse, orbit could be coplanar"
        )
    }
}

impl fmt::Display for Solve8PolyFailed {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The Aberth–Ehrlich method failed to find complex roots for the 8-order polynom."
        )
    }
}

impl fmt::Display for SpuriousRoot {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Spurious root detected for the orbital estimation")
    }
}

impl GaussObs {
    /// Initialise the struct used for the Gauss method.
    /// Use only three observations to estimate an initial orbit
    pub async fn new(
        ra: Vector3<f64>,
        dec: Vector3<f64>,
        mjd_time: Vector3<f64>,
        longitude: f64,
        latitude: f64,
        height: f64,
    ) -> GaussObs {
        let state = OutfitState::new().await;

        GaussObs {
            ra: ra,
            dec: dec,
            time: mjd_time,
            observer_position: helio_obs_pos(&mjd_time, longitude, latitude, height, &state).await,
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
            GaussObs::unit_vector(self.ra[0], self.dec[0], f64::cos(self.dec[0])),
            GaussObs::unit_vector(self.ra[1], self.dec[1], f64::cos(self.dec[1])),
            GaussObs::unit_vector(self.ra[2], self.dec[2], f64::cos(self.dec[2])),
        ])
    }

    fn gauss_prelim(
        &self,
    ) -> Result<
        (
            f64,
            f64,
            Matrix3<f64>,
            Matrix3<f64>,
            Vector3<f64>,
            Vector3<f64>,
        ),
        GaussSingMatrix,
    > {
        let tau1 = GaussGrav * (self.time[0] - self.time[1]);
        let tau3 = GaussGrav * (self.time[2] - self.time[1]);
        let tau13 = tau3 - tau1;
        let vector_a = Vector3::new(tau3 / tau13, -1.0, -(tau1 / tau13));
        let vector_b = Vector3::new(
            vector_a[0] * (tau13.powi(2) - tau3.powi(2)) / 6.0,
            0.0,
            vector_a[2] * (tau13.powi(2) - tau1.powi(2)) / 6.0,
        );

        let unit_matrix = self.unit_matrix();
        let Some(inv_unit_matrix) = unit_matrix.try_inverse() else {
            return Err(GaussSingMatrix);
        };

        Ok((tau1, tau3, unit_matrix, inv_unit_matrix, vector_a, vector_b))
    }

    fn coeff_eight_poly(
        &self,
        unit_matrix: &Matrix3<f64>,
        unit_invmatrix: &Matrix3<f64>,
        vector_a: &Vector3<f64>,
        vector_b: &Vector3<f64>,
    ) -> (f64, f64, f64) {
        let observer_position_t = self.observer_position.transpose();
        let ra = observer_position_t * vector_a;
        let rb = observer_position_t * vector_b;

        let second_row_t = unit_invmatrix.row(1).transpose();
        let a2star = second_row_t.dot(&ra);
        let b2star = second_row_t.dot(&rb);

        let r22 = self
            .observer_position
            .row(1)
            .component_mul(&self.observer_position.row(1))
            .sum();
        let s2r2 = unit_matrix
            .column(1)
            .transpose()
            .dot(&self.observer_position.row(1));

        (
            -(a2star.powi(2)) - r22 - (2.0 * a2star * s2r2),
            -(2.0 * b2star * (a2star + s2r2)),
            -(b2star.powi(2)),
        )
    }

    fn solve_8poly(
        &self,
        polynom: &[f64; 9],
        max_iterations: u32,
        aberth_epsilon: f64,
        root_acceptance_epsilon: f64,
    ) -> Result<Vec<f64>, Solve8PolyFailed> {
        let roots = aberth(&polynom, max_iterations, aberth_epsilon);
        match roots.stop_reason {
            StopReason::Converged(_) | StopReason::MaxIteration(_) => {
                return Ok(roots
                    .iter()
                    // root criteria acceptance: root must be real and real part must be positive
                    .filter(|complex| complex.re > 0. && complex.im.abs() < root_acceptance_epsilon)
                    .map(|complex| complex.re)
                    .collect::<Vec<f64>>());
            }
            StopReason::Failed(_) => {
                return Err(Solve8PolyFailed);
            }
        }
    }

    /// Compute the asteroid position vector at the time of the observation
    ///
    /// Argument
    /// --------
    /// * roots: the roots computed from solve_8poly
    /// * unit_matrix: the matrix made of unit_vector to the orbiting body
    /// * unit_matinv: the inverse of the unit matrix
    /// * vector_a and vector_b returned by the gauss_prelim function
    ///
    /// Return
    /// ------
    /// * the position vector of the asteroid at the time of the observation
    fn asteroid_position_vector(
        &self,
        unit_matrix: &Matrix3<f64>,
        unit_matinv: &Matrix3<f64>,
        vector_c: &Vector3<f64>,
    ) -> Result<Matrix3<f64>, SpuriousRoot> {
        let obs_pos_t = self.observer_position.transpose();
        let gcap = obs_pos_t * vector_c;
        let crhom = unit_matinv * gcap;
        let rho: Vector3<f64> = -crhom.component_div(&vector_c);
        if rho[1] < 0.01 {
            return Err(SpuriousRoot);
        }
        let rho_unit = Matrix3::from_columns(&[
            rho[0] * unit_matrix.column(0),
            rho[1] * unit_matrix.column(1),
            rho[2] * unit_matrix.column(2),
        ]);

        Ok(obs_pos_t + rho_unit)
    }

    /// Compute the velocity vector of the asteroid at the time of the second observation
    ///
    /// Inputs:
    ///     ast_pos_vector: Asteroid position vector at the time of the three observations
    ///     cartesian representation, unit= AU
    ///       x y z
    ///    t1 0 1 2
    ///    t2 4 5 6
    ///    t3 7 8 9
    ///    
    ///    tau1: Normalized time of the first observation
    ///    tau3: normalized time of the third observation
    ///
    /// Output:
    ///     Velocity vector of the asteroid [vx, vy, vz]
    fn gibbs_correction(
        &self,
        ast_pos_vector: &Matrix3<f64>,
        tau1: f64,
        tau3: f64,
    ) -> Vector3<f64> {
        let tau13 = tau3 - tau1;
        let r1m3 = 1. / ast_pos_vector.column(0).norm().powi(3);
        let r2m3 = 1. / ast_pos_vector.column(1).norm().powi(3);
        let r3m3 = 1. / ast_pos_vector.column(2).norm().powi(3);

        let d1 = tau3 * (r1m3 / 12. - 1. / (tau1 * tau13));
        let d2: f64 = (tau1 + tau3) * (r2m3 / 12. - 1. / (tau1 * tau3));
        let d3 = -tau1 * (r3m3 / 12. + 1. / (tau3 * tau13));
        let d_vect = Vector3::new(-d1, d2, d3);

        Vector3::new(
            GaussGrav * ast_pos_vector.row(0).dot(&d_vect.transpose()),
            GaussGrav * ast_pos_vector.row(1).dot(&d_vect.transpose()),
            GaussGrav * ast_pos_vector.row(2).dot(&d_vect.transpose()),
        )
    }

    fn accept_root(
        &self,
        root: f64,
        unit_matrix: &Matrix3<f64>,
        inv_unit_matrix: &Matrix3<f64>,
        vect_a: &Vector3<f64>,
        vect_b: &Vector3<f64>,
        tau1: f64,
        tau3: f64,
    ) -> Option<(Matrix3<f64>, Vector3<f64>)> {
        let r2m3 = 1. / root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vect_a[0] + vect_b[0] * r2m3,
            -1.,
            vect_a[2] + vect_b[2] * r2m3,
        );
        let Some(ast_pos_all_time) = self
            .asteroid_position_vector(&unit_matrix, &inv_unit_matrix, &vector_c)
            .ok()
        else {
            return None;
        };
        let ast_pos_second_time: Vector3<f64> = ast_pos_all_time.column(1).into();
        let asteroid_vel = self.gibbs_correction(&ast_pos_all_time, tau1, tau3);
        let Some((is_accepted, _, _, _)) =
            eccentricity_control(&ast_pos_second_time, &asteroid_vel, 1e3, 5.)
        else {
            return None;
        };
        if is_accepted {
            return Some((ast_pos_all_time, asteroid_vel));
        }
        return None;
    }

    fn from_position_velocity_to_orbit(
        &self,
        &asteroid_position: &Vector3<f64>,
        &asteroid_velocity: &Vector3<f64>,
    ) -> KeplerianOrbit {
        let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
        rotpn(&mut roteqec, "EQUM", "J2000", 0., "ECLM", "J2000", 0.);

        let matrix_elc_transform = Matrix3::from(roteqec).transpose();
        let ecl_pos = matrix_elc_transform * asteroid_position;
        let ecl_vel = matrix_elc_transform * asteroid_velocity;

        let ast_pos_vel: [f64; 6] = [
            ecl_pos.x, ecl_pos.y, ecl_pos.z, ecl_vel.x, ecl_vel.y, ecl_vel.z,
        ];

        let mut elem = [0.0; 6];
        let mut type_ = String::new();
        ccek1(&mut elem, &mut type_, &ast_pos_vel);
        KeplerianOrbit {
            semi_major_axis: elem[0],
            eccentricity: elem[1],
            inclination: elem[2],
            ascending_node_longitude: elem[3],
            periapsis_argument: elem[4],
            mean_anomaly: elem[5],
        }
    }

    // TODO: remplacer le return de Option<KeplerianOrbit> en enum avec 3 champs
    // un champ PrelimOrbit(KeplerienOrbit)
    // un champ CorrectedOrbit(KeplerianOrbit)
    // un champ FailOrbit
    pub fn prelim_orbit(&self) -> Option<KeplerianOrbit> {
        let (tau1, tau3, unit_matrix, inv_unit_matrix, vect_a, vect_b) =
            self.gauss_prelim().unwrap();

        let (coeff_6, coeff_3, coeff_0) =
            self.coeff_eight_poly(&unit_matrix, &inv_unit_matrix, &vect_a, &vect_b);
        let polynomial = [coeff_0, 0., 0., coeff_3, 0., 0., coeff_6, 0., 1.];

        let roots = self.solve_8poly(&polynomial, 50, 1e-6, 1e-6).unwrap();

        // get the first accepted root and return the asteroid position and asteroid velocity
        let Some((asteroid_pos_all_time, asteroid_vel)) = roots.into_iter().find_map(|root| {
            self.accept_root(
                root,
                &unit_matrix,
                &inv_unit_matrix,
                &vect_a,
                &vect_b,
                tau1,
                tau3,
            )
        }) else {
            return None;
        };

        let Some((corrected_pos, corrected_vel)) = self.pos_and_vel_correction(
            &asteroid_pos_all_time,
            &asteroid_vel,
            &unit_matrix,
            &inv_unit_matrix,
            1e3,
            5.,
            1e-10,
            50,
        ) else {
            return Some(self.from_position_velocity_to_orbit(
                &asteroid_pos_all_time.column(1).into(),
                &asteroid_vel,
            ));
        };

        Some(self.from_position_velocity_to_orbit(&corrected_pos.column(1).into(), &corrected_vel))
    }

    fn pos_and_vel_correction(
        &self,
        asteroid_position: &Matrix3<f64>,
        asteroid_velocity: &Vector3<f64>,
        unit_matrix: &Matrix3<f64>,
        inv_unit_matrix: &Matrix3<f64>,
        peri_max: f64,
        ecc_max: f64,
        err_max: f64,
        itmax: i32,
    ) -> Option<(Matrix3<f64>, Vector3<f64>)> {
        let mut previous_ast_pos = asteroid_position.clone();
        let mut previous_ast_vel = asteroid_velocity.clone();
        for _ in 0..itmax {
            let Ok((ast_vel1, f1, g1)) = velocity_correction(
                &previous_ast_pos.column(0).into(),
                &previous_ast_pos.column(1).into(),
                &previous_ast_vel,
                self.time[0] - self.time[1],
                peri_max,
                ecc_max,
            ) else {
                continue;
            };

            let Ok((ast_vel2, f2, g2)) = velocity_correction(
                &previous_ast_pos.column(2).into(),
                &previous_ast_pos.column(1).into(),
                &previous_ast_vel,
                self.time[2] - self.time[1],
                peri_max,
                ecc_max,
            ) else {
                continue;
            };

            let corrected_velocity = (ast_vel1 + ast_vel2) / 2.;
            let f_lagrange = f1 * g2 - f2 * g1;

            let c_vector = Vector3::new(g2 / f_lagrange, -1., -(g1 / f_lagrange));
            let Some(new_ast_pos) = self
                .asteroid_position_vector(&unit_matrix, &inv_unit_matrix, &c_vector)
                .ok()
            else {
                continue;
            };

            let diff_pos = new_ast_pos - previous_ast_pos;
            let diff_pos_squared = diff_pos.component_mul(&diff_pos).sum();
            let sca = new_ast_pos.component_mul(&new_ast_pos).sum();
            let ast_pos_err = (diff_pos_squared / sca).sqrt();

            let Some((is_accepted, _, _, _)) =
                eccentricity_control(&new_ast_pos.column(1).into(), &corrected_velocity, 1e3, 5.)
            else {
                return None;
            };

            if !is_accepted {
                return None;
            }

            previous_ast_pos = new_ast_pos;
            previous_ast_vel = corrected_velocity;

            if ast_pos_err <= err_max {
                break;
            }
        }
        Some((previous_ast_pos, previous_ast_vel))
    }
}

#[cfg(test)]
mod gauss_test {

    use super::*;

    #[test]
    fn test_gauss_prelim() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::zeros(),
        };

        let (tau1, tau3, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        assert_eq!(tau1, -0.35721620648079105);
        assert_eq!(tau3, 0.25342080566844405);

        let unit_mat_array: [f64; 9] = unit_matrix
            .as_slice()
            .try_into()
            .expect("Conversion failed");
        assert_eq!(
            unit_mat_array,
            [
                -0.05549934652247514,
                0.46585594034226024,
                0.8831183756345503,
                -0.06972979004485365,
                0.5827357012279389,
                0.8096646582966821,
                -0.12245931009139571,
                0.6656387438390606,
                0.7361581216507068
            ]
        );

        let inv_unit_mat_array: [f64; 9] = inv_unit_matrix
            .as_slice()
            .try_into()
            .expect("Conversion failed");
        assert_eq!(
            inv_unit_mat_array,
            [
                -18.774792915974594,
                41.814279122702025,
                -23.466669573973437,
                -8.16479071034311,
                11.489343729350427,
                -2.8418335594428186,
                4.259482782736117,
                -3.432964304649723,
                0.024345794753282718
            ]
        );

        let vect_a_array: [f64; 3] = vector_a.as_slice().try_into().expect("Conversion failed");
        assert_eq!(
            vect_a_array,
            [0.41501055557783634, -1.0, 0.5849894444221637]
        );

        let vect_b_array: [f64; 3] = vector_b.as_slice().try_into().expect("Conversion failed");
        assert_eq!(
            vect_b_array,
            [0.021349212036493866, 0.0, 0.023913797385599792]
        );
    }

    #[test]
    fn test_coeff_8poly() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::new(
                -0.26456661713915464,
                0.86893516436949503,
                0.37669962110919220,
                -0.58916318521741273,
                0.72388725167947765,
                0.31381865165245848,
                -0.77438744379695956,
                0.56128847092611645,
                0.24334971075289916,
            ),
        };

        let (tau1, tau3, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        let (coeff_6, coeff_3, coeff_0) =
            gauss.coeff_eight_poly(&unit_matrix, &inv_unit_matrix, &vector_a, &vector_b);

        assert_eq!(coeff_6, -2.615803718759013);
        assert_eq!(coeff_3, 2.0305173353541064);
        assert_eq!(coeff_0, -0.4771346939201045);
    }

    #[test]
    fn test_solving_polynom() {
        let gauss = GaussObs {
            ra: Vector3::zeros(),
            dec: Vector3::zeros(),
            time: Vector3::zeros(),
            observer_position: Matrix3::zeros(),
        };

        let polynomial = [
            -0.47713469392010482,
            0.,
            0.,
            2.0305173353541064,
            0.,
            0.,
            -2.6158037187590111,
            0.,
            1.,
        ];

        let roots = gauss.solve_8poly(&polynomial, 50, 1e-6, 1e-6).unwrap();

        assert_eq!(
            roots,
            vec![1.3856312487504954, 0.7328107254669438, 0.9540135094917113]
        );
    }

    #[test]
    fn test_asteroid_position() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::new(
                -0.26456661713915464,
                0.86893516436949503,
                0.37669962110919220,
                -0.58916318521741273,
                0.72388725167947765,
                0.31381865165245848,
                -0.77438744379695956,
                0.56128847092611645,
                0.24334971075289916,
            ),
        };

        let (_, _, unit_matrix, inv_unit_matrix, vector_a, vector_b) =
            gauss.gauss_prelim().unwrap();

        let first_root: f64 = 0.73281072546694370;
        let r2m3 = 1. / first_root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vector_a[0] + vector_b[0] * r2m3,
            -1.,
            vector_a[2] + vector_b[2] * r2m3,
        );
        let ast_pos = gauss.asteroid_position_vector(&unit_matrix, &inv_unit_matrix, &vector_c);
        assert!(ast_pos.is_err());

        let second_root: f64 = 1.3856312487504951;
        let r2m3 = 1. / second_root.powi(3);
        let vector_c: Vector3<f64> = Vector3::new(
            vector_a[0] + vector_b[0] * r2m3,
            -1.,
            vector_a[2] + vector_b[2] * r2m3,
        );
        let ast_pos = gauss
            .asteroid_position_vector(&unit_matrix, &inv_unit_matrix, &vector_c)
            .unwrap();

        let ast_pos_slice: [f64; 9] = ast_pos
            .as_slice()
            .try_into()
            .expect("test_asteroid_position result matrix into slice failed");

        assert_eq!(
            ast_pos_slice,
            [
                -0.28811969067349597,
                1.06663729794052,
                0.7514815481797275,
                -0.6235500510031637,
                1.0112601855976917,
                0.713100363506241,
                -0.8445850475187664,
                0.9428539454255418,
                0.6653391541170498
            ]
        )
    }

    #[test]
    fn test_gibbs_correction() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::zeros(),
        };
        let (tau1, tau3, _, _, _, _) = gauss.gauss_prelim().unwrap();

        let ast_pos = Matrix3::new(
            -0.28811969067349597,
            1.06663729794052,
            0.7514815481797275,
            -0.6235500510031637,
            1.0112601855976917,
            0.713100363506241,
            -0.8445850475187664,
            0.9428539454255418,
            0.6653391541170498,
        )
        .transpose();
        let asteroid_velocity = gauss.gibbs_correction(&ast_pos, tau1, tau3);

        let ast_vel_slice = asteroid_velocity.as_slice();
        assert_eq!(
            ast_vel_slice,
            [
                -0.015549845137774663,
                -0.003876936109837664,
                -0.0027014074002979886
            ]
        )
    }

    #[test]
    fn test_solve_orbit() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::new(
                -0.26456661713915464,
                0.86893516436949503,
                0.37669962110919220,
                -0.58916318521741273,
                0.72388725167947765,
                0.31381865165245848,
                -0.77438744379695956,
                0.56128847092611645,
                0.24334971075289916,
            ),
        };

        let prelim_orbit = gauss.prelim_orbit().unwrap();
        assert_eq!(
            prelim_orbit,
            KeplerianOrbit {
                semi_major_axis: 1.8015109749705496,
                eccentricity: 0.2835090890769234,
                inclination: 0.20264596920729464,
                ascending_node_longitude: 0.008079724163250196,
                periapsis_argument: 1.244926174041886,
                mean_anomaly: 0.44069722652585874
            }
        )
    }

    #[test]
    fn test_orbit_correction() {
        let gauss = GaussObs {
            ra: Vector3::new(1.6893715963476696, 1.6898894500811472, 1.7527345385664372),
            dec: Vector3::new(1.0824680373855251, 0.94358050479462163, 0.82737624078999861),
            time: Vector3::new(57028.479297592596, 57049.245147592592, 57063.977117592593),
            observer_position: Matrix3::new(
                -0.26456661713915464,
                0.86893516436949503,
                0.37669962110919220,
                -0.58916318521741273,
                0.72388725167947765,
                0.31381865165245848,
                -0.77438744379695956,
                0.56128847092611645,
                0.24334971075289916,
            ),
        };

        let ast_pos = Matrix3::new(
            -0.28811969067349602,
            1.0666372979405205,
            0.75148154817972834,
            -0.62355005100316385,
            1.0112601855976919,
            0.71310036350624140,
            -0.84458504751876651,
            0.94285394542554213,
            0.66533915411705014,
        )
        .transpose();

        let ast_vel = Vector3::new(
            -1.5549845137774663E-002,
            -3.8769361098376577E-003,
            -2.7014074002979964E-003,
        );

        let unit_matrix = Matrix3::new(
            -5.5499346522475138E-002,
            0.46585594034226024,
            0.88311837563455031,
            -6.9729790044853648E-002,
            0.58273570122793894,
            0.80966465829668210,
            -0.12245931009139571,
            0.66563874383906063,
            0.73615812165070682,
        )
        .transpose();

        let inv_unit_matrix = Matrix3::new(
            -18.774792915974604,
            41.814279122702033,
            -23.466669573973437,
            -8.1647907103431088,
            11.489343729350425,
            -2.8418335594428177,
            4.2594827827361170,
            -3.4329643046497211,
            2.4345794753281792E-002,
        )
        .transpose();

        let (new_ast_pos, new_ast_vel) = gauss
            .pos_and_vel_correction(
                &ast_pos,
                &ast_vel,
                &unit_matrix,
                &inv_unit_matrix,
                1e3,
                1.,
                1e-10,
                50,
            )
            .unwrap();

        assert_eq!(
            new_ast_pos.as_slice(),
            [
                -0.2878540141559046,
                1.06440723593647,
                0.7472540422835181,
                -0.6231216182863288,
                1.0076797497536536,
                0.7081256342111117,
                -0.8435611164802848,
                0.9372882749205874,
                0.6591838430228918
            ]
        );

        assert_eq!(
            new_ast_vel.as_slice(),
            [
                -0.015524309979972159,
                -0.003984105628190921,
                -0.0027640157742952693
            ]
        );
    }
}
