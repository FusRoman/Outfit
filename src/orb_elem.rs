use crate::orbit_type::cometary_element::CometaryElements;
use crate::orbit_type::keplerian_element::KeplerianElements;
use crate::orbit_type::OrbitalElements;

use super::constants::{DPI, GAUSS_GRAV_SQUARED};
use super::ref_system::rotmt;
use nalgebra::Vector3;

fn atan2(y: f64, x: f64) -> f64 {
    y.atan2(x)
}

fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

fn prodmm(a: &mut [[f64; 3]; 3], b: [[f64; 3]; 3], c: [[f64; 3]; 3]) {
    let w = b;
    let z = c;

    for j in 0..3 {
        for k in 0..3 {
            let mut s = 0.0;
            for (l, row) in w[j].iter().enumerate() {
                s += row * z[l][k];
            }
            a[j][k] = s;
        }
    }
}

fn prodmv(y: &mut [f64; 3], a: [[f64; 3]; 3], x: [f64; 3]) {
    let z = x;

    for (j, row) in a.iter().enumerate() {
        let mut s = 0.0;
        for (l, val) in row.iter().enumerate() {
            s += val * z[l];
        }
        y[j] = s;
    }
}

/// Converts a Cartesian state vector into orbital elements in the J2000 equatorial frame.
///
/// This function takes a 6D state vector representing position and velocity in an inertial
/// heliocentric reference frame and converts it into classical orbital elements (COEs).
///
/// The computed orbital elements depend on the type of orbit:
/// - For elliptic orbits (`1/a > 0`), the output is:  
///   `[a, e, i, Ω, ω, M]`  
///   where:
///   - `a` is the semi-major axis (AU)
///   - `e` is the eccentricity
///   - `i` is the inclination (radians)
///   - `Ω` is the longitude of ascending node (radians)
///   - `ω` is the argument of periapsis (radians)
///   - `M` is the mean anomaly (radians)
///
/// - For parabolic and hyperbolic orbits (`1/a ≤ 0`), the output is:
///   `[q, e, i, Ω, ω, ν]`  
///   where:
///   - `q` is the perihelion distance (AU)
///   - `e` is the eccentricity (=1 for parabolic)
///   - `ν` is the true anomaly (radians)
///
/// Arguments
/// ---------
/// * `xv` – A 6-element array `[x, y, z, vx, vy, vz]` in heliocentric J2000 coordinates.
///   Positions are in astronomical units (AU) and velocities in AU/day.
/// * `reference_epoch` – Epoch of the state vector in Julian Date (TDB).
///
/// Returns
/// -------
/// An `OrbitalElements` enum containing either `KeplerianElements` or `CometaryElements`.
///
/// Note
/// ----
/// This function does not account for planetary perturbations or relativistic corrections.
/// The orbital elements are computed under the assumption of a two-body Keplerian motion.
///
/// See also
/// --------
/// * [`eccentricity_control`] – Filters orbital solutions based on eccentricity and perihelion distance thresholds.
/// * [`KeplerianElements`] – Struct holding canonical orbital elements.
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Lagrange-based velocity refinement from triplet position vectors.
/// * [`GAUSS_GRAV_SQUARED`] – Gaussian gravitational constant.
pub fn ccek1(xv: &[f64; 6], reference_epoch: f64) -> OrbitalElements {
    let mut elle = [0.0; 3]; // Angular momentum vector (r × v)
    let mut elv = [0.0; 3]; // Normalized angular momentum vector
    let mut xorb = [0.0; 3]; // Position in the orbital frame
    let mut vorb = [0.0; 3]; // Velocity in the orbital frame
    let mut rot = [[0.; 3]; 3]; // Rotation matrix to orbital frame

    // === 1. Compute angular momentum vector h = r × v
    elle[0] = xv[1] * xv[5] - xv[2] * xv[4];
    elle[1] = xv[2] * xv[3] - xv[0] * xv[5];
    elle[2] = xv[0] * xv[4] - xv[1] * xv[3];

    let el2 = elle.iter().map(|&x| x * x).sum::<f64>(); // ||h||²
    let elmod = sqrt(el2); // ||h||
    for j in 0..3 {
        elv[j] = elle[j] / elmod;
    }

    // === 2. Compute inclination i and longitude of ascending node Ω
    let sini = sqrt(elv[0].powi(2) + elv[1].powi(2));
    let mut ainc = atan2(sini, elv[2]); // inclination

    // Normalize inclination to [0, 2π]
    if ainc > DPI {
        ainc -= DPI;
    }
    if ainc < 0.0 {
        ainc += DPI;
    }

    // Ascending node longitude
    let anod = if sini == 0.0 {
        0.0 // Equatorial orbit: node undefined
    } else {
        let mut anod = atan2(elv[0], -elv[1]);
        if anod > DPI {
            anod -= DPI;
        }
        if anod < 0.0 {
            anod += DPI;
        }
        anod
    };

    // === 3. Transform position and velocity to the orbital frame
    let r1 = rotmt(ainc, 0); // Rotation around X-axis (inclination)
    let r2 = rotmt(anod, 2); // Rotation around Z-axis (node)
    prodmm(&mut rot, r1.into(), r2.into()); // rot = r1 * r2

    prodmv(&mut xorb, rot, xv[0..3].try_into().unwrap()); // orbital position
    prodmv(&mut vorb, rot, xv[3..6].try_into().unwrap()); // orbital velocity

    // === 4. Scalar products needed for orbital element derivation
    let rv = xorb[0] * vorb[0] + xorb[1] * vorb[1]; // r ⋅ v
    let rs = sqrt(xorb[0].powi(2) + xorb[1].powi(2)); // |r|
    let v2 = vorb[0].powi(2) + vorb[1].powi(2); // |v|²

    // === 5. Reciprocal of semi-major axis (1/a)
    let reca = 2.0 / rs - v2 / GAUSS_GRAV_SQUARED;

    // === CASE 1: Elliptical orbit (reca > 0)
    if reca > 0.0 {
        let sma = 1.0 / reca;
        let enne = sqrt(GAUSS_GRAV_SQUARED / sma.powi(3)); // mean motion

        // Eccentricity from e⋅sinE and e⋅cosE
        let esine = rv / (enne * sma * sma);
        let ecose = v2 * rs / GAUSS_GRAV_SQUARED - 1.0;
        let ecc = sqrt(esine.powi(2) + ecose.powi(2));

        let anec = atan2(esine, ecose); // eccentric anomaly E
        let mut emme = anec - ecc * anec.sin(); // mean anomaly M
        if emme < 0.0 {
            emme += DPI;
        }
        if emme > DPI {
            emme -= DPI;
        }

        // Argument of periapsis ω
        let x1 = anec.cos() - ecc;
        let rad = 1.0 - ecc * ecc;
        let x2 = sqrt(rad) * anec.sin();
        let xm = sqrt(x1 * x1 + x2 * x2);
        let x1 = x1 / xm;
        let x2 = x2 / xm;
        let sinper = x1 * xorb[1] - x2 * xorb[0];
        let cosper = x1 * xorb[0] + x2 * xorb[1];
        let argper = atan2(sinper, cosper);

        // Output: a, e, i, Ω, ω, M
        // elem[0] = sma;
        // elem[1] = ecc;
        // elem[2] = ainc;
        // elem[3] = anod;
        // elem[4] = argper;
        // elem[5] = emme;

        OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch,
            semi_major_axis: sma,
            eccentricity: ecc,
            inclination: ainc,
            ascending_node_longitude: anod,
            periapsis_argument: argper,
            mean_anomaly: emme,
        })
    }
    // === CASE 2: Parabolic orbit (reca == 0)
    else if reca == 0.0 {
        let ecc = 1.0;
        let p = el2 / GAUSS_GRAV_SQUARED; // semi-latus rectum
        let q = p / 2.0; // perihelion distance

        // True anomaly ν
        let cosf = p / rs - 1.0;
        let sinf = rv * p / (rs * elmod);
        let effe = atan2(sinf, cosf); // ν

        // Argument of periapsis ω
        let argper = atan2(xorb[1], xorb[0]) - effe;

        // Output: q, e, i, Ω, ω, ν
        // elem[0] = q;
        // elem[1] = ecc;
        // elem[2] = ainc;
        // elem[3] = anod;
        // elem[4] = argper;
        // elem[5] = effe;

        OrbitalElements::Cometary(CometaryElements {
            reference_epoch,
            perihelion_distance: q,
            eccentricity: ecc,
            inclination: ainc,
            ascending_node_longitude: anod,
            periapsis_argument: argper,
            true_anomaly: effe,
        })
    }
    // === CASE 3: Hyperbolic orbit (reca < 0)
    else {
        let p = el2 / GAUSS_GRAV_SQUARED;
        let ecosf = p / rs - 1.0;
        let esinf = rv * p / (elmod * rs);
        let effe = atan2(esinf, ecosf); // true anomaly ν
        let ecc = sqrt(ecosf.powi(2) + esinf.powi(2));
        let q = p / (1.0 + ecc); // perihelion distance

        // Argument of periapsis ω
        let argper = atan2(xorb[1], xorb[0]) - effe;

        // Output: q, e, i, Ω, ω, ν
        // elem[0] = q;
        // elem[1] = ecc;
        // elem[2] = ainc;
        // elem[3] = anod;
        // elem[4] = argper;
        // elem[5] = effe;

        OrbitalElements::Cometary(CometaryElements {
            reference_epoch,
            perihelion_distance: q,
            eccentricity: ecc,
            inclination: ainc,
            ascending_node_longitude: anod,
            periapsis_argument: argper,
            true_anomaly: effe,
        })
    }
}

/// Checks whether an orbit is dynamically acceptable based on eccentricity and perihelion distance limits.
///
/// This function computes the orbital eccentricity, perihelion distance, and specific orbital energy
/// from a Cartesian state vector using the Lenz–Runge vector and angular momentum. It returns a boolean
/// flag indicating whether the orbit satisfies user-defined thresholds on eccentricity and perihelion.
///
/// Arguments
/// ---------
/// * `asteroid_position` – Cartesian position vector of the object in AU.
/// * `asteroid_velocity` – Cartesian velocity vector of the object in AU/day.
/// * `peri_max` – Maximum allowed perihelion distance in AU.
/// * `ecc_max` – Maximum allowed eccentricity (dimensionless).
///
/// Returns
/// -------
/// Returns `Some((is_accepted, eccentricity, perihelion_distance, specific_energy))` where:
/// * `is_accepted` – `true` if both eccentricity < `ecc_max` and perihelion < `peri_max`
/// * `eccentricity` – Computed eccentricity (norm of Lenz–Runge vector)
/// * `perihelion_distance` – Pericenter distance in AU
/// * `specific_energy` – Two-body specific orbital energy (AU²/day²)
///
/// Returns `None` if the angular momentum vector is degenerate (null), indicating a singular or invalid orbit.
///
/// See also
/// --------
/// * [`ccek1`] – Computes classical orbital elements from a Cartesian state vector.
/// * [`KeplerianElements`] – Structured representation of orbital parameters.
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Computes orbital velocity from position triplets using the Lagrange formulation.
/// * [`GAUSS_GRAV_SQUARED`] – Gaussian gravitational constant squared.
pub fn eccentricity_control(
    asteroid_position: &Vector3<f64>,
    asteroid_velocity: &Vector3<f64>,
    peri_max: f64,
    ecc_max: f64,
) -> Option<(bool, f64, f64, f64)> {
    // Compute squared magnitude of velocity
    let ast_vel_2 = asteroid_velocity.dot(asteroid_velocity);

    // Compute distance from the central body
    let distance_to_center = asteroid_position.norm();

    // Angular momentum vector h = r × v
    let angular_momentum = asteroid_position.cross(asteroid_velocity);

    // Squared norm of angular momentum
    let angmom_norm = angular_momentum.dot(&angular_momentum);

    // If angular momentum is null, the orbit is singular or undefined
    if angmom_norm.sqrt() == 0. {
        return None;
    }

    // Lenz–Runge vector: points toward perihelion, used to derive eccentricity
    let lenz_prelim = asteroid_velocity.cross(&angular_momentum) * (1. / GAUSS_GRAV_SQUARED);
    let lenz_factor = asteroid_position * (1. / distance_to_center);
    let lenz_vector = lenz_prelim - lenz_factor;

    // Eccentricity is the norm of the Lenz vector
    let eccentricity = lenz_vector.norm();

    // Perihelion distance q = h² / [μ (1 + e)]
    let perihelie = angmom_norm / (GAUSS_GRAV_SQUARED * (1. + eccentricity));

    // Specific orbital energy: ε = v²/2 - μ/r
    let energy = ast_vel_2 / 2. - GAUSS_GRAV_SQUARED / distance_to_center;

    // Return acceptability flag + diagnostics
    Some((
        eccentricity < ecc_max && perihelie < peri_max,
        eccentricity,
        perihelie,
        energy,
    ))
}

#[cfg(test)]
mod orb_elem_test {
    use crate::kepler::principal_angle;

    use super::*;
    use approx::{abs_diff_eq, assert_abs_diff_eq};
    use nalgebra::Vector3;
    use proptest::prelude::*;

    const TWO_PI: f64 = DPI;

    fn energy_from_state(xv: &[f64; 6]) -> f64 {
        // ε = v²/2 - μ/r
        let r = (xv[0].powi(2) + xv[1].powi(2) + xv[2].powi(2)).sqrt();
        let v2 = xv[3].powi(2) + xv[4].powi(2) + xv[5].powi(2);
        0.5 * v2 - GAUSS_GRAV_SQUARED / r
    }

    fn norm3(x: &[f64; 3]) -> f64 {
        (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt()
    }

    fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }

    fn wrap_0_2pi(mut ang: f64) -> f64 {
        // Normalize angle to [0, 2π)
        ang %= TWO_PI;
        if ang < 0.0 {
            ang += TWO_PI;
        }
        ang
    }

    #[test]
    fn test_elem_regression_reference_with_tolerance() {
        let xv = [
            -0.623_550_051_003_163_9,
            1.211_468_114_860_160_5,
            0.252_000_591_437_760_4,
            -1.554_984_513_777_466_3E-2,
            -4.631_577_489_268_288E-3,
            -9.363_362_126_133_925E-4,
        ];

        let orbit = ccek1(&xv, 0.0);

        let ref_elem = [
            1.815_529_716_630_423_2,
            0.289_218_264_882_582_9,
            0.204_347_857_519_529_72,
            0.007_289_013_369_044_374_5,
            1.226_373_724_947_310_3,
            0.445_547_429_557_344_05,
        ];

        let kepler_orb = orbit.as_keplerian().unwrap();
        let elem = [
            kepler_orb.semi_major_axis,
            kepler_orb.eccentricity,
            kepler_orb.inclination,
            kepler_orb.ascending_node_longitude,
            kepler_orb.periapsis_argument,
            kepler_orb.mean_anomaly,
        ];
        // Compare with a sensible tolerance for double precision.
        for (got, exp) in elem.iter().zip(ref_elem.iter()) {
            assert_abs_diff_eq!(got, exp, epsilon = 5e-13);
        }
    }

    #[test]
    fn test_kepler_energy_invariant_when_elliptic() {
        // For elliptic orbits (type_ == "KEP"), we have ε = -μ/(2a).
        let xv = [
            -0.623_550_051_003_163_9,
            1.211_468_114_860_160_5,
            0.252_000_591_437_760_4,
            -1.554_984_513_777_466_3E-2,
            -4.631_577_489_268_288E-3,
            -9.363_362_126_133_925E-4,
        ];
        let eps = energy_from_state(&xv);

        let orbit = ccek1(&xv, 0.0);

        let kepler_orb = orbit.as_keplerian().unwrap();
        let a = kepler_orb.semi_major_axis;
        let energy_kepler = -GAUSS_GRAV_SQUARED / (2.0 * a);
        assert_abs_diff_eq!(eps, energy_kepler, epsilon = 5e-12);

        let i = kepler_orb.inclination;
        assert!((0.0 - 1e-15..=std::f64::consts::PI + 1e-15).contains(&i));
        for &ang in &[
            kepler_orb.ascending_node_longitude,
            kepler_orb.periapsis_argument,
            kepler_orb.mean_anomaly,
        ] {
            let w = wrap_0_2pi(ang);
            assert!((0.0..TWO_PI + 1e-15).contains(&w));
        }
    }

    #[test]
    fn test_hyperbolic_orbit_detection_and_diagnostics() {
        // Pick r = 1 AU along x, and v > v_esc along +y to ensure hyperbolic orbit.
        // v_esc = sqrt(2μ/r), with r=1 AU.
        let r = [1.0, 0.0, 0.0];
        let vesc = (2.0 * GAUSS_GRAV_SQUARED).sqrt();
        let v = [0.0, vesc * 1.2, 0.0]; // 20% above escape

        let xv = [r[0], r[1], r[2], v[0], v[1], v[2]];
        let orbit = ccek1(&xv, 0.0);

        match orbit {
            OrbitalElements::Cometary(ce) => {
                assert!(
                    ce.perihelion_distance > 0.0,
                    "Perihelion distance must be positive, got {q}",
                    q = ce.perihelion_distance
                );
                assert!(
                    ce.eccentricity > 1.0 + 1e-12,
                    "Hyperbolic e should be > 1, got {e}",
                    e = ce.eccentricity
                );
                assert!(ce.true_anomaly.is_finite());
            }
            OrbitalElements::Keplerian(_) | OrbitalElements::Equinoctial(_) => {
                panic!("Expected Cometary elements, got {orbit:#?}")
            }
        }
    }

    #[test]
    fn test_quasi_parabolic_behaviour() {
        // Construct a nearly parabolic orbit: r = 1 AU, v ~ v_escape
        // Escape velocity at r=1 AU: v_esc = sqrt(2μ/r)
        let r = [1.0, 0.0, 0.0];
        let vesc = (2.0 * GAUSS_GRAV_SQUARED).sqrt();

        // Slightly below escape velocity → almost parabolic ellipse
        let v = [0.0, vesc * 0.999_999, 0.0];
        let xv: [f64; 6] = [r[0], r[1], r[2], v[0], v[1], v[2]];

        // For parabolic orbits, the specific orbital energy ε = 0
        // Verify that the energy is close to zero
        let eps = {
            let rnorm = (xv[0].powi(2) + xv[1].powi(2) + xv[2].powi(2)).sqrt();
            let v2 = xv[3].powi(2) + xv[4].powi(2) + xv[5].powi(2);
            0.5 * v2 - GAUSS_GRAV_SQUARED / rnorm
        };
        assert!(
            eps.abs() < 1e-8,
            "Near-parabolic specific energy should be ~0, got {eps}"
        );

        let orbit = ccek1(&xv, 0.0);

        match orbit {
            OrbitalElements::Keplerian(ke) => {
                // Eccentricity should be very close to 1
                assert!(
                    (ke.eccentricity - 1.0).abs() < 1e-5,
                    "Near-parabolic eccentricity should be ~1, got {e}",
                    e = ke.eccentricity
                );
            }

            OrbitalElements::Cometary(ce) => {
                assert!(
                    (ce.eccentricity - 1.0).abs() < 1e-5,
                    "Near-parabolic eccentricity should be ~1, got {e}",
                    e = ce.eccentricity
                );
            }

            OrbitalElements::Equinoctial(_) => {
                panic!("Expected Keplerian or Cometary elements, got {orbit:#?}")
            }
        }
    }

    #[test]
    fn test_eccentricity_control() {
        let asteroid_position = Vector3::new(
            -0.623_550_051_003_163_9,
            1.011_260_185_597_692,
            0.713_100_363_506_241_4,
        );

        let asteroid_velocity = Vector3::new(
            -1.554_984_513_777_466_3E-2,
            -3.876_936_109_837_657_7E-3,
            -2.701_407_400_297_996_4E-3,
        );

        let (accept_ecc, ecc, peri, energy) =
            eccentricity_control(&asteroid_position, &asteroid_velocity, 1e3, 2.).unwrap();

        assert!(accept_ecc);
        assert_eq!(ecc, 0.2892182648825829);
        assert_eq!(peri, 1.2904453621438048);
        assert_eq!(energy, -0.00008149473004352595);
    }

    #[test]
    fn test_eccentricity_control_degenerate_ang_momentum_returns_none() {
        // h = 0 when r and v are parallel. Example: r=[1,0,0], v=[1,0,0].
        let r = Vector3::new(1.0, 0.0, 0.0);
        let v = Vector3::new(1.0, 0.0, 0.0);

        let res = eccentricity_control(&r, &v, 1e3, 2.0);
        assert!(
            res.is_none(),
            "Degenerate angular momentum must return None"
        );
    }

    #[test]
    fn test_eccentricity_control_matches_lenz_vector_definition() {
        // Cross-check: eccentricity_control returns |e_vec| consistent with the Lenz–Runge vector.
        let r = Vector3::new(0.5, 0.2, 0.1);
        let v = Vector3::new(0.0, 0.03, -0.01);

        let (accepted, e_ctrl, q, energy) =
            eccentricity_control(&r, &v, 1e3, 2.0).expect("non-degenerate");

        // Manual Lenz–Runge computation for verification
        let h = r.cross(&v);
        let e_vec = v.cross(&h) * (1.0 / GAUSS_GRAV_SQUARED) - r * (1.0 / r.norm());
        let e_manual = e_vec.norm();

        assert!(accepted); // thresholds are huge here
        assert_abs_diff_eq!(e_ctrl, e_manual, epsilon = 1e-12);
        assert!(q > 0.0);
        assert!(energy.is_finite());
    }

    #[test]
    fn test_invariants_a_e_i_under_rotation_about_z() {
        // Rotate state about Z: a,e,i should be invariant.
        let xv = [0.8, 0.3, 0.1, 0.0, 0.02, -0.005];
        let orbit = ccek1(&xv, 0.0);

        // Build a pure Z-rotation of phi using rotmt:
        let phi = 0.73;
        let r_z = rotmt(phi, 2);

        // Apply rotation to position and velocity
        let mut rot = [[0.0; 3]; 3];
        // identity * r_z
        prodmm(
            &mut rot,
            r_z.into(),
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
        );
        let p = [xv[0], xv[1], xv[2]];
        let v = [xv[3], xv[4], xv[5]];
        let mut pr = [0.0; 3];
        let mut vr = [0.0; 3];
        prodmv(&mut pr, rot, p);
        prodmv(&mut vr, rot, v);

        let xvr = [pr[0], pr[1], pr[2], vr[0], vr[1], vr[2]];
        let orbit2 = ccek1(&xvr, 0.0);

        // Extract invariants a, e, i (Keplerian) or q, e, i (Cometary).
        let (a1, e1, i1, om1) = match &orbit {
            OrbitalElements::Keplerian(k) => (
                k.semi_major_axis,
                k.eccentricity,
                k.inclination,
                k.ascending_node_longitude,
            ),
            OrbitalElements::Cometary(c) => (
                c.perihelion_distance,
                c.eccentricity,
                c.inclination,
                c.ascending_node_longitude,
            ),
            OrbitalElements::Equinoctial(_) => panic!("Unexpected equinoctial here"),
        };
        let (a2, e2, i2, om2) = match &orbit2 {
            OrbitalElements::Keplerian(k) => (
                k.semi_major_axis,
                k.eccentricity,
                k.inclination,
                k.ascending_node_longitude,
            ),
            OrbitalElements::Cometary(c) => (
                c.perihelion_distance,
                c.eccentricity,
                c.inclination,
                c.ascending_node_longitude,
            ),
            OrbitalElements::Equinoctial(_) => panic!("Unexpected equinoctial here"),
        };

        // Invariants under active Z-rotation of (r, v):
        assert!(abs_diff_eq!(a1, a2, epsilon = 1e-12));
        assert!(abs_diff_eq!(e1, e2, epsilon = 1e-12));
        assert!(abs_diff_eq!(i1, i2, epsilon = 1e-12));

        // Optional: check Ω shift ≈ ±φ (depends on your rotation convention / order).
        // We accept either +φ or -φ modulo 2π to be robust to matrix order.
        let d_om = principal_angle(om2 - om1);
        let plus_phi = principal_angle(phi);
        let minus_phi = principal_angle(-phi);
        let ok_plus = abs_diff_eq!(d_om, plus_phi, epsilon = 1e-12);
        let ok_minus = abs_diff_eq!(d_om, minus_phi, epsilon = 1e-12);
        assert!(
            ok_plus || ok_minus,
            "Unexpected node shift: ΔΩ={d_om}, expected ±{phi}"
        );
    }

    // ---------------------------
    // Property-based tests (proptest)
    // ---------------------------

    proptest! {
        // Generate random but reasonable heliocentric states:
        // r in [0.3, 5] AU; v in [0, 0.05] AU/day
        // Avoid h ~ 0 by rejecting near-parallel r and v.
        #[test]
        fn prop_state_to_elements_basic_sanity(
            rx in 0.3f64..5.0, ry in -5.0f64..5.0, rz in -2.0f64..2.0,
            vx in -0.05f64..0.05, vy in -0.05f64..0.05, vz in -0.05f64..0.05
        ) {
            let r = [rx, ry, rz];
            let v = [vx, vy, vz];
            let rnorm = norm3(&r);
            prop_assume!(rnorm > 0.3);

            let h = cross(&r, &v);
            let hnorm = norm3(&h);
            // reject nearly singular angular momentum
            prop_assume!(hnorm > 1e-6);

            let xv = [r[0], r[1], r[2], v[0], v[1], v[2]];
            let orbit = ccek1(&xv, 0.0);

            // Compare eccentricity with eccentricity_control
            let (_, e_ctrl, q_ctrl, _) = eccentricity_control(
                &Vector3::new(r[0], r[1], r[2]),
                &Vector3::new(v[0], v[1], v[2]),
                1e4, 5.0
            ).expect("non-degenerate");
            prop_assert!(e_ctrl >= 0.0);
            prop_assert!(q_ctrl > 0.0);

            match orbit {
                OrbitalElements::Keplerian(ke) => {
                // e from both paths should match (parabolic borderline may differ by ~1e-6)
                prop_assert!((e_ctrl - ke.eccentricity).abs() < 1e-6);

                // Elliptic: a > 0, e in [0,1), M in [0,2π), ε = -μ/(2a)
                let a = ke.semi_major_axis;
                prop_assert!(a > 0.0);
                prop_assert!((0.0..1.0 + 1e-12).contains(&ke.eccentricity));

                let m = wrap_0_2pi(ke.mean_anomaly);
                prop_assert!((0.0..TWO_PI + 1e-12).contains(&m));

                let eps = energy_from_state(&xv);
                let energy_kepler = -GAUSS_GRAV_SQUARED / (2.0 * a);
                prop_assert!((eps - energy_kepler).abs() < 1e-7);

                // Angles are finite
                for &ang in &[ke.inclination, ke.ascending_node_longitude, ke.periapsis_argument, ke.mean_anomaly] {
                    let w = wrap_0_2pi(ang);
                    prop_assert!((0.0..TWO_PI + 1e-12).contains(&w));
                }

                // Inclination in [0,π]
                let inc = ke.inclination;
                prop_assert!((-1e-12..=std::f64::consts::PI + 1e-12).contains(&inc));
                }

                OrbitalElements::Cometary(ce) => {
                    // e from both paths should match (parabolic borderline may differ by ~1e-6)
                    prop_assert!((e_ctrl - ce.eccentricity).abs() < 1e-6);

                    // Cometary (parabolic/hyperbolic): e >= 1
                    prop_assert!(ce.eccentricity >= 1.0 - 1e-6);
                    prop_assert!(ce.perihelion_distance > 0.0);

                    // Angles are finite
                    for &ang in &[ce.inclination, ce.ascending_node_longitude, ce.periapsis_argument, ce.true_anomaly] {
                        let w = wrap_0_2pi(ang);
                        prop_assert!((0.0..TWO_PI + 1e-12).contains(&w));
                    }

                    // Inclination in [0,π]
                    let inc = ce.inclination;
                    prop_assert!((-1e-12..=std::f64::consts::PI + 1e-12).contains(&inc));
                }

                OrbitalElements::Equinoctial(_) => {
                    // Not tested here
                    prop_assert!(false, "Equinoctial elements not tested");
                }
            }
        }
    }
}
