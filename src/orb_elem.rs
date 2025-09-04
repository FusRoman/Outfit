use crate::orbit_type::cometary_element::CometaryElements;
use crate::orbit_type::keplerian_element::KeplerianElements;
use crate::orbit_type::OrbitalElements;

use super::constants::{DPI, GAUSS_GRAV_SQUARED};
use super::ref_system::rotmt;
use nalgebra::Vector3;

/// Convert a Cartesian heliocentric state vector into orbital elements in the J2000 equatorial frame.
///
/// This function computes the classical orbital elements from a state vector `[r, v]` expressed
/// in heliocentric J2000 coordinates. It distinguishes between **elliptical**, **parabolic**,
/// and **hyperbolic** orbits based on the reciprocal semi-major axis `1/a`.
///
/// Elliptical orbits (`1/a > 0`) return:
/// -----------------
/// * `a` – Semi-major axis (AU)  
/// * `e` – Eccentricity  
/// * `i` – Inclination (radians)  
/// * `Ω` – Longitude of ascending node (radians)  
/// * `ω` – Argument of periapsis (radians)  
/// * `M` – Mean anomaly (radians)  
///
/// Parabolic orbits (`1/a = 0`) and hyperbolic orbits (`1/a < 0`) return:
/// -----------------
/// * `q` – Perihelion distance (AU)  
/// * `e` – Eccentricity (= 1 for parabolic)  
/// * `i` – Inclination (radians)  
/// * `Ω` – Longitude of ascending node (radians)  
/// * `ω` – Argument of periapsis (radians)  
/// * `ν` – True anomaly (radians)  
///
/// Arguments
/// -----------------
/// * `position`: Position vector `[x, y, z]` in AU, heliocentric J2000.  
/// * `velocity`: Velocity vector `[vx, vy, vz]` in AU/day, heliocentric J2000.  
/// * `reference_epoch`: Epoch of the state vector in Julian Date (TDB).  
///
/// Return
/// ----------
/// * An [`OrbitalElements`] enum, which can be either:  
///   - [`KeplerianElements`] for elliptic solutions, or  
///   - [`CometaryElements`] for parabolic and hyperbolic solutions.  
///
/// Note
/// ----------
/// * The computation assumes a pure two-body Keplerian motion.  
/// * Planetary perturbations, relativistic corrections, and non-gravitational effects are **not** included.  
/// * For equatorial orbits (`i = 0`), the longitude of the ascending node is set to zero by convention.  
///
/// See also
/// ------------
/// * [`KeplerianElements`] – Canonical set of orbital elements.  
/// * [`CometaryElements`] – Orbital elements for parabolic and hyperbolic solutions.  
/// * [`eccentricity_control`] – Filters orbital solutions based on eccentricity and perihelion distance thresholds.  
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Lagrange-based velocity refinement.  
/// * [`GAUSS_GRAV_SQUARED`] – Gaussian gravitational constant squared.
pub(crate) fn ccek1(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
    reference_epoch: f64,
) -> OrbitalElements {
    // ---- Helpers & constants ------------------------------------------------
    #[inline]
    fn wrap_0_2pi(x: f64) -> f64 {
        x.rem_euclid(DPI)
    }

    #[inline]
    fn squared_norm(v: &Vector3<f64>) -> f64 {
        v.dot(v)
    }

    const EPS_EQ: f64 = 1e-15; // threshold for equatorial orbits (sin(i) ≈ 0)
    const EPS_PARAB: f64 = 1e-12; // threshold for nearly parabolic orbits (1/a ≈ 0)
    const EPS_E: f64 = 5e-15; // threshold for eccentricity ≈ 1

    let mu = GAUSS_GRAV_SQUARED;

    // ---- 1) Angular momentum and orbital plane orientation ------------------
    // Specific angular momentum vector: h = r × v
    let angular_momentum = position.cross(velocity);
    let angular_momentum_squared = squared_norm(&angular_momentum);
    let angular_momentum_norm = angular_momentum_squared.sqrt();
    let angular_momentum_unit = angular_momentum / angular_momentum_norm;

    // Inclination i = atan2(||projection in XY plane||, z-component)
    let sin_inclination = angular_momentum_unit.xy().norm();
    let mut inclination = wrap_0_2pi(sin_inclination.atan2(angular_momentum_unit.z));

    // Longitude of ascending node Ω
    let ascending_node_longitude = if sin_inclination <= EPS_EQ {
        // Equatorial orbit: conventionally set i = 0 and Ω = 0
        inclination = 0.0;
        0.0
    } else {
        wrap_0_2pi(angular_momentum_unit.x.atan2(-angular_momentum_unit.y))
    };

    // ---- 2) Rotate to orbital frame (orbital plane ≈ XY plane) --------------
    let rotation_inclination = rotmt(inclination, 0); // rotation around X-axis
    let rotation_node = rotmt(ascending_node_longitude, 2); // rotation around Z-axis
    let orbital_rotation = rotation_inclination.transpose() * rotation_node.transpose();

    let position_orbital = orbital_rotation * position;
    let velocity_orbital = orbital_rotation * velocity;

    // In-plane scalars (avoid z-component to reduce numerical noise)
    let radial_velocity_dot = position_orbital.xy().dot(&velocity_orbital.xy()); // r · v
    let radial_distance = position_orbital.xy().norm(); // |r|
    let velocity_squared = velocity_orbital.xy().norm_squared(); // |v|²

    // Reciprocal semi-major axis: 1/a = 2/|r| − |v|²/μ
    let reciprocal_semi_major_axis = 2.0 / radial_distance - velocity_squared / mu;

    // Closure: compute periapsis argument ω from true anomaly ν
    let periapsis_argument_from_true_anomaly = |true_anomaly: f64| -> f64 {
        wrap_0_2pi(position_orbital.y.atan2(position_orbital.x) - true_anomaly)
    };

    // Closure: return cometary elements in the parabolic case
    let parabolic_solution = || -> OrbitalElements {
        let eccentricity = 1.0;
        let semi_latus_rectum = angular_momentum_squared / mu;
        let perihelion_distance = semi_latus_rectum / 2.0;
        let cos_true_anomaly = semi_latus_rectum / radial_distance - 1.0;
        let sin_true_anomaly =
            radial_velocity_dot * semi_latus_rectum / (radial_distance * angular_momentum_norm);
        let true_anomaly = sin_true_anomaly.atan2(cos_true_anomaly);
        let periapsis_argument = periapsis_argument_from_true_anomaly(true_anomaly);

        OrbitalElements::Cometary(CometaryElements {
            reference_epoch,
            perihelion_distance,
            eccentricity,
            inclination,
            ascending_node_longitude,
            periapsis_argument,
            true_anomaly,
        })
    };

    // ---- 3) Orbit classification --------------------------------------------
    if reciprocal_semi_major_axis > EPS_PARAB {
        // -------- Elliptical orbit
        let semi_major_axis = 1.0 / reciprocal_semi_major_axis;
        let mean_motion = (mu / semi_major_axis.powi(3)).sqrt();

        // Components e·sinE and e·cosE
        let eccentricity_sin_ea =
            radial_velocity_dot / (mean_motion * semi_major_axis * semi_major_axis);
        let eccentricity_cos_ea = velocity_squared * radial_distance / mu - 1.0;
        let eccentricity = (eccentricity_sin_ea.powi(2) + eccentricity_cos_ea.powi(2)).sqrt();

        // Near-parabolic degeneracy
        if (eccentricity - 1.0).abs() < EPS_E {
            return parabolic_solution();
        }

        // Eccentric anomaly and mean anomaly
        let eccentric_anomaly = eccentricity_sin_ea.atan2(eccentricity_cos_ea);
        let mean_anomaly = wrap_0_2pi(eccentric_anomaly - eccentricity * eccentric_anomaly.sin());

        // Argument of periapsis ω
        let x1 = eccentric_anomaly.cos() - eccentricity;
        let rad = (1.0 - eccentricity * eccentricity).sqrt();
        let x2 = rad * eccentric_anomaly.sin();
        let norm = (x1 * x1 + x2 * x2).sqrt();
        let (x1n, x2n) = (x1 / norm, x2 / norm);

        let sin_periapsis = x1n * position_orbital.y - x2n * position_orbital.x;
        let cos_periapsis = x1n * position_orbital.x + x2n * position_orbital.y;
        let periapsis_argument = wrap_0_2pi(sin_periapsis.atan2(cos_periapsis));

        OrbitalElements::Keplerian(KeplerianElements {
            reference_epoch,
            semi_major_axis,
            eccentricity,
            inclination,
            ascending_node_longitude,
            periapsis_argument,
            mean_anomaly,
        })
    } else if reciprocal_semi_major_axis.abs() <= EPS_PARAB {
        // -------- Parabolic orbit
        parabolic_solution()
    } else {
        // -------- Hyperbolic orbit
        let semi_latus_rectum = angular_momentum_squared / mu;
        let ecc_cos_true_anomaly = semi_latus_rectum / radial_distance - 1.0;
        let ecc_sin_true_anomaly =
            radial_velocity_dot * semi_latus_rectum / (angular_momentum_norm * radial_distance);
        let true_anomaly = ecc_sin_true_anomaly.atan2(ecc_cos_true_anomaly);
        let eccentricity = (ecc_cos_true_anomaly.powi(2) + ecc_sin_true_anomaly.powi(2)).sqrt();

        // Near-parabolic degeneracy
        if (eccentricity - 1.0).abs() < EPS_E {
            return parabolic_solution();
        }

        let perihelion_distance = semi_latus_rectum / (1.0 + eccentricity);
        let periapsis_argument = periapsis_argument_from_true_anomaly(true_anomaly);

        OrbitalElements::Cometary(CometaryElements {
            reference_epoch,
            perihelion_distance,
            eccentricity,
            inclination,
            ascending_node_longitude,
            periapsis_argument,
            true_anomaly,
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
/// * [`OrbitalElements::from_orbital_state`] – Computes classical orbital elements from a Cartesian state vector.
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

    pub fn energy_from_state(position: Vector3<f64>, velocity: Vector3<f64>) -> f64 {
        let r = position.norm(); // |r|
        let v2 = velocity.norm_squared(); // |v|²
        0.5 * v2 - GAUSS_GRAV_SQUARED / r
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
        use nalgebra::Vector3;

        let position = Vector3::new(
            -0.623_550_051_003_163_9,
            1.211_468_114_860_160_5,
            0.252_000_591_437_760_4,
        );

        let velocity = Vector3::new(
            -1.554_984_513_777_466_3E-2,
            -4.631_577_489_268_288E-3,
            -9.363_362_126_133_925E-4,
        );

        let orbit = ccek1(&position, &velocity, 0.0);

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

        // Position (AU)
        let position = Vector3::new(
            -0.623_550_051_003_163_9,
            1.211_468_114_860_160_5,
            0.252_000_591_437_760_4,
        );

        // Velocity (AU/day)
        let velocity = Vector3::new(
            -1.554_984_513_777_466_3E-2,
            -4.631_577_489_268_288E-3,
            -9.363_362_126_133_925E-4,
        );

        // Specific orbital energy ε = v²/2 - μ/r
        let eps = energy_from_state(position, velocity);

        // Build orbital elements from state vector
        let orbit = ccek1(&position, &velocity, 0.0);

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
        // v_esc = sqrt(2μ/r), with r = 1 AU.
        let position: Vector3<f64> = Vector3::new(1.0, 0.0, 0.0);
        let vesc = (2.0 * GAUSS_GRAV_SQUARED / position.norm()).sqrt();
        let velocity = Vector3::new(0.0, vesc * 1.2, 0.0); // 20% above escape

        let orbit = ccek1(&position, &velocity, 0.0);

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
        let position: Vector3<f64> = Vector3::new(1.0, 0.0, 0.0);

        // Escape velocity at r=1 AU: v_esc = sqrt(2μ/r)
        let vesc: f64 = (2.0 * GAUSS_GRAV_SQUARED / position.norm()).sqrt();

        // Slightly below escape velocity → almost parabolic ellipse
        let velocity = Vector3::new(0.0, vesc * 0.999_999, 0.0);

        // For parabolic orbits, the specific orbital energy ε ≈ 0
        let eps = {
            let rnorm = position.norm();
            let v2 = velocity.norm_squared();
            0.5 * v2 - GAUSS_GRAV_SQUARED / rnorm
        };

        assert!(
            eps.abs() < 1e-8,
            "Near-parabolic specific energy should be ~0, got {eps}"
        );

        let orbit = ccek1(&position, &velocity, 0.0);

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
        // Rotate state about Z: a, e, i should be invariant.
        let position = Vector3::new(0.8, 0.3, 0.1);
        let velocity = Vector3::new(0.0, 0.02, -0.005);

        let orbit = ccek1(&position, &velocity, 0.0);

        // Build a pure Z-rotation of phi using rotmt:
        let phi = 0.73;
        let r_z = rotmt(phi, 2); // assumes rotmt returns a 3×3 rotation matrix

        // Apply rotation to position and velocity (no need for identity or manual products)
        let pr = r_z * position; // rotated position
        let vr = r_z * velocity; // rotated velocity

        let orbit2 = ccek1(&pr, &vr, 0.0);

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
        // r in [0.3, 5] AU; v in [-0.05, 0.05] AU/day
        // Avoid h ~ 0 by rejecting near-parallel r and v.
        #[test]
        fn prop_state_to_elements_basic_sanity(
            rx in 0.3f64..5.0,  ry in -5.0f64..5.0,  rz in -2.0f64..2.0,
            vx in -0.05f64..0.05, vy in -0.05f64..0.05, vz in -0.05f64..0.05
        ) {
            // State vectors
            let r = Vector3::new(rx, ry, rz);
            let v = Vector3::new(vx, vy, vz);

            // Reject too-small |r| (degenerate orbits)
            prop_assume!(r.norm() > 0.3);

            // Reject nearly singular angular momentum h = r × v
            let h = r.cross(&v);
            prop_assume!(h.norm() > 1e-6);

            // Convert state → elements
            let orbit = ccek1(&r, &v, 0.0);

            // Compare eccentricity against the control path
            let (_, e_ctrl, q_ctrl, _) = eccentricity_control(&r, &v, 1e4, 5.0).expect("non-degenerate");
            prop_assert!(e_ctrl >= 0.0);
            prop_assert!(q_ctrl > 0.0);

            // Helper to normalize any angle to [0, 2π)
            let wrap_0_2pi = |ang: f64| ang.rem_euclid(TWO_PI);

            match orbit {
                OrbitalElements::Keplerian(ke) => {
                    // Eccentricity must match the control within tolerance
                    prop_assert!((e_ctrl - ke.eccentricity).abs() < 1e-6);

                    // Elliptic: a > 0, e in [0,1), M in [0,2π), ε = -μ/(2a)
                    let a = ke.semi_major_axis;
                    prop_assert!(a > 0.0);
                    prop_assert!((0.0..1.0 + 1e-12).contains(&ke.eccentricity));

                    let m = wrap_0_2pi(ke.mean_anomaly);
                    prop_assert!((0.0..TWO_PI + 1e-12).contains(&m));

                    let eps = energy_from_state(r, v);
                    let energy_kepler = -GAUSS_GRAV_SQUARED / (2.0 * a);
                    prop_assert!((eps - energy_kepler).abs() < 1e-7);

                    // Angles are finite and normalized
                    for &ang in &[
                        ke.inclination,
                        ke.ascending_node_longitude,
                        ke.periapsis_argument,
                        ke.mean_anomaly,
                    ] {
                        let w = wrap_0_2pi(ang);
                        prop_assert!((0.0..TWO_PI + 1e-12).contains(&w));
                    }

                    // Inclination in [0, π]
                    let inc = ke.inclination;
                    prop_assert!((-1e-12..=std::f64::consts::PI + 1e-12).contains(&inc));
                }

                OrbitalElements::Cometary(ce) => {
                    // Eccentricity must match the control within tolerance
                    prop_assert!((e_ctrl - ce.eccentricity).abs() < 1e-6);

                    // Cometary (parabolic/hyperbolic): e >= 1
                    prop_assert!(ce.eccentricity >= 1.0 - 1e-6);
                    prop_assert!(ce.perihelion_distance > 0.0);

                    // Angles are finite and normalized
                    for &ang in &[
                        ce.inclination,
                        ce.ascending_node_longitude,
                        ce.periapsis_argument,
                        ce.true_anomaly,
                    ] {
                        let w = wrap_0_2pi(ang);
                        prop_assert!((0.0..TWO_PI + 1e-12).contains(&w));
                    }

                    // Inclination in [0, π]
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
