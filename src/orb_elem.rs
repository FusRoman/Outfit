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
            for l in 0..3 {
                s += w[j][l] * z[l][k];
            }
            a[j][k] = s;
        }
    }
}

fn prodmv(y: &mut [f64; 3], a: [[f64; 3]; 3], x: [f64; 3]) {
    let z = x;

    for j in 0..3 {
        let mut s = 0.0;
        for l in 0..3 {
            s += a[j][l] * z[l];
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
/// # Arguments
/// * `elem` – A mutable array to hold the resulting orbital elements (output).
/// * `type_` – A mutable string that will be set to `"KEP"` (elliptic) or `"COM"` (cometary/parabolic/hyperbolic).
/// * `xv` – A 6-element array `[x, y, z, vx, vy, vz]` in heliocentric J2000 coordinates.
///          Positions are in astronomical units (AU) and velocities in AU/day.
///
/// # Note
/// This function does not account for planetary perturbations or relativistic corrections.
/// The orbital elements are computed under the assumption of a two-body Keplerian motion.
///
/// # See also
/// * [`eccentricity_control`] – Filters orbital solutions based on eccentricity and perihelion distance thresholds.
/// * [`KeplerianElements`] – Struct holding canonical orbital elements.
/// * [`velocity_correction`] – Lagrange-based velocity refinement from triplet position vectors.
/// * [`GAUSS_GRAV_SQUARED`](crate::constants::GAUSS_GRAV_SQUARED) – Gaussian gravitational constant.
pub fn ccek1(elem: &mut [f64; 6], type_: &mut String, xv: &[f64; 6]) {
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
    prodmm(&mut rot, r1, r2); // rot = r1 * r2

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
        *type_ = String::from("KEP");
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
        elem[0] = sma;
        elem[1] = ecc;
        elem[2] = ainc;
        elem[3] = anod;
        elem[4] = argper;
        elem[5] = emme;
    }
    // === CASE 2: Parabolic orbit (reca == 0)
    else if reca == 0.0 {
        *type_ = String::from("COM");
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
        elem[0] = q;
        elem[1] = ecc;
        elem[2] = ainc;
        elem[3] = anod;
        elem[4] = argper;
        elem[5] = effe;
    }
    // === CASE 3: Hyperbolic orbit (reca < 0)
    else {
        *type_ = String::from("COM");
        let p = el2 / GAUSS_GRAV_SQUARED;
        let ecosf = p / rs - 1.0;
        let esinf = rv * p / (elmod * rs);
        let effe = atan2(esinf, ecosf); // true anomaly ν
        let ecc = sqrt(ecosf.powi(2) + esinf.powi(2));
        let q = p / (1.0 + ecc); // perihelion distance

        // Argument of periapsis ω
        let argper = atan2(xorb[1], xorb[0]) - effe;

        // Output: q, e, i, Ω, ω, ν
        elem[0] = q;
        elem[1] = ecc;
        elem[2] = ainc;
        elem[3] = anod;
        elem[4] = argper;
        elem[5] = effe;
    }
}

/// Checks whether an orbit is dynamically acceptable based on eccentricity and perihelion distance limits.
///
/// This function computes the orbital eccentricity, perihelion distance, and specific orbital energy
/// from a Cartesian state vector using the Lenz–Runge vector and angular momentum. It returns a boolean
/// flag indicating whether the orbit satisfies user-defined thresholds on eccentricity and perihelion.
///
/// # Arguments
/// * `asteroid_position` – Cartesian position vector of the object in AU.
/// * `asteroid_velocity` – Cartesian velocity vector of the object in AU/day.
/// * `peri_max` – Maximum allowed perihelion distance in AU.
/// * `ecc_max` – Maximum allowed eccentricity (dimensionless).
///
/// # Returns
/// Returns `Some((is_accepted, eccentricity, perihelion_distance, specific_energy))` where:
/// * `is_accepted` – `true` if both eccentricity < `ecc_max` and perihelion < `peri_max`
/// * `eccentricity` – Computed eccentricity (norm of Lenz–Runge vector)
/// * `perihelion_distance` – Pericenter distance in AU
/// * `specific_energy` – Two-body specific orbital energy (AU²/day²)
///
/// Returns `None` if the angular momentum vector is degenerate (null), indicating a singular or invalid orbit.
///
/// # See also
/// * [`ccek1`] – Computes classical orbital elements from a Cartesian state vector.
/// * [`KeplerianElements`] – Structured representation of orbital parameters.
/// * [`velocity_correction`] – Computes orbital velocity from position triplets using the Lagrange formulation.
/// * [`GAUSS_GRAV_SQUARED`](crate::constants::GAUSS_GRAV_SQUARED) – Gaussian gravitational constant squared.
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

    use super::*;

    #[test]
    fn test_elem() {
        let mut elem = [0.0; 6];
        let mut type_ = String::new();
        let xv = [
            -0.623_550_051_003_163_9,
            1.2114681148601605,
            0.252_000_591_437_760_4,
            -1.554_984_513_777_466_3E-2,
            -4.631_577_489_268_288E-3,
            -9.363_362_126_133_925E-4,
        ];

        ccek1(&mut elem, &mut type_, &xv);

        let ref_elem = [
            1.8155297166304232,
            0.2892182648825829,
            0.20434785751952972,
            0.0072890133690443745,
            1.2263737249473103,
            0.44554742955734405,
        ];

        assert_eq!(ref_elem, elem)
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
}
