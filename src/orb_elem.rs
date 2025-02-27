use super::constants::GaussGravSquared;
use super::ref_system::rotmt;
use nalgebra::Vector3;

use std::f64;
use std::f64::consts::PI;

const EPS: f64 = 5e-15;
const DOUBLEPI: f64 = PI * 2.;

fn atan2(y: f64, x: f64) -> f64 {
    y.atan2(x)
}

fn sqrt(x: f64) -> f64 {
    x.sqrt()
}

fn prodmm(a: &mut [[f64; 3]; 3], b: [[f64; 3]; 3], c: [[f64; 3]; 3]) {
    let mut w = b;
    let mut z = c;

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
    let z = x; // Copie du vecteur x

    for j in 0..3 {
        let mut s = 0.0;
        for l in 0..3 {
            s += a[j][l] * z[l];
        }
        y[j] = s;
    }
}

pub fn ccek1(elem: &mut [f64; 6], type_: &mut String, xv: &[f64; 6]) {
    let mut elle = [0.0; 3];
    let mut elv = [0.0; 3];
    let mut xorb = [0.0; 3];
    let mut vorb = [0.0; 3];
    let mut rot = [[0.; 3]; 3];

    // Angular momentum vector
    elle[0] = xv[1] * xv[5] - xv[2] * xv[4];
    elle[1] = xv[2] * xv[3] - xv[0] * xv[5];
    elle[2] = xv[0] * xv[4] - xv[1] * xv[3];

    // Angular momentum unit vector
    let el2 = elle.iter().map(|&x| x * x).sum::<f64>();

    let elmod = sqrt(el2);
    for j in 0..3 {
        elv[j] = elle[j] / elmod;
    }

    // Orbital inclination and longitude of the node
    let sini = sqrt(elv[0] * elv[0] + elv[1] * elv[1]);
    let mut ainc = atan2(sini, elv[2]);

    if ainc > DOUBLEPI {
        ainc -= DOUBLEPI;
    }
    if ainc < 0.0 {
        ainc += DOUBLEPI;
    }
    let mut anod = if sini == 0.0 {
        0.0
    } else {
        let mut anod = atan2(elv[0], -elv[1]);
        if anod > DOUBLEPI {
            anod -= DOUBLEPI;
        }
        if anod < 0.0 {
            anod += DOUBLEPI;
        }
        anod
    };

    // Cartesian coordinates with respect to the "orbital frame" (with X axis directed along the line of nodes)
    let mut r1 = rotmt(ainc, 0);
    let mut r2 = rotmt(anod, 2);

    prodmm(&mut rot, r1, r2);
    prodmv(&mut xorb, rot, xv[0..3].try_into().unwrap());
    prodmv(&mut vorb, rot, xv[3..6].try_into().unwrap());

    let rv = xorb[0] * vorb[0] + xorb[1] * vorb[1];

    // Reciprocal semimajor axis
    let rs = sqrt(xorb[0] * xorb[0] + xorb[1] * xorb[1]);
    let v2 = vorb[0] * vorb[0] + vorb[1] * vorb[1];

    let mut reca = 2.0 / rs - v2 / GaussGravSquared;

    // CASE 1: RECA > 0 (elliptic orbit)
    if reca > 0.0 {
        *type_ = String::from("KEP");
        let sma = 1.0 / reca;
        let enne = sqrt(GaussGravSquared / sma.powi(3));

        // Eccentricity
        let esine = rv / (enne * sma * sma);
        let ecose = v2 * rs / GaussGravSquared - 1.0;
        let ecc = sqrt(esine * esine + ecose * ecose);

        if (ecc - 1.0).abs() < EPS {
            reca = 0.0;
        } else {
            let anec = atan2(esine, ecose);
            let mut emme = anec - ecc * anec.sin();
            if emme < 0.0 {
                emme += DOUBLEPI;
            }
            if emme > DOUBLEPI {
                emme -= DOUBLEPI;
            }

            // Argument of pericenter
            let x1 = anec.cos() - ecc;
            let rad = 1.0 - ecc * ecc;
            let x2 = sqrt(rad) * anec.sin();

            let xm = sqrt(x1 * x1 + x2 * x2);
            let x1 = x1 / xm;
            let x2 = x2 / xm;
            let sinper = x1 * xorb[1] - x2 * xorb[0];
            let cosper = x1 * xorb[0] + x2 * xorb[1];
            let argper = atan2(sinper, cosper);

            elem[0] = sma;
            elem[1] = ecc;
            elem[2] = ainc;
            elem[3] = anod;
            elem[4] = argper;
            elem[5] = emme;
        }
    }
    // CASE 2: RECA = 0 (parabolic orbit)
    else if reca == 0.0 {
        *type_ = String::from("COM");
        let ecc = 1.0;

        // Semilatus rectum and pericenter distance
        let p = el2 / GaussGravSquared;
        let q = p / 2.0;

        // True anomaly
        let cosf = p / rs - 1.0;
        let sinf = rv * p / (rs * elmod);
        let effe = atan2(sinf, cosf);

        // Argument of pericenter
        let argper = atan2(xorb[1], xorb[0]) - effe;

        elem[0] = q;
        elem[1] = ecc;
        elem[2] = ainc;
        elem[3] = anod;
        elem[4] = argper;
        elem[5] = effe;
    }
    // CASE 3: RECA < 0 (hyperbolic orbit)
    else {
        *type_ = String::from("COM");

        // Semilatus rectum, true anomaly and eccentricity
        let p = el2 / GaussGravSquared;
        let ecosf = p / rs - 1.0;
        let esinf = rv * p / (elmod * rs);
        let effe = atan2(esinf, ecosf);
        let ecc = sqrt(ecosf * ecosf + esinf * esinf);
        if (ecc - 1.0).abs() < EPS {
            reca = 0.0;
        } else {
            // Pericenter distance
            let q = p / (1.0 + ecc);

            // Argument of pericenter
            let argper = atan2(xorb[1], xorb[0]) - effe;

            elem[0] = q;
            elem[1] = ecc;
            elem[2] = ainc;
            elem[3] = anod;
            elem[4] = argper;
            elem[5] = effe;
        }
    }
}

pub fn eccentricity_control(
    asteroid_position: &Vector3<f64>,
    asteroid_velocity: &Vector3<f64>,
    peri_max: f64,
    ecc_max: f64,
) -> Option<(bool, f64, f64, f64)> {
    let ast_vel_2 = asteroid_velocity.dot(&asteroid_velocity);
    let distance_to_center = asteroid_position.norm();

    let angular_momentum = asteroid_position.cross(&asteroid_velocity);

    let angmom_norm = angular_momentum.dot(&angular_momentum);
    if angmom_norm.sqrt() == 0. {
        return None;
    }

    let lenz_prelim = asteroid_velocity.cross(&angular_momentum) * (1. / GaussGravSquared);
    let lenz_factor = asteroid_position * (1. / distance_to_center);
    let lenz_vector = lenz_prelim - lenz_factor;

    let eccentricity = lenz_vector.norm();
    let perihelie = angmom_norm / (GaussGravSquared * (1. + eccentricity));
    let energy = ast_vel_2 / 2. - GaussGravSquared / distance_to_center;
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
        let mut xv = [
            -0.62355005100316385,
            1.2114681148601605,
            0.25200059143776038,
            -1.5549845137774663E-002,
            -4.6315774892682878E-003,
            -9.3633621261339246E-004,
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
            -0.62355005100316385,
            1.0112601855976919,
            0.71310036350624140,
        );

        let asteroid_velocity = Vector3::new(
            -1.5549845137774663E-002,
            -3.8769361098376577E-003,
            -2.7014074002979964E-003,
        );

        let (accept_ecc, ecc, peri, energy) =
            eccentricity_control(&asteroid_position, &asteroid_velocity, 1e3, 2.).unwrap();

        assert!(accept_ecc);
        assert_eq!(ecc, 0.2892182648825829);
        assert_eq!(peri, 1.2904453621438048);
        assert_eq!(energy, -0.00008149473004352595);
    }
}
