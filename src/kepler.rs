use super::constants::GAUSS_GRAV;
use super::constants::DPI;
use super::orb_elem::eccentricity_control;
use core::fmt;
use nalgebra::Vector3;
use std::f64::consts::PI;
use std::f64::EPSILON;

fn s_funct(psi: f64, alpha: f64) -> (f64, f64, f64, f64) {
    const JMAX: usize = 70;
    const HALFMAX: usize = 30;
    const BETACONTR: f64 = 100.0;

    let epsilon = f64::EPSILON;
    let contr = 100.0 * epsilon;
    let overfl = 1.0 / epsilon;

    let beta = alpha * psi.powi(2);
    let mut s0: f64;
    let mut s1: f64;
    let mut s2: f64;
    let mut s3: f64;

    if beta.abs() < BETACONTR {
        // Calcul par développement en série
        let mut term2 = psi.powi(2) / 2.0;
        let mut term3 = term2 * psi / 3.0;
        s2 = term2;
        s3 = term3;

        for j in 1..=JMAX {
            term2 *= beta / ((2.0 * j as f64 + 1.0) * (2.0 * j as f64 + 2.0));
            s2 += term2;
            if term2.abs() < contr || term2.abs() > overfl {
                break;
            }
        }

        for j in 1..=JMAX {
            term3 *= beta / ((2.0 * j as f64 + 2.0) * (2.0 * j as f64 + 3.0));
            s3 += term3;
            if term3.abs() < contr || term3.abs() > overfl {
                break;
            }
        }

        // Calcul de s1 et s0
        s1 = psi + alpha * s3;
        s0 = 1.0 + alpha * s2;
    } else {
        // Réduction de psi pour éviter les erreurs numériques
        let mut psi2 = psi;
        let mut nhalf = 0;

        for _ in 0..HALFMAX {
            psi2 *= 0.5;
            nhalf += 1;
            let beta_half = alpha * psi2.powi(2);
            if beta_half.abs() < BETACONTR {
                break;
            }
        }

        // Calcul initial de s0 et s1 par séries
        let mut term0 = 1.0;
        let mut term1 = psi2;
        s0 = 1.0;
        s1 = psi2;

        for j in 1..=JMAX {
            term0 *= beta / ((2 * j - 1) as f64 * (2 * j) as f64);
            s0 += term0;
            if term0.abs() < contr || term0.abs() > overfl {
                break;
            }
        }

        for j in 1..=JMAX {
            term1 *= beta / ((2 * j) as f64 * (2 * j + 1) as f64);
            s1 += term1;
            if term1.abs() < contr || term1.abs() > overfl {
                break;
            }
        }

        // Formules de duplication pour retrouver la valeur originale de psi
        for _ in 0..nhalf {
            let s02 = 2.0 * s0.powi(2) - 1.0;
            let s12 = 2.0 * s0 * s1;
            s0 = s02;
            s1 = s12;
        }

        // Calcul final de s2 et s3
        s3 = (s1 - psi) / alpha;
        s2 = (s0 - 1.0) / alpha;
    }

    (s0, s1, s2, s3)
}

/// Retourne la valeur principale d'un angle en radians dans [0, 2π].
fn principal_angle(a: f64) -> f64 {
    a.rem_euclid(DPI) // rem_euclid assure un résultat dans [0, 2π]
}

/// Retourne la différence principale entre deux angles dans [-π, π].
fn angle_diff(a: f64, b: f64) -> f64 {
    // Normalisation des angles entre [0, 2π]
    let a = principal_angle(a);
    let b = principal_angle(b);

    // Calcul de la différence
    let mut diff = a - b;

    // Ajustement pour être dans [-π, π]
    if diff > PI {
        diff -= DPI;
    } else if diff < -PI {
        diff += DPI;
    }

    diff
}

fn prelim_kepuni(
    dt: f64,
    r0: f64,
    sig0: f64,
    mu: f64,
    alpha: f64,
    e0: f64,
    contr: f64,
) -> Option<(f64, f64)> {
    const ITX: usize = 20;

    // Initialisation de psi
    let psi0;

    if alpha < 0.0 {
        // Cas elliptique
        let a0 = -mu / alpha;
        let enne = (-alpha.powi(3)).sqrt() / mu;
        let (u0, u) = if e0 < contr {
            (0.0, enne * dt)
        } else {
            let cosu0 = (1.0 - r0 / a0) / e0;
            let mut u0 = if cosu0.abs() <= 1.0 {
                cosu0.acos()
            } else if cosu0 >= 1.0 {
                0.0
            } else {
                PI
            };

            if sig0 < 0.0 {
                u0 = -u0;
            }

            u0 = principal_angle(u0);
            let ell0 = principal_angle(u0 - e0 * u0.sin());
            let mut u = PI;
            let ell = principal_angle(ell0 + enne * dt);

            for _ in 0..ITX {
                let du = -(u - e0 * u.sin() - ell) / (1.0 - e0 * u.cos());
                u += du;
                if du.abs() < contr * 1e3 {
                    break;
                }
            }
            (u0, u)
        };

        // Conversion en psi
        psi0 = angle_diff(u, u0) / (-alpha).sqrt();
    } else if alpha > 0.0 {
        // Cas hyperbolique
        let a0 = -mu / alpha;
        let enne = alpha.powi(3).sqrt() / mu;
        let coshf0 = (1.0 - r0 / a0) / e0;
        let mut f0 = if coshf0 > 1.0 {
            (coshf0 + (coshf0.powi(2) - 1.0).sqrt()).ln()
        } else {
            0.0
        };

        if sig0 < 0.0 {
            f0 = -f0;
        }

        let ell0 = e0 * f0.sinh() - f0;
        let mut f: f64 = 0.0;
        let ell = ell0 + enne * dt;

        for _ in 0..ITX {
            if f.abs() < 15.0 {
                let df = -(e0 * f.sinh() - f - ell) / (e0 * f.cosh() - 1.0);
                let ff = f + df;
                f = if f * ff < 0.0 { f / 2.0 } else { ff };
            } else {
                f /= 2.0;
            }
            if f.abs() < contr * 1e3 {
                break;
            }
        }

        // Conversion en psi
        psi0 = (f - f0) / alpha.sqrt();
    } else {
        return None; // Cas non supporté
    }

    Some((psi0, alpha))
}

/// Résout l'équation universelle de Kepler en utilisant une méthode de Newton.
/// Gère les cas elliptiques (`alpha < 0`) et hyperboliques (`alpha > 0`).
fn solve_kepuni(
    dt: f64,
    r0: f64,
    sig0: f64,
    mu: f64,
    alpha: f64,
    e0: f64,
    convergency: Option<f64>,
) -> Option<(f64, f64, f64, f64, f64)> {
    const JMAX: usize = 100;
    let contr = convergency.unwrap_or(100.0 * EPSILON);

    let Some((mut psi, alpha)) = prelim_kepuni(dt, r0, sig0, mu, alpha, e0, contr) else {
        return None;
    };

    // Méthode de Newton pour résoudre l'équation universelle de Kepler
    for _ in 0..JMAX {
        let (s0, s1, s2, s3) = s_funct(psi, alpha);

        let fun = r0 * s1 + sig0 * s2 + mu * s3 - dt;
        let funp = r0 * s0 + sig0 * s1 + mu * s2;

        let dpsi = -fun / funp;

        if s3.abs() > 1e-2 / EPSILON {
            return None; // Problème de convergence
        }

        let psi1 = psi + dpsi;
        psi = if psi1 * psi < 0.0 { psi / 2.0 } else { psi1 };

        if dpsi.abs() < contr || dpsi.abs() < contr * 10.0 * psi.abs() {
            return Some((psi, s0, s1, s2, s3));
        }
    }

    None // Convergence non atteinte
}

/// VelocityCorrectionError Error is used in case the velocity correction cannot be ended.
#[derive(Debug, Clone)]
pub struct VelocityCorrectionError;

impl fmt::Display for VelocityCorrectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "The eccentricity control detected that the asteroid as no angular momentum"
        )
    }
}

pub fn velocity_correction(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
) -> Result<(Vector3<f64>, f64, f64), VelocityCorrectionError> {
    let mu = GAUSS_GRAV.powi(2);
    let sig0 = x2.dot(&v2);
    let r2 = x2.norm();

    let Some((_, ecc, _, energy)) = eccentricity_control(&x2, &v2, peri_max, ecc_max) else {
        return Err(VelocityCorrectionError);
    };

    let eps = 1e3 * EPSILON;
    let alpha = 2. * energy;

    let Some((_, _, _, s2, s3)) = solve_kepuni(dt, r2, sig0, mu, alpha, ecc, Some(eps)) else {
        return Err(VelocityCorrectionError);
    };

    // Calcul des coefficients f et g de Lagrange
    let f = 1.0 - (mu * s2) / r2;
    let g = dt - (mu * s3);

    // Calcul de la vitesse améliorée
    let v2 = (x1 - f * x2) / g;

    Ok((v2, f, g))
}

#[cfg(test)]
mod kepler_test {

    use super::*;

    #[test]
    fn test_s_funct() {
        let psi = -15.279808141051223;
        let alpha = -1.6298946008705195e-4;

        let (s0, s1, s2, s3) = s_funct(psi, alpha);

        assert_eq!(s0, 0.9810334785583247);
        assert_eq!(s1, -15.183083836892674);
        assert_eq!(s2, 116.3665517484714);
        assert_eq!(s3, -593.4390119881925);
    }

    #[test]
    fn test_prelim_kepuni() {
        let epsilon = f64::EPSILON;
        let contr = 100.0 * epsilon;

        let dt = -20.765849999996135;
        let r0 = 1.3803870211345761;
        let sig0 = 3.7013544840038748E-003;
        let mu = 2.9591220828559115E-004;
        let alpha = -1.6421583777711407E-004;
        let e0 = 0.28359959913734450;

        let (psi, alpha) = prelim_kepuni(dt, r0, sig0, mu, alpha, e0, contr).unwrap();

        assert_eq!(psi, -15.327414893041848);
        assert_eq!(alpha, -0.00016421583777711407);

        let alpha = 1.6421583777711407E-004;
        let (psi, alpha) = prelim_kepuni(dt, r0, sig0, mu, alpha, e0, contr).unwrap();

        assert_eq!(psi, -73.1875935362658);
        assert_eq!(alpha, 0.00016421583777711407);

        let res_prelim = prelim_kepuni(dt, r0, sig0, mu, 0., e0, contr);
        assert!(res_prelim.is_none());
    }

    #[test]
    fn test_solve_kepuni() {
        let dt = -20.765849999996135;
        let r0 = 1.3803870211345761;
        let sig0 = 3.7013544840038748E-003;
        let mu = 2.9591220828559115E-004;
        let alpha = -1.6421583777711407E-004;
        let e0 = 0.28359959913734450;

        let (psi, s0, s1, s2, s3) = solve_kepuni(dt, r0, sig0, mu, alpha, e0, None).unwrap();

        assert_eq!(psi, -15.327414893041839);
        assert_eq!(s0, 0.9807723505583343);
        assert_eq!(s1, -15.229051668919967);
        assert_eq!(s2, 117.0876676813769);
        assert_eq!(s3, -598.9874390519309);

        let alpha = 1.6421583777711407E-004;
        let (psi, s0, s1, s2, s3) = solve_kepuni(dt, r0, sig0, mu, alpha, e0, None).unwrap();

        assert_eq!(psi, -15.1324122746124);
        assert_eq!(s0, 1.0188608766146905);
        assert_eq!(s1, -15.227430038021337);
        assert_eq!(s2, 114.854187452308);
        assert_eq!(s3, -578.615100072754);
    }

    #[test]
    fn test_velocity_correction() {
        let x1 = Vector3::new(
            -0.84356112612968326,
            0.93728832737077283,
            0.65918390102977664,
        );

        let x2 = Vector3::new(
            -0.62312162291738404,
            1.0076797884556383,
            0.70812568798442455,
        );

        let v2 = Vector3::new(
            -1.5524310368624056E-002,
            -3.9841041766040678E-003,
            -2.7640154361637183E-003,
        );
        let dt = 14.731970000000729;

        let (v2, f, g) = velocity_correction(&x1, &x2, &v2, dt, 1., 1.).unwrap();

        assert_eq!(f, 0.98816487709729062);
        assert_eq!(g, 14.674676076120734);
        assert_eq!(
            v2.as_slice(),
            [
                -0.015524310248562921,
                -0.003984104769239458,
                -0.0027640155187336176
            ]
        )
    }
}
