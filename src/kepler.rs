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
}
