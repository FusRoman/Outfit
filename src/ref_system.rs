use std::f64::consts::PI;

const EPS: f64 = 1e-6;
const T2000: f64 = 51544.5; // Époque J2000 en jours juliens modifiés
const RADEG: f64 = std::f64::consts::PI / 180.0; // Conversion degrés -> radians
const RADSEC: f64 = std::f64::consts::PI / 648000.0; // Conversion d'arcsecondes en radians

fn rotpn(
    rot: &mut [[f64; 3]; 3],
    rsys1: &str,
    epoch1: &str,
    date1: f64,
    rsys2: &str,
    epoch2: &str,
    date2: f64,
) {
    let mut r: [[f64; 3]; 3] = [[0.0; 3]; 3];
    let mut date: f64;
    let mut epoch = epoch1.to_string();
    let mut rsys = rsys1.to_string();
    let mut epdif: bool;

    // Vérification des systèmes de référence
    if !chkref(rsys1, epoch1) {
        eprintln!(
            "ERROR: unsupported {} reference system: RSYS = {} EPOCH = {}",
            "starting", rsys1, epoch1
        );
        panic!("**** rotpn: abnormal end ****");
    }

    if !chkref(rsys2, epoch2) {
        eprintln!(
            "ERROR: unsupported {} reference system: RSYS = {} EPOCH = {}",
            "final", rsys1, epoch1
        );
        panic!("**** rotpn: abnormal end ****");
    }

    // Déterminer la date selon l'époque
    if epoch == "J2000" {
        date = T2000;
    } else if epoch == "OFDATE" {
        date = date1;
    } else {
        panic!("**** rotpn: internal error (01) ****");
    }

    // Initialisation de la matrice de rotation (matrice identité)
    *rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];

    // Boucle pour la transformation
    let mut nit = 0;

    loop {
        // Vérifier si les époques sont les mêmes
        if epoch == epoch2 {
            if epoch == "J2000" {
                epdif = false;
            } else if epoch == "OFDATE" {
                epdif = (date - date1).abs() > EPS;
            } else {
                panic!("**** rotpn: internal error (02) ****");
            }
        } else {
            epdif = true;
        }

        if epdif {
            // Transformation passant par le système équatorial J2000
            if epoch != "J2000" {
                if rsys == "ECLM" {
                    // Transformation écliptique --> équatoriale
                    let obl = obleq(date);
                    let r = rotmt(-obl, 1);
                    matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUT" {
                    // Transformation équateur vrai --> équateur moyen
                    let mut r = rnut80(date);
                    trsp3(&mut r);
                    matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUM" {
                    // Transformation vers J2000 (précession)
                    prec(date, &mut r);
                    trsp3(&mut r);
                    matmul(&r, rot);
                    epoch = "J2000".to_string();
                    date = T2000;
                } else {
                    panic!("**** rotpn: internal error (03) ****");
                }
            } else {
                if rsys == "ECLM" {
                    // Transformation écliptique --> équatoriale
                    let obl = obleq(T2000);
                    let r = rotmt(-obl, 1);
                    matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUT" {
                    // Transformation équateur vrai --> équateur moyen
                    let mut r = rnut80(T2000);
                    trsp3(&mut r);
                    matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUM" {
                    if epoch2 == "OFDATE" {
                        prec(date2, &mut r);
                        matmul(&r, rot);
                        epoch = epoch2.to_string();
                        date = date2;
                    } else {
                        panic!("**** rotpn: internal error (04) ****");
                    }
                } else {
                    panic!("**** rotpn: internal error (05) ****");
                }
            }
        } else {
            if rsys == rsys2 {
                return;
            }

            // Transformation pour les systèmes de référence à la même époque
            if rsys == "EQUT" {
                // Transformation équateur vrai --> équateur moyen
                let mut r = rnut80(date);
                trsp3(&mut r);
                matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "ECLM" {
                // Transformation écliptique --> équatoriale
                let obl = obleq(date);
                let r = rotmt(-obl, 1);
                matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "EQUM" {
                if rsys2 == "EQUT" {
                    // Transformation équateur moyen --> équateur vrai
                    let mut r = rnut80(date);
                    matmul(&r, rot);
                    rsys = "EQUT".to_string();
                } else if rsys2 == "ECLM" {
                    // Transformation équatoriale --> écliptique
                    let obl = obleq(date);
                    let r = rotmt(obl, 1);
                    matmul(&r, rot);
                    rsys = "ECLM".to_string();
                } else {
                    panic!("**** rotpn: internal error (06) ****");
                }
            } else {
                panic!("**** rotpn: internal error (07) ****");
            }
        }

        nit += 1;
        if nit > 20 {
            panic!("**** rotpn: internal error (08) ****");
        }
    }
}

// Remarques: ces fonctions doivent être implémentées selon le même principe que dans Fortran
fn chkref(rsys: &str, epoch: &str) -> bool {
    // Flag d'erreur, commence à true
    let mut error = true;

    // Vérification du système de référence
    if rsys == "EQUM" || rsys == "EQUT" || rsys == "ECLM" {
        // Vérification de l'époque
        if epoch == "J2000" || epoch == "OFDATE" {
            error = false; // Pas d'erreur
        }
    }

    error
}

fn obleq(tjm: f64) -> f64 {
    // Coefficients d'obliquité
    let ob0 = ((23.0 * 3600.0 + 26.0 * 60.0) + 21.448) * RADSEC;
    let ob1 = -46.815 * RADSEC;
    let ob2 = -0.0006 * RADSEC;
    let ob3 = 0.00181 * RADSEC;

    let t = (tjm - T2000) / 36525.0;

    ((ob3 * t + ob2) * t + ob1) * t + ob0
}

fn rotmt(alpha: f64, k: usize) -> [[f64; 3]; 3] {
    if k < 1 || k > 3 {
        panic!("**** ROTMT: k = ??? ****");
    }

    let cosa = alpha.cos();
    let sina = alpha.sin();
    let mut r = [[0.0; 3]; 3];

    let i1 = (k - 1) % 3;
    let i2 = (i1 + 1) % 3;
    let i3 = (i2 + 1) % 3;

    r[i1][i1] = 1.0;
    r[i1][i2] = 0.0;
    r[i1][i3] = 0.0;
    r[i2][i1] = 0.0;
    r[i2][i2] = cosa;
    r[i2][i3] = sina;
    r[i3][i1] = 0.0;
    r[i3][i2] = -sina;
    r[i3][i3] = cosa;

    r
}


// arcsecond to radian conversion factor
const RS: f64 = 4.84813681109536e-6;
const P2: f64 = 2.0 * PI;
fn nutn80(tjm: f64) -> (f64, f64) {

    // computation of secular julian date
    let t1 = (tjm - T2000) / 36525.0;
    let t = t1 as f32;
    let t2 = t * t;
    let t3 = t2 * t;

    // Arguments fondamentaux (IAU 1980)
    // evaluation of the earth nutation
    // Nutation: tidal gravitational forcing on the solid Earth
    // each argument is in arcsecond converted in radians using the RS factor.

    // dl is the moon longitude representing moon mean anomaly (the moon position along its orbit)
    let dl = (485866.733 + 1717915922.633 * t1 + 31.310 * t2 as f64 + 0.064 * t3 as f64) * RS;

    // dp is the sun longitude representing sun mean anomaly (the sun position along its orbit)
    let dp = (1287099.804 + 129596581.224 * t1 - 0.577 * t2 as f64 - 0.012 * t3 as f64) * RS;

    // df is the moon latitude argument used to model moon orbit variation
    let df = (335778.877 + 1739527263.137 * t1 - 13.257 * t2 as f64 + 0.011 * t3 as f64) * RS;

    // dd is the mean distance earth-moon. 
    // it is computed as the longitude difference between sun and moon
    let dd = (1072261.307 + 1602961601.328 * t1 - 6.891 * t2 as f64 + 0.019 * t3 as f64) * RS;

    // dn is the moon ascending node, point where the moon cross the ecliptic plane from south to north.
    let dn = (450160.280 - 6962890.539 * t1 + 7.455 * t2 as f64 + 0.008 * t3 as f64) * RS;

    // Réduction des valeurs dans [0, 2π]
    let l = (dl % P2) as f32;
    let p = (dp % P2) as f32;
    let x = ((df % P2) * 2.0) as f32;
    let d = (dd % P2) as f32;
    let n = (dn % P2) as f32;

    // Calcul des valeurs trigonométriques
    let (cl, sl) = (l.cos(), l.sin());
    let (cp, sp) = (p.cos(), p.sin());
    let (cx, sx) = (x.cos(), x.sin());
    let (cd, sd) = (d.cos(), d.sin());
    let (cn, sn) = (n.cos(), n.sin());

    // Calculs auxiliaires (puissances de sinus/cosinus)
    let cp2 = 2.0 * cp * cp - 1.0;
    let sp2 = 2.0 * sp * cp;
    let cd2 = 2.0 * cd * cd - 1.0;
    let sd2 = 2.0 * sd * cd;
    let cn2 = 2.0 * cn * cn - 1.0;
    let sn2 = 2.0 * sn * cn;

    let cl2 = 2.0 * cl * cl - 1.0;
    let sl2 = 2.0 * sl * cl;

    // Séries de nutation en longitude (DPSI)
    // The series dpsi and deps correspond to a sum of 106 nutation coefficient.
    // The coefficient come from the earth nutation table established by the IAU 1980 according to the Whar theory.

    // dpsi is the ecliptic longitude of the nutation, unit = µarcsecond
    let mut dpsi = -(171996.0 + 174.2 * t as f64) * sn as f64
        + (2062.0 + 0.2 * t as f64) * sn2 as f64
        + 46.0 * (sn as f64 * cp as f64 - cn as f64 * sp as f64)
        - 11.0 * (sn as f64 * cn as f64 - cn as f64 * sn as f64)
        - (13187.0 + 1.6 * t as f64) * (sl as f64 * cn as f64 - cl as f64 * sn as f64)
        + (1426.0 - 3.4 * t as f64) * sp as f64
        - (2274.0 + 0.2 * t as f64) * (cx as f64 * cn2 as f64 - sx as f64 * sn2 as f64);

    // Séries de nutation en obliquité (DEPS)

    // deps is the ecliptic obliquity of the nutation, unit = µarcsecond
    let mut deps = (92025.0 + 8.9 * t as f64) * cn as f64 - (895.0 - 0.5 * t as f64) * cn2 as f64
        + (5736.0 - 3.1 * t as f64) * (cl as f64 * cn as f64 - sl as f64 * sn as f64)
        + (224.0 - 0.6 * t as f64) * (cx as f64 * cp as f64 - sx as f64 * sp as f64)
        - (977.0 - 0.5 * t as f64) * (cl as f64 * cn as f64 + sl as f64 * sn as f64)
        - (129.0 - 0.1 * t as f64) * (cx as f64 * cn as f64 - sx as f64 * sn as f64);

    // Conversion en arcsecondes
    dpsi *= 1e-4;
    deps *= 1e-4;

    (dpsi, deps)
}

fn rnut80(tjm: f64) -> [[f64; 3]; 3] {
    let epsm = obleq(tjm);
    let (mut dpsi, mut deps) = nutn80(tjm);
    dpsi *= RADSEC;
    let epst = epsm + deps * RADSEC;

    let r1 = rotmt(epsm, 1);
    let r2 = rotmt(-dpsi, 3);
    let r3 = rotmt(-epst, 1);

    let mut rp = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rp[i][j] = r2[i][0] * r1[0][j] + r2[i][1] * r1[1][j] + r2[i][2] * r1[2][j];
        }
    }

    let mut rnut = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rnut[i][j] = r3[i][0] * rp[0][j] + r3[i][1] * rp[1][j] + r3[i][2] * rp[2][j];
        }
    }

    rnut
}

fn trsp3(matrix: &mut [[f64; 3]; 3]) {
    for i in 0..2 {
        for j in (i + 1)..3 {
            let temp = matrix[i][j];
            matrix[i][j] = matrix[j][i];
            matrix[j][i] = temp;
        }
    }
}

fn prec(tjm: f64, rprec: &mut [[f64; 3]; 3]) {
    // Déclaration et initialisation des matrices de rotation
    let mut r1 = [[0.0; 3]; 3];
    let mut r2 = [[0.0; 3]; 3];
    let mut r3 = [[0.0; 3]; 3];
    let mut rp = [[0.0; 3]; 3];

    // Constantes de précession
    let zed = 0.6406161 * RADEG;
    let zd = 0.6406161 * RADEG;
    let thd = 0.5567530 * RADEG;

    let zedd = 0.0000839 * RADEG;
    let zdd = 0.0003041 * RADEG;
    let thdd = -0.0001185 * RADEG;

    let zeddd = 0.0000050 * RADEG;
    let zddd = 0.0000051 * RADEG;
    let thddd = -0.0000116 * RADEG;

    // Calcul des paramètres fondamentaux
    let t = (tjm - T2000) / 36525.0;
    let zeta = ((zeddd * t + zedd) * t + zed) * t;
    let z = ((zddd * t + zdd) * t + zd) * t;
    let theta = ((thddd * t + thdd) * t + thd) * t;

    // Calcul des matrices de rotation
    let r1 = rotmt(-zeta, 3);
    let r2 = rotmt(theta, 2);
    let r3 = rotmt(-z, 3);

    // Multiplication des matrices r1 et r2 pour obtenir rp
    for i in 0..3 {
        for j in 0..3 {
            rp[i][j] = r2[i][0] * r1[0][j] + r2[i][1] * r1[1][j] + r2[i][2] * r1[2][j];
        }
    }

    // Multiplication de rp par r3 pour obtenir rprec
    for i in 0..3 {
        for j in 0..3 {
            rprec[i][j] = r3[i][0] * rp[0][j] + r3[i][1] * rp[1][j] + r3[i][2] * rp[2][j];
        }
    }
}

fn matmul(a: &[[f64; 3]; 3], b: &mut [[f64; 3]; 3]) {
    // Implémentez la multiplication matricielle
    for i in 0..3 {
        for j in 0..3 {
            b[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
        }
    }
}
