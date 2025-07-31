use nalgebra::Matrix3;

use crate::{
    constants::{ArcSec, Radian, DPI, RADEG, RADSEC, T2000},
    ref_system::rotmt,
};

/// Compute the mean obliquity of the ecliptic at a given epoch (IAU 1976 model).
///
/// This function returns the mean obliquity angle ε, defined as the angle between
/// the Earth's equator and the ecliptic plane, using the standard IAU 1976 polynomial model.
/// The result is expressed in radians and is valid for dates within a few millennia
/// of the J2000 epoch.
///
/// Arguments
/// ---------
/// * `tjm`: Modified Julian Date (TT scale).
///
/// Returns
/// --------
/// * Mean obliquity of the ecliptic in radians.
///
/// Formula
/// -------
/// The obliquity ε is computed as a cubic polynomial in Julian centuries since J2000:
///
/// ```text
/// ε(t) = ε₀ + ε₁·T + ε₂·T² + ε₃·T³
/// ```
/// where:
/// - `T = (tjm - T2000) / 36525.0`,
/// - the coefficients `ε₀`, `ε₁`, `ε₂`, `ε₃` are in arcseconds and internally converted to radians.
///
/// The polynomial is evaluated using **Horner’s method** for numerical efficiency and stability:
///
/// ```text
/// ε = ((ob3 * t + ob2) * t + ob1) * t + ob0;
/// ```
///
/// # See also
/// * [`rotmt`] – constructs rotation matrices using this obliquity
/// * [`rotpn`](crate::ref_system::rotpn) – applies obliquity rotation when transforming between ecliptic and equatorial frames
pub fn obleq(tjm: f64) -> Radian {
    // Obliquity coefficients
    let ob0 = ((23.0 * 3600.0 + 26.0 * 60.0) + 21.448) * RADSEC;
    let ob1 = -46.815 * RADSEC;
    let ob2 = -0.0006 * RADSEC;
    let ob3 = 0.00181 * RADSEC;

    let t = (tjm - T2000) / 36525.0;

    ((ob3 * t + ob2) * t + ob1) * t + ob0
}

/// Compute the nutation angles in longitude and obliquity using the IAU 1980 (Wahr) model.
///
/// This function returns the nutation angles (Δψ, Δε), i.e. the periodic deviations in:
/// - ecliptic longitude (Δψ, nutation in longitude),
/// - and obliquity of the ecliptic (Δε, nutation in obliquity),
///
/// both expressed in arcseconds, using the IAU 1980 nutation theory as adopted by the IAU.
///
/// Arguments
/// ---------
/// * `tjm`: Modified Julian Date (in TT time scale).
///
/// Returns
/// --------
/// * A tuple `(Δψ, Δε)`:
///     - `Δψ`: nutation in longitude \[arcseconds\]
///     - `Δε`: nutation in obliquity \[arcseconds\]
///
/// Description
/// -----------
/// This implementation follows the IAU 1980 nutation model (Wahr), which expresses the nutation angles
/// as a sum of hundreds of periodic terms depending on five fundamental lunar and solar arguments:
/// - Mean anomaly of the Moon (l)
/// - Mean anomaly of the Sun (p)
/// - Argument of latitude of the Moon (f)
/// - Mean elongation of the Moon from the Sun (d)
/// - Longitude of the Moon's ascending node (n)
///
/// These arguments are computed as 3rd-order polynomials in time (in Julian centuries T from J2000),
/// and the nutation angles are then expressed as long trigonometric series involving various linear
/// combinations of sinusoids of these arguments.
///
/// The returned values are **in arcseconds**, as per the original IAU convention. They are typically
/// converted to radians for use in rotation matrices (via [`RADSEC`] in other modules).
///
/// # See also
/// * [`rnut80`] – uses these angles to build the nutation rotation matrix
/// * [`rotpn`](crate::ref_system::rotpn) – applies nutation when transforming between Equt and Equm systems
pub fn nutn80(tjm: f64) -> (ArcSec, ArcSec) {
    // Compute the fundamental lunar and solar arguments (in radians)
    let t1 = (tjm - T2000) / 36525.0;
    let t = t1;
    let t2 = t * t;
    let t3 = t2 * t;

    let dl = (485866.733 + 1717915922.633 * t1 + 31.310 * t2 + 0.064 * t3) * RADSEC;
    let dp = (1287099.804 + 129596581.224 * t1 - 0.577 * t2 - 0.012 * t3) * RADSEC;
    let df = (335778.877 + 1739527263.137 * t1 - 13.257 * t2 + 0.011 * t3) * RADSEC;
    let dd = (1072261.307 + 1602961601.328 * t1 - 6.891 * t2 + 0.019 * t3) * RADSEC;
    let dn = (450160.280 - 6962890.539 * t1 + 7.455 * t2 + 0.008 * t3) * RADSEC;

    let l = dl % DPI;
    let p = dp % DPI;
    let x = df % DPI * 2.0;
    let d = dd % DPI;
    let n = dn % DPI;

    // Precompute cosine and sine of fundamental arguments
    let sin_cos = |x: f64| -> (f64, f64) { (x.cos(), x.sin()) };

    let (cl, sl) = sin_cos(l);
    let (cp, sp) = sin_cos(p);
    let (cx, sx) = sin_cos(x);
    let (cd, sd) = sin_cos(d);
    let (cn, sn) = sin_cos(n);

    // Construct compound trigonometric terms used in the series expansion
    let cp2 = 2.0 * cp * cp - 1.0;

    let sp2 = 2.0 * sp * cp;
    let cd2 = 2.0 * cd * cd - 1.0;
    let sd2 = 2.0 * sd * cd;
    let cn2 = 2.0 * cn * cn - 1.0;
    let sn2 = 2.0 * sn * cn;
    let cl2 = 2.0 * cl * cl - 1.0;
    let sl2 = 2.0 * sl * cl;

    let ca = cx * cd2 + sx * sd2;
    let sa = sx * cd2 - cx * sd2;
    let cb = ca * cn - sa * sn;
    let sb = sa * cn + ca * sn;
    let cc = cb * cn - sb * sn;
    let sc = sb * cn + cb * sn;

    let cv = cx * cd2 - sx * sd2;
    let sv = sx * cd2 + cx * sd2;
    let ce = cv * cn - sv * sn;
    let se = sv * cn + cv * sn;
    let cf = ce * cn - se * sn;
    let sf = se * cn + ce * sn;

    let cg = cl * cd2 + sl * sd2;
    let sg = sl * cd2 - cl * sd2;
    let ch = cx * cn2 - sx * sn2;
    let sh = sx * cn2 + cx * sn2;
    let cj = ch * cl - sh * sl;
    let sj = sh * cl + ch * sl;

    let ck = cj * cl - sj * sl;
    let sk = sj * cl + cj * sl;
    let cm = cx * cl2 + sx * sl2;
    let sm = sx * cl2 - cx * sl2;
    let cq = cl * cd + sl * sd;
    let sq = sl * cd - cl * sd;

    let cr = 2.0 * cq * cq - 1.0;
    let sr = 2.0 * sq * cq;
    let cs = cx * cn - sx * sn;
    let ss = sx * cn + cx * sn;
    let ct = cs * cl - ss * sl;
    let st = ss * cl + cs * sl;

    let cu = cf * cl + sf * sl;
    let su = sf * cl - cf * sl;
    let cw = cp * cg - sp * sg;
    let sw = sp * cg + cp * sg;

    // Series expansion for nutation in longitude (Δψ), in 0.0001 arcseconds
    let mut dpsi =
        -(171996.0 + 174.2 * t) * sn + (2062.0 + 0.2 * t) * sn2 + 46.0 * (sm * cn + cm * sn)
            - 11.0 * sm
            - 3.0 * (sm * cn2 + cm * sn2)
            - 3.0 * (sq * cp - cq * sp)
            - 2.0 * (sb * cp2 - cb * sp2)
            + (sn * cm - cn * sm)
            - (13187.0 + 1.6 * t) * sc
            + (1426.0 - 3.4 * t) * sp
            - (517.0 - 1.2 * t) * (sc * cp + cc * sp)
            + (217.0 - 0.5 * t) * (sc * cp - cc * sp)
            + (129.0 + 0.1 * t) * sb
            + 48.0 * sr
            - 22.0 * sa
            + (17.0 - 0.1 * t) * sp2
            - 15.0 * (sp * cn + cp * sn)
            - (16.0 - 0.1 * t) * (sc * cp2 + cc * sp2)
            - 12.0 * (sn * cp - cn * sp);

    dpsi += -6.0 * (sn * cr - cn * sr) - 5.0 * (sb * cp - cb * sp)
        + 4.0 * (sr * cn + cr * sn)
        + 4.0 * (sb * cp + cb * sp)
        - 4.0 * sq
        + (sr * cp + cr * sp)
        + (sn * ca - cn * sa)
        - (sp * ca - cp * sa)
        + (sp * cn2 + cp * sn2)
        + (sn * cq - cn * sq)
        - (sp * ca + cp * sa)
        - (2274.0 + 0.2 * t) * sh
        + (712.0 + 0.1 * t) * sl
        - (386.0 + 0.4 * t) * ss
        - 301.0 * sj
        - 158.0 * sg
        + 123.0 * (sh * cl - ch * sl)
        + 63.0 * sd2
        + (63.0 + 0.1 * t) * (sl * cn + cl * sn)
        - (58.0 + 0.1 * t) * (sn * cl - cn * sl)
        - 59.0 * su
        - 51.0 * st
        - 38.0 * sf
        + 29.0 * sl2;

    dpsi += 29.0 * (sc * cl + cc * sl) - 31.0 * sk
        + 26.0 * sx
        + 21.0 * (ss * cl - cs * sl)
        + 16.0 * (sn * cg - cn * sg)
        - 13.0 * (sn * cg + cn * sg)
        - 10.0 * (se * cl - ce * sl)
        - 7.0 * (sg * cp + cg * sp)
        + 7.0 * (sh * cp + ch * sp)
        - 7.0 * (sh * cp - ch * sp)
        - 8.0 * (sf * cl + cf * sl)
        + 6.0 * (sl * cd2 + cl * sd2)
        + 6.0 * (sc * cl2 + cc * sl2)
        - 6.0 * (sn * cd2 + cn * sd2)
        - 7.0 * se
        + 6.0 * (sb * cl + cb * sl)
        - 5.0 * (sn * cd2 - cn * sd2)
        + 5.0 * (sl * cp - cl * sp)
        - 5.0 * (ss * cl2 + cs * sl2)
        - 4.0 * (sp * cd2 - cp * sd2);

    dpsi += 4.0 * (sl * cx - cl * sx) - 4.0 * sd - 3.0 * (sl * cp + cl * sp)
        + 3.0 * (sl * cx + cl * sx)
        - 3.0 * (sj * cp - cj * sp)
        - 3.0 * (su * cp - cu * sp)
        - 2.0 * (sn * cl2 - cn * sl2)
        - 3.0 * (sk * cl + ck * sl)
        - 3.0 * (sf * cp - cf * sp)
        + 2.0 * (sj * cp + cj * sp)
        - 2.0 * (sb * cl - cb * sl);

    dpsi += 2.0 * (sn * cl2 + cn * sl2) - 2.0 * (sl * cn2 + cl * sn2)
        + 2.0 * (sl * cl2 + cl * sl2)
        + 2.0 * (sh * cd + ch * sd)
        + (sn2 * cl - cn2 * sl)
        - (sg * cd2 - cg * sd2)
        + (sf * cl2 - cf * sl2)
        - 2.0 * (su * cd2 + cu * sd2)
        - (sr * cd2 - cr * sd2)
        + (sw * ch + cw * sh)
        - (sl * ce + cl * se)
        - (sf * cr - cf * sr)
        + (su * ca + cu * sa)
        + (sg * cp - cg * sp)
        + (sb * cl2 + cb * sl2)
        - (sf * cl2 + cf * sl2)
        - (st * ca - ct * sa)
        + (sc * cx + cc * sx)
        + (sj * cr + cj * sr)
        - (sg * cx + cg * sx);

    dpsi += (sp * cs + cp * ss) + (sn * cw - cn * sw)
        - (sn * cx - cn * sx)
        - (sh * cd - ch * sd)
        - (sp * cd2 + cp * sd2)
        - (sl * cv - cl * sv)
        - (ss * cp - cs * sp)
        - (sw * cn + cw * sn)
        - (sl * ca - cl * sa)
        + (sl2 * cd2 + cl2 * sd2)
        - (sf * cd2 + cf * sd2)
        + (sp * cd + cp * sd);

    // Series expansion for nutation in obliquity (Δε), in 0.0001 arcseconds
    let mut deps = (92025.0 + 8.9 * t) * cn - (895.0 - 0.5 * t) * cn2 - 24.0 * (cm * cn - sm * sn)
        + (cm * cn2 - sm * sn2)
        + (cb * cp2 + sb * sp2)
        + (5736.0 - 3.1 * t) * cc
        + (54.0 - 0.1 * t) * cp
        + (224.0 - 0.6 * t) * (cc * cp - sc * sp)
        - (95.0 - 0.3 * t) * (cc * cp + sc * sp)
        - 70.0 * cb
        + cr
        + 9.0 * (cp * cn - sp * sn)
        + 7.0 * (cc * cp2 - sc * sp2)
        + 6.0 * (cn * cp + sn * sp)
        + 3.0 * (cn * cr + sn * sr)
        + 3.0 * (cb * cp + sb * sp)
        - 2.0 * (cr * cn - sr * sn)
        - 2.0 * (cb * cp - sb * sp);

    deps += (977.0 - 0.5 * t) * ch - 7.0 * cl + 200.0 * cs + (129.0 - 0.1 * t) * cj
        - cg
        - 53.0 * (ch * cl + sh * sl)
        - 2.0 * cd2
        - 33.0 * (cl * cn - sl * sn)
        + 32.0 * (cn * cl + sn * sl)
        + 26.0 * cu
        + 27.0 * ct
        + 16.0 * cf
        - cl2
        - 12.0 * (cc * cl - sc * sl)
        + 13.0 * ck
        - cx
        - 10.0 * (cs * cl + ss * sl)
        - 8.0 * (cn * cg + sn * sg)
        + 7.0 * (cn * cg - sn * sg)
        + 5.0 * (ce * cl + se * sl)
        - 3.0 * (ch * cp - sh * sp)
        + 3.0 * (ch * cp + sh * sp)
        + 3.0 * (cf * cl - sf * sl)
        - 3.0 * (cc * cl2 - sc * sl2)
        + 3.0 * (cn * cd2 - sn * sd2)
        + 3.0 * ce
        - 3.0 * (cb * cl - sb * sl)
        + 3.0 * (cn * cd2 + sn * sd2)
        + 3.0 * (cs * cl2 - ss * sl2)
        + (cj * cp + sj * sp)
        + (cu * cp + su * sp)
        + (cn * cl2 + sn * sl2)
        + (ck * cl - sk * sl)
        + (cf * cp + sf * sp)
        - (cj * cp - sj * sp)
        + (cb * cl + sb * sl)
        - (cn * cl2 - sn * sl2)
        + (cl * cn2 - sl * sn2)
        - (ch * cd - sh * sd)
        - (cn2 * cl + sn2 * sl)
        - (cf * cl2 + sf * sl2)
        + (cu * cd2 - su * sd2)
        - (cw * ch - sw * sh)
        + (cl * ce - sl * se)
        + (cf * cr + sf * sr)
        - (cb * cl2 - sb * sl2);

    // Convert results from 0.0001 arcseconds to arcseconds
    dpsi *= 1e-4;
    deps *= 1e-4;

    (dpsi, deps)
}

/// Construct the nutation rotation matrix using the IAU 1980 nutation model.
///
/// This function returns the 3×3 rotation matrix that accounts for Earth's nutation,
/// based on the IAU 1980 theory (Wahr). It applies three successive rotations:
///
/// 1. Rotate around the X-axis by the **mean obliquity** ε (from [`obleq`]),
/// 2. Rotate around the Z-axis by the **nutation in longitude** Δψ (from [`nutn80`]),
/// 3. Rotate back around the X-axis by the **true obliquity** ε + Δε.
///
/// This yields a rotation matrix that transforms vectors from the mean equator and equinox
/// of date (Equm) to the true equator and equinox of date (Equt).
///
/// Arguments
/// ---------
/// * `tjm`: Modified Julian Date (TT scale).
///
/// Returns
/// --------
/// * A 3×3 orthonormal matrix `R` such that:
///     - `x_true = R · x_mean`
///     - where `x_mean` is a vector in the mean equatorial frame of date,
///     - and `x_true` is the same vector expressed in the true equatorial frame of date.
///
/// Notes
/// ------
/// * The obliquity angles ε and ε + Δε are in radians.
/// * The nutation angles Δψ and Δε are computed using [`nutn80`] and converted from arcseconds to radians.
/// * The rotation is expressed in the IAU 1980 nutation convention (used by OrbFit and many legacy systems).
///
/// # See also
/// * [`nutn80`] – returns the nutation angles Δψ, Δε in arcseconds
/// * [`obleq`] – computes the mean obliquity ε (radians)
/// * [`rotmt`] – builds the individual axis rotation matrices
/// * [`rotpn`](crate::ref_system::rotpn) – uses `rnut80` to transform between Equm and Equt systems
pub fn rnut80(tjm: f64) -> Matrix3<f64> {
    // Mean obliquity of the ecliptic at date (ε)
    let epsm = obleq(tjm);

    // Nutation angles in longitude (Δψ) and obliquity (Δε), in arcseconds
    let (mut dpsi, deps) = nutn80(tjm);

    // Convert nutation angles from arcseconds to radians
    dpsi *= RADSEC;
    let epst = epsm + deps * RADSEC;

    // Build individual rotation matrices:
    // R1: rotation around X by +ε (mean obliquity)
    // R2: rotation around Z by -Δψ (nutation in longitude)
    // R3: rotation around X by -ε - Δε (true obliquity)
    let r1 = rotmt(epsm, 0);
    let r2 = rotmt(-dpsi, 2);
    let r3 = rotmt(-epst, 0);

    (r1 * r2) * r3
}

/// Compute the equation of the equinoxes (nutation correction) in radians.
///
/// This term accounts for the small difference between apparent sidereal time
/// and mean sidereal time due to the nutation of Earth's rotation axis.
///
/// # Arguments
/// * `tjm` - Modified Julian Date (MJD, TT or TDB time scale)
///
/// # Returns
/// * Equation of the equinoxes in **radians**.
///
/// # Details
/// The equation of the equinoxes is computed as:
///
/// ```text
/// Eq_eq = Δψ * cos(ε)
/// ```
///
/// where:
/// * `Δψ` is the nutation in longitude (in arcseconds),
/// * `ε` is the mean obliquity of the ecliptic (in radians).
///
/// The function converts `Δψ` from arcseconds to radians using `RADSEC`.
///
/// # See also
/// * [`obleq`] – Computes the mean obliquity of the ecliptic.
/// * [`nutn80`] – Computes the 1980 IAU nutation model (Δψ and Δε).
pub fn equequ(tjm: f64) -> f64 {
    // Compute the mean obliquity of the ecliptic (ε, in radians)
    let oblm = obleq(tjm);

    // Compute nutation in longitude (Δψ) and nutation in obliquity (Δε)
    // Δψ is returned in arcseconds.
    let (dpsi, _deps) = nutn80(tjm);

    // Apply Eq_eq = Δψ * cos(ε), converting Δψ from arcseconds to radians using RADSEC
    RADSEC * dpsi * oblm.cos()
}

/// Compute the precession matrix from J2000 to the mean equator and equinox of a given epoch (IAU 1976 model).
///
/// This function constructs the 3×3 precession rotation matrix based on the IAU 1976 precession theory,
/// as formulated in the Astronomical Almanac (1987, section B18). The matrix transforms a vector
/// expressed in the J2000 mean equatorial frame into the mean equatorial frame of date.
///
/// Arguments
/// ---------
/// * `tjm`: Modified Julian Date in TT scale (epoch of transformation).
/// * `rprec`: mutable reference to a 3×3 array that will receive the resulting matrix.
///
/// Returns
/// --------
/// * The result is written in place to `rprec`, such that:
///     - `x_mean(tjm) = rprec · x_J2000`
///
/// Method
/// ------
/// The transformation is composed of three successive rotations:
/// 1. Around Z-axis by `−ζ`
/// 2. Around Y-axis by `θ`
/// 3. Around Z-axis by `−z`
///
/// where the angles are time-dependent polynomials in Julian centuries `T = (tjm - T2000) / 36525`:
///
/// ```text
/// ζ(T)     = (0.6406161 + 0.0000839·T + 0.0000050·T²) · T  [deg]
/// θ(T)     = (0.5567530 - 0.0001185·T - 0.0000116·T²) · T  [deg]
/// z(T)     = (0.6406161 + 0.0003041·T + 0.0000051·T²) · T  [deg]
/// ```
/// These angles are converted to radians internally using `RADEG`.
///
/// Remarks
/// -------
/// * This function assumes the IAU 1976 precession model, valid within a few centuries of J2000.
/// * The rotation is active (vector rotation), and the resulting matrix is orthonormal.
/// * Equivalent to the OrbFit `prec` Fortran routine used in reference frame transitions.
///
/// # See also
/// * [`rotmt`] – constructs the rotation matrices used here
/// * [`rotpn`](crate::ref_system::rotpn) – uses `prec` when converting between epochs `"OFDATE"` and `"J2000"`
pub fn prec(tjm: f64) -> Matrix3<f64> {
    // Precession polynomial coefficients (in radians)
    let zed = 0.6406161 * RADEG;
    let zd = 0.6406161 * RADEG;
    let thd = 0.5567530 * RADEG;

    let zedd = 0.0000839 * RADEG;
    let zdd = 0.0003041 * RADEG;
    let thdd = -0.0001185 * RADEG;

    let zeddd = 0.0000050 * RADEG;
    let zddd = 0.0000051 * RADEG;
    let thddd = -0.0000116 * RADEG;

    // Compute Julian centuries since J2000
    let t = (tjm - T2000) / 36525.0;

    // Compute precession angles (in radians)
    let zeta = ((zeddd * t + zedd) * t + zed) * t;
    let z = ((zddd * t + zdd) * t + zd) * t;
    let theta = ((thddd * t + thdd) * t + thd) * t;

    // Construct the three rotation matrices:
    // R1 = rotation around Z by −ζ
    // R2 = rotation around Y by θ
    // R3 = rotation around Z by −z
    let r1 = rotmt(-zeta, 2); // Z-axis
    let r2 = rotmt(theta, 1); // Y-axis
    let r3 = rotmt(-z, 2); // Z-axis

    // Compute the combined rotation matrix R3 . (R2 . R1)
    (r1 * r2) * r3
}

#[cfg(test)]
mod test_earth_oritentation {
    use super::*;

    #[test]
    fn test_obliquity() {
        let obl = obleq(T2000);
        assert_eq!(obl, 0.40909280422232897)
    }

    #[test]
    fn test_nutn80() {
        let (dpsi, deps) = nutn80(T2000);
        assert_eq!(dpsi, -13.923385169502602);
        assert_eq!(deps, -5.773808263765919);
    }

    #[test]
    fn test_rnut80() {
        let rnut = rnut80(T2000);
        let ref_rnut = [
            [
                0.9999999977217079,
                6.19323109890795e-5,
                2.6850942970991024e-5,
            ],
            [
                -6.193306258211379e-5,
                0.9999999976903892,
                2.799138089948361e-5,
            ],
            [
                -2.6849209338068913e-5,
                -2.7993043796858963e-5,
                0.9999999992477547,
            ],
        ];
        assert_eq!(rnut, ref_rnut.into());
    }

    mod tests_equequ {
        use super::*;
        use approx::assert_relative_eq;

        /// Helper: convert radians to arcseconds
        fn rad_to_arcsec(x: f64) -> f64 {
            x / RADSEC
        }

        #[test]
        fn test_equequ_at_j2000() {
            // MJD for J2000 epoch
            let tjm = T2000;
            let eqeq = equequ(tjm);

            // Expected value at J2000:
            // Δψ ≈ -13.923385 arcsec and ε ≈ 23.44°, so Δψ * cos(ε) in radians
            let expected_rad = RADSEC * (-13.923385169502602) * (obleq(tjm).cos());

            // Verify that the computed value matches the expected reference
            assert_relative_eq!(eqeq, expected_rad, epsilon = 1e-12);

            // The value is expected to be on the order of a few tens of arcseconds
            let arcsec = rad_to_arcsec(eqeq.abs());
            assert!(arcsec < 30.0, "Equation of equinoxes should be small");
        }

        #[test]
        fn test_equequ_changes_with_time() {
            let t0 = 51544.5;
            let t1 = 60000.0; // about 23 years later

            let eq0 = equequ(t0);
            let eq1 = equequ(t1);

            // The value must change over time due to nutation/precession evolution
            assert!((eq1 - eq0).abs() > 1e-7);
        }

        #[test]
        fn test_equequ_is_small() {
            let t = 58000.0;
            let val = equequ(t);

            // The equation of equinoxes must remain small (below 0.001 rad ≈ 206 arcsec)
            assert!(val.abs() < 1e-3);
        }
    }
}
