use nalgebra::Vector3;

use crate::constants::{ArcSec, Radian, VLIGHT_AU};

use super::constants::{DPI, EPS, RADEG, RADSEC, T2000};

// enum RefEpoch {
//     J2000,
//     EPOCH(f64),
// }

// enum RefSystem {
//     // Equatorial Mean, equatorial coordinates based on equator and mean equinox
//     // at a given epoch (J2000 for instance)
//     // (corrected for precession but not for nutation)
//     EQUM(RefEpoch),
//     // Equatorial True (same as EQUM but corrected for precession and nutation)
//     EQUT(RefEpoch),
//     // Ecliptic mean, ecliptic coordinates based on ecliptic and mean equinox
//     // at a given epoch (J2000 for instance)
//     ECLM(RefEpoch),
// }

/// Compute the rotation matrix between two celestial reference systems and epochs.
///
/// This function builds a composite rotation matrix that transforms coordinates
/// from a source reference system and epoch to a target system and epoch. The supported
/// systems are:
/// - `"EQUM"`: equatorial mean (precession only)
/// - `"EQUT"`: equatorial true (precession + nutation)
/// - `"ECLM"`: ecliptic mean (precession + obliquity)
///
/// Epochs can be either:
/// - `"J2000"`: standard epoch (fixed at MJD 51544.5)
/// - `"OFDATE"`: user-specified epoch (e.g., time of observation)
///
/// If the two systems differ in reference frame and/or epoch, the rotation is assembled
/// by chaining a sequence of elementary transformations: obliquity rotation, nutation,
/// and/or precession, passing if necessary through the equatorial mean J2000 frame.
///
/// Arguments
/// ---------
/// * `rot`: mutable 3×3 array to store the resulting rotation matrix.
/// * `rsys1`: name of the source reference system (`"EQUM"`, `"EQUT"`, `"ECLM"`).
/// * `epoch1`: epoch of the source system (`"J2000"` or `"OFDATE"`).
/// * `date1`: time of the source system in MJD TT (only used if `epoch1 == "OFDATE"`).
/// * `rsys2`: name of the target reference system.
/// * `epoch2`: epoch of the target system.
/// * `date2`: time of the target system in MJD TT (only used if `epoch2 == "OFDATE"`).
///
/// Output
/// -------
/// * `rot`: the rotation matrix such that `x₂ = rot · x₁`, where `x₁` is a vector
///   in the source system and `x₂` the same vector expressed in the target system.
///
/// Remarks
/// -------
/// * The rotation is built iteratively, updating the internal state until the final system/epoch is reached.
/// * Each step composes the transformation matrix `rot` from right to left.
/// * Precession uses the IAU 1976 model (`prec`), nutation uses IAU 1980 (`rnut80`), and obliquity from `obleq`.
/// * The transformation path may pass through `"EQUM", "J2000"` as an intermediate canonical frame.
///
/// Panics
/// -------
/// * If an unsupported reference system or epoch is passed.
/// * If more than 20 transformation steps are required (safety cap).
///
/// # See also
/// * [`prec`] – IAU 1976 precession matrix
/// * [`rnut80`] – IAU 1980 nutation model
/// * [`rotmt`] – rotation matrix around X/Y/Z axes
/// * [`obleq`] – mean obliquity of the ecliptic (in radians)
pub fn rotpn(
    rot: &mut [[f64; 3]; 3],
    rsys1: &str,
    epoch1: &str,
    date1: f64,
    rsys2: &str,
    epoch2: &str,
    date2: f64,
) {
    if !chkref(rsys1, epoch1) {
        panic!("ERROR: Unsupported starting reference system {rsys1} {epoch1}");
    }
    if !chkref(rsys2, epoch2) {
        panic!("ERROR: Unsupported final reference system {rsys2} {epoch2}");
    }

    let mut rsys = rsys1.to_string();
    let mut epoch = epoch1.to_string();
    let mut date = if epoch == "J2000" { T2000 } else { date1 };

    *rot = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];

    let mut nit = 0;

    loop {
        let epdif = if epoch == epoch2 {
            if epoch == "J2000" {
                false
            } else {
                (date - date1).abs() > EPS
            }
        } else {
            true
        };

        let mut r = [[0.0; 3]; 3];

        if epdif {
            if epoch != "J2000" {
                if rsys == "ECLM" {
                    let obl = obleq(date);
                    let r = rotmt(-obl, 1);
                    *rot = matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUT" {
                    let mut r = rnut80(date);
                    trsp3(&mut r);
                    *rot = matmul(&r, rot);
                    rsys = "EQUM".to_string();
                } else if rsys == "EQUM" {
                    prec(date, &mut r);
                    trsp3(&mut r);
                    *rot = matmul(&r, rot);
                    epoch = "J2000".to_string();
                    date = T2000;
                } else {
                    panic!("ERROR: Internal error (03)");
                }
            } else if rsys == "ECLM" {
                let obl = obleq(T2000);
                let r = rotmt(-obl, 1);
                *rot = matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "EQUT" {
                let mut r = rnut80(T2000);
                trsp3(&mut r);
                *rot = matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "EQUM" {
                if epoch2 == "OFDATE" {
                    prec(date2, &mut r);
                    *rot = matmul(&r, rot);
                    epoch = epoch2.to_string();
                    date = date2;
                } else {
                    panic!("ERROR: Internal error (04)");
                }
            } else {
                panic!("ERROR: Internal error (05)");
            }
        } else {
            if rsys == rsys2 {
                return;
            }

            if rsys == "EQUT" {
                let mut r = rnut80(date);
                trsp3(&mut r);
                *rot = matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "ECLM" {
                let obl = obleq(date);
                let r = rotmt(-obl, 0);
                *rot = matmul(&r, rot);
                rsys = "EQUM".to_string();
            } else if rsys == "EQUM" {
                if rsys2 == "EQUT" {
                    let r = rnut80(date);
                    *rot = matmul(&r, rot);
                    rsys = "EQUT".to_string();
                } else if rsys2 == "ECLM" {
                    let obl = obleq(date);
                    let r = rotmt(obl, 0);
                    *rot = matmul(&r, rot);
                    rsys = "ECLM".to_string();
                } else {
                    panic!("ERROR: Internal error (06)");
                }
            } else {
                panic!("ERROR: Internal error (07)");
            }
        }

        nit += 1;
        if nit > 20 {
            panic!("ERROR: Internal error (08)");
        }
    }
}

/// Check if a reference system and epoch combination is supported.
///
/// This helper function validates whether the given pair of reference system and epoch
/// identifiers corresponds to a valid configuration for transformation routines like [`rotpn`].
///
/// Supported reference systems:
/// - `"EQUM"`: Equatorial Mean
/// - `"EQUT"`: Equatorial True
/// - `"ECLM"`: Ecliptic Mean
///
/// Supported epochs:
/// - `"J2000"`: Fixed standard epoch
/// - `"OFDATE"`: Epoch of date (variable)
///
/// Arguments
/// ---------
/// * `rsys`: reference system identifier as a string slice.
/// * `epoch`: epoch descriptor as a string slice.
///
/// Returns
/// --------
/// * `true` if the combination is valid and supported, otherwise `false`.
///
/// # See also
/// * [`rotpn`] – which calls this function to validate inputs
fn chkref(rsys: &str, epoch: &str) -> bool {
    matches!(rsys, "EQUM" | "EQUT" | "ECLM") && matches!(epoch, "J2000" | "OFDATE")
}

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
/// * [`rotpn`] – applies obliquity rotation when transforming between ecliptic and equatorial frames
pub fn obleq(tjm: f64) -> Radian {
    // Obliquity coefficients
    let ob0 = ((23.0 * 3600.0 + 26.0 * 60.0) + 21.448) * RADSEC;
    let ob1 = -46.815 * RADSEC;
    let ob2 = -0.0006 * RADSEC;
    let ob3 = 0.00181 * RADSEC;

    let t = (tjm - T2000) / 36525.0;

    ((ob3 * t + ob2) * t + ob1) * t + ob0
}

/// Construct a right-handed rotation matrix around a principal axis (X, Y, or Z).
///
/// This function returns a 3×3 matrix representing an **active rotation**
/// of a vector in ℝ³ by an angle `alpha` around the axis indexed by `k`.
/// The rotation is applied in the **direct (positive/trigonometric)** sense:
/// i.e., counter-clockwise when looking along the axis toward the origin.
///
/// Arguments
/// ---------
/// * `alpha`: rotation angle in radians (positive = direct sense).
/// * `k`: index of the axis of rotation:
///     - `0` → X-axis,
///     - `1` → Y-axis,
///     - `2` → Z-axis.
///
/// Returns
/// --------
/// * A 3×3 rotation matrix `R` such that the rotated vector is `x′ = R · x`.
///
/// Remarks
/// -------
/// * This routine follows the OrbFit convention:
///     - it rotates the **vector** in a fixed reference frame,
///     - i.e., it does **not** represent a change of basis.
/// * The matrix is orthonormal and satisfies `Rᵀ = R⁻¹`.
/// * Used internally in [`prec`], [`rnut80`], and [`rotpn`] for sequential transformations.
///
/// Panics
/// -------
/// * Panics if `k > 2`, as only axes 0–2 (X, Y, Z) are valid.
///
/// # See also
/// * [`prec`] – uses `rotmt` for precession angle rotations
/// * [`rnut80`] – applies sequential `rotmt` matrices for nutation
/// * [`rotpn`] – assembles compound frame transformations
pub fn rotmt(alpha: f64, k: usize) -> [[f64; 3]; 3] {
    if k > 2 {
        panic!("**** ROTMT: k = ??? ****");
    }

    let cosa = alpha.cos();
    let sina = alpha.sin();

    let mut r = [[0.0; 3]; 3];

    let i1 = k;
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
///     - `Δψ`: nutation in longitude [arcseconds]
///     - `Δε`: nutation in obliquity [arcseconds]
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
/// Notes
/// ------
/// * The internal coefficients and series expansion are directly adapted from OrbFit (Fortran), which itself
///   follows the IAU 1980 formulation by Wahr.
/// * This implementation is compatible with OrbFit’s `nutn80` and yields identical output to within machine precision.
///
/// # See also
/// * [`rnut80`] – uses these angles to build the nutation rotation matrix
/// * [`rotpn`] – applies nutation when transforming between EQUT and EQUM systems
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
/// of date (EQUM) to the true equator and equinox of date (EQUT).
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
/// * [`rotpn`] – uses `rnut80` to transform between EQUM and EQUT systems
fn rnut80(tjm: f64) -> [[f64; 3]; 3] {
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

    // Multiply: rp = R2 · R1
    let mut rp = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            rp[i][j] = r2[i][0] * r1[0][j] + r2[i][1] * r1[1][j] + r2[i][2] * r1[2][j];
        }
    }

    // Final nutation matrix: R = R3 · (R2 · R1)
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
/// * [`rotpn`] – uses `prec` when converting between epochs `"OFDATE"` and `"J2000"`
fn prec(tjm: f64, rprec: &mut [[f64; 3]; 3]) {
    // Initialize intermediate matrix for storing R2 · R1
    let mut rp = [[0.0; 3]; 3];

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

    // Compute intermediate matrix: rp = R2 · R1
    for i in 0..3 {
        for j in 0..3 {
            rp[i][j] = r2[i][0] * r1[0][j] + r2[i][1] * r1[1][j] + r2[i][2] * r1[2][j];
        }
    }

    // Final precession matrix: rprec = R3 · (R2 · R1)
    for i in 0..3 {
        for j in 0..3 {
            rprec[i][j] = r3[i][0] * rp[0][j] + r3[i][1] * rp[1][j] + r3[i][2] * rp[2][j];
        }
    }
}

fn matmul(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
    let mut result = [[0.0; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            result[i][j] = (0..3).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    result
}

/// Apply stellar aberration correction to a relative position vector.
///
/// This function computes the apparent position of a target object by applying
/// the first-order correction for stellar aberration due to the observer's velocity.
/// It assumes the classical limit (v ≪ c), using a linear time-delay model.
///
/// Arguments
/// ---------
/// * `xrel`: relative position vector from observer to object [AU].
/// * `vrel`: velocity of the observer relative to the barycenter [AU/day].
///
/// Returns
/// --------
/// * Corrected position vector (same units and directionality as `xrel`),
///   shifted by the aberration effect.
///
/// Formula
/// -------
/// The corrected position is given by:
/// ```text
/// x_corr = xrel − (‖xrel‖ / c) · vrel
/// ```
/// where `c` is the speed of light in AU/day (`VLIGHT_AU`).
///
/// Remarks
/// -------
/// * This function does **not** normalize the output.
/// * Suitable for use in astrometric modeling or when computing apparent direction
///   of celestial objects as seen from a moving observer.
pub(crate) fn correct_aberration(xrel: Vector3<f64>, vrel: Vector3<f64>) -> Vector3<f64> {
    let norm_vector = xrel.norm();
    let dt = norm_vector / VLIGHT_AU;
    xrel - dt * vrel
}

/// Convert a 3D Cartesian position vector to right ascension and declination.
///
/// Given a position vector expressed in Cartesian coordinates (typically in an equatorial frame),
/// this function returns the corresponding right ascension (α), declination (δ), and norm (distance).
///
/// Arguments
/// ---------
/// * `cartesian_position`: 3D position vector in Cartesian coordinates [AU or any length unit].
///
/// Returns
/// --------
/// * Tuple `(α, δ, ρ)`:
///     - `α`: right ascension in radians, in the range [0, 2π).
///     - `δ`: declination in radians, in the range [−π/2, +π/2].
///     - `ρ`: Euclidean norm of the vector (distance to the origin).
///
/// Remarks
/// -------
/// * If the input vector has zero norm, the result is `(0.0, 0.0, 0.0)`.
/// * The RA computation uses `atan2` to preserve quadrant information.
/// * This function is used when converting inertial position vectors to observable angles.
///
/// # See also
/// * [`correct_aberration`] – apply aberration correction before calling this if needed
pub(crate) fn cartesian_to_radec(cartesian_position: Vector3<f64>) -> (f64, f64, f64) {
    let pos_norm = cartesian_position.norm();
    if pos_norm == 0. {
        return (0.0, 0.0, pos_norm);
    }

    let delta = (cartesian_position.z / pos_norm).asin();

    let cos_delta = delta.cos();
    if cos_delta == 0.0 {
        return (0.0, delta, pos_norm);
    }

    let cos_alpha = cartesian_position.x / (pos_norm * cos_delta);
    let sin_alpha = cartesian_position.y / (pos_norm * cos_delta);
    let alpha = sin_alpha.atan2(cos_alpha);
    let alpha = if alpha < 0.0 { alpha + DPI } else { alpha };
    (alpha, delta, pos_norm)
}

#[cfg(test)]
mod ref_system_test {

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
        assert_eq!(rnut, ref_rnut);
    }

    mod test_rotpn {
        use super::*;

        #[test]
        fn test_rotpn_equm() {
            let ref_roteqec = [
                [1.0, 0.0, 0.0],
                [0.0, 0.9174820620691818, 0.3977771559319137],
                [0.0, -0.3977771559319137, 0.9174820620691818],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "EQUM", "J2000", 0., "ECLM", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
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

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "EQUM", "J2000", 0., "EQUT", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_eclm() {
            let ref_roteqec = [
                [1.0, 0.0, 0.0],
                [0.0, 0.9174820620691818, -0.3977771559319137],
                [0.0, 0.3977771559319137, 0.9174820620691818],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "ECLM", "J2000", 0., "EQUM", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
                [
                    0.9999999977217079,
                    6.750247612406132e-5,
                    -3.3881317890172014e-21,
                ],
                [
                    -6.193306258211379e-5,
                    0.9174931942820401,
                    -0.39775147342333544,
                ],
                [
                    -2.6849209338068913e-5,
                    0.3977514725171414,
                    0.9174931963723576,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "ECLM", "J2000", 0., "EQUT", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_equt() {
            let ref_roteqec = [
                [
                    0.9999999977217079,
                    -6.193306258211379e-5,
                    -2.6849209338068913e-5,
                ],
                [
                    6.19323109890795e-5,
                    0.9999999976903892,
                    -2.7993043796858963e-5,
                ],
                [
                    2.6850942970991024e-5,
                    2.799138089948361e-5,
                    0.9999999992477547,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "EQUT", "J2000", 0., "EQUM", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
                [
                    0.9999999977217079,
                    -6.193306258211379e-5,
                    -2.6849209338068913e-5,
                ],
                [6.750247612406132e-5, 0.9174931942820401, 0.3977514725171414],
                [
                    -3.3881317890172014e-21,
                    -0.39775147342333544,
                    0.9174931963723576,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(&mut roteqec, "EQUT", "J2000", 0., "ECLM", "J2000", 0.);
            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_equt_of_date() {
            let ref_roteqec = [
                [
                    0.9999999999808916,
                    5.671879296062708e-6,
                    2.458983466038936e-6,
                ],
                [
                    -5.671991417020452e-6,
                    0.9999999989442864,
                    4.559885773575134e-5,
                ],
                [
                    -2.458724832225838e-6,
                    -4.559887168220644e-5,
                    0.9999999989573487,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "EQUT",
                "OFDATE",
                60725.5,
                "EQUM",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
                [
                    0.9999999999808916,
                    5.671879296062708e-6,
                    2.458983466038936e-6,
                ],
                [
                    -6.181974962369037e-6,
                    0.9174866172186449,
                    0.39776664916313814,
                ],
                [
                    4.235164736271502e-22,
                    -0.39776664917073884,
                    0.9174866172361764,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "EQUT",
                "OFDATE",
                60725.5,
                "ECLM",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_equm_of_date() {
            let ref_roteqec = [
                [
                    0.9999999999808916,
                    -5.671991417020452e-6,
                    -2.458724832225838e-6,
                ],
                [
                    5.671879296062708e-6,
                    0.9999999989442864,
                    -4.559887168220644e-5,
                ],
                [
                    2.458983466038936e-6,
                    4.559885773575134e-5,
                    0.9999999989573487,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "EQUM",
                "OFDATE",
                60725.5,
                "EQUT",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
                [1.0, 0.0, 0.0],
                [0.0, 0.917504753989953, 0.39772481240907737],
                [0.0, -0.39772481240907737, 0.917504753989953],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "EQUM",
                "OFDATE",
                60725.5,
                "ECLM",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_eclm_of_date() {
            let ref_roteqec = [
                [1.0, 0.0, 0.0],
                [0.0, 0.917504753989953, -0.39772481240907737],
                [0.0, 0.39772481240907737, 0.917504753989953],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "ECLM",
                "OFDATE",
                60725.5,
                "EQUM",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);

            let ref_roteqec = [
                [
                    0.9999999999808916,
                    -6.181974962369037e-6,
                    4.235164736271502e-22,
                ],
                [
                    5.671879296062708e-6,
                    0.9174866172186449,
                    -0.39776664917073884,
                ],
                [
                    2.458983466038936e-6,
                    0.39776664916313814,
                    0.9174866172361764,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            rotpn(
                &mut roteqec,
                "ECLM",
                "OFDATE",
                60725.5,
                "EQUT",
                "OFDATE",
                60730.5,
            );

            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_equt_eclm_date() {
            let ref_roteqec = [
                [
                    0.9999932036120499,
                    0.003381495004957589,
                    0.0014690885747894438,
                ],
                [
                    -0.0036868307528666357,
                    0.9174941827437706,
                    0.3977321107357815,
                ],
                [
                    -2.9510755403679666e-6,
                    -0.3977348238749929,
                    0.917500414097138,
                ],
            ];

            let mut roteqec = [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]];
            let tmjd = 57028.479297592596;
            rotpn(&mut roteqec, "EQUT", "OFDATE", tmjd, "ECLM", "J2000", 0.);

            assert_eq!(roteqec, ref_roteqec);
        }

        #[test]
        fn test_rotpn_identity_cases() {
            let mut r = [[0.; 3]; 3];
            // Identity in J2000 EQUM
            rotpn(&mut r, "EQUM", "J2000", 0., "EQUM", "J2000", 0.);
            assert_eq!(r, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);

            // Identity in OFDATE (ECLM)
            let mut r = [[0.; 3]; 3];
            rotpn(&mut r, "ECLM", "OFDATE", 60000., "ECLM", "OFDATE", 60000.);
            assert_eq!(r, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);

            // Identity in OFDATE (EQUT)
            let mut r = [[0.; 3]; 3];
            rotpn(&mut r, "EQUT", "OFDATE", 60000., "EQUT", "OFDATE", 60000.);
            assert_eq!(r, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]);
        }

        #[test]
        fn test_rotpn_inverse_transform() {
            let mut r1 = [[0.; 3]; 3];
            let mut r2 = [[0.; 3]; 3];
            rotpn(&mut r1, "EQUM", "J2000", 0., "ECLM", "J2000", 0.);
            rotpn(&mut r2, "ECLM", "J2000", 0., "EQUM", "J2000", 0.);
            let prod = matmul(&r2, &r1);
            #[allow(clippy::needless_range_loop)]
            for i in 0..3 {
                #[allow(clippy::needless_range_loop)]
                for j in 0..3 {
                    if i == j {
                        assert!((prod[i][j] - 1.0).abs() < 1e-12);
                    } else {
                        assert!(prod[i][j].abs() < 1e-12);
                    }
                }
            }
        }

        #[test]
        fn test_rotpn_large_epoch_difference() {
            let mut r = [[0.; 3]; 3];
            // From 2055 back to J2000 in EQUM
            rotpn(&mut r, "EQUM", "OFDATE", 80000., "EQUM", "J2000", 0.);
            // Just check the rotation matrix is orthonormal (basic sanity)
            let mut rt = r;
            trsp3(&mut rt);
            let prod = matmul(&r, &rt);
            #[allow(clippy::needless_range_loop)]
            for i in 0..3 {
                #[allow(clippy::needless_range_loop)]
                for j in 0..3 {
                    if i == j {
                        assert!((prod[i][j] - 1.0).abs() < 1e-12);
                    } else {
                        assert!(prod[i][j].abs() < 1e-12);
                    }
                }
            }
        }

        #[test]
        #[should_panic(expected = "ERROR: Unsupported starting reference system")]
        fn test_rotpn_invalid_reference_system() {
            let mut r = [[0.; 3]; 3];
            // Unknown system
            rotpn(&mut r, "FOO", "J2000", 0., "EQUM", "J2000", 0.);
        }

        #[test]
        fn test_rotpn_round_trip_equt_equm() {
            let mut forward = [[0.; 3]; 3];
            let mut backward = [[0.; 3]; 3];
            rotpn(
                &mut forward,
                "EQUT",
                "OFDATE",
                60725.5,
                "EQUM",
                "OFDATE",
                60730.5,
            );
            rotpn(
                &mut backward,
                "EQUM",
                "OFDATE",
                60730.5,
                "EQUT",
                "OFDATE",
                60725.5,
            );

            let prod = matmul(&backward, &forward);

            let tol = 1e-5;
            #[allow(clippy::needless_range_loop)]
            for i in 0..3 {
                #[allow(clippy::needless_range_loop)]
                for j in 0..3 {
                    if i == j {
                        assert!((prod[i][j] - 1.0).abs() < tol);
                    } else {
                        assert!(prod[i][j].abs() < tol);
                    }
                }
            }
        }
    }
}
