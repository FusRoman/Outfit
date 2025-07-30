use nalgebra::{Matrix3, Rotation3, Vector3};

use crate::earth_orientation::{obleq, prec, rnut80};

use super::constants::{EPS, T2000};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefEpoch {
    J2000,
    Epoch(f64),
}

impl RefEpoch {
    pub fn date(&self) -> f64 {
        match *self {
            RefEpoch::J2000 => T2000,
            RefEpoch::Epoch(d) => d,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefSystem {
    // Equatorial Mean, equatorial coordinates based on equator and mean equinox
    // at a given epoch (J2000 for instance)
    // (corrected for precession but not for nutation)
    Equm(RefEpoch),
    // Equatorial True (same as Equm but corrected for precession and nutation)
    Equt(RefEpoch),
    // Ecliptic mean, ecliptic coordinates based on ecliptic and mean equinox
    // at a given epoch (J2000 for instance)
    Eclm(RefEpoch),
}

impl RefSystem {
    pub fn epoch(&self) -> RefEpoch {
        match *self {
            RefSystem::Equm(e) => e,
            RefSystem::Equt(e) => e,
            RefSystem::Eclm(e) => e,
        }
    }

    /// Compare only the variant (Equm/Equt/Eclm) and ignores the epoch value.
    ///
    /// Returns `true` if both `RefSystem` values are of the same variant,
    /// regardless of the inner `RefEpoch`.
    pub fn variant_eq(&self, other: &RefSystem) -> bool {
        matches!(
            (self, other),
            (RefSystem::Equm(_), RefSystem::Equm(_))
                | (RefSystem::Equt(_), RefSystem::Equt(_))
                | (RefSystem::Eclm(_), RefSystem::Eclm(_))
        )
    }
}

/// Compute the rotation matrix between two celestial reference systems and epochs.
///
/// This function builds a composite rotation matrix that transforms coordinates
/// from a source reference system and epoch to a target system and epoch. The supported
/// systems are:
/// - `"Equm"`: equatorial mean (precession only)
/// - `"Equt"`: equatorial true (precession + nutation)
/// - `"Eclm"`: ecliptic mean (precession + obliquity)
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
/// * `rsys1`: name of the source reference system (`"Equm"`, `"Equt"`, `"Eclm"`).
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
/// * The transformation path may pass through `"Equm", "J2000"` as an intermediate canonical frame.
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
pub fn rotpn(ref_sys1: &RefSystem, ref_sys2: &RefSystem) -> Matrix3<f64> {
    let mut rsys = *ref_sys1;
    let mut epoch = ref_sys1.epoch();
    let mut date = if epoch == RefEpoch::J2000 {
        T2000
    } else {
        epoch.date()
    };

    let mut rot = Matrix3::identity();

    let mut nit = 0;

    loop {
        let epdif = if epoch == ref_sys2.epoch() {
            if epoch == RefEpoch::J2000 {
                false
            } else {
                (date - ref_sys2.epoch().date()).abs() > EPS
            }
        } else {
            true
        };

        if epdif {
            if epoch != RefEpoch::J2000 {
                if let RefSystem::Eclm(e) = rsys {
                    let obl = obleq(date);
                    let r = rotmt(-obl, 1);
                    rot *= r;
                    rsys = RefSystem::Equm(e);
                } else if let RefSystem::Equt(e) = rsys {
                    let r = rnut80(date);
                    rot *= r.transpose();
                    rsys = RefSystem::Equm(e);
                } else if let RefSystem::Equm(_) = rsys {
                    let r = prec(date);
                    rot *= r.transpose();
                    epoch = RefEpoch::J2000;
                    date = T2000;
                } else {
                    panic!("ERROR: Internal error (03)");
                }
            } else if let RefSystem::Eclm(e) = rsys {
                let obl = obleq(T2000);
                let r = rotmt(-obl, 1);
                rot *= r;
                rsys = RefSystem::Equm(e);
            } else if let RefSystem::Equt(e) = rsys {
                let r = rnut80(T2000);
                rot *= r.transpose();
                rsys = RefSystem::Equm(e);
            } else if let RefSystem::Equm(_) = rsys {
                if let RefEpoch::Epoch(_) = ref_sys2.epoch() {
                    let r = prec(ref_sys2.epoch().date());
                    rot *= r;
                    epoch = ref_sys2.epoch();
                    date = ref_sys2.epoch().date();
                } else {
                    panic!("ERROR: Internal error (04)");
                }
            } else {
                panic!("ERROR: Internal error (05)");
            }
        } else {
            if rsys.variant_eq(ref_sys2) {
                return rot;
            }

            if let RefSystem::Equt(e) = rsys {
                let r = rnut80(date);
                rot *= r.transpose();
                rsys = RefSystem::Equm(e);
            } else if let RefSystem::Eclm(e) = rsys {
                let obl = obleq(date);
                let r = rotmt(-obl, 0);
                rot *= r;
                rsys = RefSystem::Equm(e);
            } else if let RefSystem::Equm(e) = rsys {
                if let RefSystem::Equt(_) = ref_sys2 {
                    let r = rnut80(date);
                    rot *= r;
                    rsys = RefSystem::Equt(e);
                } else if let RefSystem::Eclm(_) = ref_sys2 {
                    let obl = obleq(date);
                    let r = rotmt(obl, 0);
                    rot *= r;
                    rsys = RefSystem::Eclm(e);
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

/// Construct a right-handed 3×3 rotation matrix around one of the principal axes (X, Y, or Z).
///
/// This function builds a [`nalgebra::Matrix3`] representing an **active rotation**
/// of a 3D vector by an angle `alpha` around the chosen axis.
/// The rotation follows the **direct (positive/trigonometric)** sense:
/// counter-clockwise when looking **along the axis toward the origin**.
///
/// # Arguments
///
/// * `alpha` - Rotation angle in **radians** (positive = direct/trigonometric sense).
/// * `k` - Index of the axis of rotation:
///   * `0` → X-axis
///   * `1` → Y-axis
///   * `2` → Z-axis
///
/// # Returns
///
/// A 3×3 rotation matrix `R` such that the rotated vector is `x' = R · x`.
///
/// # Remarks
///
/// * This function uses [`nalgebra::Rotation3::from_axis_angle`] internally,
///   which ensures orthonormality and numerical stability.
/// * The returned matrix is **orthonormal** and satisfies `R.transpose() == R.inverse()`.
/// * This corresponds to OrbFit's convention:
///   - the rotation is **applied to the vector** in a fixed frame,
///     and does **not** represent a change of basis.
/// * Internally used in:
///   - [`prec`] for precession
///   - [`rnut80`] for nutation
///   - [`rotpn`] for building compound frame transformations
///
/// # Panics
///
/// Panics if `k > 2`, as only axes 0–2 are valid.
///
/// # See also
/// * [`prec`] – uses `rotmt` for precession angle rotations
/// * [`rnut80`] – applies sequential `rotmt` matrices for nutation
/// * [`rotpn`] – assembles compound frame transformations
pub fn rotmt(alpha: f64, k: usize) -> Matrix3<f64> {
    let axis = match k {
        0 => Vector3::x_axis(),
        1 => Vector3::y_axis(),
        2 => Vector3::z_axis(),
        _ => panic!("**** ROTMT: invalid axis index {k} (must be 0,1,2) ****"),
    };

    Rotation3::from_axis_angle(&axis, alpha).into()
}

#[cfg(test)]
mod ref_system_test {

    use super::*;

    use approx::assert_relative_eq;

    fn assert_matrix_eq(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3], tol: f64) {
        for i in 0..3 {
            for j in 0..3 {
                assert_relative_eq!(a[i][j], b[i][j], epsilon = tol);
            }
        }
    }

    const TOLERANCE: f64 = 1e-10;

    #[test]
    fn test_rotpn_equm() {
        let ref_roteqec = [
            [1.0, 0.0, 0.0],
            [0.0, 0.9174820620691818, 0.3977771559319137],
            [0.0, -0.3977771559319137, 0.9174820620691818],
        ];

        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, ref_roteqec.into());

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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Equt(RefEpoch::J2000));
        assert_eq!(roteqec, ref_roteqec.into());
    }

    #[test]
    fn test_rotpn_eclm() {
        let ref_roteqec = [
            [1.0, 0.0, 0.0],
            [0.0, 0.9174820620691818, -0.3977771559319137],
            [0.0, 0.3977771559319137, 0.9174820620691818],
        ];

        let ref_sys1 = RefSystem::Eclm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, ref_roteqec.into());

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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Equt(RefEpoch::J2000));
        assert_eq!(roteqec, ref_roteqec.into());
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

        let ref_sys1 = RefSystem::Equt(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, ref_roteqec.into());

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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Eclm(RefEpoch::J2000));
        assert_eq!(roteqec, ref_roteqec.into());
    }

    #[test]
    fn test_rotpn_ofdate() {
        // Choose a 1000-day offset (~2.7 years) to make the precession effect visible
        let date1 = 60000.0;
        let date2 = 61000.0;

        // Compute transformation from Equm@OFDATE(date1) to Equm@OFDATE(date2)
        let ref_sys1 = RefSystem::Equm(RefEpoch::Epoch(date1));
        let ref_sys2 = RefSystem::Equm(RefEpoch::Epoch(date2));
        let rot = rotpn(&ref_sys1, &ref_sys2);

        // -------------------------------------------------------------------------
        // 1. The resulting rotation matrix should NOT be the identity matrix
        // -------------------------------------------------------------------------
        let tol = 1e-12;
        let mut is_identity = true;

        #[allow(clippy::needless_range_loop)]
        for i in 0..3 {
            #[allow(clippy::needless_range_loop)]
            for j in 0..3 {
                let expected = if i == j { 1.0 } else { 0.0 };
                if (rot[(i, j)] - expected).abs() > tol {
                    is_identity = false;
                }
            }
        }
        assert!(
            !is_identity,
            "Rotation matrix should not be identity when date1 != date2 for OFDATE"
        );

        // -------------------------------------------------------------------------
        // 2. For 1000 days, the diagonal term (rot[1][1]) must differ from 1 by at least 1e-7
        // This ensures that the precession effect has been applied.
        // -------------------------------------------------------------------------
        let delta = (1.0 - rot[(1, 1)]).abs();
        assert!(
            delta > 1e-7,
            "rot[1][1] difference too small: {delta}, expected > 1e-7"
        );
    }

    #[test]
    fn test_rotpn_equt_of_date() {
        let ref_sys1 = RefSystem::Equt(RefEpoch::Epoch(60725.5));
        let ref_sys2 = RefSystem::Equm(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        let expected = [
            [
                0.9999999999959558,
                2.6103210920298055e-6,
                1.1287777487165376e-6,
            ],
            [
                -2.610372560299571e-6,
                0.9999999989569648,
                4.559886322796942e-5,
            ],
            [
                -1.1286587198650923e-6,
                -4.559886617430879e-5,
                0.9999999989597347,
            ],
        ];

        assert_matrix_eq(&roteqec.into(), &expected, TOLERANCE);

        let ref_roteqec = [
            [
                0.9999999999959558,
                2.6103210920298055e-6,
                1.1287777487165376e-6,
            ],
            [
                -2.8439248114746454e-6,
                0.9174866295910213,
                0.3977666206629458,
            ],
            [
                2.660107394168916e-9,
                -0.3977666206645475,
                0.9174866295947346,
            ],
        ];

        let ref_sys2 = RefSystem::Eclm(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_matrix_eq(&roteqec.into(), &ref_roteqec, TOLERANCE);
    }

    #[test]
    fn test_rotpn_equm_of_date() {
        let ref_roteqec = [
            [
                0.9999999999382557,
                -1.019473782042265e-5,
                -4.422167976508847e-6,
            ],
            [
                1.0194536102237101e-5,
                0.9999999989077697,
                -4.561284900943888e-5,
            ],
            [
                4.4226329827165825e-6,
                4.561280392464384e-5,
                0.9999999989499561,
            ],
        ];

        let ref_sys1 = RefSystem::Equm(RefEpoch::Epoch(60725.5));
        let ref_sys2 = RefSystem::Equt(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_matrix_eq(&roteqec.into(), &ref_roteqec, TOLERANCE);

        let ref_roteqec = [
            [
                0.9999999999944286,
                -3.0616188567489498e-6,
                -1.330066112371995e-6,
            ],
            [
                3.3380501509251515e-6,
                0.9175047663420967,
                0.39772478390011357,
            ],
            [
                2.660299467132395e-9,
                -0.39772478390233756,
                0.9175047663472047,
            ],
        ];

        let ref_sys2 = RefSystem::Eclm(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_matrix_eq(&roteqec.into(), &ref_roteqec, TOLERANCE);
    }

    #[test]
    fn test_rotpn_eclm_of_date() {
        let ref_roteqec = [
            [
                0.9175052829851363,
                -3.0616188567489498e-6,
                0.3977235920648803,
            ],
            [
                2.809050665755966e-6,
                0.9999999999953132,
                1.2176799173935054e-6,
            ],
            [
                -0.3977235920667443,
                -2.0361171295958094e-12,
                0.9175052829894363,
            ],
        ];

        let ref_sys1 = RefSystem::Eclm(RefEpoch::Epoch(60725.5));
        let ref_sys2 = RefSystem::Equm(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_matrix_eq(&roteqec.into(), &ref_roteqec, TOLERANCE);

        let ref_roteqec = [
            [
                0.9175065127392313,
                -1.019473782042265e-5,
                0.3977207550243788,
            ],
            [
                2.7494897154247754e-5,
                0.9999999989077697,
                -3.779538585032655e-5,
            ],
            [-0.397720754204662, 4.561280392464384e-5, 0.9175065120174062],
        ];

        let ref_sys2 = RefSystem::Equt(RefEpoch::Epoch(60730.5));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_matrix_eq(&roteqec.into(), &ref_roteqec, TOLERANCE);
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

        let tmjd = 57028.479297592596;
        let ref_sys1 = RefSystem::Equt(RefEpoch::Epoch(tmjd));
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        assert_eq!(roteqec, ref_roteqec.into());
    }

    #[test]
    fn test_rotpn_identity_cases() {
        // Identity in J2000 Equm
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());

        // Identity in OFDATE (Eclm)
        let ref_sys1 = RefSystem::Eclm(RefEpoch::Epoch(60000.));
        let ref_sys2 = RefSystem::Eclm(RefEpoch::Epoch(60000.));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());

        // Identity in OFDATE (Equt)
        let ref_sys1 = RefSystem::Equt(RefEpoch::Epoch(60000.));
        let ref_sys2 = RefSystem::Equt(RefEpoch::Epoch(60000.));
        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());
    }

    #[test]
    fn test_rotpn_inverse_transform() {
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);

        let roteqec = rotpn(&ref_sys1, &ref_sys2);
        let roteqec2 = rotpn(&ref_sys2, &ref_sys1);
        let prod = roteqec2 * roteqec;
        #[allow(clippy::needless_range_loop)]
        for i in 0..3 {
            #[allow(clippy::needless_range_loop)]
            for j in 0..3 {
                if i == j {
                    assert!((prod[(i, j)] - 1.0).abs() < 1e-12);
                } else {
                    assert!(prod[(i, j)].abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_rotpn_large_epoch_difference() {
        // From 2055 back to J2000 in Equm

        let ref_sys1 = RefSystem::Equm(RefEpoch::Epoch(80000.));
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2);

        // Just check the rotation matrix is orthonormal (basic sanity)
        let prod = roteqec * roteqec.transpose();
        #[allow(clippy::needless_range_loop)]
        for i in 0..3 {
            #[allow(clippy::needless_range_loop)]
            for j in 0..3 {
                if i == j {
                    assert!((prod[(i, j)] - 1.0).abs() < 1e-12);
                } else {
                    assert!(prod[(i, j)].abs() < 1e-12);
                }
            }
        }
    }

    #[test]
    fn test_rotpn_round_trip_equt_equm() {
        let ref_sys1 = RefSystem::Equt(RefEpoch::Epoch(60725.5));
        let ref_sys2 = RefSystem::Equm(RefEpoch::Epoch(60730.5));
        let forward = rotpn(&ref_sys1, &ref_sys2);
        let backward = rotpn(&ref_sys2, &ref_sys1);

        let prod = backward * forward;

        let tol = 1e-5;
        #[allow(clippy::needless_range_loop)]
        for i in 0..3 {
            #[allow(clippy::needless_range_loop)]
            for j in 0..3 {
                if i == j {
                    assert!((prod[(i, j)] - 1.0).abs() < tol);
                } else {
                    assert!(prod[(i, j)].abs() < tol);
                }
            }
        }
    }
}
