//! # Reference-frame transformations (equatorial/ecliptic, mean/true, epoch-aware)
//!
//! This module defines the types and routines used to build **3×3 rotation
//! matrices** between common celestial reference systems:
//! - **Equatorial mean** (precession only),
//! - **Equatorial true** (precession + nutation),
//! - **Ecliptic mean** (mean obliquity),
//! each optionally tied to an **epoch** (J2000 or a specific Modified Julian Date, TT).
//!
//! The central entry point is [`rotpn`](crate::ref_system::rotpn), which composes the necessary rotations
//! (precession, nutation, obliquity) to transform vectors between any two supported
//! frames and epochs. All rotations are **active**, right‑handed and orthonormal.
//!
//! ## Coordinate systems & epochs
//!
//! - [`RefSystem::Equm(e)`] — *Equatorial mean* at epoch `e` (precession only).
//! - [`RefSystem::Equt(e)`] — *Equatorial true* at epoch `e` (precession + nutation).
//! - [`RefSystem::Eclm(e)`] — *Ecliptic mean* at epoch `e` (mean obliquity).
//!
//! Epochs are represented by [`RefEpoch`](crate::ref_system::RefEpoch):
//! - [`RefEpoch::J2000`](crate::ref_system::RefEpoch::J2000) — Fixed epoch at **MJD 51544.5 (TT)**.
//! - [`RefEpoch::Epoch(d)`] — “of-date” epoch at MJD `d` (TT).
//!
//! **Units & conventions**
//! -----------------------
//! - Angles are in **radians**.
//! - Dates are **Modified Julian Date** in **Terrestrial Time (TT)**.
//! - Rotations are **active**: the matrix `R` rotates vectors (`x_dst = R · x_src`).
//! - All matrices are **orthonormal**: `R.transpose() == R.inverse()`.
//! - Right‑handed axes (X,Y,Z). See [`rotmt`](crate::ref_system::rotmt) for axis indexing.
//!
//! ## Mathematical models
//!
//! - **Precession**: IAU 1976, via [`prec`](crate::earth_orientation::prec) (matrix from/to J2000 or of‑date).
//! - **Nutation**: IAU 1980, via [`rnut80`](crate::earth_orientation::rnut80) (true ↔ mean equator/equinox).
//! - **Mean obliquity**: via [`obleq`](crate::earth_orientation::obleq) (ecliptic ↔ equatorial).
//!
//! The composition strategy used by [`rotpn`](crate::ref_system::rotpn) is intentionally explicit and
//! mirror‑friendly to classic astrometry codes (e.g., OrbFit): whenever helpful,
//! transformations are routed through **Equatorial Mean J2000** as an intermediate,
//! ensuring determinism and easing verification against reference implementations.
//!
//! ## Typical usage
//!
//! ```rust, no_run
//! use outfit::ref_system::{RefEpoch, RefSystem, rotpn};
//! use nalgebra::{Matrix3, Vector3};
//!
//! // Equatorial mean J2000 → Ecliptic mean J2000
//! let src = RefSystem::Equm(RefEpoch::J2000);
//! let dst = RefSystem::Eclm(RefEpoch::J2000);
//! let r_eq_to_ec: Matrix3<f64> = rotpn(&src, &dst)?;
//!
//! // Apply to a position vector (active rotation):
//! let x_eq = Vector3::new(1.0, 0.0, 0.0);
//! let x_ec = r_eq_to_ec * x_eq;
//!
//! // “Of-date” example: Equatorial true @ t1 → Equatorial mean @ t2
//! let t1 = RefEpoch::Epoch(60725.5);
//! let t2 = RefEpoch::Epoch(60730.5);
//! let r_t1_to_t2 = rotpn(&RefSystem::Equt(t1), &RefSystem::Equm(t2))?;
//! # Ok::<(), outfit::outfit_errors::OutfitError>(())
//! ```
//!
//! ## API overview
//!
//! - [`RefEpoch`](crate::ref_system::RefEpoch) — Epoch tagging (J2000 vs. MJD(TT) “of‑date”).
//! - [`RefSystem`](crate::ref_system::RefSystem) — Frame tagging (equatorial/ecliptic, mean/true) + epoch.
//! - [`rotpn`](crate::ref_system::rotpn) — Build the **full rotation** from any source to any destination frame.
//! - [`rotmt`](crate::ref_system::rotmt) — Primitive axis‑angle rotation used by precession/nutation steps.
//!
//! Internally, [`RefSystem`](crate::ref_system::RefSystem) offers two helpers used by [`rotpn`](crate::ref_system::rotpn):
//! - `transform_to_equm_date` — Normalize to **Equatorial Mean** at a target epoch
//!   (handles precession, removes nutation/obliquity as needed).
//! - `transform_to_target_system` — Switch **system type** at fixed epoch
//!   (apply/remove nutation or obliquity).
//!
//! ## Numerical behavior & guarantees
//!
//! - **Round‑trip** stability: forward/backward transforms compose to (near) identity;
//!   test coverage verifies this for representative cases.
//! - **Orthogonality** is preserved within floating‑point tolerance.
//! - **Epoch equality** uses a small tolerance [`EPS`](crate::constants::EPS) to avoid spurious precession steps.
//! - A safety cap (20 iterations) protects against non‑convergent sequences; exceeding it
//!   yields an [`OutfitError::InvalidRefSystem`](crate::outfit_errors::OutfitError::InvalidRefSystem).
//!
//! ## Common pitfalls
//!
//! - **Active vs. passive**: results here are **active** rotations (rotate vectors).
//!   If you need change‑of‑basis, use the transpose/inverse accordingly.
//! - **Time scale**: epochs must be in **TT**; mixing with UTC/TDB without conversion
//!   will introduce arcsecond‑level errors over long spans.
//! - **Axis index**: [`rotmt`](crate::ref_system::rotmt) panics for `k > 2` (valid axes: 0→X, 1→Y, 2→Z).
//!
//! ## Testing hints
//!
//! - Identity cases when source and destination frames match (including of‑date).
//! - Inverse consistency: `rotpn(dst, src) * rotpn(src, dst) ≈ I`.
//! - Large epoch spans: check orthonormality and non‑identity behavior.
//! - Sensitivity: small but non‑zero precession deltas over ~years (diagonal drift).
//!
//! ## See also
//! ------------
//! * [`rotpn`](crate::ref_system::rotpn) – Build compound frame transformations (precession/nutation/obliquity).
//! * [`rotmt`](crate::ref_system::rotmt) – Right‑handed axis rotations used by the low‑level models.
//! * [`prec`](crate::earth_orientation::prec) – IAU 1976 precession matrix builder.
//! * [`rnut80`](crate::earth_orientation::rnut80) – IAU 1980 nutation rotation.
//! * [`obleq`](crate::earth_orientation::obleq) – Mean obliquity (radians) used for ecliptic ↔ equatorial.
use nalgebra::{Matrix3, Rotation3, Vector3};

use crate::{
    earth_orientation::{obleq, prec, rnut80},
    outfit_errors::OutfitError,
};

use super::constants::{EPS, T2000};

/// Represents the reference epoch associated with a celestial reference frame.
///
/// The epoch defines the time at which the frame (mean equator, true equator,
/// ecliptic, etc.) is valid. It influences the precession and nutation corrections
/// used to transform coordinates.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefEpoch {
    /// Standard J2000 epoch, corresponding to MJD 51544.5 (TT).
    J2000,

    /// A custom epoch, expressed as a Modified Julian Date (TT).
    ///
    /// For example, `RefEpoch::Epoch(60000.0)` corresponds to
    /// MJD = 60000.0 (TT), which can be used for "of-date" transformations.
    Epoch(f64),
}

impl RefEpoch {
    /// Return the epoch as a Modified Julian Date (TT).
    ///
    /// # Returns
    /// * For [`RefEpoch::J2000`]: returns the constant [`T2000`].
    /// * For [`RefEpoch::Epoch(d)`]: returns the stored `d`.
    pub fn date(&self) -> f64 {
        match *self {
            RefEpoch::J2000 => T2000,
            RefEpoch::Epoch(d) => d,
        }
    }
}

/// Represents a celestial reference system, including both
/// the **coordinate basis** (equatorial/ecliptic) and the **epoch**.
///
/// # Variants
///
/// * `Equm(RefEpoch)` – **Equatorial mean**:
///   Based on the mean equator and equinox (precession only).
///
/// * `Equt(RefEpoch)` – **Equatorial true**:
///   Based on the true equator and equinox (precession + nutation).
///
/// * `Eclm(RefEpoch)` – **Ecliptic mean**:
///   Based on the mean ecliptic and mean equinox (precession + obliquity).
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RefSystem {
    /// Equatorial mean (precession only, no nutation)
    Equm(RefEpoch),

    /// Equatorial true (precession + nutation)
    Equt(RefEpoch),

    /// Ecliptic mean (precession + obliquity)
    Eclm(RefEpoch),
}

impl RefSystem {
    /// Returns the [`RefEpoch`] associated with this reference system.
    pub fn epoch(&self) -> RefEpoch {
        match *self {
            RefSystem::Equm(e) => e,
            RefSystem::Equt(e) => e,
            RefSystem::Eclm(e) => e,
        }
    }

    /// Compare only the type of system (Equm/Equt/Eclm) while ignoring the epoch.
    ///
    /// # Returns
    /// `true` if `self` and `other` are of the same variant, regardless of epoch.
    ///
    /// # Example
    ///
    /// ```rust, ignore
    /// assert!(RefSystem::Equm(RefEpoch::J2000)
    ///     .variant_eq(&RefSystem::Equm(RefEpoch::Epoch(60000.0))));
    /// ```
    pub fn variant_eq(&self, other: &RefSystem) -> bool {
        matches!(
            (self, other),
            (RefSystem::Equm(_), RefSystem::Equm(_))
                | (RefSystem::Equt(_), RefSystem::Equt(_))
                | (RefSystem::Eclm(_), RefSystem::Eclm(_))
        )
    }

    /// Convert the current reference system to an **equatorial mean (Equm)** system
    /// at a specified target epoch.
    ///
    /// This function is typically the **first step** in a reference frame transformation.
    /// It normalizes the current frame so that:
    ///
    /// * The reference plane becomes **equatorial** (removing ecliptic obliquity if needed).
    /// * Nutation effects are removed (if the current frame is true equator/equinox).
    /// * The epoch is converted to the requested target epoch (if necessary) using precession.
    ///
    /// After this transformation, the result is always an `Equm` system, ready for
    /// further transformations (e.g., conversion to another system with `transform_to_target_system`).
    ///
    /// # Behavior
    ///
    /// Depending on the **current epoch** and **frame type**, the following actions are performed:
    ///
    /// ### If the current epoch is `J2000`
    /// - `Eclm(J2000)` → `Equm(J2000)` using an obliquity rotation.
    /// - `Equt(J2000)` → `Equm(J2000)` using the **inverse** nutation matrix.
    /// - `Equm(J2000)` → `Equm(target_epoch)` using a **precession** matrix.
    ///
    /// ### If the current epoch is a custom epoch `Epoch(d)`
    /// - `Eclm(d)` → `Equm(d)` using an obliquity rotation.
    /// - `Equt(d)` → `Equm(d)` using the **inverse** nutation matrix.
    /// - `Equm(d)` → `Equm(J2000)` using the **transpose** of the precession matrix.
    ///
    /// Note that when the starting epoch is not J2000, this function may bring
    /// the frame **back to J2000** before applying other transformations.
    ///
    /// # Returns
    ///
    /// A tuple `(next_ref, rot)` where:
    /// * `next_ref` – The new reference system, always `RefSystem::Equm`.
    /// * `rot` – The 3×3 rotation matrix to convert a vector from the original frame to `next_ref`.
    ///
    /// # Errors
    ///
    /// Returns an `Err(OutfitError)` if the conversion request is redundant, for example
    /// `Equm(J2000)` to `Equm(J2000)`.
    ///
    /// # See also
    ///
    /// * [`RefSystem::transform_to_target_system`] – For changing between Equm/Equt/Eclm once the epochs match.
    /// * [`prec`] – Precession matrix (IAU 1976).
    /// * [`rnut80`] – Nutation model (IAU 1980).
    /// * [`obleq`] – Mean obliquity of the ecliptic.
    fn transform_to_equm_date(
        &self,
        target_epoch: RefEpoch,
    ) -> Result<(RefSystem, Matrix3<f64>), OutfitError> {
        match self.epoch() {
            RefEpoch::J2000 => match self {
                RefSystem::Eclm(e) => Ok((RefSystem::Equm(*e), rotmt(-obleq(T2000), 1))),
                RefSystem::Equt(e) => Ok((RefSystem::Equm(*e), rnut80(T2000).transpose())),
                RefSystem::Equm(_) => match target_epoch {
                    RefEpoch::J2000 => Err(OutfitError::InvalidRefSystem(
                        "Cannot convert Equm to Equm J2000, already in that system".into(),
                    )),
                    RefEpoch::Epoch(target_date) => Ok((
                        RefSystem::Equm(RefEpoch::Epoch(target_date)),
                        prec(target_date),
                    )),
                },
            },
            RefEpoch::Epoch(_) => match self {
                RefSystem::Eclm(e) => Ok((RefSystem::Equm(*e), rotmt(-obleq(e.date()), 1))),
                RefSystem::Equt(e) => Ok((RefSystem::Equm(*e), rnut80(e.date()).transpose())),
                RefSystem::Equm(e) => {
                    Ok((RefSystem::Equm(RefEpoch::J2000), prec(e.date()).transpose()))
                }
            },
        }
    }

    /// Apply a **system-type transformation** (nutation or obliquity) without changing the epoch.
    ///
    /// This function converts between the three reference system variants
    /// (`Equm`, `Equt`, `Eclm`) at a fixed epoch. Depending on the current system:
    ///
    /// * `Equt(e)` → `Equm(e)`: removes nutation (applies the inverse nutation matrix).
    /// * `Eclm(e)` → `Equm(e)`: removes obliquity (rotates back from ecliptic to equatorial).
    /// * `Equm(e)` → `Equt(e)`: applies nutation.
    /// * `Equm(e)` → `Eclm(e)`: applies obliquity rotation.
    ///
    /// This function **does not** handle precession or epoch changes; it only
    /// switches between system types at the same epoch.
    ///
    /// # Returns
    ///
    /// A tuple `(next_ref, rot)` where:
    /// * `next_ref` – The resulting reference system after applying the transformation.
    /// * `rot` – A 3×3 rotation matrix that converts coordinates from `self` to `next_ref`.
    ///
    /// # Errors
    ///
    /// Returns an `Err(OutfitError)` if no transformation is needed
    /// (e.g. Equm → Equm).
    fn transform_to_target_system(
        &self,
        target_system: RefSystem,
    ) -> Result<(RefSystem, Matrix3<f64>), OutfitError> {
        match self {
            RefSystem::Equt(e) => Ok((RefSystem::Equm(*e), rnut80(e.date()).transpose())),
            RefSystem::Eclm(e) => Ok((RefSystem::Equm(*e), rotmt(-obleq(e.date()), 0))),
            RefSystem::Equm(e) => match target_system {
                RefSystem::Equt(_) => Ok((RefSystem::Equt(*e), rnut80(e.date()))),
                RefSystem::Eclm(_) => Ok((RefSystem::Eclm(*e), rotmt(obleq(e.date()), 0))),
                RefSystem::Equm(_) => Err(OutfitError::InvalidRefSystem(
                    "Cannot convert Equm to Equm, already in that system".into(),
                )),
            },
        }
    }
}

/// Compute the 3×3 rotation matrix between two celestial reference systems and epochs.
///
/// This function builds a composite rotation matrix that transforms a vector from a **source**
/// reference system (with its epoch) to a **destination** reference system (with its epoch).
///
/// The supported reference system variants are:
/// - [`RefSystem::Equm`] – Equatorial mean (precession only),
/// - [`RefSystem::Equt`] – Equatorial true (precession + nutation),
/// - [`RefSystem::Eclm`] – Ecliptic mean (precession + obliquity).
///
/// Supported epochs:
/// - [`RefEpoch::J2000`] – Standard J2000 epoch (fixed at MJD 51544.5),
/// - [`RefEpoch::Epoch`] – Epoch tied to a specific date (e.g. time of observation).
///
/// ## Transformation procedure
///
/// If the two frames differ by epoch and/or system, the resulting rotation matrix is constructed
/// as a product of simpler transformations:
///
/// 1. **Epoch alignment** – Applies precession to match the epoch of the destination frame.
/// 2. **System alignment** – Applies nutation and/or obliquity rotation to convert between
///    equatorial/ecliptic and mean/true frames.
///
/// Internally, transformations may pass through the canonical frame
/// **Equatorial Mean J2000** as an intermediate step.
///
/// The resulting rotation satisfies:
///
/// ```text
/// x_dst = R · x_src
/// ```
///
/// where `x_src` are coordinates in the source frame, and `x_dst` in the destination frame.
///
/// # Arguments
///
/// * `src` – The source reference frame, including its epoch.
/// * `dst` – The destination reference frame, including its epoch.
///
/// # Returns
///
/// A [`Matrix3<f64>`] that rotates vectors from `src` to `dst`.
///
/// # Algorithm details
///
/// * Precession model: IAU 1976 (`prec`)
/// * Nutation model: IAU 1980 (`rnut80`)
/// * Mean obliquity: [`obleq`] (in radians)
///
/// The algorithm proceeds iteratively, modifying an internal reference frame (`current`) until
/// it matches the target frame. The cumulative product of all stepwise rotations is returned.
///
/// # Errors
///
/// Returns an [`OutfitError`] if:
/// * An unsupported reference system or epoch is encountered,
/// * The transformation does not converge within 20 iterations (safety cap).
///
/// # See also
///
/// * [`prec`] – IAU 1976 precession matrix
/// * [`rnut80`] – IAU 1980 nutation model
/// * [`rotmt`] – Rotation around a principal axis
/// * [`obleq`] – Mean obliquity of the ecliptic
pub fn rotpn(src: &RefSystem, dst: &RefSystem) -> Result<Matrix3<f64>, OutfitError> {
    let mut current = *src;
    let mut rotation = Matrix3::identity();

    for _ in 0..20 {
        let epochs_equal = match (current.epoch(), dst.epoch()) {
            (RefEpoch::J2000, RefEpoch::J2000) => true,
            (e1, e2) => (e1.date() - e2.date()).abs() <= EPS,
        };

        if !epochs_equal {
            // Step 1: adjust the epoch
            let (next_ref, step_rot) = current.transform_to_equm_date(dst.epoch())?;
            rotation *= step_rot;
            current = next_ref;
            continue;
        }

        // Epochs are now aligned, check if reference systems match
        if current.variant_eq(dst) {
            return Ok(rotation);
        }

        // Step 2: align the reference system (nutation / obliquity)
        let (next_ref, step_rot) = current.transform_to_target_system(*dst)?;
        rotation *= step_rot;
        current = next_ref;
    }

    Err(OutfitError::InvalidRefSystem(
        "Transformation did not converge in 20 iterations".into(),
    ))
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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Equt(RefEpoch::J2000)).unwrap();
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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Equt(RefEpoch::J2000)).unwrap();
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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
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

        let roteqec = rotpn(&ref_sys1, &RefSystem::Eclm(RefEpoch::J2000)).unwrap();
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
        let rot = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

        assert_eq!(roteqec, ref_roteqec.into());
    }

    #[test]
    fn test_rotpn_identity_cases() {
        // Identity in J2000 Equm
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Equm(RefEpoch::J2000);
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());

        // Identity in OFDATE (Eclm)
        let ref_sys1 = RefSystem::Eclm(RefEpoch::Epoch(60000.));
        let ref_sys2 = RefSystem::Eclm(RefEpoch::Epoch(60000.));
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());

        // Identity in OFDATE (Equt)
        let ref_sys1 = RefSystem::Equt(RefEpoch::Epoch(60000.));
        let ref_sys2 = RefSystem::Equt(RefEpoch::Epoch(60000.));
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
        assert_eq!(roteqec, [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]].into());
    }

    #[test]
    fn test_rotpn_inverse_transform() {
        let ref_sys1 = RefSystem::Equm(RefEpoch::J2000);
        let ref_sys2 = RefSystem::Eclm(RefEpoch::J2000);

        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();
        let roteqec2 = rotpn(&ref_sys2, &ref_sys1).unwrap();
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
        let roteqec = rotpn(&ref_sys1, &ref_sys2).unwrap();

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
        let forward = rotpn(&ref_sys1, &ref_sys2).unwrap();
        let backward = rotpn(&ref_sys2, &ref_sys1).unwrap();

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
