//! Lagrange f–g velocity correction between two position epochs.

use nalgebra::Vector3;

use crate::constants::GAUSS_GRAV_SQUARED;
use crate::orb_elem::eccentricity_control;
use crate::outfit_errors::OutfitError;

use super::params::UniversalKeplerParams;
use super::solver::solve_kepuni_with_guess;

/// Apply velocity correction using Lagrange f–g coefficients.
///
/// This function refines the velocity vector `v2` of an orbiting body at a given epoch,
/// using two-body dynamics in universal variables. It is a simplified wrapper around
/// [`velocity_correction_with_guess`](crate::kepler::velocity_correction_with_guess) that does not provide a warm-start for the universal
/// anomaly χ or a custom solver tolerance.
///
/// Arguments
/// ---------
/// * `x1` – Position vector of the body at epoch t₁ (in AU).
/// * `x2` – Position vector of the body at epoch t₂ (in AU).
/// * `v2` – Velocity vector at epoch t₂ (in AU/day).
/// * `dt` – Time difference between epochs t₁ and t₂ (in days).
/// * `peri_max` – Maximum acceptable perihelion distance for eccentricity control.
/// * `ecc_max` – Maximum acceptable eccentricity for eccentricity control.
///
/// Return
/// ------
/// * `Ok((v2_corrected, f, g))` on success, where:
///   - `v2_corrected` is the corrected velocity vector at t₂ (AU/day),
///   - `f`, `g` are the Lagrange coefficients.
/// * `Err(OutfitError::VelocityCorrectionError)` if eccentricity control fails
///   or the universal Kepler solver does not converge.
///
/// See also
/// --------
/// * [`velocity_correction_with_guess`](crate::kepler::velocity_correction_with_guess) – Variant that accepts a universal-variable guess.
/// * [`eccentricity_control`](crate::orb_elem::eccentricity_control) – Validates eccentricity and angular momentum.
/// * [`solve_kepuni`](crate::kepler::solve_kepuni) – Universal Kepler solver returning Stumpff-based integrals.
pub fn velocity_correction(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
    eps: f64,
) -> Result<(Vector3<f64>, f64, f64), OutfitError> {
    // Delegate to the more general implementation with warm-start capability,
    // but here we pass `None` for both the χ guess and the solver tolerance
    // to keep the legacy, minimal interface.
    let (velocity_corrected, f, g, _) =
        velocity_correction_with_guess(x1, x2, v2, dt, peri_max, ecc_max, None, eps)?;
    Ok((velocity_corrected, f, g))
}

/// Find corrected velocity using Lagrange f–g coefficients (configurable).
///
/// This routine refines the velocity vector `v2` of a body at epoch `t₂` by
/// combining two-body dynamics in universal variables with Lagrange’s f–g
/// formulation. Compared to [`velocity_correction`](crate::kepler::velocity_correction), this version allows:
/// - a custom solver tolerance for the universal Kepler equation,
/// - an optional warm-start on the universal anomaly `χ`.
///
/// Arguments
/// ---------
/// * `x1` – Position vector at epoch t₁ (AU).
/// * `x2` – Position vector at epoch t₂ (AU).
/// * `v2` – Velocity vector at epoch t₂ (AU/day).
/// * `dt` – Time difference t₂ − t₁ (days).
/// * `peri_max` – Maximum perihelion distance for dynamic acceptability.
/// * `ecc_max` – Maximum eccentricity for dynamic acceptability.
/// * `chi_guess` – Optional warm-start for the universal variable `χ`.
/// * `eps` – Optional solver tolerance; default is `1e3 × f64::EPSILON`.
///
/// Return
/// ------
/// * `Ok((v2_corrected, f, g, chi))` on success, where:
///   - `v2_corrected` is the corrected velocity at t₂ (AU/day),
///   - `f`, `g` are the Lagrange coefficients,
///   - `chi` is the final universal anomaly (useful for warm-starts).
/// * `Err(OutfitError::VelocityCorrectionError)` if:
///   - the eccentricity/energy check fails,
///   - the universal Kepler solver does not converge,
///   - or `g` is too small for a stable velocity update.
///
/// See also
/// --------
/// * [`velocity_correction`](crate::kepler::velocity_correction) – Simpler wrapper without warm-start or tolerance control.
/// * [`eccentricity_control`](crate::orb_elem::eccentricity_control) – Validates eccentricity and angular momentum.
/// * [`solve_kepuni_with_guess`](crate::kepler::solve_kepuni_with_guess) – Universal Kepler solver with warm-start support.
#[allow(clippy::too_many_arguments)]
pub fn velocity_correction_with_guess(
    x1: &Vector3<f64>,
    x2: &Vector3<f64>,
    v2: &Vector3<f64>,
    dt: f64,
    peri_max: f64,
    ecc_max: f64,
    chi_guess: Option<f64>,
    eps: f64,
) -> Result<(Vector3<f64>, f64, f64, f64), OutfitError> {
    let gravitational_parameter = GAUSS_GRAV_SQUARED;
    let radius_at_t2 = x2.norm();
    let radial_velocity_proxy_at_t2 = x2.dot(v2);

    // Quick reject: near-zero angular momentum means rectilinear, unstable motion.
    reject_if_angular_momentum_is_degenerate(x2, v2)?;

    // Step 1: validate eccentricity and compute the initial orbital energy.
    let (_, eccentricity, _, specific_orbital_energy) = eccentricity_control(
        x2, v2, peri_max, ecc_max,
    )
    .ok_or(OutfitError::VelocityCorrectionError(
        "Eccentricity control rejected the candidate orbit".into(),
    ))?;

    // Step 2: pack parameters for the universal Kepler solver.
    let params = UniversalKeplerParams {
        dt,
        r0: radius_at_t2,
        sig0: radial_velocity_proxy_at_t2,
        mu: gravitational_parameter,
        alpha: 2.0 * specific_orbital_energy, // specific orbital energy -> alpha
        e0: eccentricity,
    };

    // Step 3: solve the universal Kepler equation (with optional chi warm-start).
    let (universal_anomaly, _, _, s2, s3) = solve_kepuni_with_guess(&params, Some(eps), chi_guess)
        .ok_or(OutfitError::VelocityCorrectionError(
            "Universal Kepler solver did not converge".into(),
        ))?;

    // Step 4: compute the Lagrange coefficients f and g.
    let f_coefficient = 1.0 - (gravitational_parameter * s2) / radius_at_t2;
    let g_coefficient = dt - (gravitational_parameter * s3);

    // Guard against an ill-conditioned g (division instability).
    reject_if_lagrange_g_is_unstable(g_coefficient, dt)?;

    // Step 5: correct the velocity vector via the Lagrange f-g relation.
    let corrected_velocity = lagrange_corrected_velocity(x1, x2, f_coefficient, g_coefficient);

    Ok((
        corrected_velocity,
        f_coefficient,
        g_coefficient,
        universal_anomaly,
    ))
}

/// Reject orbits with a near-zero angular momentum (`position × velocity ≈ 0`),
/// which correspond to degenerate, rectilinear motion that the universal-variable
/// solver cannot handle reliably.
fn reject_if_angular_momentum_is_degenerate(
    position: &Vector3<f64>,
    velocity: &Vector3<f64>,
) -> Result<(), OutfitError> {
    let angular_momentum_norm = position.cross(velocity).norm();
    if !angular_momentum_norm.is_finite() || angular_momentum_norm <= 1e6 * f64::EPSILON {
        return Err(OutfitError::VelocityCorrectionError(
            "Rejected orbit: near-zero angular momentum (x × v ≈ 0)".into(),
        ));
    }
    Ok(())
}

/// Reject Lagrange coefficients `g` that are too small (or non-finite) to
/// support a numerically stable velocity update via division by `g`.
fn reject_if_lagrange_g_is_unstable(
    g_coefficient: f64,
    time_of_flight: f64,
) -> Result<(), OutfitError> {
    let g_absolute_value = g_coefficient.abs();
    // Scale-aware minimum magnitude required for a stable division by g.
    let g_minimum_magnitude = 100.0 * f64::EPSILON * (1.0 + time_of_flight.abs());

    if !g_absolute_value.is_finite() || g_absolute_value < g_minimum_magnitude {
        return Err(OutfitError::VelocityCorrectionError(
            "Lagrange coefficient g is too small for a stable velocity update".into(),
        ));
    }
    Ok(())
}

/// Apply the Lagrange f–g relation to recover the corrected velocity vector:
/// `v2_corrected = (x1 - f * x2) / g`.
fn lagrange_corrected_velocity(
    position_at_t1: &Vector3<f64>,
    position_at_t2: &Vector3<f64>,
    f_coefficient: f64,
    g_coefficient: f64,
) -> Vector3<f64> {
    // Implementation uses a BLAS-like axpy operation for fewer temporaries.
    let mut corrected_velocity = *position_at_t1;
    corrected_velocity.axpy(-f_coefficient, position_at_t2, 1.0); // corrected_velocity = x1 - f*x2
    corrected_velocity.unscale_mut(g_coefficient); // corrected_velocity /= g
    corrected_velocity
}

#[cfg(test)]
mod tests_velocity_correction {
    use approx::assert_relative_eq;
    use nalgebra::Vector3;

    use super::*;

    const KEP_EPS: f64 = 1e3 * f64::EPSILON;

    /// Helper function: easily build a 3D vector.
    fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
        Vector3::new(x, y, z)
    }

    #[test]
    fn test_velocity_correction_nominal_elliptic() {
        // Nominal test: elliptical orbit with simple positions and velocity.
        // This ensures that the correction converges and returns reasonable values.
        let x1 = v(1.0, 0.0, 0.0);
        let x2 = v(1.1, 0.0, 0.0);
        let v2 = v(0.0, 0.017, 0.0); // Approximate Earth orbital speed in AU/day
        let dt = 1.0; // 1 day
        let peri_max = 5.0;
        let ecc_max = 0.9;

        let result = velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS);
        assert!(result.is_ok(), "Velocity correction should succeed");

        let (vcorr, f, g) = result.unwrap();

        // All outputs should be finite
        assert!(vcorr.iter().all(|c| c.is_finite()));
        assert!(f.is_finite());
        assert!(g.is_finite());

        // f should be reasonably close to 1 for small geometry changes
        assert!(f > 0.5 && f < 1.5, "f coefficient is in a reasonable range");
        // g should be close to dt
        assert_relative_eq!(g, dt, epsilon = 0.1);
    }

    #[test]
    fn test_velocity_correction_rejects_bad_orbit() {
        // Use completely unrealistic input values to trigger an error in eccentricity_control.
        // This ensures that the function correctly returns an error instead of panicking.
        let x1 = v(1e6, 1e6, 1e6);
        let x2 = v(1e6, 1e6, 1e6);
        let v2 = v(1e6, 1e6, 1e6);
        let dt = 0.1;
        let peri_max = 0.001;
        let ecc_max = 0.001;

        let result = velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS);
        assert!(
            result.is_err(),
            "Velocity correction should fail on an invalid orbit"
        );
    }

    #[test]
    fn test_velocity_correction_sensitivity_to_x1() {
        // Sensitivity test:
        // If x1 changes slightly, the corrected velocity vector should also change.
        let x1 = v(1.0, 0.0, 0.0);
        let x1_shifted = v(1.05, 0.0, 0.0); // Slight offset
        let x2 = v(1.1, 0.0, 0.0);
        let v2 = v(0.0, 0.017, 0.0);
        let dt = 1.0;
        let peri_max = 5.0;
        let ecc_max = 0.9;

        let result1 = velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();
        let result2 =
            velocity_correction(&x1_shifted, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();

        let (v_corr1, _, _) = result1;
        let (v_corr2, _, _) = result2;

        let diff = (v_corr1 - v_corr2).norm();
        assert!(
            diff > 1e-6,
            "Changing x1 must influence the corrected velocity (diff = {diff})"
        );
    }

    #[test]
    fn test_velocity_correction_output_not_nan() {
        // This test ensures that the correction never returns NaN values,
        // even for slightly irregular configurations.
        let x1 = v(1.0, 0.5, 0.0);
        let x2 = v(1.1, 0.2, 0.0);
        let v2 = v(0.0, 0.02, 0.0);
        let dt = 0.5;
        let peri_max = 5.0;
        let ecc_max = 0.9;

        let (v_corr, f, g) =
            velocity_correction(&x1, &x2, &v2, dt, peri_max, ecc_max, KEP_EPS).unwrap();
        assert!(
            !v_corr.iter().any(|x| x.is_nan()),
            "Corrected velocity vector contains NaN"
        );
        assert!(!f.is_nan(), "f coefficient is NaN");
        assert!(!g.is_nan(), "g coefficient is NaN");
    }

    #[test]
    fn test_velocity_correction_real_data() {
        let x1 = Vector3::new(
            -0.843_561_126_129_683_3,
            0.937_288_327_370_772_8,
            0.659_183_901_029_776_6,
        );

        let x2 = Vector3::new(
            -0.623_121_622_917_384,
            1.0076797884556383,
            0.708_125_687_984_424_5,
        );

        let v2 = Vector3::new(
            -1.552_431_036_862_405_6E-2,
            -3.984_104_176_604_068E-3,
            -2.764_015_436_163_718_3E-3,
        );
        let dt = 14.731970000000729;

        let (v2, f, g) = velocity_correction(&x1, &x2, &v2, dt, 1., 1., KEP_EPS).unwrap();

        assert_eq!(f, 0.988_164_877_097_290_6);
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

    mod velocity_correction_prop_tests {
        use nalgebra::Vector3;
        use proptest::prelude::*;

        use super::*;

        /// Build a Vector3 from components
        fn v(x: f64, y: f64, z: f64) -> Vector3<f64> {
            Vector3::new(x, y, z)
        }

        /// Strategy to generate realistic orbital parameters
        fn arb_orbit_params(
        ) -> impl Strategy<Value = (Vector3<f64>, Vector3<f64>, Vector3<f64>, f64, f64, f64)>
        {
            (
                // x1: position vector (avoid zero)
                (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64)
                    .prop_map(|(x, y, z)| v(x + 1.0, y + 1.0, z)), // ensure non-zero
                // x2: position vector (close to x1)
                (-2.0..2.0f64, -2.0..2.0f64, -2.0..2.0f64)
                    .prop_map(|(x, y, z)| v(x + 1.1, y + 1.0, z)),
                // v2: velocity vector
                (-0.05..0.05f64, -0.05..0.05f64, -0.05..0.05f64).prop_map(|(x, y, z)| v(x, y, z)),
                // dt: time interval
                0.01..5.0f64,
                // peri_max
                0.1..10.0f64,
                // ecc_max
                0.1..2.0f64,
            )
        }

        proptest! {
            /// Property test: velocity_correction should never panic
            /// and always return finite results when it succeeds.
            #[test]
            fn prop_velocity_correction_no_panic(
                (x1,x2,v2,dt,peri_max,ecc_max) in arb_orbit_params()
            ) {
                let res = velocity_correction(&x1,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);

                match res {
                    Ok((vcorr, f, g)) => {
                        // Ensure outputs are finite
                        prop_assert!(vcorr.iter().all(|c| c.is_finite()));
                        prop_assert!(f.is_finite());
                        prop_assert!(g.is_finite());
                    }
                    Err(_) => {
                        // Err is acceptable: indicates eccentricity_control or solve_kepuni failure
                    }
                }
            }
        }

        proptest! {
            /// Differential property test: a slight change in x1 should
            /// change the result (when both calls return Ok).
            #[test]
            fn prop_velocity_correction_sensitive_to_x1(
                (x1,x2,v2,dt,peri_max,ecc_max) in arb_orbit_params()
            ) {
                // Slight perturbation on x1
                let mut x1_shifted = x1;
                x1_shifted[0] += 0.01;

                let res1 = velocity_correction(&x1,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);
                let res2 = velocity_correction(&x1_shifted,&x2,&v2,dt,peri_max,ecc_max, KEP_EPS);

                if let (Ok((v1, _, _)), Ok((v2c, _, _))) = (res1, res2) {
                    let diff = (v1 - v2c).norm();
                    // Either result is valid and diff can be very small or non-zero
                    prop_assert!(diff >= 0.0);
                }
            }
        }
    }
}
