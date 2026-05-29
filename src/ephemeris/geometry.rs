//! Geometric quantities derived from the apparent-position pipeline.
//!
//! This module sits inside the ephemeris pipeline and is responsible for
//! computing the five observational quantities stored in [`BodyGeometry`]:
//! phase angle, solar elongation, radial velocity, and the apparent angular
//! rates in right ascension and declination.  All quantities are derived from
//! the same [`PropagatedState`] that feeds [`ApparentPosition`], so the two
//! types share a single orbit propagation when computed together.
//!
//! # Role in the ephemeris pipeline
//!
//! ```text
//! OrbitalElements
//!       │
//!       ▼  propagate  (TwoBody or NBody)
//! PropagatedState  [equatorial mean J2000]
//!       │
//!       ├──────────────────────────────────────────────────┐
//!       ▼  apparent_position::assemble_apparent_position   ▼  geometry::compute_geometry
//! ApparentPosition                                   BodyGeometry
//! { coord, geocentric_dist,               { phase_angle, solar_elongation,
//!   heliocentric_dist }                     radial_velocity, d_ra_dt, d_dec_dt }
//! ```
//!
//! # Coordinate frame
//!
//! All vectors consumed by this module — body position, observer position, and
//! their velocities — are expressed in the **equatorial mean J2000** frame
//! (the same frame used for the final [`ApparentPosition`] output).
//! Positions are in **AU**, velocities in **AU/day**, and angles in **radians**.
//!
//! # Relationship to [`ApparentPosition`]
//!
//! Both types are computed from the same [`PropagatedState`].  When both are
//! needed, use [`OrbitalElements::apparent_position_and_geometry`] to avoid
//! propagating the orbit twice.  When only one is needed, use
//! [`OrbitalElements::apparent_position`] or
//! [`OrbitalElements::body_geometry`] respectively.
//!
//! # Bulk API
//!
//! Four methods on [`OrbitalElements`] expose this module at the public level:
//!
//! | Method | Time input | Returns |
//! |---|---|---|
//! | [`body_geometry_range`](OrbitalElements::body_geometry_range) | uniform range | [`GeometryTable`] |
//! | [`body_geometry_at`](OrbitalElements::body_geometry_at) | arbitrary epochs | [`GeometryTable`] |
//! | [`apparent_position_and_geometry_range`](OrbitalElements::apparent_position_and_geometry_range) | uniform range | [`ApparentPositionAndGeometryTable`] |
//! | [`apparent_position_and_geometry_at`](OrbitalElements::apparent_position_and_geometry_at) | arbitrary epochs | [`ApparentPositionAndGeometryTable`] |
//!
//! Errors are collected per epoch — a failure at one instant does not abort
//! computation for the remaining epochs.

use nalgebra::Vector3;

use crate::{EquinoctialElements, OutfitError};

use super::{apparent_position::PropagatedState, AberrationOrder};

// ---------------------------------------------------------------------------
// BodyGeometry
// ---------------------------------------------------------------------------

/// Geometric quantities for a solar system body at a given epoch.
///
/// Returned by [`crate::OrbitalElements::body_geometry`] and as the second element of
/// the tuple from [`crate::OrbitalElements::apparent_position_and_geometry`].
///
/// All angles are in **radians**, velocities in **AU/day**.
///
/// # Physical definitions
///
/// | Field | Definition |
/// |---|---|
/// | [`phase_angle`](Self::phase_angle) | Angle at the body between the Sun and the observer: Sun–body–observer |
/// | [`solar_elongation`](Self::solar_elongation) | Angle at the observer between the Sun and the body: Sun–observer–body |
/// | [`radial_velocity`](Self::radial_velocity) | Rate of change of the observer–body distance |
/// | [`d_ra_dt`](Self::d_ra_dt) | Apparent angular rate in right ascension |
/// | [`d_dec_dt`](Self::d_dec_dt) | Apparent angular rate in declination |
///
/// # Usage
///
/// The example below uses [`crate::OrbitalElements::body_geometry_at`] to compute
/// geometry at a set of discrete epochs and then reads all five fields:
///
/// ```rust,ignore
/// use hifitime::Epoch;
///
/// // `elements`, `observer`, `jpl`, `ut1`, `config` are assumed to be
/// // already constructed.
/// let times = vec![
///     Epoch::from_mjd_tt(60310.0),
///     Epoch::from_mjd_tt(60320.0),
///     Epoch::from_mjd_tt(60330.0),
/// ];
///
/// let table = elements.body_geometry_at(times, &observer, &jpl, &ut1, &config);
///
/// for (epoch, result) in table {
///     if let Ok(geo) = result {
///         println!(
///             "{epoch}: \
///              phase = {:.4} rad, \
///              elong = {:.4} rad, \
///              rv = {:.6} AU/day, \
///              dRA/dt = {:.6} rad/day, \
///              dDec/dt = {:.6} rad/day",
///             geo.phase_angle,
///             geo.solar_elongation,
///             geo.radial_velocity,
///             geo.d_ra_dt,
///             geo.d_dec_dt,
///         );
///     }
/// }
/// ```
///
/// When both sky coordinates and geometry are needed for the same epochs,
/// prefer [`crate::OrbitalElements::apparent_position_and_geometry_at`] to avoid
/// propagating the orbit twice:
///
/// ```rust,ignore
/// let table = elements.apparent_position_and_geometry_at(
///     times, &observer, &jpl, &ut1, &config,
/// );
/// for (epoch, result) in table {
///     if let Ok((pos, geo)) = result {
///         println!("{epoch}: RA = {:.4}, phase = {:.4}", pos.coord.ra, geo.phase_angle);
///     }
/// }
/// ```
///
/// # See also
///
/// - [`crate::OrbitalElements::body_geometry`] — single-epoch entry point.
/// - [`crate::OrbitalElements::apparent_position_and_geometry`] — combined single-epoch
///   entry point (one propagation, two outputs).
/// - [`crate::ephemeris::GeometryTable`] — bulk return type for geometry-only queries.
/// - [`crate::ephemeris::ApparentPositionAndGeometryTable`] — bulk return type for combined queries.
#[derive(Debug, Clone, PartialEq)]
pub struct BodyGeometry {
    /// Phase angle — Sun–body–observer \[rad\], $\phi \in [0, \pi]$.
    ///
    /// $$\phi = \arccos\left(\frac{\mathbf{r}\_\text{body} \cdot \mathbf{d}}{r\_\text{helio}\\,\rho}\right)$$
    ///
    /// where $\mathbf{d}$ is the aberration-corrected topocentric vector
    /// and $\rho = |\mathbf{d}|$.
    ///
    /// - $\phi = 0$: body is in opposition (fully illuminated as seen by the observer).
    /// - $\phi = \pi$: body is at superior conjunction (dark side facing the observer).
    pub phase_angle: f64,

    /// Solar elongation — Sun–observer–body \[rad\], $\varepsilon \in [0, \pi]$.
    ///
    /// $$\varepsilon = \arccos\left(\frac{-\mathbf{r}\_\text{obs} \cdot \mathbf{d}}{|\mathbf{r}\_\text{obs}|\\,\rho}\right)$$
    ///
    /// Small values indicate the body is close to the Sun on the sky and may
    /// be difficult to observe.
    pub solar_elongation: f64,

    /// Observer-relative radial velocity \[AU/day\].
    ///
    /// $$\dot{\rho} = \frac{\mathbf{d} \cdot \mathbf{v}\_\text{topo}}{\rho}$$
    ///
    /// where $\mathbf{v}\_\text{topo} = \mathbf{v}\_\text{body} - \mathbf{v}\_\text{obs}$
    /// is the true topocentric velocity (observer velocity from
    /// [`ResolvedObserver::pvobs`](crate::observer_extension::ResolvedObserver::pvobs)).
    ///
    /// - Positive: body is receding from the observer.
    /// - Negative: body is approaching the observer.
    pub radial_velocity: f64,

    /// Apparent angular rate in right ascension \[rad/day\].
    ///
    /// $$\dot{\alpha} = \frac{\partial\alpha}{\partial\mathbf{r}} \cdot \mathbf{v}\_\text{topo}$$
    ///
    /// Positive eastward.  Includes the true topocentric velocity of the observer.
    pub d_ra_dt: f64,

    /// Apparent angular rate in declination \[rad/day\].
    ///
    /// $$\dot{\delta} = \frac{\partial\delta}{\partial\mathbf{r}} \cdot \mathbf{v}\_\text{topo}$$
    ///
    /// Positive northward.
    pub d_dec_dt: f64,
}

// ---------------------------------------------------------------------------
// Computation
// ---------------------------------------------------------------------------

/// Compute the [`BodyGeometry`] from a [`PropagatedState`].
///
/// This function is the single entry point for geometry computation inside the
/// ephemeris pipeline.  It executes the following steps in order:
///
/// 1. **Aberration correction** — builds the raw topocentric vector
///    $\mathbf{d}_0 = \mathbf{r}_\text{body} - \mathbf{r}_\text{obs}$
///    and applies either the first-order or the second-order stellar-aberration
///    model (controlled by `aberration`) to obtain the corrected line-of-sight
///    $\mathbf{d}$.
/// 2. **Norms** — computes the topocentric distance $\rho = |\mathbf{d}|$,
///    the heliocentric distance $r_\text{helio} = |\mathbf{r}_\text{body}|$,
///    and the observer–Sun distance $r_\text{obs} = |\mathbf{r}_\text{obs}|$.
/// 3. **Phase angle** $\phi$ — via [`phase_angle`] (Sun–body–observer,
///    dot product clamped before `acos`).
/// 4. **Solar elongation** $\varepsilon$ — via [`solar_elongation`]
///    (Sun–observer–body, same clamping strategy).
/// 5. **True topocentric velocity** —
///    $\mathbf{v}_\text{topo} = \mathbf{v}_\text{body} - \mathbf{v}_\text{obs}$.
/// 6. **Radial velocity** $\dot{\rho}$ — via [`radial_velocity`], the
///    projection of $\mathbf{v}_\text{topo}$ onto the line of sight.
/// 7. **Angular rates** $(\dot{\alpha}, \dot{\delta})$ — via
///    [`angular_rates`], using the geometric Jacobian of the spherical-coordinate
///    transform.
///
/// All input vectors are expected in the **equatorial mean J2000** frame with
/// positions in AU and velocities in AU/day.  The only additional cost over
/// computing [`ApparentPosition`] alone is one `acos` per angle and a few
/// dot products.
///
/// # Arguments
///
/// - `state`      – Propagated body + observer state at the observation epoch.
/// - `elements`   – Orbital elements (needed by second-order aberration).
/// - `aberration` – Aberration correction order.
///
/// # Panics
///
/// This function does not panic under any input.  Potential division-by-zero
/// conditions (zero topocentric distance, body on the celestial pole) are
/// guarded by [`radial_velocity`] and [`angular_rates`].
///
/// # Errors
///
/// Returns [`OutfitError`] if the second-order aberration back-propagation
/// fails.  The first-order path is always infallible.
pub(crate) fn compute_geometry(
    state: &PropagatedState,
    elements: &EquinoctialElements,
    aberration: &AberrationOrder,
) -> Result<BodyGeometry, OutfitError> {
    // Aberration-corrected topocentric vector and its norm.
    let topo_raw = state.ast_pos_equ - state.obs_pos_equ;
    let topo = aberration.apply_aberration(
        topo_raw,
        state.ast_vel_equ,
        elements,
        state.obs_time_mjd,
        state.obs_pos_equ,
    )?;
    let rho = topo.norm();

    let r_helio = state.ast_pos_equ.norm();
    let r_obs = state.obs_pos_equ.norm();

    let phase_angle = phase_angle(state.ast_pos_equ, topo, r_helio, rho);
    let solar_elongation = solar_elongation(state.obs_pos_equ, topo, r_obs, rho);

    // True topocentric velocity: body velocity minus observer velocity.
    let v_topo = state.ast_vel_equ - state.obs_vel_equ;

    let radial_velocity = radial_velocity(topo, v_topo, rho);
    let (d_ra_dt, d_dec_dt) = angular_rates(topo, v_topo, rho);

    Ok(BodyGeometry {
        phase_angle,
        solar_elongation,
        radial_velocity,
        d_ra_dt,
        d_dec_dt,
    })
}

// ---------------------------------------------------------------------------
// Private helpers
// ---------------------------------------------------------------------------

/// Compute the phase angle (Sun–body–observer) \[rad\].
///
/// The phase angle $\phi$ is the angle at the body between the direction
/// toward the Sun and the direction toward the observer.  It determines the
/// fraction of the illuminated disk visible to the observer: $\phi = 0$
/// corresponds to full illumination (opposition), $\phi = \pi$ to a dark
/// face (superior conjunction).
///
/// $$\phi = \arccos\left(\frac{\mathbf{r}_\text{body} \cdot \mathbf{d}}{r_\text{helio}\,\rho}\right)$$
///
/// The cosine argument is clamped to $[-1, 1]$ before calling `acos` to
/// guard against floating-point rounding that would otherwise produce `NaN`
/// when the body is exactly at opposition or conjunction.
#[inline]
fn phase_angle(ast_pos: Vector3<f64>, topo: Vector3<f64>, r_helio: f64, rho: f64) -> f64 {
    let cos_phi = (ast_pos.dot(&topo) / (r_helio * rho)).clamp(-1.0, 1.0);
    cos_phi.acos()
}

/// Compute the solar elongation (Sun–observer–body) \[rad\].
///
/// The solar elongation $\varepsilon$ is the angle at the observer
/// between the direction toward the Sun and the direction toward the body.
/// Small values ($\varepsilon \lesssim 20°$) indicate that the body is
/// close to the Sun on the sky and may be unobservable from the ground.
///
/// $$\varepsilon = \arccos\left(\frac{-\mathbf{r}_\text{obs} \cdot \mathbf{d}}{|\mathbf{r}_\text{obs}|\,\rho}\right)$$
///
/// The negation of $\mathbf{r}_\text{obs}$ converts from the
/// observer-to-Sun direction to the Sun-to-observer direction, giving the
/// correct sign convention.  As with [`phase_angle`], the cosine argument is
/// clamped to $[-1, 1]$ to prevent `NaN` from rounding.
#[inline]
fn solar_elongation(obs_pos: Vector3<f64>, topo: Vector3<f64>, r_obs: f64, rho: f64) -> f64 {
    let cos_eps = (-obs_pos.dot(&topo) / (r_obs * rho)).clamp(-1.0, 1.0);
    cos_eps.acos()
}

/// Compute the observer-relative radial velocity \[AU/day\].
///
/// The radial velocity $\dot{\rho}$ is the rate of change of the
/// topocentric distance: positive when the body is receding from the observer,
/// negative when it is approaching.  It is the component of the true
/// topocentric velocity $\mathbf{v}_\text{topo}$ projected onto the
/// unit line-of-sight vector.
///
/// $$\dot{\rho} = \frac{\mathbf{d} \cdot \mathbf{v}_\text{topo}}{\rho}$$
///
/// No special edge-case handling is needed here: $\rho = 0$ would mean
/// the observer is at the body, which is physically impossible in practice.
#[inline]
fn radial_velocity(topo: Vector3<f64>, v_topo: Vector3<f64>, rho: f64) -> f64 {
    topo.dot(&v_topo) / rho
}

/// Compute the apparent angular rates $(\dot{\alpha}, \dot{\delta})$
/// \[rad/day\].
///
/// The angular rates give how fast the body moves across the sky as seen by
/// the observer: $\dot{\alpha}$ is positive eastward and $\dot{\delta}$
/// is positive northward.  Both include the contribution from the observer's
/// own velocity (diurnal and annual motion).
///
/// Uses the geometric Jacobians
///
/// $$\frac{\partial\alpha}{\partial\mathbf{r}} = \frac{1}{d_x^2+d_y^2}\begin{pmatrix}-d_y\\d_x\\0\end{pmatrix}$$
///
/// $$\frac{\partial\delta}{\partial\mathbf{r}} = \frac{1}{\rho^2\sqrt{d_x^2+d_y^2}}\begin{pmatrix}-d_z d_x\\-d_z d_y\\d_x^2+d_y^2\end{pmatrix}$$
///
/// where $\mathbf{d} = (d_x, d_y, d_z)$ is the aberration-corrected
/// topocentric vector.
///
/// **Degenerate pole case**: when the body lies on (or very close to) the
/// celestial pole, $d_x^2 + d_y^2 \approx 0$ and the Jacobians are
/// singular.  The function detects this condition — specifically when
/// $\sqrt{d_x^2+d_y^2} < \varepsilon_\text{machine} \cdot \rho$ — and
/// returns $(0, 0)$ in that case.
///
/// Returns $(\dot{\alpha}, \dot{\delta})$ in rad/day.
fn angular_rates(topo: Vector3<f64>, v_topo: Vector3<f64>, rho: f64) -> (f64, f64) {
    let dx = topo[0];
    let dy = topo[1];
    let dz = topo[2];

    let dxy2 = dx * dx + dy * dy; // ‖(dx, dy)‖²
    let dxy = dxy2.sqrt();

    // Guard against degenerate case (body on the pole, dxy ≈ 0).
    if dxy < f64::EPSILON * rho {
        return (0.0, 0.0);
    }

    // ∂α/∂r · v_topo
    let d_ra_dt = (-dy * v_topo[0] + dx * v_topo[1]) / dxy2;

    // ∂δ/∂r · v_topo
    let d_dec_dt =
        (-dz * dx * v_topo[0] - dz * dy * v_topo[1] + dxy2 * v_topo[2]) / (rho * rho * dxy);

    (d_ra_dt, d_dec_dt)
}
