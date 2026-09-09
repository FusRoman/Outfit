//! Time-dependent perturber positions for the N-body propagator.
//!
//! The DOP853 right-hand side (`NBodyOde`) needs the heliocentric position of
//! every perturbing body at the *current* integration time, but querying the JPL
//! ephemeris on every call to `diff` is far too slow and cannot report errors.
//! This module pre-computes, once per observation arc, a **piecewise Chebyshev
//! interpolation** of each perturber's heliocentric position (ecliptic J2000,
//! AU) that `diff` can evaluate allocation-free and infallibly.
//!
//! ## Design
//!
//! - The validity window is split into contiguous, equal-length **panels**. Each
//!   panel carries one truncated Chebyshev series per Cartesian axis, fitted from
//!   samples taken at the Chebyshev–Gauss nodes of that panel.
//! - Panelling keeps the required degree low and independent of the arc length,
//!   so fast bodies (Mercury) and long backward arcs are handled without tuning.
//! - The **Sun** is a special case: its heliocentric position is identically
//!   zero, so it carries no series and `position_at` returns exactly
//!   `Vector3::zeros()`, bit-identical to the previous frozen snapshot behaviour.
//!
//! Everything that can fail (missing GM, ephemeris lookup) happens in the `build`
//! constructors, which return [`OutfitError`]; the hot-path accessors never fail.

use std::f64::consts::PI;

use hifitime::{Epoch, TimeScale};
use nalgebra::Vector3;

use crate::jpl_ephem::naif::naif_ids::{solar_system_bary::SolarSystemBary, NaifIds};
use crate::jpl_ephem::JPLEphem;
use crate::outfit_errors::OutfitError;
use crate::propagator::{planet_gm::gm_au3_day2, NBodyConfig};

// ───────────────────────────────────────────────────────────────────────────
// Pure numerical helpers (no I/O, unit- and property-testable in isolation)
// ───────────────────────────────────────────────────────────────────────────

/// Returns the `degree + 1` Chebyshev–Gauss nodes on `[-1, 1]`.
///
/// The nodes are `x_k = cos(π (k + ½) / (degree + 1))` for `k = 0 ..= degree`;
/// they lie strictly inside `(-1, 1)`, are symmetric about `0`, and are listed
/// in strictly decreasing order.
///
/// # Arguments
///
/// * `degree` – degree of the Chebyshev series the nodes will be used to fit.
///
/// # Returns
///
/// A `Vec<f64>` of length `degree + 1` holding the nodes.
fn chebyshev_gauss_nodes(degree: usize) -> Vec<f64> {
    let n_plus_1 = degree + 1;
    (0..n_plus_1)
        .map(|k| (PI * (k as f64 + 0.5) / n_plus_1 as f64).cos())
        .collect()
}

/// Fits a truncated Chebyshev series to values sampled at the Chebyshev–Gauss
/// nodes.
///
/// Uses the discrete orthogonality of `T_j` on those nodes:
/// `c_j = (2 / (N+1)) Σ_k f_k cos(π j (k + ½) / (N+1))`. The returned vector
/// stores `c_0 / 2` in slot `0` and `c_j` elsewhere, so the interpolant is the
/// plain sum `Σ_j coeffs[j] T_j(x)` evaluated by [`clenshaw`] without any
/// special first-term handling.
///
/// # Arguments
///
/// * `samples_at_nodes` – values `f(x_k)`, in the same order as
///   [`chebyshev_gauss_nodes`]; length must equal `degree + 1`.
/// * `degree` – degree of the fitted series.
///
/// # Returns
///
/// A `Vec<f64>` of `degree + 1` Chebyshev coefficients (first term pre-halved).
fn fit_chebyshev_series(samples_at_nodes: &[f64], degree: usize) -> Vec<f64> {
    let n_plus_1 = degree + 1;
    debug_assert_eq!(samples_at_nodes.len(), n_plus_1);

    let scale = 2.0 / n_plus_1 as f64;
    let mut coeffs = vec![0.0_f64; n_plus_1];
    for (j, coeff) in coeffs.iter_mut().enumerate() {
        let acc: f64 = samples_at_nodes
            .iter()
            .enumerate()
            .map(|(k, &f_k)| f_k * (PI * j as f64 * (k as f64 + 0.5) / n_plus_1 as f64).cos())
            .sum();
        *coeff = scale * acc;
    }
    coeffs[0] *= 0.5;
    coeffs
}

/// Evaluates `Σ_j coeffs[j] T_j(tau)` by the Clenshaw recurrence.
///
/// `coeffs[0]` is used as-is (callers such as [`fit_chebyshev_series`] pre-halve
/// the constant term). The routine performs no heap allocation.
///
/// # Arguments
///
/// * `coeffs` – Chebyshev coefficients, lowest order first.
/// * `tau` – evaluation point, expected in `[-1, 1]`.
///
/// # Returns
///
/// The value of the Chebyshev sum at `tau`.
#[inline]
fn clenshaw(coeffs: &[f64], tau: f64) -> f64 {
    match coeffs {
        [] => 0.0,
        [c0] => *c0,
        [c0, rest @ ..] => {
            let mut b_kp1 = 0.0_f64;
            let mut b_kp2 = 0.0_f64;
            for &c in rest.iter().rev() {
                let b_k = 2.0 * tau * b_kp1 - b_kp2 + c;
                b_kp2 = b_kp1;
                b_kp1 = b_k;
            }
            tau * b_kp1 - b_kp2 + c0
        }
    }
}

/// Maps a time `t` onto the panel-local coordinate `(t - mid) / radius`.
///
/// # Arguments
///
/// * `t` – time to map (MJD TT).
/// * `mid` – panel midpoint (MJD TT).
/// * `radius` – panel half-length in days; must be strictly positive.
///
/// # Returns
///
/// The normalised coordinate; it lies in `[-1, 1]` when `t` is inside the panel.
#[inline]
fn map_to_unit(t: f64, mid: f64, radius: f64) -> f64 {
    (t - mid) / radius
}

/// Computes the uniform panel decomposition of a time window.
///
/// The window `[t_start, t_end]` is split into `ceil((t_end - t_start) /
/// panel_days)` equal panels (at least one). A non-finite or non-positive
/// `panel_days` falls back to a single panel spanning the whole window.
///
/// # Arguments
///
/// * `t_start` – window start (MJD TT); must be strictly less than `t_end`.
/// * `t_end` – window end (MJD TT).
/// * `panel_days` – target maximum panel length in days.
///
/// # Returns
///
/// `(n_panels, radius)` where `radius` is the common panel half-length in days.
fn panel_layout(t_start: f64, t_end: f64, panel_days: f64) -> (usize, f64) {
    let span = t_end - t_start;
    let target = if panel_days.is_finite() && panel_days > 0.0 {
        panel_days
    } else {
        span
    };
    let n_panels = (span / target).ceil().max(1.0) as usize;
    let radius = span / (2.0 * n_panels as f64);
    (n_panels, radius)
}

/// Returns the midpoint (MJD TT) of panel `index` in a uniform decomposition.
///
/// # Arguments
///
/// * `t_start` – window start (MJD TT).
/// * `radius` – common panel half-length in days.
/// * `index` – zero-based panel index.
///
/// # Returns
///
/// `t_start + radius · (2·index + 1)`.
#[inline]
fn panel_mid(t_start: f64, radius: f64, index: usize) -> f64 {
    t_start + radius * (2 * index + 1) as f64
}

/// Selects the panel index whose interval contains `t`, clamped to the valid
/// range.
///
/// Times before the window map to panel `0` and times after it to the last
/// panel, so a small extrapolation margin (e.g. a rejected DOP853 step) stays
/// well-defined.
///
/// # Arguments
///
/// * `t` – query time (MJD TT).
/// * `t_start` – window start (MJD TT).
/// * `radius` – common panel half-length in days; must be strictly positive.
/// * `n_panels` – number of panels; must be at least `1`.
///
/// # Returns
///
/// A panel index in `0 ..= n_panels - 1`.
#[inline]
fn select_panel_index(t: f64, t_start: f64, radius: f64, n_panels: usize) -> usize {
    let raw = ((t - t_start) / (2.0 * radius)).floor();
    if raw < 0.0 {
        0
    } else {
        (raw as usize).min(n_panels - 1)
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Ephemeris sampling (the only I/O in this module)
// ───────────────────────────────────────────────────────────────────────────

/// Samples one perturber's heliocentric position at each of `node_epochs`.
///
/// Only the position is read; `JPLEphem::body_ephemeris` velocities are not used
/// here (they carry a backend-dependent scaling that is irrelevant to a
/// position-only Chebyshev fit).
///
/// # Arguments
///
/// * `body` – perturbing body to sample.
/// * `jpl` – opened JPL ephemeris.
/// * `node_epochs` – epochs at which to sample, one per Chebyshev–Gauss node.
///
/// # Returns
///
/// `[xs, ys, zs]`, three `Vec<f64>` of AU components in ecliptic J2000, aligned
/// with `node_epochs`.
///
/// # Errors
///
/// Propagates any [`OutfitError`] returned by [`JPLEphem::body_ephemeris`].
fn sample_component_grids(
    body: NaifIds,
    jpl: &JPLEphem,
    node_epochs: &[Epoch],
) -> Result<[Vec<f64>; 3], OutfitError> {
    let mut xs = Vec::with_capacity(node_epochs.len());
    let mut ys = Vec::with_capacity(node_epochs.len());
    let mut zs = Vec::with_capacity(node_epochs.len());
    for epoch in node_epochs {
        let (position, _velocity) = jpl.body_ephemeris(body, epoch)?;
        xs.push(position[0]);
        ys.push(position[1]);
        zs.push(position[2]);
    }
    Ok([xs, ys, zs])
}

// ───────────────────────────────────────────────────────────────────────────
// Public types
// ───────────────────────────────────────────────────────────────────────────

/// One Chebyshev panel: a truncated series per Cartesian axis valid on
/// `[mid - radius, mid + radius]` (MJD TT), returning AU in ecliptic J2000.
#[derive(Debug, Clone)]
struct ChebPanel {
    /// Panel midpoint (MJD TT).
    mid: f64,
    /// Panel half-length (days); strictly positive.
    radius: f64,
    /// Chebyshev coefficients for the `x`, `y`, `z` components (first term
    /// pre-halved, see [`fit_chebyshev_series`]).
    coeffs: [Vec<f64>; 3],
}

/// Piecewise Chebyshev interpolation of a single perturber's heliocentric
/// position over a time window, plus its gravitational parameter.
#[derive(Debug, Clone)]
pub(crate) struct PerturberEphemeris {
    /// Contiguous panels ordered by increasing `mid`. Empty when `is_sun`.
    panels: Vec<ChebPanel>,
    /// Gravitational parameter GM (AU³/day²).
    gravitational_parameter: f64,
    /// `true` for the Sun, whose heliocentric position is identically zero.
    is_sun: bool,
}

impl PerturberEphemeris {
    /// Builds the interpolation for `body` over `[t_start_mjd_tt, t_end_mjd_tt]`.
    ///
    /// For the Sun no sampling is done. For any other body the window is split
    /// into panels (see [`panel_layout`]) and each panel's per-axis series is
    /// fitted from `degree + 1` ephemeris samples taken at the Chebyshev–Gauss
    /// nodes.
    ///
    /// # Arguments
    ///
    /// * `body` – perturbing body.
    /// * `jpl` – opened JPL ephemeris.
    /// * `t_start_mjd_tt` / `t_end_mjd_tt` – validity window (MJD TT), `t_start < t_end`.
    /// * `degree` – Chebyshev degree per panel; values below `1` are treated as `1`.
    /// * `panel_days` – target maximum panel length in days.
    ///
    /// # Returns
    ///
    /// A ready-to-evaluate [`PerturberEphemeris`].
    ///
    /// # Errors
    ///
    /// - [`OutfitError::EphemerisBodyNotSupported`] if `body` has no tabulated GM.
    /// - Any [`OutfitError`] from [`JPLEphem::body_ephemeris`] while sampling.
    pub(crate) fn build(
        body: NaifIds,
        jpl: &JPLEphem,
        t_start_mjd_tt: f64,
        t_end_mjd_tt: f64,
        degree: usize,
        panel_days: f64,
    ) -> Result<Self, OutfitError> {
        let gravitational_parameter = gm_au3_day2(body).ok_or_else(|| {
            OutfitError::EphemerisBodyNotSupported(format!(
                "No GM available for perturber {body:?}"
            ))
        })?;

        if matches!(body, NaifIds::SSB(SolarSystemBary::Sun)) {
            return Ok(Self {
                panels: Vec::new(),
                gravitational_parameter,
                is_sun: true,
            });
        }

        let degree = degree.max(1);
        let (n_panels, radius) = panel_layout(t_start_mjd_tt, t_end_mjd_tt, panel_days);
        let nodes = chebyshev_gauss_nodes(degree);

        let mut panels = Vec::with_capacity(n_panels);
        for index in 0..n_panels {
            let mid = panel_mid(t_start_mjd_tt, radius, index);
            let node_epochs: Vec<Epoch> = nodes
                .iter()
                .map(|&x| Epoch::from_mjd_in_time_scale(mid + radius * x, TimeScale::TT))
                .collect();
            let grids = sample_component_grids(body, jpl, &node_epochs)?;
            panels.push(ChebPanel {
                mid,
                radius,
                coeffs: [
                    fit_chebyshev_series(&grids[0], degree),
                    fit_chebyshev_series(&grids[1], degree),
                    fit_chebyshev_series(&grids[2], degree),
                ],
            });
        }

        Ok(Self {
            panels,
            gravitational_parameter,
            is_sun: false,
        })
    }

    /// Evaluates the perturber's heliocentric position at `mjd_tt`.
    ///
    /// Infallible and allocation-free: this is called on every DOP853 stage.
    /// The Sun returns exactly `Vector3::zeros()`; other bodies select the
    /// enclosing panel and evaluate its series with the argument clamped to the
    /// panel interval.
    ///
    /// # Arguments
    ///
    /// * `mjd_tt` – evaluation epoch (MJD TT).
    ///
    /// # Returns
    ///
    /// The heliocentric position in AU, ecliptic J2000.
    #[inline]
    pub(crate) fn position_at(&self, mjd_tt: f64) -> Vector3<f64> {
        if self.is_sun {
            return Vector3::zeros();
        }
        let Some(first) = self.panels.first() else {
            debug_assert!(false, "non-Sun perturber must have at least one panel");
            return Vector3::zeros();
        };
        let radius = first.radius;
        let window_start = first.mid - radius;
        let index = select_panel_index(mjd_tt, window_start, radius, self.panels.len());
        let panel = &self.panels[index];
        let tau = map_to_unit(mjd_tt, panel.mid, panel.radius).clamp(-1.0, 1.0);
        Vector3::new(
            clenshaw(&panel.coeffs[0], tau),
            clenshaw(&panel.coeffs[1], tau),
            clenshaw(&panel.coeffs[2], tau),
        )
    }

    /// Returns the perturber's gravitational parameter GM (AU³/day²).
    #[inline]
    pub(crate) fn gm(&self) -> f64 {
        self.gravitational_parameter
    }
}

/// A full set of per-body position interpolations covering one observation arc,
/// built once and shared by every propagation of that trajectory.
///
/// It appears in the public `EquinoctialElements::propagate_nbody` and
/// `single_iteration` signatures so callers can hand in a pre-built table; the
/// differential-correction driver builds one per trajectory automatically.
#[derive(Debug, Clone)]
pub struct PerturberEphemerisSet {
    perturbers: Vec<PerturberEphemeris>,
    t_start: f64,
    t_end: f64,
}

impl PerturberEphemerisSet {
    /// Builds one per-body position interpolation for every entry of
    /// `config.perturbing_bodies` over `[t_start_mjd_tt, t_end_mjd_tt]`.
    ///
    /// # Arguments
    ///
    /// * `config` – N-body configuration (body list, interpolation degree, panel length).
    /// * `jpl` – opened JPL ephemeris.
    /// * `t_start_mjd_tt` / `t_end_mjd_tt` – validity window (MJD TT); `t_start < t_end`.
    ///
    /// # Returns
    ///
    /// A [`PerturberEphemerisSet`] whose entries match `config.perturbing_bodies`
    /// in order.
    ///
    /// # Errors
    ///
    /// - [`OutfitError::NBodyPropagationFailed`] if the window is empty or invalid.
    /// - [`OutfitError::EphemerisBodyNotSupported`] or any ephemeris-lookup error
    ///   raised while fitting one body's series.
    pub fn build(
        config: &NBodyConfig,
        jpl: &JPLEphem,
        t_start_mjd_tt: f64,
        t_end_mjd_tt: f64,
    ) -> Result<Self, OutfitError> {
        let span = t_end_mjd_tt - t_start_mjd_tt;
        if !span.is_finite() || span <= 0.0 {
            return Err(OutfitError::NBodyPropagationFailed(format!(
                "invalid perturber-ephemeris window: t_start={t_start_mjd_tt}, t_end={t_end_mjd_tt}"
            )));
        }

        let perturbers = config
            .perturbing_bodies
            .iter()
            .map(|&body| {
                PerturberEphemeris::build(
                    body,
                    jpl,
                    t_start_mjd_tt,
                    t_end_mjd_tt,
                    config.perturber_interp_degree,
                    config.perturber_panel_days,
                )
            })
            .collect::<Result<Vec<_>, _>>()?;

        Ok(Self {
            perturbers,
            t_start: t_start_mjd_tt,
            t_end: t_end_mjd_tt,
        })
    }

    /// Returns `true` if `[t_lo, t_hi]` is fully inside the built window.
    ///
    /// # Arguments
    ///
    /// * `t_lo` / `t_hi` – closed interval to test (MJD TT).
    ///
    /// # Returns
    ///
    /// Whether every point of `[t_lo, t_hi]` can be evaluated without extrapolation.
    #[inline]
    pub fn covers(&self, t_lo: f64, t_hi: f64) -> bool {
        t_lo >= self.t_start && t_hi <= self.t_end
    }

    /// Returns the per-body interpolators in `config.perturbing_bodies` order.
    #[inline]
    pub(crate) fn perturbers(&self) -> &[PerturberEphemeris] {
        &self.perturbers
    }
}

// ───────────────────────────────────────────────────────────────────────────
// Tests
// ───────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod perturber_ephemeris_tests {
    use super::*;
    use crate::jpl_ephem::naif::naif_ids::planet_bary::PlanetaryBary;
    use crate::test_fixture::JPL_EPHEM_HORIZON;
    use approx::assert_abs_diff_eq;
    use proptest::prelude::*;

    // ── Oracle for the Chebyshev basis (naive recurrence) ──────────────────

    /// Evaluates `Σ_j coeffs[j] T_j(tau)` directly from the `T_j` recurrence.
    fn naive_cheb_sum(coeffs: &[f64], tau: f64) -> f64 {
        if coeffs.is_empty() {
            return 0.0;
        }
        let mut t_prev = 1.0_f64; // T_0
        let mut acc = coeffs[0] * t_prev;
        if coeffs.len() == 1 {
            return acc;
        }
        let mut t_curr = tau; // T_1
        acc += coeffs[1] * t_curr;
        for &c in &coeffs[2..] {
            let t_next = 2.0 * tau * t_curr - t_prev;
            acc += c * t_next;
            t_prev = t_curr;
            t_curr = t_next;
        }
        acc
    }

    /// Horner evaluation of a monomial polynomial (highest order last).
    fn poly_eval(monomials: &[f64], x: f64) -> f64 {
        monomials.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    // ── Unit tests: pure helpers ──────────────────────────────────────────

    #[test]
    fn nodes_have_expected_shape() {
        let nodes = chebyshev_gauss_nodes(12);
        assert_eq!(nodes.len(), 13);
        for &x in &nodes {
            assert!(x > -1.0 && x < 1.0);
        }
        // strictly decreasing
        for pair in nodes.windows(2) {
            assert!(pair[0] > pair[1]);
        }
        // symmetric about 0
        let n = nodes.len();
        for k in 0..n {
            assert_abs_diff_eq!(nodes[k], -nodes[n - 1 - k], epsilon = 1e-14);
        }
    }

    #[test]
    fn clenshaw_matches_low_order_chebyshev() {
        // T_0..T_3 at a handful of points.
        for &tau in &[-1.0, -0.5, 0.0, 0.25, 1.0] {
            assert_abs_diff_eq!(clenshaw(&[1.0], tau), 1.0, epsilon = 1e-15);
            assert_abs_diff_eq!(clenshaw(&[0.0, 1.0], tau), tau, epsilon = 1e-15);
            assert_abs_diff_eq!(
                clenshaw(&[0.0, 0.0, 1.0], tau),
                2.0 * tau * tau - 1.0,
                epsilon = 1e-15
            );
            assert_abs_diff_eq!(
                clenshaw(&[0.0, 0.0, 0.0, 1.0], tau),
                4.0 * tau.powi(3) - 3.0 * tau,
                epsilon = 1e-14
            );
        }
    }

    #[test]
    fn fit_recovers_a_constant_with_halved_first_term() {
        let degree = 6;
        let samples = vec![3.5_f64; degree + 1];
        let coeffs = fit_chebyshev_series(&samples, degree);
        assert_abs_diff_eq!(coeffs[0], 3.5, epsilon = 1e-12);
        for c in &coeffs[1..] {
            assert_abs_diff_eq!(*c, 0.0, epsilon = 1e-12);
        }
    }

    #[test]
    fn map_to_unit_is_the_inverse_of_the_node_mapping() {
        let (mid, radius) = (59_000.0, 16.0);
        for &x in &[-1.0, -0.3, 0.0, 0.75, 1.0] {
            assert_abs_diff_eq!(
                map_to_unit(mid + radius * x, mid, radius),
                x,
                epsilon = 1e-12
            );
        }
    }

    #[test]
    fn panel_layout_covers_the_window() {
        let (t_start, t_end) = (59_000.0, 59_400.0);
        let (n_panels, radius) = panel_layout(t_start, t_end, 32.0);
        assert!(n_panels >= 1);
        assert_abs_diff_eq!(
            2.0 * radius * n_panels as f64,
            t_end - t_start,
            epsilon = 1e-9
        );
        // first panel starts at t_start, last ends at t_end
        assert_abs_diff_eq!(
            panel_mid(t_start, radius, 0) - radius,
            t_start,
            epsilon = 1e-9
        );
        assert_abs_diff_eq!(
            panel_mid(t_start, radius, n_panels - 1) + radius,
            t_end,
            epsilon = 1e-9
        );
    }

    #[test]
    fn select_panel_index_clamps_outside_the_window() {
        let (t_start, t_end) = (59_000.0, 59_400.0);
        let (n_panels, radius) = panel_layout(t_start, t_end, 32.0);
        assert_eq!(
            select_panel_index(t_start - 100.0, t_start, radius, n_panels),
            0
        );
        assert_eq!(
            select_panel_index(t_end + 100.0, t_start, radius, n_panels),
            n_panels - 1
        );
    }

    #[test]
    fn build_rejects_body_without_gm() {
        let err = PerturberEphemeris::build(
            NaifIds::SSB(SolarSystemBary::SSB),
            &JPL_EPHEM_HORIZON,
            59_000.0,
            59_100.0,
            12,
            32.0,
        )
        .unwrap_err();
        assert!(matches!(err, OutfitError::EphemerisBodyNotSupported(_)));
    }

    #[test]
    fn set_build_rejects_empty_window() {
        let cfg = NBodyConfig::default();
        let err =
            PerturberEphemerisSet::build(&cfg, &JPL_EPHEM_HORIZON, 59_100.0, 59_000.0).unwrap_err();
        assert!(matches!(err, OutfitError::NBodyPropagationFailed(_)));
    }

    #[test]
    fn covers_is_inclusive_at_the_bounds() {
        let cfg = NBodyConfig::default();
        let set =
            PerturberEphemerisSet::build(&cfg, &JPL_EPHEM_HORIZON, 59_000.0, 59_400.0).unwrap();
        assert!(set.covers(59_000.0, 59_400.0));
        assert!(set.covers(59_100.0, 59_200.0));
        assert!(!set.covers(58_999.0, 59_200.0));
        assert!(!set.covers(59_200.0, 59_401.0));
    }

    // ── Oracle test: interpolation vs the ephemeris itself ─────────────────

    fn max_position_error_over_grid(body: NaifIds, t0: f64, span: f64) -> f64 {
        let (t_start, t_end) = if span >= 0.0 {
            (t0, t0 + span)
        } else {
            (t0 + span, t0)
        };
        let cfg = NBodyConfig {
            perturbing_bodies: vec![body],
            ..NBodyConfig::default()
        };
        let set =
            PerturberEphemerisSet::build(&cfg, &JPL_EPHEM_HORIZON, t_start - 5.0, t_end + 5.0)
                .unwrap();
        let interp = &set.perturbers()[0];

        let mut max_err = 0.0_f64;
        let steps = 200;
        for i in 0..=steps {
            let t = t_start + (t_end - t_start) * (i as f64 / steps as f64);
            let epoch = Epoch::from_mjd_in_time_scale(t, TimeScale::TT);
            let (truth, _) = JPL_EPHEM_HORIZON.body_ephemeris(body, &epoch).unwrap();
            let err = (interp.position_at(t) - truth).norm();
            max_err = max_err.max(err);
        }
        max_err
    }

    #[test]
    fn interpolation_tracks_the_ephemeris_forward_and_backward() {
        let t0 = 60_000.0;
        for body in [
            NaifIds::PB(PlanetaryBary::Mercury),
            NaifIds::PB(PlanetaryBary::EarthMoon),
            NaifIds::PB(PlanetaryBary::Mars),
            NaifIds::PB(PlanetaryBary::Jupiter),
        ] {
            for span in [400.0_f64, -400.0_f64] {
                let err = max_position_error_over_grid(body, t0, span);
                assert!(
                    err < 1e-9,
                    "{body:?} span {span}: max position error {err:e} AU exceeds 1e-9"
                );
            }
        }
    }

    #[test]
    fn sun_position_is_exactly_zero() {
        let interp = PerturberEphemeris::build(
            NaifIds::SSB(SolarSystemBary::Sun),
            &JPL_EPHEM_HORIZON,
            59_000.0,
            59_400.0,
            12,
            32.0,
        )
        .unwrap();
        assert_eq!(interp.position_at(59_123.4), Vector3::zeros());
    }

    // ── Property-based tests ──────────────────────────────────────────────

    proptest! {
        /// Clenshaw agrees with the naive `T_j` recurrence.
        #[test]
        fn clenshaw_agrees_with_naive_sum(
            coeffs in prop::collection::vec(-10.0_f64..10.0, 1..20),
            tau in -1.0_f64..1.0,
        ) {
            prop_assert!((clenshaw(&coeffs, tau) - naive_cheb_sum(&coeffs, tau)).abs() < 1e-9);
        }

        /// Fitting then evaluating reproduces any polynomial of degree ≤ `degree`.
        #[test]
        fn fit_is_exact_for_low_degree_polynomials(
            monomials in prop::collection::vec(-5.0_f64..5.0, 1..9),
            x in -1.0_f64..1.0,
        ) {
            let degree = monomials.len().max(2) + 2; // headroom above the polynomial degree
            let nodes = chebyshev_gauss_nodes(degree);
            let samples: Vec<f64> = nodes.iter().map(|&n| poly_eval(&monomials, n)).collect();
            let coeffs = fit_chebyshev_series(&samples, degree);
            let got = clenshaw(&coeffs, x);
            let want = poly_eval(&monomials, x);
            prop_assert!((got - want).abs() < 1e-8, "got {got}, want {want}");
        }

        /// `map_to_unit` inverts the node placement for any positive radius.
        #[test]
        fn map_to_unit_round_trips(
            mid in 40_000.0_f64..70_000.0,
            radius in 0.5_f64..200.0,
            x in -1.0_f64..1.0,
        ) {
            prop_assert!((map_to_unit(mid + radius * x, mid, radius) - x).abs() < 1e-9);
        }

        /// Every time in the window lands in a panel whose interval contains it.
        #[test]
        fn panel_selection_is_a_covering_partition(
            span in 5.0_f64..800.0,
            panel_days in 1.0_f64..64.0,
            frac in 0.0_f64..1.0,
        ) {
            let t_start = 59_000.0;
            let t_end = t_start + span;
            let (n_panels, radius) = panel_layout(t_start, t_end, panel_days);
            let window_start = panel_mid(t_start, radius, 0) - radius;
            let t = t_start + frac * span;
            let idx = select_panel_index(t, window_start, radius, n_panels);
            prop_assert!(idx < n_panels);
            let mid = panel_mid(t_start, radius, idx);
            prop_assert!((t - mid).abs() <= radius + 1e-6);
        }
    }
}
