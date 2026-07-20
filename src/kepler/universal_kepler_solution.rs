/// Solution of the universal Kepler equation.
///
/// Returned by [`solve`](crate::kepler::UniversalKeplerParams::solve) on convergence, this struct bundles
/// the universal anomaly $\psi$ with the four Stumpff auxiliary values
/// $(s_0, s_1, s_2, s_3)$ evaluated at $(\psi, \alpha)$.
///
/// # Physical meaning
///
/// ## Universal anomaly $\psi$
///
/// The universal anomaly generalizes the classical orbital anomalies across all
/// conic regimes:
///
/// | Regime | Relation |
/// |---|---|
/// | Elliptic  | $\psi = \sqrt{a} \cdot E$ where $E$ is the eccentric anomaly |
/// | Hyperbolic | $\psi = \sqrt{-a} \cdot H$ where $H$ is the hyperbolic anomaly |
///
/// A single equation covers all regimes:
///
/// $$r_0 \cdot s_1(\psi, \alpha) + \sigma_0 \cdot s_2(\psi, \alpha) + s_3(\psi, \alpha) = \sqrt{\mu} \cdot \Delta t$$
///
/// where `alpha` is the reciprocal semi-major-axis convention
/// ( $\alpha = -1/a$ ), not the raw vis-viva $2E$ — see
/// [`UniversalKeplerParams`](crate::kepler::UniversalKeplerParams).
///
/// ## Stumpff auxiliary functions $s_n(\psi, \alpha)$
///
/// These functions unify trigonometric and hyperbolic functions through the
/// energy parameter $\alpha$ via the substitution $\beta = \alpha \psi^2$:
///
/// $$s_n(\psi, \alpha) = \psi^n \cdot c_n(\alpha \psi^2)$$
///
/// where $c_n$ are the Stumpff functions defined by:
///
/// $$c_n(x) = \sum_{k=0}^{\infty} \frac{(-x)^k}{(n + 2k)!}$$
///
/// Their explicit forms depend on the orbital regime ($\beta = \alpha\psi^2$):
///
/// | | Elliptic $(\beta < 0)$ | Hyperbolic $(\beta > 0)$ |
/// |---|---|---|
/// | $s_0$ | $\cos\!\sqrt{-\beta}$ | $\cosh\!\sqrt{\beta}$ |
/// | $s_1$ | $\sin\!\sqrt{-\beta}\,/\sqrt{-\alpha}$ | $\sinh\!\sqrt{\beta}\,/\sqrt{\alpha}$ |
/// | $s_2$ | $(1 - \cos\!\sqrt{-\beta})\,/(-\alpha)$ | $(\cosh\!\sqrt{\beta}-1)\,/\alpha$ |
/// | $s_3$ | $(\sqrt{-\beta}-\sin\!\sqrt{-\beta})\,/(-\alpha)^{3/2}$ | $(\sinh\!\sqrt{\beta}-\sqrt{\beta})\,/\alpha^{3/2}$ |
///
/// The differentiation chain with respect to $\psi$:
///
/// $$\frac{d s_3}{d\psi} = s_2, \qquad \frac{d s_2}{d\psi} = s_1, \qquad \frac{d s_1}{d\psi} = s_0$$
///
/// which directly yields the residual derivative:
///
/// $$f'(\psi) = r_0 \cdot s_0 + \sigma_0 \cdot s_1 + s_2$$
///
/// ## Role of each $s_n$ in orbit propagation
///
/// The Stumpff functions appear in the **Lagrange coefficients** $f$ and $g$,
/// which propagate position and velocity from epoch $t_0$ to $t_0 + \Delta t$:
///
/// $$\vec{r}(t) = f \cdot \vec{r}_0 + g \cdot \vec{v}_0$$
/// $$\vec{v}(t) = \dot{f} \cdot \vec{r}_0 + \dot{g} \cdot \vec{v}_0$$
///
/// | Field | Stumpff analog | Role in propagation |
/// |---|---|---|
/// | `cos_like` ($s_0$) | $\cos$ / $\cosh$ | Lagrange coefficient $f$ and $\dot{g}$ |
/// | `sin_like` ($s_1$) | $\sin/\omega$ / $\sinh/\omega$ | Lagrange coefficient $g$ (units of time) |
/// | `one_minus_cos_like` ($s_2$) | $(1-\cos)/\alpha$ | Lagrange coefficient $\dot{f}$ |
/// | `time_integral` ($s_3$) | time-of-flight kernel | Carries $\Delta t$ in the Kepler equation |
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UniversalKeplerSolution {
    /// Universal anomaly $\psi$ (rad, generalized units).
    ///
    /// The root of the universal Kepler equation:
    ///
    /// $$f(\psi) = r_0 \cdot s_1 + \sigma_0 \cdot s_2 + s_3 - \sqrt{\mu} \cdot \Delta t = 0$$
    pub universal_anomaly: f64,

    /// $s_0(\psi, \alpha)$ — cosine-like Stumpff function.
    ///
    /// Reduces to $\cos\!\sqrt{-\alpha}\,\psi$ (elliptic) or
    /// $\cosh\!\sqrt{\alpha}\,\psi$ (hyperbolic).
    ///
    /// Appears in the Lagrange coefficients $f$ and $\dot{g}$.
    pub cos_like: f64,

    /// $s_1(\psi, \alpha)$ — sine-like Stumpff function (units: $\text{length}^{1/2}$).
    ///
    /// Reduces to $\sin\!\sqrt{-\alpha}\,\psi\,/\sqrt{-\alpha}$ (elliptic) or
    /// $\sinh\!\sqrt{\alpha}\,\psi\,/\sqrt{\alpha}$ (hyperbolic).
    ///
    /// Appears in the Lagrange coefficient $g$.
    pub sin_like: f64,

    /// $s_2(\psi, \alpha)$ — $(1 - \cos)$-like Stumpff function (units: $\text{length}$).
    ///
    /// Reduces to $(1 - \cos\!\sqrt{-\alpha}\,\psi)\,/(-\alpha)$ (elliptic) or
    /// $(\cosh\!\sqrt{\alpha}\,\psi - 1)\,/\alpha$ (hyperbolic).
    ///
    /// Appears in the Lagrange coefficient $\dot{f}$.
    pub one_minus_cos_like: f64,

    /// $s_3(\psi, \alpha)$ — time-of-flight kernel Stumpff function
    /// (units: $\text{length}^{3/2} \cdot \text{time}^{-1}$).
    ///
    /// Reduces to $(\sqrt{-\alpha}\,\psi - \sin\!\sqrt{-\alpha}\,\psi)\,/(-\alpha)^{3/2}$
    /// (elliptic) or $(\sinh\!\sqrt{\alpha}\,\psi - \sqrt{\alpha}\,\psi)\,/\alpha^{3/2}$
    /// (hyperbolic).
    ///
    /// This is the dominant term carrying the time of flight $\Delta t$ in the
    /// Kepler equation.
    pub time_integral: f64,
}

impl UniversalKeplerSolution {
    /// Construct a [`UniversalKeplerSolution`] from the raw tuple returned by
    /// [`s_funct`](crate::kepler::s_funct).
    ///
    /// Arguments
    /// ---------
    /// * `psi` – Converged universal anomaly.
    /// * `s` – Tuple `(s0, s1, s2, s3)` from `s_funct(psi, alpha)`.
    #[inline]
    pub fn from_raw(psi: f64, s: (f64, f64, f64, f64)) -> Self {
        Self {
            universal_anomaly: psi,
            cos_like: s.0,
            sin_like: s.1,
            one_minus_cos_like: s.2,
            time_integral: s.3,
        }
    }

    /// Return the raw tuple `(s0, s1, s2, s3)` for use in downstream
    /// computations that expect the legacy positional form.
    #[inline]
    pub fn as_raw_stumpff(&self) -> (f64, f64, f64, f64) {
        (
            self.cos_like,
            self.sin_like,
            self.one_minus_cos_like,
            self.time_integral,
        )
    }
}
