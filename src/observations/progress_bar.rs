//! Lightweight iteration timing utilities.
//!
//! This module provides helpers to measure and report iteration times
//! in long-running loops, e.g. when combined with a progress bar
//! (see the `progress` feature).
//!
//! Components
//! -----------------
//! * [`IterTimer`] – Tracks per-iteration durations and computes a
//!   smoothed **exponential moving average** (EMA).
//!   Useful to get a stable estimate of iteration time even when
//!   individual steps fluctuate.
//!
//! * [`fmt_dur`] – Human-readable formatter for [`Duration`] values,
//!   producing strings like `"253µs"`, `"42ms"`, or `"3.14s"` depending
//!   on the scale.
//!
//! Usage
//! -----------------
//! Typical workflow inside a loop:
//!
//! ```rust, no_run
//! use std::time::Duration;
//! use your_crate::iter_timer::{IterTimer, fmt_dur};
//!
//! let mut timer = IterTimer::new(0.2); // smoothing factor α = 0.2
//!
//! for i in 0..10 {
//!     // ... some expensive work ...
//!
//!     let dt = timer.tick();
//!     println!(
//!         "iter {i} took {}, EMA = {}",
//!         fmt_dur(dt),
//!         fmt_dur(timer.avg())
//!     );
//! }
//! ```
//!
//! Design notes
//! -----------------
//! * The EMA update rule is:  
//!   `ema ← α·dt + (1–α)·ema`  
//!   with `α ∈ (0,1]`.  
//!   - `α = 1.0` → no smoothing (EMA = last sample).  
//!   - small `α` → stronger smoothing, slower adaptation.
//!
//! * [`IterTimer::tick`] must be called at each iteration boundary.
//!   The first tick initializes the average to the first duration.
//!
//! * [`IterTimer::avg`] returns the smoothed duration as a [`Duration`].
//!
//! * This module is enabled only with the `progress` feature.
#[cfg(feature = "progress")]
use std::time::{Duration, Instant};

pub struct IterTimer {
    last: Instant,
    ema_ns: f64,
    alpha: f64,
    count: u64,
}

impl IterTimer {
    pub fn new(alpha: f64) -> Self {
        Self {
            last: Instant::now(),
            ema_ns: 0.0,
            alpha,
            count: 0,
        }
    }

    #[inline]
    pub fn tick(&mut self) -> Duration {
        let now = Instant::now();
        let dt = now.duration_since(self.last);
        self.last = now;
        self.count += 1;

        let dt_ns = dt.as_nanos() as f64;
        self.ema_ns = if self.count == 1 {
            dt_ns
        } else {
            self.alpha * dt_ns + (1.0 - self.alpha) * self.ema_ns
        };

        dt
    }

    #[inline]
    pub fn avg(&self) -> Duration {
        if self.count == 0 {
            Duration::from_nanos(0)
        } else {
            Duration::from_nanos(self.ema_ns as u64)
        }
    }
}

#[inline]
pub fn fmt_dur(d: Duration) -> String {
    let us = d.as_micros();
    if us < 1_000 {
        format!("{us}µs")
    } else {
        let ms = d.as_millis();
        if ms < 1_000 {
            format!("{ms}ms")
        } else {
            let s = d.as_secs_f32();
            format!("{s:.2}s")
        }
    }
}
