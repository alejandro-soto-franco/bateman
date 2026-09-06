//! Exposure integrators: scalar summaries of a concentration-time profile.

use bateman_core::units::TimeHr;

/// A scalar summary of a concentration-time profile.
pub trait ExposureIntegrator {
    /// Reduce a profile to one number. `t` and `c` are the same length and `t`
    /// ascends.
    fn apply(&self, t: &[TimeHr], c: &[f64]) -> f64;
}

/// Area under the concentration-time curve, by the trapezoid rule.
#[derive(Debug, Clone, Copy, Default)]
pub struct Auc;

impl ExposureIntegrator for Auc {
    fn apply(&self, t: &[TimeHr], c: &[f64]) -> f64 {
        t.windows(2)
            .zip(c.windows(2))
            .map(|(ts, cs)| 0.5 * (cs[0] + cs[1]) * (ts[1].0 - ts[0].0))
            .sum()
    }
}

/// Peak concentration over the profile.
#[derive(Debug, Clone, Copy, Default)]
pub struct Cmax;

impl ExposureIntegrator for Cmax {
    fn apply(&self, _t: &[TimeHr], c: &[f64]) -> f64 {
        c.iter().copied().fold(f64::NEG_INFINITY, f64::max)
    }
}

/// Average concentration: the area under the curve divided by its span.
#[derive(Debug, Clone, Copy, Default)]
pub struct Cavg;

impl ExposureIntegrator for Cavg {
    fn apply(&self, t: &[TimeHr], c: &[f64]) -> f64 {
        match (t.first(), t.last()) {
            (Some(a), Some(b)) if b.0 > a.0 => Auc.apply(t, c) / (b.0 - a.0),
            _ => 0.0,
        }
    }
}

/// Time spent at or above a threshold concentration, in hours.
///
/// Crossings are located by linear interpolation between output times, so the
/// answer does not quantise to the output grid.
#[derive(Debug, Clone, Copy)]
pub struct TimeAboveThreshold(
    /// The threshold concentration.
    pub f64,
);

impl ExposureIntegrator for TimeAboveThreshold {
    fn apply(&self, t: &[TimeHr], c: &[f64]) -> f64 {
        let thr = self.0;
        let mut total = 0.0;
        for (ts, cs) in t.windows(2).zip(c.windows(2)) {
            let (t0, t1) = (ts[0].0, ts[1].0);
            let (c0, c1) = (cs[0], cs[1]);
            let dt = t1 - t0;
            match (c0 >= thr, c1 >= thr) {
                (true, true) => total += dt,
                (false, false) => {}
                // One crossing inside the interval: interpolate to find it.
                (true, false) => total += dt * (c0 - thr) / (c0 - c1),
                (false, true) => total += dt * (c1 - thr) / (c1 - c0),
            }
        }
        total
    }
}
