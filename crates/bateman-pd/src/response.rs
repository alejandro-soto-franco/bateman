//! Response functions: exposure to effect.

/// A map from a scalar exposure to a scalar response.
pub trait ResponseFn {
    /// The response at this exposure.
    fn response(&self, exposure: f64) -> f64;
}

/// Sigmoid Hill response with an adjustable coefficient.
///
/// `n = 1` is the Emax model. A negative or zero exposure gives zero response,
/// since the Hill form is undefined below zero for a fractional coefficient.
#[derive(Debug, Clone, Copy)]
pub struct Hill {
    /// Maximum attainable response.
    pub emax: f64,
    /// Exposure producing half of `emax`.
    pub ec50: f64,
    /// Hill coefficient.
    pub n: f64,
}

impl ResponseFn for Hill {
    fn response(&self, exposure: f64) -> f64 {
        if exposure <= 0.0 {
            return 0.0;
        }
        let e = exposure.powf(self.n);
        let half = self.ec50.powf(self.n);
        self.emax * e / (half + e)
    }
}

/// Hyperbolic Emax response, the Hill model at unit coefficient.
#[derive(Debug, Clone, Copy)]
pub struct Emax {
    /// Maximum attainable response.
    pub emax: f64,
    /// Exposure producing half of `emax`.
    pub ec50: f64,
}

impl ResponseFn for Emax {
    fn response(&self, exposure: f64) -> f64 {
        Hill {
            emax: self.emax,
            ec50: self.ec50,
            n: 1.0,
        }
        .response(exposure)
    }
}

/// Straight-line response.
#[derive(Debug, Clone, Copy)]
pub struct Linear {
    /// Response per unit exposure.
    pub slope: f64,
    /// Response at zero exposure.
    pub intercept: f64,
}

impl ResponseFn for Linear {
    fn response(&self, exposure: f64) -> f64 {
        self.slope * exposure + self.intercept
    }
}
