//! Exposure and response on top of a [`bateman_core`] trajectory.
//!
//! A pharmacodynamic statement is two steps: reduce a concentration-time
//! profile to a scalar exposure, then map that exposure to an effect.
//! [`PdLink`] composes the two so the choice of summary and the choice of
//! response curve stay separable.
//!
//! ```
//! use bateman_core::prelude::*;
//! use bateman_pd::{Auc, Hill, PdLink};
//!
//! let graph = presets::murine_gut()?;
//! let stomach = graph.compartment_id("stomach").unwrap();
//! let schedule = Schedule::single_gavage(
//!     DoseSpec::NmolAbsolute(1000.0),
//!     BodyWeightG(25.0),
//!     stomach,
//!     presets::even_prodrug_mix(&graph),
//!     TimeHr(24.0),
//! );
//! let trajectory = Solver::default().run(&graph, &schedule, &ModelParameters::new())?;
//!
//! let plasma = graph.compartment_id("plasma").unwrap();
//! let product = graph.species_id("free_product").unwrap();
//! let index = graph.index_of(plasma, product).unwrap();
//! let volume = graph.compartment("plasma").unwrap().volume;
//!
//! let link = PdLink::new(Auc, Hill { emax: 1.0, ec50: 20.0, n: 2.0 });
//! let effect = link.apply(&trajectory, index, volume);
//! assert!((0.0..=1.0).contains(&effect));
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod exposure;
pub mod response;

pub use exposure::{Auc, Cavg, Cmax, ExposureIntegrator, TimeAboveThreshold};
pub use response::{Emax, Hill, Linear, ResponseFn};

use bateman_core::solver::Trajectory;
use bateman_core::units::VolumeML;

/// An exposure summary composed with a response curve.
#[derive(Debug, Clone, Copy)]
pub struct PdLink<E, R> {
    /// How the profile is reduced to a scalar.
    pub exposure: E,
    /// How that scalar becomes an effect.
    pub response: R,
}

impl<E: ExposureIntegrator, R: ResponseFn> PdLink<E, R> {
    /// Compose an exposure integrator with a response function.
    pub fn new(exposure: E, response: R) -> Self {
        Self { exposure, response }
    }

    /// The exposure of one state entry, as a concentration summary.
    pub fn exposure_of(&self, trajectory: &Trajectory, index: usize, volume: VolumeML) -> f64 {
        let c = trajectory.concentration_series(index, volume);
        self.exposure.apply(&trajectory.t, &c)
    }

    /// The response driven by that exposure.
    pub fn apply(&self, trajectory: &Trajectory, index: usize, volume: VolumeML) -> f64 {
        self.response
            .response(self.exposure_of(trajectory, index, volume))
    }
}
