//! Compartmental pharmacokinetics for orally dosed, site-of-release-controlled
//! prodrugs.
//!
//! A model is a [`graph::CompartmentGraph`]: compartments that admit species,
//! first-order transit between them, release kernels converting a prodrug into
//! its product where the chemistry says it should, absorption into a systemic
//! compartment, and elimination from it. A [`dose::Schedule`] says what is
//! given and when. [`solver::Solver`] integrates the two.
//!
//! ```
//! use bateman_core::prelude::*;
//! use std::collections::BTreeMap;
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
//!
//! let trajectory = Solver::default().run(&graph, &schedule, &ModelParameters::new())?;
//! assert!(trajectory.diagnostics.relative_error < 1e-6);
//! # Ok::<(), Box<dyn std::error::Error>>(())
//! ```
//!
//! # Conventions
//!
//! Every rate is per hour and every time is in hours. The state vector is in
//! nanomoles; concentrations are derived at output time by dividing by a
//! compartment volume. Release is one-to-one in molar terms, so one nanomole of
//! prodrug becomes one nanomole of product.
//!
//! Nothing leaves the state vector. Faecal loss and systemic elimination both
//! move material into terminal compartments, which is what makes the
//! mass-balance check in [`solver::MassBalanceReport`] a real check on the
//! integration rather than an accounting identity.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

pub mod dose;
pub mod error;
pub mod graph;
pub mod kernel;
pub mod params;
pub mod presets;
pub mod solver;
pub mod units;

/// The types most models need.
pub mod prelude {
    pub use crate::dose::{Dose, DoseSpec, Schedule};
    pub use crate::error::{GraphError, SolverError};
    pub use crate::graph::{CompartmentGraph, CompartmentId, SpeciesId};
    pub use crate::kernel::{
        ChargeDependent, EnzymeDriven, KernelCtx, ReleaseKernel, UserSupplied,
    };
    pub use crate::params::ModelParameters;
    pub use crate::presets;
    pub use crate::solver::{MassBalanceReport, Solver, Trajectory};
    pub use crate::units::{AmountNmol, BodyWeightG, ConcentrationUM, RateHrInv, TimeHr, VolumeML};
}
