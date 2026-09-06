//! Errors raised while building a graph or running a solve.

use crate::graph::CompartmentId;

/// A graph that cannot describe a well-posed model.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum GraphError {
    /// Two compartments share a name.
    #[error("duplicate compartment name: {0}")]
    DuplicateCompartment(
        /// The repeated name.
        String,
    ),
    /// Two species share a name.
    #[error("duplicate species name: {0}")]
    DuplicateSpecies(
        /// The repeated name.
        String,
    ),
    /// An edge names a compartment that was never declared.
    #[error("unknown compartment: {0}")]
    UnknownCompartment(
        /// The name that matched no compartment.
        String,
    ),
    /// An edge names a species that was never declared.
    #[error("unknown species: {0}")]
    UnknownSpecies(
        /// The name that matched no species.
        String,
    ),
    /// An edge moves a species through a compartment that does not admit it.
    #[error("compartment {0:?} does not admit the species on this edge")]
    SpeciesNotInCompartment(
        /// The compartment that does not admit it.
        CompartmentId,
    ),
    /// No compartment was marked as a dose site.
    #[error("the graph has no dose site")]
    NoDoseSite,
    /// A compartment volume is zero or negative, so no concentration is defined.
    #[error("compartment {0} has non-positive volume {1}")]
    NonPositiveVolume(
        /// Compartment name.
        String,
        /// The volume given.
        f64,
    ),
    /// A rate constant is negative.
    #[error("rate on edge {0} is negative: {1}")]
    NegativeRate(
        /// Canonical name of the edge.
        String,
        /// The rate given.
        f64,
    ),
}

/// A solve that cannot be completed.
#[derive(Debug, thiserror::Error)]
pub enum SolverError {
    /// The integrator failed on an interval.
    #[error("integration failed between t={t0} and t={t1}: {message}")]
    Integration {
        /// Start of the interval.
        t0: f64,
        /// End of the interval.
        t1: f64,
        /// What the integrator reported.
        message: String,
    },
    /// A dose is scheduled before zero or past the horizon.
    #[error("a dose at t={0} falls outside the schedule horizon")]
    DoseOutsideHorizon(
        /// The offending dose time, in hours.
        f64,
    ),
    /// The schedule horizon is zero or negative.
    #[error("schedule horizon must be positive, got {0}")]
    NonPositiveHorizon(
        /// The horizon that was given.
        f64,
    ),
    /// An override names a parameter the graph does not have.
    #[error("unknown parameter override: {0}")]
    UnknownParameter(
        /// The name that matched nothing.
        String,
    ),
    /// The integration lost or gained material beyond tolerance.
    #[error("mass balance violated: relative error {0:e} exceeds tolerance {1:e}")]
    MassBalance(
        /// Relative error observed.
        f64,
        /// Tolerance it exceeded.
        f64,
    ),
}
