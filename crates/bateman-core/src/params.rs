//! Named parameter overrides applied to a graph at solve time.
//!
//! A [`crate::graph::CompartmentGraph`] has a default value on every edge, so a graph
//! alone is enough to simulate. Fitting needs to vary a subset of those values
//! without rebuilding the topology, so every edge also has a canonical name and
//! [`ModelParameters`] maps a name to a replacement value.
//!
//! The names are:
//!
//! | Edge | Name |
//! |---|---|
//! | transit | `k:{from}->{to}:{species}` |
//! | absorption | `ka:{from}->{to}:{species}` |
//! | elimination | `ke:{compartment}:{species}` |
//! | release kernel | `kernel:{compartment}:{field}` |
//!
//! [`ModelParameters::names_in`] lists every name a given graph accepts, so a
//! caller can discover the parameter space rather than guess at it.

use std::collections::BTreeMap;

/// Replacement values for named edges, applied at the start of a solve.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ModelParameters {
    values: BTreeMap<String, f64>,
}

impl ModelParameters {
    /// An empty set of overrides: every edge keeps its graph default.
    pub fn new() -> Self {
        Self::default()
    }

    /// Override one named edge.
    pub fn set(&mut self, name: impl Into<String>, value: f64) -> &mut Self {
        self.values.insert(name.into(), value);
        self
    }

    /// Builder form of [`ModelParameters::set`].
    pub fn with(mut self, name: impl Into<String>, value: f64) -> Self {
        self.set(name, value);
        self
    }

    /// The override for `name`, when one is set.
    pub fn get(&self, name: &str) -> Option<f64> {
        self.values.get(name).copied()
    }

    /// The override for `name`, falling back to `default`.
    pub fn get_or(&self, name: &str, default: f64) -> f64 {
        self.get(name).unwrap_or(default)
    }

    /// Every name this set overrides, in sorted order.
    pub fn names(&self) -> impl Iterator<Item = &str> {
        self.values.keys().map(String::as_str)
    }

    /// True when no override is set.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }
}
