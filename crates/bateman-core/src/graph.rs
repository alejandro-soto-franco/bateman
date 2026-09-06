//! The compartment graph: compartments, the species each admits, and the edges
//! along which material moves.
//!
//! A graph is built through [`GraphBuilder`] and validated on [`GraphBuilder::build`],
//! so an invalid topology is rejected before any integration runs.

use std::collections::BTreeMap;

use crate::error::GraphError;
use crate::kernel::ReleaseKernel;
use crate::params::ModelParameters;
use crate::units::{RateHrInv, VolumeML};

/// Index of a compartment within its graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CompartmentId(
    /// Position of the compartment in its graph's compartment list.
    pub u32,
);

/// Index of a species within its graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpeciesId(
    /// Position of the species in its graph's species list.
    pub u32,
);

/// A well-mixed volume admitting a fixed set of species.
#[derive(Debug, Clone)]
pub struct Compartment {
    /// This compartment's index in its graph.
    pub id: CompartmentId,
    /// Name used in override strings and in output.
    pub name: String,
    /// Species this compartment admits. Anything else is rejected at build time.
    pub species: Vec<SpeciesId>,
    /// Volume used to turn an amount into a concentration.
    pub volume: VolumeML,
}

/// First-order movement of one species between two compartments.
#[derive(Debug, Clone)]
pub struct Transfer {
    /// Source compartment.
    pub from: CompartmentId,
    /// Destination compartment.
    pub to: CompartmentId,
    /// The species that moves.
    pub species: SpeciesId,
    /// First-order rate constant.
    pub rate: RateHrInv,
}

/// Conversion of a prodrug into its product, in place, at a kernel-set rate.
#[derive(Debug)]
pub struct ReleaseSite {
    /// Where the conversion happens.
    pub compartment: CompartmentId,
    /// The species consumed.
    pub prodrug: SpeciesId,
    /// The species produced, one for one in molar terms.
    pub product: SpeciesId,
    /// The rate constant, which may depend on time, site and state.
    pub kernel: Box<dyn ReleaseKernel>,
}

/// First-order absorption of one species into a systemic compartment.
#[derive(Debug, Clone)]
pub struct AbsorptionSink {
    /// The compartment absorbed from.
    pub from_compartment: CompartmentId,
    /// The systemic compartment absorbed into.
    pub to_compartment: CompartmentId,
    /// The species absorbed.
    pub species: SpeciesId,
    /// First-order absorption rate constant.
    pub rate: RateHrInv,
}

/// First-order removal of one species from the system.
#[derive(Debug, Clone)]
pub struct Elimination {
    /// The compartment eliminated from.
    pub compartment: CompartmentId,
    /// The species eliminated.
    pub species: SpeciesId,
    /// First-order elimination rate constant.
    pub rate: RateHrInv,
    /// Terminal compartment the removed material is accounted into, so nothing
    /// leaves the state vector and mass balance stays checkable.
    pub sink: CompartmentId,
}

/// A validated topology.
#[derive(Debug)]
pub struct CompartmentGraph {
    compartments: Vec<Compartment>,
    species_names: Vec<String>,
    transfers: Vec<Transfer>,
    releases: Vec<ReleaseSite>,
    absorptions: Vec<AbsorptionSink>,
    eliminations: Vec<Elimination>,
    dose_sites: Vec<CompartmentId>,
    /// State index for each (compartment, species) pair that exists.
    index: BTreeMap<(CompartmentId, SpeciesId), usize>,
    n_state: usize,
}

impl CompartmentGraph {
    /// Start building a graph.
    pub fn builder() -> GraphBuilder {
        GraphBuilder::default()
    }

    /// Number of entries in the state vector.
    pub fn state_len(&self) -> usize {
        self.n_state
    }

    /// State index of a (compartment, species) pair, when the compartment
    /// admits that species.
    pub fn index_of(&self, c: CompartmentId, s: SpeciesId) -> Option<usize> {
        self.index.get(&(c, s)).copied()
    }

    /// Look a compartment up by name.
    pub fn compartment(&self, name: &str) -> Option<&Compartment> {
        self.compartments.iter().find(|c| c.name == name)
    }

    /// Look a compartment id up by name.
    pub fn compartment_id(&self, name: &str) -> Option<CompartmentId> {
        self.compartment(name).map(|c| c.id)
    }

    /// Look a species id up by name.
    pub fn species_id(&self, name: &str) -> Option<SpeciesId> {
        self.species_names
            .iter()
            .position(|n| n == name)
            .map(|i| SpeciesId(i as u32))
    }

    /// The name of a compartment.
    pub fn compartment_name(&self, c: CompartmentId) -> &str {
        &self.compartments[c.0 as usize].name
    }

    /// The name of a species.
    pub fn species_name(&self, s: SpeciesId) -> &str {
        &self.species_names[s.0 as usize]
    }

    /// Every compartment, in id order.
    pub fn compartments(&self) -> &[Compartment] {
        &self.compartments
    }

    pub(crate) fn transfers(&self) -> &[Transfer] {
        &self.transfers
    }

    pub(crate) fn releases(&self) -> &[ReleaseSite] {
        &self.releases
    }

    pub(crate) fn absorptions(&self) -> &[AbsorptionSink] {
        &self.absorptions
    }

    pub(crate) fn eliminations(&self) -> &[Elimination] {
        &self.eliminations
    }

    /// The compartments a dose may be placed in.
    pub fn dose_sites(&self) -> &[CompartmentId] {
        &self.dose_sites
    }

    /// Canonical override name for a transit edge.
    pub fn transit_name(&self, t: &Transfer) -> String {
        format!(
            "k:{}->{}:{}",
            self.compartment_name(t.from),
            self.compartment_name(t.to),
            self.species_name(t.species)
        )
    }

    /// Canonical override name for an absorption edge.
    pub fn absorption_name(&self, a: &AbsorptionSink) -> String {
        format!(
            "ka:{}->{}:{}",
            self.compartment_name(a.from_compartment),
            self.compartment_name(a.to_compartment),
            self.species_name(a.species)
        )
    }

    /// Canonical override name for an elimination edge.
    pub fn elimination_name(&self, e: &Elimination) -> String {
        format!(
            "ke:{}:{}",
            self.compartment_name(e.compartment),
            self.species_name(e.species)
        )
    }

    /// Canonical override prefix for a release site's kernel fields.
    pub fn kernel_prefix(&self, r: &ReleaseSite) -> String {
        format!("kernel:{}:", self.compartment_name(r.compartment))
    }

    /// Every override name this graph accepts, in sorted order.
    ///
    /// Use it to discover the parameter space instead of guessing at names.
    pub fn parameter_names(&self) -> Vec<String> {
        let mut out = Vec::new();
        for t in &self.transfers {
            out.push(self.transit_name(t));
        }
        for a in &self.absorptions {
            out.push(self.absorption_name(a));
        }
        for e in &self.eliminations {
            out.push(self.elimination_name(e));
        }
        for r in &self.releases {
            let prefix = self.kernel_prefix(r);
            for field in r.kernel.fields().keys() {
                out.push(format!("{prefix}{field}"));
            }
        }
        out.sort();
        out
    }

    /// Reject any override name this graph does not accept.
    ///
    /// A misspelt parameter name would otherwise be applied to nothing and the
    /// solve would silently run on defaults.
    pub fn validate_parameters(&self, params: &ModelParameters) -> Result<(), String> {
        let known = self.parameter_names();
        for name in params.names() {
            if !known.iter().any(|k| k == name) {
                return Err(name.to_string());
            }
        }
        Ok(())
    }
}

/// Incremental builder for a [`CompartmentGraph`].
#[derive(Debug, Default)]
pub struct GraphBuilder {
    compartments: Vec<(String, VolumeML, Vec<String>)>,
    species: Vec<String>,
    transfers: Vec<(String, String, String, RateHrInv)>,
    absorptions: Vec<(String, String, String, RateHrInv)>,
    eliminations: Vec<(String, String, RateHrInv, String)>,
    releases: Vec<(String, String, String, Box<dyn ReleaseKernel>)>,
    dose_sites: Vec<String>,
}

impl GraphBuilder {
    /// Declare a species by name.
    pub fn species(mut self, name: &str) -> Self {
        self.species.push(name.to_string());
        self
    }

    /// Declare a compartment, its volume, and the species it admits.
    pub fn compartment(mut self, name: &str, volume: VolumeML, species: &[&str]) -> Self {
        self.compartments.push((
            name.to_string(),
            volume,
            species.iter().map(|s| s.to_string()).collect(),
        ));
        self
    }

    /// Add a first-order transit edge.
    pub fn transfer(mut self, from: &str, to: &str, species: &str, rate: RateHrInv) -> Self {
        self.transfers
            .push((from.to_string(), to.to_string(), species.to_string(), rate));
        self
    }

    /// Add a first-order absorption edge into a systemic compartment.
    pub fn absorption(mut self, from: &str, to: &str, species: &str, rate: RateHrInv) -> Self {
        self.absorptions
            .push((from.to_string(), to.to_string(), species.to_string(), rate));
        self
    }

    /// Add first-order elimination, accounted into a terminal sink compartment.
    pub fn elimination(mut self, from: &str, species: &str, rate: RateHrInv, sink: &str) -> Self {
        self.eliminations.push((
            from.to_string(),
            species.to_string(),
            rate,
            sink.to_string(),
        ));
        self
    }

    /// Attach a release kernel converting `prodrug` into `product` at a site.
    pub fn release(
        mut self,
        compartment: &str,
        prodrug: &str,
        product: &str,
        kernel: Box<dyn ReleaseKernel>,
    ) -> Self {
        self.releases.push((
            compartment.to_string(),
            prodrug.to_string(),
            product.to_string(),
            kernel,
        ));
        self
    }

    /// Mark a compartment as somewhere a dose may be placed.
    pub fn dose_site(mut self, name: &str) -> Self {
        self.dose_sites.push(name.to_string());
        self
    }

    /// Validate and freeze the topology.
    pub fn build(self) -> Result<CompartmentGraph, GraphError> {
        let mut species_names: Vec<String> = Vec::new();
        for s in &self.species {
            if species_names.contains(s) {
                return Err(GraphError::DuplicateSpecies(s.clone()));
            }
            species_names.push(s.clone());
        }
        let sid = |name: &str| -> Result<SpeciesId, GraphError> {
            species_names
                .iter()
                .position(|n| n == name)
                .map(|i| SpeciesId(i as u32))
                .ok_or_else(|| GraphError::UnknownSpecies(name.to_string()))
        };

        let mut names: Vec<String> = Vec::new();
        for (n, v, _) in &self.compartments {
            if names.contains(n) {
                return Err(GraphError::DuplicateCompartment(n.clone()));
            }
            if v.0 <= 0.0 {
                return Err(GraphError::NonPositiveVolume(n.clone(), v.0));
            }
            names.push(n.clone());
        }
        let cid = |name: &str| -> Result<CompartmentId, GraphError> {
            names
                .iter()
                .position(|n| n == name)
                .map(|i| CompartmentId(i as u32))
                .ok_or_else(|| GraphError::UnknownCompartment(name.to_string()))
        };

        let mut compartments = Vec::with_capacity(self.compartments.len());
        for (i, (name, volume, spec)) in self.compartments.iter().enumerate() {
            let mut ids = Vec::with_capacity(spec.len());
            for s in spec {
                ids.push(sid(s)?);
            }
            compartments.push(Compartment {
                id: CompartmentId(i as u32),
                name: name.clone(),
                species: ids,
                volume: *volume,
            });
        }

        let admits =
            |c: CompartmentId, s: SpeciesId| compartments[c.0 as usize].species.contains(&s);

        let mut transfers = Vec::with_capacity(self.transfers.len());
        for (from, to, species, rate) in &self.transfers {
            let (f, t, s) = (cid(from)?, cid(to)?, sid(species)?);
            if rate.0 < 0.0 {
                return Err(GraphError::NegativeRate(
                    format!("k:{from}->{to}:{species}"),
                    rate.0,
                ));
            }
            if !admits(f, s) {
                return Err(GraphError::SpeciesNotInCompartment(f));
            }
            if !admits(t, s) {
                return Err(GraphError::SpeciesNotInCompartment(t));
            }
            transfers.push(Transfer {
                from: f,
                to: t,
                species: s,
                rate: *rate,
            });
        }

        let mut absorptions = Vec::with_capacity(self.absorptions.len());
        for (from, to, species, rate) in &self.absorptions {
            let (f, t, s) = (cid(from)?, cid(to)?, sid(species)?);
            if rate.0 < 0.0 {
                return Err(GraphError::NegativeRate(
                    format!("ka:{from}->{to}:{species}"),
                    rate.0,
                ));
            }
            if !admits(f, s) || !admits(t, s) {
                return Err(GraphError::SpeciesNotInCompartment(f));
            }
            absorptions.push(AbsorptionSink {
                from_compartment: f,
                to_compartment: t,
                species: s,
                rate: *rate,
            });
        }

        let mut eliminations = Vec::with_capacity(self.eliminations.len());
        for (from, species, rate, sink) in &self.eliminations {
            let (f, s, k) = (cid(from)?, sid(species)?, cid(sink)?);
            if rate.0 < 0.0 {
                return Err(GraphError::NegativeRate(
                    format!("ke:{from}:{species}"),
                    rate.0,
                ));
            }
            if !admits(f, s) || !admits(k, s) {
                return Err(GraphError::SpeciesNotInCompartment(f));
            }
            eliminations.push(Elimination {
                compartment: f,
                species: s,
                rate: *rate,
                sink: k,
            });
        }

        let mut releases = Vec::with_capacity(self.releases.len());
        for (compartment, prodrug, product, kernel) in self.releases {
            let (c, p, q) = (cid(&compartment)?, sid(&prodrug)?, sid(&product)?);
            if !admits(c, p) || !admits(c, q) {
                return Err(GraphError::SpeciesNotInCompartment(c));
            }
            releases.push(ReleaseSite {
                compartment: c,
                prodrug: p,
                product: q,
                kernel,
            });
        }

        if self.dose_sites.is_empty() {
            return Err(GraphError::NoDoseSite);
        }
        let mut dose_sites = Vec::with_capacity(self.dose_sites.len());
        for name in &self.dose_sites {
            dose_sites.push(cid(name)?);
        }

        let mut index = BTreeMap::new();
        let mut n = 0usize;
        for c in &compartments {
            for s in &c.species {
                index.insert((c.id, *s), n);
                n += 1;
            }
        }

        Ok(CompartmentGraph {
            compartments,
            species_names,
            transfers,
            releases,
            absorptions,
            eliminations,
            dose_sites,
            index,
            n_state: n,
        })
    }
}
