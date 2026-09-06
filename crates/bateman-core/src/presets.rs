//! Ready-made topologies.
//!
//! The rate values here are illustrative defaults chosen to give a
//! well-conditioned test problem. They are not fitted to any published dataset,
//! and nothing in this crate is parameterised to a particular platform. Fit the
//! rates to your own data through [`crate::params::ModelParameters`] before
//! reading anything quantitative into a simulation.

use std::collections::BTreeMap;

use crate::error::GraphError;
use crate::graph::CompartmentGraph;
use crate::kernel::ChargeDependent;
use crate::units::{RateHrInv, VolumeML};

/// Names of the transit compartments, mouth-ward to anus-ward.
pub const GUT_SEGMENTS: [&str; 7] = [
    "stomach",
    "duodenum",
    "jejunum",
    "ileum",
    "cecum",
    "proximal_colon",
    "distal_colon",
];

/// A murine gastrointestinal topology with a plasma compartment.
///
/// Eight compartments (seven gut segments plus `plasma`) and two terminal
/// sinks, `feces_out` and `elim_out`. Three species: a neutral prodrug, an
/// anionic prodrug, and the free product they both release.
///
/// The neutral prodrug releases in the small intestine and the anionic one in
/// the colon, which is the qualitative behaviour a charge-controlled conjugate
/// is designed for. The free product is absorbed into plasma from every gut
/// segment below the stomach and eliminated from plasma.
pub fn murine_gut() -> Result<CompartmentGraph, GraphError> {
    const PRODRUG_NEUTRAL: &str = "prodrug_neutral";
    const PRODRUG_ANIONIC: &str = "prodrug_anionic";
    const PRODUCT: &str = "free_product";

    let all = [PRODRUG_NEUTRAL, PRODRUG_ANIONIC, PRODUCT];
    let volumes: [f64; 7] = [0.5, 0.15, 0.35, 0.35, 0.6, 0.25, 0.25];

    let mut g = CompartmentGraph::builder()
        .species(PRODRUG_NEUTRAL)
        .species(PRODRUG_ANIONIC)
        .species(PRODUCT);

    for (segment, volume) in GUT_SEGMENTS.iter().zip(volumes) {
        g = g.compartment(segment, VolumeML(volume), &all);
    }
    g = g
        .compartment("plasma", VolumeML(2.0), &[PRODUCT])
        .compartment("feces_out", VolumeML(1.0), &all)
        .compartment("elim_out", VolumeML(1.0), &[PRODUCT]);

    // Transit along the tract, every species moving together with the bulk.
    let transit: [f64; 7] = [1.2, 3.0, 1.5, 1.0, 0.35, 0.30, 0.25];
    for (i, rate) in transit.iter().enumerate() {
        let from = GUT_SEGMENTS[i];
        let to = if i + 1 < GUT_SEGMENTS.len() {
            GUT_SEGMENTS[i + 1]
        } else {
            "feces_out"
        };
        for species in all {
            g = g.transfer(from, to, species, RateHrInv(*rate));
        }
    }

    // Absorption of the free product into plasma, below the stomach only.
    for segment in &GUT_SEGMENTS[1..] {
        g = g.absorption(segment, "plasma", PRODUCT, RateHrInv(0.4));
    }

    g = g.elimination("plasma", PRODUCT, RateHrInv(1.8), "elim_out");

    // Charge-controlled release: neutral in the small intestine, anionic in the colon.
    g = g
        .release(
            "ileum",
            PRODRUG_NEUTRAL,
            PRODUCT,
            Box::new(ChargeDependent {
                k_gastric: RateHrInv(0.02),
                k_small: RateHrInv(1.6),
                k_colon: RateHrInv(0.4),
            }),
        )
        .release(
            "cecum",
            PRODRUG_ANIONIC,
            PRODUCT,
            Box::new(ChargeDependent {
                k_gastric: RateHrInv(0.01),
                k_small: RateHrInv(0.05),
                k_colon: RateHrInv(1.4),
            }),
        );

    g.dose_site("stomach").build()
}

/// An even split of a dose between the two prodrug species of [`murine_gut`].
pub fn even_prodrug_mix(graph: &CompartmentGraph) -> BTreeMap<crate::graph::SpeciesId, f64> {
    let mut mix = BTreeMap::new();
    if let Some(s) = graph.species_id("prodrug_neutral") {
        mix.insert(s, 0.5);
    }
    if let Some(s) = graph.species_id("prodrug_anionic") {
        mix.insert(s, 0.5);
    }
    mix
}
