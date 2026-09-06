# bateman

Compartmental pharmacokinetics and pharmacodynamics for orally dosed,
site-of-release-controlled prodrugs.

A prodrug designed to release its payload in a particular part of the gut needs
a model that knows where the release happens. `bateman` builds that model from a
compartment graph: compartments admitting species, first-order transit between
them, release kernels converting a prodrug into its product where the chemistry
says it should, absorption into a systemic compartment, and elimination from it.

The name is the Bateman function, the closed-form solution of the two-step
absorption and elimination cascade that the second integration test checks
against.

## Crates

| Crate | Contents |
|---|---|
| `bateman-core` | compartment graph, release kernels, dosing schedules, adaptive ODE solver, mass-balance diagnostics |
| `bateman-pd` | exposure integrators (AUC, Cmax, Cavg, time above threshold) and response functions (Hill, Emax, linear), composed through `PdLink` |

## Install

```toml
[dependencies]
bateman-core = "0.1"
bateman-pd = "0.1"
```

## Example

```rust
use bateman_core::prelude::*;
use bateman_pd::{Auc, Hill, PdLink};

let graph = presets::murine_gut()?;
let stomach = graph.compartment_id("stomach").unwrap();

let schedule = Schedule::single_gavage(
    DoseSpec::NmolAbsolute(1000.0),
    BodyWeightG(25.0),
    stomach,
    presets::even_prodrug_mix(&graph),
    TimeHr(24.0),
);

let trajectory = Solver::default().run(&graph, &schedule, &ModelParameters::new())?;
assert!(trajectory.diagnostics.relative_error < 1e-6);

let plasma = graph.compartment_id("plasma").unwrap();
let product = graph.species_id("free_product").unwrap();
let index = graph.index_of(plasma, product).unwrap();
let volume = graph.compartment("plasma").unwrap().volume;

let effect = PdLink::new(Auc, Hill { emax: 1.0, ec50: 20.0, n: 2.0 })
    .apply(&trajectory, index, volume);
# Ok::<(), Box<dyn std::error::Error>>(())
```

## Conventions

- Every rate is per hour and every time is in hours. The unit is fixed at the
  type level, so seconds and hours cannot be mixed by accident.
- The state vector is in nanomoles. Concentrations are derived at output time by
  dividing by a compartment volume.
- Release is one for one in molar terms: one nanomole of prodrug becomes one
  nanomole of product.
- Nothing leaves the state vector. Faecal loss and systemic elimination both
  move material into terminal compartments, which is what makes the
  mass-balance report a check on the integration rather than an accounting
  identity. A solve whose relative error exceeds the tolerance fails rather
  than returning a quietly wrong trajectory.

## Fitting parameters

A graph has a default value on every edge, and every edge also has a
canonical name, so a subset of values can be varied without rebuilding the
topology:

| Edge | Name |
|---|---|
| transit | `k:{from}->{to}:{species}` |
| absorption | `ka:{from}->{to}:{species}` |
| elimination | `ke:{compartment}:{species}` |
| release kernel | `kernel:{compartment}:{field}` |

`CompartmentGraph::parameter_names` lists every name a graph accepts. A name the
graph does not accept is rejected before the solve runs, so a misspelt parameter
cannot silently leave the model on its defaults.

## Status

Implemented and tested: the compartment graph and its validation, the three
release kernels, dosing schedules, the deterministic solver with mass-balance
diagnostics, named parameter overrides, and the exposure and response layer.
Correctness is checked against closed forms: single-compartment exponential
decay and the two-step Bateman cascade, both to 1e-6.

Not yet implemented: the stochastic transit-time solver, Bayesian parameter
inference, and the Python bindings.

The rate values in `presets` are illustrative defaults chosen to give a
well-conditioned test problem. They are fitted to no published dataset, and no
platform-specific parameterisation lives in either crate.

## Licence

Apache-2.0. See `LICENSE-APACHE` and `NOTICE`.
