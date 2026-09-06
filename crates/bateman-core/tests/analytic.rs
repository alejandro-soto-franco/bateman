//! Correctness checks against closed-form solutions and conservation laws.

use bateman_core::prelude::*;
use std::collections::BTreeMap;

/// A single compartment eliminating into a sink at rate `k`.
fn one_compartment(k: f64) -> CompartmentGraph {
    CompartmentGraph::builder()
        .species("drug")
        .compartment("central", VolumeML(1.0), &["drug"])
        .compartment("out", VolumeML(1.0), &["drug"])
        .elimination("central", "drug", RateHrInv(k), "out")
        .dose_site("central")
        .build()
        .expect("valid graph")
}

#[test]
fn first_order_elimination_matches_the_exponential() {
    let k = 0.7;
    let graph = one_compartment(k);
    let central = graph.compartment_id("central").unwrap();
    let drug = graph.species_id("drug").unwrap();
    let idx = graph.index_of(central, drug).unwrap();

    let dose = 1000.0;
    let schedule = Schedule::single_gavage(
        DoseSpec::NmolAbsolute(dose),
        BodyWeightG(25.0),
        central,
        BTreeMap::from([(drug, 1.0)]),
        TimeHr(10.0),
    );

    let traj = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .expect("solve");

    for (t, row) in traj.t.iter().zip(&traj.y) {
        let want = dose * (-k * t.0).exp();
        assert!(
            (row[idx] - want).abs() <= 1e-6 * dose.max(1.0),
            "at t={}: got {}, want {}",
            t.0,
            row[idx],
            want
        );
    }
}

#[test]
fn two_step_transfer_matches_the_bateman_function() {
    // A -> B -> sink with distinct rates is the classic Bateman solution, which
    // is where this crate takes its name.
    let (ka, ke) = (1.3, 0.35);
    let graph = CompartmentGraph::builder()
        .species("drug")
        .compartment("gut", VolumeML(1.0), &["drug"])
        .compartment("plasma", VolumeML(1.0), &["drug"])
        .compartment("out", VolumeML(1.0), &["drug"])
        .absorption("gut", "plasma", "drug", RateHrInv(ka))
        .elimination("plasma", "drug", RateHrInv(ke), "out")
        .dose_site("gut")
        .build()
        .expect("valid graph");

    let gut = graph.compartment_id("gut").unwrap();
    let plasma = graph.compartment_id("plasma").unwrap();
    let drug = graph.species_id("drug").unwrap();
    let idx = graph.index_of(plasma, drug).unwrap();

    let dose = 500.0;
    let schedule = Schedule::single_gavage(
        DoseSpec::NmolAbsolute(dose),
        BodyWeightG(25.0),
        gut,
        BTreeMap::from([(drug, 1.0)]),
        TimeHr(12.0),
    );
    let traj = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .expect("solve");

    for (t, row) in traj.t.iter().zip(&traj.y) {
        let want = dose * ka / (ka - ke) * ((-ke * t.0).exp() - (-ka * t.0).exp());
        assert!(
            (row[idx] - want).abs() <= 1e-6 * dose,
            "at t={}: got {}, want {}",
            t.0,
            row[idx],
            want
        );
    }
}

#[test]
fn release_conserves_moles_one_for_one() {
    let graph = CompartmentGraph::builder()
        .species("prodrug")
        .species("product")
        .compartment("site", VolumeML(1.0), &["prodrug", "product"])
        .release(
            "site",
            "prodrug",
            "product",
            Box::new(ChargeDependent {
                k_gastric: RateHrInv(0.5),
                k_small: RateHrInv(0.5),
                k_colon: RateHrInv(0.5),
            }),
        )
        .dose_site("site")
        .build()
        .expect("valid graph");

    let site = graph.compartment_id("site").unwrap();
    let prodrug = graph.species_id("prodrug").unwrap();
    let product = graph.species_id("product").unwrap();

    let schedule = Schedule::single_gavage(
        DoseSpec::NmolAbsolute(100.0),
        BodyWeightG(25.0),
        site,
        BTreeMap::from([(prodrug, 1.0)]),
        TimeHr(8.0),
    );
    let traj = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .expect("solve");

    let (i, j) = (
        graph.index_of(site, prodrug).unwrap(),
        graph.index_of(site, product).unwrap(),
    );
    for row in &traj.y {
        assert!((row[i] + row[j] - 100.0).abs() < 1e-8);
    }
    // The colon rate is the one that applies at an unnamed site, and it is 0.5.
    let last = traj.y.last().unwrap();
    assert!((last[i] - 100.0 * (-0.5f64 * 8.0).exp()).abs() < 1e-6);
}

#[test]
fn the_preset_conserves_mass_and_reaches_plasma() {
    let graph = presets::murine_gut().expect("preset builds");
    let stomach = graph.compartment_id("stomach").unwrap();
    let plasma = graph.compartment_id("plasma").unwrap();
    let product = graph.species_id("free_product").unwrap();
    let idx = graph.index_of(plasma, product).unwrap();

    let schedule = Schedule::single_gavage(
        DoseSpec::NmolAbsolute(1000.0),
        BodyWeightG(25.0),
        stomach,
        presets::even_prodrug_mix(&graph),
        TimeHr(24.0),
    );
    let traj = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .expect("solve");

    assert!(traj.diagnostics.relative_error < 1e-6);
    let peak = traj.y.iter().map(|r| r[idx]).fold(0.0, f64::max);
    assert!(peak > 0.0, "the product must reach plasma");
}

#[test]
fn repeated_dosing_accumulates_every_dose() {
    let graph = one_compartment(0.1);
    let central = graph.compartment_id("central").unwrap();
    let drug = graph.species_id("drug").unwrap();

    let schedule = Schedule::qd(
        DoseSpec::NmolAbsolute(100.0),
        BodyWeightG(25.0),
        central,
        BTreeMap::from([(drug, 1.0)]),
        3.0,
    );
    assert_eq!(schedule.doses.len(), 3);
    assert!((schedule.total_nmol() - 300.0).abs() < 1e-12);

    let traj = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .expect("solve");
    assert!(traj.diagnostics.relative_error < 1e-6);
}

#[test]
fn an_override_changes_the_answer_and_a_typo_is_rejected() {
    let graph = one_compartment(0.7);
    let central = graph.compartment_id("central").unwrap();
    let drug = graph.species_id("drug").unwrap();
    let idx = graph.index_of(central, drug).unwrap();
    let schedule = Schedule::single_gavage(
        DoseSpec::NmolAbsolute(1000.0),
        BodyWeightG(25.0),
        central,
        BTreeMap::from([(drug, 1.0)]),
        TimeHr(5.0),
    );

    assert!(
        graph
            .parameter_names()
            .contains(&"ke:central:drug".to_string())
    );

    let base = Solver::default()
        .run(&graph, &schedule, &ModelParameters::new())
        .unwrap();
    let faster = Solver::default()
        .run(
            &graph,
            &schedule,
            &ModelParameters::new().with("ke:central:drug", 2.0),
        )
        .unwrap();
    assert!(faster.y.last().unwrap()[idx] < base.y.last().unwrap()[idx]);

    let typo = Solver::default().run(
        &graph,
        &schedule,
        &ModelParameters::new().with("ke:centrl:drug", 2.0),
    );
    assert!(matches!(typo, Err(SolverError::UnknownParameter(_))));
}

#[test]
fn an_invalid_graph_is_rejected_at_build_time() {
    let no_dose = CompartmentGraph::builder()
        .species("drug")
        .compartment("a", VolumeML(1.0), &["drug"])
        .build();
    assert_eq!(no_dose.unwrap_err(), GraphError::NoDoseSite);

    let bad_volume = CompartmentGraph::builder()
        .species("drug")
        .compartment("a", VolumeML(0.0), &["drug"])
        .dose_site("a")
        .build();
    assert!(matches!(
        bad_volume.unwrap_err(),
        GraphError::NonPositiveVolume(_, _)
    ));

    let unknown_species = CompartmentGraph::builder()
        .species("drug")
        .compartment("a", VolumeML(1.0), &["ghost"])
        .dose_site("a")
        .build();
    assert!(matches!(
        unknown_species.unwrap_err(),
        GraphError::UnknownSpecies(_)
    ));

    let negative = CompartmentGraph::builder()
        .species("drug")
        .compartment("a", VolumeML(1.0), &["drug"])
        .compartment("b", VolumeML(1.0), &["drug"])
        .transfer("a", "b", "drug", RateHrInv(-1.0))
        .dose_site("a")
        .build();
    assert!(matches!(
        negative.unwrap_err(),
        GraphError::NegativeRate(_, _)
    ));
}

#[test]
fn dose_specs_convert_to_nanomoles() {
    // 10 mg/kg into a 25 g mouse is 0.25 mg; at 100 g/mol that is 2500 nmol.
    let d = DoseSpec::MgPerKg {
        mg_per_kg: 10.0,
        molar_mass_g_per_mol: 100.0,
    };
    assert!((d.nmol(BodyWeightG(25.0)) - 2500.0).abs() < 1e-9);

    // 50 ug at 250 g/mol is 200 nmol.
    let d = DoseSpec::UgAbsolute {
        ug: 50.0,
        molar_mass_g_per_mol: 250.0,
    };
    assert!((d.nmol(BodyWeightG(25.0)) - 200.0).abs() < 1e-9);
}
