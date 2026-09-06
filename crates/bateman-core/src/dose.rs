//! Doses and dosing schedules.

use std::collections::BTreeMap;

use crate::graph::{CompartmentId, SpeciesId};
use crate::units::{BodyWeightG, TimeHr};

/// How much material a dose delivers.
///
/// The mass forms need a molar mass, since the state vector is in nanomoles and
/// nothing else in the model knows the compound's weight.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DoseSpec {
    /// An absolute amount, in nanomoles.
    NmolAbsolute(f64),
    /// A weight-scaled dose, in milligrams per kilogram of body weight.
    MgPerKg {
        /// Dose per kilogram of body weight, in milligrams.
        mg_per_kg: f64,
        /// Molar mass of the dosed compound, in grams per mole.
        molar_mass_g_per_mol: f64,
    },
    /// An absolute mass, in micrograms.
    UgAbsolute {
        /// Absolute mass, in micrograms.
        ug: f64,
        /// Molar mass of the dosed compound, in grams per mole.
        molar_mass_g_per_mol: f64,
    },
}

impl DoseSpec {
    /// The amount delivered, in nanomoles.
    pub fn nmol(&self, body_weight: BodyWeightG) -> f64 {
        match *self {
            DoseSpec::NmolAbsolute(n) => n,
            DoseSpec::MgPerKg {
                mg_per_kg,
                molar_mass_g_per_mol,
            } => {
                let mg = mg_per_kg * (body_weight.0 / 1000.0);
                mg * 1.0e6 / molar_mass_g_per_mol
            }
            DoseSpec::UgAbsolute {
                ug,
                molar_mass_g_per_mol,
            } => ug * 1.0e3 / molar_mass_g_per_mol,
        }
    }
}

/// One administration, split across species by a mixing ratio.
#[derive(Debug, Clone)]
pub struct Dose {
    /// How much is given.
    pub amount: DoseSpec,
    /// Body weight, used by the weight-scaled dose forms.
    pub body_weight: BodyWeightG,
    /// Compartment the dose lands in.
    pub site: CompartmentId,
    /// When it is given, in hours from the start of the schedule.
    pub time: TimeHr,
    /// Fraction of the dose delivered as each species. Fractions are used as
    /// given, so a mix summing to something other than one is a scaled dose.
    pub species_mix: BTreeMap<SpeciesId, f64>,
}

impl Dose {
    /// A single administration of one species.
    pub fn single(
        amount: DoseSpec,
        body_weight: BodyWeightG,
        site: CompartmentId,
        time: TimeHr,
        species: SpeciesId,
    ) -> Self {
        Self {
            amount,
            body_weight,
            site,
            time,
            species_mix: BTreeMap::from([(species, 1.0)]),
        }
    }

    /// Amount of each species this dose delivers, in nanomoles.
    pub fn amounts_nmol(&self) -> BTreeMap<SpeciesId, f64> {
        let total = self.amount.nmol(self.body_weight);
        self.species_mix
            .iter()
            .map(|(s, f)| (*s, total * f))
            .collect()
    }
}

/// A dosing schedule over a finite horizon.
#[derive(Debug, Clone)]
pub struct Schedule {
    /// Every administration, in any order.
    pub doses: Vec<Dose>,
    /// End of the simulation horizon, in hours.
    pub t_end: TimeHr,
}

impl Schedule {
    /// One dose at t = 0, simulated to `t_end`.
    pub fn single_gavage(
        amount: DoseSpec,
        body_weight: BodyWeightG,
        site: CompartmentId,
        species_mix: BTreeMap<SpeciesId, f64>,
        t_end: TimeHr,
    ) -> Self {
        Self {
            doses: vec![Dose {
                amount,
                body_weight,
                site,
                time: TimeHr(0.0),
                species_mix,
            }],
            t_end,
        }
    }

    /// Repeated dosing at a fixed interval, starting at t = 0.
    pub fn repeated(
        amount: DoseSpec,
        body_weight: BodyWeightG,
        site: CompartmentId,
        species_mix: BTreeMap<SpeciesId, f64>,
        interval_hr: f64,
        duration_days: f64,
    ) -> Self {
        let t_end = duration_days * 24.0;
        let mut doses = Vec::new();
        let mut t = 0.0;
        while t < t_end {
            doses.push(Dose {
                amount,
                body_weight,
                site,
                time: TimeHr(t),
                species_mix: species_mix.clone(),
            });
            t += interval_hr;
        }
        Self {
            doses,
            t_end: TimeHr(t_end),
        }
    }

    /// Once daily for `duration_days`.
    pub fn qd(
        amount: DoseSpec,
        body_weight: BodyWeightG,
        site: CompartmentId,
        species_mix: BTreeMap<SpeciesId, f64>,
        duration_days: f64,
    ) -> Self {
        Self::repeated(amount, body_weight, site, species_mix, 24.0, duration_days)
    }

    /// Twice daily for `duration_days`.
    pub fn bid(
        amount: DoseSpec,
        body_weight: BodyWeightG,
        site: CompartmentId,
        species_mix: BTreeMap<SpeciesId, f64>,
        duration_days: f64,
    ) -> Self {
        Self::repeated(amount, body_weight, site, species_mix, 12.0, duration_days)
    }

    /// Total amount delivered across every dose, in nanomoles.
    pub fn total_nmol(&self) -> f64 {
        self.doses
            .iter()
            .map(|d| d.amounts_nmol().values().sum::<f64>())
            .sum()
    }
}
