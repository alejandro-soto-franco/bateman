//! Newtypes for the quantities the model carries.
//!
//! Every rate in this crate is per hour and every time is in hours. Mixing
//! seconds and hours is the classic silent error in a compartmental model, so
//! the unit is fixed at the type level rather than left to a comment.

/// An amount of substance, in nanomoles. The state vector stores amounts.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct AmountNmol(
    /// Nanomoles.
    pub f64,
);

/// A concentration, in micromolar. Derived at output time as amount / volume.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct ConcentrationUM(
    /// Micromolar.
    pub f64,
);

/// A compartment volume, in millilitres.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct VolumeML(
    /// Millilitres.
    pub f64,
);

/// A first-order rate constant, per hour.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct RateHrInv(
    /// Per hour.
    pub f64,
);

/// A time, in hours.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct TimeHr(
    /// Hours.
    pub f64,
);

/// A body weight, in grams.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Default)]
pub struct BodyWeightG(
    /// Grams.
    pub f64,
);

impl AmountNmol {
    /// Concentration of this amount distributed through `volume`.
    ///
    /// One nmol in one mL is one micromolar, so the conversion is a division
    /// with no scale factor.
    pub fn concentration(self, volume: VolumeML) -> ConcentrationUM {
        ConcentrationUM(self.0 / volume.0)
    }
}
