//! Release kernels: the rate at which a prodrug converts to its product.
//!
//! A kernel answers one question, at one compartment and one time: what
//! first-order rate constant applies to the prodrug sitting there. The solver
//! multiplies that constant by the prodrug amount, so a kernel returning a
//! constant is ordinary first-order release and a kernel reading the clock or
//! the state is anything else.

use std::collections::BTreeMap;
use std::fmt;
use std::sync::Arc;

use crate::params::ModelParameters;
use crate::units::{RateHrInv, TimeHr};

/// What a kernel is allowed to see when it answers.
pub struct KernelCtx<'a> {
    /// Time since the start of the schedule.
    pub t: TimeHr,
    /// The name of the compartment the release is happening in.
    pub compartment: &'a str,
    /// Amount of prodrug currently in that compartment, in nmol.
    pub prodrug_amount: f64,
}

/// A first-order release rate constant, possibly varying in time or state.
pub trait ReleaseKernel: Send + Sync + fmt::Debug {
    /// The rate constant that applies right now.
    fn rate(&self, ctx: &KernelCtx<'_>) -> RateHrInv;

    /// The kernel's own parameters, as `(field, value)` pairs.
    ///
    /// The solver prefixes each field with `kernel:{compartment}:` to build the
    /// override name, so a field named `k_small` at the ileum is overridden by
    /// `kernel:ileum:k_small`.
    fn fields(&self) -> BTreeMap<String, f64>;

    /// A copy of this kernel with `params` applied, given the name prefix that
    /// identifies its site.
    fn with_overrides(&self, prefix: &str, params: &ModelParameters) -> Box<dyn ReleaseKernel>;
}

/// Release whose rate depends on where in the gut the prodrug sits.
///
/// This is the charge-dependent behaviour the ButM platform reports: a
/// negatively charged conjugate that survives the stomach and releases in the
/// colon, against a neutral one that releases in the small intestine.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChargeDependent {
    /// Rate in the stomach and duodenum.
    pub k_gastric: RateHrInv,
    /// Rate through the rest of the small intestine.
    pub k_small: RateHrInv,
    /// Rate from the cecum onward.
    pub k_colon: RateHrInv,
}

impl ChargeDependent {
    fn region(compartment: &str) -> Region {
        match compartment {
            "stomach" | "duodenum" => Region::Gastric,
            "jejunum" | "ileum" => Region::Small,
            _ => Region::Colon,
        }
    }
}

enum Region {
    Gastric,
    Small,
    Colon,
}

impl ReleaseKernel for ChargeDependent {
    fn rate(&self, ctx: &KernelCtx<'_>) -> RateHrInv {
        match Self::region(ctx.compartment) {
            Region::Gastric => self.k_gastric,
            Region::Small => self.k_small,
            Region::Colon => self.k_colon,
        }
    }

    fn fields(&self) -> BTreeMap<String, f64> {
        BTreeMap::from([
            ("k_gastric".to_string(), self.k_gastric.0),
            ("k_small".to_string(), self.k_small.0),
            ("k_colon".to_string(), self.k_colon.0),
        ])
    }

    fn with_overrides(&self, prefix: &str, params: &ModelParameters) -> Box<dyn ReleaseKernel> {
        Box::new(Self {
            k_gastric: RateHrInv(params.get_or(&format!("{prefix}k_gastric"), self.k_gastric.0)),
            k_small: RateHrInv(params.get_or(&format!("{prefix}k_small"), self.k_small.0)),
            k_colon: RateHrInv(params.get_or(&format!("{prefix}k_colon"), self.k_colon.0)),
        })
    }
}

/// Release driven by local enzyme activity, scaled per compartment.
///
/// The rate is `k_base` times the activity recorded for the compartment, and
/// zero where no activity is recorded, so a compartment absent from the map
/// does not release at all.
#[derive(Debug, Clone, PartialEq)]
pub struct EnzymeDriven {
    /// Rate at unit activity.
    pub k_base: RateHrInv,
    /// Relative activity per compartment name.
    pub activity: BTreeMap<String, f64>,
}

impl ReleaseKernel for EnzymeDriven {
    fn rate(&self, ctx: &KernelCtx<'_>) -> RateHrInv {
        let a = self.activity.get(ctx.compartment).copied().unwrap_or(0.0);
        RateHrInv(self.k_base.0 * a)
    }

    fn fields(&self) -> BTreeMap<String, f64> {
        let mut out = BTreeMap::from([("k_base".to_string(), self.k_base.0)]);
        for (site, a) in &self.activity {
            out.insert(format!("activity_{site}"), *a);
        }
        out
    }

    fn with_overrides(&self, prefix: &str, params: &ModelParameters) -> Box<dyn ReleaseKernel> {
        let activity = self
            .activity
            .iter()
            .map(|(site, a)| {
                (
                    site.clone(),
                    params.get_or(&format!("{prefix}activity_{site}"), *a),
                )
            })
            .collect();
        Box::new(Self {
            k_base: RateHrInv(params.get_or(&format!("{prefix}k_base"), self.k_base.0)),
            activity,
        })
    }
}

/// A kernel supplied by the caller.
///
/// The closure sees the same context as any other kernel. It has no named
/// fields, so [`ModelParameters`] cannot reach inside it; a caller who needs a
/// fitted parameter in a custom kernel should implement [`ReleaseKernel`]
/// directly.
#[derive(Clone)]
pub struct UserSupplied(
    /// The rate the caller computes.
    pub Arc<dyn Fn(&KernelCtx<'_>) -> RateHrInv + Send + Sync>,
);

impl fmt::Debug for UserSupplied {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("UserSupplied(<closure>)")
    }
}

impl ReleaseKernel for UserSupplied {
    fn rate(&self, ctx: &KernelCtx<'_>) -> RateHrInv {
        (self.0)(ctx)
    }

    fn fields(&self) -> BTreeMap<String, f64> {
        BTreeMap::new()
    }

    fn with_overrides(&self, _prefix: &str, _params: &ModelParameters) -> Box<dyn ReleaseKernel> {
        Box::new(self.clone())
    }
}
