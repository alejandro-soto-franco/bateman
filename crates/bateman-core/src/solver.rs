//! The deterministic solver: an adaptive Dormand-Prince integration of the
//! graph's mass-action right-hand side, restarted at every dose time.

use ode_solvers::dopri5::Dopri5;
use ode_solvers::{DVector, System};

use crate::dose::Schedule;
use crate::error::SolverError;
use crate::graph::CompartmentGraph;
use crate::kernel::{KernelCtx, ReleaseKernel};
use crate::params::ModelParameters;
use crate::units::{AmountNmol, TimeHr, VolumeML};

type SVec = DVector<f64>;

/// A solved trajectory, plus the mass-balance check taken over it.
#[derive(Debug, Clone)]
pub struct Trajectory {
    /// Output times, in hours, ascending.
    pub t: Vec<TimeHr>,
    /// State at each output time, in nanomoles, indexed as the graph indexes.
    pub y: Vec<Vec<f64>>,
    /// Mass-balance check taken across the whole run.
    pub diagnostics: MassBalanceReport,
}

impl Trajectory {
    /// The amount of one state entry across every output time.
    pub fn series(&self, index: usize) -> Vec<AmountNmol> {
        self.y.iter().map(|row| AmountNmol(row[index])).collect()
    }

    /// The concentration of one state entry across every output time.
    pub fn concentration_series(&self, index: usize, volume: VolumeML) -> Vec<f64> {
        self.y
            .iter()
            .map(|row| AmountNmol(row[index]).concentration(volume).0)
            .collect()
    }
}

/// Total in, total accounted for, and the gap between them.
///
/// Every removal route in this crate moves material into a terminal
/// compartment rather than deleting it, so a correct solve conserves the total
/// exactly and the residual measures integration error alone.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassBalanceReport {
    /// Everything the schedule delivered.
    pub input_total: AmountNmol,
    /// Everything present in the state vector at the end of the run.
    pub final_total: AmountNmol,
    /// Input minus final. Zero for an exact integration.
    pub residual: AmountNmol,
    /// Residual as a fraction of input.
    pub relative_error: f64,
}

/// Integration settings.
#[derive(Debug, Clone, Copy)]
pub struct Solver {
    /// Relative tolerance passed to the integrator.
    pub rtol: f64,
    /// Absolute tolerance passed to the integrator.
    pub atol: f64,
    /// Spacing of the output grid, in hours.
    pub dt_out: f64,
    /// Largest relative mass-balance error accepted before the solve is an error.
    pub mass_tol: f64,
}

impl Default for Solver {
    fn default() -> Self {
        Self {
            rtol: 1e-8,
            atol: 1e-10,
            dt_out: 0.05,
            mass_tol: 1e-6,
        }
    }
}

struct Rhs<'a> {
    transfers: Vec<(usize, usize, f64)>,
    absorptions: Vec<(usize, usize, f64)>,
    eliminations: Vec<(usize, usize, f64)>,
    releases: Vec<(usize, usize, &'a str, Box<dyn ReleaseKernel>)>,
}

impl System<f64, SVec> for Rhs<'_> {
    fn system(&self, t: f64, y: &SVec, dy: &mut SVec) {
        dy.fill(0.0);
        for &(from, to, rate) in &self.transfers {
            let flux = rate * y[from];
            dy[from] -= flux;
            dy[to] += flux;
        }
        for &(from, to, rate) in &self.absorptions {
            let flux = rate * y[from];
            dy[from] -= flux;
            dy[to] += flux;
        }
        for &(from, sink, rate) in &self.eliminations {
            let flux = rate * y[from];
            dy[from] -= flux;
            dy[sink] += flux;
        }
        for (prodrug, product, name, kernel) in &self.releases {
            let amount = y[*prodrug];
            let ctx = KernelCtx {
                t: TimeHr(t),
                compartment: name,
                prodrug_amount: amount,
            };
            let flux = kernel.rate(&ctx).0 * amount;
            dy[*prodrug] -= flux;
            dy[*product] += flux;
        }
    }
}

impl System<f64, SVec> for &Rhs<'_> {
    fn system(&self, t: f64, y: &SVec, dy: &mut SVec) {
        (**self).system(t, y, dy)
    }
}

impl Solver {
    /// Integrate `schedule` through `graph` with `params` applied.
    pub fn run(
        &self,
        graph: &CompartmentGraph,
        schedule: &Schedule,
        params: &ModelParameters,
    ) -> Result<Trajectory, SolverError> {
        if schedule.t_end.0 <= 0.0 {
            return Err(SolverError::NonPositiveHorizon(schedule.t_end.0));
        }
        if let Err(name) = graph.validate_parameters(params) {
            return Err(SolverError::UnknownParameter(name));
        }
        for d in &schedule.doses {
            if d.time.0 < 0.0 || d.time.0 > schedule.t_end.0 {
                return Err(SolverError::DoseOutsideHorizon(d.time.0));
            }
        }

        let n = graph.state_len();
        let rhs = self.build_rhs(graph, params);

        // Dose times, ascending and deduplicated, always including t = 0.
        let mut breakpoints: Vec<f64> = schedule.doses.iter().map(|d| d.time.0).collect();
        breakpoints.push(0.0);
        breakpoints.sort_by(|a, b| a.partial_cmp(b).unwrap());
        breakpoints.dedup_by(|a, b| (*a - *b).abs() < 1e-12);

        let mut y = vec![0.0f64; n];
        let mut out_t: Vec<f64> = Vec::new();
        let mut out_y: Vec<Vec<f64>> = Vec::new();

        for (i, &t0) in breakpoints.iter().enumerate() {
            // Apply every dose landing at this breakpoint before integrating on.
            for d in schedule
                .doses
                .iter()
                .filter(|d| (d.time.0 - t0).abs() < 1e-12)
            {
                for (species, amount) in d.amounts_nmol() {
                    let idx = graph.index_of(d.site, species).ok_or_else(|| {
                        SolverError::UnknownParameter(format!(
                            "dose site {} does not admit species {}",
                            graph.compartment_name(d.site),
                            graph.species_name(species)
                        ))
                    })?;
                    y[idx] += amount;
                }
            }

            let t1 = breakpoints.get(i + 1).copied().unwrap_or(schedule.t_end.0);
            if t1 <= t0 {
                continue;
            }

            let y0 = SVec::from_vec(y.clone());
            let mut stepper = Dopri5::new(&rhs, t0, t1, self.dt_out, y0, self.rtol, self.atol);
            stepper.integrate().map_err(|e| SolverError::Integration {
                t0,
                t1,
                message: e.to_string(),
            })?;

            let xs = stepper.x_out();
            let ys = stepper.y_out();
            for (x, state) in xs.iter().zip(ys.iter()) {
                // Skip the duplicated interval boundary, keeping the post-dose value.
                if out_t.last().is_some_and(|last| (last - x).abs() < 1e-12) {
                    continue;
                }
                out_t.push(*x);
                out_y.push(state.iter().copied().collect());
            }
            if let Some(last) = ys.last() {
                y = last.iter().copied().collect();
            }
        }

        let input_total = schedule.total_nmol();
        let final_total: f64 = y.iter().sum();
        let residual = input_total - final_total;
        let relative_error = if input_total.abs() > 0.0 {
            (residual / input_total).abs()
        } else {
            residual.abs()
        };
        if relative_error > self.mass_tol {
            return Err(SolverError::MassBalance(relative_error, self.mass_tol));
        }

        Ok(Trajectory {
            t: out_t.into_iter().map(TimeHr).collect(),
            y: out_y,
            diagnostics: MassBalanceReport {
                input_total: AmountNmol(input_total),
                final_total: AmountNmol(final_total),
                residual: AmountNmol(residual),
                relative_error,
            },
        })
    }

    fn build_rhs<'a>(&self, graph: &'a CompartmentGraph, params: &ModelParameters) -> Rhs<'a> {
        let mut transfers = Vec::new();
        for t in graph.transfers() {
            let rate = params.get_or(&graph.transit_name(t), t.rate.0);
            let (from, to) = (
                graph.index_of(t.from, t.species).expect("validated"),
                graph.index_of(t.to, t.species).expect("validated"),
            );
            transfers.push((from, to, rate));
        }

        let mut absorptions = Vec::new();
        for a in graph.absorptions() {
            let rate = params.get_or(&graph.absorption_name(a), a.rate.0);
            let (from, to) = (
                graph
                    .index_of(a.from_compartment, a.species)
                    .expect("validated"),
                graph
                    .index_of(a.to_compartment, a.species)
                    .expect("validated"),
            );
            absorptions.push((from, to, rate));
        }

        let mut eliminations = Vec::new();
        for e in graph.eliminations() {
            let rate = params.get_or(&graph.elimination_name(e), e.rate.0);
            let (from, sink) = (
                graph.index_of(e.compartment, e.species).expect("validated"),
                graph.index_of(e.sink, e.species).expect("validated"),
            );
            eliminations.push((from, sink, rate));
        }

        let mut releases = Vec::new();
        for r in graph.releases() {
            let prefix = graph.kernel_prefix(r);
            let kernel = r.kernel.with_overrides(&prefix, params);
            let (prodrug, product) = (
                graph.index_of(r.compartment, r.prodrug).expect("validated"),
                graph.index_of(r.compartment, r.product).expect("validated"),
            );
            releases.push((
                prodrug,
                product,
                graph.compartment_name(r.compartment),
                kernel,
            ));
        }

        Rhs {
            transfers,
            absorptions,
            eliminations,
            releases,
        }
    }
}
