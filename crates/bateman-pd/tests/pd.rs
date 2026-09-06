use bateman_core::units::TimeHr;
use bateman_pd::{
    Auc, Cavg, Cmax, ExposureIntegrator, Hill, Linear, ResponseFn, TimeAboveThreshold,
};

fn grid(n: usize, dt: f64) -> Vec<TimeHr> {
    (0..n).map(|i| TimeHr(i as f64 * dt)).collect()
}

#[test]
fn auc_of_a_constant_is_the_rectangle() {
    let t = grid(11, 0.5);
    let c = vec![4.0; 11];
    assert!((Auc.apply(&t, &c) - 20.0).abs() < 1e-12);
}

#[test]
fn auc_of_a_ramp_is_the_triangle() {
    let t = grid(11, 1.0);
    let c: Vec<f64> = (0..11).map(|i| i as f64).collect();
    // The trapezoid rule is exact on a straight line.
    assert!((Auc.apply(&t, &c) - 50.0).abs() < 1e-12);
}

#[test]
fn cmax_and_cavg_read_the_profile() {
    let t = grid(5, 1.0);
    let c = vec![0.0, 3.0, 9.0, 3.0, 0.0];
    assert!((Cmax.apply(&t, &c) - 9.0).abs() < 1e-12);
    assert!((Cavg.apply(&t, &c) - Auc.apply(&t, &c) / 4.0).abs() < 1e-12);
}

#[test]
fn time_above_threshold_interpolates_the_crossings() {
    // Straight up from 0 to 10 over one hour, straight back down over the next.
    let t = vec![TimeHr(0.0), TimeHr(1.0), TimeHr(2.0)];
    let c = vec![0.0, 10.0, 0.0];
    // At a threshold of 5 the profile is above for the top half of each leg.
    assert!((TimeAboveThreshold(5.0).apply(&t, &c) - 1.0).abs() < 1e-12);
    assert!((TimeAboveThreshold(0.0).apply(&t, &c) - 2.0).abs() < 1e-12);
    assert!(TimeAboveThreshold(20.0).apply(&t, &c).abs() < 1e-12);
}

#[test]
fn hill_is_half_maximal_at_ec50_and_saturates() {
    let h = Hill {
        emax: 2.0,
        ec50: 5.0,
        n: 2.0,
    };
    assert!((h.response(5.0) - 1.0).abs() < 1e-12);
    assert!(h.response(0.0).abs() < 1e-12);
    assert!(h.response(-1.0).abs() < 1e-12);
    assert!(h.response(1e6) > 1.999);
}

#[test]
fn linear_response_is_a_line() {
    let l = Linear {
        slope: 3.0,
        intercept: -1.0,
    };
    assert!((l.response(2.0) - 5.0).abs() < 1e-12);
}
