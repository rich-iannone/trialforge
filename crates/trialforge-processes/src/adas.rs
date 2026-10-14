//! ADAS-Cog item responses from a latent cognitive-impairment trajectory, layered on a
//! fixed-duration course (FT findings at the assessment visits).
//!
//! Each subject has a latent impairment θ (higher = worse): θ(t) = θ₀ + (β + u)·t + s_arm·g(t),
//! with θ₀ ~ N(μ, 1), a subject slope deviation u, a drug effect s_arm that builds up with
//! g(t) = (1 − e^(−t/τ)) / (1 − e^(−T/τ)) over the treatment period T, and visit-to-visit noise.
//! Items follow a graded response model: P(score ≥ k | θ) = logistic(a·(θ − δ_j − τ_k)), with
//! one uniform per item and visit, so scores are monotone in θ.
//!
//! The controls are in ADAS points and are calibrated to the latent scale deterministically
//! (quadrature over θ₀): μ gives the target mean baseline total, β the placebo change over T,
//! and each s_arm that arm's change.

use std::sync::OnceLock;

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject};

use crate::reference::{AdasItem, adas_items};
use crate::timeline::{Finding, FindingDomain, Timeline};

/// Item discrimination.
const A: f64 = 1.7;
/// Spread of an item's thresholds.
const SPREAD: f64 = 1.0;

/// ADAS-Cog parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct AdasParams {
    /// Mean total at baseline.
    pub baseline_mean: f64,
    /// Mean change in the total over the treatment period, per arm.
    pub change_by_arm: Vec<f64>,
    /// Study days of the assessment visits (scheduled visits on these days, plus early
    /// termination).
    pub visit_days: Vec<i32>,
    /// Treatment period length in days (T).
    pub treatment_days: i32,
    /// Onset time constant of the drug effect, in days (τ).
    pub onset_days: f64,
    /// SD of the subject slope deviation (latent units per year).
    pub slope_sd: f64,
    /// SD of visit-to-visit latent noise.
    pub visit_sd: f64,
    /// Probability an item is not done.
    pub item_missing: f64,
}

/// Calibrated latent parameters: (μ, β per year, shift per arm).
#[derive(Clone, Debug, PartialEq)]
pub struct Calibration {
    /// Mean baseline θ.
    pub mu: f64,
    /// Placebo-arm drift (latent units per year).
    pub beta: f64,
    /// Drug shift at the end of treatment, per arm.
    pub shift: Vec<f64>,
}

/// Item difficulties δ_j: E[score_j | θ = 0] = baseline_fraction · max.
fn deltas() -> &'static [f64] {
    static D: OnceLock<Vec<f64>> = OnceLock::new();
    D.get_or_init(|| {
        adas_items()
            .iter()
            .map(|it| {
                let target = it.baseline_fraction * f64::from(it.max);
                bisect(|d| expected_item(it, d, 0.0), -10.0, 10.0, target, true)
            })
            .collect()
    })
}

fn logistic(x: f64) -> f64 {
    1.0 / (1.0 + math::exp(-x))
}

/// Threshold offsets τ_k (k = 1..=max), spread symmetrically.
fn tau(it: &AdasItem, k: u32) -> f64 {
    SPREAD * math::norm_inv(f64::from(k) / f64::from(it.max + 1))
}

fn p_ge(it: &AdasItem, delta: f64, theta: f64, k: u32) -> f64 {
    logistic(A * (theta - delta - tau(it, k)))
}

fn expected_item(it: &AdasItem, delta: f64, theta: f64) -> f64 {
    (1..=it.max).map(|k| p_ge(it, delta, theta, k)).sum()
}

/// Expected total at latent θ.
pub fn expected_total(theta: f64) -> f64 {
    adas_items()
        .iter()
        .zip(deltas())
        .map(|(it, d)| expected_item(it, *d, theta))
        .sum()
}

/// E over θ₀ ~ N(mu, 1) of f(θ₀), by trapezoidal quadrature on ±8 SD.
fn over_population(mu: f64, f: impl Fn(f64) -> f64) -> f64 {
    let n = 320;
    let (lo, hi) = (-8.0, 8.0);
    let h = (hi - lo) / f64::from(n);
    let mut sum = 0.0;
    for i in 0..=n {
        let z = lo + h * f64::from(i);
        let w = if i == 0 || i == n { 0.5 } else { 1.0 };
        sum += w * math::exp(-z * z / 2.0) * f(mu + z);
    }
    sum * h / math::sqrt(2.0 * core::f64::consts::PI)
}

/// Bisection for an increasing (or, with `increasing = false`, decreasing) function.
fn bisect(f: impl Fn(f64) -> f64, mut lo: f64, mut hi: f64, target: f64, decreasing: bool) -> f64 {
    for _ in 0..100 {
        let mid = (lo + hi) / 2.0;
        let above = f(mid) > target;
        if above != decreasing {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    (lo + hi) / 2.0
}

/// Calibrates the latent model to the targets.
pub fn calibrate(p: &AdasParams) -> Calibration {
    let mu = bisect(
        |m| over_population(m, expected_total),
        -10.0,
        10.0,
        p.baseline_mean,
        false,
    );
    let base = over_population(mu, expected_total);
    let change = |shift: f64| over_population(mu, |t| expected_total(t + shift)) - base;
    let placebo = p.change_by_arm.first().copied().unwrap_or(0.0);
    let years = f64::from(p.treatment_days) / 365.25;
    let drift = bisect(change, -10.0, 10.0, placebo, false);
    let shift = p
        .change_by_arm
        .iter()
        .map(|c| bisect(change, -10.0, 10.0, *c, false) - drift)
        .collect();
    Calibration {
        mu,
        beta: drift / years,
        shift,
    }
}

/// Adds ADAS-Cog item findings to `timelines` (in subject order).
pub fn simulate(design: &Design, p: &AdasParams, subjects: &[Subject], timelines: &mut [Timeline]) {
    let cal = calibrate(p);
    let schedule = design
        .schedule
        .as_ref()
        .expect("an ADAS study has a schedule");
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let stream = design.subject_stream(s.index).label("adas");
        let mut d = Draws::new(stream.label("latent"));
        // Draw order is part of the output contract; append new draws at the end.
        let theta0 = d.normal(cal.mu, 1.0);
        let slope = cal.beta + d.normal(0.0, p.slope_sd);
        let shift = cal.shift.get(s.arm).copied().unwrap_or(0.0);
        let first = t.first_dose();
        let years = |date: Date| f64::from(date.days_since(first)) / 365.25;
        let ramp = |date: Date| {
            let days = f64::from(date.days_since(first)).max(0.0);
            let full = 1.0 - math::exp(-f64::from(p.treatment_days) / p.onset_days);
            ((1.0 - math::exp(-days / p.onset_days)) / full).min(1.0)
        };
        let assessed: Vec<(i64, Date)> = t
            .visits
            .iter()
            .filter(|v| {
                schedule.visit(v.visitnum).is_some_and(|pv| {
                    pv.day.is_some_and(|d| p.visit_days.contains(&d))
                        || pv.kind == trialforge_model::VisitKind::EarlyTermination
                })
            })
            .map(|v| (v.visitnum, v.date))
            .collect();
        for (k, (visitnum, date)) in assessed.into_iter().enumerate() {
            let mut v = Draws::new(stream.label("visit").index(k as u64));
            let theta =
                theta0 + slope * years(date) + shift * ramp(date) + v.normal(0.0, p.visit_sd);
            for (it, delta) in adas_items().iter().zip(deltas()) {
                let u = v.uniform();
                let missing = v.uniform() < p.item_missing;
                let score = (1..=it.max)
                    .filter(|&kk| u < p_ge(it, *delta, theta, kk))
                    .count();
                t.findings.push(Finding {
                    domain: FindingDomain::Ft,
                    testcd: it.testcd.clone(),
                    visitnum,
                    date,
                    value: if missing { f64::NAN } else { score as f64 },
                    laterality: None,
                    dp: None,
                });
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calibration_hits_its_targets() {
        let p = AdasParams {
            baseline_mean: 24.0,
            change_by_arm: vec![2.5, 1.0],
            visit_days: vec![1, 57, 113, 169],
            treatment_days: 168,
            onset_days: 28.0,
            slope_sd: 0.3,
            visit_sd: 0.25,
            item_missing: 0.0,
        };
        let c = calibrate(&p);
        let base = over_population(c.mu, expected_total);
        assert!((base - 24.0).abs() < 1e-6, "{base}");
        let drift = c.beta * 168.0 / 365.25;
        for (arm, want) in [(0usize, 2.5), (1, 1.0)] {
            let got = over_population(c.mu, |t| expected_total(t + drift + c.shift[arm])) - base;
            assert!((got - want).abs() < 1e-6, "{got} vs {want}");
        }
        // Items are harder than the share at θ = 0 for higher θ.
        assert!(expected_total(1.0) > expected_total(0.0));
    }
}
