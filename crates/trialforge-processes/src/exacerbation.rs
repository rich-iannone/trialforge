//! Disease exacerbations (CE) and lung function (RE) layered on a fixed-duration course, for a
//! respiratory study.
//!
//! - **Exacerbations**: a recurrent-event process on the subject's at-risk time from first dose
//!   to the end of participation. The annual rate is the arm's rate × the stratum's multiplier ×
//!   a gamma frailty with mean 1 and variance θ, so counts over a fixed period are negative
//!   binomial (Var = μ + θμ²). Each event lasts a log-normal time; the subject is not at risk
//!   while it lasts nor for 7 days after it resolves. A share of events is severe
//!   (hospitalised), the rest moderate.
//! - **FEV1** (trough, L) at every attended visit: the subject's level + the arm's effect
//!   building up as 1 − e^(−t/τ) (normalised to the full effect at the end of treatment) − a
//!   yearly decline + measurement noise, reported to 0.01 L. The effect is set so that the
//!   mean change at the end of treatment, decline included, is the arm's target.

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject};

use crate::timeline::{ClinicalEvent, Finding, FindingDomain, Severity, Timeline};

/// Days after an exacerbation resolves before another one can start.
pub const EVENT_GAP_DAYS: i32 = 7;

/// Exacerbation and lung-function parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct ExacerbationParams {
    /// Preferred Term of an exacerbation (`CEDECOD`), e.g. the disease.
    pub pt: String,
    /// Reported term (`CETERM`).
    pub term: String,
    /// Annual exacerbation rate per arm (stratum multiplier 1).
    pub rate_by_arm: Vec<f64>,
    /// Rate multiplier per stratum.
    pub stratum_rr: Vec<f64>,
    /// Gamma frailty variance θ (negative binomial overdispersion).
    pub dispersion: f64,
    /// Share of exacerbations that are severe.
    pub severe_fraction: f64,
    /// Median exacerbation duration (days).
    pub median_days: f64,
    /// FEV1 at baseline: population (mean, SD), L.
    pub fev1_baseline: (f64, f64),
    /// Mean change in FEV1 at the end of treatment, per arm (L), including the decline.
    pub fev1_change_by_arm: Vec<f64>,
    /// FEV1 decline (L per year).
    pub fev1_decline: f64,
    /// Onset time constant of the FEV1 effect (days).
    pub onset_days: f64,
    /// Treatment period length (days).
    pub treatment_days: i32,
}

/// Adds exacerbations and FEV1 to `timelines` (in subject order).
pub fn simulate(
    design: &Design,
    p: &ExacerbationParams,
    subjects: &[Subject],
    timelines: &mut [Timeline],
) {
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let stream = design.subject_stream(s.index).label("respiratory");
        let mut d = Draws::new(stream.label("latent"));
        // Draw order is part of the output contract; append new draws at the end.
        let frailty = if p.dispersion > 0.0 {
            d.gamma(1.0 / p.dispersion, p.dispersion)
        } else {
            d.uniform();
            1.0
        };
        let fev1 = d.normal(p.fev1_baseline.0, p.fev1_baseline.1).max(0.5);
        let response = d.normal(0.0, 0.05);

        // Exacerbations.
        let rate = p.rate_by_arm.get(s.arm).copied().unwrap_or(0.0)
            * p.stratum_rr.get(s.stratum).copied().unwrap_or(1.0)
            * frailty;
        let mut e = Draws::new(stream.label("exacerbations"));
        let first = t.first_dose();
        let mut at_risk = first;
        let mut events = Vec::new();
        // No draws at a zero rate.
        if rate > 0.0 {
            loop {
                let gap_years = e.exponential(rate);
                let duration = e.lognormal(p.median_days, 0.5).round().clamp(2.0, 90.0) as i32;
                let severe = e.uniform() < p.severe_fraction;
                // Compare in floating point first: tiny rates give gaps beyond any day count.
                if gap_years * 365.25 > f64::from(t.last_contact.days_since(at_risk)) {
                    break;
                }
                let start = at_risk.add_days((gap_years * 365.25).floor() as i32);
                if start > t.last_contact {
                    break;
                }
                let end = start.add_days(duration - 1);
                events.push(ClinicalEvent {
                    pt: p.pt.clone(),
                    term: p.term.clone(),
                    category: "EXACERBATION",
                    date: start,
                    end: (end <= t.last_contact).then_some(end),
                    severity: Some(if severe {
                        Severity::Severe
                    } else {
                        Severity::Moderate
                    }),
                });
                at_risk = end.add_days(EVENT_GAP_DAYS + 1);
            }
        }
        t.clinical_events = events;

        // FEV1 at every attended visit.
        let full = 1.0 - math::exp(-f64::from(p.treatment_days) / p.onset_days);
        let decline_at_end = p.fev1_decline * f64::from(p.treatment_days) / 365.25;
        let effect =
            p.fev1_change_by_arm.get(s.arm).copied().unwrap_or(0.0) + decline_at_end + response;
        let visits: Vec<(i64, Date)> = t.visits.iter().map(|v| (v.visitnum, v.date)).collect();
        for (visitnum, date) in visits {
            let mut v = Draws::new(stream.label("fev1").index(visitnum as u64));
            let days = f64::from(date.days_since(first)).max(0.0);
            let ramp = ((1.0 - math::exp(-days / p.onset_days)) / full).min(1.0);
            let value = fev1 + effect * ramp - p.fev1_decline * days / 365.25 + v.normal(0.0, 0.1);
            t.findings.push(Finding {
                domain: FindingDomain::Re,
                testcd: "FEV1".into(),
                visitnum,
                date,
                value: math::round_dp(value.max(0.3), 2),
                laterality: None,
                dp: Some(2),
            });
        }
    }
}
