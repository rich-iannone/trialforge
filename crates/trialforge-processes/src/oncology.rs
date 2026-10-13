//! A simplified solid-tumour course, after `clinical-trials-info/data-raw/sim-onco.R`.
//!
//! Latent times are drawn in continuous time and then observed on the study calendar:
//! progression is seen at the first tumour assessment on or after the latent time, and death,
//! withdrawal and progression are only observed before the data cutoff and before the
//! subject leaves the study. Every subject draws the same uniforms in the same order, whatever
//! their arm, so arms differ only through the parameters that transform those uniforms
//! (common random numbers).

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject};

/// Days per month, as used for medians given in months.
pub const DAYS_PER_MONTH: f64 = 30.4375;

/// Parameters of the course.
#[derive(Clone, Debug, PartialEq)]
pub struct OncologyParams {
    /// Latent median time to progression, in months, per arm.
    pub median_pd_months: Vec<f64>,
    /// Latent median post-progression survival, in months, per arm.
    pub median_pps_months: Vec<f64>,
    /// Hazard multiplier per stratification level (applies to progression and death).
    pub stratum_hr: Vec<f64>,
    /// Monthly hazard of death before progression.
    pub death_before_pd_rate: f64,
    /// Probability of early withdrawal from the study.
    pub withdrawal_prob: f64,
    /// Window (days from randomisation) in which withdrawal happens.
    pub withdrawal_window: (f64, f64),
    /// Withdrawal reasons (`DSDECOD`) and their weights.
    pub withdrawal_reasons: Vec<(String, f64)>,
    /// Probability of stopping treatment for an adverse event.
    pub ae_stop_prob: f64,
    /// Window (days from first dose) in which AE-driven treatment stops happen.
    pub ae_stop_window: (f64, f64),
    /// Days between tumour assessments.
    pub assessment_interval_days: i32,
    /// Maximum number of treatment cycles.
    pub max_cycles: u32,
    /// For subjects alive and in follow-up at cutoff, last contact is up to this many days
    /// before the cutoff.
    pub last_contact_lag_days: i32,
}

/// Why treatment stopped.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TreatmentStop {
    /// Still on treatment at the data cutoff.
    Ongoing,
    /// Progressive disease.
    Progression,
    /// Death.
    Death,
    /// Withdrawal from the study, with its `DSDECOD`.
    Withdrawal(String),
    /// Adverse event.
    AdverseEvent,
}

/// One subject's observed course.
#[derive(Clone, Debug, PartialEq)]
pub struct Course {
    /// Date progression was observed (a tumour assessment date).
    pub progression: Option<Date>,
    /// Date of death, if observed.
    pub death: Option<Date>,
    /// Date of withdrawal from the study, if observed.
    pub withdrawal: Option<Date>,
    /// Reason for withdrawal (`DSDECOD`), when `withdrawal` is set.
    pub withdrawal_reason: Option<String>,
    /// Date treatment stopped (or the cutoff, if ongoing).
    pub treatment_stop: Date,
    /// Why treatment stopped.
    pub treatment_stop_reason: TreatmentStop,
    /// Dose administration dates, in order. Never empty.
    pub doses: Vec<Date>,
    /// Last contact date.
    pub last_contact: Date,
}

/// Simulates the course for every subject.
pub fn simulate(design: &Design, params: &OncologyParams, subjects: &[Subject]) -> Vec<Course> {
    subjects
        .iter()
        .map(|s| simulate_one(design, params, s))
        .collect()
}

fn simulate_one(design: &Design, p: &OncologyParams, s: &Subject) -> Course {
    let mut d = Draws::new(design.subject_stream(s.index).label("oncology"));
    // Draw order is part of the output contract; append new draws at the end.
    let u_pd = d.uniform_pos();
    let u_death = d.uniform_pos();
    let u_pps = d.uniform_pos();
    let withdraws = d.bernoulli(p.withdrawal_prob);
    let u_withdrawal = d.uniform();
    let ae_stop = d.bernoulli(p.ae_stop_prob);
    let u_ae = d.uniform();
    let reason_weights: Vec<f64> = p.withdrawal_reasons.iter().map(|(_, w)| *w).collect();
    let reason = p.withdrawal_reasons[d.weighted(&reason_weights)].0.clone();
    let contact_lag = d.int_range(0, i64::from(p.last_contact_lag_days)) as i32;

    let hr = p.stratum_hr[s.stratum];
    let ln2 = core::f64::consts::LN_2;
    let exp_time = |u: f64, rate: f64| -math::ln(u) / rate;
    let pd_days = exp_time(
        u_pd,
        ln2 / (p.median_pd_months[s.arm] * DAYS_PER_MONTH) * hr,
    );
    let death_pre_pd = exp_time(u_death, p.death_before_pd_rate / DAYS_PER_MONTH * hr);
    let pps_days = exp_time(
        u_pps,
        ln2 / (p.median_pps_months[s.arm] * DAYS_PER_MONTH) * hr,
    );

    let rand = s.randomization;
    let first_dose_day = s.treatment_start.days_since(rand);
    let day = |t: f64| rand.add_days(t.floor() as i32);

    // Nobody dies before their first dose (keeps EX and DM consistent).
    let death_latent = day(death_pre_pd
        .min(pd_days + pps_days)
        .max(f64::from(first_dose_day + 1)));
    let withdrawal_latent = withdraws.then(|| {
        day(p.withdrawal_window.0 + u_withdrawal * (p.withdrawal_window.1 - p.withdrawal_window.0))
    });

    let cutoff = design.cutoff;
    let death = (death_latent <= cutoff && withdrawal_latent.is_none_or(|w| death_latent < w))
        .then_some(death_latent);
    let withdrawal = withdrawal_latent.filter(|&w| w <= cutoff && death.is_none_or(|dd| w < dd));

    // Progression is seen at the first assessment on or after the latent time, if the subject
    // is still in the study then.
    let interval = f64::from(p.assessment_interval_days);
    let pd_assessment = rand.add_days(((pd_days / interval).ceil() * interval) as i32);
    let progression = (pd_assessment <= cutoff
        && pd_assessment < death_latent
        && withdrawal.is_none_or(|w| pd_assessment < w))
    .then_some(pd_assessment);

    let ae_latent = ae_stop.then(|| {
        s.treatment_start.add_days(
            (p.ae_stop_window.0 + u_ae * (p.ae_stop_window.1 - p.ae_stop_window.0)) as i32,
        )
    });

    let candidates = [
        (progression, TreatmentStop::Progression),
        (death, TreatmentStop::Death),
        (withdrawal, TreatmentStop::Withdrawal(reason.clone())),
        (ae_latent, TreatmentStop::AdverseEvent),
    ];
    let (treatment_stop, treatment_stop_reason) = candidates
        .into_iter()
        .filter_map(|(date, why)| date.filter(|&dd| dd <= cutoff).map(|dd| (dd, why)))
        .min_by_key(|(dd, _)| *dd)
        .unwrap_or((cutoff, TreatmentStop::Ongoing));

    let cycle = design.arms[s.arm].treatment.cycle_days;
    let doses: Vec<Date> = (0..p.max_cycles as i32)
        .map(|k| s.treatment_start.add_days(k * cycle))
        .take_while(|&dd| dd == s.treatment_start || dd < treatment_stop)
        .collect();
    let last_dose = *doses.last().expect("the first dose is always given");

    let last_contact = death
        .or(withdrawal)
        .unwrap_or_else(|| cutoff.add_days(-contact_lag))
        .max(last_dose)
        .max(s.treatment_start);

    Course {
        progression,
        death,
        withdrawal,
        withdrawal_reason: withdrawal.map(|_| reason),
        treatment_stop,
        treatment_stop_reason,
        doses,
        last_contact,
    }
}
