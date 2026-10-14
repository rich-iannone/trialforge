//! Pharmacokinetics: timed oral dosing and plasma concentrations, layered on a subject's
//! course from another engine (the fixed-duration engine for phase 1 studies).
//!
//! The model is a population model with first-order absorption and one or two compartments:
//!
//! - Apparent clearance, volumes and absorption rate are log-normal around the population
//!   values, with weight scaling (CL^0.75, V^1 at 70 kg) and a dose-proportionality exponent β
//!   (exposure ∝ dose^β, through CL ∝ (dose / reference)^(1 − β)).
//! - Concentrations from repeated doses superpose; observations carry proportional and additive
//!   error and are reported to 3 significant digits, or as BLQ below the LLOQ.
//! - Sampling: a pre-dose sample before the first dose; a full profile after a single dose or
//!   the last of several; a profile to τ after the first of several; and a trough (τ after the
//!   dose, just before the next) after every other dose.
//!
//! Each sample references the most recent dose before it, and its `VISITNUM` is that dose's
//! visit (`DAY k`), so profiles that run into the next days stay in one visit.

use trialforge_core::{Date, DateTime, Draws, math};
use trialforge_model::{Design, Subject};

use crate::timeline::{Exposure, FindingDomain, PkSample, Timeline};

/// Pharmacokinetic parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct PkParams {
    /// Analyte (the treatment name of active arms).
    pub analyte: String,
    /// 1 or 2 compartments.
    pub compartments: u8,
    /// Apparent clearance CL/F (L/h) at 70 kg and the reference dose.
    pub cl: f64,
    /// Apparent central volume V/F (L) at 70 kg.
    pub v: f64,
    /// Absorption rate constant (1/h).
    pub ka: f64,
    /// Apparent intercompartmental clearance Q/F (L/h), two compartments only.
    pub q: f64,
    /// Apparent peripheral volume V2/F (L), two compartments only.
    pub v2: f64,
    /// Between-subject coefficients of variation of CL, V and ka.
    pub cv: (f64, f64, f64),
    /// Dose-proportionality exponent β (1 = linear).
    pub proportionality: f64,
    /// Reference dose (mg) for β.
    pub reference_dose: f64,
    /// Proportional residual error (CV).
    pub proportional_error: f64,
    /// Additive residual error SD (ng/mL).
    pub additive_error: f64,
    /// Lower limit of quantitation (ng/mL).
    pub lloq: f64,
    /// Nominal post-dose sampling hours of a full profile.
    pub profile_hours: Vec<f64>,
    /// Nominal post-dose sampling hours of a first-dose profile in a multiple-dose study (the
    /// last is τ, taken just before the next dose).
    pub tau_hours: Vec<f64>,
    /// Dosing interval τ (h) for repeated doses.
    pub tau: f64,
    /// Planned number of doses.
    pub planned_doses: u32,
    /// Nominal dosing time of day (hour, minute).
    pub dose_time: (u32, u32),
    /// Study days of single administrations in successive crossover periods (period k gets
    /// the k-th day); `None` doses on every day of the exposure the course recorded.
    pub dose_days: Option<Vec<i32>>,
    /// Relative bioavailability by treatment name or crossover element (`ETCD`); default 1.
    pub relative_f: Vec<(String, f64)>,
    /// Absorption-rate multiplier by treatment name or crossover element (`ETCD`); default 1.
    pub relative_ka: Vec<(String, f64)>,
    /// Within-subject (inter-occasion) CV of clearance between administrations.
    pub iov_cv: f64,
}

/// Formats a nominal time point as `PCTPT` (`Pre-dose`, `0.5h Post-dose`, `24h Post-dose`).
pub fn time_point(hours: f64) -> String {
    if hours < 0.0 {
        "Pre-dose".into()
    } else {
        format!("{}h Post-dose", fmt_hours(hours))
    }
}

/// Hours without a trailing `.0` (e.g. `0.5`, `24`).
pub fn fmt_hours(h: f64) -> String {
    let s = format!("{}", math::round_dp(h, 3));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

/// Individual parameters.
#[derive(Clone, Copy)]
struct Individual {
    cl: f64,
    v: f64,
    ka: f64,
    q: f64,
    v2: f64,
}

/// Concentration (ng/mL) at `t` hours after a single oral dose of `dose` mg.
fn single_dose(p: &PkParams, i: &Individual, dose: f64, t: f64) -> f64 {
    if t <= 0.0 || dose <= 0.0 {
        return 0.0;
    }
    let ka = i.ka;
    let mg_per_l = if p.compartments == 1 {
        let k = i.cl / i.v;
        let ka = if (ka - k).abs() < 1e-9 {
            ka * (1.0 + 1e-6)
        } else {
            ka
        };
        dose * ka / (i.v * (ka - k)) * (math::exp(-k * t) - math::exp(-ka * t))
    } else {
        let k10 = i.cl / i.v;
        let k12 = i.q / i.v;
        let k21 = i.q / i.v2;
        let sum = k10 + k12 + k21;
        let disc = math::sqrt(sum * sum - 4.0 * k10 * k21);
        let alpha = (sum + disc) / 2.0;
        let beta = (sum - disc) / 2.0;
        let ka = if (ka - alpha).abs() < 1e-9 || (ka - beta).abs() < 1e-9 {
            ka * (1.0 + 1e-6)
        } else {
            ka
        };
        dose * ka / i.v
            * ((k21 - alpha) / ((ka - alpha) * (beta - alpha)) * math::exp(-alpha * t)
                + (k21 - beta) / ((ka - beta) * (alpha - beta)) * math::exp(-beta * t)
                + (k21 - ka) / ((alpha - ka) * (beta - ka)) * math::exp(-ka * t))
    };
    (mg_per_l * 1000.0).max(0.0)
}

/// Adds timed dosing, PK samples and multi-day visit spans to `timelines` (in subject order).
pub fn simulate(design: &Design, p: &PkParams, subjects: &[Subject], timelines: &mut [Timeline]) {
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        simulate_one(design, p, s, t);
    }
}

fn simulate_one(design: &Design, p: &PkParams, s: &Subject, t: &mut Timeline) {
    let stream = design.subject_stream(s.index).label("pk");
    let mut d = Draws::new(stream.label("individual"));
    // Draw order is part of the output contract; append new draws at the end.
    let omega = |cv: f64| math::sqrt(math::ln(1.0 + cv * cv));
    let eta = [
        d.normal(0.0, omega(p.cv.0)),
        d.normal(0.0, omega(p.cv.1)),
        d.normal(0.0, omega(p.cv.2)),
    ];
    let dose = design.arms[s.arm].treatment.dose;
    let weight = t
        .findings
        .iter()
        .find(|f| f.domain == FindingDomain::Vs && f.testcd == "WEIGHT")
        .map_or(70.0, |f| f.value);
    let wt = weight / 70.0;
    let dose_effect = if dose > 0.0 {
        math::pow(dose / p.reference_dose, 1.0 - p.proportionality)
    } else {
        1.0
    };
    let ind = Individual {
        cl: p.cl * math::pow(wt, 0.75) * dose_effect * math::exp(eta[0]),
        v: p.v * wt * math::exp(eta[1]),
        ka: p.ka * math::exp(eta[2]),
        q: p.q * math::pow(wt, 0.75),
        v2: p.v2 * wt,
    };

    // Administrations: every dosing day of the course, or one per crossover period that
    // starts before treatment ends.
    let first = t.first_dose();
    let arm = &design.arms[s.arm];
    let admin: Vec<(Date, Option<u8>)> = match &p.dose_days {
        Some(days) => days
            .iter()
            .enumerate()
            .map(|(k, day)| (first.add_days(day - 1), Some(k as u8 + 1)))
            .filter(|(date, _)| *date <= t.treatment_end)
            .collect(),
        None => t
            .exposure
            .iter()
            .flat_map(|e| (0..=e.end.days_since(e.start)).map(move |k| (e.start.add_days(k), None)))
            .collect(),
    };
    // A dose is given only at an attended visit (a subject who leaves on a dosing day is not
    // dosed that day).
    let planned: Vec<(i32, i64)> = design
        .schedule
        .as_ref()
        .map(|sch| {
            sch.visits
                .iter()
                .filter_map(|v| Some((v.day?, v.visitnum)))
                .collect()
        })
        .unwrap_or_default();
    let admin: Vec<(Date, Option<u8>)> = admin
        .into_iter()
        .filter(|(date, _)| {
            let day = date.days_since(first) + 1;
            planned
                .iter()
                .filter(|(d, _)| *d == day)
                .all(|(_, num)| t.visits.iter().any(|v| v.visitnum == *num))
        })
        .collect();
    let treatment_of = |period: Option<u8>| match period {
        Some(k) if !arm.periods.is_empty() => {
            let p = &arm.periods[usize::from(k) - 1];
            (&p.treatment, p.etcd.as_str())
        }
        _ => (&arm.treatment, ""),
    };
    let lookup = |table: &[(String, f64)], name: &str| {
        table
            .iter()
            .find(|(n, _)| n == name)
            .map_or(1.0, |(_, x)| *x)
    };
    // (time, study day, amount, individual parameters, relative F)
    let doses: Vec<(DateTime, i64, f64, Individual, f64)> = admin
        .iter()
        .enumerate()
        .map(|(k, &(day, period))| {
            let mut dd = Draws::new(stream.label("dose").index(k as u64));
            let minutes = dd.int_range(0, 15);
            let iov = dd.normal(0.0, omega(p.iov_cv));
            let (trt, etcd) = treatment_of(period);
            let key = |table: &[(String, f64)]| lookup(table, &trt.name) * lookup(table, etcd);
            let visit_day = i64::from(day.days_since(first)) + 1;
            (
                DateTime::new(day, p.dose_time.0, p.dose_time.1).add_minutes(minutes),
                visit_day,
                trt.dose,
                Individual {
                    cl: ind.cl * math::exp(iov),
                    ka: ind.ka * key(&p.relative_ka),
                    ..ind
                },
                key(&p.relative_f),
            )
        })
        .collect();
    t.exposure = admin
        .iter()
        .zip(&doses)
        .map(|(&(day, period), (time, _, amount, _, _))| Exposure {
            start: day,
            end: day,
            dose: *amount,
            time: Some(*time),
            period,
            location: None,
            group: None,
            visitnum: None,
            switched: false,
        })
        .collect();
    let periodic = p.dose_days.is_some();

    let schedule = design
        .schedule
        .as_ref()
        .expect("phase 1 designs have a schedule");
    let visit_of = |day: i64| {
        schedule
            .visits
            .iter()
            .find(|v| v.day == Some(day as i32))
            .map(|v| v.visitnum)
    };
    let last_contact = DateTime::new(t.last_contact, 23, 59);
    let n_planned = p.planned_doses as usize;
    let conc_at = |time: DateTime| -> f64 {
        doses
            .iter()
            .filter(|(dt, ..)| *dt < time)
            .map(|(dt, _, amount, i, f)| f * single_dose(p, i, *amount, time.hours_since(*dt)))
            .sum()
    };
    let mut samples: Vec<PkSample> = Vec::new();
    let mut push = |d: &mut Draws, k: usize, nominal: f64, time: DateTime| {
        let (reference, day) = (doses[k].0, doses[k].1);
        let true_conc = conc_at(time);
        let e1 = d.normal(0.0, 1.0);
        let e2 = d.normal(0.0, 1.0);
        let observed = true_conc * (1.0 + p.proportional_error * e1) + p.additive_error * e2;
        // Reported to 3 significant digits; BLQ when the reported value is below the LLOQ.
        let reported = math::signif(observed, 3);
        let conc = (true_conc > 0.0 && reported >= p.lloq).then_some(reported);
        if time > last_contact {
            return;
        }
        let Some(visitnum) = visit_of(day) else {
            return;
        };
        samples.push(PkSample {
            analyte: p.analyte.clone(),
            visitnum,
            tpt: time_point(nominal),
            tptnum: nominal,
            reference,
            time,
            conc,
            lloq: p.lloq,
        });
    };
    for k in 0..doses.len() {
        let dose_time = doses[k].0;
        let next = doses.get(k + 1).map(|x| x.0);
        // Each dose's samples draw from their own stream.
        let mut d = Draws::new(stream.label("samples").index(k as u64));
        if k == 0 || periodic {
            let lead = d.int_range(5, 60);
            push(&mut d, k, -0.5, dose_time.add_minutes(-lead));
        }
        let hours: &[f64] = if periodic || n_planned <= 1 || k + 1 == n_planned {
            &p.profile_hours
        } else if k == 0 {
            &p.tau_hours
        } else {
            &p.tau_hours[p.tau_hours.len() - 1..]
        };
        for &h in hours {
            let jitter = if h <= 1.0 {
                2
            } else if h <= 12.0 {
                5
            } else {
                15
            };
            let j = d.int_range(-jitter, jitter);
            let early = d.int_range(2, 10);
            let mut time = dose_time.add_minutes((h * 60.0).round() as i64 + j);
            // A sample due at the next dose is taken just before it.
            if let Some(nt) = next.filter(|nt| time >= *nt) {
                time = nt.add_minutes(-early);
            }
            push(&mut d, k, h, time);
        }
    }

    // Visits that hold a profile end on the day of its last sample.
    for v in &mut t.visits {
        let end = samples
            .iter()
            .filter(|x| x.visitnum == v.visitnum)
            .map(|x| x.time.date())
            .max();
        if let Some(e) = end.filter(|e| *e > v.date) {
            v.end = Some(e);
        }
    }
    // Keep only samples whose visit was attended.
    samples.retain(|x| t.visits.iter().any(|v| v.visitnum == x.visitnum));
    t.pk = samples;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn params(compartments: u8) -> PkParams {
        PkParams {
            analyte: "X".into(),
            compartments,
            cl: 10.0,
            v: 50.0,
            ka: 1.0,
            q: 5.0,
            v2: 100.0,
            cv: (0.3, 0.2, 0.4),
            proportionality: 1.0,
            reference_dose: 100.0,
            proportional_error: 0.1,
            additive_error: 0.5,
            lloq: 1.0,
            profile_hours: vec![],
            tau_hours: vec![],
            tau: 24.0,
            planned_doses: 1,
            dose_time: (8, 0),
            dose_days: None,
            relative_f: Vec::new(),
            relative_ka: Vec::new(),
            iov_cv: 0.0,
        }
    }

    fn ind() -> Individual {
        Individual {
            cl: 10.0,
            v: 50.0,
            ka: 1.0,
            q: 5.0,
            v2: 100.0,
        }
    }

    #[test]
    fn one_compartment_auc_is_dose_over_clearance() {
        // Trapezoid AUC over a long fine grid ≈ dose / CL (mg·h/L → ng·h/mL × 1000).
        let p = params(1);
        let mut auc = 0.0;
        let dt = 0.01;
        for k in 0..20_000 {
            let t = k as f64 * dt;
            auc += (single_dose(&p, &ind(), 100.0, t) + single_dose(&p, &ind(), 100.0, t + dt))
                / 2.0
                * dt;
        }
        assert!(
            (auc - 100.0 / 10.0 * 1000.0).abs() / 10_000.0 < 1e-3,
            "{auc}"
        );
    }

    #[test]
    fn two_compartment_auc_is_dose_over_clearance() {
        let p = params(2);
        let mut auc = 0.0;
        let dt = 0.02;
        for k in 0..50_000 {
            let t = k as f64 * dt;
            auc += (single_dose(&p, &ind(), 100.0, t) + single_dose(&p, &ind(), 100.0, t + dt))
                / 2.0
                * dt;
        }
        assert!((auc - 10_000.0).abs() / 10_000.0 < 2e-3, "{auc}");
    }

    #[test]
    fn time_points_parse_as_admiral_expects() {
        assert_eq!(time_point(-0.5), "Pre-dose");
        assert_eq!(time_point(0.5), "0.5h Post-dose");
        assert_eq!(time_point(24.0), "24h Post-dose");
    }
}
