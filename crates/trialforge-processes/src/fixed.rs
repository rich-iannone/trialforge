//! The course of a fixed-duration study with a visit schedule (e.g. a parallel-group trial
//! with a continuous endpoint).
//!
//! Per subject, in this order:
//!
//! 1. **Adverse events** over the whole potential observation window, each Preferred Term on
//!    its own stream: a Poisson process with the background rate × a subject frailty, times the
//!    arm's relative risk while on treatment.
//! 2. **Disposition**: the earliest of AE-driven discontinuation, dropout and death ends
//!    treatment early; otherwise the subject completes treatment and follow-up.
//! 3. **Visits** on schedule within their windows (some missed), plus an early-termination
//!    visit for subjects who leave early (unless lost to follow-up or dead).
//! 4. **Dosing** daily from Day 1 to the last dose, interrupted by some AEs and split into one
//!    exposure interval per dispensing visit.
//! 5. **Medical history**, **medications** (chronic, for ongoing history, and acute, for
//!    AEs) and **findings** (labs and vital signs at attended visits).
//!
//! Every subject draws from its own streams in a fixed order, so arms differ only through the
//! parameters applied to those draws.

use std::collections::HashMap;

use trialforge_core::{Date, Draws, Stream, math};
use trialforge_model::{Design, Epoch, PlannedVisit, Schedule, Subject, VisitKind};

use crate::reference::{self, Scale, TreatmentKind};
use crate::timeline::*;

const DAYS_PER_YEAR: f64 = 365.25;

/// The expected change in an analyte by the end of treatment, per arm.
#[derive(Clone, Debug, PartialEq)]
pub struct EndpointEffect {
    /// `LBTESTCD` or `VSTESTCD` (must be a linear-scale analyte).
    pub testcd: String,
    /// Mean change from baseline at the end of treatment, per arm (standard units).
    pub change_by_arm: Vec<f64>,
    /// Time constant (days) of the approach to the full effect.
    pub onset_days: f64,
}

/// Overrides for an analyte's population distribution.
#[derive(Clone, Debug, PartialEq)]
pub struct AnalyteOverride {
    /// `LBTESTCD` or `VSTESTCD`.
    pub testcd: String,
    /// Population mean for (female, male), standard units.
    pub mean: Option<(f64, f64)>,
    /// Between-subject SD on the model scale.
    pub sd_between: Option<f64>,
    /// Inclusion range for the screening value (subjects outside are not enrolled).
    pub inclusion: Option<(f64, f64)>,
}

/// Parameters of a fixed-duration course.
#[derive(Clone, Debug, PartialEq)]
pub struct CourseParams {
    /// Annual hazard of dropping out (other than for AEs or death).
    pub dropout_annual_rate: f64,
    /// Dropout reasons (`DSDECOD`) and weights.
    pub dropout_reasons: Vec<(String, f64)>,
    /// Annual hazard of death.
    pub death_annual_rate: f64,
    /// Probability of missing an optional visit.
    pub missed_visit_prob: f64,
    /// Multiplier on every background AE rate.
    pub ae_rate_multiplier: f64,
    /// Per arm: Preferred Term → relative risk while on treatment.
    pub ae_relative_risk: Vec<HashMap<String, f64>>,
    /// Probability an on-treatment AE leads to discontinuation, by severity.
    pub ae_discontinuation: [f64; 3],
    /// Probability an on-treatment AE interrupts dosing, by severity.
    pub ae_interruption: [f64; 3],
    /// Medical-history conditions every subject has (e.g. the disease under study).
    pub required_history: Vec<String>,
    /// (indication, ingredient, probability): background therapy for a required condition.
    pub background_medication: Vec<(String, String, f64)>,
    /// Analyte distribution overrides.
    pub analyte_overrides: Vec<AnalyteOverride>,
    /// Treatment effects on analytes.
    pub endpoints: Vec<EndpointEffect>,
    /// Baseline BMI (mean, SD).
    pub bmi: (f64, f64),
    /// Dose only on the first this-many days of treatment (e.g. 1 for a single-dose study);
    /// `None` doses every day of the treatment period.
    pub dosing_days: Option<i32>,
}

/// Why a course could not be simulated.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CourseError(pub String);

impl core::fmt::Display for CourseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for CourseError {}

/// Simulates every subject's course.
pub fn simulate(
    design: &Design,
    params: &CourseParams,
    subjects: &[Subject],
) -> Result<Vec<Timeline>, CourseError> {
    let schedule = design
        .schedule
        .as_ref()
        .ok_or_else(|| CourseError("a fixed-duration course needs a visit schedule".into()))?;
    validate(design, params)?;
    Ok(subjects
        .iter()
        .map(|s| simulate_one(design, schedule, params, s))
        .collect())
}

fn validate(design: &Design, p: &CourseParams) -> Result<(), CourseError> {
    let err = |m: String| Err(CourseError(m));
    if p.ae_relative_risk.len() != design.arms.len() {
        return err("`ae_relative_risk` needs one entry per arm".into());
    }
    for e in &p.endpoints {
        if e.change_by_arm.len() != design.arms.len() {
            return err(format!("endpoint {} needs one change per arm", e.testcd));
        }
        let linear = reference::lab_panel()
            .iter()
            .any(|a| a.testcd == e.testcd && a.scale == Scale::Linear)
            || reference::vital_signs()
                .iter()
                .any(|v| v.testcd == e.testcd);
        if !linear {
            return err(format!(
                "endpoint {} must be a linear-scale lab analyte or a vital sign",
                e.testcd
            ));
        }
    }
    let schedule = design.schedule.as_ref().expect("checked by the caller");
    for kind in [VisitKind::Baseline, VisitKind::EndOfTreatment] {
        if schedule.visits.iter().filter(|v| v.kind == kind).count() != 1 {
            return err(format!(
                "a fixed-duration schedule needs exactly one {kind:?} visit"
            ));
        }
    }
    for o in &p.analyte_overrides {
        let log_lab = reference::lab_panel()
            .iter()
            .any(|a| a.testcd == o.testcd && a.scale == Scale::Log);
        if o.inclusion.is_some() && log_lab {
            return err(format!(
                "inclusion ranges are only supported on linear-scale analytes ({})",
                o.testcd
            ));
        }
    }
    for pt in &p.required_history {
        if !reference::history_conditions().iter().any(|h| &h.pt == pt) {
            return err(format!(
                "required history {pt:?} is not in the medical-history reference"
            ));
        }
    }
    for (rate, name) in [
        (p.dropout_annual_rate, "dropout"),
        (p.death_annual_rate, "death"),
    ] {
        if !(rate >= 0.0 && rate.is_finite()) {
            return err(format!("the {name} rate must be a non-negative number"));
        }
    }
    Ok(())
}

/// A candidate AE before disposition is known.
struct AeDraft {
    pt: String,
    start: Date,
    duration: i32,
    severity: Severity,
    serious: bool,
    u_disc: f64,
    u_interrupt: f64,
    u_related: [f64; 2],
    u_outcome: f64,
    u_treat: f64,
    u_med: f64,
}

/// Why treatment ended early.
#[derive(Clone, Debug)]
enum Reason {
    AdverseEvent(usize),
    Dropout(String),
    Death,
}

fn severity_index(s: Severity) -> usize {
    match s {
        Severity::Mild => 0,
        Severity::Moderate => 1,
        Severity::Severe => 2,
    }
}

fn simulate_one(design: &Design, schedule: &Schedule, p: &CourseParams, s: &Subject) -> Timeline {
    let stream = design.subject_stream(s.index);
    let first_dose = s.treatment_start;
    let day = |d: i32| first_dose.add_days(d - 1);
    let eot_visit = schedule
        .visits
        .iter()
        .find(|v| v.kind == VisitKind::EndOfTreatment)
        .expect("validated");
    let fu_visit = schedule
        .visits
        .iter()
        .find(|v| v.kind == VisitKind::FollowUp);

    // Visit dates as planned (with window jitter); whether each is attended is decided later.
    let mut vd = Draws::new(stream.label("visits"));
    let mut planned_dates: Vec<(usize, Date, bool)> = Vec::new(); // (index, date, missed)
    let mut last = s.screening;
    for (i, v) in schedule.visits.iter().enumerate() {
        let jitter = vd.int_range(-i64::from(v.window.0), i64::from(v.window.1)) as i32;
        let missed_draw = vd.uniform();
        let date = match v.kind {
            VisitKind::Screening => s.screening,
            VisitKind::Baseline => first_dose,
            VisitKind::EarlyTermination => continue,
            VisitKind::FollowUp => continue, // placed after the actual end-of-treatment date
            _ => day(v.day.expect("scheduled visits have a day")).add_days(jitter),
        };
        let date = if v.kind == VisitKind::Screening {
            date
        } else {
            date.max(last.add_days(1))
        };
        last = date;
        let optional = matches!(v.kind, VisitKind::Treatment);
        planned_dates.push((i, date, optional && missed_draw < p.missed_visit_prob));
    }
    let eot_date = planned_dates
        .iter()
        .find(|(i, _, _)| schedule.visits[*i].kind == VisitKind::EndOfTreatment)
        .map(|(_, d, _)| *d)
        .expect("validated");
    let fu_date = fu_visit.map(|fu| {
        let gap = fu.day.unwrap_or(eot_visit.day.unwrap_or(0)) - eot_visit.day.unwrap_or(0);
        let jitter = vd.int_range(-i64::from(fu.window.0), i64::from(fu.window.1)) as i32;
        eot_date.add_days((gap + jitter).max(1))
    });
    let window_end = fu_date.unwrap_or(eot_date);

    // 1. Adverse events.
    let rr = &p.ae_relative_risk[s.arm];
    let frailty = Draws::new(stream.label("ae-frailty")).gamma(2.0, 0.5);
    let mut drafts: Vec<AeDraft> = Vec::new();
    for bg in reference::ae_background() {
        let mut d = Draws::new(stream.label("ae").label(&bg.pt));
        let base = bg.rate_per_year * p.ae_rate_multiplier * frailty / DAYS_PER_YEAR;
        let on = base * rr.get(&bg.pt).copied().unwrap_or(1.0);
        let mut t = s.consent;
        loop {
            // Piecewise-constant hazard: background before first dose and after planned EOT.
            let rate = if t >= first_dose && t <= eot_date {
                on
            } else {
                base
            };
            if rate <= 0.0 {
                break;
            }
            let gap = d.exponential(rate).ceil().max(1.0);
            // Compare in floating point first: tiny rates give gaps far beyond any i32 day count.
            if gap > f64::from(window_end.days_since(t)) {
                break;
            }
            let next = t.add_days(gap as i32);
            // The hazard is re-evaluated at each event, which approximates the piecewise rate
            // closely enough for daily data.
            t = next;
            let sev_i = d.weighted(&bg.severity);
            let severity = [Severity::Mild, Severity::Moderate, Severity::Severe][sev_i];
            let serious = d.bernoulli((bg.p_serious * [1.0, 1.5, 3.0][sev_i]).min(1.0));
            let duration = d
                .lognormal(bg.median_days * [1.0, 1.5, 2.5][sev_i], 0.8)
                .round()
                .max(1.0) as i32;
            drafts.push(AeDraft {
                pt: bg.pt.clone(),
                start: t,
                duration,
                severity,
                serious,
                u_disc: d.uniform(),
                u_interrupt: d.uniform(),
                u_related: [d.uniform(), d.uniform()],
                u_outcome: d.uniform(),
                u_treat: d.uniform(),
                u_med: d.uniform(),
            });
        }
    }
    drafts.sort_by(|a, b| (a.start, &a.pt).cmp(&(b.start, &b.pt)));

    // 2. Disposition.
    let mut cd = Draws::new(stream.label("course"));
    let dropout_days = cd.exponential(p.dropout_annual_rate.max(1e-12) / DAYS_PER_YEAR);
    let reason_w: Vec<f64> = p.dropout_reasons.iter().map(|(_, w)| *w).collect();
    let dropout_reason = p.dropout_reasons[cd.weighted(&reason_w)].0.clone();
    let death_days = cd.exponential(p.death_annual_rate.max(1e-12) / DAYS_PER_YEAR);
    let fatal_pick = cd.uniform();
    let fatal_lead = cd.int_range(0, 5) as i32;

    let ae_stop = drafts.iter().enumerate().find_map(|(i, a)| {
        let on_treatment = a.start >= first_dose && a.start < eot_date;
        (on_treatment && a.u_disc < p.ae_discontinuation[severity_index(a.severity)]).then(|| {
            (
                a.start.add_days(2.min(a.duration - 1).max(0)),
                Reason::AdverseEvent(i),
            )
        })
    });
    let dropout = (
        first_dose.add_days(dropout_days.floor().min(1e6) as i32),
        Reason::Dropout(dropout_reason),
    );
    let death_latent = first_dose.add_days(death_days.floor().min(1e6) as i32);
    let candidates = [ae_stop, Some(dropout), Some((death_latent, Reason::Death))];
    let early = candidates
        .into_iter()
        .flatten()
        .filter(|(d, _)| *d < eot_date)
        .min_by_key(|(d, _)| *d);

    let (treatment_end, last_dose, last_contact, eot_reason, eos_reason, death) = match &early {
        Some((d, Reason::Death)) => (
            *d,
            d.add_days(-1).max(first_dose),
            *d,
            "DEATH".to_string(),
            "DEATH".to_string(),
            Some(*d),
        ),
        Some((d, r)) => {
            let reason = match r {
                Reason::AdverseEvent(_) => "ADVERSE EVENT".to_string(),
                Reason::Dropout(why) => why.clone(),
                Reason::Death => unreachable!(),
            };
            (
                *d,
                d.add_days(-1).max(first_dose),
                *d,
                reason.clone(),
                reason,
                None,
            )
        }
        None => {
            let end = fu_date.unwrap_or(eot_date);
            if death_latent > eot_date && death_latent < end {
                (
                    eot_date,
                    eot_date.add_days(-1),
                    death_latent,
                    "COMPLETED".into(),
                    "DEATH".into(),
                    Some(death_latent),
                )
            } else {
                (
                    eot_date,
                    eot_date.add_days(-1),
                    end,
                    "COMPLETED".into(),
                    "COMPLETED".into(),
                    None,
                )
            }
        }
    };
    let last_dose = last_dose.max(first_dose);
    let last_dose = p.dosing_days.map_or(last_dose, |n| {
        last_dose.min(first_dose.add_days(n.max(1) - 1))
    });
    let discontinuing_ae = match &early {
        Some((_, Reason::AdverseEvent(i))) => Some(*i),
        _ => None,
    };

    // AEs as observed: start within participation; ends after the last contact are ongoing.
    let mut interruptions: Vec<(Date, Date)> = Vec::new();
    let mut adverse_events = Vec::new();
    let mut acute_meds = Vec::new();
    let active = rr.values().any(|&r| r > 1.0);
    for (i, a) in drafts.iter().enumerate() {
        if a.start > last_contact {
            continue;
        }
        let end = a.start.add_days(a.duration - 1);
        let end = (end <= last_contact).then_some(end);
        let on_drug = a.start >= first_dose && a.start <= last_dose;
        let sev_i = severity_index(a.severity);
        let action = if !on_drug {
            "NOT APPLICABLE"
        } else if Some(i) == discontinuing_ae {
            "DRUG WITHDRAWN"
        } else if a.u_interrupt < p.ae_interruption[sev_i]
            && a.start > first_dose
            && a.start < last_dose
        {
            let stop = end
                .unwrap_or(last_contact)
                .min(a.start.add_days(6))
                .min(last_dose.add_days(-1));
            if stop >= a.start {
                interruptions.push((a.start, stop));
            }
            "DRUG INTERRUPTED"
        } else {
            "DOSE NOT CHANGED"
        };
        let term_rr = rr.get(&a.pt).copied().unwrap_or(1.0);
        let p_related = if a.start < first_dose {
            0.0
        } else {
            0.10 + if active && term_rr > 1.0 {
                0.5 * (1.0 - 1.0 / term_rr)
            } else {
                0.0
            }
        };
        let causality = if a.u_related[0] < p_related {
            ["POSSIBLY RELATED", "PROBABLY RELATED", "RELATED"][(a.u_related[1] * 3.0) as usize % 3]
        } else if a.u_related[1] < 0.3 && a.start >= first_dose {
            "UNLIKELY RELATED"
        } else {
            "NOT RELATED"
        };
        let outcome = match end {
            Some(_) if a.severity == Severity::Severe && a.u_outcome < 0.1 => {
                "RECOVERED/RESOLVED WITH SEQUELAE"
            }
            Some(_) => "RECOVERED/RESOLVED",
            None if a.u_outcome < 0.7 => "NOT RECOVERED/NOT RESOLVED",
            None => "RECOVERING/RESOLVING",
        };
        let serious = a.serious.then_some(SeriousCriteria {
            hospitalization: true,
            life_threatening: a.severity == Severity::Severe && a.u_outcome < 0.15,
            ..SeriousCriteria::default()
        });
        let mut treated = false;
        if let Some(ind) = reference::indication(&a.pt, TreatmentKind::Acute) {
            if a.u_treat < (ind.p_treat * [0.7, 1.2, 1.5][sev_i]).min(1.0) {
                treated = true;
                let w: Vec<f64> = ind.ingredients.iter().map(|(_, w)| *w).collect();
                let pick = pick_weighted(&w, a.u_med);
                acute_meds.push(Medication {
                    ingredient: ind.ingredients[pick].0.clone(),
                    indication: a.pt.clone(),
                    start: a.start.to_string(),
                    end,
                });
            }
        }
        adverse_events.push(AdverseEvent {
            pt: a.pt.clone(),
            start: a.start,
            end,
            severity: a.severity,
            serious,
            causality,
            action,
            outcome,
            treated,
        });
    }
    if let Some(dd) = death {
        // The fatal event behind the death.
        let pt = [
            "Myocardial infarction",
            "Pneumonia",
            "Acute kidney injury",
            "Atrial fibrillation",
        ][pick_weighted(&[0.45, 0.3, 0.15, 0.1], fatal_pick)];
        let start = dd.add_days(-fatal_lead).max(s.consent);
        adverse_events.push(AdverseEvent {
            pt: pt.to_string(),
            start,
            end: Some(dd),
            severity: Severity::Severe,
            serious: Some(SeriousCriteria {
                death: true,
                hospitalization: true,
                ..SeriousCriteria::default()
            }),
            causality: "NOT RELATED",
            action: if start >= first_dose && start <= last_dose {
                "DRUG WITHDRAWN"
            } else {
                "NOT APPLICABLE"
            },
            outcome: "FATAL",
            treated: false,
        });
    }
    adverse_events.sort_by(|a, b| (a.start, &a.pt).cmp(&(b.start, &b.pt)));

    // 3. Visits actually attended.
    let mut visits: Vec<VisitOccurrence> = Vec::new();
    for (i, date, missed) in &planned_dates {
        let v = &schedule.visits[*i];
        let attended = match v.kind {
            VisitKind::Screening | VisitKind::Baseline => true,
            VisitKind::EndOfTreatment => early.is_none(),
            _ => !missed && *date < treatment_end,
        };
        if attended {
            visits.push(VisitOccurrence {
                visitnum: v.visitnum,
                date: *date,
                end: None,
            });
        }
    }
    if let (Some(fu), Some(fd)) = (fu_visit, fu_date) {
        if early.is_none() && death.is_none() {
            visits.push(VisitOccurrence {
                visitnum: fu.visitnum,
                date: fd,
                end: None,
            });
        }
    }
    let et_visit = schedule
        .visits
        .iter()
        .find(|v| v.kind == VisitKind::EarlyTermination);
    if let (Some((d, reason)), Some(et)) = (&early, et_visit) {
        let attends = match reason {
            Reason::Death => false,
            Reason::Dropout(why) => why != "LOST TO FOLLOW-UP",
            Reason::AdverseEvent(_) => true,
        };
        if attends && visits.last().is_none_or(|v| v.date < *d) {
            visits.push(VisitOccurrence {
                visitnum: et.visitnum,
                date: *d,
                end: None,
            });
        }
    }

    // 4. Dosing.
    let dispensing: Vec<Date> = visits
        .iter()
        .filter(|v| v.date > first_dose && v.date <= last_dose)
        .map(|v| v.date)
        .collect();
    let exposure = exposure_intervals(
        first_dose,
        last_dose,
        &interruptions,
        &dispensing,
        design.arms[s.arm].treatment.dose,
    );

    // Elements.
    let mut elements = vec![
        ElementSpan {
            epoch: Epoch::Screening,
            start: s.consent,
            end: first_dose,
        },
        ElementSpan {
            epoch: Epoch::Treatment,
            start: first_dose,
            end: treatment_end,
        },
    ];
    if last_contact > treatment_end {
        elements.push(ElementSpan {
            epoch: Epoch::FollowUp,
            start: treatment_end,
            end: last_contact,
        });
    }

    // Disposition records.
    let disposition = vec![
        DispositionEvent {
            decod: "INFORMED CONSENT OBTAINED".into(),
            kind: DispositionKind::Milestone,
            date: s.consent,
        },
        DispositionEvent {
            decod: "RANDOMIZED".into(),
            kind: DispositionKind::Milestone,
            date: s.randomization,
        },
        DispositionEvent {
            decod: eot_reason,
            kind: DispositionKind::EndOfTreatment,
            date: treatment_end,
        },
        DispositionEvent {
            decod: eos_reason,
            kind: DispositionKind::EndOfStudy,
            date: last_contact,
        },
    ];

    // 5. History, medications, findings.
    let history = medical_history(stream, design, p, s);
    let mut medications = chronic_medications(stream, p, s, &history);
    for m in acute_meds {
        if !medications
            .iter()
            .any(|x| x.ingredient == m.ingredient && x.start == m.start)
        {
            medications.push(m);
        }
    }
    let female = design.population.sex.levels[s.sex].0 == "F";
    let findings = findings(
        stream,
        schedule,
        p,
        s,
        female,
        &visits,
        first_dose,
        day(eot_visit.day.unwrap_or(1)),
        treatment_end,
    );

    Timeline {
        exposure,
        elements,
        visits,
        disposition,
        death,
        last_contact,
        treatment_end,
        adverse_events,
        medications,
        history,
        findings,
        lesions: Vec::new(),
        tumor: Vec::new(),
        pk: Vec::new(),
        clinical_events: Vec::new(),
    }
}

fn pick_weighted(weights: &[f64], u: f64) -> usize {
    let total: f64 = weights.iter().sum();
    let mut target = u * total;
    for (i, w) in weights.iter().enumerate() {
        if target < *w {
            return i;
        }
        target -= w;
    }
    weights.len() - 1
}

/// Daily dosing from `first` to `last`, minus interruptions, split at dispensing visits.
fn exposure_intervals(
    first: Date,
    last: Date,
    interruptions: &[(Date, Date)],
    splits: &[Date],
    dose: f64,
) -> Vec<Exposure> {
    let mut out: Vec<Exposure> = Vec::new();
    let mut start: Option<Date> = None;
    let mut d = first;
    while d <= last {
        let interrupted = interruptions.iter().any(|(a, b)| *a <= d && d <= *b);
        let split_here = splits.contains(&d);
        if let Some(st) = start {
            if interrupted || split_here {
                out.push(Exposure {
                    start: st,
                    end: d.add_days(-1),
                    dose,
                    time: None,
                    period: None,
                });
                start = None;
            }
        }
        if !interrupted && start.is_none() {
            start = Some(d);
        }
        d = d.add_days(1);
    }
    if let Some(st) = start {
        out.push(Exposure {
            start: st,
            end: last,
            dose,
            time: None,
            period: None,
        });
    }
    out
}

fn partial_date(d: Date, with_month: bool) -> String {
    let (y, m, _) = d.ymd();
    if with_month {
        format!("{y:04}-{m:02}")
    } else {
        format!("{y:04}")
    }
}

fn medical_history(
    stream: Stream,
    design: &Design,
    p: &CourseParams,
    s: &Subject,
) -> Vec<HistoryItem> {
    let female = design.population.sex.levels[s.sex].0 == "F";
    let mut out = Vec::new();
    for c in reference::history_conditions() {
        let mut d = Draws::new(stream.label("mh").label(&c.pt));
        let u = d.uniform();
        let years = d.uniform_range(0.5, c.max_years.min((s.age - 17) as f64).max(0.6));
        let with_month = d.bernoulli(0.6);
        let ongoing = d.bernoulli(if c.chronic { 0.92 } else { 0.0 });
        let required = p.required_history.contains(&c.pt);
        let prev = if female {
            c.prevalence.0
        } else {
            c.prevalence.1
        };
        let present = required || {
            prev > 0.0 && {
                let logit = math::ln(prev / (1.0 - prev))
                    + c.log_odds_per_decade * (s.age as f64 - 55.0) / 10.0;
                u < 1.0 / (1.0 + math::exp(-logit))
            }
        };
        if present {
            let onset = s
                .consent
                .add_days(-((years * DAYS_PER_YEAR) as i32).max(30));
            out.push(HistoryItem {
                pt: c.pt.clone(),
                category: c.category.clone(),
                start: partial_date(onset, with_month),
                ongoing: ongoing || required,
            });
        }
    }
    out
}

fn chronic_medications(
    stream: Stream,
    p: &CourseParams,
    s: &Subject,
    history: &[HistoryItem],
) -> Vec<Medication> {
    let mut out: Vec<Medication> = Vec::new();
    for h in history.iter().filter(|h| h.ongoing) {
        let mut d = Draws::new(stream.label("cm").label(&h.pt));
        let u = d.uniform();
        let u_pick = d.uniform();
        let lag_years = d.uniform_range(0.0, 1.0);
        let background = p
            .background_medication
            .iter()
            .find(|(ind, _, _)| *ind == h.pt);
        let chosen = match (
            background,
            reference::indication(&h.pt, TreatmentKind::Chronic),
        ) {
            (Some((_, ing, prob)), _) => (u < *prob).then(|| ing.clone()),
            (None, Some(ind)) if u < ind.p_treat => {
                let w: Vec<f64> = ind.ingredients.iter().map(|(_, w)| *w).collect();
                Some(ind.ingredients[pick_weighted(&w, u_pick)].0.clone())
            }
            _ => None,
        };
        let Some(ingredient) = chosen else { continue };
        if out.iter().any(|m| m.ingredient == ingredient) {
            continue;
        }
        // Started some time after the condition's onset, before consent.
        let onset_year: i32 = h.start[..4].parse().expect("history start has a year");
        let (consent_year, _, _) = s.consent.ymd();
        let start_year = onset_year + ((consent_year - onset_year) as f64 * lag_years) as i32;
        out.push(Medication {
            ingredient,
            indication: h.pt.clone(),
            start: format!("{:04}", start_year.min(consent_year)),
            end: None,
        });
    }
    out
}

/// A subject's analyte: everything needed to produce its value on a date.
struct SubjectAnalyte {
    testcd: String,
    domain: FindingDomain,
    scale: Scale,
    level: f64,    // subject mean on the model scale
    pop_mean: f64, // population mean on the model scale
    sd_between: f64,
    sd_within: f64,
    effect: Option<f64>,
    onset: f64,
    floor: f64,
}

#[allow(clippy::too_many_arguments)]
fn findings(
    stream: Stream,
    schedule: &Schedule,
    p: &CourseParams,
    s: &Subject,
    female: bool,
    visits: &[VisitOccurrence],
    first_dose: Date,
    planned_eot: Date,
    treatment_end: Date,
) -> Vec<Finding> {
    let fd = stream.label("findings");
    let group_z = |g: &str| Draws::new(fd.label("group").label(g)).normal(0.0, 1.0);
    let rho: f64 = 0.6;
    let mut analytes: Vec<SubjectAnalyte> = Vec::new();

    let effect_for = |testcd: &str| p.endpoints.iter().find(|e| e.testcd == testcd);
    let override_for = |testcd: &str| p.analyte_overrides.iter().find(|o| o.testcd == testcd);

    for a in reference::lab_panel() {
        let mut d = Draws::new(fd.label(&a.testcd));
        let ov = override_for(&a.testcd);
        let mean_pair = ov.and_then(|o| o.mean).unwrap_or(a.mean);
        let mean = if female { mean_pair.0 } else { mean_pair.1 };
        let sd_b = ov.and_then(|o| o.sd_between).unwrap_or(a.sd_between);
        let zg = a.group.as_deref().map_or(0.0, group_z);
        let zo = d.normal(0.0, 1.0);
        let z = if a.group.is_some() {
            math::sqrt(rho) * zg + math::sqrt(1.0 - rho) * zo
        } else {
            zo
        };
        let pop_mean = match a.scale {
            Scale::Linear => mean,
            Scale::Log => math::ln(mean),
        };
        let level = pop_mean + sd_b * z;
        let e = effect_for(&a.testcd);
        analytes.push(SubjectAnalyte {
            testcd: a.testcd.clone(),
            domain: FindingDomain::Lb,
            scale: a.scale,
            level,
            pop_mean,
            sd_between: sd_b,
            sd_within: a.sd_within,
            effect: e.map(|e| e.change_by_arm[s.arm]),
            onset: e.map_or(1.0, |e| e.onset_days),
            floor: 0.05 * mean,
        });
    }
    // Vital signs; weight follows from BMI and height.
    let mut vd = Draws::new(fd.label("anthropometry"));
    let bmi = vd.normal(p.bmi.0, p.bmi.1).clamp(18.0, 55.0);
    let mut height = 170.0;
    for v in reference::vital_signs() {
        let mut d = Draws::new(fd.label(&v.testcd));
        let ov = override_for(&v.testcd);
        let mean_pair = ov.and_then(|o| o.mean).unwrap_or(v.mean);
        let mean = if female { mean_pair.0 } else { mean_pair.1 };
        let sd_b = ov.and_then(|o| o.sd_between).unwrap_or(v.sd_between);
        let zg = v.group.as_deref().map_or(0.0, group_z);
        let zo = d.normal(0.0, 1.0);
        let z = if v.group.is_some() {
            math::sqrt(rho) * zg + math::sqrt(1.0 - rho) * zo
        } else {
            zo
        };
        let level = match v.testcd.as_str() {
            "HEIGHT" => {
                height = (mean + sd_b * z).clamp(140.0, 210.0);
                height
            }
            "WEIGHT" => bmi * (height / 100.0) * (height / 100.0),
            _ => mean + sd_b * z,
        };
        let e = effect_for(&v.testcd);
        analytes.push(SubjectAnalyte {
            testcd: v.testcd.clone(),
            domain: FindingDomain::Vs,
            scale: Scale::Linear,
            level,
            pop_mean: mean,
            sd_between: sd_b,
            sd_within: v.sd_within,
            effect: e.map(|e| e.change_by_arm[s.arm]),
            onset: e.map_or(1.0, |e| e.onset_days),
            floor: 0.3 * mean.max(1.0),
        });
    }

    // Inclusion criteria on screening values: re-draw the subject level until it qualifies.
    let screening = visits.iter().find(|v| {
        schedule
            .visit(v.visitnum)
            .is_some_and(|pv| pv.kind == VisitKind::Screening)
    });
    for a in analytes.iter_mut() {
        let Some((lo, hi)) = override_for(&a.testcd).and_then(|o| o.inclusion) else {
            continue;
        };
        let Some(scr) = screening else { continue };
        // The same draw the screening measurement uses below.
        let e_scr =
            Draws::new(fd.label(&a.testcd).index(scr.visitnum as u64)).normal(0.0, a.sd_within);
        let mut d = Draws::new(fd.label("inclusion").label(&a.testcd));
        for _ in 0..200 {
            if (lo..=hi).contains(&(a.level + e_scr)) {
                break;
            }
            a.level = a.pop_mean + a.sd_between * d.normal(0.0, 1.0);
        }
        a.level = a.level.clamp(lo - e_scr, hi - e_scr);
    }

    let total = f64::from(planned_eot.days_since(first_dose).max(1));
    let mut out = Vec::new();
    for v in visits {
        let pv: &PlannedVisit = schedule
            .visit(v.visitnum)
            .expect("attended visits are planned");
        let t = f64::from(v.date.days_since(first_dose).max(0));
        for a in &analytes {
            let wanted = match a.domain {
                FindingDomain::Lb => pv.labs || (pv.endpoint && a.effect.is_some()),
                FindingDomain::Vs => {
                    (pv.vitals && a.testcd != "HEIGHT")
                        || (a.testcd == "HEIGHT" && pv.kind == VisitKind::Screening)
                }
                FindingDomain::Ft => false,
            };
            if !wanted {
                continue;
            }
            let e =
                Draws::new(fd.label(&a.testcd).index(v.visitnum as u64)).normal(0.0, a.sd_within);
            let effect = a.effect.map_or(0.0, |delta| {
                let shape = |x: f64| 1.0 - math::exp(-x / a.onset);
                let full = shape(total);
                if v.date <= treatment_end {
                    delta * shape(t) / full
                } else {
                    let since = f64::from(v.date.days_since(treatment_end));
                    delta * shape(f64::from(treatment_end.days_since(first_dose))) / full
                        * math::exp(-since / a.onset)
                }
            });
            let value = match a.scale {
                Scale::Linear => a.level + effect + e,
                Scale::Log => math::exp(a.level + e),
            };
            let value = if a.testcd == "HEIGHT" {
                a.level
            } else {
                value.max(a.floor)
            };
            out.push(Finding {
                domain: a.domain,
                testcd: a.testcd.clone(),
                visitnum: v.visitnum,
                date: v.date,
                value,
            });
        }
        // Diastolic stays below systolic.
        let n = out.len();
        let sys = out[..n]
            .iter()
            .rposition(|f| f.testcd == "SYSBP" && f.date == v.date);
        let dia = out[..n]
            .iter()
            .rposition(|f| f.testcd == "DIABP" && f.date == v.date);
        if let (Some(si), Some(di)) = (sys, dia) {
            out[di].value = out[di].value.min(out[si].value - 10.0);
        }
    }
    out
}
