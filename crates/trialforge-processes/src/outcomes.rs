//! An event-driven cardiovascular outcomes course.
//!
//! Per subject, from randomisation: recurrent non-fatal events (e.g. myocardial infarction,
//! ischaemic stroke) and terminal events (cardiovascular death, and non-cardiovascular death as
//! a competing risk), each a Poisson process with an annual rate × the arm's hazard ratio × the
//! stratum's multiplier × a subject frailty (log-normal, shared across event types), plus
//! withdrawal and permanent treatment discontinuation.
//!
//! The data cutoff is event-driven: the date of the K-th first MACE (any event flagged MACE)
//! across all subjects, or the planned maximum when fewer occur. Subjects stay in follow-up after
//! stopping treatment (intention to treat), with visits every few months.

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Epoch, Subject};

use crate::timeline::*;

const DAYS_PER_YEAR: f64 = 365.25;

/// One clinical event type.
#[derive(Clone, Debug, PartialEq)]
pub struct EventType {
    /// Preferred Term (`CEDECOD`).
    pub pt: String,
    /// Annual rate in the reference arm (stratum multiplier 1, frailty 1).
    pub rate: f64,
    /// Hazard ratio per arm.
    pub hr: Vec<f64>,
    /// Whether the event ends follow-up (death).
    pub terminal: bool,
    /// Whether the event counts toward the MACE composite (and the event-driven cutoff).
    pub mace: bool,
}

/// Parameters of the outcomes course.
#[derive(Clone, Debug, PartialEq)]
pub struct OutcomesParams {
    /// Event types.
    pub events: Vec<EventType>,
    /// Hazard multiplier per stratum.
    pub stratum_hr: Vec<f64>,
    /// SD of the log frailty.
    pub frailty_sd: f64,
    /// Annual hazard of withdrawing (or being lost to follow-up).
    pub withdrawal_rate: f64,
    /// Withdrawal reasons (`DSDECOD`) and weights.
    pub withdrawal_reasons: Vec<(String, f64)>,
    /// Annual hazard of stopping study treatment for good.
    pub discontinuation_rate: f64,
    /// Days between scheduled visits.
    pub visit_interval: i32,
    /// First MACE events that trigger the cutoff.
    pub target_events: usize,
    /// An interim cutoff at this many first MACE events (a snapshot of the same study, which
    /// may fall during enrollment).
    pub interim_events: Option<usize>,
    /// Longest possible follow-up from the first randomisation, in days.
    pub max_days: i32,
}

/// A subject's latent course.
struct Latent {
    events: Vec<(Date, usize)>,
    withdrawal: Option<Date>,
    stop: Option<Date>,
    reason: String,
    contact_lag: i32,
}

fn latent(design: &Design, p: &OutcomesParams, s: &Subject) -> Latent {
    let stream = design.subject_stream(s.index).label("outcomes");
    let mut d = Draws::new(stream.label("subject"));
    // Draw order is part of the output contract; append new draws at the end.
    let frailty = math::exp(d.normal(0.0, p.frailty_sd) - p.frailty_sd * p.frailty_sd / 2.0);
    let withdraw_years = d.exponential(p.withdrawal_rate.max(1e-12));
    let stop_years = d.exponential(p.discontinuation_rate.max(1e-12));
    let w: Vec<f64> = p.withdrawal_reasons.iter().map(|(_, w)| *w).collect();
    let reason = p.withdrawal_reasons[d.weighted(&w)].0.clone();
    let contact_lag = d.int_range(0, 30) as i32;
    let horizon = f64::from(p.max_days) / DAYS_PER_YEAR;
    let rand = s.randomization;
    let day = |years: f64| rand.add_days((years * DAYS_PER_YEAR).floor().min(1e6) as i32);
    let mut events: Vec<(Date, usize)> = Vec::new();
    for (k, e) in p.events.iter().enumerate() {
        let mut de = Draws::new(stream.label("event").label(&e.pt));
        let rate = e.rate * e.hr[s.arm] * p.stratum_hr[s.stratum] * frailty;
        let mut t = 0.0;
        loop {
            t += de.exponential(rate.max(1e-12));
            if t > horizon {
                break;
            }
            let date = day(t).max(rand.add_days(1));
            // A recurrence on the same day is the same clinical event.
            if events.last().is_none_or(|(d, kk)| *kk != k || *d < date) {
                events.push((date, k));
            }
            if e.terminal {
                break;
            }
        }
    }
    events.sort();
    // Nothing after the first death.
    if let Some(i) = events.iter().position(|(_, k)| p.events[*k].terminal) {
        events.truncate(i + 1);
    }
    Latent {
        events,
        withdrawal: (withdraw_years < horizon).then(|| day(withdraw_years).max(rand.add_days(1))),
        stop: (stop_years < horizon).then(|| day(stop_years).max(s.treatment_start)),
        reason,
        contact_lag,
    }
}

/// The event-driven cutoff: the date of the `target_events`-th first MACE (before each subject's
/// death or withdrawal), or `fallback` when fewer occur.
fn cutoff(p: &OutcomesParams, latents: &[Latent], fallback: Date) -> Date {
    let mut first: Vec<Date> = latents
        .iter()
        .filter_map(|l| {
            let end = l.withdrawal.unwrap_or(Date::from_days(i32::MAX / 2));
            l.events
                .iter()
                .find(|(d, k)| p.events[*k].mace && *d <= end)
                .map(|(d, _)| *d)
        })
        .collect();
    first.sort();
    match first.get(p.target_events.saturating_sub(1)) {
        Some(&d) if p.target_events > 0 => d.min(fallback),
        _ => fallback,
    }
}

/// Simulates every subject; returns the timelines and the data cutoff.
pub fn simulate(
    design: &Design,
    p: &OutcomesParams,
    subjects: &[Subject],
) -> (Vec<Timeline>, Date) {
    let latents: Vec<Latent> = subjects.iter().map(|s| latent(design, p, s)).collect();
    let last_rand = subjects
        .iter()
        .map(|s| s.randomization)
        .max()
        .unwrap_or(design.first_randomization);
    let cut = match p.interim_events {
        Some(k) => cutoff(
            &OutcomesParams {
                target_events: k,
                ..p.clone()
            },
            &latents,
            design.cutoff,
        ),
        None => cutoff(p, &latents, design.cutoff).max(last_rand.add_days(1)),
    };
    let timelines = subjects
        .iter()
        .zip(&latents)
        .map(|(s, l)| observe(design, p, s, l, cut))
        .collect();
    (timelines, cut)
}

fn observe(design: &Design, p: &OutcomesParams, s: &Subject, l: &Latent, cutoff: Date) -> Timeline {
    let death = l
        .events
        .iter()
        .find(|(d, k)| p.events[*k].terminal && *d <= cutoff && l.withdrawal.is_none_or(|w| *d < w))
        .map(|(d, _)| *d);
    let withdrawal = l
        .withdrawal
        .filter(|&w| w <= cutoff && death.is_none_or(|dd| w < dd));
    // Events up to the cutoff are ascertained at close-out, so the last contact is no earlier
    // than the last of them.
    let last_event = l
        .events
        .iter()
        .map(|(d, _)| *d)
        .filter(|d| *d <= cutoff && withdrawal.is_none_or(|w| *d <= w))
        .max();
    let last_contact = death.or(withdrawal).unwrap_or_else(|| {
        cutoff
            .add_days(-l.contact_lag)
            .max(s.treatment_start)
            .max(last_event.unwrap_or(s.treatment_start))
    });
    // Treatment stops at discontinuation, death or withdrawal; otherwise it is ongoing.
    let stop = [l.stop, death, withdrawal]
        .into_iter()
        .flatten()
        .filter(|d| *d <= last_contact)
        .min();
    let first_dose = s.treatment_start;
    let treatment_end = stop.unwrap_or(last_contact).max(first_dose);

    // Visits: screening, randomisation (Day 1), then every `visit_interval` days.
    let schedule = design
        .schedule
        .as_ref()
        .expect("outcomes designs have a schedule");
    let mut visits = Vec::new();
    for pv in &schedule.visits {
        let Some(day) = pv.day else { continue };
        let date = if pv.kind == trialforge_model::VisitKind::Screening {
            s.screening
        } else {
            first_dose.add_days(day - 1)
        };
        if date <= last_contact {
            visits.push(VisitOccurrence {
                visitnum: pv.visitnum,
                date,
                end: None,
            });
        }
    }

    // Dispensing intervals between visits while on treatment.
    let mut exposure = Vec::new();
    let mut starts: Vec<Date> = visits
        .iter()
        .map(|v| v.date)
        .filter(|d| *d >= first_dose && *d < treatment_end)
        .collect();
    starts.dedup();
    if starts.first() != Some(&first_dose) {
        starts.insert(0, first_dose);
    }
    let dose = design.arms[s.arm].treatment.dose;
    for (k, &st) in starts.iter().enumerate() {
        let end = starts
            .get(k + 1)
            .map_or(treatment_end, |n| n.add_days(-1))
            .max(st);
        exposure.push(Exposure {
            start: st,
            end,
            dose,
            time: None,
            period: None,
            location: None,
            group: None,
            visitnum: None,
            switched: false,
        });
    }

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
    let mut disposition = vec![
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
    ];
    if let Some(st) = stop {
        let why = if Some(st) == death {
            "DEATH".to_string()
        } else if Some(st) == withdrawal {
            l.reason.clone()
        } else {
            "ADVERSE EVENT".to_string()
        };
        disposition.push(DispositionEvent {
            decod: why,
            kind: DispositionKind::EndOfTreatment,
            date: st,
        });
    }
    match (death, withdrawal) {
        (Some(dd), _) => disposition.push(DispositionEvent {
            decod: "DEATH".into(),
            kind: DispositionKind::EndOfStudy,
            date: dd,
        }),
        (None, Some(w)) => disposition.push(DispositionEvent {
            decod: l.reason.clone(),
            kind: DispositionKind::EndOfStudy,
            date: w,
        }),
        (None, None) => {}
    }
    let clinical_events = l
        .events
        .iter()
        .filter(|(d, _)| *d <= last_contact)
        .map(|(d, k)| ClinicalEvent {
            pt: p.events[*k].pt.clone(),
            term: p.events[*k].pt.clone(),
            category: "ADJUDICATED ENDPOINT EVENT",
            date: *d,
            end: None,
            severity: None,
        })
        .collect();

    Timeline {
        exposure,
        elements,
        visits,
        disposition,
        death,
        last_contact,
        treatment_end,
        adverse_events: Vec::new(),
        medications: Vec::new(),
        history: Vec::new(),
        findings: Vec::new(),
        lesions: Vec::new(),
        tumor: Vec::new(),
        pk: Vec::new(),
        clinical_events,
        study_eye: None,
        birth_date: None,
        diaries: Vec::new(),
    }
}
