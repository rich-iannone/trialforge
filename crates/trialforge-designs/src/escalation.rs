//! `dose_escalation`: an open-label, first-in-patient dose-escalation study of Drug A in
//! advanced solid tumours. Sequential cohorts are treated at dose levels chosen by a 3+3,
//! BOIN or mTPI-2 design from the dose-limiting toxicities (DLTs) in the first 28-day cycle;
//! `n` is the maximum sample size. DLTs are AEs flagged in SUPPAE (`AEDLT`); ADSL has the
//! DLT-evaluable population and ADAE the DLT flag.
//!
//! The escalation is decided from each subject's latent DLT tolerance on their own stream
//! (`trialforge_processes::escalation`), so it is part of the design and reproducible.

use serde::{Deserialize, Serialize};
use trialforge_model::{Cohort, Design, Epoch, PlannedVisit, Schedule, VisitKind};
use trialforge_processes::escalation::{EscalationParams, Method, Trace, escalate, tolerance};
use trialforge_processes::fixed::CourseParams;

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `dose_escalation`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// `3+3`, `boin` or `mtpi2`.
    pub method: String,
    /// Dose levels (mg once daily), in escalation order.
    pub doses: Vec<f64>,
    /// True probability of a DLT at each dose level.
    pub dlt_rates: Vec<f64>,
    /// Target DLT rate (BOIN, mTPI-2).
    pub target: f64,
    /// Subjects per cohort (BOIN, mTPI-2; the 3+3 treats 3).
    pub cohort_size: u32,
    /// Treatment duration in weeks (a multiple of 4; 28-day cycles).
    pub weeks: u32,
    /// Annual dropout hazard (reasons other than AEs and death).
    pub dropout_rate: f64,
    /// Multiplier on every background AE rate.
    pub ae_rate_multiplier: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            method: "boin".into(),
            doses: vec![10.0, 20.0, 40.0, 80.0, 160.0, 320.0],
            dlt_rates: vec![0.02, 0.05, 0.10, 0.25, 0.40, 0.55],
            target: 0.25,
            cohort_size: 3,
            weeks: 8,
            dropout_rate: 0.1,
            ae_rate_multiplier: 1.5,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new("method", "Dose-finding design: 3+3, boin or mtpi2"),
        ControlDoc::new("doses", "Dose levels in mg once daily, in escalation order"),
        ControlDoc::new(
            "dlt_rates",
            "True probability of a DLT in cycle 1 at each dose level",
        ),
        ControlDoc::new("target", "Target DLT rate (BOIN, mTPI-2)"),
        ControlDoc::new(
            "cohort_size",
            "Subjects per cohort (BOIN, mTPI-2; the 3+3 treats 3); n is the maximum sample size",
        ),
        ControlDoc::new(
            "weeks",
            "Treatment duration in weeks (a multiple of 4; 28-day cycles)",
        ),
        ControlDoc::new(
            "dropout_rate",
            "Annual dropout hazard (reasons other than adverse events and death)",
        ),
        ControlDoc::new(
            "ae_rate_multiplier",
            "Multiplier on every background adverse event rate",
        ),
    ]
}

fn method(c: &Controls) -> Option<Method> {
    match c.method.as_str() {
        "3+3" => Some(Method::ThreeThree),
        "boin" => Some(Method::Boin),
        "mtpi2" => Some(Method::Mtpi2),
        _ => None,
    }
}

pub(crate) fn validate(c: &Controls, n: usize) -> Result<(), String> {
    if method(c).is_none() {
        return Err("`method` must be \"3+3\", \"boin\" or \"mtpi2\"".into());
    }
    if c.doses.is_empty() || c.doses.len() > 10 {
        return Err("`doses` needs 1 to 10 dose levels".into());
    }
    if !c.doses.windows(2).all(|w| w[0] < w[1])
        || !c.doses.iter().all(|d| d.is_finite() && *d > 0.0)
    {
        return Err("`doses` must be positive and increasing".into());
    }
    if c.dlt_rates.len() != c.doses.len() || !c.dlt_rates.iter().all(|p| (0.0..=1.0).contains(p)) {
        return Err("`dlt_rates` needs one probability in [0, 1] per dose".into());
    }
    if !c.dlt_rates.windows(2).all(|w| w[0] <= w[1]) {
        return Err("`dlt_rates` must not decrease with dose".into());
    }
    if !(0.05..=0.5).contains(&c.target) {
        return Err("`target` must be between 0.05 and 0.5".into());
    }
    if !(1..=6).contains(&c.cohort_size) {
        return Err("`cohort_size` must be between 1 and 6".into());
    }
    if !(3..=300).contains(&n) {
        return Err("`n` (the maximum sample size) must be between 3 and 300".into());
    }
    if !(4..=52).contains(&c.weeks) || !c.weeks.is_multiple_of(4) {
        return Err("`weeks` must be a multiple of 4 between 4 and 52".into());
    }
    if !(0.0..=5.0).contains(&c.dropout_rate) || !(0.0..=20.0).contains(&c.ae_rate_multiplier) {
        return Err("`dropout_rate` in [0, 5] and `ae_rate_multiplier` in [0, 20]".into());
    }
    Ok(())
}

/// Weekly visits in cycle 1 (DLT monitoring), day 1 of each later cycle, the end of treatment
/// and a safety follow-up 30 days later.
fn schedule(weeks: u32) -> Schedule {
    let visit = |num: i64, name: String, day: i32, window, epoch, kind, labs| PlannedVisit {
        visitnum: num,
        name,
        day: Some(day),
        window,
        epoch,
        kind,
        labs,
        vitals: true,
        endpoint: true,
    };
    let mut visits = vec![
        visit(
            1,
            "SCREENING".into(),
            -14,
            (14, 0),
            Epoch::Screening,
            VisitKind::Screening,
            true,
        ),
        visit(
            2,
            "CYCLE 1 DAY 1".into(),
            1,
            (0, 0),
            Epoch::Treatment,
            VisitKind::Baseline,
            true,
        ),
    ];
    for d in [8, 15, 22] {
        let n = visits.len() as i64 + 1;
        visits.push(visit(
            n,
            format!("CYCLE 1 DAY {d}"),
            d,
            (1, 1),
            Epoch::Treatment,
            VisitKind::Treatment,
            true,
        ));
    }
    let cycles = weeks / 4;
    for k in 2..=cycles {
        let n = visits.len() as i64 + 1;
        visits.push(visit(
            n,
            format!("CYCLE {k} DAY 1"),
            (28 * (k - 1) + 1) as i32,
            (2, 2),
            Epoch::Treatment,
            VisitKind::Treatment,
            true,
        ));
    }
    let eot = (weeks * 7 + 1) as i32;
    let n = visits.len() as i64 + 1;
    visits.push(visit(
        n,
        "END OF TREATMENT".into(),
        eot,
        (3, 3),
        Epoch::Treatment,
        VisitKind::EndOfTreatment,
        true,
    ));
    visits.push(visit(
        n + 1,
        "SAFETY FOLLOW-UP".into(),
        eot + 30,
        (5, 5),
        Epoch::FollowUp,
        VisitKind::FollowUp,
        true,
    ));
    visits.push(PlannedVisit {
        visitnum: 99,
        name: "EARLY TERMINATION".into(),
        day: None,
        window: (0, 0),
        epoch: Epoch::Treatment,
        kind: VisitKind::EarlyTermination,
        labs: true,
        vitals: true,
        endpoint: true,
    });
    Schedule {
        visits,
        treatment_days: eot - 1,
        followup_days: 30,
    }
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

pub(crate) fn resolve(
    opts: &StudyOptions,
    c: &Controls,
) -> (Design, CourseParams, EscalationParams, Trace) {
    let arms: Vec<ArmSpec> = c
        .doses
        .iter()
        .enumerate()
        .map(|(k, &d)| ArmSpec {
            armcd: format!("DL{}", k + 1),
            arm: format!("Drug A {} mg", fmt_dose(d)),
            treatment: "DRUG A".into(),
            dose: d,
            ratio: 1,
            change: 0.0,
            ae_share: 0.0,
        })
        .collect();
    let (mut design, mut course) = build(
        opts,
        &Common {
            indication: "advanced_solid_tumors",
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: 0.02,
            studyid: "SIMESC01",
            phase: "PHASE I TRIAL",
            control: "NONE",
            title: None,
        },
        &arms,
    );
    let p = EscalationParams {
        method: method(c).expect("validated"),
        dlt_rates: c.dlt_rates.clone(),
        target: c.target,
        cohort_size: c.cohort_size as usize,
        max_n: opts.n,
    };
    let trace = escalate(&p, |i, level| tolerance(&design, i) < p.dlt_rates[level]);
    // Each cohort enrols after the previous one has completed its DLT window.
    design.cohorts = trace
        .cohorts
        .iter()
        .enumerate()
        .map(|(k, &(level, size))| Cohort {
            arms: vec![(level, size as u32)],
            start_day: 35 * k as i32,
            enrollment_days: 7,
        })
        .collect();
    design.n = trace.cohorts.iter().map(|c| c.1).sum();
    design.enrollment_days = 35 * trace.cohorts.len() as i32;
    design.cutoff = design
        .first_randomization
        .add_days(design.enrollment_days + (c.weeks * 7) as i32 + 120);
    design.schedule = Some(schedule(c.weeks));
    design.population.age_min = 18;
    design.population.age_max = 85;
    let info = &mut design.info;
    info.blinding = "OPEN LABEL".into();
    info.trial_types = vec!["SAFETY".into(), "DOSE FINDING".into()];
    info.intervention_model = "SEQUENTIAL".into();
    course.death_annual_rate = 0.05;
    (design, course, p, trace)
}
