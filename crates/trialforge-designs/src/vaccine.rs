//! `vaccine`: a randomised, observer-blind, placebo-controlled study of a two-dose (or
//! single-dose) vaccine in healthy adults. Neutralising antibody titres (IS) are measured
//! before each dose, 28 days after the last dose and at the end of follow-up; solicited local
//! and systemic reactions are recorded in a 7-day diary after each vaccination (FACE, with
//! temperatures in VS and a summary per reaction in CE). ADIS, ADCE and ADFACE follow
//! admiralvaccine's templates.
//!
//! Vaccinations, diaries and titres come from `trialforge_processes::vaccine`.

use serde::{Deserialize, Serialize};
use trialforge_model::{Categorical, Design, Epoch, PlannedVisit, Schedule, VisitKind};
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::vaccine::{Solicited, VaccineParams};

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `vaccine`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Vaccinations in the primary series (1 or 2).
    pub doses: u32,
    /// Days between vaccinations.
    pub interval_days: u32,
    /// Geometric mean neutralising titre 28 days after the last dose (vaccine arm).
    pub gmt: f64,
    /// Share of subjects seropositive at baseline.
    pub seropositive: f64,
    /// Share of the log-titre rise reached after the first of two doses.
    pub first_dose_fraction: f64,
    /// SD of individual log10 titres after vaccination.
    pub titer_sd: f64,
    /// Multiplier on the vaccine arm's solicited reaction probabilities.
    pub reactogenicity: f64,
    /// Probability of fever (≥ 38.0 °C) after a vaccination (vaccine arm).
    pub fever_rate: f64,
    /// Follow-up after the last dose, in days.
    pub followup_days: u32,
    /// Age range in years: [youngest, oldest].
    pub age_range: [i64; 2],
    /// Allocation ratio: [placebo, vaccine].
    pub ratio: [u32; 2],
    /// Annual dropout hazard (reasons other than AEs and death).
    pub dropout_rate: f64,
    /// Multiplier on every background (unsolicited) AE rate.
    pub ae_rate_multiplier: f64,
    /// Probability an optional visit is missed.
    pub missed_visit: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            doses: 2,
            interval_days: 28,
            gmt: 640.0,
            seropositive: 0.1,
            first_dose_fraction: 0.6,
            titer_sd: 0.45,
            reactogenicity: 1.0,
            fever_rate: 0.08,
            followup_days: 180,
            age_range: [18, 64],
            ratio: [1, 1],
            dropout_rate: 0.05,
            ae_rate_multiplier: 1.0,
            missed_visit: 0.02,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new("doses", "Vaccinations in the primary series (1 or 2)"),
        ControlDoc::new("interval_days", "Days between vaccinations (21 to 84)"),
        ControlDoc::new(
            "gmt",
            "Geometric mean neutralising titre 28 days after the last dose (vaccine arm)",
        ),
        ControlDoc::new("seropositive", "Share of subjects seropositive at baseline"),
        ControlDoc::new(
            "first_dose_fraction",
            "Share of the log-titre rise reached after the first of two doses",
        ),
        ControlDoc::new(
            "titer_sd",
            "SD of individual log10 titres after vaccination",
        ),
        ControlDoc::new(
            "reactogenicity",
            "Multiplier on the vaccine arm's solicited reaction probabilities",
        ),
        ControlDoc::new(
            "fever_rate",
            "Probability of fever (38.0 C or more) after a vaccination, vaccine arm",
        ),
        ControlDoc::new("followup_days", "Follow-up after the last dose, in days"),
        ControlDoc::new("age_range", "Age range in years: [youngest, oldest]"),
        ControlDoc::new("ratio", "Allocation ratio: [placebo, vaccine]"),
        ControlDoc::new(
            "dropout_rate",
            "Annual dropout hazard (reasons other than adverse events and death)",
        ),
        ControlDoc::new(
            "ae_rate_multiplier",
            "Multiplier on every background (unsolicited) adverse event rate",
        ),
        ControlDoc::new("missed_visit", "Probability an optional visit is missed"),
    ]
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    if !(1..=2).contains(&c.doses) {
        return Err("`doses` must be 1 or 2".into());
    }
    if !(21..=84).contains(&c.interval_days) {
        return Err("`interval_days` must be between 21 and 84".into());
    }
    if !(c.gmt.is_finite() && (10.0..=10_000.0).contains(&c.gmt)) {
        return Err("`gmt` must be between 10 and 10000".into());
    }
    let unit = |x: f64| (0.0..=1.0).contains(&x);
    if !unit(c.seropositive) || !unit(c.first_dose_fraction) || !unit(c.fever_rate) {
        return Err(
            "`seropositive`, `first_dose_fraction` and `fever_rate` must be in [0, 1]".into(),
        );
    }
    if !(0.05..=1.5).contains(&c.titer_sd) || !(0.0..=3.0).contains(&c.reactogenicity) {
        return Err("`titer_sd` must be in [0.05, 1.5] and `reactogenicity` in [0, 3]".into());
    }
    if !(56..=730).contains(&c.followup_days) {
        return Err("`followup_days` must be between 56 and 730".into());
    }
    if c.age_range[0] < 18 || c.age_range[1] > 90 || c.age_range[0] > c.age_range[1] {
        return Err("`age_range` must be [youngest, oldest] within 18 to 90 years".into());
    }
    if c.ratio.contains(&0) || c.ratio.iter().any(|&r| r > 10) {
        return Err("`ratio` entries must be between 1 and 10".into());
    }
    if !(0.0..=5.0).contains(&c.dropout_rate)
        || !(0.0..=20.0).contains(&c.ae_rate_multiplier)
        || !(0.0..=0.5).contains(&c.missed_visit)
    {
        return Err(
            "`dropout_rate` in [0, 5], `ae_rate_multiplier` in [0, 20], `missed_visit` in [0, 0.5]"
                .into(),
        );
    }
    Ok(())
}

/// Planned study days: (vaccinations, end of the vaccination phase, end of follow-up).
fn days(c: &Controls) -> (Vec<i32>, i32, i32) {
    let vax: Vec<i32> = (0..c.doses)
        .map(|k| 1 + (k * c.interval_days) as i32)
        .collect();
    let last = *vax.last().expect("at least one dose");
    (vax, last + 28, last + c.followup_days as i32)
}

/// Screening, the vaccination visits, 28 days after the last dose (end of the vaccination
/// phase) and the end of follow-up.
fn schedule(c: &Controls) -> Schedule {
    let (vax, eot, fu) = days(c);
    let visit = |num: i64, day: i32, window, epoch, kind, labs| PlannedVisit {
        visitnum: num,
        name: if kind == VisitKind::Screening {
            "SCREENING".into()
        } else {
            format!("DAY {day}")
        },
        day: Some(day),
        window,
        epoch,
        kind,
        labs,
        vitals: true,
        endpoint: true,
    };
    let mut visits = vec![visit(
        1,
        -7,
        (7, 0),
        Epoch::Screening,
        VisitKind::Screening,
        true,
    )];
    for (k, &d) in vax.iter().enumerate() {
        let kind = if k == 0 {
            VisitKind::Baseline
        } else {
            VisitKind::Treatment
        };
        visits.push(visit(
            2 + k as i64,
            d,
            if k == 0 { (0, 0) } else { (3, 3) },
            Epoch::Treatment,
            kind,
            false,
        ));
    }
    let n = visits.len() as i64;
    visits.push(visit(
        n + 1,
        eot,
        (3, 3),
        Epoch::Treatment,
        VisitKind::EndOfTreatment,
        true,
    ));
    visits.push(visit(
        n + 2,
        fu,
        (14, 14),
        Epoch::FollowUp,
        VisitKind::FollowUp,
        false,
    ));
    visits.push(PlannedVisit {
        visitnum: 99,
        name: "EARLY TERMINATION".into(),
        day: None,
        window: (0, 0),
        epoch: Epoch::Treatment,
        kind: VisitKind::EarlyTermination,
        labs: false,
        vitals: true,
        endpoint: true,
    });
    Schedule {
        visits,
        treatment_days: eot - 1,
        followup_days: fu - eot,
    }
}

/// Solicited reactions: (term, Preferred Term, site, measured, placebo and vaccine rates).
const SOLICITED: [(&str, &str, bool, bool, f64, f64); 10] = [
    (
        "PAIN AT INJECTION SITE",
        "Injection site pain",
        true,
        false,
        0.20,
        0.75,
    ),
    ("REDNESS", "Injection site erythema", true, true, 0.03, 0.12),
    (
        "SWELLING",
        "Injection site swelling",
        true,
        true,
        0.02,
        0.12,
    ),
    ("FATIGUE", "Fatigue", false, false, 0.20, 0.42),
    ("HEADACHE", "Headache", false, false, 0.20, 0.38),
    (
        "NEW OR WORSENED MUSCLE PAIN",
        "Myalgia",
        false,
        false,
        0.08,
        0.28,
    ),
    (
        "NEW OR WORSENED JOINT PAIN",
        "Arthralgia",
        false,
        false,
        0.06,
        0.16,
    ),
    ("CHILLS", "Chills", false, false, 0.04, 0.14),
    ("VOMITING", "Vomiting", false, false, 0.01, 0.02),
    ("DIARRHEA", "Diarrhoea", false, false, 0.05, 0.06),
];

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, VaccineParams) {
    let arms = [
        ArmSpec {
            armcd: "PBO".into(),
            arm: "Placebo".into(),
            treatment: "PLACEBO".into(),
            dose: 0.5,
            ratio: c.ratio[0],
            change: 0.0,
            ae_share: 0.0,
        },
        ArmSpec {
            armcd: "VACA".into(),
            arm: "Vaccine A".into(),
            treatment: "VACCINE A".into(),
            dose: 0.5,
            ratio: c.ratio[1],
            change: 0.0,
            ae_share: 0.0,
        },
    ];
    let (mut design, mut course) = build(
        opts,
        &Common {
            indication: "vaccine",
            weeks: 1,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMVAC01",
            phase: "PHASE III TRIAL",
            control: "PLACEBO",
            title: None,
        },
        &arms,
    );
    let (vax, eot, fu) = days(c);
    design.schedule = Some(schedule(c));
    design.cutoff = design
        .first_randomization
        .add_days(design.enrollment_days + fu + 60);
    let pop = &mut design.population;
    pop.age_min = c.age_range[0];
    pop.age_max = c.age_range[1];
    pop.age_mean = (c.age_range[0] + c.age_range[1]) as f64 / 2.0;
    pop.age_sd = ((c.age_range[1] - c.age_range[0]) as f64 / 4.0).max(1.0);
    pop.sex = Categorical::new([("F", 0.5), ("M", 0.5)]);
    for arm in &mut design.arms {
        let t = &mut arm.treatment;
        t.unit = "mL".into();
        t.form = "INJECTION".into();
        t.route = "INTRAMUSCULAR".into();
        t.frequency = "ONCE".into();
        t.cycle_days = c.interval_days as i32;
    }
    let info = &mut design.info;
    info.trial_types = vec![
        "IMMUNOGENICITY".into(),
        "REACTOGENICITY".into(),
        "SAFETY".into(),
    ];
    info.healthy_subjects = true;
    info.intervention_type = "VACCINE".into();
    course.death_annual_rate = 0.001;
    let solicited = SOLICITED
        .iter()
        .map(|&(term, pt, site, measured, placebo, vaccine)| Solicited {
            term,
            pt,
            site,
            measured,
            rate: vec![placebo, (vaccine * c.reactogenicity).min(0.95)],
        })
        .collect();
    let mut assays = vax.clone();
    assays.extend([eot, fu]);
    let p = VaccineParams {
        vaccination_days: vax,
        assay_days: assays,
        dose: 0.5,
        solicited,
        fever_rate: vec![0.005, c.fever_rate],
        later_dose_multiplier: vec![1.0, 1.3],
        seropositive: c.seropositive,
        gmt_by_arm: vec![10.0, c.gmt],
        responds: vec![false, true],
        first_dose_fraction: c.first_dose_fraction,
        log10_sd: c.titer_sd,
        waning_per_year: 0.4,
        limits: (10.0, 20480.0),
    };
    (design, course, p)
}
