//! `cv_outcomes`: an event-driven, randomised, double-blind, placebo-controlled cardiovascular
//! outcomes trial. The primary endpoint is time to first MACE (cardiovascular death, non-fatal
//! myocardial infarction or non-fatal ischaemic stroke); non-cardiovascular death competes.
//! MI and stroke can recur. The study ends when a target number of first MACE events occur.

use serde::{Deserialize, Serialize};
use trialforge_core::Date;
use trialforge_model::{
    Arm, Categorical, Design, Epoch, PlannedVisit, Population, Schedule, Site, Timing, Treatment,
    TrialInfo, VisitKind,
};
use trialforge_processes::outcomes::{EventType, OutcomesParams};

use crate::{ControlDoc, StudyOptions};

/// Controls for `cv_outcomes`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// First MACE events that end the study.
    pub target_events: usize,
    /// Annual rate of first MACE on placebo (for a subject without prior MI).
    pub mace_rate: f64,
    /// Hazard ratio of each MACE component, Drug A vs placebo.
    pub hr: f64,
    /// Shares of MACE among [CV death, MI, stroke].
    pub components: [f64; 3],
    /// Annual rate of non-cardiovascular death (competing risk, unaffected by treatment).
    pub non_cv_death_rate: f64,
    /// Hazard multiplier for subjects with a prior MI (the stratification factor).
    pub prior_mi_hr: f64,
    /// SD of the log frailty shared by a subject's events.
    pub frailty_sd: f64,
    /// Enrollment period in months.
    pub enrollment_months: f64,
    /// Longest possible study duration in months (the cutoff if the target is not reached).
    pub max_months: f64,
    /// Annual hazard of permanently stopping study treatment.
    pub discontinuation_rate: f64,
    /// Annual hazard of withdrawing from follow-up.
    pub withdrawal_rate: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            target_events: 150,
            mace_rate: 0.045,
            hr: 0.8,
            components: [0.3, 0.45, 0.25],
            non_cv_death_rate: 0.01,
            prior_mi_hr: 1.5,
            frailty_sd: 0.5,
            enrollment_months: 18.0,
            max_months: 60.0,
            discontinuation_rate: 0.06,
            withdrawal_rate: 0.02,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "target_events",
            "First MACE events that end the study (event-driven cutoff)",
        ),
        ControlDoc::new(
            "mace_rate",
            "Annual rate of first MACE on placebo (no prior MI)",
        ),
        ControlDoc::new(
            "hr",
            "Hazard ratio of each MACE component, Drug A vs placebo",
        ),
        ControlDoc::new("components", "Shares of MACE among [CV death, MI, stroke]"),
        ControlDoc::new(
            "non_cv_death_rate",
            "Annual rate of non-cardiovascular death (competing risk, unaffected by treatment)",
        ),
        ControlDoc::new(
            "prior_mi_hr",
            "Hazard multiplier for subjects with a prior MI (stratum)",
        ),
        ControlDoc::new(
            "frailty_sd",
            "SD of the log frailty shared by a subject's events",
        ),
        ControlDoc::new("enrollment_months", "Enrollment period in months"),
        ControlDoc::new(
            "max_months",
            "Longest study duration in months, the cutoff if the target is not reached",
        ),
        ControlDoc::new(
            "discontinuation_rate",
            "Annual hazard of permanently stopping study treatment",
        ),
        ControlDoc::new(
            "withdrawal_rate",
            "Annual hazard of withdrawing from follow-up",
        ),
    ]
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    let pos = |x: f64| x.is_finite() && x > 0.0;
    if !(1..=100_000).contains(&c.target_events) {
        return Err("`target_events` must be between 1 and 100000".into());
    }
    for (name, x) in [
        ("mace_rate", c.mace_rate),
        ("hr", c.hr),
        ("prior_mi_hr", c.prior_mi_hr),
        ("enrollment_months", c.enrollment_months),
    ] {
        if !pos(x) {
            return Err(format!("`{name}` must be positive"));
        }
    }
    if c.mace_rate > 2.0 || c.hr > 5.0 || c.prior_mi_hr > 10.0 {
        return Err(
            "`mace_rate` must be at most 2, `hr` at most 5 and `prior_mi_hr` at most 10".into(),
        );
    }
    if !c.components.iter().all(|x| x.is_finite() && *x >= 0.0)
        || c.components.iter().sum::<f64>() <= 0.0
    {
        return Err("`components` must be non-negative with a positive sum".into());
    }
    for (name, x) in [
        ("non_cv_death_rate", c.non_cv_death_rate),
        ("discontinuation_rate", c.discontinuation_rate),
        ("withdrawal_rate", c.withdrawal_rate),
    ] {
        if !(0.0..=2.0).contains(&x) {
            return Err(format!("`{name}` must be between 0 and 2 per year"));
        }
    }
    if !(0.0..=2.0).contains(&c.frailty_sd) {
        return Err("`frailty_sd` must be between 0 and 2".into());
    }
    if c.max_months <= c.enrollment_months || c.max_months > 180.0 {
        return Err("`max_months` must exceed `enrollment_months` and be at most 180".into());
    }
    Ok(())
}

fn month(m: f64) -> i32 {
    (m * 30.4375).round() as i32
}

/// Visits: screening, Day 1, then every 4 months up to the longest possible follow-up.
fn schedule(c: &Controls) -> Schedule {
    let visit = |num: i64, name: String, day: i32, kind| PlannedVisit {
        visitnum: num,
        name,
        day: Some(day),
        window: (0, 0),
        epoch: if kind == VisitKind::Screening {
            Epoch::Screening
        } else {
            Epoch::Treatment
        },
        kind,
        labs: false,
        vitals: false,
        endpoint: true,
    };
    let mut visits = vec![
        visit(1, "SCREENING".into(), -14, VisitKind::Screening),
        visit(2, "DAY 1".into(), 1, VisitKind::Baseline),
    ];
    let mut k = 1;
    while k * 120 < month(c.max_months) {
        visits.push(visit(
            2 + i64::from(k),
            format!("MONTH {}", 4 * k),
            1 + k * 120,
            VisitKind::Treatment,
        ));
        k += 1;
    }
    Schedule {
        visits,
        treatment_days: month(c.max_months),
        followup_days: 0,
    }
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, OutcomesParams) {
    let treatment = |name: &str, dose: f64| Treatment {
        name: name.into(),
        dose,
        unit: "mg".into(),
        form: "TABLET".into(),
        route: "ORAL".into(),
        frequency: "QD".into(),
        cycle_days: 1,
    };
    let countries = [
        ("USA", 10),
        ("CAN", 3),
        ("GBR", 3),
        ("DEU", 4),
        ("POL", 4),
        ("BRA", 3),
        ("CHN", 3),
    ];
    let sites = countries
        .iter()
        .flat_map(|&(c, k)| std::iter::repeat_n(c, k))
        .enumerate()
        .map(|(i, c)| Site {
            siteid: (301 + i).to_string(),
            country: c.into(),
        })
        .collect();
    let first = Date::from_ymd(2023, 6, 1).expect("valid date");
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| "SIMCVO01".into()),
        seed: opts.seed,
        n: opts.n,
        arms: vec![
            Arm { armcd: "DRUGA".into(), arm: "Drug A".into(), ratio: 1, treatment: treatment("DRUG A", 20.0), periods: Vec::new() },
            Arm { armcd: "PBO".into(), arm: "Placebo".into(), ratio: 1, treatment: treatment("PLACEBO", 0.0), periods: Vec::new() },
        ],
        sites,
        first_randomization: first,
        enrollment_days: month(c.enrollment_months),
        cutoff: first.add_days(month(c.max_months)),
        block_size: 4,
        population: Population {
            age_mean: 66.0,
            age_sd: 8.0,
            age_min: 45,
            age_max: 90,
            sex: Categorical::new([("M", 0.7), ("F", 0.3)]),
            race: Categorical::new([
                ("WHITE", 0.75),
                ("BLACK OR AFRICAN AMERICAN", 0.08),
                ("ASIAN", 0.14),
                ("MULTIPLE", 0.03),
            ]),
            ethnic: Categorical::new([("NOT HISPANIC OR LATINO", 0.88), ("HISPANIC OR LATINO", 0.12)]),
            strata: Categorical::new([("NO PRIOR MI", 0.6), ("PRIOR MI", 0.4)]),
        },
        timing: Timing { screening_lead: (3, 21), consent_lead: (0, 7), dose_lag: (0, 1) },
        schedule: Some(schedule(c)),
        cohorts: Vec::new(),
        info: TrialInfo {
            title: "A Randomized, Double-Blind, Placebo-Controlled, Event-Driven Cardiovascular Outcomes Study of Drug A".into(),
            phase: "PHASE III TRIAL".into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: "Atherosclerotic cardiovascular disease".into(),
            therapeutic_area: "CARDIOVASCULAR".into(),
            primary_objective: "To compare the time to first major adverse cardiovascular event (CV death, non-fatal MI or non-fatal ischaemic stroke) between Drug A and placebo".into(),
            blinding: "DOUBLE BLIND".into(),
            control: "PLACEBO".into(),
            add_on: true,
            trial_types: vec!["EFFICACY".into(), "SAFETY".into()],
            intervention_model: "PARALLEL".into(),
            healthy_subjects: false,
            criteria: vec![
                ("INCL01".into(), "INCLUSION".into(), "Adults aged 45 or older with established atherosclerotic cardiovascular disease".into()),
                ("EXCL01".into(), "EXCLUSION".into(), "Myocardial infarction or stroke within 30 days before randomization".into()),
            ],
        },
    };
    let share = c.components.iter().sum::<f64>();
    let rate = |k: usize| c.mace_rate * c.components[k] / share;
    let event = |pt: &str, rate: f64, hr: f64, terminal: bool, mace: bool| EventType {
        pt: pt.into(),
        rate,
        hr: vec![hr, 1.0],
        terminal,
        mace,
    };
    let params = OutcomesParams {
        events: vec![
            event("Cardiac death", rate(0), c.hr, true, true),
            event("Myocardial infarction", rate(1), c.hr, false, true),
            event("Ischaemic stroke", rate(2), c.hr, false, true),
            event("Death", c.non_cv_death_rate, 1.0, true, false),
        ],
        stratum_hr: vec![1.0, c.prior_mi_hr],
        frailty_sd: c.frailty_sd,
        withdrawal_rate: c.withdrawal_rate,
        withdrawal_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.7),
            ("LOST TO FOLLOW-UP".into(), 0.3),
        ],
        discontinuation_rate: c.discontinuation_rate,
        visit_interval: 120,
        target_events: c.target_events,
        max_days: month(c.max_months),
    };
    (design, params)
}
