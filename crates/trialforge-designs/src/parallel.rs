//! `parallel_continuous`: a randomised, double-blind, placebo-controlled parallel-group study
//! with a continuous primary endpoint measured at scheduled visits.
//!
//! Indication presets set the endpoint, the population and the safety profile:
//! - `type_2_diabetes` (default): HbA1c (%), screening HbA1c 7.0–10.5%, background
//!   metformin, GI-heavy AE profile for the active drug (GLP-1-like).
//! - `hypertension`: sitting systolic blood pressure (mmHg), screening SBP 140–180 mmHg,
//!   dizziness/oedema signal for the active drug.
//! - `obesity`: body weight (kg), BMI around 36, GI-heavy AE profile.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use trialforge_core::Date;
use trialforge_model::{
    Arm, Categorical, Design, Epoch, PlannedVisit, Population, Schedule, Site, Timing, Treatment,
    TrialInfo, VisitKind,
};
use trialforge_processes::fixed::{AnalyteOverride, CourseParams, EndpointEffect};

use crate::{ControlDoc, StudyOptions};

/// Controls for `parallel_continuous`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// `type_2_diabetes`, `hypertension` or `obesity`.
    pub indication: String,
    /// Treatment duration in weeks (visits every 4 weeks after Week 4).
    pub weeks: u32,
    /// Mean change from baseline at the end of treatment: [placebo, active]. `null` uses the
    /// indication's default.
    pub change: Option<[f64; 2]>,
    /// Allocation ratio: [placebo, active].
    pub ratio: [u32; 2],
    /// Active dose (mg, once daily).
    pub dose: f64,
    /// Annual dropout hazard (reasons other than AEs and death).
    pub dropout_rate: f64,
    /// Multiplier on every background AE rate.
    pub ae_rate_multiplier: f64,
    /// Probability an optional visit is missed.
    pub missed_visit: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            indication: "type_2_diabetes".into(),
            weeks: 24,
            change: None,
            ratio: [1, 1],
            dose: 10.0,
            dropout_rate: 0.15,
            ae_rate_multiplier: 1.5,
            missed_visit: 0.04,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "indication",
            "Indication preset: type_2_diabetes (HbA1c), hypertension (SYSBP) or obesity (WEIGHT)",
        ),
        ControlDoc::new(
            "weeks",
            "Treatment duration in weeks; visits at Weeks 2, 4 and every 4 weeks after",
        ),
        ControlDoc::new(
            "change",
            "Mean change from baseline at end of treatment: [placebo, active]; null for the preset",
        ),
        ControlDoc::new("ratio", "Allocation ratio: [placebo, active]"),
        ControlDoc::new("dose", "Active dose in mg, once daily"),
        ControlDoc::new(
            "dropout_rate",
            "Annual dropout hazard (reasons other than adverse events and death)",
        ),
        ControlDoc::new(
            "ae_rate_multiplier",
            "Multiplier on every background adverse event rate",
        ),
        ControlDoc::new(
            "missed_visit",
            "Probability an optional on-treatment visit is missed",
        ),
    ]
}

/// What an indication preset fixes.
struct Preset {
    endpoint: &'static str,
    change: [f64; 2],
    onset_days: f64,
    title: &'static str,
    disease: &'static str,
    area: &'static str,
    objective: &'static str,
    age: (f64, f64),
    bmi: (f64, f64),
    required_history: &'static [&'static str],
    background: &'static [(&'static str, &'static str, f64)],
    overrides: Vec<AnalyteOverride>,
    ae_risk: &'static [(&'static str, f64)],
    criteria: &'static [(&'static str, &'static str, &'static str)],
}

fn preset(indication: &str) -> Option<Preset> {
    let ov = |testcd: &str,
              mean: Option<(f64, f64)>,
              sd: Option<f64>,
              inclusion: Option<(f64, f64)>| AnalyteOverride {
        testcd: testcd.into(),
        mean,
        sd_between: sd,
        inclusion,
    };
    const GI: &[(&str, f64)] = &[
        ("Nausea", 3.0),
        ("Vomiting", 2.5),
        ("Diarrhoea", 2.0),
        ("Decreased appetite", 2.5),
        ("Constipation", 1.8),
        ("Dyspepsia", 1.8),
        ("Headache", 1.2),
    ];
    match indication {
        "type_2_diabetes" => Some(Preset {
            endpoint: "HBA1C",
            change: [-0.2, -1.0],
            onset_days: 35.0,
            title: "A Randomized, Double-Blind, Placebo-Controlled Study of Drug A in Adults with Type 2 Diabetes Mellitus",
            disease: "Type 2 diabetes mellitus",
            area: "ENDOCRINOLOGY",
            objective: "To compare the change from baseline in HbA1c at the end of treatment between Drug A and placebo",
            age: (57.0, 10.0),
            bmi: (32.0, 5.0),
            required_history: &["Type 2 diabetes mellitus"],
            background: &[("Type 2 diabetes mellitus", "METFORMIN", 0.85)],
            overrides: vec![
                ov("HBA1C", Some((8.2, 8.2)), Some(0.8), Some((7.0, 10.5))),
                ov("GLUC", Some((9.5, 9.6)), Some(0.25), None),
            ],
            ae_risk: &[
                ("Nausea", 3.0),
                ("Vomiting", 2.5),
                ("Diarrhoea", 2.0),
                ("Decreased appetite", 2.5),
                ("Constipation", 1.8),
                ("Dyspepsia", 1.8),
                ("Headache", 1.2),
                ("Hypoglycaemia", 3.0),
            ],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 18 to 80 years with type 2 diabetes mellitus",
                ),
                ("INCL02", "INCLUSION", "HbA1c of 7.0% to 10.5% at screening"),
                (
                    "INCL03",
                    "INCLUSION",
                    "Stable metformin dose, or no glucose-lowering therapy, for at least 90 days",
                ),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Type 1 diabetes mellitus or a history of diabetic ketoacidosis",
                ),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "eGFR below 30 mL/min/1.73 m2 at screening",
                ),
            ],
        }),
        "hypertension" => Some(Preset {
            endpoint: "SYSBP",
            change: [-4.0, -12.0],
            onset_days: 14.0,
            title: "A Randomized, Double-Blind, Placebo-Controlled Study of Drug A in Adults with Essential Hypertension",
            disease: "Essential hypertension",
            area: "CARDIOVASCULAR",
            objective: "To compare the change from baseline in sitting systolic blood pressure at the end of treatment between Drug A and placebo",
            age: (58.0, 10.0),
            bmi: (30.0, 5.0),
            required_history: &["Hypertension"],
            background: &[("Hypertension", "LISINOPRIL", 0.0)],
            overrides: vec![
                ov(
                    "SYSBP",
                    Some((152.0, 154.0)),
                    Some(9.0),
                    Some((140.0, 180.0)),
                ),
                ov("DIABP", Some((91.0, 93.0)), Some(7.0), None),
            ],
            ae_risk: &[
                ("Dizziness", 1.8),
                ("Headache", 1.3),
                ("Oedema peripheral", 2.0),
                ("Fatigue", 1.3),
            ],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 18 to 80 years with essential hypertension",
                ),
                (
                    "INCL02",
                    "INCLUSION",
                    "Mean sitting systolic blood pressure of 140 to 180 mmHg at screening",
                ),
                ("EXCL01", "EXCLUSION", "Secondary hypertension"),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "Myocardial infarction or stroke within 6 months before screening",
                ),
            ],
        }),
        "obesity" => Some(Preset {
            endpoint: "WEIGHT",
            change: [-1.5, -8.0],
            onset_days: 60.0,
            title: "A Randomized, Double-Blind, Placebo-Controlled Study of Drug A for Weight Management in Adults with Obesity",
            disease: "Obesity",
            area: "ENDOCRINOLOGY",
            objective: "To compare the change from baseline in body weight at the end of treatment between Drug A and placebo",
            age: (48.0, 11.0),
            bmi: (36.0, 4.0),
            required_history: &["Obesity"],
            background: &[],
            overrides: vec![],
            ae_risk: GI,
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 18 to 75 years with a BMI of at least 30 kg/m2",
                ),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Weight change of more than 5 kg in the 90 days before screening",
                ),
            ],
        }),
        _ => None,
    }
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    if preset(&c.indication).is_none() {
        return Err(format!(
            "unknown indication {:?}; available: type_2_diabetes, hypertension, obesity",
            c.indication
        ));
    }
    if !(4..=104).contains(&c.weeks) {
        return Err("`weeks` must be between 4 and 104".into());
    }
    if c.ratio.contains(&0) || c.ratio.iter().any(|&r| r > 10) {
        return Err("`ratio` entries must be between 1 and 10".into());
    }
    if !(c.dose.is_finite() && c.dose > 0.0) {
        return Err("`dose` must be positive".into());
    }
    if let Some(ch) = c.change {
        if !ch.iter().all(|x| x.is_finite()) {
            return Err("`change` must be finite".into());
        }
    }
    let p = |x: f64| x.is_finite() && (0.0..=1.0).contains(&x);
    if !(c.dropout_rate.is_finite() && c.dropout_rate >= 0.0 && c.dropout_rate <= 5.0) {
        return Err("`dropout_rate` must be between 0 and 5 per year".into());
    }
    if !(c.ae_rate_multiplier.is_finite()
        && c.ae_rate_multiplier >= 0.0
        && c.ae_rate_multiplier <= 20.0)
    {
        return Err("`ae_rate_multiplier` must be between 0 and 20".into());
    }
    if !p(c.missed_visit) || c.missed_visit > 0.5 {
        return Err("`missed_visit` must be between 0 and 0.5".into());
    }
    Ok(())
}

/// The visit schedule for `weeks` of treatment plus a 2-week safety follow-up.
fn schedule(weeks: u32) -> Schedule {
    let mut visits = Vec::new();
    let visit = |num: i64,
                 name: String,
                 day: Option<i32>,
                 window: (i32, i32),
                 epoch,
                 kind,
                 labs,
                 vitals| PlannedVisit {
        visitnum: num,
        name,
        day,
        window,
        epoch,
        kind,
        labs,
        vitals,
        endpoint: true,
    };
    visits.push(visit(
        1,
        "SCREENING".into(),
        Some(-14),
        (14, 0),
        Epoch::Screening,
        VisitKind::Screening,
        true,
        true,
    ));
    visits.push(visit(
        2,
        "BASELINE".into(),
        Some(1),
        (0, 0),
        Epoch::Treatment,
        VisitKind::Baseline,
        true,
        true,
    ));
    let mut wks: Vec<u32> = vec![2, 4];
    wks.extend((8..=weeks).step_by(4));
    if *wks.last().expect("non-empty") != weeks {
        wks.push(weeks);
    }
    wks.retain(|&w| w <= weeks);
    for (k, w) in wks.iter().enumerate() {
        let last = *w == weeks;
        let window = if *w <= 4 { (3, 3) } else { (7, 7) };
        visits.push(visit(
            3 + k as i64,
            format!("WEEK {w}"),
            Some((w * 7 + 1) as i32),
            window,
            Epoch::Treatment,
            if last {
                VisitKind::EndOfTreatment
            } else {
                VisitKind::Treatment
            },
            last || *w == 4 || *w == weeks / 2,
            true,
        ));
    }
    let fu_day = (weeks * 7 + 1 + 14) as i32;
    visits.push(visit(
        90,
        "FOLLOW-UP".into(),
        Some(fu_day),
        (3, 3),
        Epoch::FollowUp,
        VisitKind::FollowUp,
        true,
        true,
    ));
    visits.push(visit(
        99,
        "EARLY TERMINATION".into(),
        None,
        (0, 0),
        Epoch::Treatment,
        VisitKind::EarlyTermination,
        true,
        true,
    ));
    Schedule {
        visits,
        treatment_days: (weeks * 7) as i32,
        followup_days: 14,
    }
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams) {
    let pre = preset(&c.indication).expect("validated");
    let change = c.change.unwrap_or(pre.change);
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
        ("USA", 12),
        ("CAN", 3),
        ("MEX", 3),
        ("POL", 4),
        ("DEU", 3),
        ("ESP", 3),
        ("JPN", 2),
    ];
    let sites = countries
        .iter()
        .flat_map(|&(c, k)| std::iter::repeat_n(c, k))
        .enumerate()
        .map(|(i, c)| Site {
            siteid: (201 + i).to_string(),
            country: c.into(),
        })
        .collect();
    let first = Date::from_ymd(2024, 2, 1).expect("valid date");
    let enrollment_days = 365;
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| "SIMPAR01".into()),
        seed: opts.seed,
        n: opts.n,
        arms: vec![
            Arm {
                armcd: "PBO".into(),
                arm: "Placebo".into(),
                ratio: c.ratio[0],
                treatment: treatment("PLACEBO", 0.0),
            },
            Arm {
                armcd: "DRUGA".into(),
                arm: format!("Drug A {} mg", c.dose),
                ratio: c.ratio[1],
                treatment: treatment("DRUG A", c.dose),
            },
        ],
        sites,
        first_randomization: first,
        enrollment_days,
        cutoff: first.add_days(enrollment_days + (c.weeks * 7) as i32 + 120),
        block_size: (c.ratio[0] + c.ratio[1]) * 2,
        population: Population {
            age_mean: pre.age.0,
            age_sd: pre.age.1,
            age_min: 18,
            age_max: 80,
            sex: Categorical::new([("F", 0.48), ("M", 0.52)]),
            race: Categorical::new([
                ("WHITE", 0.70),
                ("BLACK OR AFRICAN AMERICAN", 0.12),
                ("ASIAN", 0.12),
                ("AMERICAN INDIAN OR ALASKA NATIVE", 0.03),
                ("MULTIPLE", 0.03),
            ]),
            ethnic: Categorical::new([
                ("NOT HISPANIC OR LATINO", 0.82),
                ("HISPANIC OR LATINO", 0.18),
            ]),
            strata: Categorical::new([("ALL", 1.0)]),
        },
        timing: Timing {
            screening_lead: (7, 21),
            consent_lead: (0, 7),
            dose_lag: (0, 0),
        },
        schedule: Some(schedule(c.weeks)),
        info: TrialInfo {
            title: pre.title.into(),
            phase: "PHASE III TRIAL".into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: pre.disease.into(),
            therapeutic_area: pre.area.into(),
            primary_objective: pre.objective.into(),
            blinding: "DOUBLE BLIND".into(),
            control: "PLACEBO".into(),
            add_on: !pre.background.is_empty() && pre.background.iter().any(|b| b.2 > 0.0),
            criteria: pre
                .criteria
                .iter()
                .map(|(a, b, c)| (a.to_string(), b.to_string(), c.to_string()))
                .collect(),
        },
    };
    let active_risk: HashMap<String, f64> = pre
        .ae_risk
        .iter()
        .map(|(pt, rr)| (pt.to_string(), *rr))
        .collect();
    let params = CourseParams {
        dropout_annual_rate: c.dropout_rate,
        dropout_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.55),
            ("LOST TO FOLLOW-UP".into(), 0.2),
            ("PHYSICIAN DECISION".into(), 0.1),
            ("LACK OF EFFICACY".into(), 0.1),
            ("PROTOCOL DEVIATION".into(), 0.05),
        ],
        death_annual_rate: 0.004,
        missed_visit_prob: c.missed_visit,
        ae_rate_multiplier: c.ae_rate_multiplier,
        ae_relative_risk: vec![HashMap::new(), active_risk],
        ae_discontinuation: [0.004, 0.04, 0.25],
        ae_interruption: [0.0, 0.05, 0.3],
        required_history: pre.required_history.iter().map(|s| s.to_string()).collect(),
        background_medication: pre
            .background
            .iter()
            .map(|(a, b, p)| (a.to_string(), b.to_string(), *p))
            .collect(),
        analyte_overrides: pre.overrides,
        endpoints: vec![EndpointEffect {
            testcd: pre.endpoint.into(),
            change_by_arm: change.to_vec(),
            onset_days: pre.onset_days,
        }],
        bmi: pre.bmi,
    };
    (design, params)
}
