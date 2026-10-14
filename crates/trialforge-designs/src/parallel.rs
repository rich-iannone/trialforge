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
    /// `placebo` or `active` (an active comparator, for non-inferiority).
    pub comparator: String,
    /// The active comparator's dose (mg, once daily).
    pub comparator_dose: f64,
    /// Non-inferiority margin for Drug A minus the comparator (endpoint units; lower is
    /// better), or `null`.
    pub ni_margin: Option<f64>,
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
            comparator: "placebo".into(),
            comparator_dose: 10.0,
            ni_margin: None,
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
        ControlDoc::new(
            "comparator",
            "placebo, or active (an active comparator: change defaults to equal effects)",
        ),
        ControlDoc::new(
            "comparator_dose",
            "Active comparator dose in mg, once daily",
        ),
        ControlDoc::new(
            "ni_margin",
            "Non-inferiority margin for Drug A minus comparator (lower is better); null for none",
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
        // For cns_alzheimers: the endpoint is ADAS-Cog (FT), not a lab or vital sign.
        "alzheimers" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 28.0,
            title: "A Randomized, Double-Blind, Placebo-Controlled Study of Drug A in Mild to Moderate Alzheimer's Disease",
            disease: "Alzheimer's disease",
            area: "NEUROLOGY",
            objective: "To compare the change from baseline in ADAS-Cog (11) at the end of treatment between Drug A and placebo",
            age: (75.0, 8.0),
            bmi: (25.0, 4.0),
            required_history: &["Dementia Alzheimer's type"],
            background: &[("Dementia Alzheimer's type", "DONEPEZIL", 0.0)],
            overrides: vec![],
            ae_risk: &[
                ("Nausea", 2.5),
                ("Vomiting", 2.5),
                ("Diarrhoea", 2.0),
                ("Dizziness", 1.8),
                ("Decreased appetite", 1.8),
                ("Syncope", 2.0),
                ("Headache", 1.2),
            ],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 50 to 89 years with probable Alzheimer's disease (NINCDS-ADRDA)",
                ),
                (
                    "INCL02",
                    "INCLUSION",
                    "Mini-Mental State Examination score of 10 to 23 at screening",
                ),
                (
                    "INCL03",
                    "INCLUSION",
                    "A caregiver who sees the participant at least 3 days a week",
                ),
                ("EXCL01", "EXCLUSION", "Any other cause of dementia"),
            ],
        }),
        // For ophthalmology: the endpoint is visual acuity (OE), not a lab or vital sign.
        "neovascular_amd" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 21.0,
            title: "A Randomized, Double-Masked, Active-Controlled Study of Intravitreal Drug A in Neovascular Age-Related Macular Degeneration",
            disease: "Neovascular age-related macular degeneration",
            area: "OPHTHALMOLOGY",
            objective: "To show that intravitreal Drug A is non-inferior to the comparator in the change from baseline in best corrected visual acuity of the study eye at the end of treatment",
            age: (77.0, 8.0),
            bmi: (27.0, 4.5),
            required_history: &["Neovascular age-related macular degeneration"],
            background: &[],
            overrides: vec![],
            ae_risk: &[],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 50 or older with active choroidal neovascularisation secondary to age-related macular degeneration in the study eye",
                ),
                (
                    "INCL02",
                    "INCLUSION",
                    "Best corrected visual acuity of 78 to 24 ETDRS letters in the study eye at screening",
                ),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Prior anti-VEGF treatment of the study eye",
                ),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "Active intraocular or periocular infection in either eye",
                ),
            ],
        }),
        // For pediatric_growth: the endpoint is height velocity, from the growth layer.
        "growth_hormone_deficiency" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 30.0,
            title: "A Randomized, Open-Label, Active-Controlled Study of Once-Weekly Drug A Versus Daily Somatropin in Prepubertal Children with Growth Hormone Deficiency",
            disease: "Growth hormone deficiency",
            area: "ENDOCRINOLOGY",
            objective: "To show that once-weekly Drug A is non-inferior to daily somatropin in annualized height velocity at the end of treatment",
            age: (7.0, 2.3),
            bmi: (16.5, 2.0),
            required_history: &["Growth hormone deficiency"],
            background: &[],
            overrides: vec![
                ov("SYSBP", Some((100.0, 102.0)), Some(8.0), None),
                ov("DIABP", Some((62.0, 63.0)), Some(6.0), None),
                ov("PULSE", Some((92.0, 90.0)), Some(10.0), None),
                ov("RESP", Some((20.0, 20.0)), Some(2.5), None),
            ],
            ae_risk: &[],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Prepubertal children aged 3 to 11 years with growth hormone deficiency confirmed by two stimulation tests",
                ),
                (
                    "INCL02",
                    "INCLUSION",
                    "Height at least 2.0 standard deviations below the mean for age and sex (CDC 2000) at screening",
                ),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Prior treatment with growth hormone or IGF-1",
                ),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "Chromosomal abnormalities or syndromes associated with short stature",
                ),
            ],
        }),
        // For respiratory_exacerbation: exacerbations (CE) and FEV1 (RE) come from their own layer.
        "copd" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 14.0,
            title: "A Randomized, Double-Blind, Placebo-Controlled Study of Drug A in Adults with Chronic Obstructive Pulmonary Disease and a History of Exacerbations",
            disease: "Chronic obstructive pulmonary disease",
            area: "RESPIRATORY",
            objective: "To compare the annualized rate of moderate or severe COPD exacerbations between Drug A and placebo",
            age: (65.0, 8.0),
            bmi: (26.0, 5.0),
            required_history: &["Chronic obstructive pulmonary disease"],
            background: &[("Chronic obstructive pulmonary disease", "TIOTROPIUM", 0.9)],
            overrides: vec![],
            ae_risk: &[],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults aged 40 to 85 years with COPD and a post-bronchodilator FEV1/FVC below 0.70",
                ),
                (
                    "INCL02",
                    "INCLUSION",
                    "At least one moderate or severe COPD exacerbation in the 12 months before screening",
                ),
                ("EXCL01", "EXCLUSION", "A current diagnosis of asthma"),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "A COPD exacerbation within 4 weeks before screening",
                ),
            ],
        }),
        // For vaccine: healthy adults; immunogenicity and reactogenicity come from their layer.
        "vaccine" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 14.0,
            title: "A Randomized, Observer-Blind, Placebo-Controlled Study of the Immunogenicity, Reactogenicity and Safety of Vaccine A in Healthy Adults",
            disease: "Prevention of Virus X infection",
            area: "VACCINES",
            objective: "To compare the neutralizing antibody geometric mean titer 28 days after the last vaccination between Vaccine A and placebo",
            age: (40.0, 12.0),
            bmi: (26.0, 4.5),
            required_history: &[],
            background: &[],
            overrides: vec![],
            ae_risk: &[],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Healthy adults, or adults with stable chronic conditions, within the protocol age range",
                ),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Receipt of any vaccine within 28 days before the first vaccination",
                ),
                (
                    "EXCL02",
                    "EXCLUSION",
                    "Known immunodeficiency or immunosuppressive therapy",
                ),
            ],
        }),
        // For dose_escalation: phase 1 in advanced solid tumours; DLTs come from their layer.
        "advanced_solid_tumors" => Some(Preset {
            endpoint: "",
            change: [0.0, 0.0],
            onset_days: 14.0,
            title: "An Open-Label, Dose-Escalation Study of the Safety and Tolerability of Drug A in Adults with Advanced Solid Tumors",
            disease: "Advanced solid tumors",
            area: "ONCOLOGY",
            objective: "To determine the maximum tolerated dose of Drug A from dose-limiting toxicities in the first 28-day cycle",
            age: (60.0, 10.0),
            bmi: (26.0, 4.5),
            required_history: &["Neoplasm malignant"],
            background: &[],
            overrides: vec![],
            ae_risk: &[],
            criteria: &[
                (
                    "INCL01",
                    "INCLUSION",
                    "Adults with an advanced or metastatic solid tumor for which no standard therapy is available",
                ),
                ("INCL02", "INCLUSION", "ECOG performance status of 0 or 1"),
                (
                    "EXCL01",
                    "EXCLUSION",
                    "Anticancer therapy within 28 days before the first dose",
                ),
            ],
        }),
        _ => None,
    }
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    if preset(&c.indication).is_none()
        || [
            "alzheimers",
            "neovascular_amd",
            "growth_hormone_deficiency",
            "copd",
            "vaccine",
            "advanced_solid_tumors",
        ]
        .contains(&c.indication.as_str())
    {
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
    if c.comparator != "placebo" && c.comparator != "active" {
        return Err("`comparator` must be \"placebo\" or \"active\"".into());
    }
    if !(c.comparator_dose.is_finite() && c.comparator_dose > 0.0) {
        return Err("`comparator_dose` must be positive".into());
    }
    if c.ni_margin.is_some_and(|m| !(m.is_finite() && m > 0.0)) {
        return Err("`ni_margin` must be positive".into());
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

/// One arm of a parallel-group design.
pub(crate) struct ArmSpec {
    pub armcd: String,
    pub arm: String,
    /// `EXTRT`.
    pub treatment: String,
    pub dose: f64,
    pub ratio: u32,
    /// Mean change from baseline in the endpoint at the end of treatment.
    pub change: f64,
    /// Share of the indication's active-drug AE signal (0 = none, 1 = full).
    pub ae_share: f64,
}

/// What a parallel-group archetype fixes beyond its arms.
pub(crate) struct Common<'a> {
    pub indication: &'a str,
    pub weeks: u32,
    pub dropout_rate: f64,
    pub ae_rate_multiplier: f64,
    pub missed_visit: f64,
    pub studyid: &'a str,
    pub phase: &'a str,
    /// `TCNTRL`: `PLACEBO` or `ACTIVE`.
    pub control: &'a str,
    /// Replaces the preset's title and primary objective.
    pub title: Option<(String, String)>,
}

/// The indication preset's [placebo, active] change at the end of treatment.
pub(crate) fn preset_change(indication: &str) -> Option<[f64; 2]> {
    preset(indication).map(|p| p.change)
}

/// A parallel-group design and course with the given arms.
pub(crate) fn build(
    opts: &StudyOptions,
    c: &Common<'_>,
    arms: &[ArmSpec],
) -> (Design, CourseParams) {
    let pre = preset(c.indication).expect("validated");
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
    let (title, objective) = c
        .title
        .clone()
        .unwrap_or_else(|| (pre.title.to_string(), pre.objective.to_string()));
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| c.studyid.into()),
        seed: opts.seed,
        n: opts.n,
        arms: arms
            .iter()
            .map(|a| Arm {
                armcd: a.armcd.clone(),
                arm: a.arm.clone(),
                ratio: a.ratio,
                treatment: treatment(&a.treatment, a.dose),
                periods: Vec::new(),
                switch_to: None,
            })
            .collect(),
        sites,
        first_randomization: first,
        enrollment_days,
        cutoff: first.add_days(enrollment_days + (c.weeks * 7) as i32 + 120),
        interim: false,
        block_size: arms.iter().map(|a| a.ratio).sum::<u32>() * 2,
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
        cohorts: Vec::new(),
        info: TrialInfo {
            title,
            phase: c.phase.into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: pre.disease.into(),
            therapeutic_area: pre.area.into(),
            primary_objective: objective,
            blinding: "DOUBLE BLIND".into(),
            control: c.control.into(),
            add_on: !pre.background.is_empty() && pre.background.iter().any(|b| b.2 > 0.0),
            trial_types: vec!["EFFICACY".into()],
            intervention_model: "PARALLEL".into(),
            healthy_subjects: false,
            intervention_type: "DRUG".into(),
            criteria: pre
                .criteria
                .iter()
                .map(|(a, b, c)| (a.to_string(), b.to_string(), c.to_string()))
                .collect(),
        },
    };
    let ae_relative_risk: Vec<HashMap<String, f64>> = arms
        .iter()
        .map(|a| {
            if a.ae_share <= 0.0 {
                HashMap::new()
            } else {
                pre.ae_risk
                    .iter()
                    .map(|(pt, rr)| {
                        let r = if a.ae_share >= 1.0 {
                            *rr
                        } else {
                            1.0 + (rr - 1.0) * a.ae_share
                        };
                        (pt.to_string(), r)
                    })
                    .collect()
            }
        })
        .collect();
    let params = CourseParams {
        dosing_days: None,
        extra_ae: Vec::new(),
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
        ae_relative_risk,
        ae_discontinuation: [0.004, 0.04, 0.25],
        ae_interruption: [0.0, 0.05, 0.3],
        required_history: pre.required_history.iter().map(|s| s.to_string()).collect(),
        background_medication: pre
            .background
            .iter()
            .map(|(a, b, p)| (a.to_string(), b.to_string(), *p))
            .collect(),
        analyte_overrides: pre.overrides,
        endpoints: if pre.endpoint.is_empty() {
            Vec::new()
        } else {
            vec![EndpointEffect {
                testcd: pre.endpoint.into(),
                change_by_arm: arms.iter().map(|a| a.change).collect(),
                onset_days: pre.onset_days,
            }]
        },
        bmi: pre.bmi,
    };
    (design, params)
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams) {
    let pre = preset(&c.indication).expect("validated");
    let active = c.comparator == "active";
    // Against an active comparator the default is two equally effective drugs.
    let change = c.change.unwrap_or(if active {
        [pre.change[1]; 2]
    } else {
        pre.change
    });
    let control = if active {
        ArmSpec {
            armcd: "ACTCTRL".into(),
            arm: format!("Comparator {} mg", c.comparator_dose),
            treatment: "COMPARATOR".into(),
            dose: c.comparator_dose,
            ratio: c.ratio[0],
            change: change[0],
            ae_share: 1.0,
        }
    } else {
        ArmSpec {
            armcd: "PBO".into(),
            arm: "Placebo".into(),
            treatment: "PLACEBO".into(),
            dose: 0.0,
            ratio: c.ratio[0],
            change: change[0],
            ae_share: 0.0,
        }
    };
    let arms = [
        control,
        ArmSpec {
            armcd: "DRUGA".into(),
            arm: format!("Drug A {} mg", c.dose),
            treatment: "DRUG A".into(),
            dose: c.dose,
            ratio: c.ratio[1],
            change: change[1],
            ae_share: 1.0,
        },
    ];
    let title = active.then(|| {
        (
            format!("A Randomized, Double-Blind, Active-Controlled, Non-Inferiority Study of Drug A in Adults with {}", pre.disease),
            format!("To show that Drug A is non-inferior to the comparator in the change from baseline in {} at the end of treatment", pre.endpoint),
        )
    });
    build(
        opts,
        &Common {
            indication: &c.indication,
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMPAR01",
            phase: "PHASE III TRIAL",
            control: if active { "ACTIVE" } else { "PLACEBO" },
            title,
        },
        &arms,
    )
}
