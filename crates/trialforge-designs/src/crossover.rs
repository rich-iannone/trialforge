//! `crossover_be`: a bioequivalence (or food-effect) crossover study in healthy volunteers.
//! Each subject receives a single oral dose in every period, with intensive PK sampling over a
//! 3-day in-house stay and a washout between doses.
//!
//! - `2x2`: Test and Reference, sequences TR and RT.
//! - `2x2x4`: full replicate, sequences TRTR and RTRT.
//! - `food_effect`: the same product fed and fasted, sequences FED-FASTED and FASTED-FED.
//!
//! The Test (or fed) product differs from the Reference (or fasted) by a relative
//! bioavailability (`gmr`, the true geometric mean ratio of AUC) and an absorption-rate ratio
//! (`ka_ratio`, which moves Cmax and Tmax); clearance varies between periods within a subject
//! (`intra_cv`).

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use trialforge_core::Date;
use trialforge_model::{
    Arm, Categorical, Design, Epoch, Period, PlannedVisit, Population, Schedule, Site, Timing,
    Treatment, TrialInfo, VisitKind,
};
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::pk::PkParams;

use crate::phase1::DRUG;
use crate::{ControlDoc, StudyOptions};

/// Days after a dose before discharge (72 h of sampling).
pub const IN_HOUSE_DAYS: i32 = 3;

/// Controls for `crossover_be`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// `2x2`, `2x2x4` or `food_effect`.
    pub design: String,
    /// Dose (mg) in every period.
    pub dose: f64,
    /// True Test/Reference (or fed/fasted) geometric mean ratio of AUC.
    pub gmr: f64,
    /// Test/Reference (or fed/fasted) ratio of the absorption rate constant.
    pub ka_ratio: f64,
    /// Within-subject CV of exposure between periods.
    pub intra_cv: f64,
    /// Between-subject CV of clearance.
    pub inter_cv: f64,
    /// Days from one period's dose to the next (at least 5).
    pub days_between_doses: u32,
    /// Apparent clearance CL/F (L/h) of a 70 kg subject on the Reference product.
    pub cl: f64,
    /// Apparent volume V/F (L) of a 70 kg subject.
    pub v: f64,
    /// Absorption rate constant (1/h) of the Reference product.
    pub ka: f64,
    /// Lower limit of quantitation (ng/mL).
    pub lloq: f64,
    /// Annual dropout hazard.
    pub dropout_rate: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            design: "2x2".into(),
            dose: 100.0,
            gmr: 0.95,
            ka_ratio: 1.0,
            intra_cv: 0.2,
            inter_cv: 0.3,
            days_between_doses: 8,
            cl: 12.0,
            v: 90.0,
            ka: 1.2,
            lloq: 0.5,
            dropout_rate: 1.0,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "design",
            "2x2 (TR/RT), 2x2x4 (TRTR/RTRT) or food_effect (fed/fasted 2x2)",
        ),
        ControlDoc::new("dose", "Dose in mg given in every period"),
        ControlDoc::new(
            "gmr",
            "True Test/Reference (or fed/fasted) geometric mean ratio of AUC",
        ),
        ControlDoc::new(
            "ka_ratio",
            "Test/Reference (or fed/fasted) ratio of the absorption rate, which moves Cmax and Tmax",
        ),
        ControlDoc::new("intra_cv", "Within-subject CV of exposure between periods"),
        ControlDoc::new("inter_cv", "Between-subject CV of clearance"),
        ControlDoc::new(
            "days_between_doses",
            "Days from one period's dose to the next (at least 5)",
        ),
        ControlDoc::new(
            "cl",
            "Apparent clearance CL/F (L/h) of a 70 kg subject on the Reference",
        ),
        ControlDoc::new("v", "Apparent volume V/F (L) of a 70 kg subject"),
        ControlDoc::new("ka", "Absorption rate constant (1/h) of the Reference"),
        ControlDoc::new("lloq", "Lower limit of quantitation (ng/mL)"),
        ControlDoc::new("dropout_rate", "Annual dropout hazard"),
    ]
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if !matches!(c.design.as_str(), "2x2" | "2x2x4" | "food_effect") {
        return Err("`design` must be \"2x2\", \"2x2x4\" or \"food_effect\"".into());
    }
    for (name, x) in [
        ("dose", c.dose),
        ("cl", c.cl),
        ("v", c.v),
        ("ka", c.ka),
        ("lloq", c.lloq),
    ] {
        if !positive(x) {
            return Err(format!("`{name}` must be positive"));
        }
    }
    if !(0.2..=5.0).contains(&c.gmr) || !(0.1..=10.0).contains(&c.ka_ratio) {
        return Err("`gmr` must be in [0.2, 5] and `ka_ratio` in [0.1, 10]".into());
    }
    if !(0.0..=1.0).contains(&c.intra_cv) || !(0.0..=1.5).contains(&c.inter_cv) {
        return Err("`intra_cv` must be in [0, 1] and `inter_cv` in [0, 1.5]".into());
    }
    if !(5..=60).contains(&c.days_between_doses) {
        return Err("`days_between_doses` must be between 5 and 60".into());
    }
    if !(0.0..=10.0).contains(&c.dropout_rate) {
        return Err("`dropout_rate` must be between 0 and 10".into());
    }
    Ok(())
}

/// (sequence codes, treatment codes and labels) of a design.
fn layout(c: &Controls) -> (Vec<&'static str>, [(&'static str, &'static str); 2]) {
    match c.design.as_str() {
        "2x2x4" => (vec!["TRTR", "RTRT"], [("T", "test"), ("R", "reference")]),
        "food_effect" => (vec!["AB", "BA"], [("FED", "fed"), ("FAST", "fasted")]),
        _ => (vec!["TR", "RT"], [("T", "test"), ("R", "reference")]),
    }
}

/// Study day of each period's dose.
pub fn period_days(c: &Controls) -> Vec<i32> {
    let n = layout(c).0[0].len() as i32;
    (0..n)
        .map(|k| 1 + k * c.days_between_doses as i32)
        .collect()
}

fn schedule(c: &Controls) -> Schedule {
    let days = period_days(c);
    let n = days.len();
    let visit = |num: i64, name: String, day: i32, epoch, kind, labs| PlannedVisit {
        visitnum: num,
        name,
        day: Some(day),
        window: (0, 0),
        epoch,
        kind,
        labs,
        vitals: true,
        endpoint: false,
    };
    let mut visits = vec![PlannedVisit {
        window: (7, 0),
        ..visit(
            1,
            "SCREENING".into(),
            -21,
            Epoch::Screening,
            VisitKind::Screening,
            true,
        )
    }];
    for (k, &d) in days.iter().enumerate() {
        let p = k as i64 + 1;
        let epoch = Epoch::Period(p as u8);
        let kind = if k == 0 {
            VisitKind::Baseline
        } else {
            VisitKind::Treatment
        };
        visits.push(visit(
            10 * p,
            format!("PERIOD {p} DAY 1"),
            d,
            epoch,
            kind,
            k == 0,
        ));
        let last = k + 1 == n;
        visits.push(visit(
            10 * p + 1,
            format!("PERIOD {p} DAY {}", IN_HOUSE_DAYS + 1),
            d + IN_HOUSE_DAYS,
            epoch,
            if last {
                VisitKind::EndOfTreatment
            } else {
                VisitKind::Treatment
            },
            last,
        ));
    }
    let end = days[n - 1] + IN_HOUSE_DAYS;
    visits.push(PlannedVisit {
        window: (2, 2),
        ..visit(
            90,
            "FOLLOW-UP".into(),
            end + 7,
            Epoch::FollowUp,
            VisitKind::FollowUp,
            true,
        )
    });
    visits.push(PlannedVisit {
        day: None,
        ..visit(
            99,
            "EARLY TERMINATION".into(),
            0,
            Epoch::Period(1),
            VisitKind::EarlyTermination,
            true,
        )
    });
    Schedule {
        visits,
        treatment_days: end,
        followup_days: 7,
    }
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, PkParams) {
    let (sequences, codes) = layout(c);
    let food = c.design == "food_effect";
    let treatment = |code: &str| {
        let form = if food || code == "T" {
            "TABLET"
        } else {
            "CAPSULE"
        };
        Treatment {
            name: DRUG.to_string(),
            dose: c.dose,
            unit: "mg".into(),
            form: form.into(),
            route: "ORAL".into(),
            frequency: "ONCE".into(),
            cycle_days: 1,
        }
    };
    let label = |code: &str| {
        codes
            .iter()
            .find(|(cd, _)| *cd == code)
            .map(|(_, l)| *l)
            .unwrap_or("")
    };
    let element = |code: &str| {
        let t = treatment(code);
        format!(
            "{DRUG} {} mg {} ({})",
            fmt_dose(c.dose),
            t.form.to_lowercase(),
            label(code)
        )
    };
    let arms: Vec<Arm> = sequences
        .iter()
        .map(|seq| {
            let steps: Vec<&str> = if food {
                seq.chars()
                    .map(|ch| if ch == 'A' { "FED" } else { "FAST" })
                    .collect()
            } else {
                seq.chars()
                    .map(|ch| if ch == 'T' { "T" } else { "R" })
                    .collect()
            };
            let periods: Vec<Period> = steps
                .iter()
                .map(|code| Period {
                    etcd: (*code).to_string(),
                    element: element(code),
                    treatment: treatment(code),
                })
                .collect();
            let arm = steps
                .iter()
                .map(|code| capitalise(label(code)))
                .collect::<Vec<_>>()
                .join(" - ");
            Arm {
                armcd: if food {
                    steps.join("-")
                } else {
                    (*seq).to_string()
                },
                arm,
                ratio: 1,
                treatment: periods[0].treatment.clone(),
                periods,
            }
        })
        .collect();
    let first = Date::from_ymd(2024, 4, 8).expect("valid date");
    let enrollment_days = 120;
    let sch = schedule(c);
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| "SIMBE01".into()),
        seed: opts.seed,
        n: opts.n,
        arms,
        sites: vec![Site {
            siteid: "101".into(),
            country: "USA".into(),
        }],
        first_randomization: first,
        enrollment_days,
        cutoff: first.add_days(enrollment_days + sch.treatment_days + sch.followup_days + 30),
        block_size: 4,
        population: Population {
            age_mean: 33.0,
            age_sd: 9.0,
            age_min: 18,
            age_max: 55,
            sex: Categorical::new([("M", 0.55), ("F", 0.45)]),
            race: Categorical::new([
                ("WHITE", 0.62),
                ("BLACK OR AFRICAN AMERICAN", 0.25),
                ("ASIAN", 0.08),
                ("MULTIPLE", 0.05),
            ]),
            ethnic: Categorical::new([
                ("NOT HISPANIC OR LATINO", 0.8),
                ("HISPANIC OR LATINO", 0.2),
            ]),
            strata: Categorical::new([("ALL", 1.0)]),
        },
        timing: Timing {
            screening_lead: (2, 21),
            consent_lead: (0, 7),
            dose_lag: (0, 0),
        },
        schedule: Some(sch),
        cohorts: Vec::new(),
        info: TrialInfo {
            title: if food {
                format!(
                    "A Randomized, Open-Label, Two-Period Crossover Study of the Effect of Food on the Pharmacokinetics of {DRUG} in Healthy Adults"
                )
            } else {
                format!(
                    "A Randomized, Open-Label, Crossover Bioequivalence Study of Two Formulations of {DRUG} in Healthy Adults"
                )
            },
            phase: "PHASE I TRIAL".into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: "Healthy volunteers".into(),
            therapeutic_area: "CLINICAL PHARMACOLOGY".into(),
            primary_objective: if food {
                format!(
                    "To assess the effect of food on the pharmacokinetics of a single oral dose of {DRUG}"
                )
            } else {
                format!(
                    "To demonstrate bioequivalence of the test and reference formulations of {DRUG}"
                )
            },
            blinding: "OPEN LABEL".into(),
            control: "NONE".into(),
            add_on: false,
            trial_types: vec!["PHARMACOKINETIC".into()],
            intervention_model: "CROSSOVER".into(),
            healthy_subjects: true,
            criteria: vec![
                (
                    "INCL01".into(),
                    "INCLUSION".into(),
                    "Healthy adult aged 18 to 55 years".into(),
                ),
                (
                    "INCL02".into(),
                    "INCLUSION".into(),
                    "Body mass index 18.5 to 32 kg/m2".into(),
                ),
                (
                    "EXCL01".into(),
                    "EXCLUSION".into(),
                    "Use of any medication within 14 days of dosing".into(),
                ),
            ],
        },
    };
    let course = CourseParams {
        dropout_annual_rate: c.dropout_rate,
        dropout_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.7),
            ("PHYSICIAN DECISION".into(), 0.2),
            ("PROTOCOL DEVIATION".into(), 0.1),
        ],
        death_annual_rate: 0.0,
        missed_visit_prob: 0.0,
        ae_rate_multiplier: 1.0,
        ae_relative_risk: vec![HashMap::new(); design.arms.len()],
        ae_discontinuation: [0.0, 0.02, 0.3],
        ae_interruption: [0.0, 0.0, 0.0],
        required_history: Vec::new(),
        background_medication: Vec::new(),
        analyte_overrides: Vec::new(),
        endpoints: Vec::new(),
        bmi: (24.5, 3.0),
        dosing_days: None,
    };
    // The Test (or fed) element's relative F and absorption rate.
    let test = codes[0].0.to_string();
    let (relative_f, relative_ka) = (vec![(test.clone(), c.gmr)], vec![(test, c.ka_ratio)]);
    let pk = PkParams {
        analyte: DRUG.into(),
        compartments: 1,
        cl: c.cl,
        v: c.v,
        ka: c.ka,
        q: 1.0,
        v2: 1.0,
        cv: (c.inter_cv, 0.2, 0.3),
        proportionality: 1.0,
        reference_dose: c.dose,
        proportional_error: 0.08,
        additive_error: c.lloq * 0.2,
        lloq: c.lloq,
        profile_hours: vec![
            0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0,
        ],
        tau_hours: Vec::new(),
        tau: 24.0,
        planned_doses: period_days(c).len() as u32,
        dose_time: (8, 0),
        dose_days: Some(period_days(c)),
        relative_f,
        relative_ka,
        iov_cv: c.intra_cv,
    };
    (design, course, pk)
}

fn capitalise(s: &str) -> String {
    let mut c = s.chars();
    c.next().map_or(String::new(), |f| {
        f.to_uppercase().collect::<String>() + c.as_str()
    })
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}
