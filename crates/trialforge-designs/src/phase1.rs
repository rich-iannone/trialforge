//! `phase1_sad_mad`: a first-in-human study in healthy volunteers. Sequential dose cohorts
//! (6 active : 2 placebo by default) receive a single ascending dose (`part = "sad"`) or once
//! daily multiple ascending doses (`part = "mad"`), with intensive PK sampling, safety labs
//! and vital signs during an in-house period, a discharge visit and a follow-up visit.
//!
//! PK follows a population model (see `trialforge_processes::pk`) with controls in its own
//! terms: apparent clearance, volume and absorption rate, their between-subject variability,
//! and a dose-proportionality exponent.

use std::collections::HashMap;

use serde::{Deserialize, Serialize};
use trialforge_core::Date;
use trialforge_model::{
    Arm, Categorical, Cohort, Design, Epoch, PlannedVisit, Population, Schedule, Site, Timing,
    Treatment, TrialInfo, VisitKind,
};
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::pk::PkParams;

use crate::{ControlDoc, StudyOptions};

/// The investigational product.
pub const DRUG: &str = "TFX-101";

/// Controls for `phase1_sad_mad`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// `sad` (single ascending dose) or `mad` (multiple ascending dose, once daily).
    pub part: String,
    /// Dose (mg) of each cohort, in escalation order; `n` is split evenly across cohorts.
    pub doses: Vec<f64>,
    /// Share of each cohort on placebo.
    pub placebo_fraction: f64,
    /// Days of once-daily dosing in the MAD part.
    pub mad_days: u32,
    /// 1 or 2 compartments.
    pub compartments: u8,
    /// Apparent clearance CL/F (L/h) of a 70 kg subject at the reference (first) dose.
    pub cl: f64,
    /// Apparent central volume V/F (L) of a 70 kg subject.
    pub v: f64,
    /// Absorption rate constant (1/h).
    pub ka: f64,
    /// Apparent intercompartmental clearance Q/F (L/h); two compartments only.
    pub q: f64,
    /// Apparent peripheral volume V2/F (L); two compartments only.
    pub v2: f64,
    /// Between-subject CV of [CL, V, ka].
    pub cv: [f64; 3],
    /// Dose-proportionality exponent β: AUC ∝ dose^β (through clearance).
    pub proportionality: f64,
    /// Proportional residual error (CV).
    pub proportional_error: f64,
    /// Lower limit of quantitation (ng/mL).
    pub lloq: f64,
    /// Multiplier on every background AE rate.
    pub ae_rate_multiplier: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            part: "sad".into(),
            doses: vec![10.0, 30.0, 100.0, 300.0, 600.0],
            placebo_fraction: 0.25,
            mad_days: 7,
            compartments: 1,
            cl: 12.0,
            v: 90.0,
            ka: 1.2,
            q: 6.0,
            v2: 150.0,
            cv: [0.3, 0.25, 0.4],
            proportionality: 1.0,
            proportional_error: 0.12,
            lloq: 0.5,
            ae_rate_multiplier: 1.0,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "part",
            "sad (single ascending dose) or mad (multiple ascending dose, once daily)",
        ),
        ControlDoc::new(
            "doses",
            "Dose in mg of each cohort, in escalation order; n is split evenly across cohorts",
        ),
        ControlDoc::new(
            "placebo_fraction",
            "Share of each cohort randomised to placebo",
        ),
        ControlDoc::new("mad_days", "Days of once-daily dosing in the MAD part"),
        ControlDoc::new(
            "compartments",
            "PK model: 1 or 2 compartments, first-order absorption",
        ),
        ControlDoc::new(
            "cl",
            "Apparent clearance CL/F (L/h) of a 70 kg subject at the first dose",
        ),
        ControlDoc::new("v", "Apparent central volume V/F (L) of a 70 kg subject"),
        ControlDoc::new("ka", "Absorption rate constant (1/h)"),
        ControlDoc::new(
            "q",
            "Apparent intercompartmental clearance Q/F (L/h); 2 compartments only",
        ),
        ControlDoc::new(
            "v2",
            "Apparent peripheral volume V2/F (L); 2 compartments only",
        ),
        ControlDoc::new(
            "cv",
            "Between-subject coefficients of variation of [CL, V, ka]",
        ),
        ControlDoc::new(
            "proportionality",
            "Dose-proportionality exponent: AUC grows as dose^proportionality (through clearance)",
        ),
        ControlDoc::new(
            "proportional_error",
            "Proportional residual error (CV) of concentrations",
        ),
        ControlDoc::new("lloq", "Lower limit of quantitation (ng/mL)"),
        ControlDoc::new(
            "ae_rate_multiplier",
            "Multiplier on every background adverse event rate",
        ),
    ]
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    if c.part != "sad" && c.part != "mad" {
        return Err("`part` must be \"sad\" or \"mad\"".into());
    }
    if c.doses.is_empty() || c.doses.len() > 12 || !c.doses.iter().all(|d| positive(*d)) {
        return Err("`doses` needs 1–12 positive doses".into());
    }
    if !(0.0..0.9).contains(&c.placebo_fraction) {
        return Err("`placebo_fraction` must be in [0, 0.9)".into());
    }
    if !(2..=28).contains(&c.mad_days) {
        return Err("`mad_days` must be between 2 and 28".into());
    }
    if c.compartments != 1 && c.compartments != 2 {
        return Err("`compartments` must be 1 or 2".into());
    }
    for (name, x) in [
        ("cl", c.cl),
        ("v", c.v),
        ("ka", c.ka),
        ("q", c.q),
        ("v2", c.v2),
        ("lloq", c.lloq),
    ] {
        if !positive(x) {
            return Err(format!("`{name}` must be positive"));
        }
    }
    if !c.cv.iter().all(|x| (0.0..=1.5).contains(x)) {
        return Err("`cv` values must be between 0 and 1.5".into());
    }
    if !(0.3..=1.7).contains(&c.proportionality) {
        return Err("`proportionality` must be between 0.3 and 1.7".into());
    }
    if !(0.0..=0.5).contains(&c.proportional_error) {
        return Err("`proportional_error` must be between 0 and 0.5".into());
    }
    if !(0.0..=10.0).contains(&c.ae_rate_multiplier) {
        return Err("`ae_rate_multiplier` must be between 0 and 10".into());
    }
    Ok(())
}

/// Planned doses.
fn planned_doses(c: &Controls) -> u32 {
    if c.part == "sad" { 1 } else { c.mad_days }
}

fn schedule(c: &Controls) -> Schedule {
    let n_doses = planned_doses(c) as i32;
    let discharge = n_doses + 3;
    let visit =
        |num: i64, name: String, day: Option<i32>, window, epoch, kind, labs| PlannedVisit {
            visitnum: num,
            name,
            day,
            window,
            epoch,
            kind,
            labs,
            vitals: true,
            endpoint: false,
        };
    let mut visits = vec![
        visit(
            1,
            "SCREENING".into(),
            Some(-21),
            (7, 0),
            Epoch::Screening,
            VisitKind::Screening,
            true,
        ),
        visit(
            2,
            "DAY 1".into(),
            Some(1),
            (0, 0),
            Epoch::Treatment,
            VisitKind::Baseline,
            true,
        ),
    ];
    for d in 2..=n_doses {
        let labs = d == n_doses || d % 7 == 0;
        visits.push(visit(
            1 + i64::from(d),
            format!("DAY {d}"),
            Some(d),
            (0, 0),
            Epoch::Treatment,
            VisitKind::Treatment,
            labs,
        ));
    }
    visits.push(visit(
        1 + i64::from(discharge),
        format!("DAY {discharge}"),
        Some(discharge),
        (0, 0),
        Epoch::Treatment,
        VisitKind::EndOfTreatment,
        true,
    ));
    visits.push(visit(
        90,
        "FOLLOW-UP".into(),
        Some(discharge + 7),
        (2, 2),
        Epoch::FollowUp,
        VisitKind::FollowUp,
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
    ));
    Schedule {
        visits,
        treatment_days: discharge,
        followup_days: 7,
    }
}

/// Cohort sizes: `n` split evenly across dose levels (earlier cohorts take the remainder).
fn cohort_sizes(n: usize, levels: usize) -> Vec<usize> {
    (0..levels)
        .map(|k| n / levels + usize::from(k < n % levels))
        .collect()
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, PkParams) {
    let mad = c.part == "mad";
    let frequency = if mad { "QD" } else { "ONCE" };
    let treatment = |name: &str, dose: f64| Treatment {
        name: name.into(),
        dose,
        unit: "mg".into(),
        form: "TABLET".into(),
        route: "ORAL".into(),
        frequency: frequency.into(),
        cycle_days: 1,
    };
    let mut arms: Vec<Arm> = c
        .doses
        .iter()
        .enumerate()
        .map(|(k, d)| Arm {
            armcd: format!("C{}", k + 1),
            arm: format!("{DRUG} {} mg", fmt_dose(*d)),
            ratio: 1,
            treatment: treatment(DRUG, *d),
            periods: Vec::new(),
            switch_to: None,
        })
        .collect();
    let placebo = arms.len();
    arms.push(Arm {
        armcd: "PBO".into(),
        arm: "Placebo".into(),
        ratio: 1,
        treatment: treatment("PLACEBO", 0.0),
        periods: Vec::new(),
        switch_to: None,
    });
    let spacing = if mad { 35 } else { 21 };
    let cohorts: Vec<Cohort> = cohort_sizes(opts.n, c.doses.len())
        .into_iter()
        .enumerate()
        .filter(|(_, size)| *size > 0)
        .map(|(k, size)| {
            let pbo =
                ((size as f64 * c.placebo_fraction).round() as usize).min(size.saturating_sub(1));
            let mut arms = vec![(k, (size - pbo) as u32)];
            if pbo > 0 {
                arms.push((placebo, pbo as u32));
            }
            Cohort {
                arms,
                start_day: spacing * k as i32,
                enrollment_days: 3,
            }
        })
        .collect();
    let first = Date::from_ymd(2024, 1, 8).expect("valid date");
    let enrollment_days = spacing * c.doses.len() as i32;
    let sch = schedule(c);
    let study_days = sch.treatment_days + sch.followup_days + 30;
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| "SIMPK01".into()),
        seed: opts.seed,
        n: opts.n,
        arms,
        sites: vec![Site {
            siteid: "101".into(),
            country: "USA".into(),
        }],
        first_randomization: first,
        enrollment_days,
        cutoff: first.add_days(enrollment_days + study_days),
        interim: false,
        block_size: (c.doses.len() + 1) as u32,
        population: Population {
            age_mean: 32.0,
            age_sd: 8.0,
            age_min: 18,
            age_max: 55,
            sex: Categorical::new([("M", 0.6), ("F", 0.4)]),
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
        cohorts,
        info: TrialInfo {
            title: format!(
                "A Randomized, Double-Blind, Placebo-Controlled, {} Ascending Dose Study of the Safety, Tolerability and Pharmacokinetics of {DRUG} in Healthy Adults",
                if mad { "Multiple" } else { "Single" }
            ),
            phase: "PHASE I TRIAL".into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: "Healthy volunteers".into(),
            therapeutic_area: "CLINICAL PHARMACOLOGY".into(),
            primary_objective: format!(
                "To evaluate the safety, tolerability and pharmacokinetics of {} doses of {DRUG}",
                if mad { "multiple" } else { "single" }
            ),
            blinding: "DOUBLE BLIND".into(),
            control: "PLACEBO".into(),
            add_on: false,
            trial_types: vec!["SAFETY".into(), "PHARMACOKINETIC".into()],
            intervention_model: "SEQUENTIAL".into(),
            healthy_subjects: true,
            intervention_type: "DRUG".into(),
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
                    "Clinically significant abnormal laboratory value at screening".into(),
                ),
                (
                    "EXCL02".into(),
                    "EXCLUSION".into(),
                    "Use of any medication within 14 days of dosing".into(),
                ),
            ],
        },
    };
    // Dose-related tolerability signal for the active arms.
    let reference = c.doses[0];
    let ae_relative_risk: Vec<HashMap<String, f64>> = design
        .arms
        .iter()
        .map(|a| {
            if a.treatment.dose > 0.0 {
                let scale = 1.0
                    + 0.5 * trialforge_core::math::ln(1.0 + a.treatment.dose / reference).min(4.0);
                [("Headache", 1.5), ("Nausea", 2.0), ("Dizziness", 1.5)]
                    .iter()
                    .map(|(pt, r)| ((*pt).to_string(), 1.0 + (r - 1.0) * scale))
                    .collect()
            } else {
                HashMap::new()
            }
        })
        .collect();
    let course = CourseParams {
        dropout_annual_rate: 0.8,
        dropout_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.7),
            ("PHYSICIAN DECISION".into(), 0.2),
            ("PROTOCOL DEVIATION".into(), 0.1),
        ],
        death_annual_rate: 0.0,
        missed_visit_prob: 0.0,
        ae_rate_multiplier: c.ae_rate_multiplier,
        ae_relative_risk,
        ae_discontinuation: [0.0, 0.02, 0.3],
        ae_interruption: [0.0, 0.0, 0.0],
        required_history: Vec::new(),
        background_medication: Vec::new(),
        analyte_overrides: Vec::new(),
        endpoints: Vec::new(),
        bmi: (24.5, 3.0),
        dosing_days: Some(planned_doses(c) as i32),
        extra_ae: Vec::new(),
    };
    let pk = PkParams {
        analyte: DRUG.into(),
        compartments: c.compartments,
        cl: c.cl,
        v: c.v,
        ka: c.ka,
        q: c.q,
        v2: c.v2,
        cv: (c.cv[0], c.cv[1], c.cv[2]),
        proportionality: c.proportionality,
        reference_dose: reference,
        proportional_error: c.proportional_error,
        additive_error: c.lloq * 0.2,
        lloq: c.lloq,
        profile_hours: vec![
            0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0, 36.0, 48.0, 72.0,
        ],
        tau_hours: vec![0.25, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0, 12.0, 24.0],
        tau: 24.0,
        planned_doses: planned_doses(c),
        dose_time: (8, 0),
        dose_days: None,
        relative_f: Vec::new(),
        relative_ka: Vec::new(),
        iov_cv: 0.0,
    };
    (design, course, pk)
}

/// A dose without a trailing `.0`.
fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}
