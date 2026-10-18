//! `pediatric_growth`: a randomised, open-label, active-controlled non-inferiority study of
//! once-weekly Drug A versus daily somatropin in prepubertal children with growth hormone
//! deficiency. Height, weight and BMI are measured every 13 weeks; ADVS adds height-, weight-
//! and BMI-for-age z-scores and percentiles (CDC 2000 charts, admiralpeds) and the annualized
//! height velocity, the primary endpoint.
//!
//! Growth comes from `trialforge_processes::growth`; the controls are in cm/yr and
//! height-for-age z-scores.

use serde::{Deserialize, Serialize};
use trialforge_model::{Design, Epoch, PlannedVisit, Schedule, VisitKind};
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::growth::GrowthParams;
use trialforge_processes::reference::AeBackground;

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `pediatric_growth`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Treatment duration in weeks (a multiple of 13; visits every 13 weeks).
    pub weeks: u32,
    /// Age range at enrolment in years: [youngest, oldest].
    pub age_range: [i64; 2],
    /// Drug A dose (mg/kg, once weekly).
    pub dose: f64,
    /// Somatropin dose (mg/kg, once daily).
    pub comparator_dose: f64,
    /// Mean annualized height velocity on treatment: [somatropin, Drug A] (cm/yr).
    pub height_velocity: [f64; 2],
    /// SD of children's height velocity around their arm's mean (cm/yr).
    pub velocity_sd: f64,
    /// Mean height-for-age z-score at baseline.
    pub height_sds: f64,
    /// Inclusion: the highest screening height-for-age z-score.
    pub max_screening_sds: f64,
    /// Height velocity before treatment (cm/yr).
    pub untreated_velocity: f64,
    /// Non-inferiority margin for height velocity, Drug A minus somatropin (cm/yr; higher is
    /// better), or `null`.
    pub ni_margin: Option<f64>,
    /// Allocation ratio: [somatropin, Drug A].
    pub ratio: [u32; 2],
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
            weeks: 52,
            age_range: [3, 11],
            dose: 0.24,
            comparator_dose: 0.034,
            height_velocity: [10.3, 11.2],
            velocity_sd: 2.0,
            height_sds: -2.9,
            max_screening_sds: -2.0,
            untreated_velocity: 4.0,
            ni_margin: Some(2.0),
            ratio: [1, 2],
            dropout_rate: 0.05,
            ae_rate_multiplier: 1.0,
            missed_visit: 0.02,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "weeks",
            "Treatment duration in weeks (a multiple of 13; visits every 13 weeks)",
        ),
        ControlDoc::new(
            "age_range",
            "Age range at enrolment in years (3 to 15): [youngest, oldest]",
        ),
        ControlDoc::new("dose", "Drug A dose in mg/kg, once weekly"),
        ControlDoc::new("comparator_dose", "Somatropin dose in mg/kg, once daily"),
        ControlDoc::new(
            "height_velocity",
            "Mean annualized height velocity on treatment: [somatropin, Drug A] (cm/yr)",
        ),
        ControlDoc::new(
            "velocity_sd",
            "SD of children's height velocity around their arm's mean (cm/yr)",
        ),
        ControlDoc::new(
            "height_sds",
            "Mean height-for-age z-score at baseline (CDC 2000 charts)",
        ),
        ControlDoc::new(
            "max_screening_sds",
            "Inclusion: the highest screening height-for-age z-score",
        ),
        ControlDoc::new(
            "untreated_velocity",
            "Height velocity before treatment (cm/yr)",
        ),
        ControlDoc::new(
            "ni_margin",
            "Non-inferiority margin for height velocity, Drug A minus somatropin (cm/yr); null for none",
        ),
        ControlDoc::new("ratio", "Allocation ratio: [somatropin, Drug A]"),
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

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    if !(26..=104).contains(&c.weeks) || !c.weeks.is_multiple_of(13) {
        return Err("`weeks` must be a multiple of 13 between 26 and 104".into());
    }
    if c.age_range[0] < 3 || c.age_range[1] > 15 || c.age_range[0] > c.age_range[1] {
        return Err("`age_range` must be [youngest, oldest] within 3 to 15 years".into());
    }
    let pos = |x: f64| x.is_finite() && x > 0.0;
    if !pos(c.dose) || !pos(c.comparator_dose) {
        return Err("`dose` and `comparator_dose` must be positive".into());
    }
    if !c
        .height_velocity
        .iter()
        .all(|v| v.is_finite() && (1.0..=20.0).contains(v))
    {
        return Err("`height_velocity` entries must be between 1 and 20 cm/yr".into());
    }
    if !(0.0..=5.0).contains(&c.velocity_sd) {
        return Err("`velocity_sd` must be between 0 and 5 cm/yr".into());
    }
    if !(-4.5..=-1.0).contains(&c.height_sds) || !(-4.0..=0.0).contains(&c.max_screening_sds) {
        return Err("`height_sds` must be in [-4.5, -1] and `max_screening_sds` in [-4, 0]".into());
    }
    if !(0.5..=10.0).contains(&c.untreated_velocity) {
        return Err("`untreated_velocity` must be between 0.5 and 10 cm/yr".into());
    }
    if c.ni_margin.is_some_and(|m| !pos(m)) {
        return Err("`ni_margin` must be positive".into());
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

/// Screening, baseline, then a visit every 13 weeks to the end of treatment.
fn schedule(weeks: u32) -> Schedule {
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
            endpoint: true,
        };
    let mut visits = vec![
        visit(
            1,
            "SCREENING".into(),
            Some(-28),
            (14, 0),
            Epoch::Screening,
            VisitKind::Screening,
            true,
        ),
        visit(
            2,
            "BASELINE".into(),
            Some(1),
            (0, 0),
            Epoch::Treatment,
            VisitKind::Baseline,
            true,
        ),
    ];
    for (k, w) in (13..=weeks).step_by(13).enumerate() {
        let last = w == weeks;
        visits.push(visit(
            3 + k as i64,
            format!("WEEK {w}"),
            Some((w * 7 + 1) as i32),
            (7, 7),
            Epoch::Treatment,
            if last {
                VisitKind::EndOfTreatment
            } else {
                VisitKind::Treatment
            },
            last || w == 26,
        ));
    }
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
        treatment_days: (weeks * 7) as i32,
        followup_days: 0,
    }
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, GrowthParams) {
    let arms = [
        ArmSpec {
            armcd: "SOMA".into(),
            arm: format!("Somatropin {} mg/kg/day", fmt_dose(c.comparator_dose)),
            treatment: "SOMATROPIN".into(),
            dose: c.comparator_dose,
            ratio: c.ratio[0],
            change: 0.0,
            ae_share: 0.0,
        },
        ArmSpec {
            armcd: "DRUGA".into(),
            arm: format!("Drug A {} mg/kg/week", fmt_dose(c.dose)),
            treatment: "DRUG A".into(),
            dose: c.dose,
            ratio: c.ratio[1],
            change: 0.0,
            ae_share: 0.0,
        },
    ];
    let (mut design, mut course) = build(
        opts,
        &Common {
            indication: "growth_hormone_deficiency",
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMPED01",
            phase: "PHASE III TRIAL",
            control: "ACTIVE",
            title: None,
        },
        &arms,
    );
    design.schedule = Some(schedule(c.weeks));
    design.info.blinding = "OPEN LABEL".into();
    let pop = &mut design.population;
    pop.age_min = c.age_range[0];
    pop.age_max = c.age_range[1];
    pop.age_mean = (c.age_range[0] + c.age_range[1]) as f64 / 2.0;
    pop.age_sd = ((c.age_range[1] - c.age_range[0]) as f64 / 3.5).max(0.5);
    pop.sex = trialforge_model::Categorical::new([("F", 0.35), ("M", 0.65)]);
    for (arm, frequency) in design.arms.iter_mut().zip(["QD", "EVERY WEEK"]) {
        let t = &mut arm.treatment;
        t.unit = "mg/kg".into();
        t.form = "INJECTION".into();
        t.route = "SUBCUTANEOUS".into();
        t.frequency = frequency.into();
        t.cycle_days = if frequency == "QD" { 1 } else { 7 };
    }
    course.death_annual_rate = 0.0002;
    course.extra_ae = vec![
        AeBackground {
            pt: "Injection site reaction".into(),
            rate_per_year: 0.4,
            median_days: 2.0,
            p_serious: 0.0,
            severity: [0.85, 0.14, 0.01],
        },
        AeBackground {
            pt: "Hypothyroidism".into(),
            rate_per_year: 0.03,
            median_days: 300.0,
            p_serious: 0.0,
            severity: [0.6, 0.38, 0.02],
        },
    ];
    let growth = GrowthParams {
        height_sds: (c.height_sds, 0.6),
        max_screening_sds: c.max_screening_sds,
        velocity_by_arm: c.height_velocity.to_vec(),
        velocity_sd: c.velocity_sd,
        untreated_velocity: c.untreated_velocity,
    };
    (design, course, growth)
}
