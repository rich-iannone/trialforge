//! `cns_alzheimers`: a randomised, double-blind, placebo-controlled study of two doses of Drug A
//! in mild to moderate Alzheimer's disease (CDISC pilot-like), with ADAS-Cog (11) assessed at
//! baseline and every 8 weeks.
//!
//! ADAS-Cog items are generated from a latent cognitive trajectory with a graded response
//! model (see `trialforge_processes::adas`); the controls are in ADAS-Cog points.

use serde::{Deserialize, Serialize};
use trialforge_model::Design;
use trialforge_processes::adas::AdasParams;
use trialforge_processes::fixed::CourseParams;

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `cns_alzheimers`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Treatment duration in weeks (a multiple of 8).
    pub weeks: u32,
    /// Doses (mg, once daily) of the active arms.
    pub doses: Vec<f64>,
    /// Mean ADAS-Cog (11) total at baseline (0–70, higher is worse).
    pub baseline_mean: f64,
    /// Mean change in the total at the end of treatment: [placebo, then each dose].
    pub change: Vec<f64>,
    /// SD of each subject's rate of decline around the mean (ADAS-Cog points per year at the
    /// population mean).
    pub decline_sd: f64,
    /// Probability an item is not done.
    pub item_missing: f64,
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
            weeks: 24,
            doses: vec![5.0, 10.0],
            baseline_mean: 24.0,
            change: vec![2.5, 1.2, 0.3],
            decline_sd: 3.0,
            item_missing: 0.01,
            dropout_rate: 0.3,
            ae_rate_multiplier: 1.2,
            missed_visit: 0.04,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "weeks",
            "Treatment duration in weeks (a multiple of 8; ADAS-Cog every 8 weeks)",
        ),
        ControlDoc::new("doses", "Doses in mg (once daily) of the active arms"),
        ControlDoc::new(
            "baseline_mean",
            "Mean ADAS-Cog (11) total at baseline (0-70, higher is worse)",
        ),
        ControlDoc::new(
            "change",
            "Mean change in the ADAS-Cog total at the end of treatment: [placebo, each dose]",
        ),
        ControlDoc::new(
            "decline_sd",
            "SD of subjects' rate of decline (ADAS-Cog points per year)",
        ),
        ControlDoc::new("item_missing", "Probability an ADAS-Cog item is not done"),
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
    if !(8..=104).contains(&c.weeks) || c.weeks % 8 != 0 {
        return Err("`weeks` must be a multiple of 8 between 8 and 104".into());
    }
    if c.doses.is_empty() || c.doses.len() > 4 || !c.doses.iter().all(|d| d.is_finite() && *d > 0.0)
    {
        return Err("`doses` needs 1 to 4 positive doses".into());
    }
    if c.change.len() != c.doses.len() + 1
        || !c.change.iter().all(|x| x.is_finite() && x.abs() <= 30.0)
    {
        return Err("`change` needs one value per arm (placebo first), each within ±30".into());
    }
    if !(5.0..=60.0).contains(&c.baseline_mean) {
        return Err("`baseline_mean` must be between 5 and 60".into());
    }
    if !(0.0..=20.0).contains(&c.decline_sd) || !(0.0..=0.3).contains(&c.item_missing) {
        return Err("`decline_sd` must be in [0, 20] and `item_missing` in [0, 0.3]".into());
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

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, AdasParams) {
    let top = c.doses.iter().copied().fold(0.0, f64::max);
    let mut arms = vec![ArmSpec {
        armcd: "PBO".into(),
        arm: "Placebo".into(),
        treatment: "PLACEBO".into(),
        dose: 0.0,
        ratio: 1,
        change: 0.0,
        ae_share: 0.0,
    }];
    for (k, &d) in c.doses.iter().enumerate() {
        let label = fmt_dose(d);
        arms.push(ArmSpec {
            armcd: format!("D{}", k + 1),
            arm: format!("Drug A {label} mg"),
            treatment: "DRUG A".into(),
            dose: d,
            ratio: 1,
            change: 0.0,
            ae_share: d / top,
        });
    }
    let (design, course) = build(
        opts,
        &Common {
            indication: "alzheimers",
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMCNS01",
            phase: "PHASE II TRIAL",
            control: "PLACEBO",
            title: None,
        },
        &arms,
    );
    let treatment_days = (c.weeks * 7) as i32;
    // ADAS-Cog at baseline and every 8 weeks (the schedule has a visit on each of these days).
    let visit_days: Vec<i32> = std::iter::once(1)
        .chain((1..=c.weeks / 8).map(|k| (k * 56 + 1) as i32))
        .collect();
    // Decline SD in points per year → latent units via the total's slope near the mean
    // (≈ 9 points per latent unit for these items).
    let adas = AdasParams {
        baseline_mean: c.baseline_mean,
        change_by_arm: c.change.clone(),
        visit_days,
        treatment_days,
        onset_days: 28.0,
        slope_sd: c.decline_sd / 9.0,
        visit_sd: 0.25,
        item_missing: c.item_missing,
    };
    (design, course, adas)
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}
