//! `dose_ranging_ph2`: a randomised, double-blind, placebo-controlled phase 2 study of several
//! doses with a continuous endpoint, for dose-response modelling (MCP-Mod style).
//!
//! The mean change at the end of treatment follows a sigmoid Emax model of dose:
//! `placebo_change + emax · d^hill / (ed50^hill + d^hill)`. The indication presets, schedule
//! and safety profile are those of `parallel_continuous`; the active-drug AE signal grows with
//! dose.

use serde::{Deserialize, Serialize};
use trialforge_model::Design;
use trialforge_processes::fixed::CourseParams;

use crate::parallel::{ArmSpec, Common, build, preset_change};
use crate::{ControlDoc, StudyOptions};

/// Controls for `dose_ranging_ph2`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// `type_2_diabetes`, `hypertension` or `obesity`.
    pub indication: String,
    /// Treatment duration in weeks.
    pub weeks: u32,
    /// Doses (mg) including 0 for placebo, in increasing order.
    pub doses: Vec<f64>,
    /// Mean change from baseline on placebo; `null` uses the indication's.
    pub placebo_change: Option<f64>,
    /// Maximum effect beyond placebo (endpoint units); `null` uses 1.25 times the indication's
    /// active-minus-placebo difference.
    pub emax: Option<f64>,
    /// Dose (mg) giving half the maximum effect.
    pub ed50: f64,
    /// Hill coefficient (steepness).
    pub hill: f64,
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
            weeks: 12,
            doses: vec![0.0, 2.5, 5.0, 10.0, 20.0, 40.0],
            placebo_change: None,
            emax: None,
            ed50: 7.5,
            hill: 1.0,
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
        ControlDoc::new("weeks", "Treatment duration in weeks"),
        ControlDoc::new("doses", "Doses in mg, increasing, including 0 for placebo"),
        ControlDoc::new(
            "placebo_change",
            "Mean change from baseline on placebo; null for the indication's",
        ),
        ControlDoc::new(
            "emax",
            "Maximum effect beyond placebo (endpoint units); null for 1.25 times the indication's difference",
        ),
        ControlDoc::new("ed50", "Dose (mg) with half the maximum effect"),
        ControlDoc::new("hill", "Hill coefficient of the sigmoid Emax model"),
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
    if preset_change(&c.indication).is_none() || c.indication == "alzheimers" {
        return Err(format!(
            "unknown indication {:?}; available: type_2_diabetes, hypertension, obesity",
            c.indication
        ));
    }
    if !(4..=104).contains(&c.weeks) {
        return Err("`weeks` must be between 4 and 104".into());
    }
    if c.doses.len() < 2 || c.doses.len() > 8 {
        return Err("`doses` needs 2 to 8 doses".into());
    }
    if !c.doses.iter().all(|d| d.is_finite() && *d >= 0.0)
        || !c.doses.windows(2).all(|w| w[0] < w[1])
    {
        return Err("`doses` must be non-negative and increasing".into());
    }
    for (name, x) in [("ed50", c.ed50), ("hill", c.hill)] {
        if !(x.is_finite() && x > 0.0) {
            return Err(format!("`{name}` must be positive"));
        }
    }
    if c.hill > 10.0 {
        return Err("`hill` must be at most 10".into());
    }
    if [c.placebo_change, c.emax]
        .iter()
        .flatten()
        .any(|x| !x.is_finite())
    {
        return Err("`placebo_change` and `emax` must be finite".into());
    }
    if !(0.0..=5.0).contains(&c.dropout_rate) {
        return Err("`dropout_rate` must be between 0 and 5 per year".into());
    }
    if !(0.0..=20.0).contains(&c.ae_rate_multiplier) {
        return Err("`ae_rate_multiplier` must be between 0 and 20".into());
    }
    if !(0.0..=0.5).contains(&c.missed_visit) {
        return Err("`missed_visit` must be between 0 and 0.5".into());
    }
    Ok(())
}

/// (E0, Emax, ED50, Hill) of the dose-response model.
pub fn model(c: &Controls) -> (f64, f64, f64, f64) {
    let pre = preset_change(&c.indication).expect("validated");
    let e0 = c.placebo_change.unwrap_or(pre[0]);
    let emax = c.emax.unwrap_or(1.25 * (pre[1] - pre[0]));
    (e0, emax, c.ed50, c.hill)
}

/// Mean change at dose `d`.
pub fn response(model: (f64, f64, f64, f64), d: f64) -> f64 {
    let (e0, emax, ed50, hill) = model;
    if d <= 0.0 {
        return e0;
    }
    let dh = trialforge_core::math::pow(d, hill);
    e0 + emax * dh / (trialforge_core::math::pow(ed50, hill) + dh)
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams) {
    let m = model(c);
    let top = c.doses.iter().copied().fold(0.0, f64::max);
    let arms: Vec<ArmSpec> = c
        .doses
        .iter()
        .map(|&d| {
            let label = fmt_dose(d);
            if d <= 0.0 {
                ArmSpec {
                    armcd: "PBO".into(),
                    arm: "Placebo".into(),
                    treatment: "PLACEBO".into(),
                    dose: 0.0,
                    ratio: 1,
                    change: response(m, 0.0),
                    ae_share: 0.0,
                }
            } else {
                ArmSpec {
                    armcd: format!("D{}", label.replace('.', "P"))
                        .chars()
                        .take(8)
                        .collect(),
                    arm: format!("Drug A {label} mg"),
                    treatment: "DRUG A".into(),
                    dose: d,
                    ratio: 1,
                    change: response(m, d),
                    ae_share: d / top,
                }
            }
        })
        .collect();
    build(
        opts,
        &Common {
            indication: &c.indication,
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMDR01",
            phase: "PHASE II TRIAL",
            control: "PLACEBO",
            title: Some((
                "A Randomized, Double-Blind, Placebo-Controlled, Dose-Ranging Study of Drug A".into(),
                "To characterise the dose-response of Drug A on the change from baseline in the primary endpoint at the end of treatment".into(),
            )),
        },
        &arms,
    )
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}
