//! `respiratory_exacerbation`: a randomised, double-blind, placebo-controlled study of Drug A
//! (subcutaneous, every 4 weeks) added to inhaled maintenance therapy in COPD with a history of
//! exacerbations. The primary endpoint is the annualized rate of moderate or severe
//! exacerbations (negative binomial counts, CE); trough FEV1 (RE) is a secondary endpoint.
//! Randomisation is stratified by exacerbations in the prior year (1 or ≥2).
//!
//! Exacerbations and FEV1 come from `trialforge_processes::exacerbation`.

use serde::{Deserialize, Serialize};
use trialforge_model::{Categorical, Design, Epoch, PlannedVisit, Schedule, VisitKind};
use trialforge_processes::exacerbation::ExacerbationParams;
use trialforge_processes::fixed::CourseParams;

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `respiratory_exacerbation`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Treatment duration in weeks (a multiple of 4).
    pub weeks: u32,
    /// Drug A dose (mg every 4 weeks).
    pub dose: f64,
    /// Annual rate of moderate or severe exacerbations on placebo (one prior exacerbation).
    pub exacerbation_rate: f64,
    /// Exacerbation rate ratio, Drug A vs placebo.
    pub rate_ratio: f64,
    /// Negative binomial overdispersion θ (Var = μ + θμ²).
    pub dispersion: f64,
    /// Share of exacerbations that are severe (hospitalised).
    pub severe_fraction: f64,
    /// Rate multiplier for subjects with two or more exacerbations in the prior year.
    pub prior_exacerbations_rr: f64,
    /// Mean trough FEV1 at baseline (L).
    pub fev1_baseline: f64,
    /// Mean change in trough FEV1 at the end of treatment: [placebo, Drug A] (L).
    pub fev1_change: [f64; 2],
    /// Allocation ratio: [placebo, Drug A].
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
            dose: 300.0,
            exacerbation_rate: 1.0,
            rate_ratio: 0.7,
            dispersion: 0.6,
            severe_fraction: 0.2,
            prior_exacerbations_rr: 1.6,
            fev1_baseline: 1.35,
            fev1_change: [-0.02, 0.08],
            ratio: [1, 1],
            dropout_rate: 0.12,
            ae_rate_multiplier: 1.0,
            missed_visit: 0.03,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new("weeks", "Treatment duration in weeks (a multiple of 4)"),
        ControlDoc::new("dose", "Drug A dose in mg, subcutaneous every 4 weeks"),
        ControlDoc::new(
            "exacerbation_rate",
            "Annual rate of moderate or severe exacerbations on placebo (one prior exacerbation)",
        ),
        ControlDoc::new("rate_ratio", "Exacerbation rate ratio, Drug A vs placebo"),
        ControlDoc::new(
            "dispersion",
            "Negative binomial overdispersion theta (Var = mu + theta mu^2)",
        ),
        ControlDoc::new(
            "severe_fraction",
            "Share of exacerbations that are severe (hospitalised)",
        ),
        ControlDoc::new(
            "prior_exacerbations_rr",
            "Rate multiplier for subjects with two or more exacerbations in the prior year (stratum)",
        ),
        ControlDoc::new("fev1_baseline", "Mean trough FEV1 at baseline (L)"),
        ControlDoc::new(
            "fev1_change",
            "Mean change in trough FEV1 at the end of treatment: [placebo, Drug A] (L)",
        ),
        ControlDoc::new("ratio", "Allocation ratio: [placebo, Drug A]"),
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
    if !(12..=104).contains(&c.weeks) || !c.weeks.is_multiple_of(4) {
        return Err("`weeks` must be a multiple of 4 between 12 and 104".into());
    }
    let pos = |x: f64| x.is_finite() && x > 0.0;
    if !pos(c.dose) {
        return Err("`dose` must be positive".into());
    }
    if !(c.exacerbation_rate.is_finite() && (0.0..=10.0).contains(&c.exacerbation_rate)) {
        return Err("`exacerbation_rate` must be between 0 and 10 per year".into());
    }
    if !pos(c.rate_ratio) || c.rate_ratio > 5.0 {
        return Err("`rate_ratio` must be positive and at most 5".into());
    }
    if !(0.0..=5.0).contains(&c.dispersion) || !(0.0..=1.0).contains(&c.severe_fraction) {
        return Err("`dispersion` must be in [0, 5] and `severe_fraction` in [0, 1]".into());
    }
    if !pos(c.prior_exacerbations_rr) || c.prior_exacerbations_rr > 10.0 {
        return Err("`prior_exacerbations_rr` must be positive and at most 10".into());
    }
    if !(0.5..=4.0).contains(&c.fev1_baseline)
        || !c
            .fev1_change
            .iter()
            .all(|x| x.is_finite() && x.abs() <= 1.0)
    {
        return Err("`fev1_baseline` must be in [0.5, 4] L and `fev1_change` within ±1 L".into());
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

/// Screening, baseline, Weeks 4, 8 and 12, then every 8 weeks to the end of treatment.
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
            Some(-14),
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
    let mut wks: Vec<u32> = vec![4, 8, 12];
    wks.extend((20..weeks).step_by(8));
    wks.retain(|&w| w < weeks);
    wks.push(weeks);
    for (k, w) in wks.iter().enumerate() {
        let last = *w == weeks;
        visits.push(visit(
            3 + k as i64,
            format!("WEEK {w}"),
            Some((w * 7 + 1) as i32),
            (5, 5),
            Epoch::Treatment,
            if last {
                VisitKind::EndOfTreatment
            } else {
                VisitKind::Treatment
            },
            last || *w == 12,
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

pub(crate) fn resolve(
    opts: &StudyOptions,
    c: &Controls,
) -> (Design, CourseParams, ExacerbationParams) {
    let dose = trialforge_core::math::round_dp(c.dose, 4);
    let arms = [
        ArmSpec {
            armcd: "PBO".into(),
            arm: "Placebo".into(),
            treatment: "PLACEBO".into(),
            dose: 0.0,
            ratio: c.ratio[0],
            change: 0.0,
            ae_share: 0.0,
        },
        ArmSpec {
            armcd: "DRUGA".into(),
            arm: format!("Drug A {dose} mg"),
            treatment: "DRUG A".into(),
            dose: c.dose,
            ratio: c.ratio[1],
            change: 0.0,
            ae_share: 0.0,
        },
    ];
    let (mut design, course) = build(
        opts,
        &Common {
            indication: "copd",
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMRES01",
            phase: "PHASE III TRIAL",
            control: "PLACEBO",
            title: None,
        },
        &arms,
    );
    design.schedule = Some(schedule(c.weeks));
    let pop = &mut design.population;
    pop.age_min = 40;
    pop.age_max = 85;
    pop.sex = Categorical::new([("F", 0.35), ("M", 0.65)]);
    pop.strata = Categorical::new([("1 EXACERBATION", 0.6), (">=2 EXACERBATIONS", 0.4)]);
    for arm in &mut design.arms {
        let t = &mut arm.treatment;
        t.form = "INJECTION".into();
        t.route = "SUBCUTANEOUS".into();
        t.frequency = "EVERY 4 WEEKS".into();
        t.cycle_days = 28;
    }
    let p = ExacerbationParams {
        pt: "Chronic obstructive pulmonary disease".into(),
        term: "COPD EXACERBATION".into(),
        rate_by_arm: vec![c.exacerbation_rate, c.exacerbation_rate * c.rate_ratio],
        stratum_rr: vec![1.0, c.prior_exacerbations_rr],
        dispersion: c.dispersion,
        severe_fraction: c.severe_fraction,
        median_days: 10.0,
        fev1_baseline: (c.fev1_baseline, 0.45),
        fev1_change_by_arm: c.fev1_change.to_vec(),
        fev1_decline: 0.04,
        onset_days: 14.0,
        treatment_days: (c.weeks * 7) as i32,
    };
    (design, course, p)
}
