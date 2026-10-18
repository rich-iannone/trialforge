//! `ophthalmology`: a randomised, double-masked, active-controlled non-inferiority study of
//! intravitreal Drug A in neovascular age-related macular degeneration. One eye per subject is
//! the study eye; both eyes are examined every 4 weeks (best corrected visual acuity in ETDRS
//! letters and central subfield thickness by SD-OCT, in OE), and the study eye is injected
//! monthly for a loading phase and then every few weeks.
//!
//! Ophthalmic findings come from `trialforge_processes::ophtha`; the controls are in ETDRS
//! letters and micrometres.

use serde::{Deserialize, Serialize};
use trialforge_model::{Design, Epoch, PlannedVisit, Schedule, VisitKind};
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::ophtha::OphthaParams;
use trialforge_processes::reference::AeBackground;

use crate::parallel::{ArmSpec, Common, build};
use crate::{ControlDoc, StudyOptions};

/// Controls for `ophthalmology`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Treatment duration in weeks (a multiple of 4; visits every 4 weeks).
    pub weeks: u32,
    /// Drug A dose (mg per injection).
    pub dose: f64,
    /// Comparator dose (mg per injection).
    pub comparator_dose: f64,
    /// Monthly loading injections before the maintenance interval.
    pub loading_doses: u32,
    /// Maintenance injection interval in weeks: [comparator, Drug A].
    pub interval_weeks: [u32; 2],
    /// Mean study-eye BCVA at baseline (ETDRS letters).
    pub bcva_baseline: f64,
    /// Mean change in study-eye BCVA at the end of treatment: [comparator, Drug A] (letters).
    pub bcva_change: [f64; 2],
    /// SD of subjects' BCVA change around their arm's mean (letters).
    pub response_sd: f64,
    /// Mean change in study-eye central subfield thickness at the end of treatment:
    /// [comparator, Drug A] (µm).
    pub cst_change: [f64; 2],
    /// Non-inferiority margin for the BCVA change, Drug A minus comparator (letters; higher is
    /// better), or `null`.
    pub ni_margin: Option<f64>,
    /// Allocation ratio: [comparator, Drug A].
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
            weeks: 48,
            dose: 6.0,
            comparator_dose: 2.0,
            loading_doses: 3,
            interval_weeks: [8, 8],
            bcva_baseline: 55.0,
            bcva_change: [8.0, 8.0],
            response_sd: 9.0,
            cst_change: [-140.0, -160.0],
            ni_margin: Some(4.0),
            ratio: [1, 1],
            dropout_rate: 0.1,
            ae_rate_multiplier: 1.0,
            missed_visit: 0.03,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "weeks",
            "Treatment duration in weeks (a multiple of 4; visits every 4 weeks)",
        ),
        ControlDoc::new("dose", "Drug A dose in mg per intravitreal injection"),
        ControlDoc::new(
            "comparator_dose",
            "Comparator dose in mg per intravitreal injection",
        ),
        ControlDoc::new(
            "loading_doses",
            "Monthly loading injections before the maintenance interval",
        ),
        ControlDoc::new(
            "interval_weeks",
            "Maintenance injection interval in weeks (4, 8, 12 or 16): [comparator, Drug A]",
        ),
        ControlDoc::new(
            "bcva_baseline",
            "Mean study-eye best corrected visual acuity at baseline (ETDRS letters)",
        ),
        ControlDoc::new(
            "bcva_change",
            "Mean change in study-eye BCVA at the end of treatment: [comparator, Drug A] (letters)",
        ),
        ControlDoc::new(
            "response_sd",
            "SD of subjects' BCVA change around their arm's mean (letters)",
        ),
        ControlDoc::new(
            "cst_change",
            "Mean change in study-eye central subfield thickness at the end of treatment: [comparator, Drug A] (um)",
        ),
        ControlDoc::new(
            "ni_margin",
            "Non-inferiority margin for the BCVA change, Drug A minus comparator (letters); null for none",
        ),
        ControlDoc::new("ratio", "Allocation ratio: [comparator, Drug A]"),
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
    if !(c.dose.is_finite()
        && c.dose > 0.0
        && c.comparator_dose.is_finite()
        && c.comparator_dose > 0.0)
    {
        return Err("`dose` and `comparator_dose` must be positive".into());
    }
    if !(1..=6).contains(&c.loading_doses) || c.loading_doses * 4 > c.weeks {
        return Err(
            "`loading_doses` must be between 1 and 6 and fit in the treatment period".into(),
        );
    }
    if !c.interval_weeks.iter().all(|w| [4, 8, 12, 16].contains(w)) {
        return Err("`interval_weeks` entries must be 4, 8, 12 or 16".into());
    }
    if !(30.0..=75.0).contains(&c.bcva_baseline) {
        return Err("`bcva_baseline` must be between 30 and 75 letters".into());
    }
    if !c
        .bcva_change
        .iter()
        .all(|x| x.is_finite() && x.abs() <= 30.0)
    {
        return Err("`bcva_change` entries must be within ±30 letters".into());
    }
    if !(0.0..=25.0).contains(&c.response_sd) {
        return Err("`response_sd` must be between 0 and 25 letters".into());
    }
    if !c
        .cst_change
        .iter()
        .all(|x| x.is_finite() && (-300.0..=100.0).contains(x))
    {
        return Err("`cst_change` entries must be between -300 and 100 um".into());
    }
    if c.ni_margin.is_some_and(|m| !(m.is_finite() && m > 0.0)) {
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

/// Screening, baseline, then a visit every 4 weeks to the end of treatment.
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
    for (k, w) in (4..=weeks).step_by(4).enumerate() {
        let last = w == weeks;
        visits.push(visit(
            3 + k as i64,
            format!("WEEK {w}"),
            Some((w * 7 + 1) as i32),
            (3, 3),
            Epoch::Treatment,
            if last {
                VisitKind::EndOfTreatment
            } else {
                VisitKind::Treatment
            },
            last || w == weeks / 2,
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

/// Planned injection study days: monthly loading doses, then every `interval` weeks, before the
/// end-of-treatment visit.
pub(crate) fn injection_days(c: &Controls, interval: u32) -> Vec<i32> {
    let mut weeks: Vec<u32> = (0..c.loading_doses).map(|k| 4 * k).collect();
    let mut w = 4 * (c.loading_doses - 1) + interval;
    while w < c.weeks {
        weeks.push(w);
        w += interval;
    }
    weeks.into_iter().map(|w| (w * 7 + 1) as i32).collect()
}

fn fmt_dose(d: f64) -> String {
    let s = format!("{}", trialforge_core::math::round_dp(d, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, CourseParams, OphthaParams) {
    let arms = [
        ArmSpec {
            armcd: "ACTCTRL".into(),
            arm: format!("Comparator {} mg", fmt_dose(c.comparator_dose)),
            treatment: "COMPARATOR".into(),
            dose: c.comparator_dose,
            ratio: c.ratio[0],
            change: 0.0,
            ae_share: 0.0,
        },
        ArmSpec {
            armcd: "DRUGA".into(),
            arm: format!("Drug A {} mg", fmt_dose(c.dose)),
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
            indication: "neovascular_amd",
            weeks: c.weeks,
            dropout_rate: c.dropout_rate,
            ae_rate_multiplier: c.ae_rate_multiplier,
            missed_visit: c.missed_visit,
            studyid: "SIMOPH01",
            phase: "PHASE III TRIAL",
            control: "ACTIVE",
            title: None,
        },
        &arms,
    );
    design.schedule = Some(schedule(c.weeks));
    for (arm, weeks) in design.arms.iter_mut().zip(c.interval_weeks) {
        let t = &mut arm.treatment;
        t.form = "INJECTION".into();
        t.route = "INTRAVITREAL".into();
        t.frequency = format!("EVERY {weeks} WEEKS");
        t.cycle_days = (weeks * 7) as i32;
    }
    let ae = |pt: &str, rate: f64, days: f64, serious: f64, severity: [f64; 3]| AeBackground {
        pt: pt.into(),
        rate_per_year: rate,
        median_days: days,
        p_serious: serious,
        severity,
    };
    // Injection-related ocular events (both arms).
    course.extra_ae = vec![
        ae(
            "Conjunctival haemorrhage",
            0.15,
            10.0,
            0.0,
            [0.85, 0.14, 0.01],
        ),
        ae("Eye pain", 0.10, 2.0, 0.0, [0.70, 0.27, 0.03]),
        ae(
            "Intraocular pressure increased",
            0.06,
            3.0,
            0.0,
            [0.75, 0.22, 0.03],
        ),
        ae("Vitreous floaters", 0.06, 30.0, 0.0, [0.80, 0.19, 0.01]),
        ae("Endophthalmitis", 0.003, 21.0, 1.0, [0.0, 0.3, 0.7]),
    ];
    let ophtha = OphthaParams {
        bcva_baseline: (c.bcva_baseline, 11.0),
        bcva_inclusion: (24.0, 78.0),
        bcva_change_by_arm: c.bcva_change.to_vec(),
        response_sd: c.response_sd,
        bcva_noise_sd: 4.0,
        cst_change_by_arm: c.cst_change.to_vec(),
        onset_days: 21.0,
        treatment_days: (c.weeks * 7) as i32,
        injection_days: c
            .interval_weeks
            .iter()
            .map(|&w| injection_days(c, w))
            .collect(),
    };
    (design, course, ophtha)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn injections_load_monthly_then_follow_the_interval() {
        let c = Controls::default();
        assert_eq!(injection_days(&c, 8), vec![1, 29, 57, 113, 169, 225, 281]);
        assert_eq!(injection_days(&c, 4).len(), 12);
    }
}
