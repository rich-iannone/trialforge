//! Calibration of the oncology course to clinical targets.

use trialforge_core::math;
use trialforge_model::{Design, enroll};
use trialforge_processes::oncology::{OncologyParams, latent, responded, simulate_one};

/// Clinical targets for the oncology course.
#[derive(Clone, Debug, PartialEq)]
pub struct OncologyTargets {
    /// Population median overall survival of the reference (first control) arm, in days.
    pub os_median_days: f64,
    /// OS hazard ratio per arm versus the reference arm (reference = 1).
    pub os_hr: Vec<f64>,
    /// Population median PFS (latent: min(progression, death)) of the reference arm, in days.
    pub pfs_median_days: f64,
    /// Progression hazard ratio per arm versus the reference arm.
    pub pd_hr: Vec<f64>,
    /// Confirmed objective response rate per arm.
    pub orr: Vec<f64>,
    /// Index of the reference arm.
    pub reference_arm: usize,
}

/// What calibration achieved on the calibration cohort.
#[derive(Clone, Debug, PartialEq)]
pub struct Calibration {
    /// The calibrated parameters.
    pub params: OncologyParams,
    /// Confirmed ORR per arm achieved on the calibration cohort.
    pub achieved_orr: Vec<f64>,
    /// Whether every ORR target was reachable (a target above what the schedule allows, e.g.
    /// because subjects die before confirmation, is clamped).
    pub orr_reachable: bool,
}

/// Size of the calibration cohort (per study; split across arms by the allocation ratio). The
/// Monte Carlo error of a calibrated rate is about sqrt(p(1-p)/(COHORT/arms)): 0.3 points for
/// an ORR of 20% with two arms.
pub const COHORT: usize = 40_000;

/// Calibrates `base` (whose structural parameters are kept) to `targets` for `design`.
pub fn calibrate(design: &Design, base: &OncologyParams, targets: &OncologyTargets) -> Calibration {
    let mut p = base.clone();
    let weights: Vec<f64> = design
        .population
        .strata
        .levels
        .iter()
        .map(|(_, w)| *w)
        .collect();
    let r = targets.reference_arm;

    // OS: exact Weibull scale for the reference arm; other arms by proportional hazards.
    let os_ref = crate::median_scale(targets.os_median_days, p.os_shape, &weights, &p.stratum_hr);
    p.os_scale_days = targets
        .os_hr
        .iter()
        .map(|hr| os_ref * math::pow(*hr, -1.0 / p.os_shape))
        .collect();

    // A synthetic cohort from the same design (its own seed stream), for Monte Carlo.
    let mut cohort_design = design.clone();
    cohort_design.n = COHORT;
    cohort_design.seed = design.stream().label("calibration").seed() >> 11;
    let cohort = enroll(&cohort_design).expect("the design was validated");

    // PFS: bisect the reference arm's progression scale on the latent median of min(PD, death).
    // Progression times scale linearly with the Weibull scale, so draw them once at scale 1.
    let mut unit = p.clone();
    unit.pd_scale_days = vec![1.0; design.arms.len()];
    let reference: Vec<(f64, f64)> = cohort
        .iter()
        .filter(|s| s.arm == r)
        .map(|s| {
            let (l, _) = latent(&cohort_design, &unit, s);
            (l.pd_days, l.os_days)
        })
        .collect();
    let pfs_median = |scale: f64| -> f64 {
        let mut t: Vec<f64> = reference
            .iter()
            .map(|(pd, os)| (pd * scale).min(*os))
            .collect();
        t.sort_by(f64::total_cmp);
        t[t.len() / 2]
    };
    let pd_ref = crate::bisect(pfs_median, 1.0, 1e5, targets.pfs_median_days, 60);
    p.pd_scale_days = targets
        .pd_hr
        .iter()
        .map(|hr| pd_ref * math::pow(*hr, -1.0 / p.pd_shape))
        .collect();

    // Confirmed ORR: a subject responds when its latent score exceeds 1 - response_prob, and
    // nothing else in its course depends on response_prob. So simulate each subject once as a
    // responder and once not, then bisect the probability on those outcomes.
    let confirmed = |q: &OncologyParams, s: &trialforge_model::Subject| {
        let t = simulate_one(&cohort_design, q, s);
        let resp: Vec<_> = t
            .tumor
            .iter()
            .filter_map(|a| a.response.as_ref().map(|r| (a.date, r.2.as_str())))
            .collect();
        responded(q.criteria, &resp)
    };
    let mut achieved = Vec::new();
    let mut reachable = true;
    for (arm, &target) in targets.orr.iter().enumerate() {
        let (mut always, mut never) = (p.clone(), p.clone());
        always.response_prob[arm] = 1.0;
        never.response_prob[arm] = 0.0;
        // (score, confirmed as a responder, confirmed as a non-responder)
        let outcomes: Vec<(f64, bool, bool)> = cohort
            .iter()
            .filter(|s| s.arm == arm)
            .map(|s| {
                let score = latent(&cohort_design, &p, s).0.response_score;
                (score, confirmed(&always, s), confirmed(&never, s))
            })
            .collect();
        let orr = |prob: f64| -> f64 {
            let k = outcomes
                .iter()
                .filter(|(score, yes, no)| if *score > 1.0 - prob { *yes } else { *no })
                .count();
            k as f64 / outcomes.len().max(1) as f64
        };
        let ceiling = orr(1.0);
        let prob = if target >= ceiling {
            reachable &= target <= ceiling + 1e-9;
            1.0
        } else {
            crate::bisect(orr, 0.0, 1.0, target, 40)
        };
        achieved.push(orr(prob));
        p.response_prob[arm] = prob;
    }
    Calibration {
        params: p,
        achieved_orr: achieved,
        orr_reachable: reachable,
    }
}
