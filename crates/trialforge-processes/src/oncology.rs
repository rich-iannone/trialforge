//! A randomised solid-tumour course with RECIST 1.1 tumour assessments, after
//! `clinical-trials-info/data-raw/sim-onco.R` (SIMONC01).
//!
//! **Latent survival.** Overall survival and the latent progression time have Weibull margins
//! with exact proportional hazards: arm and ECOG stratum multiply the hazard. A Gaussian copula
//! (correlation `os_pd_correlation`) links the two, so progression and death are dependent while
//! the OS hazard ratio stays a well-defined, exact estimand. Progression only happens if it
//! precedes death.
//!
//! **Response.** A subject is a latent responder with probability `response_prob[arm]`; the
//! response draw is copula-linked to the progression time, so responders progress later.
//! Responders are CR or PR; some non-responders have a single unconfirmed PR. The archetype
//! calibrates `response_prob` so the *confirmed* response rate hits its target.
//!
//! **Tumour assessments** (scans every `assessment_interval_days`, ±7 days) follow SIMONC01's
//! lesion model: 1–5 target lesions (at most 2 per organ, lymph nodes by short axis), 0–3
//! non-target lesions, shrinkage by response class, occasional non-evaluable lesions, and at
//! progression a target, new-lesion or non-target progression. Visit responses come from the
//! measurements by RECIST 1.1 ([`recist_visit`]), so RS is consistent with TR by construction.
//! Scans stop at progression, death, withdrawal or the data cutoff.
//!
//! Every subject draws the same latent uniforms in the same order whatever its arm (common
//! random numbers); tumour measurements use a separate stream.

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Epoch, Subject};

use crate::timeline::{
    DispositionEvent, DispositionKind, ElementSpan, Exposure, Lesion, LesionKind, Timeline,
    TumorAssessment, VisitOccurrence,
};

/// Days per month, as used for medians given in months.
pub const DAYS_PER_MONTH: f64 = 30.4375;

/// Parameters of the course (after calibration).
#[derive(Clone, Debug, PartialEq)]
pub struct OncologyParams {
    /// Weibull scale (days) of overall survival per arm, for a stratum hazard multiplier of 1.
    pub os_scale_days: Vec<f64>,
    /// Weibull shape of overall survival.
    pub os_shape: f64,
    /// Weibull scale (days) of the latent progression time per arm.
    pub pd_scale_days: Vec<f64>,
    /// Weibull shape of the latent progression time.
    pub pd_shape: f64,
    /// Hazard multiplier per stratification level (OS and progression).
    pub stratum_hr: Vec<f64>,
    /// Gaussian-copula correlation between OS and the progression time.
    pub os_pd_correlation: f64,
    /// Probability of being a latent responder, per arm.
    pub response_prob: Vec<f64>,
    /// Gaussian-copula correlation between response and the progression time.
    pub response_pd_correlation: f64,
    /// Share of responders who reach a complete response.
    pub cr_fraction: f64,
    /// Probability a non-responder shows a single, unconfirmed PR.
    pub unconfirmed_prob: f64,
    /// Probability of early withdrawal from the study.
    pub withdrawal_prob: f64,
    /// Window (days from randomisation) in which withdrawal happens.
    pub withdrawal_window: (f64, f64),
    /// Withdrawal reasons (`DSDECOD`) and their weights.
    pub withdrawal_reasons: Vec<(String, f64)>,
    /// Probability of stopping treatment for an adverse event.
    pub ae_stop_prob: f64,
    /// Window (days from first dose) in which AE-driven treatment stops happen.
    pub ae_stop_window: (f64, f64),
    /// Days between tumour assessments.
    pub assessment_interval_days: i32,
    /// Maximum number of treatment cycles.
    pub max_cycles: u32,
    /// For subjects alive and in follow-up at cutoff, last contact is up to this many days
    /// before the cutoff.
    pub last_contact_lag_days: i32,
    /// The response criteria of the tumour assessments.
    pub criteria: Criteria,
    /// Crossover on progression for arms with a switch treatment: (probability that a subject
    /// who progresses switches, hazard ratio of death after switching relative to the arm's).
    pub crossover: Option<(f64, f64)>,
    /// An interim data cutoff before the design's (a snapshot of the same study).
    pub data_cutoff: Option<Date>,
}

/// Tumour response criteria.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Criteria {
    /// RECIST 1.1 (solid tumours).
    Recist11,
    /// The Lugano 2014 classification (lymphoma), CT- and PET-CT-based.
    Lugano,
}

/// Why treatment stopped.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TreatmentStop {
    /// Still on treatment at the data cutoff.
    Ongoing,
    /// Progressive disease.
    Progression,
    /// Death.
    Death,
    /// Withdrawal from the study, with its `DSDECOD`.
    Withdrawal(String),
    /// Adverse event.
    AdverseEvent,
}

/// The latent draws of one subject (exposed for calibration).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Latent {
    /// Latent overall survival, days from randomisation.
    pub os_days: f64,
    /// Latent progression time, days from randomisation.
    pub pd_days: f64,
    /// Whether the subject is a latent responder.
    pub responder: bool,
    /// The response score Φ(z): the subject responds when it exceeds `1 - response_prob`.
    pub response_score: f64,
}

/// The latent survival and response draws for one subject.
pub fn latent(design: &Design, p: &OncologyParams, s: &Subject) -> (Latent, [f64; 8]) {
    let mut d = Draws::new(design.subject_stream(s.index).label("oncology"));
    // Draw order is part of the output contract; append new draws at the end.
    let z_os = d.normal(0.0, 1.0);
    let z_pd_own = d.normal(0.0, 1.0);
    let z_resp_own = d.normal(0.0, 1.0);
    let rest = [
        d.uniform(),
        d.uniform(),
        d.uniform(),
        d.uniform(),
        d.uniform(),
        d.uniform(),
        d.uniform(),
        d.uniform(),
    ];
    let rho = p.os_pd_correlation;
    let z_pd = rho * z_os + math::sqrt(1.0 - rho * rho) * z_pd_own;
    let rr = p.response_pd_correlation;
    let z_resp = rr * z_pd + math::sqrt(1.0 - rr * rr) * z_resp_own;
    let hr = p.stratum_hr[s.stratum];
    // S(t) = exp(-hr (t / scale)^shape) = Φ(-z): larger z → longer time.
    let time = |z: f64, scale: f64, shape: f64| {
        scale * math::pow(-math::ln(math::norm_cdf(-z)) / hr, 1.0 / shape)
    };
    let os_days = time(z_os, p.os_scale_days[s.arm], p.os_shape);
    let pd_days = time(z_pd, p.pd_scale_days[s.arm], p.pd_shape);
    let response_score = math::norm_cdf(z_resp);
    let responder = response_score > 1.0 - p.response_prob[s.arm];
    (
        Latent {
            os_days,
            pd_days,
            responder,
            response_score,
        },
        rest,
    )
}

/// RECIST 1.1 visit-level response. `ld`: target measurements (`None` = not evaluable);
/// `node`: lymph-node flags; `base_sum`, `nadir`: baseline sum and smallest sum so far;
/// `ntl`: non-target states; `new_les`: whether a new lesion is present.
/// Returns (target, non-target, overall) responses and the sum when every lesion was measured.
pub fn recist_visit(
    ld: &[Option<f64>],
    node: &[bool],
    base_sum: f64,
    nadir: f64,
    ntl: &[String],
    new_les: bool,
) -> (String, Option<String>, String, Option<f64>) {
    let all_measured = ld.iter().all(Option::is_some);
    let s: f64 = ld.iter().flatten().sum();
    // ≥ 20% and ≥ 5 mm over the nadir; a partial sum can still establish PD.
    let tl = if s >= 1.2 * nadir && s - nadir >= 5.0 {
        "PD"
    } else if !all_measured {
        "NE"
    } else if ld.iter().zip(node).all(|(v, &n)| {
        if n {
            v.unwrap_or(0.0) < 10.0
        } else {
            v.unwrap_or(0.0) == 0.0
        }
    }) {
        "CR"
    } else if s <= 0.7 * base_sum {
        "PR"
    } else {
        "SD"
    };
    let nt = if ntl.is_empty() {
        None
    } else if ntl.iter().any(|x| x == "UNEQUIVOCAL PROGRESSION") {
        Some("PD")
    } else if ntl.iter().all(|x| x == "ABSENT") {
        Some("CR")
    } else {
        Some("NON-CR/NON-PD")
    };
    let ovr = if new_les || tl == "PD" || nt == Some("PD") {
        "PD"
    } else if tl == "NE" {
        "NE"
    } else if tl == "CR" {
        if nt.is_none_or(|n| n == "CR") {
            "CR"
        } else {
            "PR"
        }
    } else {
        tl
    };
    (
        tl.to_string(),
        nt.map(str::to_string),
        ovr.to_string(),
        all_measured.then_some(s),
    )
}

/// Response class of a subject.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Class {
    Cr,
    Pr,
    Unconfirmed,
    Stable,
}

const TL_ORGANS: [&str; 6] = [
    "LUNG",
    "LIVER",
    "LYMPH NODE",
    "ADRENAL GLAND",
    "PERITONEUM",
    "KIDNEY",
];
const TL_WEIGHTS: [f64; 6] = [0.3, 0.25, 0.25, 0.08, 0.07, 0.05];
const NTL_ORGANS: [&str; 5] = ["BONE", "LUNG", "LIVER", "LYMPH NODE", "PLEURA"];
const NEW_ORGANS: [&str; 5] = ["LIVER", "LUNG", "BONE", "BRAIN", "LYMPH NODE"];
const PAIRED: [&str; 4] = ["LUNG", "ADRENAL GLAND", "KIDNEY", "PLEURA"];

/// One subject's tumour course: lesions and assessments, and when progression was seen.
struct TumorCourse {
    lesions: Vec<Lesion>,
    assessments: Vec<TumorAssessment>,
    progression: Option<Date>,
}

#[allow(clippy::too_many_arguments)]
fn tumor_course(
    design: &Design,
    p: &OncologyParams,
    s: &Subject,
    class: Class,
    onset: usize,
    depth: f64,
    sd_level: f64,
    pd_days: f64,
    end_day: f64,
) -> TumorCourse {
    let mut d = Draws::new(design.subject_stream(s.index).label("tumour"));
    let lat = |d: &mut Draws, loc: &str| {
        PAIRED
            .contains(&loc)
            .then(|| ["LEFT", "RIGHT"][d.below(2) as usize].to_string())
    };

    // Baseline lesions: at most 2 target lesions per organ.
    let n_tl = 1 + d.weighted(&[0.15, 0.25, 0.25, 0.2, 0.15]);
    let mut slots: Vec<usize> = (0..TL_ORGANS.len()).flat_map(|o| [o, o]).collect();
    let mut tl_loc: Vec<&str> = Vec::new();
    for _ in 0..n_tl {
        let w: Vec<f64> = slots.iter().map(|&o| TL_WEIGHTS[o]).collect();
        let k = d.weighted(&w);
        tl_loc.push(TL_ORGANS[slots.swap_remove(k)]);
    }
    let node: Vec<bool> = tl_loc.iter().map(|l| *l == "LYMPH NODE").collect();
    let base: Vec<f64> = node
        .iter()
        .map(|&n| {
            if n {
                d.int_range(15, 35) as f64
            } else {
                d.lognormal(24.0, 0.4).round().clamp(10.0, 90.0)
            }
        })
        .collect();
    let base_sum: f64 = base.iter().sum();
    let n_ntl = d.weighted(&[0.25, 0.35, 0.25, 0.15]);
    let ntl_loc: Vec<&str> = (0..n_ntl)
        .map(|_| NTL_ORGANS[d.below(NTL_ORGANS.len() as u64) as usize])
        .collect();
    let pd_type = match d.weighted(&[0.65, 0.2, 0.15]) {
        0 => "TARGET",
        2 if n_ntl > 0 => "NTL",
        _ => "NEW",
    };

    let screening = s.screening;
    let mut lesions = Vec::new();
    for (k, loc) in tl_loc.iter().enumerate() {
        lesions.push(Lesion {
            id: format!("T{:02}", k + 1),
            kind: LesionKind::Target,
            location: loc.to_string(),
            laterality: lat(&mut d, loc),
            node: node[k],
            visitnum: 1,
            date: screening,
        });
    }
    for (k, loc) in ntl_loc.iter().enumerate() {
        lesions.push(Lesion {
            id: format!("NT{:02}", k + 1),
            kind: LesionKind::NonTarget,
            location: loc.to_string(),
            laterality: lat(&mut d, loc),
            node: *loc == "LYMPH NODE",
            visitnum: 1,
            date: screening,
        });
    }
    let tl_ids: Vec<String> = (1..=n_tl).map(|k| format!("T{k:02}")).collect();
    let ntl_ids: Vec<String> = (1..=n_ntl).map(|k| format!("NT{k:02}")).collect();
    let mut assessments = vec![TumorAssessment {
        visitnum: 1,
        date: screening,
        target: tl_ids
            .iter()
            .cloned()
            .zip(base.iter().map(|&b| Some(b)))
            .collect(),
        sum: Some(base_sum),
        non_target: ntl_ids
            .iter()
            .map(|id| (id.clone(), "PRESENT".to_string()))
            .collect(),
        new_lesions: Vec::new(),
        response: None,
        perpendicular: Vec::new(),
        deauville: None,
    }];

    let interval = p.assessment_interval_days;
    let mut nadir = base_sum;
    let mut n_new = 0;
    let mut progression = None;
    for k in 1usize.. {
        let jitter = d.normal(0.0, 3.0).round().clamp(-7.0, 7.0) as i32;
        let t_k = interval * k as i32 + jitter;
        if f64::from(t_k) > end_day {
            break;
        }
        let pd_now = f64::from(t_k) >= pd_days;
        let mut cr_now = class == Class::Cr && k > onset;
        let pc = match class {
            Class::Cr => {
                if k < onset {
                    -0.2
                } else {
                    -0.55
                }
            }
            Class::Pr => {
                if k < onset {
                    (depth / 2.0).max(-0.25)
                } else {
                    depth
                }
            }
            Class::Unconfirmed => {
                if k < onset {
                    -0.15
                } else if k == onset {
                    -0.38
                } else {
                    -0.2
                }
            }
            Class::Stable => sd_level * (k as f64 / 2.0).min(1.0),
        } + d.normal(0.0, 0.03);
        let noise: Vec<f64> = (0..n_tl).map(|_| d.normal(0.0, 0.05)).collect();
        let node_cr: Vec<f64> = (0..n_tl).map(|_| d.int_range(4, 9) as f64).collect();
        let u_ne = d.uniform();
        let ne_pick = d.below(n_tl as u64) as usize;
        let pd_factor = d.uniform_range(1.3, 1.6);
        let ntl_pick = if n_ntl > 0 {
            d.below(n_ntl as u64) as usize
        } else {
            0
        };
        let new_loc = NEW_ORGANS[d.below(NEW_ORGANS.len() as u64) as usize];

        let mut new_les = false;
        let mut ntl_state: Vec<String> =
            vec![if cr_now { "ABSENT" } else { "PRESENT" }.to_string(); n_ntl];
        let mut design_sum = base_sum * (1.0 + pc);
        if pd_now && pd_type == "TARGET" {
            design_sum = (nadir * pd_factor).max(nadir + 6.0);
            cr_now = false;
        } else {
            design_sum = design_sum.min(nadir * 1.1); // no progression before the latent time
            if pd_now && pd_type == "NEW" {
                new_les = true;
            }
            if pd_now && pd_type == "NTL" {
                ntl_state[ntl_pick] = "UNEQUIVOCAL PROGRESSION".into();
            }
        }
        let mut ld: Vec<Option<f64>> = base
            .iter()
            .zip(&noise)
            .map(|(b, e)| Some((b * design_sum / base_sum * math::exp(*e)).round()))
            .collect();
        if cr_now {
            ld = node
                .iter()
                .zip(&node_cr)
                .map(|(&n, &v)| Some(if n { v } else { 0.0 }))
                .collect();
        }
        // Occasionally a target lesion cannot be evaluated (not at a progression).
        if !pd_now && u_ne < 0.05 {
            ld[ne_pick] = None;
        }
        let (tl, nt, ovr, sum) = recist_visit(&ld, &node, base_sum, nadir, &ntl_state, new_les);
        if let Some(sum) = sum {
            nadir = nadir.min(sum);
        }
        let visitnum = k as i64 + 1;
        let date = s.randomization.add_days(t_k);
        let mut new_ids = Vec::new();
        if new_les {
            n_new += 1;
            let id = format!("NEW{n_new:02}");
            lesions.push(Lesion {
                id: id.clone(),
                kind: LesionKind::New,
                location: new_loc.to_string(),
                laterality: lat(&mut d, new_loc),
                node: new_loc == "LYMPH NODE",
                visitnum,
                date,
            });
            new_ids.push(id);
        }
        let is_pd = ovr == "PD";
        assessments.push(TumorAssessment {
            visitnum,
            date,
            target: tl_ids.iter().cloned().zip(ld).collect(),
            sum,
            non_target: ntl_ids.iter().cloned().zip(ntl_state).collect(),
            new_lesions: new_ids,
            response: Some((tl, nt, ovr)),
            perpendicular: Vec::new(),
            deauville: None,
        });
        if is_pd {
            progression = Some(date);
            break;
        }
    }
    TumorCourse {
        lesions,
        assessments,
        progression,
    }
}

/// A target lesion's measurement at one Lugano assessment: (longest diameter, perpendicular),
/// in mm; `None` if not evaluable.
pub type Ppd = Option<(f64, f64)>;

/// Lugano 2014 visit-level responses: (CT-based, PET-CT-based) and the sum of the products of
/// the perpendicular diameters (SPD, mm²) when every target lesion was measured.
///
/// - CT: PD for a new lesion, unequivocal non-target progression, or a target lesion with a
///   longest diameter over 15 mm whose product grew by at least 50% and whose longest diameter
///   grew by at least 5 mm over that lesion's nadir; NE when a target lesion is not measured;
///   CR when every nodal target lesion is 15 mm or less, every extranodal one is gone and no
///   non-target lesion is present; PR for an SPD at most half the baseline; else SD.
/// - PET-CT: PMD on CT progression; NE without a Deauville score; CMR for a Deauville score of
///   3 or less; PMR for 4–5 with an SPD at most half the baseline (the residual mass of a
///   metabolic partial response); NE for 4–5 when the SPD is unknown; else SMD.
///
/// `nadir`: per target lesion, the smallest (product, longest diameter) so far.
pub fn lugano_visit(
    m: &[Ppd],
    node: &[bool],
    base_spd: f64,
    nadir: &[(f64, f64)],
    ntl: &[String],
    new_les: bool,
    deauville: Option<u8>,
) -> (String, String, Option<f64>) {
    let all_measured = m.iter().all(Option::is_some);
    let spd: f64 = m.iter().flatten().map(|(l, s)| l * s).sum();
    let lesion_pd = m.iter().zip(nadir).any(|(v, &(np, nl))| {
        v.is_some_and(|(l, s)| l > 15.0 && l * s >= 1.5 * np && l - nl >= 5.0)
    });
    let ntl_pd = ntl.iter().any(|x| x == "UNEQUIVOCAL PROGRESSION");
    let ct = if new_les || ntl_pd || lesion_pd {
        "PD"
    } else if !all_measured {
        "NE"
    } else if m.iter().zip(node).all(|(v, &n)| {
        let l = v.map_or(0.0, |x| x.0);
        if n { l <= 15.0 } else { l == 0.0 }
    }) && ntl.iter().all(|x| x == "ABSENT")
    {
        "CR"
    } else if spd <= 0.5 * base_spd {
        "PR"
    } else {
        "SD"
    };
    let pet = match (ct, deauville) {
        ("PD", _) => "PMD",
        (_, None) => "NE",
        (_, Some(d)) if d <= 3 => "CMR",
        ("NE", Some(_)) => "NE",
        (_, Some(_)) if spd <= 0.5 * base_spd => "PMR",
        _ => "SMD",
    };
    (ct.to_string(), pet.to_string(), all_measured.then_some(spd))
}

/// Nodal and extranodal sites of lymphoma lesions.
const LYMPHOMA_SITES: [&str; 4] = ["LYMPH NODE", "SPLEEN", "LIVER", "LUNG"];

/// One subject's lymphoma course under Lugano 2014: up to 6 target lesions (mostly lymph
/// nodes) measured by longest and perpendicular diameter, 0–3 non-target lesions, and an
/// FDG-PET Deauville score at every assessment, with the same response classes and latent
/// progression as the solid-tumour course.
#[allow(clippy::too_many_arguments)]
fn lymphoma_course(
    design: &Design,
    p: &OncologyParams,
    s: &Subject,
    class: Class,
    onset: usize,
    depth: f64,
    sd_level: f64,
    pd_days: f64,
    end_day: f64,
) -> TumorCourse {
    let mut d = Draws::new(design.subject_stream(s.index).label("lymphoma"));
    let lat = |d: &mut Draws, loc: &str| {
        matches!(loc, "LYMPH NODE" | "LUNG")
            .then(|| ["LEFT", "RIGHT"][d.below(2) as usize].to_string())
    };
    // Draw order is part of the output contract; append new draws at the end.
    let n_tl = 1 + d.weighted(&[0.15, 0.2, 0.25, 0.2, 0.1, 0.1]);
    let sites: Vec<&str> = (0..n_tl)
        .map(|_| LYMPHOMA_SITES[d.weighted(&[0.8, 0.08, 0.07, 0.05])])
        .collect();
    let node: Vec<bool> = sites.iter().map(|l| *l == "LYMPH NODE").collect();
    let base: Vec<(f64, f64)> = sites
        .iter()
        .map(|_| {
            let l = d.lognormal(32.0, 0.35).round().clamp(16.0, 100.0);
            let w = (l * d.uniform_range(0.5, 0.9)).round().max(10.0);
            (l, w)
        })
        .collect();
    let base_spd: f64 = base.iter().map(|(l, w)| l * w).sum();
    let n_ntl = d.weighted(&[0.4, 0.35, 0.15, 0.1]);
    let pd_type = match d.weighted(&[0.65, 0.2, 0.15]) {
        0 => "TARGET",
        2 if n_ntl > 0 => "NTL",
        _ => "NEW",
    };
    let screening = s.screening;
    let mut lesions = Vec::new();
    for (k, loc) in sites.iter().enumerate() {
        lesions.push(Lesion {
            id: format!("T{:02}", k + 1),
            kind: LesionKind::Target,
            location: loc.to_string(),
            laterality: lat(&mut d, loc),
            node: node[k],
            visitnum: 1,
            date: screening,
        });
    }
    for k in 0..n_ntl {
        lesions.push(Lesion {
            id: format!("NT{:02}", k + 1),
            kind: LesionKind::NonTarget,
            location: "LYMPH NODE".into(),
            laterality: lat(&mut d, "LYMPH NODE"),
            node: true,
            visitnum: 1,
            date: screening,
        });
    }
    let tl_ids: Vec<String> = (1..=n_tl).map(|k| format!("T{k:02}")).collect();
    let ntl_ids: Vec<String> = (1..=n_ntl).map(|k| format!("NT{k:02}")).collect();
    let base_deauville = 4 + u8::from(d.uniform() < 0.6);
    let mut assessments = vec![TumorAssessment {
        visitnum: 1,
        date: screening,
        target: tl_ids
            .iter()
            .cloned()
            .zip(base.iter().map(|b| Some(b.0)))
            .collect(),
        sum: Some(base_spd),
        non_target: ntl_ids
            .iter()
            .map(|id| (id.clone(), "PRESENT".to_string()))
            .collect(),
        new_lesions: Vec::new(),
        response: None,
        perpendicular: tl_ids
            .iter()
            .cloned()
            .zip(base.iter().map(|b| Some(b.1)))
            .collect(),
        deauville: Some(base_deauville),
    }];
    let mut nadir: Vec<(f64, f64)> = base.iter().map(|(l, w)| (l * w, *l)).collect();
    let interval = p.assessment_interval_days;
    let mut n_new = 0;
    let mut progression = None;
    for k in 1usize.. {
        let jitter = d.normal(0.0, 3.0).round().clamp(-7.0, 7.0) as i32;
        let t_k = interval * k as i32 + jitter;
        if f64::from(t_k) > end_day {
            break;
        }
        let pd_now = f64::from(t_k) >= pd_days;
        let cr_now = class == Class::Cr && k > onset && !pd_now;
        // Fractional change in each diameter (an area change of (1 + pc)²).
        let pc = match class {
            Class::Cr => {
                if k < onset {
                    -0.3
                } else {
                    -0.6
                }
            }
            Class::Pr => {
                if k < onset {
                    (depth / 2.0).max(-0.3)
                } else {
                    depth
                }
            }
            Class::Unconfirmed => {
                if k == onset {
                    -0.4
                } else {
                    -0.2
                }
            }
            Class::Stable => sd_level * (k as f64 / 2.0).min(1.0) / 2.0,
        } + d.normal(0.0, 0.03);
        let noise: Vec<f64> = (0..n_tl).map(|_| d.normal(0.0, 0.04)).collect();
        let node_cr: Vec<f64> = (0..n_tl).map(|_| d.int_range(6, 14) as f64).collect();
        let u_ne = d.uniform();
        let ne_pick = d.below(n_tl as u64) as usize;
        let pd_pick = d.below(n_tl as u64) as usize;
        let pd_factor = d.uniform_range(1.35, 1.6);
        let ntl_pick = if n_ntl > 0 {
            d.below(n_ntl as u64) as usize
        } else {
            0
        };
        let u_pet = d.uniform();
        let u_pet_done = d.uniform();
        let mut m: Vec<Ppd> = base
            .iter()
            .zip(&noise)
            .map(|((l, w), e)| {
                let f = ((1.0 + pc) * math::exp(*e)).max(0.05);
                Some(((l * f).round(), (w * f).round().max(1.0)))
            })
            .collect();
        if cr_now {
            m = node
                .iter()
                .zip(&node_cr)
                .map(|(&n, &v)| {
                    Some(if n {
                        (v, (v * 0.7).round())
                    } else {
                        (0.0, 0.0)
                    })
                })
                .collect();
        }
        let mut new_les = false;
        let mut ntl_state: Vec<String> =
            vec![if cr_now { "ABSENT" } else { "PRESENT" }.to_string(); n_ntl];
        if pd_now {
            match pd_type {
                "TARGET" => {
                    // Grow one lesion well past its nadir.
                    let (np, nl) = nadir[pd_pick];
                    let l = (nl.max(16.0) * pd_factor).round().max(nl + 6.0).max(17.0);
                    let w = ((np * pd_factor * pd_factor / l).round() + 1.0).max(10.0);
                    m[pd_pick] = Some((l, w));
                }
                "NTL" => ntl_state[ntl_pick] = "UNEQUIVOCAL PROGRESSION".into(),
                _ => new_les = true,
            }
        } else {
            // No lesion-level progression before the latent time: cap growth below the rule.
            for (v, &(np, nl)) in m.iter_mut().zip(&nadir) {
                if let Some((l, w)) = v
                    && (*l * *w >= 1.4 * np || *l - nl >= 5.0)
                {
                    let cap = (nl + 4.0).min(*l);
                    *w = (*w).min((1.4 * np / cap.max(1.0)).floor()).max(1.0);
                    *l = cap;
                }
            }
        }
        if !pd_now && u_ne < 0.04 {
            m[ne_pick] = None;
        }
        let deauville = (u_pet_done >= 0.03).then(|| match (pd_now, class) {
            (true, _) => 5,
            (_, Class::Cr) if k > onset => 1 + (u_pet * 3.0) as u8,
            (_, Class::Cr) | (_, Class::Pr) => 4,
            (_, Class::Unconfirmed) if k == onset => 4,
            _ => 4 + u8::from(u_pet < 0.5),
        });
        let (ct, pet, spd) =
            lugano_visit(&m, &node, base_spd, &nadir, &ntl_state, new_les, deauville);
        for (v, n) in m.iter().zip(nadir.iter_mut()) {
            if let Some((l, w)) = v {
                n.0 = n.0.min(l * w);
                n.1 = n.1.min(*l);
            }
        }
        let visitnum = k as i64 + 1;
        let date = s.randomization.add_days(t_k);
        let mut new_ids = Vec::new();
        if new_les {
            n_new += 1;
            let id = format!("NEW{n_new:02}");
            lesions.push(Lesion {
                id: id.clone(),
                kind: LesionKind::New,
                location: "LYMPH NODE".into(),
                laterality: lat(&mut d, "LYMPH NODE"),
                node: true,
                visitnum,
                date,
            });
            new_ids.push(id);
        }
        let is_pd = pet == "PMD";
        assessments.push(TumorAssessment {
            visitnum,
            date,
            target: tl_ids
                .iter()
                .cloned()
                .zip(m.iter().map(|v| v.map(|x| x.0)))
                .collect(),
            sum: spd,
            non_target: ntl_ids.iter().cloned().zip(ntl_state).collect(),
            new_lesions: new_ids,
            response: Some((ct, None, pet)),
            perpendicular: tl_ids
                .iter()
                .cloned()
                .zip(m.iter().map(|v| v.map(|x| x.1)))
                .collect(),
            deauville,
        });
        if is_pd {
            progression = Some(date);
            break;
        }
    }
    TumorCourse {
        lesions,
        assessments,
        progression,
    }
}

/// Whether a list of responses (date, overall) shows an objective response: confirmed CR/PR
/// under RECIST 1.1 ([`confirmed_response`]), any CMR or PMR under Lugano (no confirmation).
pub fn responded(criteria: Criteria, responses: &[(Date, &str)]) -> bool {
    match criteria {
        Criteria::Recist11 => confirmed_response(responses),
        Criteria::Lugano => responses.iter().any(|(_, r)| matches!(*r, "CMR" | "PMR")),
    }
}

/// Simulates every subject's course as a [`Timeline`].
pub fn simulate(design: &Design, params: &OncologyParams, subjects: &[Subject]) -> Vec<Timeline> {
    subjects
        .iter()
        .map(|s| simulate_one(design, params, s))
        .collect()
}

/// One subject's course.
pub fn simulate_one(design: &Design, p: &OncologyParams, s: &Subject) -> Timeline {
    let (lat, u) = latent(design, p, s);
    let [
        u_class,
        u_onset,
        u_depth,
        u_sd,
        u_withdraw,
        u_wd_time,
        u_ae,
        u_ae_time,
    ] = u;
    let mut d = Draws::new(design.subject_stream(s.index).label("oncology-disposition"));
    let reason_w: Vec<f64> = p.withdrawal_reasons.iter().map(|(_, w)| *w).collect();
    let reason = p.withdrawal_reasons[d.weighted(&reason_w)].0.clone();
    let contact_lag = d.int_range(0, i64::from(p.last_contact_lag_days)) as i32;

    let class = if lat.responder {
        if u_class < p.cr_fraction {
            Class::Cr
        } else {
            Class::Pr
        }
    } else if u_class < p.unconfirmed_prob {
        Class::Unconfirmed
    } else {
        Class::Stable
    };
    let onset = 1 + usize::from(u_onset >= 0.6);
    let depth = -0.75 + u_depth * (0.75 - 0.36);
    let sd_level = -0.25 + u_sd * 0.37;

    let rand = s.randomization;
    let first_dose_day = s.treatment_start.days_since(rand);
    let day = |t: f64| rand.add_days(t.floor().min(1e6) as i32);
    let cutoff = p.data_cutoff.unwrap_or(design.cutoff);
    let death_latent = day(lat.os_days.max(f64::from(first_dose_day + 1)));
    let withdrawal_latent = (u_withdraw < p.withdrawal_prob).then(|| {
        day(p.withdrawal_window.0 + u_wd_time * (p.withdrawal_window.1 - p.withdrawal_window.0))
    });
    let mut death = (death_latent <= cutoff && withdrawal_latent.is_none_or(|w| death_latent < w))
        .then_some(death_latent);
    let mut withdrawal =
        withdrawal_latent.filter(|&w| w <= cutoff && death.is_none_or(|dd| w < dd));

    // Scans stop at death, withdrawal or cutoff; progression only if it precedes death.
    let end = [Some(death_latent), withdrawal_latent, Some(cutoff)]
        .into_iter()
        .flatten()
        .min()
        .expect("cutoff");
    let end_day = f64::from(end.days_since(rand));
    let pd_days = if lat.pd_days < lat.os_days {
        lat.pd_days
    } else {
        f64::INFINITY
    };
    let tumor = match p.criteria {
        Criteria::Recist11 => tumor_course(
            design, p, s, class, onset, depth, sd_level, pd_days, end_day,
        ),
        Criteria::Lugano => lymphoma_course(
            design, p, s, class, onset, depth, sd_level, pd_days, end_day,
        ),
    };
    let progression = tumor.progression;

    // Crossover on progression: switching changes the hazard of death from then on (Weibull
    // proportional hazards on the same latent uniform), so the death time moves.
    let mut crossover: Option<Date> = None;
    if let (Some((prob, hr)), Some(pd), Some(_)) = (
        p.crossover,
        progression,
        design.arms[s.arm].switch_to.as_ref(),
    ) {
        let mut x = Draws::new(design.subject_stream(s.index).label("crossover"));
        let (u_switch, delay) = (x.uniform(), x.int_range(7, 28) as i32);
        let date = pd.add_days(delay);
        let alive = date < death_latent && withdrawal_latent.is_none_or(|w| date < w);
        if u_switch < prob && alive && date <= cutoff {
            let k = p.os_shape;
            let tx = f64::from(date.days_since(rand));
            let os = math::pow(
                math::pow(tx, k) + (math::pow(lat.os_days, k) - math::pow(tx, k)) / hr,
                1.0 / k,
            );
            let latent_death = day(os).max(date.add_days(1));
            death = (latent_death <= cutoff && withdrawal_latent.is_none_or(|w| latent_death < w))
                .then_some(latent_death);
            withdrawal =
                withdrawal_latent.filter(|&w| w <= cutoff && death.is_none_or(|dd| w < dd));
            crossover = Some(date);
        }
    }

    let ae_latent = (u_ae < p.ae_stop_prob).then(|| {
        s.treatment_start.add_days(
            (p.ae_stop_window.0 + u_ae_time * (p.ae_stop_window.1 - p.ae_stop_window.0)) as i32,
        )
    });
    let candidates = [
        (progression, TreatmentStop::Progression),
        (death, TreatmentStop::Death),
        (withdrawal, TreatmentStop::Withdrawal(reason.clone())),
        (ae_latent, TreatmentStop::AdverseEvent),
    ];
    let (treatment_stop, stop_reason) = candidates
        .into_iter()
        .filter_map(|(date, why)| date.filter(|&dd| dd <= cutoff).map(|dd| (dd, why)))
        .min_by_key(|(dd, _)| *dd)
        .unwrap_or((cutoff, TreatmentStop::Ongoing));

    let cycle = design.arms[s.arm].treatment.cycle_days;
    let doses: Vec<Date> = (0..p.max_cycles as i32)
        .map(|k| s.treatment_start.add_days(k * cycle))
        .take_while(|&dd| dd == s.treatment_start || dd < treatment_stop)
        .collect();
    let last_dose = *doses.last().expect("the first dose is always given");
    // Switch-treatment doses from the crossover until death, withdrawal or the cutoff.
    let switch_doses: Vec<Date> = crossover.map_or_else(Vec::new, |x| {
        let stop = [death, withdrawal, Some(cutoff.add_days(1))]
            .into_iter()
            .flatten()
            .min()
            .expect("cutoff");
        let cycle = design.arms[s.arm]
            .switch_to
            .as_ref()
            .map_or(21, |t| t.0.cycle_days);
        (0..p.max_cycles as i32)
            .map(|k| x.add_days(k * cycle))
            .take_while(|&dd| dd == x || dd < stop)
            .collect()
    });
    let last_dose = switch_doses.last().copied().unwrap_or(last_dose);
    let last_scan = tumor.assessments.last().map_or(s.screening, |a| a.date);
    let stopped = stop_reason != TreatmentStop::Ongoing;
    let last_contact = death
        .or(withdrawal)
        .unwrap_or_else(|| cutoff.add_days(-contact_lag))
        .max(last_dose)
        .max(last_scan)
        .max(if stopped {
            treatment_stop
        } else {
            s.treatment_start
        });

    let dose = design.arms[s.arm].treatment.dose;
    let switch_dose = design.arms[s.arm]
        .switch_to
        .as_ref()
        .map_or(dose, |t| t.0.dose);
    let exposure = doses
        .iter()
        .map(|&dd| (dd, false))
        .chain(switch_doses.iter().map(|&dd| (dd, true)))
        .map(|(dd, switched)| Exposure {
            start: dd,
            end: dd,
            dose: if switched { switch_dose } else { dose },
            time: None,
            period: None,
            location: None,
            group: None,
            visitnum: None,
            switched,
        })
        .collect();
    let first_dose = s.treatment_start;
    let treatment_end = if stopped {
        treatment_stop
    } else {
        last_contact
    };
    let mut elements = vec![
        ElementSpan {
            epoch: Epoch::Screening,
            start: s.consent,
            end: first_dose,
        },
        ElementSpan {
            epoch: Epoch::Treatment,
            start: first_dose,
            end: treatment_end,
        },
    ];
    match crossover {
        Some(x) => {
            if x > treatment_end {
                elements.push(ElementSpan {
                    epoch: Epoch::FollowUp,
                    start: treatment_end,
                    end: x,
                });
            }
            elements.push(ElementSpan {
                epoch: Epoch::Continuation,
                start: x,
                end: last_contact,
            });
        }
        None if last_contact > treatment_end => elements.push(ElementSpan {
            epoch: Epoch::FollowUp,
            start: treatment_end,
            end: last_contact,
        }),
        None => {}
    }
    let mut disposition = vec![
        DispositionEvent {
            decod: "INFORMED CONSENT OBTAINED".into(),
            kind: DispositionKind::Milestone,
            date: s.consent,
        },
        DispositionEvent {
            decod: "RANDOMIZED".into(),
            kind: DispositionKind::Milestone,
            date: rand,
        },
    ];
    let treatment_term = match &stop_reason {
        TreatmentStop::Ongoing => None,
        TreatmentStop::Progression => Some("PROGRESSIVE DISEASE".to_string()),
        TreatmentStop::Death => Some("DEATH".to_string()),
        TreatmentStop::Withdrawal(r) => Some(r.clone()),
        TreatmentStop::AdverseEvent => Some("ADVERSE EVENT".to_string()),
    };
    if let Some(decod) = treatment_term {
        disposition.push(DispositionEvent {
            decod,
            kind: DispositionKind::EndOfTreatment,
            date: treatment_stop,
        });
    }
    match (death, withdrawal) {
        (Some(dd), _) => disposition.push(DispositionEvent {
            decod: "DEATH".into(),
            kind: DispositionKind::EndOfStudy,
            date: dd,
        }),
        (None, Some(w)) => disposition.push(DispositionEvent {
            decod: reason,
            kind: DispositionKind::EndOfStudy,
            date: w,
        }),
        (None, None) => {}
    }
    let visits = tumor
        .assessments
        .iter()
        .map(|a| VisitOccurrence {
            visitnum: a.visitnum,
            date: a.date,
            end: None,
        })
        .collect();

    Timeline {
        exposure,
        elements,
        visits,
        disposition,
        death,
        last_contact,
        treatment_end,
        adverse_events: Vec::new(),
        medications: Vec::new(),
        history: Vec::new(),
        findings: Vec::new(),
        lesions: tumor.lesions,
        tumor: tumor.assessments,
        pk: Vec::new(),
        clinical_events: Vec::new(),
        study_eye: None,
        birth_date: None,
        diaries: Vec::new(),
    }
}

/// Whether a subject's overall responses contain a confirmed response (CR or PR confirmed by a
/// CR/PR at least 28 days later, with at most one NE in between and no PR after a CR), as
/// admiralonco's `crsp_y_cr` / `crsp_y_pr`. Used to calibrate response rates.
pub fn confirmed_response(responses: &[(Date, &str)]) -> bool {
    for (i, &(di, ri)) in responses.iter().enumerate() {
        if ri != "CR" && ri != "PR" {
            continue;
        }
        // The joined records run from i+1 up to the first qualifying confirmation.
        let Some(j) = responses[i + 1..].iter().position(|&(dj, rj)| {
            (if ri == "CR" {
                rj == "CR"
            } else {
                rj == "CR" || rj == "PR"
            }) && dj.days_since(di) >= 28
        }) else {
            continue;
        };
        let joined = &responses[i + 1..=i + 1 + j];
        let allowed = |r: &str| {
            if ri == "CR" {
                r == "CR" || r == "NE"
            } else {
                matches!(r, "CR" | "PR" | "NE")
            }
        };
        let ne = joined.iter().filter(|(_, r)| *r == "NE").count();
        let ok = joined.iter().all(|(_, r)| allowed(r)) && ne <= 1;
        let pr_after_cr = ri == "PR" && {
            let first_cr = joined.iter().find(|(_, r)| *r == "CR").map(|(d, _)| *d);
            let last_pr = joined
                .iter()
                .rev()
                .find(|(_, r)| *r == "PR")
                .map(|(d, _)| *d);
            matches!((first_cr, last_pr), (Some(c), Some(p)) if c <= p)
        };
        if ok && !pr_after_cr {
            return true;
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ntl(v: &[&str]) -> Vec<String> {
        v.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn recist_rules() {
        let node = [false, true];
        // PR: 30% down from baseline 100.
        let r = recist_visit(
            &[Some(50.0), Some(15.0)],
            &node,
            100.0,
            100.0,
            &ntl(&["PRESENT"]),
            false,
        );
        assert_eq!((r.0.as_str(), r.2.as_str()), ("PR", "PR"));
        // CR needs non-nodes at 0 and nodes < 10 mm; NTL present downgrades to PR.
        let r = recist_visit(
            &[Some(0.0), Some(8.0)],
            &node,
            100.0,
            60.0,
            &ntl(&["PRESENT"]),
            false,
        );
        assert_eq!((r.0.as_str(), r.2.as_str()), ("CR", "PR"));
        let r = recist_visit(
            &[Some(0.0), Some(8.0)],
            &node,
            100.0,
            60.0,
            &ntl(&["ABSENT"]),
            false,
        );
        assert_eq!(r.2, "CR");
        // PD: +20% and +5 mm over the nadir, even with a lesion missing.
        let r = recist_visit(&[Some(70.0), None], &node, 100.0, 55.0, &[], false);
        assert_eq!((r.0.as_str(), r.3), ("PD", None));
        // NE when not all measured and no PD.
        let r = recist_visit(&[Some(50.0), None], &node, 100.0, 60.0, &[], false);
        assert_eq!(r.2, "NE");
        // New lesion → PD.
        assert_eq!(
            recist_visit(&[Some(90.0), Some(20.0)], &node, 100.0, 100.0, &[], true).2,
            "PD"
        );
    }

    #[test]
    fn confirmation_rules() {
        let d = |n: i32| Date::from_days(n);
        assert!(confirmed_response(&[(d(0), "PR"), (d(56), "PR")]));
        assert!(!confirmed_response(&[(d(0), "PR"), (d(20), "PR")]));
        assert!(confirmed_response(&[
            (d(0), "PR"),
            (d(20), "NE"),
            (d(56), "PR")
        ]));
        assert!(!confirmed_response(&[
            (d(0), "PR"),
            (d(20), "NE"),
            (d(30), "NE"),
            (d(56), "PR")
        ]));
        assert!(!confirmed_response(&[(d(0), "PR"), (d(56), "SD")]));
        assert!(confirmed_response(&[(d(0), "CR"), (d(56), "CR")]));
        assert!(!confirmed_response(&[(d(0), "CR"), (d(56), "PR")]));
        assert!(confirmed_response(&[(d(0), "PR"), (d(56), "CR")]));
        // A PR after a CR inside the window does not confirm the first PR.
        assert!(!confirmed_response(&[
            (d(0), "PR"),
            (d(10), "CR"),
            (d(40), "PR")
        ]));
    }
}
