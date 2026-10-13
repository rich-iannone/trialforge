//! `oncology_solid_tumor`: SIMONC01's design with clinical-target controls. 1:1 Drug A vs
//! Control, permuted blocks of 4 stratified by ECOG (0/1), 18 months of enrollment from
//! 2023-03-01, data cutoff 30 months after the first randomisation, Q3W dosing and RECIST 1.1
//! tumour assessments every 8 weeks.
//!
//! Targets are population quantities: OS and PFS medians of the control arm, hazard ratios for
//! Drug A (OS exactly proportional; the progression hazard ratio applies to the latent
//! progression time), and confirmed ORR per arm. They are calibrated to latent parameters once
//! per control set (independently of the study's seed and size) and cached.

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};

use serde::{Deserialize, Serialize};
use trialforge_calibrate::oncology::{Calibration, OncologyTargets, calibrate};
use trialforge_core::Date;
use trialforge_model::{
    Arm, Categorical, Design, Epoch, PlannedVisit, Population, Schedule, Site, Timing, Treatment,
    TrialInfo, VisitKind,
};
use trialforge_processes::oncology::{DAYS_PER_MONTH, OncologyParams};

use crate::{ControlDoc, StudyOptions};

/// Controls for `oncology_solid_tumor`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Median overall survival of the Control arm, in months.
    pub os_median_months: f64,
    /// OS hazard ratio, Drug A vs Control.
    pub os_hr: f64,
    /// Median progression-free survival of the Control arm, in months.
    pub pfs_median_months: f64,
    /// Progression hazard ratio, Drug A vs Control.
    pub pfs_hr: f64,
    /// Confirmed objective response rate: [Drug A, Control].
    pub orr: [f64; 2],
    /// Weibull shape of overall survival (1 = exponential).
    pub os_shape: f64,
    /// Weibull shape of the progression time.
    pub pfs_shape: f64,
    /// Correlation between progression and death (Gaussian copula).
    pub os_pfs_correlation: f64,
    /// Probability of early withdrawal.
    pub withdrawal_prob: f64,
    /// Probability of stopping treatment for an adverse event.
    pub ae_stop_prob: f64,
    /// Enrollment period in months.
    pub enrollment_months: f64,
    /// Data cutoff, months after the first randomisation.
    pub cutoff_months: f64,
}

impl Default for Controls {
    fn default() -> Self {
        Self {
            os_median_months: 15.0,
            os_hr: 0.75,
            pfs_median_months: 4.5,
            pfs_hr: 0.65,
            orr: [0.35, 0.20],
            os_shape: 1.2,
            pfs_shape: 1.1,
            os_pfs_correlation: 0.6,
            withdrawal_prob: 0.08,
            ae_stop_prob: 0.10,
            enrollment_months: 18.0,
            cutoff_months: 30.0,
        }
    }
}

pub(crate) fn docs() -> Vec<ControlDoc> {
    vec![
        ControlDoc::new(
            "os_median_months",
            "Median overall survival of the Control arm, in months",
        ),
        ControlDoc::new(
            "os_hr",
            "Overall-survival hazard ratio, Drug A vs Control (exact proportional hazards)",
        ),
        ControlDoc::new(
            "pfs_median_months",
            "Median progression-free survival of the Control arm, in months (of the latent time to \
             progression or death; observed KM medians fall on the 8-weekly assessment dates)",
        ),
        ControlDoc::new("pfs_hr", "Progression hazard ratio, Drug A vs Control"),
        ControlDoc::new(
            "orr",
            "Confirmed objective response rate (RECIST 1.1): [Drug A, Control]",
        ),
        ControlDoc::new(
            "os_shape",
            "Weibull shape of overall survival (1 = exponential)",
        ),
        ControlDoc::new("pfs_shape", "Weibull shape of the progression time"),
        ControlDoc::new(
            "os_pfs_correlation",
            "Correlation between progression and death (Gaussian copula)",
        ),
        ControlDoc::new("withdrawal_prob", "Probability a subject withdraws early"),
        ControlDoc::new(
            "ae_stop_prob",
            "Probability a subject stops treatment for an adverse event",
        ),
        ControlDoc::new(
            "enrollment_months",
            "Length of the enrollment period in months",
        ),
        ControlDoc::new(
            "cutoff_months",
            "Data cutoff, in months after the first randomisation",
        ),
    ]
}

pub(crate) fn validate(c: &Controls) -> Result<(), String> {
    let positive = |x: f64| x.is_finite() && x > 0.0;
    for (name, x) in [
        ("os_median_months", c.os_median_months),
        ("pfs_median_months", c.pfs_median_months),
        ("os_hr", c.os_hr),
        ("pfs_hr", c.pfs_hr),
    ] {
        if !positive(x) {
            return Err(format!("`{name}` must be positive"));
        }
    }
    if c.pfs_median_months >= c.os_median_months {
        return Err("`pfs_median_months` must be shorter than `os_median_months`".into());
    }
    if !(0.2..=5.0).contains(&c.os_shape) || !(0.2..=5.0).contains(&c.pfs_shape) {
        return Err("Weibull shapes must be between 0.2 and 5".into());
    }
    if !(0.0..=0.95).contains(&c.os_pfs_correlation) {
        return Err("`os_pfs_correlation` must be between 0 and 0.95".into());
    }
    if !c.orr.iter().all(|r| (0.0..=0.9).contains(r)) {
        return Err("`orr` values must be between 0 and 0.9".into());
    }
    if !(0.0..=1.0).contains(&c.withdrawal_prob) || !(0.0..=1.0).contains(&c.ae_stop_prob) {
        return Err("probabilities must be in [0, 1]".into());
    }
    if !positive(c.enrollment_months)
        || c.cutoff_months <= c.enrollment_months
        || c.cutoff_months > 120.0
    {
        return Err(
            "`cutoff_months` must exceed a positive `enrollment_months` (and be at most 120)"
                .into(),
        );
    }
    Ok(())
}

fn month(m: f64) -> i32 {
    (m * DAYS_PER_MONTH).round() as i32
}

/// Tumour assessments every 8 weeks, planned until the data cutoff.
fn schedule(c: &Controls) -> Schedule {
    let interval = 56;
    let k_max = month(c.cutoff_months) / interval + 1;
    let mut visits = vec![PlannedVisit {
        visitnum: 1,
        name: "SCREENING".into(),
        day: Some(-14),
        window: (14, 0),
        epoch: Epoch::Screening,
        kind: VisitKind::Screening,
        labs: false,
        vitals: false,
        endpoint: true,
    }];
    for k in 1..=k_max {
        visits.push(PlannedVisit {
            visitnum: i64::from(k) + 1,
            name: format!("WEEK {}", 8 * k),
            day: Some(interval * k + 1),
            window: (7, 7),
            epoch: Epoch::Treatment,
            kind: VisitKind::Treatment,
            labs: false,
            vitals: false,
            endpoint: true,
        });
    }
    Schedule {
        visits,
        treatment_days: month(c.cutoff_months),
        followup_days: 0,
    }
}

fn base_params() -> OncologyParams {
    OncologyParams {
        os_scale_days: vec![1.0, 1.0],
        os_shape: 1.0,
        pd_scale_days: vec![1.0, 1.0],
        pd_shape: 1.0,
        stratum_hr: vec![0.8, 1.25],
        os_pd_correlation: 0.6,
        response_prob: vec![0.3, 0.2],
        response_pd_correlation: 0.6,
        cr_fraction: 0.15,
        unconfirmed_prob: 0.06,
        withdrawal_prob: 0.08,
        withdrawal_window: (15.0, 400.0),
        withdrawal_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.7),
            ("LOST TO FOLLOW-UP".into(), 0.3),
        ],
        ae_stop_prob: 0.10,
        ae_stop_window: (20.0, 300.0),
        assessment_interval_days: 56,
        max_cycles: 61,
        last_contact_lag_days: 28,
    }
}

/// Seed of the calibration cohort (fixed, so calibration depends on the controls alone).
const CALIBRATION_SEED: u64 = 20_261_013;

/// Calibrated parameters for a control set, computed once per process.
fn calibrated(c: &Controls, design: &Design) -> Calibration {
    static CACHE: OnceLock<Mutex<HashMap<String, Calibration>>> = OnceLock::new();
    let key = serde_json::to_string(c).expect("controls serialize");
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(hit) = cache.lock().expect("calibration cache").get(&key) {
        return hit.clone();
    }
    let mut base = base_params();
    base.os_shape = c.os_shape;
    base.pd_shape = c.pfs_shape;
    base.os_pd_correlation = c.os_pfs_correlation;
    base.withdrawal_prob = c.withdrawal_prob;
    base.ae_stop_prob = c.ae_stop_prob;
    // Calibration depends on the targets and design structure only, not the study's seed or n.
    let mut cal_design = design.clone();
    cal_design.seed = CALIBRATION_SEED;
    let targets = OncologyTargets {
        os_median_days: c.os_median_months * DAYS_PER_MONTH,
        os_hr: vec![c.os_hr, 1.0],
        pfs_median_days: c.pfs_median_months * DAYS_PER_MONTH,
        pd_hr: vec![c.pfs_hr, 1.0],
        orr: c.orr.to_vec(),
        reference_arm: 1,
    };
    let result = calibrate(&cal_design, &base, &targets);
    cache
        .lock()
        .expect("calibration cache")
        .insert(key, result.clone());
    result
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, OncologyParams, Calibration) {
    let first = Date::from_ymd(2023, 3, 1).expect("valid date");
    let treatment = |name: &str, dose: f64, unit: &str| Treatment {
        name: name.into(),
        dose,
        unit: unit.into(),
        form: "SOLUTION".into(),
        route: "INTRAVENOUS".into(),
        frequency: "EVERY 3 WEEKS".into(),
        cycle_days: 21,
    };
    let countries = [
        ("USA", 9),
        ("CAN", 3),
        ("GBR", 5),
        ("DEU", 5),
        ("FRA", 4),
        ("ESP", 4),
    ];
    let sites = countries
        .iter()
        .flat_map(|&(c, k)| std::iter::repeat_n(c, k))
        .enumerate()
        .map(|(i, c)| Site {
            siteid: (101 + i).to_string(),
            country: c.into(),
        })
        .collect();
    let design = Design {
        studyid: opts.studyid.clone().unwrap_or_else(|| "SIMONC02".into()),
        seed: opts.seed,
        n: opts.n,
        arms: vec![
            Arm {
                armcd: "A".into(),
                arm: "Drug A".into(),
                ratio: 1,
                treatment: treatment("DRUG A", 200.0, "mg"),
            },
            Arm {
                armcd: "C".into(),
                arm: "Control".into(),
                ratio: 1,
                treatment: treatment("CONTROL", 75.0, "mg/m2"),
            },
        ],
        sites,
        first_randomization: first,
        enrollment_days: month(c.enrollment_months),
        cutoff: first.add_days(month(c.cutoff_months)),
        block_size: 4,
        population: Population {
            age_mean: 62.0,
            age_sd: 10.0,
            age_min: 28,
            age_max: 85,
            sex: Categorical::new([("M", 0.52), ("F", 0.48)]),
            race: Categorical::new([
                ("WHITE", 0.74),
                ("BLACK OR AFRICAN AMERICAN", 0.10),
                ("ASIAN", 0.12),
                ("AMERICAN INDIAN OR ALASKA NATIVE", 0.02),
                ("MULTIPLE", 0.02),
            ]),
            ethnic: Categorical::new([
                ("NOT HISPANIC OR LATINO", 0.88),
                ("HISPANIC OR LATINO", 0.12),
            ]),
            strata: Categorical::new([("ECOG 0", 0.45), ("ECOG 1", 0.55)]),
        },
        timing: Timing {
            screening_lead: (1, 14),
            consent_lead: (0, 14),
            dose_lag: (0, 3),
        },
        schedule: Some(schedule(c)),
        info: TrialInfo {
            title:
                "A Randomized, Open-Label Study of Drug A versus Control in Advanced Solid Tumors"
                    .into(),
            phase: "PHASE III TRIAL".into(),
            sponsor: "trialforge (synthetic)".into(),
            indication: "Advanced solid tumor".into(),
            therapeutic_area: "ONCOLOGY".into(),
            primary_objective: "To compare progression-free survival between Drug A and Control"
                .into(),
            blinding: "OPEN LABEL".into(),
            control: "ACTIVE".into(),
            add_on: false,
            criteria: vec![
                (
                    "INCL01".into(),
                    "INCLUSION".into(),
                    "Histologically confirmed advanced solid tumor".into(),
                ),
                (
                    "INCL02".into(),
                    "INCLUSION".into(),
                    "At least one measurable lesion per RECIST 1.1".into(),
                ),
                (
                    "INCL03".into(),
                    "INCLUSION".into(),
                    "ECOG performance status 0 or 1".into(),
                ),
                (
                    "EXCL01".into(),
                    "EXCLUSION".into(),
                    "Untreated central nervous system metastases".into(),
                ),
            ],
        },
    };
    let cal = calibrated(c, &design);
    (design, cal.params.clone(), cal)
}
