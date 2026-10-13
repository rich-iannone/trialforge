//! `oncology_solid_tumor`: SIMONC01's design. 1:1 Drug A vs Control, permuted blocks of 4
//! stratified by ECOG (0/1), 18 months of enrollment from 2023-03-01, data cutoff 30 months
//! after the first randomisation, Q3W dosing and tumour assessments every 8 weeks.
//! Calibration to clinical targets comes in Phase 3, so latent medians are set directly.

use serde::{Deserialize, Serialize};
use trialforge_core::Date;
use trialforge_model::{Arm, Categorical, Design, Population, Site, Timing, Treatment, TrialInfo};
use trialforge_processes::oncology::OncologyParams;

use crate::{ControlDoc, StudyOptions};

/// Controls for `oncology_solid_tumor`.
#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Controls {
    /// Latent median time to progression in months: [Drug A, Control].
    pub median_pd_months: [f64; 2],
    /// Latent median post-progression survival in months: [Drug A, Control].
    pub median_pps_months: [f64; 2],
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
            median_pd_months: [7.0, 4.5],
            median_pps_months: [13.0, 13.0],
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
            "median_pd_months",
            "Latent median time to progression in months: [Drug A, Control]",
        ),
        ControlDoc::new(
            "median_pps_months",
            "Latent median post-progression survival in months: [Drug A, Control]",
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
    if !c
        .median_pd_months
        .iter()
        .chain(&c.median_pps_months)
        .all(|&m| positive(m))
    {
        return Err("medians must be positive".into());
    }
    if !(0.0..=1.0).contains(&c.withdrawal_prob) || !(0.0..=1.0).contains(&c.ae_stop_prob) {
        return Err("probabilities must be in [0, 1]".into());
    }
    if !positive(c.enrollment_months) || c.cutoff_months <= c.enrollment_months {
        return Err("`cutoff_months` must exceed a positive `enrollment_months`".into());
    }
    Ok(())
}

pub(crate) fn resolve(opts: &StudyOptions, c: &Controls) -> (Design, OncologyParams) {
    let first = Date::from_ymd(2023, 3, 1).expect("valid date");
    let month = |m: f64| (m * trialforge_processes::oncology::DAYS_PER_MONTH).round() as i32;
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
        schedule: None,
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
    let params = OncologyParams {
        median_pd_months: c.median_pd_months.to_vec(),
        median_pps_months: c.median_pps_months.to_vec(),
        stratum_hr: vec![0.8, 1.25],
        death_before_pd_rate: 0.008,
        withdrawal_prob: c.withdrawal_prob,
        withdrawal_window: (15.0, 400.0),
        withdrawal_reasons: vec![
            ("WITHDRAWAL BY SUBJECT".into(), 0.7),
            ("LOST TO FOLLOW-UP".into(), 0.3),
        ],
        ae_stop_prob: c.ae_stop_prob,
        ae_stop_window: (20.0, 300.0),
        assessment_interval_days: 56,
        max_cycles: 61,
        last_contact_lag_days: 28,
    };
    (design, params)
}
