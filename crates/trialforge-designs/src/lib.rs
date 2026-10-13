//! Study-type archetypes: default specs, type-specific controls and calibration recipes.
//!
//! The spike has one archetype, `oncology_solid_tumor`, with SIMONC01's design: 1:1 Drug A vs
//! Control, permuted blocks of 4 stratified by ECOG (0/1), 18 months of enrollment from
//! 2023-03-01, data cutoff 30 months after the first randomisation, Q3W dosing and tumour
//! assessments every 8 weeks. Calibration to clinical targets comes in Phase 3, so the latent
//! medians here are set directly.

use trialforge_core::Date;
use trialforge_model::{Arm, Categorical, Design, Population, Site, Treatment};
use trialforge_processes::oncology::OncologyParams;

/// The archetypes this build knows, with one-line descriptions.
pub const ARCHETYPES: &[(&str, &str)] = &[(
    "oncology_solid_tumor",
    "Randomised solid-tumour study (SIMONC01-like): DM, EX, DS and ADSL in this preview",
)];

/// The largest seed accepted: 2^53, the largest integer every binding represents exactly.
pub const MAX_SEED: u64 = 1 << 53;

/// User-facing options shared by all archetypes.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StudyOptions {
    /// Number of randomised subjects.
    pub n: usize,
    /// Root seed. At most [`MAX_SEED`], so it survives R doubles and JSON numbers.
    pub seed: u64,
    /// `STUDYID`; each archetype has a default.
    pub studyid: Option<String>,
}

/// A resolved oncology study: its design and course parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct OncologyStudy {
    /// The design.
    pub design: Design,
    /// The course parameters.
    pub params: OncologyParams,
}

/// An archetype name this build does not know.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnknownArchetype(pub String);

impl core::fmt::Display for UnknownArchetype {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        let known: Vec<&str> = ARCHETYPES.iter().map(|(n, _)| *n).collect();
        write!(
            f,
            "unknown archetype {:?}; available: {}",
            self.0,
            known.join(", ")
        )
    }
}

impl std::error::Error for UnknownArchetype {}

/// Resolves an archetype name and options to a concrete study.
pub fn resolve(archetype: &str, opts: &StudyOptions) -> Result<OncologyStudy, UnknownArchetype> {
    match archetype {
        "oncology_solid_tumor" => Ok(oncology_solid_tumor(opts)),
        other => Err(UnknownArchetype(other.to_string())),
    }
}

fn oncology_solid_tumor(opts: &StudyOptions) -> OncologyStudy {
    let first = Date::from_ymd(2023, 3, 1).expect("valid date");
    let treatment = |name: &str, dose: f64, unit: &str| Treatment {
        name: name.into(),
        dose,
        unit: unit.into(),
        form: "SOLUTION".into(),
        route: "INTRAVENOUS".into(),
        frequency: "Q3W".into(),
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
        enrollment_days: 548,
        cutoff: first.add_days(912),
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
    };
    let params = OncologyParams {
        median_pd_months: vec![7.0, 4.5],
        median_pps_months: vec![13.0, 13.0],
        stratum_hr: vec![0.8, 1.25],
        death_before_pd_rate: 0.008,
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
    };
    OncologyStudy { design, params }
}
