//! Study-type archetypes: default designs, typed **controls** with defaults, and the course
//! parameters they resolve to.
//!
//! Controls arrive as a JSON object (from Python keyword arguments, an R list, or a study-spec
//! file) and are merged over the archetype's defaults; unknown names are an error.

pub mod cns;
pub mod crossover;
pub mod cv;
pub mod dose_ranging;
pub mod escalation;
pub mod lymphoma;
pub mod oncology;
pub mod ophthalmology;
pub mod parallel;
pub mod pediatric;
pub mod phase1;
pub mod respiratory;
pub mod vaccine;

use serde::Serialize;
use serde_json::Value;
use trialforge_model::Design;
use trialforge_processes::adas::AdasParams;
use trialforge_processes::escalation::{EscalationParams, Trace};
use trialforge_processes::exacerbation::ExacerbationParams;
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::growth::GrowthParams;
use trialforge_processes::oncology::OncologyParams;
use trialforge_processes::ophtha::OphthaParams;
use trialforge_processes::outcomes::OutcomesParams;
use trialforge_processes::pk::PkParams;
use trialforge_processes::vaccine::VaccineParams;

/// The archetypes this build knows, with one-line descriptions.
pub const ARCHETYPES: &[(&str, &str)] = &[
    (
        "parallel_continuous",
        "Randomised, double-blind, placebo-controlled parallel study with a continuous endpoint (HbA1c, SBP or weight)",
    ),
    (
        "heme_lymphoma",
        "Randomised study in relapsed or refractory DLBCL: Lugano 2014 CT and PET-CT responses (TU, TR, RS), ADRS, PFS and OS",
    ),
    (
        "oncology_solid_tumor",
        "Randomised solid-tumour study (SIMONC01-like): RECIST 1.1 tumour assessments, response, progression and survival",
    ),
    (
        "cns_alzheimers",
        "Placebo-controlled study in mild to moderate Alzheimer's disease with ADAS-Cog (11) items (FT) and LOCF analysis",
    ),
    (
        "crossover_be",
        "Crossover bioequivalence or food-effect study (2x2, 2x2x4 replicate or fed/fasted) with PK and NCA",
    ),
    (
        "cv_outcomes",
        "Event-driven cardiovascular outcomes trial: recurrent MI and stroke, CV and non-CV death, time to first MACE",
    ),
    (
        "dose_escalation",
        "Phase 1 dose escalation in advanced solid tumours (3+3, BOIN or mTPI-2): DLTs in SUPPAE, DLT-evaluable population, MTD",
    ),
    (
        "dose_ranging_ph2",
        "Phase 2 dose-ranging study (placebo + several doses) with a sigmoid Emax dose-response",
    ),
    (
        "ophthalmology",
        "Active-controlled non-inferiority study of intravitreal Drug A in neovascular AMD: study and fellow eye BCVA and retinal thickness (OE)",
    ),
    (
        "pediatric_growth",
        "Active-controlled study of once-weekly Drug A vs daily somatropin in children with growth hormone deficiency: height velocity, CDC growth z-scores",
    ),
    (
        "phase1_sad_mad",
        "First-in-human single or multiple ascending dose study in healthy volunteers, with intensive PK and NCA",
    ),
    (
        "vaccine",
        "Placebo-controlled vaccine study: neutralising antibody titres (IS), 7-day reactogenicity diaries (CE, FACE), ADIS/ADCE/ADFACE",
    ),
    (
        "respiratory_exacerbation",
        "Placebo-controlled COPD study: moderate and severe exacerbations (CE, negative binomial rate ratio) and trough FEV1 (RE)",
    ),
];

/// The largest seed accepted: 2^53, the largest integer every binding represents exactly.
pub const MAX_SEED: u64 = 1 << 53;

/// User-facing options shared by all archetypes.
#[derive(Clone, Debug, PartialEq)]
pub struct StudyOptions {
    /// Number of randomised subjects.
    pub n: usize,
    /// Root seed. At most [`MAX_SEED`], so it survives R doubles and JSON numbers.
    pub seed: u64,
    /// `STUDYID`; each archetype has a default.
    pub studyid: Option<String>,
    /// Archetype controls as a JSON object (`Value::Null` or `{}` for the defaults).
    pub controls: Value,
}

impl StudyOptions {
    /// Options with default controls.
    pub fn new(n: usize, seed: u64) -> Self {
        Self {
            n,
            seed,
            studyid: None,
            controls: Value::Null,
        }
    }
}

/// One documented control.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct ControlDoc {
    /// Control name.
    pub name: &'static str,
    /// What it does.
    pub description: &'static str,
    /// Default value.
    pub default: Value,
}

impl ControlDoc {
    pub(crate) fn new(name: &'static str, description: &'static str) -> Self {
        Self {
            name,
            description,
            default: Value::Null,
        }
    }
}

/// What an archetype resolves to.
#[derive(Clone, Debug, PartialEq)]
// One value per study: boxing the larger variant would buy nothing.
#[allow(clippy::large_enum_variant)]
pub enum Course {
    /// A fixed-duration study with a visit schedule.
    Fixed(CourseParams),
    /// An event-driven oncology course, with the calibration that produced its parameters.
    Oncology(OncologyParams, trialforge_calibrate::oncology::Calibration),
    /// A fixed-duration course with timed dosing and pharmacokinetics layered on it.
    Pk(CourseParams, PkParams),
    /// A crossover: the fixed-duration course split into periods (each period's dose study
    /// day), with pharmacokinetics.
    Crossover(CourseParams, PkParams, Vec<i32>),
    /// A fixed-duration course with ADAS-Cog assessments layered on it.
    Adas(CourseParams, AdasParams),
    /// An event-driven outcomes course.
    Outcomes(OutcomesParams),
    /// A fixed-duration course with ophthalmic examinations and intravitreal injections.
    Ophtha(CourseParams, OphthaParams),
    /// A fixed-duration paediatric course with growth measurements.
    Growth(CourseParams, GrowthParams),
    /// A fixed-duration respiratory course with exacerbations and lung function.
    Exacerbation(CourseParams, ExacerbationParams),
    /// A fixed-duration vaccine course: vaccinations, reactogenicity and immunogenicity.
    Vaccine(CourseParams, VaccineParams),
    /// A dose-escalation course: the cohorts the design treated and their DLTs.
    Escalation(CourseParams, EscalationParams, Trace),
}

/// A resolved study: design, course parameters and the controls actually used.
#[derive(Clone, Debug, PartialEq)]
pub struct Resolved {
    /// The design.
    pub design: Design,
    /// The course parameters.
    pub course: Course,
    /// The full control set used (defaults merged with the user's), as JSON.
    pub controls: Value,
}

/// Bad archetype name or controls.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ResolveError(pub String);

impl core::fmt::Display for ResolveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for ResolveError {}

fn unknown(archetype: &str) -> ResolveError {
    let known: Vec<&str> = ARCHETYPES.iter().map(|(n, _)| *n).collect();
    ResolveError(format!(
        "unknown archetype {archetype:?}; available: {}",
        known.join(", ")
    ))
}

/// Merges `user` over the default controls `T` and parses the result.
fn controls<T>(user: &Value) -> Result<(T, Value), ResolveError>
where
    T: Default + Serialize + serde::de::DeserializeOwned,
{
    let mut merged = serde_json::to_value(T::default()).expect("controls serialize");
    match user {
        Value::Null => {}
        Value::Object(map) => {
            let target = merged.as_object_mut().expect("controls are an object");
            for (k, v) in map {
                if !target.contains_key(k) {
                    let names: Vec<&str> = target.keys().map(String::as_str).collect();
                    return Err(ResolveError(format!(
                        "unknown control {k:?}; available: {}",
                        names.join(", ")
                    )));
                }
                target.insert(k.clone(), v.clone());
            }
        }
        _ => return Err(ResolveError("controls must be a JSON object".into())),
    }
    let parsed: T = serde_json::from_value(merged.clone())
        .map_err(|e| ResolveError(format!("invalid controls: {e}")))?;
    Ok((parsed, merged))
}

/// The documented controls of `archetype`, with their defaults.
pub fn control_docs(archetype: &str) -> Result<Vec<ControlDoc>, ResolveError> {
    let (mut docs, defaults) = match archetype {
        "parallel_continuous" => (
            parallel::docs(),
            serde_json::to_value(parallel::Controls::default()),
        ),
        "oncology_solid_tumor" => (
            oncology::docs(),
            serde_json::to_value(oncology::Controls::default()),
        ),
        "cns_alzheimers" => (cns::docs(), serde_json::to_value(cns::Controls::default())),
        "cv_outcomes" => (cv::docs(), serde_json::to_value(cv::Controls::default())),
        "dose_ranging_ph2" => (
            dose_ranging::docs(),
            serde_json::to_value(dose_ranging::Controls::default()),
        ),
        "crossover_be" => (
            crossover::docs(),
            serde_json::to_value(crossover::Controls::default()),
        ),
        "ophthalmology" => (
            ophthalmology::docs(),
            serde_json::to_value(ophthalmology::Controls::default()),
        ),
        "pediatric_growth" => (
            pediatric::docs(),
            serde_json::to_value(pediatric::Controls::default()),
        ),
        "phase1_sad_mad" => (
            phase1::docs(),
            serde_json::to_value(phase1::Controls::default()),
        ),
        "heme_lymphoma" => (
            lymphoma::docs(),
            serde_json::to_value(lymphoma::Controls::default()),
        ),
        "dose_escalation" => (
            escalation::docs(),
            serde_json::to_value(escalation::Controls::default()),
        ),
        "vaccine" => (
            vaccine::docs(),
            serde_json::to_value(vaccine::Controls::default()),
        ),
        "respiratory_exacerbation" => (
            respiratory::docs(),
            serde_json::to_value(respiratory::Controls::default()),
        ),
        other => return Err(unknown(other)),
    };
    let defaults = defaults.expect("controls serialize");
    for d in &mut docs {
        d.default = defaults.get(d.name).cloned().unwrap_or(Value::Null);
    }
    Ok(docs)
}

/// Resolves an archetype name and options to a concrete study.
pub fn resolve(archetype: &str, opts: &StudyOptions) -> Result<Resolved, ResolveError> {
    match archetype {
        "parallel_continuous" => {
            let (c, json) = controls::<parallel::Controls>(&opts.controls)?;
            parallel::validate(&c).map_err(ResolveError)?;
            let (design, params) = parallel::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Fixed(params),
                controls: json,
            })
        }
        "oncology_solid_tumor" => {
            let (c, json) = controls::<oncology::Controls>(&opts.controls)?;
            oncology::validate(&c).map_err(ResolveError)?;
            let (design, params, calibration) = oncology::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Oncology(params, calibration),
                controls: json,
            })
        }
        "cns_alzheimers" => {
            let (c, json) = controls::<cns::Controls>(&opts.controls)?;
            cns::validate(&c).map_err(ResolveError)?;
            let (design, course, adas) = cns::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Adas(course, adas),
                controls: json,
            })
        }
        "cv_outcomes" => {
            let (c, json) = controls::<cv::Controls>(&opts.controls)?;
            cv::validate(&c).map_err(ResolveError)?;
            let (design, params) = cv::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Outcomes(params),
                controls: json,
            })
        }
        "dose_ranging_ph2" => {
            let (c, json) = controls::<dose_ranging::Controls>(&opts.controls)?;
            dose_ranging::validate(&c).map_err(ResolveError)?;
            let (design, params) = dose_ranging::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Fixed(params),
                controls: json,
            })
        }
        "crossover_be" => {
            let (c, json) = controls::<crossover::Controls>(&opts.controls)?;
            crossover::validate(&c).map_err(ResolveError)?;
            let (design, course, pk) = crossover::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Crossover(course, pk, crossover::period_days(&c)),
                controls: json,
            })
        }
        "ophthalmology" => {
            let (c, json) = controls::<ophthalmology::Controls>(&opts.controls)?;
            ophthalmology::validate(&c).map_err(ResolveError)?;
            let (design, course, ophtha) = ophthalmology::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Ophtha(course, ophtha),
                controls: json,
            })
        }
        "pediatric_growth" => {
            let (c, json) = controls::<pediatric::Controls>(&opts.controls)?;
            pediatric::validate(&c).map_err(ResolveError)?;
            let (design, course, growth) = pediatric::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Growth(course, growth),
                controls: json,
            })
        }
        "phase1_sad_mad" => {
            let (c, json) = controls::<phase1::Controls>(&opts.controls)?;
            phase1::validate(&c).map_err(ResolveError)?;
            let (design, course, pk) = phase1::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Pk(course, pk),
                controls: json,
            })
        }
        "heme_lymphoma" => {
            let (c, json) = controls::<lymphoma::Controls>(&opts.controls)?;
            oncology::validate(&c.0).map_err(ResolveError)?;
            let (design, params, calibration) = lymphoma::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Oncology(params, calibration),
                controls: json,
            })
        }
        "dose_escalation" => {
            let (c, json) = controls::<escalation::Controls>(&opts.controls)?;
            escalation::validate(&c, opts.n).map_err(ResolveError)?;
            let (design, course, p, trace) = escalation::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Escalation(course, p, trace),
                controls: json,
            })
        }
        "vaccine" => {
            let (c, json) = controls::<vaccine::Controls>(&opts.controls)?;
            vaccine::validate(&c).map_err(ResolveError)?;
            let (design, course, p) = vaccine::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Vaccine(course, p),
                controls: json,
            })
        }
        "respiratory_exacerbation" => {
            let (c, json) = controls::<respiratory::Controls>(&opts.controls)?;
            respiratory::validate(&c).map_err(ResolveError)?;
            let (design, course, p) = respiratory::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Exacerbation(course, p),
                controls: json,
            })
        }
        other => Err(unknown(other)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_archetype_resolves_with_defaults() {
        for (name, _) in ARCHETYPES {
            let r = resolve(name, &StudyOptions::new(10, 1)).unwrap();
            r.design.validate().unwrap();
            assert!(!control_docs(name).unwrap().is_empty());
        }
    }

    #[test]
    fn controls_merge_and_reject_unknowns() {
        let mut o = StudyOptions::new(10, 1);
        o.controls = serde_json::json!({"weeks": 12, "indication": "hypertension"});
        let r = resolve("parallel_continuous", &o).unwrap();
        assert_eq!(r.controls["weeks"], 12);
        assert_eq!(r.design.schedule.as_ref().unwrap().treatment_days, 84);
        o.controls = serde_json::json!({"wekks": 12});
        assert!(
            resolve("parallel_continuous", &o)
                .unwrap_err()
                .0
                .contains("unknown control")
        );
        o.controls = serde_json::json!({"weeks": "twelve"});
        assert!(
            resolve("parallel_continuous", &o)
                .unwrap_err()
                .0
                .contains("invalid controls")
        );
        o.controls = serde_json::json!({"indication": "gout"});
        assert!(
            resolve("parallel_continuous", &o)
                .unwrap_err()
                .0
                .contains("unknown indication")
        );
    }
}
