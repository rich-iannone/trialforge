//! Study-type archetypes: default designs, typed **controls** with defaults, and the course
//! parameters they resolve to.
//!
//! Controls arrive as a JSON object (from Python keyword arguments, an R list, or a study-spec
//! file) and are merged over the archetype's defaults; unknown names are an error.

pub mod oncology;
pub mod parallel;

use serde::Serialize;
use serde_json::Value;
use trialforge_model::Design;
use trialforge_processes::fixed::CourseParams;
use trialforge_processes::oncology::OncologyParams;

/// The archetypes this build knows, with one-line descriptions.
pub const ARCHETYPES: &[(&str, &str)] = &[
    (
        "parallel_continuous",
        "Randomised, double-blind, placebo-controlled parallel study with a continuous endpoint (HbA1c, SBP or weight)",
    ),
    (
        "oncology_solid_tumor",
        "Randomised solid-tumour study (SIMONC01-like): progression, death and disposition; no tumour data yet",
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
pub enum Course {
    /// A fixed-duration study with a visit schedule.
    Fixed(CourseParams),
    /// An event-driven oncology course.
    Oncology(OncologyParams),
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
            let (design, params) = oncology::resolve(opts, &c);
            Ok(Resolved {
                design,
                course: Course::Oncology(params),
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
