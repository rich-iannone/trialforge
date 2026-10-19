//! Study-spec files (`spec/study-spec.schema.json`): JSON or YAML.
//!
//! The engine currently consumes `study.archetype`, `study.seed`, `study.studyid`,
//! `population.n`, `controls`, `data_quality.defects` and `outputs`. Sections the schema
//! sketches for later phases (`design`, `endpoints`, `safety`, `dropout`) are rejected with a
//! pointer to `controls`, rather than silently ignored.

use serde::Deserialize;
use serde_json::Value;

use crate::{Error, StudyOptions, defects::DefectSpec};

/// A parsed study spec.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Spec {
    /// Spec format version (0).
    #[serde(default)]
    pub spec_version: Option<u32>,
    /// The study section.
    pub study: StudySection,
    /// The population section.
    #[serde(default)]
    pub population: PopulationSection,
    /// Archetype controls (see `controls(archetype)`).
    #[serde(default)]
    pub controls: Value,
    /// Data-quality defects.
    #[serde(default)]
    pub data_quality: DataQualitySection,
    /// What to write.
    #[serde(default)]
    pub outputs: OutputSection,
    #[serde(default)]
    design: Option<Value>,
    #[serde(default)]
    endpoints: Option<Value>,
    #[serde(default)]
    safety: Option<Value>,
    #[serde(default)]
    dropout: Option<Value>,
}

/// `study:`.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct StudySection {
    /// Archetype name.
    pub archetype: String,
    /// Root seed (required: specs are for reproducible studies).
    pub seed: u64,
    /// `STUDYID`.
    #[serde(default)]
    pub studyid: Option<String>,
    /// Standards versions; only SDTMIG 3.4 / ADaMIG 1.3 / Define 2.1 / CT 2024-09-27 for now.
    #[serde(default)]
    pub standards: Option<Standards>,
    /// Not supported yet (each archetype fixes its calendar).
    #[serde(default)]
    start_date: Option<Value>,
    /// Accepted for forward compatibility; no modifiers are implemented yet.
    #[serde(default)]
    modifiers: Vec<String>,
}

/// `study.standards:`.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Standards {
    /// SDTMIG version.
    #[serde(default)]
    pub sdtmig: Option<String>,
    /// ADaMIG version.
    #[serde(default)]
    pub adamig: Option<String>,
    /// Define-XML version.
    #[serde(default)]
    pub define: Option<String>,
    /// CT package date.
    #[serde(default)]
    pub ct: Option<String>,
}

/// `population:`.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PopulationSection {
    /// Number of randomised subjects.
    #[serde(default = "default_n")]
    pub n: usize,
}

impl Default for PopulationSection {
    fn default() -> Self {
        Self { n: default_n() }
    }
}

fn default_n() -> usize {
    300
}

/// `data_quality:`.
#[derive(Clone, Debug, Default, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct DataQualitySection {
    /// Defects to inject into the clean study: `none`, or a map of defect id to a count (an
    /// integer) or a rate (see `trialforge_defects::CATALOGUE`).
    #[serde(default)]
    pub defects: Option<DefectSpec>,
    /// Not supported yet.
    #[serde(default)]
    partial_dates: Option<Value>,
    /// Not supported yet.
    #[serde(default)]
    missed_visit: Option<Value>,
}

/// `outputs:`.
#[derive(Clone, Debug, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputSection {
    /// Formats for the flat layout (see `FORMATS`).
    #[serde(default = "default_formats")]
    pub formats: Vec<String>,
    /// `flat` (one folder, every format) or `submission` (`sdtm/` and `adam/` with XPT and
    /// define.xml).
    #[serde(default = "default_layout")]
    pub layout: String,
}

impl Default for OutputSection {
    fn default() -> Self {
        Self {
            formats: default_formats(),
            layout: default_layout(),
        }
    }
}

fn default_formats() -> Vec<String> {
    vec!["xpt".into()]
}

fn default_layout() -> String {
    "flat".into()
}

impl Spec {
    /// Parses a JSON or YAML spec (JSON if it starts with `{`).
    pub fn parse(text: &str) -> Result<Self, Error> {
        let value: Value = if text.trim_start().starts_with('{') {
            serde_json::from_str(text)
                .map_err(|e| Error::Invalid(format!("spec is not valid JSON: {e}")))?
        } else {
            serde_norway::from_str(text)
                .map_err(|e| Error::Invalid(format!("spec is not valid YAML: {e}")))?
        };
        Self::from_value(value)
    }

    /// Builds a spec from an already-parsed JSON value.
    pub fn from_value(value: Value) -> Result<Self, Error> {
        let spec: Spec = serde_json::from_value(value)
            .map_err(|e| Error::Invalid(format!("invalid spec: {e}")))?;
        spec.validate()?;
        Ok(spec)
    }

    fn validate(&self) -> Result<(), Error> {
        let later = |name: &str| {
            Error::Invalid(format!(
                "spec section `{name}` is not supported yet; set the equivalent archetype controls under `controls` \
                 (see `controls(\"{}\")`)",
                self.study.archetype
            ))
        };
        for (name, v) in [
            ("design", &self.design),
            ("endpoints", &self.endpoints),
            ("safety", &self.safety),
            ("dropout", &self.dropout),
            ("study.start_date", &self.study.start_date),
        ] {
            if v.is_some() {
                return Err(later(name));
            }
        }
        for (name, v) in [
            ("partial_dates", &self.data_quality.partial_dates),
            ("missed_visit", &self.data_quality.missed_visit),
        ] {
            if v.is_some() {
                return Err(Error::Invalid(format!(
                    "`data_quality.{name}` is not supported yet; `data_quality.defects` is"
                )));
            }
        }
        if !self.study.modifiers.is_empty() {
            return Err(Error::Invalid(
                "study modifiers are not implemented yet".into(),
            ));
        }
        if let Some(v) = self.spec_version.filter(|&v| v != 0) {
            return Err(Error::Invalid(format!(
                "unsupported spec_version {v}; this build reads version 0"
            )));
        }
        if let Some(s) = &self.study.standards {
            for (name, got, want) in [
                ("sdtmig", &s.sdtmig, "3.4"),
                ("adamig", &s.adamig, "1.3"),
                ("define", &s.define, "2.1"),
                ("ct", &s.ct, "2024-09-27"),
            ] {
                if let Some(g) = got.as_deref().filter(|g| *g != want) {
                    return Err(Error::Invalid(format!(
                        "standards.{name} {g:?} is not available; this build supports {want}"
                    )));
                }
            }
        }
        if !matches!(self.outputs.layout.as_str(), "flat" | "submission") {
            return Err(Error::Invalid(format!(
                "outputs.layout must be \"flat\" or \"submission\", not {:?}",
                self.outputs.layout
            )));
        }
        if let Some(bad) = self
            .outputs
            .formats
            .iter()
            .find(|f| !crate::FORMATS.contains(&f.as_str()))
        {
            return Err(Error::Invalid(format!(
                "unknown output format {bad:?}; available: {}",
                crate::FORMATS.join(", ")
            )));
        }
        Ok(())
    }

    /// The options to generate this study with.
    pub fn options(&self) -> StudyOptions {
        StudyOptions {
            n: self.population.n,
            seed: self.study.seed,
            studyid: self.study.studyid.clone(),
            controls: self.controls.clone(),
            defects: self.data_quality.defects.clone(),
        }
    }
}

/// Parses a spec and generates its study.
pub fn study_from_spec(text: &str) -> Result<(crate::Study, Spec), Error> {
    let spec = Spec::parse(text)?;
    let study = crate::study(&spec.study.archetype, &spec.options())?;
    Ok((study, spec))
}

#[cfg(test)]
mod tests {
    use super::*;

    const YAML: &str = "
study:
  archetype: parallel_continuous
  studyid: DEMO01
  seed: 42
  standards: { sdtmig: \"3.4\", adamig: \"1.3\" }
population:
  n: 20
controls:
  indication: hypertension
  weeks: 12
outputs:
  formats: [xpt, csv]
  layout: submission
";

    #[test]
    fn yaml_and_json_agree() {
        let y = Spec::parse(YAML).unwrap();
        let as_json =
            serde_json::to_string(&serde_norway::from_str::<Value>(YAML).unwrap()).unwrap();
        assert_eq!(y, Spec::parse(&as_json).unwrap());
        assert_eq!(y.options().n, 20);
        assert_eq!(y.controls["weeks"], 12);
        assert_eq!(y.outputs.layout, "submission");
    }

    #[test]
    fn generates_and_rejects_unsupported_sections() {
        let (study, _) = study_from_spec(YAML).unwrap();
        assert_eq!(study.studyid, "DEMO01");
        let err = |text: &str| Spec::parse(text).unwrap_err().to_string();
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1}\ndesign: {arms: []}")
                .contains("controls")
        );
        assert!(err("study: {archetype: parallel_continuous}").contains("seed"));
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1}\nfoo: 1")
                .contains("unknown field")
        );
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1, standards: {sdtmig: \"3.3\"}}")
                .contains("3.4")
        );
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1}\noutputs: {layout: zip}")
                .contains("layout")
        );
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1}\ndata_quality: {defects: {typo: 1}}")
                .contains("unknown defect")
        );
        assert!(
            err("study: {archetype: parallel_continuous, seed: 1}\ndata_quality: {missed_visit: 0.1}")
                .contains("not supported yet")
        );
    }

    #[test]
    fn reads_defects() {
        let base = "study: {archetype: parallel_continuous, seed: 1}\ndata_quality:\n  defects: ";
        let spec = Spec::parse(&format!("{base}{{out_of_ct: 2, duplicate_seq: 0.5}}")).unwrap();
        let d = spec.options().defects.unwrap();
        assert_eq!(d.defects.len(), 2);
        let none = Spec::parse(&format!("{base}none")).unwrap();
        assert!(none.options().defects.unwrap().is_empty());
    }
}
