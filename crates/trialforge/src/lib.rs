//! trialforge: synthetic CDISC clinical-trial data from declarative study specs.
//!
//! This facade re-exports the workspace crates and runs the pipeline:
//! archetype + controls → design → enrollment → subject timelines → SDTM → ADaM → checks →
//! truth manifest, then, when asked for, seeded data-quality defects. The Python and R
//! bindings and the CLI depend on this crate alone, so the internal crate layout can change
//! without touching them.

pub mod spec;
mod truth;

use std::collections::HashMap;
use std::path::Path;

pub use trialforge_adam as adam;
pub use trialforge_calibrate as calibrate;
pub use trialforge_cdisc as cdisc;
pub use trialforge_check as check;
pub use trialforge_core as core;
pub use trialforge_defects as defects;
pub use trialforge_designs as designs;
pub use trialforge_hooks as hooks;
pub use trialforge_io as io;
pub use trialforge_model as model;
pub use trialforge_processes as processes;
pub use trialforge_sdtm as sdtm;

pub use spec::{Spec, study_from_spec};
pub use trialforge_core::stats::km_median;
pub use trialforge_core::{Column, Dataset, Date, Values};
pub use trialforge_designs::{ARCHETYPES, MAX_SEED, StudyOptions};

/// The version of the trialforge core.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Output formats [`Study::write`] understands in this build.
#[cfg(feature = "parquet")]
pub const FORMATS: &[&str] = &["xpt", "csv", "dataset-json", "parquet"];
/// Output formats [`Study::write`] understands in this build.
#[cfg(not(feature = "parquet"))]
pub const FORMATS: &[&str] = &["xpt", "csv", "dataset-json"];

/// A generated study.
#[derive(Clone, Debug, PartialEq)]
pub struct Study {
    /// `STUDYID`.
    pub studyid: String,
    /// The archetype it was generated from.
    pub archetype: String,
    /// The root seed (rerun with it to reproduce the study).
    pub seed: u64,
    /// The data cutoff (event-driven studies) or the last contact date (fixed-duration).
    pub cutoff: Date,
    /// The trial title (TS `TITLE`).
    pub title: String,
    /// The full control set used, as JSON.
    pub controls: serde_json::Value,
    /// SDTM datasets, in standard order.
    pub sdtm: Vec<Dataset>,
    /// ADaM datasets, in standard order.
    pub adam: Vec<Dataset>,
    /// Targets vs realised values, as JSON (of the clean study, before any defects).
    pub truth: serde_json::Value,
    /// The defects injected into `sdtm` and `adam`, if any were requested.
    pub defects: Option<defects::Manifest>,
}

/// Anything that can go wrong generating or writing a study.
#[derive(Debug)]
pub enum Error {
    /// Bad user input (unknown archetype, invalid options, unknown format, …).
    Invalid(String),
    /// Generated data broke its own metadata or derivation contract. This is a bug.
    Internal(String),
    /// Writing output failed.
    Io(std::io::Error),
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Error::Invalid(m) => f.write_str(m),
            Error::Internal(m) => write!(f, "internal error (please report): {m}"),
            Error::Io(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for Error {}

/// A fresh seed in `[0, MAX_SEED]`, for callers that did not supply one. The seed is recorded
/// on the [`Study`], so the result is still reproducible.
pub fn random_seed() -> u64 {
    use std::hash::{BuildHasher, Hasher};
    // RandomState is seeded from OS randomness; this is not an output path, only seed choice.
    let mut h = std::collections::hash_map::RandomState::new().build_hasher();
    h.write_u128(
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map_or(0, |d| d.as_nanos()),
    );
    h.finish() % (MAX_SEED + 1)
}

/// The documented controls of `archetype` (name, description, default) as JSON.
pub fn controls(archetype: &str) -> Result<serde_json::Value, Error> {
    let docs = designs::control_docs(archetype).map_err(|e| Error::Invalid(e.to_string()))?;
    Ok(serde_json::to_value(docs).expect("control docs serialize"))
}

/// Generates a study from an archetype.
pub fn study(archetype: &str, opts: &StudyOptions) -> Result<Study, Error> {
    if opts.seed > MAX_SEED {
        return Err(Error::Invalid(format!(
            "`seed` must be at most 2^53 ({MAX_SEED})"
        )));
    }
    let mut resolved =
        designs::resolve(archetype, opts).map_err(|e| Error::Invalid(e.to_string()))?;
    let design = &resolved.design;
    let mut subjects = model::enroll(design).map_err(|e| Error::Invalid(e.to_string()))?;
    let internal = |e: &dyn std::fmt::Display| Error::Internal(e.to_string());

    let (timelines, strata, cutoff, interim) = match &resolved.course {
        designs::Course::Fixed(p) => {
            let t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            (t, HashMap::new(), None, false)
        }
        designs::Course::Pk(p, pk) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::pk::simulate(design, pk, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Crossover(p, pk, days) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            for x in &mut t {
                processes::crossover::apply(x, days, designs::crossover::IN_HOUSE_DAYS);
            }
            processes::pk::simulate(design, pk, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Adas(p, adas) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::adas::simulate(design, adas, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Ophtha(p, ophtha) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::ophtha::simulate(design, ophtha, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Growth(p, growth) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::growth::simulate(design, growth, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Exacerbation(p, ex) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::exacerbation::simulate(design, ex, &subjects, &mut t);
            (t, strata_of(design, &subjects), None, false)
        }
        designs::Course::Vaccine(p, vx) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::vaccine::simulate(design, vx, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Escalation(p, esc, _) => {
            let mut t = processes::fixed::simulate(design, p, &subjects)
                .map_err(|e| Error::Invalid(e.to_string()))?;
            processes::escalation::simulate(design, esc, &subjects, &mut t);
            (t, HashMap::new(), None, false)
        }
        designs::Course::Outcomes(p) => {
            let (t, cut) = processes::outcomes::simulate(design, p, &subjects);
            let interim = p.interim_events.is_some();
            (t, strata_of(design, &subjects), Some(cut), interim)
        }
        designs::Course::Oncology(p, _) => {
            let t = processes::oncology::simulate(design, p, &subjects);
            let cut = p.data_cutoff.unwrap_or(design.cutoff);
            (
                t,
                strata_of(design, &subjects),
                Some(cut),
                p.data_cutoff.is_some(),
            )
        }
    };
    // An interim snapshot reports its own cutoff (TS DCUTDTC, DCUTDESC).
    if let (true, Some(cut)) = (interim, cutoff) {
        resolved.design.cutoff = cut;
        resolved.design.interim = true;
    }
    let design = &resolved.design;
    // A data cutoff excludes subjects randomised after it (an interim snapshot can fall during
    // enrollment); subject IDs are those of the full study.
    let mut timelines = timelines;
    if let Some(cut) = cutoff {
        let keep: Vec<bool> = subjects.iter().map(|s| s.randomization <= cut).collect();
        if !keep.iter().any(|k| *k) {
            return Err(Error::Invalid(format!(
                "no subject is randomised by the data cutoff {cut}"
            )));
        }
        let mut k = keep.iter();
        subjects.retain(|_| *k.next().expect("one flag per subject"));
        let mut k = keep.iter();
        timelines.retain(|_| *k.next().expect("one flag per timeline"));
    }
    let ctx = sdtm::Context::new(design, &subjects, &timelines);
    let sdtm = sdtm::emit(&ctx).map_err(|e| internal(&e))?;
    let adam =
        adam::derive(&sdtm, &adam::AdslInputs { strata, cutoff }).map_err(|e| internal(&e))?;
    // Verification levels 0–2 are part of generation: a failure is a trialforge bug.
    let report = check::run(&sdtm, &adam);
    if !report.passed() {
        return Err(Error::Internal(format!(
            "generated data failed verification:\n{}",
            report.summary()
        )));
    }
    let truth = truth::manifest(&resolved, archetype, &sdtm, &adam);
    // Defects go into the verified datasets, after the truth manifest has described them.
    let (mut sdtm, mut adam) = (sdtm, adam);
    let defects =
        opts.defects.as_ref().filter(|d| !d.is_empty()).map(|spec| {
            defects::inject(design.stream().label("defects"), spec, &mut sdtm, &mut adam)
        });
    let last_contact = timelines
        .iter()
        .map(|t| t.last_contact)
        .max()
        .unwrap_or(design.cutoff);

    Ok(Study {
        studyid: design.studyid.clone(),
        archetype: archetype.to_string(),
        seed: design.seed,
        cutoff: cutoff.unwrap_or(last_contact),
        title: design.info.title.clone(),
        controls: resolved.controls.clone(),
        sdtm,
        adam,
        truth,
        defects,
    })
}

/// `USUBJID` → randomised stratum, as the randomisation system records it.
fn strata_of(design: &model::Design, subjects: &[model::Subject]) -> HashMap<String, String> {
    subjects
        .iter()
        .map(|s| {
            (
                s.usubjid.clone(),
                design.population.strata.levels[s.stratum].0.clone(),
            )
        })
        .collect()
}

/// Re-runs verification levels 0–2 on a study (for example after user post-processing, or to
/// find its injected defects). Every finding is kept.
pub fn verify(study: &Study) -> check::Report {
    check::run_all(&study.sdtm, &study.adam)
}

impl Study {
    /// The define.xml for the SDTM or ADaM datasets.
    pub fn define_xml(&self, standard: io::define::Standard) -> Result<String, Error> {
        let datasets = match standard {
            io::define::Standard::Sdtm => &self.sdtm,
            io::define::Standard::Adam => &self.adam,
        };
        let info = io::define::DefineInfo {
            studyid: &self.studyid,
            description: &self.title,
            created: self.cutoff,
            standard,
        };
        io::define::to_string(datasets, &info).map_err(|e| Error::Internal(e.to_string()))
    }

    /// Writes a submission-style package: `dir/sdtm/` and `dir/adam/`, each with one XPT file
    /// per dataset and a `define.xml`, plus `dir/truth.json` (and `dir/defects.json` when the
    /// study has defects). Returns the paths written.
    pub fn write_submission(&self, dir: &Path) -> Result<Vec<std::path::PathBuf>, Error> {
        let mut written = Vec::new();
        for (sub, sets, standard) in [
            ("sdtm", &self.sdtm, io::define::Standard::Sdtm),
            ("adam", &self.adam, io::define::Standard::Adam),
        ] {
            let folder = dir.join(sub);
            std::fs::create_dir_all(&folder).map_err(Error::Io)?;
            for ds in sets {
                let path = folder.join(format!("{}.xpt", ds.name.to_ascii_lowercase()));
                io::xpt::write_file(ds, &path, self.cutoff).map_err(|e| match e {
                    io::xpt::XptError::Io(e) => Error::Io(e),
                    io::xpt::XptError::Invalid(m) => Error::Internal(m),
                })?;
                written.push(path);
            }
            let define = folder.join("define.xml");
            std::fs::write(&define, self.define_xml(standard)?).map_err(Error::Io)?;
            written.push(define);
        }
        written.extend(self.write_manifests(dir)?);
        Ok(written)
    }

    /// The dataset called `name` (case-insensitive), from SDTM or ADaM.
    pub fn dataset(&self, name: &str) -> Option<&Dataset> {
        self.sdtm
            .iter()
            .chain(&self.adam)
            .find(|d| d.name.eq_ignore_ascii_case(name))
    }

    /// Writes every dataset to `dir` (created if needed) in each of `formats`, plus
    /// `truth.json` (and `defects.json` when the study has defects). Files are named by
    /// lower-case dataset name, e.g. `dm.xpt`. Returns the paths written.
    pub fn write(&self, dir: &Path, formats: &[&str]) -> Result<Vec<std::path::PathBuf>, Error> {
        if let Some(bad) = formats.iter().find(|f| !FORMATS.contains(f)) {
            return Err(Error::Invalid(format!(
                "unknown format {bad:?}; available: {}",
                FORMATS.join(", ")
            )));
        }
        std::fs::create_dir_all(dir).map_err(Error::Io)?;
        let mut written = Vec::new();
        let is_adam = |ds: &Dataset| self.adam.iter().any(|a| a.name == ds.name);
        for ds in self.sdtm.iter().chain(&self.adam) {
            for &format in formats {
                let ext = if format == "dataset-json" {
                    "json"
                } else {
                    format
                };
                let path = dir.join(format!("{}.{ext}", ds.name.to_ascii_lowercase()));
                match format {
                    "xpt" => io::xpt::write_file(ds, &path, self.cutoff).map_err(|e| match e {
                        io::xpt::XptError::Io(e) => Error::Io(e),
                        io::xpt::XptError::Invalid(m) => Error::Internal(m),
                    })?,
                    "csv" => io::csv::write_file(ds, &path).map_err(Error::Io)?,
                    "dataset-json" => {
                        let info = io::define::DefineInfo {
                            studyid: &self.studyid,
                            description: &self.title,
                            created: self.cutoff,
                            standard: if is_adam(ds) {
                                io::define::Standard::Adam
                            } else {
                                io::define::Standard::Sdtm
                            },
                        };
                        std::fs::write(&path, io::dataset_json::to_string(ds, &info, VERSION))
                            .map_err(Error::Io)?;
                    }
                    #[cfg(feature = "parquet")]
                    "parquet" => io::parquet::write_file(ds, &path).map_err(|e| match e {
                        io::parquet::ParquetError::Io(e) => Error::Io(e),
                        io::parquet::ParquetError::Encode(m) => Error::Internal(m),
                    })?,
                    _ => unreachable!("formats validated above"),
                }
                written.push(path);
            }
        }
        written.extend(self.write_manifests(dir)?);
        Ok(written)
    }

    /// The defect manifest as written to `defects.json`: the request (the profile, if one was
    /// named, and the amount per type), the seed and every defect, with zero-based rows.
    pub fn defects_json(&self) -> Option<serde_json::Value> {
        let m = self.defects.as_ref()?;
        Some(serde_json::json!({
            "studyid": self.studyid,
            "seed": self.seed,
            "rows": "zero-based",
            "profile": m.spec.profile,
            "spec": m.spec,
            "defects": m.defects,
        }))
    }

    /// The injected defects as a table, one row per defect; `expected` lists the expected
    /// findings as `RULE DATASET row`, separated by `; `. Rows are numbered from `first_row`
    /// (0 in Rust and Python, 1 in R).
    pub fn defects_table(&self, first_row: usize) -> Option<Dataset> {
        let m = self.defects.as_ref()?;
        let number = |r: usize| (r + first_row) as i64;
        let char_col = |name: &str, label: &str, f: &dyn Fn(&defects::Defect) -> Option<String>| {
            Column::new(name, label, Values::Char(m.defects.iter().map(f).collect()))
        };
        let columns = vec![
            char_col("id", "Defect type", &|d| Some(d.id.to_string())),
            char_col("dataset", "Dataset", &|d| Some(d.dataset.clone())),
            Column::new(
                "row",
                "Row",
                Values::Int(m.defects.iter().map(|d| Some(number(d.row))).collect()),
            ),
            Column::new(
                "paired_row",
                "Paired row",
                Values::Int(m.defects.iter().map(|d| d.paired_row.map(number)).collect()),
            ),
            char_col("usubjid", "Unique Subject Identifier", &|d| {
                d.usubjid.clone()
            }),
            char_col("key", "Record key", &|d| d.key.clone()),
            char_col("variable", "Variable", &|d| Some(d.variable.clone())),
            char_col("original", "Original value", &|d| d.original.clone()),
            char_col("injected", "Injected value", &|d| d.injected.clone()),
            char_col("expected", "Expected findings", &|d| {
                Some(
                    d.expected
                        .iter()
                        .map(|e| match e.row {
                            Some(r) => format!("{} {} {}", e.rule, e.dataset, number(r)),
                            None => format!("{} {}", e.rule, e.dataset),
                        })
                        .collect::<Vec<_>>()
                        .join("; "),
                )
            }),
        ];
        Some(Dataset::new("DEFECTS", "Injected defects", columns).expect("equal-length columns"))
    }

    /// Scores a validator's findings against the injected defects (precision and recall, per
    /// defect type). Fails if the study has no defects.
    pub fn score(&self, findings: &[defects::Observed<'_>]) -> Result<defects::Score, Error> {
        let m = self.defects.as_ref().ok_or_else(|| {
            Error::Invalid("this study has no injected defects; generate it with `defects`".into())
        })?;
        Ok(defects::score(&m.defects, findings))
    }

    /// Writes `truth.json` and, when the study has defects, `defects.json` to `dir`.
    fn write_manifests(&self, dir: &Path) -> Result<Vec<std::path::PathBuf>, Error> {
        let mut written = Vec::new();
        let files = [
            ("truth.json", Some(self.truth.clone())),
            ("defects.json", self.defects_json()),
        ];
        for (name, value) in files {
            let Some(value) = value else { continue };
            let path = dir.join(name);
            std::fs::write(
                &path,
                serde_json::to_string_pretty(&value).expect("manifests serialize") + "\n",
            )
            .map_err(Error::Io)?;
            written.push(path);
        }
        Ok(written)
    }
}
