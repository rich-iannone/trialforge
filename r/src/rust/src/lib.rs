//! R bindings for trialforge. Keep this layer thin: conversion and error mapping only.
//! Generation logic belongs in the Rust crates so Python gets it too.

use extendr_api::prelude::*;
use extendr_api::{Error, Result};

/// The Rust core's version string.
/// @export
#[extendr]
fn core_version() -> &'static str {
    trialforge::VERSION
}

/// Seed of an RNG stream
///
/// Returns the seed of the RNG stream at `path` under `seed`. Whole numbers in `path` are
/// index components and strings are label components, e.g.
/// `stream_seed(1, list("subject", 17, "ae"))`. Used by the cross-language determinism tests
/// and by user hooks.
///
/// R has no unsigned 64-bit integer type, so the seed is returned as a decimal string. It
/// matches the integer returned by the Python package exactly.
///
/// @param seed A whole number in `[0, 2^53]`.
/// @param path A list of strings and whole numbers.
/// @return A string holding the unsigned 64-bit stream seed.
/// @export
#[extendr]
fn stream_seed(seed: f64, path: List) -> Result<String> {
    let mut s = trialforge::core::Stream::root(as_u64(seed, "`seed`")?);
    for (i, (_, part)) in path.iter().enumerate() {
        s = if let Some(label) = part.as_str() {
            s.label(label)
        } else if let Some(n) = part.as_integer() {
            s.index(as_u64(f64::from(n), &format!("`path[[{}]]`", i + 1))?)
        } else if let Some(x) = part.as_real() {
            s.index(as_u64(x, &format!("`path[[{}]]`", i + 1))?)
        } else {
            return Err(Error::Other(format!(
                "`path[[{}]]` must be a string or a whole number",
                i + 1
            )));
        };
    }
    Ok(s.seed().to_string())
}

/// Convert an R number to `u64`, refusing values doubles cannot represent exactly.
fn as_u64(x: f64, what: &str) -> Result<u64> {
    const MAX_EXACT: f64 = 9_007_199_254_740_992.0; // 2^53
    if x.is_finite() && x >= 0.0 && x.fract() == 0.0 && x <= MAX_EXACT {
        Ok(x as u64)
    } else {
        Err(Error::Other(format!(
            "{what} must be a whole number in [0, 2^53]"
        )))
    }
}

/// (name, description) pairs for every archetype, as a two-column list.
/// @noRd
#[extendr]
fn archetypes_impl() -> List {
    let names: Vec<&str> = trialforge::ARCHETYPES.iter().map(|(n, _)| *n).collect();
    let descriptions: Vec<&str> = trialforge::ARCHETYPES.iter().map(|(_, d)| *d).collect();
    list!(name = names, description = descriptions)
}

/// The documented controls of an archetype, as a JSON string.
/// @noRd
#[extendr]
fn controls_impl(archetype: &str) -> Result<String> {
    trialforge::controls(archetype)
        .map(|v| v.to_string())
        .map_err(|e| Error::Other(e.to_string()))
}

/// Output formats `write_study()` supports.
/// @noRd
#[extendr]
fn formats_impl() -> Vec<String> {
    trialforge::FORMATS.iter().map(|s| s.to_string()).collect()
}

/// Generates a study. `seed`, `studyid`, `controls` and `defects` are `NA` when not supplied;
/// `controls` is a JSON object string and `defects` a JSON string (`"none"` or an object of
/// defect id to count or rate).
/// @noRd
#[extendr]
fn study_impl(
    archetype: &str,
    n: f64,
    seed: Option<f64>,
    studyid: Option<String>,
    controls: Option<String>,
    defects: Option<String>,
) -> Result<RStudy> {
    let n = as_u64(n, "`n`")? as usize;
    let seed = match seed {
        Some(s) => as_u64(s, "`seed`")?,
        None => trialforge::random_seed(),
    };
    let controls = match controls {
        Some(c) => serde_json::from_str(&c)
            .map_err(|e| Error::Other(format!("controls are not valid JSON: {e}")))?,
        None => serde_json::Value::Null,
    };
    let defects = match defects {
        Some(d) => {
            let v = serde_json::from_str(&d)
                .map_err(|e| Error::Other(format!("defects are not valid JSON: {e}")))?;
            Some(trialforge::defects::DefectSpec::from_value(&v).map_err(Error::Other)?)
        }
        None => None,
    };
    let opts = trialforge::StudyOptions {
        n,
        seed,
        studyid,
        controls,
        defects,
    };
    trialforge::study(archetype, &opts)
        .map(|inner| RStudy { inner })
        .map_err(|e| Error::Other(e.to_string()))
}

/// Generates a study from JSON or YAML spec text.
/// @noRd
#[extendr]
fn study_from_spec_impl(text: &str) -> Result<RStudy> {
    trialforge::study_from_spec(text)
        .map(|(inner, _)| RStudy { inner })
        .map_err(|e| Error::Other(e.to_string()))
}

/// A generated study, held in Rust.
#[extendr]
struct RStudy {
    inner: trialforge::Study,
}

/// @noRd
#[extendr]
impl RStudy {
    fn studyid(&self) -> String {
        self.inner.studyid.clone()
    }

    fn archetype(&self) -> String {
        self.inner.archetype.clone()
    }

    /// Seeds are at most 2^53, so a double holds them exactly.
    fn seed(&self) -> f64 {
        self.inner.seed as f64
    }

    fn cutoff_days(&self) -> i32 {
        self.inner.cutoff.days()
    }

    fn controls_json(&self) -> String {
        self.inner.controls.to_string()
    }

    fn truth_json(&self) -> String {
        self.inner.truth.to_string()
    }

    /// The injected defects as a data.frame (rows numbered from 1), or `NULL` for a clean
    /// study.
    fn defects_table(&self) -> Result<Robj> {
        match self.inner.defects_table(1) {
            Some(ds) => data_frame(&ds),
            None => Ok(Robj::from(())),
        }
    }

    /// Scores findings (parallel vectors; `row` is 1-based, `NA` for none) against the
    /// injected defects. Returns the score as a JSON string, whose `missed` and `unmatched`
    /// indices are 0-based.
    fn score(
        &self,
        level: Vec<String>,
        rule: Vec<String>,
        dataset: Vec<String>,
        row: Integers,
    ) -> Result<String> {
        let rows: Vec<Option<usize>> = row
            .iter()
            .map(|r| match Option::<i32>::from(r) {
                None => Ok(None),
                Some(r) => usize::try_from(r - 1)
                    .map(Some)
                    .map_err(|_| Error::Other("finding rows must be positive".into())),
            })
            .collect::<Result<_>>()?;
        if [rule.len(), dataset.len(), rows.len()] != [level.len(); 3] {
            return Err(Error::Other(
                "`level`, `rule`, `dataset` and `row` must have equal lengths".into(),
            ));
        }
        let observed: Vec<trialforge::defects::Observed<'_>> = (0..level.len())
            .map(|i| trialforge::defects::Observed {
                level: &level[i],
                rule: &rule[i],
                dataset: &dataset[i],
                row: rows[i],
            })
            .collect();
        let score = self
            .inner
            .score(&observed)
            .map_err(|e| Error::Other(e.to_string()))?;
        Ok(serde_json::to_string(&score).expect("scores serialize"))
    }

    fn write_submission(&self, path: &str) -> Result<Vec<String>> {
        self.inner
            .write_submission(std::path::Path::new(path))
            .map(|paths| {
                paths
                    .iter()
                    .map(|p| p.to_string_lossy().into_owned())
                    .collect()
            })
            .map_err(|e| Error::Other(e.to_string()))
    }

    fn define_xml(&self, standard: &str) -> Result<String> {
        let std = match standard {
            "sdtm" => trialforge::io::define::Standard::Sdtm,
            "adam" => trialforge::io::define::Standard::Adam,
            other => {
                return Err(Error::Other(format!(
                    "`standard` must be \"sdtm\" or \"adam\", not {other:?}"
                )));
            }
        };
        self.inner
            .define_xml(std)
            .map_err(|e| Error::Other(e.to_string()))
    }

    /// Verification levels 0–2: a list with `passed` and a `findings` data frame.
    fn verify(&self) -> Result<Robj> {
        let r = trialforge::verify(&self.inner);
        let f = &r.findings;
        let col = |g: &dyn Fn(&trialforge::check::Finding) -> Option<String>| {
            Strings::from_values(
                f.iter()
                    .map(|x| g(x).map_or_else(Rstr::na, |s| Rstr::from(s.as_str()))),
            )
        };
        let row = Integers::from_values(
            f.iter()
                .map(|x| x.row.map_or_else(Rint::na, |r| Rint::from(r as i32 + 1))),
        );
        let mut findings: Robj = List::from_names_and_values(
            ["level", "rule", "dataset", "variable", "row", "message"],
            [
                col(&|x| Some(format!("{:?}", x.level))).into(),
                col(&|x| Some(x.rule.to_string())).into(),
                col(&|x| Some(x.dataset.clone())).into(),
                col(&|x| x.variable.clone()).into(),
                Robj::from(row),
                col(&|x| Some(x.message.clone())).into(),
            ],
        )?
        .into();
        findings.set_class(["data.frame"])?;
        findings.set_attrib(
            row_names_symbol(),
            [Rint::na(), Rint::from(-(f.len() as i32))]
                .into_iter()
                .collect::<Integers>(),
        )?;
        Ok(list!(passed = r.passed(), findings = findings).into())
    }

    fn sdtm_names(&self) -> Vec<String> {
        self.inner.sdtm.iter().map(|d| d.name.clone()).collect()
    }

    fn adam_names(&self) -> Vec<String> {
        self.inner.adam.iter().map(|d| d.name.clone()).collect()
    }

    /// The dataset called `name` as a data.frame with haven-style `label` attributes.
    fn table(&self, name: &str) -> Result<Robj> {
        let ds = self
            .inner
            .dataset(name)
            .ok_or_else(|| Error::Other(format!("no dataset named {name:?}")))?;
        data_frame(ds)
    }

    fn write(&self, path: &str, formats: Vec<String>) -> Result<Vec<String>> {
        let formats: Vec<&str> = formats.iter().map(String::as_str).collect();
        self.inner
            .write(std::path::Path::new(path), &formats)
            .map(|paths| {
                paths
                    .iter()
                    .map(|p| p.to_string_lossy().into_owned())
                    .collect()
            })
            .map_err(|e| Error::Other(e.to_string()))
    }
}

/// A dataset as a data.frame. Character missing is `NA`; dates are `Date`; date-times are
/// `POSIXct` in UTC (the study clock has no time zone; UTC keeps the clock time); whole numbers are
/// integer; every column and the frame carry a `label` attribute.
fn data_frame(ds: &trialforge::Dataset) -> Result<Robj> {
    let mut names = Vec::with_capacity(ds.columns().len());
    let mut values = Vec::with_capacity(ds.columns().len());
    for c in ds.columns() {
        let mut col: Robj = match &c.values {
            trialforge::Values::Char(v) => Strings::from_values(
                v.iter()
                    .map(|x| x.as_deref().map_or_else(Rstr::na, Rstr::from)),
            )
            .into(),
            trialforge::Values::Num(v) => {
                Doubles::from_values(v.iter().map(|x| x.map_or_else(Rfloat::na, Rfloat::from)))
                    .into()
            }
            trialforge::Values::Int(v) => Integers::from_values(v.iter().map(|x| match x {
                Some(n) => i32::try_from(*n).map_or_else(|_| Rint::na(), Rint::from),
                None => Rint::na(),
            }))
            .into(),
            trialforge::Values::Date(v) => {
                let mut d: Robj = Doubles::from_values(
                    v.iter()
                        .map(|x| x.map_or_else(Rfloat::na, |d| Rfloat::from(f64::from(d.days())))),
                )
                .into();
                d.set_class(["Date"])?;
                d
            }
            trialforge::Values::DateTime(v) => {
                let mut d: Robj = Doubles::from_values(
                    v.iter()
                        .map(|x| x.map_or_else(Rfloat::na, |d| Rfloat::from(d.seconds() as f64))),
                )
                .into();
                d.set_class(["POSIXct", "POSIXt"])?;
                d.set_attrib("tzone", "UTC")?;
                d
            }
        };
        col.set_attrib("label", c.label.as_str())?;
        names.push(c.name.as_str());
        values.push(col);
    }
    let mut df: Robj = List::from_names_and_values(names, values)?.into();
    df.set_class(["data.frame"])?;
    let n = i32::try_from(ds.n_rows()).map_err(|_| Error::Other("too many rows for R".into()))?;
    df.set_attrib(
        row_names_symbol(),
        [Rint::na(), Rint::from(-n)]
            .into_iter()
            .collect::<Integers>(),
    )?;
    df.set_attrib("label", ds.label.as_str())?;
    Ok(df)
}

// Macro to generate exports. This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod trialforge;
    fn core_version;
    fn stream_seed;
    fn archetypes_impl;
    fn controls_impl;
    fn formats_impl;
    fn study_impl;
    fn study_from_spec_impl;
    impl RStudy;
}
