//! Python bindings for trialforge (`trialforge._core`). Keep this layer thin: conversion and
//! error mapping only. Generation logic belongs in the Rust crates so R gets it too.

use std::path::PathBuf;

use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{RecordBatch, RecordBatchIterator};
use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use trialforge::StudyOptions;

/// The Rust core's version string.
#[pyfunction]
fn core_version() -> &'static str {
    trialforge::VERSION
}

/// Seed of the RNG stream at `path` under `seed`. Integers in `path` are index components and
/// strings are label components, e.g. `stream_seed(1, ["subject", 17, "ae"])`. Used by the
/// cross-language determinism tests and by user hooks.
#[pyfunction]
fn stream_seed(seed: u64, path: Vec<StreamPart>) -> u64 {
    let mut s = trialforge::core::Stream::root(seed);
    for part in path {
        s = match part {
            StreamPart::Index(i) => s.index(i),
            StreamPart::Label(l) => s.label(&l),
        };
    }
    s.seed()
}

#[derive(FromPyObject)]
enum StreamPart {
    Index(u64),
    Label(String),
}

/// (name, description) for every archetype.
#[pyfunction]
fn archetypes() -> Vec<(&'static str, &'static str)> {
    trialforge::ARCHETYPES.to_vec()
}

/// The documented controls of an archetype, as a JSON string.
#[pyfunction]
fn controls(archetype: &str) -> PyResult<String> {
    trialforge::controls(archetype)
        .map(|v| v.to_string())
        .map_err(to_py_err)
}

/// Output formats `Study.write` understands.
#[pyfunction]
fn formats() -> Vec<&'static str> {
    trialforge::FORMATS.to_vec()
}

/// Generates a study. `seed=None` draws a fresh seed (recorded on the result); `controls` is a
/// JSON object string of archetype controls.
#[pyfunction]
#[pyo3(signature = (archetype, n, seed=None, studyid=None, controls=None))]
fn study(
    py: Python<'_>,
    archetype: &str,
    n: usize,
    seed: Option<u64>,
    studyid: Option<String>,
    controls: Option<&str>,
) -> PyResult<Study> {
    let controls = match controls {
        Some(c) => serde_json::from_str(c)
            .map_err(|e| PyValueError::new_err(format!("controls are not valid JSON: {e}")))?,
        None => serde_json::Value::Null,
    };
    let opts = StudyOptions {
        n,
        seed: seed.unwrap_or_else(trialforge::random_seed),
        studyid,
        controls,
    };
    let inner = py
        .detach(|| trialforge::study(archetype, &opts))
        .map_err(to_py_err)?;
    Ok(Study { inner })
}

/// Generates a study from a JSON or YAML study spec (the text, not a path).
#[pyfunction]
fn study_from_spec(py: Python<'_>, text: &str) -> PyResult<Study> {
    let text = text.to_string();
    let (inner, _) = py
        .detach(|| trialforge::study_from_spec(&text))
        .map_err(to_py_err)?;
    Ok(Study { inner })
}

fn to_py_err(e: trialforge::Error) -> PyErr {
    match e {
        trialforge::Error::Invalid(m) => PyValueError::new_err(m),
        e @ trialforge::Error::Internal(_) => PyRuntimeError::new_err(e.to_string()),
        trialforge::Error::Io(e) => PyOSError::new_err(e.to_string()),
    }
}

/// A generated study, held in Rust.
#[pyclass(frozen, module = "trialforge._core")]
struct Study {
    inner: trialforge::Study,
}

#[pymethods]
impl Study {
    #[getter]
    fn studyid(&self) -> &str {
        &self.inner.studyid
    }

    #[getter]
    fn archetype(&self) -> &str {
        &self.inner.archetype
    }

    #[getter]
    fn seed(&self) -> u64 {
        self.inner.seed
    }

    /// The data cutoff (or last contact) date as days since 1970-01-01.
    #[getter]
    fn cutoff_days(&self) -> i32 {
        self.inner.cutoff.days()
    }

    /// The controls used, as a JSON string.
    #[getter]
    fn controls_json(&self) -> String {
        self.inner.controls.to_string()
    }

    /// The truth manifest, as a JSON string.
    #[getter]
    fn truth_json(&self) -> String {
        self.inner.truth.to_string()
    }

    #[getter]
    fn sdtm_names(&self) -> Vec<String> {
        self.inner.sdtm.iter().map(|d| d.name.clone()).collect()
    }

    #[getter]
    fn adam_names(&self) -> Vec<String> {
        self.inner.adam.iter().map(|d| d.name.clone()).collect()
    }

    /// The dataset called `name` as an Arrow table (Arrow PyCapsule interface).
    fn table(&self, name: &str) -> PyResult<ArrowTable> {
        let ds = self
            .inner
            .dataset(name)
            .ok_or_else(|| PyValueError::new_err(format!("no dataset named {name:?}")))?;
        let batch = trialforge::io::arrow::to_record_batch(ds)
            .map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        Ok(ArrowTable { batch })
    }

    /// Writes every dataset to `path` in each of `formats`, plus `truth.json`; returns the
    /// paths written.
    fn write(&self, py: Python<'_>, path: PathBuf, formats: Vec<String>) -> PyResult<Vec<PathBuf>> {
        let formats: Vec<&str> = formats.iter().map(String::as_str).collect();
        py.detach(|| self.inner.write(&path, &formats))
            .map_err(to_py_err)
    }

    /// Writes `path/sdtm/` and `path/adam/` (XPT + define.xml each) and `path/truth.json`;
    /// returns the paths written.
    fn write_submission(&self, py: Python<'_>, path: PathBuf) -> PyResult<Vec<PathBuf>> {
        py.detach(|| self.inner.write_submission(&path))
            .map_err(to_py_err)
    }

    /// The define.xml for `standard` ("sdtm" or "adam").
    fn define_xml(&self, standard: &str) -> PyResult<String> {
        let std = match standard {
            "sdtm" => trialforge::io::define::Standard::Sdtm,
            "adam" => trialforge::io::define::Standard::Adam,
            other => {
                return Err(PyValueError::new_err(format!(
                    "`standard` must be \"sdtm\" or \"adam\", not {other:?}"
                )));
            }
        };
        self.inner.define_xml(std).map_err(to_py_err)
    }

    /// Re-runs verification levels 0–2: (passed, [(level, rule, dataset, variable, row, message)]).
    #[allow(clippy::type_complexity)]
    fn verify(
        &self,
    ) -> (
        bool,
        Vec<(
            String,
            String,
            String,
            Option<String>,
            Option<usize>,
            String,
        )>,
    ) {
        let r = trialforge::verify(&self.inner);
        let findings = r
            .findings
            .iter()
            .map(|f| {
                (
                    format!("{:?}", f.level),
                    f.rule.to_string(),
                    f.dataset.clone(),
                    f.variable.clone(),
                    f.row,
                    f.message.clone(),
                )
            })
            .collect();
        (r.passed(), findings)
    }
}

/// An Arrow table exported through the Arrow PyCapsule interface (`__arrow_c_stream__`), which
/// polars, pyarrow, pandas (via pyarrow), DuckDB and others consume without copying.
#[pyclass(frozen, module = "trialforge._core")]
struct ArrowTable {
    batch: RecordBatch,
}

#[pymethods]
impl ArrowTable {
    #[getter]
    fn num_rows(&self) -> usize {
        self.batch.num_rows()
    }

    /// Exports the table as an `ArrowArrayStream` capsule. `requested_schema` is accepted for
    /// protocol compatibility and ignored: the table has a single fixed schema.
    #[pyo3(signature = (requested_schema=None))]
    fn __arrow_c_stream__<'py>(
        &self,
        py: Python<'py>,
        requested_schema: Option<Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, PyCapsule>> {
        let _ = requested_schema;
        let reader = RecordBatchIterator::new([Ok(self.batch.clone())], self.batch.schema());
        let stream = FFI_ArrowArrayStream::new(Box::new(reader));
        // The capsule owns the stream. A consumer moves it out (leaving a released stream
        // behind); otherwise dropping it with the capsule releases it.
        PyCapsule::new_with_value(py, stream, c"arrow_array_stream")
    }
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(stream_seed, m)?)?;
    m.add_function(wrap_pyfunction!(archetypes, m)?)?;
    m.add_function(wrap_pyfunction!(study, m)?)?;
    m.add_function(wrap_pyfunction!(controls, m)?)?;
    m.add_function(wrap_pyfunction!(study_from_spec, m)?)?;
    m.add_function(wrap_pyfunction!(formats, m)?)?;
    m.add_class::<Study>()?;
    m.add_class::<ArrowTable>()?;
    m.add("MAX_SEED", trialforge::MAX_SEED)?;
    Ok(())
}
