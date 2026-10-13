//! Python bindings for trialforge (`trialforge._core`). Keep this layer thin: conversion and
//! error mapping only. Generation logic belongs in the Rust crates so R gets it too.

use std::collections::HashMap;
use std::path::PathBuf;
use std::sync::Arc;

use arrow_array::ffi_stream::FFI_ArrowArrayStream;
use arrow_array::{
    ArrayRef, Date32Array, Float64Array, Int64Array, RecordBatch, RecordBatchIterator, StringArray,
};
use arrow_schema::{DataType, Field, Schema};
use pyo3::exceptions::{PyOSError, PyRuntimeError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::PyCapsule;
use trialforge::{Dataset, StudyOptions, Values};

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

/// Generates a study. `seed=None` draws a fresh seed (recorded on the result).
#[pyfunction]
#[pyo3(signature = (archetype, n, seed=None, studyid=None))]
fn study(
    py: Python<'_>,
    archetype: &str,
    n: usize,
    seed: Option<u64>,
    studyid: Option<String>,
) -> PyResult<Study> {
    let opts = StudyOptions {
        n,
        seed: seed.unwrap_or_else(trialforge::random_seed),
        studyid,
    };
    let inner = py
        .detach(|| trialforge::study(archetype, &opts))
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

    /// The data cutoff date as days since 1970-01-01.
    #[getter]
    fn cutoff_days(&self) -> i32 {
        self.inner.cutoff.days()
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
        to_arrow(ds)
    }

    /// Writes every dataset to `path` in each of `formats`; returns the paths written.
    fn write(&self, py: Python<'_>, path: PathBuf, formats: Vec<String>) -> PyResult<Vec<PathBuf>> {
        let formats: Vec<&str> = formats.iter().map(String::as_str).collect();
        py.detach(|| self.inner.write(&path, &formats))
            .map_err(to_py_err)
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

/// A dataset as an Arrow table. Labels travel as `label` metadata on each field and on the
/// schema; dates are `Date32`.
fn to_arrow(ds: &Dataset) -> PyResult<ArrowTable> {
    let mut fields = Vec::with_capacity(ds.columns().len());
    let mut arrays: Vec<ArrayRef> = Vec::with_capacity(ds.columns().len());
    for c in ds.columns() {
        let (dtype, array): (DataType, ArrayRef) = match &c.values {
            Values::Char(v) => (
                DataType::Utf8,
                Arc::new(StringArray::from_iter(v.iter().map(Option::as_deref))),
            ),
            Values::Num(v) => (DataType::Float64, Arc::new(Float64Array::from(v.clone()))),
            Values::Int(v) => (DataType::Int64, Arc::new(Int64Array::from(v.clone()))),
            Values::Date(v) => (
                DataType::Date32,
                Arc::new(Date32Array::from(
                    v.iter().map(|d| d.map(|d| d.days())).collect::<Vec<_>>(),
                )),
            ),
        };
        let metadata = HashMap::from([("label".to_string(), c.label.clone())]);
        fields.push(Field::new(&c.name, dtype, true).with_metadata(metadata));
        arrays.push(array);
    }
    let schema = Arc::new(Schema::new_with_metadata(
        fields,
        HashMap::from([
            ("name".to_string(), ds.name.clone()),
            ("label".to_string(), ds.label.clone()),
        ]),
    ));
    let batch =
        RecordBatch::try_new(schema, arrays).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
    Ok(ArrowTable { batch })
}

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(stream_seed, m)?)?;
    m.add_function(wrap_pyfunction!(archetypes, m)?)?;
    m.add_function(wrap_pyfunction!(study, m)?)?;
    m.add_class::<Study>()?;
    m.add_class::<ArrowTable>()?;
    m.add("MAX_SEED", trialforge::MAX_SEED)?;
    Ok(())
}
