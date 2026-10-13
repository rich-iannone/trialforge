//! Python bindings for trialforge (`trialforge._core`). Keep this layer thin: conversion and
//! error mapping only. Generation logic belongs in the Rust crates so R gets it too.

use pyo3::prelude::*;

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

#[pymodule]
fn _core(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(core_version, m)?)?;
    m.add_function(wrap_pyfunction!(stream_seed, m)?)?;
    Ok(())
}
