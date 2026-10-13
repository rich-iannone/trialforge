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

// Macro to generate exports. This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod trialforge;
    fn core_version;
    fn stream_seed;
}
