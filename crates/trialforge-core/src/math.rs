//! Portable maths for output paths. These go through the pure-Rust `libm` crate, so results are
//! bit-identical on every platform. Use them instead of `f64::ln`, `f64::exp` and friends
//! anywhere the result can reach generated data.

/// Natural logarithm.
pub fn ln(x: f64) -> f64 {
    libm::log(x)
}

/// Exponential function.
pub fn exp(x: f64) -> f64 {
    libm::exp(x)
}

/// Square root.
pub fn sqrt(x: f64) -> f64 {
    libm::sqrt(x)
}
