//! Foundations shared by every other crate: deterministic RNG streams (see [`rng`]),
//! distributions and copulas, the study calendar and ISO 8601 (including partial) dates.

pub mod rng;

pub use rng::Stream;
