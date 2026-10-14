//! Foundations shared by every other crate: deterministic RNG streams ([`rng`]), portable
//! random draws ([`dist`]), portable maths ([`math`]), the study calendar ([`date`]) and
//! in-memory datasets ([`table`]) and small survival statistics ([`stats`]).

pub mod date;
pub mod dist;
pub mod math;
pub mod nca;
pub mod rng;
pub mod stats;
pub mod table;

pub use date::{Date, DateTime, study_day};
pub use dist::Draws;
pub use rng::Stream;
pub use table::{Column, Dataset, TableError, Values};
