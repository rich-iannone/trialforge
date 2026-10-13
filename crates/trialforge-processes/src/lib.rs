//! Per-subject process engines that build each subject's [`timeline::Timeline`]: dosing, visits,
//! adverse events, medications, medical history, labs and vital signs, time-to-event and more.
//!
//! - [`fixed`]: fixed-duration studies with a visit schedule (e.g. `parallel_continuous`).
//! - [`oncology`]: an event-driven solid-tumour course (SIMONC01-like).
//! - [`reference`]: the reference data both draw on.

pub mod fixed;
pub mod oncology;
pub mod reference;
pub mod timeline;
