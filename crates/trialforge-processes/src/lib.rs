//! Per-subject process engines that build each subject's [`timeline::Timeline`]: dosing, visits,
//! adverse events, medications, medical history, labs and vital signs, time-to-event and more.
//!
//! - [`adas`]: ADAS-Cog item responses from a latent cognitive trajectory (FT findings).
//! - [`crossover`]: the fixed-duration course split into crossover periods and washouts.
//! - [`fixed`]: fixed-duration studies with a visit schedule (e.g. `parallel_continuous`).
//! - [`outcomes`]: an event-driven cardiovascular outcomes course (CE events, competing death).
//! - [`oncology`]: an event-driven solid-tumour course (SIMONC01-like).
//! - [`pk`]: timed dosing and plasma concentrations layered on another course (phase 1).
//! - [`reference`]: the reference data they draw on.

pub mod adas;
pub mod crossover;
pub mod fixed;
pub mod oncology;
pub mod outcomes;
pub mod pk;
pub mod reference;
pub mod timeline;
