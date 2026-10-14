//! Per-subject process engines that build each subject's [`timeline::Timeline`]: dosing, visits,
//! adverse events, medications, medical history, labs and vital signs, time-to-event and more.
//!
//! - [`adas`]: ADAS-Cog item responses from a latent cognitive trajectory (FT findings).
//! - [`crossover`]: the fixed-duration course split into crossover periods and washouts.
//! - [`escalation`]: 3+3, BOIN and mTPI-2 dose escalation and dose-limiting toxicities.
//! - [`exacerbation`]: respiratory exacerbations (CE, negative binomial counts) and FEV1 (RE).
//! - [`fixed`]: fixed-duration studies with a visit schedule (e.g. `parallel_continuous`).
//! - [`ophtha`]: ophthalmic examinations of both eyes and intravitreal injections (OE).
//! - [`growth`]: paediatric birth dates and growth (height, weight, BMI) on CDC charts.
//! - [`outcomes`]: an event-driven cardiovascular outcomes course (CE events, competing death).
//! - [`oncology`]: an event-driven solid-tumour course (SIMONC01-like).
//! - [`pk`]: timed dosing and plasma concentrations layered on another course (phase 1).
//! - [`vaccine`]: vaccinations, reactogenicity diaries (CE, FACE, VS) and antibody titres (IS).
//! - [`reference`]: the reference data they draw on.

pub mod adas;
pub mod crossover;
pub mod escalation;
pub mod exacerbation;
pub mod fixed;
pub mod growth;
pub mod oncology;
pub mod ophtha;
pub mod outcomes;
pub mod pk;
pub mod reference;
pub mod timeline;
pub mod vaccine;
