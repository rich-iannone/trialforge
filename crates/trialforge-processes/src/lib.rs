//! Per-subject process engines that build each subject's timeline: dosing, visits, adverse
//! events, con meds, labs, vitals, time-to-event, tumour, PK and more.
//!
//! The spike has one engine, [`oncology`]: a simplified SIMONC01 course (progression, death,
//! withdrawal, discontinuation for AEs, Q-cycle dosing) without lesion-level tumour data.

pub mod oncology;
