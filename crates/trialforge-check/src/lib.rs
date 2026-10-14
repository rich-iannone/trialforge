//! Verification levels 0–2, run on every generation.
//!
//! - **L0 structural** ([`structural`]): metadata conformance, controlled terminology, keys,
//!   ISO 8601.
//! - **L1 referential & temporal** ([`referential`]): subjects, dates, visits, epochs,
//!   exposure, disposition and supplemental-qualifier links agree across domains; RS responses
//!   follow RECIST 1.1 from TR ([`onco`]); PC is consistent with EX ([`pk`]).
//! - **L2 derivational** ([`derivational`]): ADaM is recomputed from SDTM by an independent
//!   implementation (this crate does not depend on `trialforge-adam`) and compared, including
//!   admiralonco's response and time-to-event derivations.
//!
//! Every finding is an error: trialforge promises its output passes all three levels.

mod adas;
mod derivational;
mod escalation;
mod exacerbation;
mod growth;
mod lugano;
mod onco;
mod ophtha;
mod outcomes;
mod periods;
mod pk;
mod referential;
mod structural;
mod switching;
mod table;
mod vaccine;

use std::collections::BTreeMap;

use trialforge_core::Dataset;

/// A verification level.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Level {
    /// Structural.
    L0,
    /// Referential and temporal.
    L1,
    /// Derivational.
    L2,
}

/// One failed check.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Finding {
    /// Level.
    pub level: Level,
    /// Rule identifier, e.g. `L1.DEATH`.
    pub rule: &'static str,
    /// Dataset.
    pub dataset: String,
    /// Variable, if the finding is about one.
    pub variable: Option<String>,
    /// Zero-based row, if the finding is about one.
    pub row: Option<usize>,
    /// What is wrong.
    pub message: String,
}

impl std::fmt::Display for Finding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "[{}] {}", self.rule, self.dataset)?;
        if let Some(v) = &self.variable {
            write!(f, ".{v}")?;
        }
        if let Some(r) = self.row {
            write!(f, " row {}", r + 1)?;
        }
        write!(f, ": {}", self.message)
    }
}

/// The result of a verification run.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Report {
    /// Findings (at most [`MAX_PER_RULE`] per rule and dataset).
    pub findings: Vec<Finding>,
    /// Total failures per rule, including those not kept in `findings`.
    pub counts: BTreeMap<&'static str, usize>,
}

/// Findings kept per (rule, dataset); the rest are only counted.
pub const MAX_PER_RULE: usize = 20;

impl Report {
    /// Whether nothing failed.
    pub fn passed(&self) -> bool {
        self.counts.is_empty()
    }

    pub(crate) fn add(
        &mut self,
        level: Level,
        rule: &'static str,
        dataset: &str,
        variable: Option<&str>,
        row: Option<usize>,
        message: String,
    ) {
        let total = self.counts.entry(rule).or_insert(0);
        *total += 1;
        let kept = self
            .findings
            .iter()
            .filter(|f| f.rule == rule && f.dataset == dataset)
            .count();
        if kept < MAX_PER_RULE {
            self.findings.push(Finding {
                level,
                rule,
                dataset: dataset.to_string(),
                variable: variable.map(str::to_string),
                row,
                message,
            });
        }
    }

    /// A short multi-line summary.
    pub fn summary(&self) -> String {
        if self.passed() {
            return "all L0–L2 checks passed".into();
        }
        let mut s = String::new();
        for (rule, n) in &self.counts {
            s.push_str(&format!("{rule}: {n} failure(s)\n"));
        }
        for f in self.findings.iter().take(10) {
            s.push_str(&format!("  {f}\n"));
        }
        s
    }
}

/// Runs L0–L2 over a study's SDTM and ADaM datasets.
pub fn run(sdtm: &[Dataset], adam: &[Dataset]) -> Report {
    let mut r = Report::default();
    structural::run(sdtm, adam, &mut r);
    referential::run(sdtm, adam, &mut r);
    derivational::run(sdtm, adam, &mut r);
    r
}
