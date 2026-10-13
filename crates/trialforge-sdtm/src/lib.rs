//! Observation layer: turns subject [`Timeline`]s into SDTM datasets.
//!
//! - Trial design: TA, TE, TV, TI, TS ([`tdm`]).
//! - Subject level: DM, SE, SV, EX, CM, AE (+SUPPAE), DS, MH, LB, VS ([`subject`]).
//!
//! Conventions shared by every domain:
//! - Records are sorted by `USUBJID` and the domain's natural key, then `--SEQ` is assigned.
//! - `--DY` is relative to `RFSTDTC` (first dose), with no day 0.
//! - `EPOCH` comes from the subject's elements ([`Timeline::epoch_of`]).
//! - Datasets with no records are not emitted (DM always is).

pub mod subject;
pub mod tdm;

use trialforge_cdisc::MetadataError;
use trialforge_core::{Dataset, Date, Values, study_day};
use trialforge_model::{Design, Subject};
use trialforge_processes::timeline::Timeline;

/// Everything the emitters read.
pub struct Context<'a> {
    /// The design.
    pub design: &'a Design,
    /// Subjects, in enrollment order.
    pub subjects: &'a [Subject],
    /// One timeline per subject.
    pub timelines: &'a [Timeline],
    /// Subject indices sorted by `USUBJID`.
    pub order: Vec<usize>,
}

impl<'a> Context<'a> {
    /// A context over parallel `subjects` and `timelines`.
    pub fn new(design: &'a Design, subjects: &'a [Subject], timelines: &'a [Timeline]) -> Self {
        assert_eq!(subjects.len(), timelines.len(), "one timeline per subject");
        let mut order: Vec<usize> = (0..subjects.len()).collect();
        order.sort_by(|&a, &b| subjects[a].usubjid.cmp(&subjects[b].usubjid));
        Self {
            design,
            subjects,
            timelines,
            order,
        }
    }

    /// (subject, timeline) pairs in `USUBJID` order.
    pub fn sorted(&self) -> impl Iterator<Item = (&'a Subject, &'a Timeline)> + '_ {
        self.order
            .iter()
            .map(|&i| (&self.subjects[i], &self.timelines[i]))
    }
}

/// Emits every SDTM dataset that has records, in standard order.
pub fn emit(ctx: &Context<'_>) -> Result<Vec<Dataset>, MetadataError> {
    let mut out = tdm::emit(ctx)?;
    out.extend(subject::emit(ctx)?);
    Ok(out
        .into_iter()
        .filter(|d| d.n_rows() > 0 || d.name == "DM")
        .collect())
}

// Column helpers: build a column from rows with an accessor.

pub(crate) fn chars<T>(rows: &[T], f: impl Fn(&T) -> Option<String>) -> Values {
    Values::Char(rows.iter().map(f).collect())
}

pub(crate) fn ints<T>(rows: &[T], f: impl Fn(&T) -> Option<i64>) -> Values {
    Values::Int(rows.iter().map(f).collect())
}

pub(crate) fn nums<T>(rows: &[T], f: impl Fn(&T) -> Option<f64>) -> Values {
    Values::Num(rows.iter().map(f).collect())
}

pub(crate) fn constant<T>(rows: &[T], v: &str) -> Values {
    Values::Char(vec![Some(v.to_string()); rows.len()])
}

pub(crate) fn dy(d: Date, reference: Date) -> Option<i64> {
    Some(i64::from(study_day(d, reference)))
}

pub(crate) fn yn(b: bool) -> Option<String> {
    Some(if b { "Y" } else { "N" }.to_string())
}

/// `x` rounded to `dp` decimals and formatted, without a negative zero.
pub(crate) fn fmt_dp(x: f64, dp: u32) -> String {
    let r = trialforge_core::math::round_dp(x, dp);
    let s = format!("{:.*}", dp as usize, r);
    if s.starts_with('-') && s[1..].chars().all(|c| c == '0' || c == '.') {
        s[1..].to_string()
    } else {
        s
    }
}
