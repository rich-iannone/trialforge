//! Derivation layer: SDTM to ADaM (ADSL, ADAE, ADLB, ADVS, ADRS, ADTTE, ADTR, ADPC, ADPP),
//! following
//! admiral and admiralonco conventions.
//!
//! Derivations read SDTM *datasets* by variable name, never the simulation's internal state,
//! so what ADaM says is exactly what the SDTM supports. Inputs a real study also takes from
//! outside SDTM (randomisation strata, the data cutoff) are passed explicitly.

mod adae;
mod adas;
mod adsl;
mod bds;
mod onco;
mod outcomes;
mod pk;

use std::collections::HashMap;

use trialforge_cdisc::MetadataError;
use trialforge_core::{Dataset, Date, Values};

pub use adae::adae;
pub use adas::adadas;
pub use adsl::{AdslInputs, adsl};
pub use bds::{adlb, advs};
pub use onco::{adrs, adtr, adtte};
pub use outcomes::adtte_outcomes;
pub use pk::{adpc, adpp};

/// Errors from deriving ADaM.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum DeriveError {
    /// An input dataset lacks a variable or has it with the wrong type.
    Input(String),
    /// The output does not fit its metadata.
    Metadata(MetadataError),
}

impl core::fmt::Display for DeriveError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            DeriveError::Input(m) => f.write_str(m),
            DeriveError::Metadata(e) => e.fmt(f),
        }
    }
}

impl std::error::Error for DeriveError {}

impl From<MetadataError> for DeriveError {
    fn from(e: MetadataError) -> Self {
        DeriveError::Metadata(e)
    }
}

/// Typed access to a dataset's columns by name.
pub(crate) struct Frame<'a> {
    pub ds: &'a Dataset,
}

impl<'a> Frame<'a> {
    pub fn new(ds: &'a Dataset) -> Self {
        Self { ds }
    }

    fn column(&self, name: &str) -> Result<&'a trialforge_core::Column, DeriveError> {
        self.ds
            .column(name)
            .ok_or_else(|| DeriveError::Input(format!("{} has no {name} variable", self.ds.name)))
    }

    pub fn has(&self, name: &str) -> bool {
        self.ds.column(name).is_some()
    }

    pub fn chars(&self, name: &str) -> Result<&'a [Option<String>], DeriveError> {
        self.column(name)?
            .as_char()
            .ok_or_else(|| DeriveError::Input(format!("{}.{name} must be character", self.ds.name)))
    }

    /// Numeric values as `f64` (integer and float columns alike).
    pub fn nums(&self, name: &str) -> Result<Vec<Option<f64>>, DeriveError> {
        match &self.column(name)?.values {
            trialforge_core::Values::Num(v) => Ok(v.clone()),
            trialforge_core::Values::Int(v) => Ok(v.iter().map(|x| x.map(|n| n as f64)).collect()),
            _ => Err(DeriveError::Input(format!(
                "{}.{name} must be numeric",
                self.ds.name
            ))),
        }
    }

    pub fn ints(&self, name: &str) -> Result<Vec<Option<i64>>, DeriveError> {
        match &self.column(name)?.values {
            trialforge_core::Values::Int(v) => Ok(v.clone()),
            _ => Err(DeriveError::Input(format!(
                "{}.{name} must be integer",
                self.ds.name
            ))),
        }
    }

    pub fn dates(&self, name: &str) -> Result<Vec<Option<Date>>, DeriveError> {
        match &self.column(name)?.values {
            trialforge_core::Values::Date(v) => Ok(v.clone()),
            _ => Err(DeriveError::Input(format!(
                "{}.{name} must be a date",
                self.ds.name
            ))),
        }
    }

    pub fn n(&self) -> usize {
        self.ds.n_rows()
    }
}

/// The dataset called `name`, if present.
pub(crate) fn find<'a>(datasets: &'a [Dataset], name: &str) -> Option<&'a Dataset> {
    datasets.iter().find(|d| d.name == name)
}

/// The dataset called `name`, or an input error.
pub(crate) fn require<'a>(datasets: &'a [Dataset], name: &str) -> Result<&'a Dataset, DeriveError> {
    find(datasets, name).ok_or_else(|| DeriveError::Input(format!("SDTM has no {name} dataset")))
}

/// The date of a complete ISO 8601 date or date-time. Partial dates map to `None` (no
/// imputation yet).
pub(crate) fn dt(v: &Option<String>) -> Option<Date> {
    v.as_deref().and_then(|s| {
        Date::parse_iso(s).or_else(|| trialforge_core::DateTime::parse_iso(s).map(|d| d.date()))
    })
}

/// A complete ISO 8601 date-time.
pub(crate) fn dtm(v: &Option<String>) -> Option<trialforge_core::DateTime> {
    v.as_deref().and_then(trialforge_core::DateTime::parse_iso)
}

/// ADaM relative day: no day 0; the reference date is day 1.
pub(crate) fn rel_day(d: Option<Date>, reference: Option<Date>) -> Option<i64> {
    match (d, reference) {
        (Some(d), Some(r)) => Some(i64::from(trialforge_core::study_day(d, r))),
        _ => None,
    }
}

/// Crossover periods from ADSL (`APxxSDT`, `APxxEDT`, `TRTxxP/PN/A/AN`); empty otherwise.
pub(crate) struct Periods {
    /// Per ADSL row: (period, start, end, TRTxxP, TRTxxPN, TRTxxA, TRTxxAN).
    rows: Vec<Vec<PeriodRow>>,
}

/// One subject's period: (period, start, end, planned, planned (N), actual, actual (N)).
pub(crate) type PeriodRow = (
    i64,
    Date,
    Date,
    Option<String>,
    Option<f64>,
    Option<String>,
    Option<f64>,
);

impl Periods {
    pub fn new(adsl: &Dataset) -> Result<Self, DeriveError> {
        let sl = Frame::new(adsl);
        let mut rows = vec![Vec::new(); sl.n()];
        for k in 1..=4 {
            let (sv, ev) = (format!("AP{k:02}SDT"), format!("AP{k:02}EDT"));
            if !sl.has(&sv) {
                break;
            }
            let (st, en) = (sl.dates(&sv)?, sl.dates(&ev)?);
            let p = sl.chars(&format!("TRT{k:02}P"))?;
            let a = sl.chars(&format!("TRT{k:02}A"))?;
            let pn = sl.nums(&format!("TRT{k:02}PN"))?;
            let an = sl.nums(&format!("TRT{k:02}AN"))?;
            for (r, out) in rows.iter_mut().enumerate() {
                if let (Some(s), Some(e)) = (st[r], en[r]) {
                    out.push((k, s, e, p[r].clone(), pn[r], a[r].clone(), an[r]));
                }
            }
        }
        Ok(Self { rows })
    }

    /// Whether the study has periods.
    pub fn any(&self) -> bool {
        self.rows.iter().any(|r| !r.is_empty())
    }

    /// The period of ADSL row `r` containing `date`.
    pub fn of(&self, r: usize, date: Option<Date>) -> Option<&PeriodRow> {
        let d = date?;
        self.rows[r].iter().find(|p| p.1 <= d && d <= p.2)
    }

    /// The treatment attribution of `date`: the latest period started on or before it, or the
    /// first period for dates before any (washout and follow-up records keep the treatment of
    /// the period before them; screening records the first period's).
    pub fn attributed(&self, r: usize, date: Option<Date>) -> Option<&PeriodRow> {
        let p = &self.rows[r];
        date.and_then(|d| p.iter().rev().find(|x| x.1 <= d))
            .or(p.first())
    }
}

/// For a crossover, replaces `TRTP`/`TRTPN`/`TRTA`/`TRTAN` (those present) with the values of
/// the period each record is attributed to (see [`Periods::attributed`]) and adds `APERIOD`
/// (only for records within a period).
pub(crate) fn by_period(
    cols: &mut Vec<(&str, Values)>,
    periods: &Periods,
    rows: &[(usize, Option<Date>)],
) {
    if !periods.any() {
        return;
    }
    let p: Vec<Option<&PeriodRow>> = rows
        .iter()
        .map(|(r, d)| periods.attributed(*r, *d))
        .collect();
    let within: Vec<Option<i64>> = rows
        .iter()
        .map(|(r, d)| periods.of(*r, *d).map(|x| x.0))
        .collect();
    for (name, values) in cols.iter_mut() {
        *values = match *name {
            "TRTP" => Values::Char(p.iter().map(|x| x.and_then(|x| x.3.clone())).collect()),
            "TRTPN" => Values::Num(p.iter().map(|x| x.and_then(|x| x.4)).collect()),
            "TRTA" => Values::Char(p.iter().map(|x| x.and_then(|x| x.5.clone())).collect()),
            "TRTAN" => Values::Num(p.iter().map(|x| x.and_then(|x| x.6)).collect()),
            _ => continue,
        };
    }
    cols.push(("APERIOD", Values::Int(within)));
}

/// Row indices of `usubjid` values, for joining.
pub(crate) fn index(usubjid: &[Option<String>]) -> HashMap<&str, usize> {
    usubjid
        .iter()
        .enumerate()
        .filter_map(|(i, u)| u.as_deref().map(|u| (u, i)))
        .collect()
}

/// Derives ADSL, and ADAE, ADLB, ADVS, ADRS, ADTTE and ADTR where their SDTM sources exist.
pub fn derive(sdtm: &[Dataset], inputs: &AdslInputs) -> Result<Vec<Dataset>, DeriveError> {
    let adsl = adsl(sdtm, inputs)?;
    let tv = find(sdtm, "TV");
    let mut out = Vec::new();
    if let Some(ae) = find(sdtm, "AE") {
        out.push(adae(ae, &adsl)?);
    }
    if let Some(lb) = find(sdtm, "LB") {
        out.push(adlb(lb, &adsl, tv)?);
    }
    if let Some(vs) = find(sdtm, "VS") {
        out.push(advs(vs, &adsl, tv)?);
    }
    if let Some(rs) = find(sdtm, "RS") {
        let adrs = adrs(rs, find(sdtm, "TU"), &adsl)?;
        out.push(adtte(&adrs, &adsl)?);
        out.insert(out.len() - 1, adrs);
    }
    if let (Some(pc), Some(ex)) = (find(sdtm, "PC"), find(sdtm, "EX")) {
        out.push(adpc(pc, ex, &adsl)?);
    }
    if let Some(ce) = find(sdtm, "CE") {
        out.push(adtte_outcomes(ce, &adsl)?);
    }
    if let Some(ft) = find(sdtm, "FT") {
        out.push(adadas(ft, &adsl, tv)?);
    }
    if let Some(pp) = find(sdtm, "PP") {
        out.push(adpp(pp, &adsl)?);
    }
    if let (Some(tr), Some(tu), Some(rs)) = (find(sdtm, "TR"), find(sdtm, "TU"), find(sdtm, "RS")) {
        out.push(adtr(tr, tu, rs, &adsl)?);
    }
    out.insert(0, adsl);
    Ok(out)
}
