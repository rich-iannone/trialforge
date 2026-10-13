//! Derivation layer: SDTM to ADaM (ADSL, ADAE, ADLB, ADVS), following admiral conventions.
//!
//! Derivations read SDTM *datasets* by variable name, never the simulation's internal state,
//! so what ADaM says is exactly what the SDTM supports. Inputs a real study also takes from
//! outside SDTM (randomisation strata, the data cutoff) are passed explicitly.

mod adae;
mod adsl;
mod bds;

use std::collections::HashMap;

use trialforge_cdisc::MetadataError;
use trialforge_core::{Dataset, Date};

pub use adae::adae;
pub use adsl::{AdslInputs, adsl};
pub use bds::{adlb, advs};

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

/// Parses a complete ISO 8601 date. Partial dates map to `None` (no imputation yet).
pub(crate) fn dt(v: &Option<String>) -> Option<Date> {
    v.as_deref().and_then(Date::parse_iso)
}

/// ADaM relative day: no day 0; the reference date is day 1.
pub(crate) fn rel_day(d: Option<Date>, reference: Option<Date>) -> Option<i64> {
    match (d, reference) {
        (Some(d), Some(r)) => Some(i64::from(trialforge_core::study_day(d, r))),
        _ => None,
    }
}

/// Row indices of `usubjid` values, for joining.
pub(crate) fn index(usubjid: &[Option<String>]) -> HashMap<&str, usize> {
    usubjid
        .iter()
        .enumerate()
        .filter_map(|(i, u)| u.as_deref().map(|u| (u, i)))
        .collect()
}

/// Derives ADSL, and ADAE, ADLB and ADVS where their SDTM sources exist.
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
    out.insert(0, adsl);
    Ok(out)
}
