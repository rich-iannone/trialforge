//! Derivation layer: SDTM to ADaM (ADSL, OCCDS, BDS, TTE, ADRS), following admiral
//! conventions documented in `spec/derivations/`.
//!
//! Derivations read SDTM *datasets* by variable name, never the simulation's internal state,
//! so what ADaM says is exactly what the SDTM supports. The only inputs from outside SDTM are
//! the ones a real study also takes from elsewhere: the stratification factor (from the
//! randomisation system, as in SIMONC01) and the data cutoff date.

use std::collections::HashMap;

use trialforge_cdisc::{ADSL, MetadataError};
use trialforge_core::{Column, Dataset, Date, Values};

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

fn column<'a>(ds: &'a Dataset, name: &str) -> Result<&'a Column, DeriveError> {
    ds.column(name)
        .ok_or_else(|| DeriveError::Input(format!("{} has no {name} variable", ds.name)))
}

fn chars<'a>(ds: &'a Dataset, name: &str) -> Result<&'a [Option<String>], DeriveError> {
    column(ds, name)?
        .as_char()
        .ok_or_else(|| DeriveError::Input(format!("{}.{name} must be character", ds.name)))
}

fn ints<'a>(ds: &'a Dataset, name: &str) -> Result<&'a [Option<i64>], DeriveError> {
    column(ds, name)?
        .as_int()
        .ok_or_else(|| DeriveError::Input(format!("{}.{name} must be numeric", ds.name)))
}

/// Parses a complete ISO 8601 date. Partial dates map to `None` (no imputation in the spike).
fn dt(v: &Option<String>) -> Option<Date> {
    v.as_deref().and_then(Date::parse_iso)
}

/// Per-subject facts gathered from EX and DS.
#[derive(Default)]
struct FromEvents {
    first_dose: Option<Date>,
    last_dose: Option<Date>,
    randomized: Option<Date>,
    treatment_reason: Option<String>,
    study_reason: Option<String>,
    study_end: Option<Date>,
}

/// Derives ADSL from DM, EX and DS.
///
/// `strata` maps `USUBJID` to the randomised stratum (`STRAT1R`); `cutoff` is `DCUTDT`.
pub fn adsl(
    dm: &Dataset,
    ex: &Dataset,
    ds: &Dataset,
    strata: &HashMap<String, String>,
    cutoff: Date,
) -> Result<Dataset, DeriveError> {
    let mut facts: HashMap<&str, FromEvents> = HashMap::new();

    // EX: first and last exposure (admiral: derive_vars_merged on EXSTDTC/EXENDTC).
    for ((u, st), en) in chars(ex, "USUBJID")?
        .iter()
        .zip(chars(ex, "EXSTDTC")?)
        .zip(chars(ex, "EXENDTC")?)
    {
        let Some(u) = u.as_deref() else { continue };
        let f = facts.entry(u).or_default();
        if let Some(d) = dt(st) {
            f.first_dose = Some(f.first_dose.map_or(d, |x| x.min(d)));
        }
        if let Some(d) = dt(en) {
            f.last_dose = Some(f.last_dose.map_or(d, |x| x.max(d)));
        }
    }

    // DS: randomisation date and the two disposition events.
    let ds_cols = (
        chars(ds, "USUBJID")?,
        chars(ds, "DSDECOD")?,
        chars(ds, "DSCAT")?,
        chars(ds, "EPOCH")?,
        chars(ds, "DSSTDTC")?,
    );
    for i in 0..ds.n_rows() {
        let Some(u) = ds_cols.0[i].as_deref() else {
            continue;
        };
        let f = facts.entry(u).or_default();
        let (decod, cat, epoch, date) = (
            &ds_cols.1[i],
            ds_cols.2[i].as_deref(),
            ds_cols.3[i].as_deref(),
            dt(&ds_cols.4[i]),
        );
        match (decod.as_deref(), cat, epoch) {
            (Some("RANDOMIZED"), _, _) => f.randomized = date,
            (_, Some("DISPOSITION EVENT"), Some("TREATMENT")) => f.treatment_reason = decod.clone(),
            (_, Some("DISPOSITION EVENT"), Some("FOLLOW-UP")) => {
                f.study_reason = decod.clone();
                f.study_end = date;
            }
            _ => {}
        }
    }

    let usubjid = chars(dm, "USUBJID")?;
    let none = FromEvents::default();
    let per: Vec<&FromEvents> = usubjid
        .iter()
        .map(|u| u.as_deref().and_then(|u| facts.get(u)).unwrap_or(&none))
        .collect();
    let copy = |name: &str| -> Result<Values, DeriveError> { Ok(column(dm, name)?.values.clone()) };
    let age = ints(dm, "AGE")?;
    let arm = chars(dm, "ARM")?;
    let actarm = chars(dm, "ACTARM")?;
    let dthdtc = chars(dm, "DTHDTC")?;
    let rfpendtc = chars(dm, "RFPENDTC")?;
    let n = dm.n_rows();

    let status = |reason: &Option<String>| {
        Some(
            if reason.is_some() {
                "DISCONTINUED"
            } else {
                "ONGOING"
            }
            .to_string(),
        )
    };
    let flag = |b: bool| Some(if b { "Y" } else { "N" }.to_string());

    Ok(ADSL.build(vec![
        ("STUDYID", copy("STUDYID")?),
        ("USUBJID", copy("USUBJID")?),
        ("SUBJID", copy("SUBJID")?),
        ("SITEID", copy("SITEID")?),
        ("COUNTRY", copy("COUNTRY")?),
        ("AGE", copy("AGE")?),
        ("AGEU", copy("AGEU")?),
        (
            "AGEGR1",
            Values::Char(
                age.iter()
                    .map(|a| a.map(|a| if a < 65 { "<65" } else { ">=65" }.to_string()))
                    .collect(),
            ),
        ),
        ("SEX", copy("SEX")?),
        ("RACE", copy("RACE")?),
        ("ETHNIC", copy("ETHNIC")?),
        ("ARMCD", copy("ARMCD")?),
        ("ARM", copy("ARM")?),
        ("ACTARMCD", copy("ACTARMCD")?),
        ("ACTARM", copy("ACTARM")?),
        ("TRT01P", Values::Char(arm.to_vec())),
        ("TRT01A", Values::Char(actarm.to_vec())),
        (
            "STRAT1R",
            Values::Char(
                usubjid
                    .iter()
                    .map(|u| u.as_ref().and_then(|u| strata.get(u).cloned()))
                    .collect(),
            ),
        ),
        (
            "ITTFL",
            Values::Char(per.iter().map(|f| flag(f.randomized.is_some())).collect()),
        ),
        (
            "SAFFL",
            Values::Char(per.iter().map(|f| flag(f.first_dose.is_some())).collect()),
        ),
        (
            "RANDDT",
            Values::Date(per.iter().map(|f| f.randomized).collect()),
        ),
        (
            "TRTSDT",
            Values::Date(per.iter().map(|f| f.first_dose).collect()),
        ),
        (
            "TRTEDT",
            Values::Date(per.iter().map(|f| f.last_dose).collect()),
        ),
        (
            "EOTSTT",
            Values::Char(per.iter().map(|f| status(&f.treatment_reason)).collect()),
        ),
        (
            "DCTREAS",
            Values::Char(per.iter().map(|f| f.treatment_reason.clone()).collect()),
        ),
        (
            "EOSSTT",
            Values::Char(per.iter().map(|f| status(&f.study_reason)).collect()),
        ),
        (
            "DCSREAS",
            Values::Char(per.iter().map(|f| f.study_reason.clone()).collect()),
        ),
        (
            "EOSDT",
            Values::Date(per.iter().map(|f| f.study_end).collect()),
        ),
        ("DTHFL", copy("DTHFL")?),
        ("DTHDT", Values::Date(dthdtc.iter().map(dt).collect())),
        // Last known alive: the later of the last dose and the end of participation.
        (
            "LSTALVDT",
            Values::Date(
                per.iter()
                    .zip(rfpendtc)
                    .map(|(f, p)| f.last_dose.max(dt(p)))
                    .collect(),
            ),
        ),
        ("DCUTDT", Values::Date(vec![Some(cutoff); n])),
    ])?)
}
