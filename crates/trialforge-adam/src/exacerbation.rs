//! Exacerbation analysis datasets from CE (`CECAT = EXACERBATION`) and ADSL:
//!
//! - ADEXAC, one record per subject per parameter, for rate analyses (negative binomial with a
//!   log follow-up offset): `EXACN` moderate or severe exacerbations and `SEVEXACN` severe ones
//!   starting from `TRTSDT` to the end of observation (`LSTALVDT`), `FUPYR` the observation
//!   period in years (`(LSTALVDT − TRTSDT + 1) / 365.25`) and `EXACRATE` the annualized rate
//!   (`EXACN / FUPYR`). Counts follow admiral `derive_var_merged_summary()` semantics.
//! - ADTTE `TTFEXAC`: time from randomisation to the first exacerbation, censored at
//!   `LSTALVDT` (admiral `derive_param_tte()`, as for the outcomes endpoints).

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::onco::{ADSL_KEYS, Source, adsl_columns, tte};
use crate::outcomes::{TteRec, tte_dataset};
use crate::{DeriveError, Frame, dt, index};

/// `CECAT` of exacerbations.
pub(crate) const CATEGORY: &str = "EXACERBATION";

/// (PARAMCD, PARAM).
const PARAMS: [(&str, &str); 4] = [
    ("EXACN", "Number of Moderate or Severe Exacerbations"),
    ("SEVEXACN", "Number of Severe Exacerbations"),
    ("FUPYR", "Observation Period (years)"),
    ("EXACRATE", "Annualized Rate of Exacerbations"),
];

/// An exacerbation onset: (start, severe, CESEQ).
type Onset = (Date, bool, Option<i64>);

/// Exacerbation onsets per ADSL row, in date order.
fn onsets(
    ce: &Frame<'_>,
    sl_idx: &HashMap<&str, usize>,
) -> Result<HashMap<usize, Vec<Onset>>, DeriveError> {
    let (u, cat, sev, dtc, seq) = (
        ce.chars("USUBJID")?,
        ce.chars("CECAT")?,
        ce.chars("CESEV")?,
        ce.chars("CESTDTC")?,
        ce.ints("CESEQ")?,
    );
    let mut out: HashMap<usize, Vec<(Date, bool, Option<i64>)>> = HashMap::new();
    for i in 0..ce.n() {
        if cat[i].as_deref() != Some(CATEGORY) {
            continue;
        }
        let (Some(r), Some(d)) = (u[i].as_deref().and_then(|x| sl_idx.get(x)), dt(&dtc[i])) else {
            continue;
        };
        out.entry(*r)
            .or_default()
            .push((d, sev[i].as_deref() == Some("SEVERE"), seq[i]));
    }
    for v in out.values_mut() {
        v.sort();
    }
    Ok(out)
}

/// Derives ADEXAC.
pub fn adexac(ce: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let ce = Frame::new(ce);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let (trtsdt, lstalvdt) = (sl.dates("TRTSDT")?, sl.dates("LSTALVDT")?);
    let events = onsets(&ce, &sl_idx)?;
    let mut rows: Vec<(usize, usize, f64)> = Vec::new();
    for r in 0..sl.n() {
        let (Some(start), Some(end)) = (trtsdt[r], lstalvdt[r]) else {
            continue;
        };
        let within: Vec<&(Date, bool, Option<i64>)> = events
            .get(&r)
            .into_iter()
            .flatten()
            .filter(|e| e.0 >= start && e.0 <= end)
            .collect();
        let n = within.len() as f64;
        let severe = within.iter().filter(|e| e.1).count() as f64;
        let years = f64::from(end.days_since(start) + 1) / 365.25;
        for (k, v) in [n, severe, years, n / years].into_iter().enumerate() {
            rows.push((r, k, v));
        }
    }
    let sl_rows: Vec<usize> = rows.iter().map(|x| x.0).collect();
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_KEYS)?;
    cols.extend(adsl_columns(
        &sl,
        &sl_rows,
        &["AGE", "SEX", "TRTSDT", "LSTALVDT"],
    )?);
    cols.extend([
        (
            "TRTP",
            adsl_columns(&sl, &sl_rows, &["TRT01P"])?.remove(0).1,
        ),
        (
            "TRTA",
            adsl_columns(&sl, &sl_rows, &["TRT01A"])?.remove(0).1,
        ),
        (
            "PARAMCD",
            Values::Char(
                rows.iter()
                    .map(|x| Some(PARAMS[x.1].0.to_string()))
                    .collect(),
            ),
        ),
        (
            "PARAM",
            Values::Char(
                rows.iter()
                    .map(|x| Some(PARAMS[x.1].1.to_string()))
                    .collect(),
            ),
        ),
        (
            "PARAMN",
            Values::Num(rows.iter().map(|x| Some(x.1 as f64 + 1.0)).collect()),
        ),
        (
            "AVAL",
            Values::Num(rows.iter().map(|x| Some(x.2)).collect()),
        ),
    ]);
    Ok(build("ADEXAC", cols)?)
}

/// Derives ADTTE for an exacerbation study (`TTFEXAC`).
pub fn adtte_exacerbation(ce: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let ce = Frame::new(ce);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let (randdt, lstalvdt) = (sl.dates("RANDDT")?, sl.dates("LSTALVDT")?);
    let events = onsets(&ce, &sl_idx)?;
    let mut recs = Vec::new();
    for r in 0..sl.n() {
        let Some(start) = randdt[r] else { continue };
        let first = events
            .get(&r)
            .and_then(|v| v.first())
            .map(|&(date, _, seq)| Source {
                date,
                event: true,
                desc: "Exacerbation",
                cnsdesc: None,
                srcdom: "CE",
                srcvar: "CESTDTC",
                srcseq: seq,
            });
        let censor = |d: Option<Date>, desc: &'static str, cns: &'static str, var: &'static str| {
            d.map(|date| Source {
                date,
                event: false,
                desc,
                cnsdesc: Some(cns),
                srcdom: "ADSL",
                srcvar: var,
                srcseq: None,
            })
        };
        let censors = [
            censor(lstalvdt[r], "Alive", "Alive During Study", "LSTALVDT"),
            censor(randdt[r], "Randomization", "Randomization", "RANDDT"),
        ];
        if let Some(s) = tte(&[first], &censors) {
            recs.push(TteRec {
                sl: r,
                paramcd: "TTFEXAC",
                param: "Time to First Moderate or Severe Exacerbation",
                start,
                s,
            });
        }
    }
    tte_dataset(&sl, &recs)
}
