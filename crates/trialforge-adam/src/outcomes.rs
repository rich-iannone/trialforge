//! ADTTE for an outcomes study, from CE and ADSL, with admiral `derive_param_tte()` semantics
//! (the earliest event; otherwise the latest censoring date, `LSTALVDT` or `RANDDT`):
//!
//! - `MACE`: time to the first of cardiac death, myocardial infarction or ischaemic stroke;
//! - `CVDEATH`, `MI`, `STROKE`: each component (non-cardiovascular death censors at death,
//!   through `LSTALVDT`);
//! - `DEATH`: all-cause death (`ADSL.DTHDT`).

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::onco::{ADSL_KEYS, Source, adsl_columns, tte};
use crate::{DeriveError, Frame, dt, index};

/// (PARAMCD, PARAM, event PTs, event description).
const PARAMS: [(&str, &str, &[&str]); 5] = [
    (
        "CVDEATH",
        "Time to Cardiovascular Death",
        &["Cardiac death"],
    ),
    ("DEATH", "Time to Death from Any Cause", &[]),
    (
        "MACE",
        "Time to First MACE",
        &["Cardiac death", "Myocardial infarction", "Ischaemic stroke"],
    ),
    (
        "MI",
        "Time to First Myocardial Infarction",
        &["Myocardial infarction"],
    ),
    (
        "STROKE",
        "Time to First Ischaemic Stroke",
        &["Ischaemic stroke"],
    ),
];

/// `CECAT` of the adjudicated endpoint events.
pub(crate) const CATEGORY: &str = "ADJUDICATED ENDPOINT EVENT";

/// The event description of a Preferred Term.
fn desc(pt: &str) -> &'static str {
    match pt {
        "Cardiac death" => "Cardiovascular death",
        "Myocardial infarction" => "Myocardial infarction",
        "Ischaemic stroke" => "Ischaemic stroke",
        _ => "Death",
    }
}

/// Derives ADTTE.
pub fn adtte_outcomes(ce: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let ce = Frame::new(ce);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let (randdt, lstalvdt, dthdt) = (
        sl.dates("RANDDT")?,
        sl.dates("LSTALVDT")?,
        sl.dates("DTHDT")?,
    );
    // First event of each PT per subject: (date, CESEQ).
    let mut first: HashMap<(usize, &str), (Date, Option<i64>)> = HashMap::new();
    let (u, decod, dtc, seq, cat) = (
        ce.chars("USUBJID")?,
        ce.chars("CEDECOD")?,
        ce.chars("CESTDTC")?,
        ce.ints("CESEQ")?,
        ce.chars("CECAT")?,
    );
    for i in 0..ce.n() {
        if cat[i].as_deref() != Some(CATEGORY) {
            continue;
        }
        let (Some(x), Some(pt), Some(d)) = (u[i].as_deref(), decod[i].as_deref(), dt(&dtc[i]))
        else {
            continue;
        };
        let Some(&r) = sl_idx.get(x) else { continue };
        let e = first.entry((r, pt)).or_insert((d, seq[i]));
        if (d, seq[i]) < *e {
            *e = (d, seq[i]);
        }
    }
    let mut recs = Vec::new();
    for r in 0..sl.n() {
        let Some(start) = randdt[r] else { continue };
        let censor = |d: Option<Date>, desc: &'static str, var: &'static str| {
            d.map(|date| Source {
                date,
                event: false,
                desc,
                cnsdesc: Some(if var == "LSTALVDT" {
                    "Alive During Study"
                } else {
                    "Randomization"
                }),
                srcdom: "ADSL",
                srcvar: var,
                srcseq: None,
            })
        };
        let censors = [
            censor(lstalvdt[r], "Alive", "LSTALVDT"),
            censor(randdt[r], "Randomization", "RANDDT"),
        ];
        for (paramcd, param, pts) in PARAMS {
            let events: Vec<Option<Source>> = if pts.is_empty() {
                vec![dthdt[r].map(|date| Source {
                    date,
                    event: true,
                    desc: "Death",
                    cnsdesc: None,
                    srcdom: "ADSL",
                    srcvar: "DTHDT",
                    srcseq: None,
                })]
            } else {
                pts.iter()
                    .map(|pt| {
                        first.get(&(r, *pt)).map(|&(date, seq)| Source {
                            date,
                            event: true,
                            desc: desc(pt),
                            cnsdesc: None,
                            srcdom: "CE",
                            srcvar: "CESTDTC",
                            srcseq: seq,
                        })
                    })
                    .collect()
            };
            if let Some(s) = tte(&events, &censors) {
                recs.push(TteRec {
                    sl: r,
                    paramcd,
                    param,
                    start,
                    s,
                });
            }
        }
    }
    tte_dataset(&sl, &recs)
}

/// One ADTTE record: the ADSL row, parameter, start date and the chosen event or censoring.
pub(crate) struct TteRec {
    pub sl: usize,
    pub paramcd: &'static str,
    pub param: &'static str,
    pub start: Date,
    pub s: Source,
}

/// Builds ADTTE from records in ADSL-row then parameter order.
pub(crate) fn tte_dataset(sl: &Frame<'_>, recs: &[TteRec]) -> Result<Dataset, DeriveError> {
    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let mut aseq: Vec<Option<i64>> = Vec::with_capacity(recs.len());
    for (i, r) in recs.iter().enumerate() {
        let k = if i > 0 && recs[i - 1].sl == r.sl {
            aseq[i - 1].unwrap_or(0) + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }
    let adt: Vec<Date> = recs.iter().map(|r| r.s.date.max(r.start)).collect();
    let mut cols = adsl_columns(sl, &sl_rows, &ADSL_KEYS)?;
    cols.extend(adsl_columns(sl, &sl_rows, &["AGE", "SEX", "RANDDT"])?);
    cols.extend(vec![
        ("ASEQ", Values::Int(aseq)),
        (
            "PARAMCD",
            Values::Char(recs.iter().map(|r| Some(r.paramcd.to_string())).collect()),
        ),
        (
            "PARAM",
            Values::Char(recs.iter().map(|r| Some(r.param.to_string())).collect()),
        ),
        (
            "AVAL",
            Values::Num(
                recs.iter()
                    .zip(&adt)
                    .map(|(r, d)| Some(f64::from(d.days_since(r.start) + 1)))
                    .collect(),
            ),
        ),
        (
            "AVALU",
            Values::Char(vec![Some("DAYS".to_string()); recs.len()]),
        ),
        (
            "STARTDT",
            Values::Date(recs.iter().map(|r| Some(r.start)).collect()),
        ),
        ("ADT", Values::Date(adt.iter().map(|d| Some(*d)).collect())),
        (
            "CNSR",
            Values::Int(recs.iter().map(|r| Some(i64::from(!r.s.event))).collect()),
        ),
        (
            "EVNTDESC",
            Values::Char(recs.iter().map(|r| Some(r.s.desc.to_string())).collect()),
        ),
        (
            "CNSDTDSC",
            Values::Char(
                recs.iter()
                    .map(|r| r.s.cnsdesc.map(str::to_string))
                    .collect(),
            ),
        ),
        (
            "SRCDOM",
            Values::Char(recs.iter().map(|r| Some(r.s.srcdom.to_string())).collect()),
        ),
        (
            "SRCVAR",
            Values::Char(recs.iter().map(|r| Some(r.s.srcvar.to_string())).collect()),
        ),
        (
            "SRCSEQ",
            Values::Int(recs.iter().map(|r| r.s.srcseq).collect()),
        ),
    ]);
    Ok(build("ADTTE", cols)?)
}
