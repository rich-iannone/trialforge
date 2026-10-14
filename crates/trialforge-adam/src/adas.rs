//! ADADAS (ADAS-Cog analysis) from FT, ADSL and TV, in admiral style:
//!
//! - Item records (`PARAMCD` = `FTTESTCD`) and a total, `ACTOT`, per visit: the sum of the 11
//!   items, prorated as sum × 11 / n when 8 to 10 are present (missing with fewer).
//! - Analysis windows by `ADY`: Baseline (≤ 1), then each planned assessment, bounded by the
//!   midpoints between planned days; in each window the record closest to the planned day (the
//!   earlier on ties) is the analysis record (`ANL01FL`).
//! - `ABLFL`: the last record with `ADY` ≤ 1; `BASE`, `CHG` on post-baseline records.
//! - LOCF for `ACTOT` (admiral `derive_locf_records()`): every planned visit without an
//!   analysis record takes the last earlier analysis value (`DTYPE = LOCF`).

use std::collections::BTreeMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::{DeriveError, Frame, dt, index, rel_day};

const ITEMS: usize = 11;

#[derive(Clone, Debug, Default)]
struct Rec {
    sl: usize,
    paramcd: String,
    param: String,
    paramn: f64,
    aval: Option<f64>,
    adt: Option<Date>,
    ady: Option<i64>,
    avisitn: Option<i64>,
    target: Option<i64>,
    ablfl: bool,
    anl01: bool,
    locf: bool,
    srcseq: Option<i64>,
}

/// Planned assessment days after baseline (from TV VISITDY of the visits FT was collected at).
fn planned(ft: &Frame<'_>, tv: Option<&Dataset>) -> Result<Vec<i64>, DeriveError> {
    let visits: std::collections::BTreeSet<i64> =
        ft.ints("VISITDY")?.into_iter().flatten().collect();
    let mut days: Vec<i64> = match tv.map(Frame::new) {
        Some(tv) => tv
            .ints("VISITDY")?
            .into_iter()
            .flatten()
            .filter(|d| visits.contains(d))
            .collect(),
        None => visits.into_iter().collect(),
    };
    days.retain(|d| *d > 1);
    days.sort_unstable();
    days.dedup();
    Ok(days)
}

/// Derives ADADAS.
pub fn adadas(ft: &Dataset, adsl: &Dataset, tv: Option<&Dataset>) -> Result<Dataset, DeriveError> {
    let ft = Frame::new(ft);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;
    let targets = planned(&ft, tv)?;
    let week = |day: i64| (day - 1) / 7;
    // Window of an analysis day: (AVISITN, target day).
    let window = |ady: i64| -> (i64, Option<i64>) {
        if ady <= 1 {
            return (0, None);
        }
        let mut lo = 1;
        for (k, &t) in targets.iter().enumerate() {
            let hi = targets.get(k + 1).map_or(i64::MAX, |n| (t + n) / 2);
            if ady > lo && ady <= hi {
                return (week(t), Some(t));
            }
            lo = hi;
        }
        (week(*targets.last().unwrap_or(&1)), targets.last().copied())
    };

    let (u, testcd, test, stresn, dtc, seq, visitnum) = (
        ft.chars("USUBJID")?,
        ft.chars("FTTESTCD")?,
        ft.chars("FTTEST")?,
        ft.nums("FTSTRESN")?,
        ft.chars("FTDTC")?,
        ft.ints("FTSEQ")?,
        ft.ints("VISITNUM")?,
    );
    let order: Vec<&str> = {
        let mut o: Vec<&str> = Vec::new();
        for t in testcd.iter().flatten() {
            if !o.contains(&t.as_str()) {
                o.push(t);
            }
        }
        o
    };
    let mut recs: Vec<Rec> = Vec::new();
    // Totals per (subject, visit, date): (sum, n).
    let mut totals: BTreeMap<(usize, i64, Date), (f64, usize)> = BTreeMap::new();
    for i in 0..ft.n() {
        let (Some(x), Some(v)) = (u[i].as_deref(), stresn[i]) else {
            continue;
        };
        let Some(&r) = sl_idx.get(x) else { continue };
        let adt = dt(&dtc[i]);
        let ady = rel_day(adt, trtsdt[r]);
        let (avisitn, target) = ady.map_or((None, None), |d| {
            let (w, t) = window(d);
            (Some(w), t)
        });
        let code = testcd[i].clone().unwrap_or_default();
        recs.push(Rec {
            sl: r,
            paramn: order
                .iter()
                .position(|o| *o == code)
                .map_or(0.0, |p| p as f64 + 1.0),
            paramcd: code,
            param: test[i].clone().unwrap_or_default(),
            aval: Some(v),
            adt,
            ady,
            avisitn,
            target,
            srcseq: seq[i],
            ..Rec::default()
        });
        if let (Some(vn), Some(d)) = (visitnum[i], adt) {
            let e = totals.entry((r, vn, d)).or_insert((0.0, 0));
            e.0 += v;
            e.1 += 1;
        }
    }
    for ((r, _, d), (sum, n)) in totals {
        let aval = if n == ITEMS {
            Some(sum)
        } else if n >= 8 {
            Some(sum * ITEMS as f64 / n as f64)
        } else {
            None
        };
        let Some(aval) = aval else { continue };
        let ady = rel_day(Some(d), trtsdt[r]);
        let (avisitn, target) = ady.map_or((None, None), |dd| {
            let (w, t) = window(dd);
            (Some(w), t)
        });
        recs.push(Rec {
            sl: r,
            paramcd: "ACTOT".into(),
            param: "ADAS-Cog(11) Subscore".into(),
            paramn: 12.0,
            aval: Some(aval),
            adt: Some(d),
            ady,
            avisitn,
            target,
            ..Rec::default()
        });
    }

    // Flags per subject and parameter.
    let mut groups: BTreeMap<(usize, String), Vec<usize>> = BTreeMap::new();
    for (k, x) in recs.iter().enumerate() {
        groups.entry((x.sl, x.paramcd.clone())).or_default().push(k);
    }
    let mut locf: Vec<Rec> = Vec::new();
    for ((sl_row, paramcd), idx) in &groups {
        let base = idx
            .iter()
            .copied()
            .filter(|&k| recs[k].ady.is_some_and(|d| d <= 1))
            .max_by_key(|&k| (recs[k].adt, recs[k].srcseq));
        if let Some(b) = base {
            recs[b].ablfl = true;
            recs[b].anl01 = true;
        }
        let mut best: BTreeMap<i64, usize> = BTreeMap::new();
        for &k in idx {
            let (Some(w), Some(t)) = (recs[k].avisitn, recs[k].target) else {
                continue;
            };
            let key = |k: usize| {
                (
                    (recs[k].ady.unwrap_or(0) - t).abs(),
                    recs[k].ady,
                    recs[k].srcseq,
                )
            };
            match best.get(&w) {
                Some(&b) if key(b) <= key(k) => {}
                _ => {
                    best.insert(w, k);
                }
            }
        }
        for &k in best.values() {
            recs[k].anl01 = true;
        }
        // LOCF for the total over the planned visits.
        if paramcd == "ACTOT" {
            let mut avail: BTreeMap<i64, f64> = BTreeMap::new();
            for &k in idx.iter().filter(|&&k| recs[k].anl01) {
                if let (Some(w), Some(v)) = (recs[k].avisitn, recs[k].aval) {
                    avail.insert(w, v);
                }
            }
            for &t in &targets {
                let w = week(t);
                if avail.contains_key(&w) {
                    continue;
                }
                if let Some((_, &v)) = avail.range(..w).next_back() {
                    locf.push(Rec {
                        sl: *sl_row,
                        paramcd: "ACTOT".into(),
                        param: "ADAS-Cog(11) Subscore".into(),
                        paramn: 12.0,
                        aval: Some(v),
                        avisitn: Some(w),
                        anl01: true,
                        locf: true,
                        ..Rec::default()
                    });
                    avail.insert(w, v);
                }
            }
        }
    }
    recs.extend(locf);
    let mut bases: BTreeMap<(usize, String), f64> = BTreeMap::new();
    for x in recs.iter().filter(|x| x.ablfl) {
        if let Some(v) = x.aval {
            bases.insert((x.sl, x.paramcd.clone()), v);
        }
    }
    recs.sort_by(|a, b| {
        (
            a.sl,
            a.paramn as i64,
            a.avisitn,
            a.adt.is_none(),
            a.adt,
            a.srcseq,
        )
            .cmp(&(
                b.sl,
                b.paramn as i64,
                b.avisitn,
                b.adt.is_none(),
                b.adt,
                b.srcseq,
            ))
    });

    let sl_rows: Vec<usize> = recs.iter().map(|x| x.sl).collect();
    let from_sl = |name: &str| -> Result<Values, DeriveError> {
        Ok(
            match &sl
                .ds
                .column(name)
                .ok_or_else(|| DeriveError::Input(format!("ADSL has no {name}")))?
                .values
            {
                Values::Char(v) => Values::Char(sl_rows.iter().map(|&r| v[r].clone()).collect()),
                Values::Num(v) => Values::Num(sl_rows.iter().map(|&r| v[r]).collect()),
                Values::Int(v) => Values::Int(sl_rows.iter().map(|&r| v[r]).collect()),
                Values::Date(v) => Values::Date(sl_rows.iter().map(|&r| v[r]).collect()),
                Values::DateTime(v) => Values::DateTime(sl_rows.iter().map(|&r| v[r]).collect()),
            },
        )
    };
    let base: Vec<Option<f64>> = recs
        .iter()
        .map(|x| bases.get(&(x.sl, x.paramcd.clone())).copied())
        .collect();
    let y = |b: bool| b.then(|| "Y".to_string());
    let label = |w: Option<i64>| {
        w.map(|w| {
            if w == 0 {
                "Baseline".to_string()
            } else {
                format!("Week {w}")
            }
        })
    };
    Ok(build(
        "ADADAS",
        vec![
            ("STUDYID", from_sl("STUDYID")?),
            ("USUBJID", from_sl("USUBJID")?),
            ("SUBJID", from_sl("SUBJID")?),
            ("SITEID", from_sl("SITEID")?),
            ("AGE", from_sl("AGE")?),
            ("SEX", from_sl("SEX")?),
            ("RACE", from_sl("RACE")?),
            ("SAFFL", from_sl("SAFFL")?),
            ("ITTFL", from_sl("ITTFL")?),
            ("TRTP", from_sl("TRT01P")?),
            ("TRTPN", from_sl("TRT01PN")?),
            ("TRTA", from_sl("TRT01A")?),
            ("TRTAN", from_sl("TRT01AN")?),
            (
                "PARAMCD",
                Values::Char(recs.iter().map(|x| Some(x.paramcd.clone())).collect()),
            ),
            (
                "PARAM",
                Values::Char(recs.iter().map(|x| Some(x.param.clone())).collect()),
            ),
            (
                "PARAMN",
                Values::Num(recs.iter().map(|x| Some(x.paramn)).collect()),
            ),
            ("AVAL", Values::Num(recs.iter().map(|x| x.aval).collect())),
            ("BASE", Values::Num(base.clone())),
            (
                "CHG",
                Values::Num(
                    recs.iter()
                        .zip(&base)
                        .map(|(x, b)| match (x.avisitn, x.aval, b) {
                            (Some(w), Some(a), Some(b)) if w > 0 => Some(a - b),
                            _ => None,
                        })
                        .collect(),
                ),
            ),
            ("ADT", Values::Date(recs.iter().map(|x| x.adt).collect())),
            ("ADY", Values::Int(recs.iter().map(|x| x.ady).collect())),
            (
                "AVISIT",
                Values::Char(recs.iter().map(|x| label(x.avisitn)).collect()),
            ),
            (
                "AVISITN",
                Values::Num(recs.iter().map(|x| x.avisitn.map(|w| w as f64)).collect()),
            ),
            (
                "ABLFL",
                Values::Char(recs.iter().map(|x| y(x.ablfl)).collect()),
            ),
            (
                "ANL01FL",
                Values::Char(recs.iter().map(|x| y(x.anl01)).collect()),
            ),
            (
                "DTYPE",
                Values::Char(
                    recs.iter()
                        .map(|x| x.locf.then(|| "LOCF".to_string()))
                        .collect(),
                ),
            ),
            (
                "SRCDOM",
                Values::Char(
                    recs.iter()
                        .map(|x| x.srcseq.map(|_| "FT".to_string()))
                        .collect(),
                ),
            ),
            (
                "SRCVAR",
                Values::Char(
                    recs.iter()
                        .map(|x| x.srcseq.map(|_| "FTSTRESN".to_string()))
                        .collect(),
                ),
            ),
            (
                "SRCSEQ",
                Values::Int(recs.iter().map(|x| x.srcseq).collect()),
            ),
        ],
    )?)
}
