//! BDS datasets: ADLB from LB and ADVS from VS (with ADSL and TV).
//!
//! Following admiral's `ad_adlb.R` / `ad_advs.R` templates:
//! - `ABLFL`: the last non-missing record on or before `TRTSDT`, per parameter.
//! - `BASE`, `BNRIND` from that record; `CHG`, `PCHG` on post-baseline records only.
//! - `AVISIT`: the planned visit; early-termination records are windowed to the nearest
//!   planned on-treatment visit (by planned day from TV, with midpoint window bounds).
//! - `ANL01FL`: one record per subject, parameter and analysis visit: baseline, then for each
//!   post-baseline window the record closest to the target day (ties: the later one).

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values, math};

use crate::{DeriveError, Frame, dt, index, rel_day};

/// Derives ADLB.
pub fn adlb(lb: &Dataset, adsl: &Dataset, tv: Option<&Dataset>) -> Result<Dataset, DeriveError> {
    bds("ADLB", "LB", lb, adsl, tv)
}

/// Derives ADVS.
pub fn advs(vs: &Dataset, adsl: &Dataset, tv: Option<&Dataset>) -> Result<Dataset, DeriveError> {
    bds("ADVS", "VS", vs, adsl, tv)
}

/// Planned post-baseline visits from TV: (VISITNUM, VISIT, VISITDY).
fn planned(tv: Option<&Dataset>) -> Result<Vec<(i64, String, i64)>, DeriveError> {
    let Some(tv) = tv.map(Frame::new) else {
        return Ok(Vec::new());
    };
    let (num, name, day) = (
        tv.ints("VISITNUM")?,
        tv.chars("VISIT")?,
        tv.ints("VISITDY")?,
    );
    let mut out: Vec<(i64, String, i64)> = (0..tv.n())
        .filter_map(|i| match (num[i], &name[i], day[i]) {
            (Some(n), Some(v), Some(d)) if d > 1 => Some((n, v.clone(), d)),
            _ => None,
        })
        .collect();
    out.sort_by_key(|v| v.2);
    Ok(out)
}

/// Decimals in a character result (for rounding derived differences).
fn decimals(s: &Option<String>) -> u32 {
    s.as_deref()
        .and_then(|s| s.split_once('.'))
        .map_or(0, |(_, f)| f.len() as u32)
}

fn bds(
    name: &str,
    src: &str,
    ds: &Dataset,
    adsl: &Dataset,
    tv: Option<&Dataset>,
) -> Result<Dataset, DeriveError> {
    let f = Frame::new(ds);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let col = |v: &str| format!("{src}{v}");

    let usubjid = f.chars("USUBJID")?;
    let testcd = f.chars(&col("TESTCD"))?;
    let test = f.chars(&col("TEST"))?;
    let stresn = f.nums(&col("STRESN"))?;
    let stresc = f.chars(&col("STRESC"))?;
    let stresu = f.chars(&col("STRESU"))?;
    let seq = f.ints(&col("SEQ"))?;
    let dtc = f.chars(&col("DTC"))?;
    let visitnum = f.ints("VISITNUM")?;
    let visit = f.chars("VISIT")?;
    let has_ranges = f.has(&col("STNRLO"));
    let (lo, hi) = if has_ranges {
        (f.nums(&col("STNRLO"))?, f.nums(&col("STNRHI"))?)
    } else {
        (vec![None; f.n()], vec![None; f.n()])
    };
    let cat = if src == "LB" {
        Some(f.chars("LBCAT")?)
    } else {
        None
    };

    let trtsdt = sl.dates("TRTSDT")?;
    let trtedt = sl.dates("TRTEDT")?;
    let windows = planned(tv)?;

    // Records with a value, sorted by subject, parameter order, date, sequence.
    let mut param_order: Vec<&str> = Vec::new();
    for t in testcd.iter().flatten() {
        if !param_order.contains(&t.as_str()) {
            param_order.push(t);
        }
    }
    let pn = |t: &Option<String>| {
        t.as_deref()
            .and_then(|t| param_order.iter().position(|p| *p == t))
            .unwrap_or(usize::MAX)
    };
    let mut rows: Vec<usize> = (0..f.n()).filter(|&i| stresn[i].is_some()).collect();
    rows.sort_by(|&a, &b| {
        (&usubjid[a], pn(&testcd[a]), dt(&dtc[a]), seq[a]).cmp(&(
            &usubjid[b],
            pn(&testcd[b]),
            dt(&dtc[b]),
            seq[b],
        ))
    });

    struct Rec {
        i: usize,
        sl: Option<usize>,
        adt: Option<Date>,
        avisitn: Option<i64>,
        avisit: Option<String>,
        target: Option<i64>,
        ablfl: bool,
        base: Option<f64>,
        bnrind: Option<String>,
        anrind: Option<String>,
        anl: bool,
        post: bool,
    }
    let nrind = |v: Option<f64>, lo: Option<f64>, hi: Option<f64>| match (v, lo, hi) {
        (Some(v), Some(lo), Some(hi)) => Some(
            if v < lo {
                "LOW"
            } else if v > hi {
                "HIGH"
            } else {
                "NORMAL"
            }
            .to_string(),
        ),
        _ => None,
    };

    let mut recs: Vec<Rec> = rows
        .iter()
        .map(|&i| {
            let sl_row = usubjid[i].as_deref().and_then(|u| sl_idx.get(u).copied());
            let adt = dt(&dtc[i]);
            let sdt = sl_row.and_then(|r| trtsdt[r]);
            let post = matches!((adt, sdt), (Some(a), Some(s)) if a > s);
            let ady = rel_day(adt, sdt);
            // Scheduled visits keep their own window; others go to the nearest planned visit.
            let scheduled = windows.iter().find(|w| Some(w.0) == visitnum[i]);
            let (avisitn, avisit, target) = match (scheduled, post, ady) {
                (Some(w), _, _) => (Some(w.0), Some(w.1.clone()), Some(w.2)),
                (None, true, Some(d)) if !windows.is_empty() => {
                    let w = windows
                        .iter()
                        .min_by_key(|w| ((w.2 - d).abs(), -w.2))
                        .expect("non-empty");
                    (Some(w.0), Some(w.1.clone()), Some(w.2))
                }
                _ => (visitnum[i], visit[i].clone(), None),
            };
            Rec {
                i,
                sl: sl_row,
                adt,
                avisitn,
                avisit,
                target,
                ablfl: false,
                base: None,
                bnrind: None,
                anrind: nrind(stresn[i], lo[i], hi[i]),
                anl: false,
                post,
            }
        })
        .collect();

    // Baseline and analysis flags, per subject and parameter.
    let mut start = 0;
    while start < recs.len() {
        let key = (&usubjid[recs[start].i], &testcd[recs[start].i]);
        let mut end = start;
        while end < recs.len() && (&usubjid[recs[end].i], &testcd[recs[end].i]) == key {
            end += 1;
        }
        let group = &mut recs[start..end];
        let base_pos = group.iter().rposition(|r| !r.post && r.adt.is_some());
        if let Some(b) = base_pos {
            group[b].ablfl = true;
            group[b].anl = true;
            let (bv, bn) = (stresn[group[b].i], group[b].anrind.clone());
            for r in group.iter_mut() {
                r.base = bv;
                r.bnrind = bn.clone();
            }
        }
        // Closest to target per post-baseline window; unplanned-window records compete too.
        let mut best: HashMap<i64, (i64, usize)> = HashMap::new();
        for (k, r) in group.iter().enumerate() {
            if !r.post {
                continue;
            }
            let (Some(av), Some(adt)) = (r.avisitn, r.adt) else {
                continue;
            };
            let sdt = r.sl.and_then(|s| trtsdt[s]);
            let ady = rel_day(Some(adt), sdt).unwrap_or(0);
            let dist = r.target.map_or(0, |t| (ady - t).abs());
            match best.get(&av) {
                Some((d, _)) if *d < dist => {}
                _ => {
                    best.insert(av, (dist, k));
                }
            }
        }
        for (_, k) in best.values() {
            group[*k].anl = true;
        }
        start = end;
    }

    let n = recs.len();
    let sl_char = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.chars(name)?;
        Ok(Values::Char(
            recs.iter()
                .map(|r| r.sl.and_then(|s| v[s].clone()))
                .collect(),
        ))
    };
    let sl_num = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.nums(name)?;
        Ok(Values::Num(
            recs.iter().map(|r| r.sl.and_then(|s| v[s])).collect(),
        ))
    };
    let chg: Vec<Option<f64>> = recs
        .iter()
        .map(|r| match (r.post, stresn[r.i], r.base) {
            (true, Some(v), Some(b)) => Some(math::round_dp(v - b, decimals(&stresc[r.i]))),
            _ => None,
        })
        .collect();
    let y = |b: bool| b.then(|| "Y".to_string());

    let mut cols: Vec<(&str, Values)> = vec![
            ("STUDYID", sl_char("STUDYID")?),
            ("USUBJID", Values::Char(recs.iter().map(|r| usubjid[r.i].clone()).collect())),
            ("SUBJID", sl_char("SUBJID")?),
            ("SITEID", sl_char("SITEID")?),
            ("AGE", sl_num("AGE")?),
            ("SEX", sl_char("SEX")?),
            ("RACE", sl_char("RACE")?),
            ("SAFFL", sl_char("SAFFL")?),
            ("TRTP", sl_char("TRT01P")?),
            ("TRTPN", sl_num("TRT01PN")?),
            ("TRTA", sl_char("TRT01A")?),
            ("TRTAN", sl_num("TRT01AN")?),
            ("PARAMCD", Values::Char(recs.iter().map(|r| testcd[r.i].clone()).collect())),
            (
                "PARAM",
                Values::Char(
                    recs.iter()
                        .map(|r| match (&test[r.i], &stresu[r.i]) {
                            (Some(t), Some(u)) => Some(format!("{t} ({u})")),
                            (t, _) => t.clone(),
                        })
                        .collect(),
                ),
            ),
            ("PARAMN", Values::Num(recs.iter().map(|r| Some(pn(&testcd[r.i]) as f64 + 1.0)).collect())),
            ("AVAL", Values::Num(recs.iter().map(|r| stresn[r.i]).collect())),
            ("BASE", Values::Num(recs.iter().map(|r| r.base).collect())),
            ("CHG", Values::Num(chg.clone())),
            (
                "PCHG",
                Values::Num(
                    // As admiral's derive_var_pchg: (AVAL - BASE) / |BASE| * 100, unrounded.
                    recs.iter()
                        .map(|r| match (r.post, stresn[r.i], r.base) {
                            (true, Some(v), Some(b)) if b != 0.0 => Some((v - b) / b.abs() * 100.0),
                            _ => None,
                        })
                        .collect(),
                ),
            ),
            ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
            ("ADY", Values::Int(recs.iter().map(|r| rel_day(r.adt, r.sl.and_then(|s| trtsdt[s]))).collect())),
            ("AVISIT", Values::Char(recs.iter().map(|r| r.avisit.clone()).collect())),
            ("AVISITN", Values::Num(recs.iter().map(|r| r.avisitn.map(|v| v as f64)).collect())),
            ("ABLFL", Values::Char(recs.iter().map(|r| y(r.ablfl)).collect())),
            (
                "ONTRTFL",
                Values::Char(
                    recs.iter()
                        .map(|r| {
                            let (s, e) = (r.sl.and_then(|s| trtsdt[s]), r.sl.and_then(|s| trtedt[s]));
                            y(matches!((r.adt, s, e), (Some(a), Some(s), Some(e)) if a > s && a <= e.add_days(1)))
                        })
                        .collect(),
                ),
            ),
            ("ANL01FL", Values::Char(recs.iter().map(|r| y(r.anl)).collect())),
            ("SRCDOM", Values::Char(vec![Some(src.to_string()); n])),
            ("SRCVAR", Values::Char(vec![Some(col("STRESN")); n])),
            ("SRCSEQ", Values::Int(recs.iter().map(|r| seq[r.i]).collect())),
    ];
    // Category and reference-range variables only exist where the source has them (LB, not VS).
    if let Some(cat) = cat {
        cols.push((
            "PARCAT1",
            Values::Char(recs.iter().map(|r| cat[r.i].clone()).collect()),
        ));
    }
    if has_ranges {
        cols.extend(vec![
            ("ANRLO", Values::Num(recs.iter().map(|r| lo[r.i]).collect())),
            ("ANRHI", Values::Num(recs.iter().map(|r| hi[r.i]).collect())),
            (
                "ANRIND",
                Values::Char(recs.iter().map(|r| r.anrind.clone()).collect()),
            ),
            (
                "BNRIND",
                Values::Char(recs.iter().map(|r| r.bnrind.clone()).collect()),
            ),
            (
                "SHIFT1",
                Values::Char(
                    recs.iter()
                        .map(|r| match (r.post, &r.bnrind, &r.anrind) {
                            (true, Some(b), Some(a)) => Some(format!("{b} to {a}")),
                            _ => None,
                        })
                        .collect(),
                ),
            ),
        ]);
    }
    Ok(build(name, cols)?)
}
