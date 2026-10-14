//! Pharmacokinetic BDS datasets, following admiral 1.5's templates:
//!
//! - ADPC (`ad_adpc.R`) from PC, EX and ADSL: concentration records plus one `DOSE` record per
//!   administration; nominal and actual times relative to the first and to the reference dose;
//!   BLQ values imputed at half the LLOQ (`DTYPE = HALFLLOQ`, added records); samples taken at
//!   the next dose's nominal time copied as that dose's pre-dose record (`DTYPE = COPY`); and
//!   pre-dose baselines per reference dose (`BASETYPE`). Adaptation: a pre-dose record is a
//!   baseline only if it precedes its reference dose (`ARRLT` <= 0), which matters only for
//!   crossovers, where every period has a pre-dose sample.
//! - ADPP (`ad_adpp.R`) from PP and ADSL. Adaptation: `PARAMCD`/`PARAM` are `PPTESTCD`/`PPTEST`
//!   for every parameter (the template's lookup lists only some) and the date comes from
//!   `PPRFTDTC` (the template reads the older `PPRFDTC`).

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, DateTime, Values, math};

use crate::{DeriveError, Frame, dt, index, rel_day};

/// Hours of a planned time point, as admiral's `convert_xxtpt_to_hours()` for the forms
/// trialforge writes (`Pre-dose`, `0.5h Post-dose`, `24h Post-dose`).
fn tpt_hours(tpt: &str) -> Option<f64> {
    let t = tpt.trim();
    let lower = t.to_ascii_lowercase();
    if lower == "pre-dose" || lower == "predose" {
        return Some(0.0);
    }
    let num = lower.strip_suffix(" post-dose")?.strip_suffix('h')?;
    num.trim().parse().ok()
}

/// R's `prettyNum()` of a number (default 7 significant digits; scientific when shorter).
fn pretty(x: f64) -> String {
    let fixed = {
        let s = format!("{}", math::signif(x, 7));
        s.strip_suffix(".0").map_or(s.clone(), str::to_string)
    };
    let e = decade(x.abs());
    let mantissa = math::signif(x / math::pow(10.0, f64::from(e)), 7);
    let m = format!("{mantissa}");
    let m = m.strip_suffix(".0").map_or(m.clone(), str::to_string);
    let sci = format!("{m}e{}{:02}", if e < 0 { "-" } else { "+" }, e.abs());
    if x != 0.0 && sci.len() < fixed.len() {
        sci
    } else {
        fixed
    }
}

/// ⌊log10 x⌋ for x > 0 (0 for 0), by repeated scaling (exact for decimal powers).
fn decade(x: f64) -> i32 {
    if x == 0.0 {
        return 0;
    }
    let mut e = 0i32;
    let mut v = x;
    while v >= 10.0 {
        v /= 10.0;
        e += 1;
    }
    while v < 1.0 {
        v *= 10.0;
        e -= 1;
    }
    e
}

/// One ADPC record.
#[derive(Clone, Debug, Default)]
struct Rec {
    sl: usize,
    evid: u8,
    paramcd: String,
    param: String,
    parcat1: Option<String>,
    aval: Option<f64>,
    avalu: Option<String>,
    avalcat1: Option<String>,
    adtm: Option<DateTime>,
    nfrlt: Option<f64>,
    avisitn: Option<f64>,
    atpt: Option<String>,
    atptn: Option<f64>,
    atptref: Option<String>,
    basetype: Option<String>,
    base: Option<f64>,
    dtype: Option<&'static str>,
    ablfl: bool,
    afrlt: Option<f64>,
    arrlt: Option<f64>,
    nrrlt: Option<f64>,
    alloq: Option<f64>,
    dosea: Option<f64>,
    pcrftdtm: Option<DateTime>,
    fanldtm: Option<DateTime>,
    srcdom: &'static str,
    srcseq: Option<i64>,
    // For the COPY records.
    nxrlt: Option<f64>,
    axrlt: Option<f64>,
    next_visit: Option<String>,
    next_adtm: Option<DateTime>,
    next_dose: Option<f64>,
    blq: bool,
}

/// A dose (EX with EXDOSE > 0).
#[derive(Clone)]
struct Dose {
    adtm: DateTime,
    nfrlt: f64,
    dose: f64,
    unit: Option<String>,
    drug: String,
    seq: Option<i64>,
}

fn avisit(nfrlt: f64) -> (f64, String) {
    let n = (nfrlt / 24.0).floor() + 1.0;
    (n, format!("Day {n}"))
}

/// Derives ADPC.
pub fn adpc(pc: &Dataset, ex: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let pc = Frame::new(pc);
    let ex = Frame::new(ex);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;

    // Doses per subject, in time order.
    let mut doses: HashMap<&str, Vec<Dose>> = HashMap::new();
    let (eu, edose, edtc, etrt, eunit, evisitdy, eseq) = (
        ex.chars("USUBJID")?,
        ex.nums("EXDOSE")?,
        ex.chars("EXSTDTC")?,
        ex.chars("EXTRT")?,
        ex.chars("EXDOSU")?,
        ex.ints("VISITDY")?,
        ex.ints("EXSEQ")?,
    );
    for i in 0..ex.n() {
        let (Some(u), Some(d), Some(adtm)) = (eu[i].as_deref(), edose[i], crate::dtm(&edtc[i]))
        else {
            continue;
        };
        if d <= 0.0 {
            continue;
        }
        doses.entry(u).or_default().push(Dose {
            adtm,
            nfrlt: evisitdy[i].map_or(0.0, |v| (v - 1) as f64 * 24.0),
            dose: d,
            unit: eunit[i].clone(),
            drug: etrt[i].clone().unwrap_or_default(),
            seq: eseq[i],
        });
    }
    for v in doses.values_mut() {
        v.sort_by_key(|d| (d.adtm, d.seq));
    }

    let (
        pu,
        ptestcd,
        ptest,
        pstresn,
        pstresc,
        pstresu,
        pspec,
        plloq,
        pdtc,
        ptpt,
        ptptnum,
        pvisitdy,
        pseq,
    ) = (
        pc.chars("USUBJID")?,
        pc.chars("PCTESTCD")?,
        pc.chars("PCTEST")?,
        pc.nums("PCSTRESN")?,
        pc.chars("PCSTRESC")?,
        pc.chars("PCSTRESU")?,
        pc.chars("PCSPEC")?,
        pc.nums("PCLLOQ")?,
        pc.chars("PCDTC")?,
        pc.chars("PCTPT")?,
        pc.nums("PCTPTNUM")?,
        pc.ints("VISITDY")?,
        pc.ints("PCSEQ")?,
    );

    let mut per: Vec<Vec<Rec>> = vec![Vec::new(); sl.n()];
    for i in 0..pc.n() {
        let Some(u) = pu[i].as_deref() else { continue };
        let Some(&r) = sl_idx.get(u) else { continue };
        let drug = ptest[i].clone().unwrap_or_default();
        let Some(ds) = doses.get(u) else { continue };
        // First dose of this analyte.
        let Some(first) = ds.iter().filter(|d| d.drug == drug).map(|d| d.adtm).min() else {
            continue;
        };
        let adtm = crate::dtm(&pdtc[i]);
        let nfrlt = match (pvisitdy[i], ptpt[i].as_deref().and_then(tpt_hours)) {
            (Some(v), Some(h)) => Some((v - 1) as f64 * 24.0 + h),
            _ => None,
        };
        let prev = adtm.and_then(|a| ds.iter().filter(|d| d.adtm < a).max_by_key(|d| d.adtm));
        let next = adtm.and_then(|a| ds.iter().filter(|d| d.adtm >= a).min_by_key(|d| d.adtm));
        let nom_prev = nfrlt.and_then(|n| {
            ds.iter()
                .map(|d| d.nfrlt)
                .filter(|x| *x < n)
                .fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x))))
        });
        let nom_next = nfrlt.and_then(|n| {
            ds.iter()
                .map(|d| d.nfrlt)
                .filter(|x| *x >= n)
                .fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.min(x))))
        });
        let hours = |a: DateTime, b: DateTime| a.hours_since(b);
        let blq = pstresc[i].as_deref() == Some("<BLQ") && pstresn[i].is_none();
        let reference = prev.or(next);
        per[r].push(Rec {
            sl: r,
            evid: 0,
            paramcd: ptestcd[i].clone().unwrap_or_default(),
            param: format!("{} concentration", drug),
            parcat1: pspec[i].clone(),
            aval: pstresn[i],
            avalu: pstresu[i].clone(),
            avalcat1: if pstresc[i].as_deref() == Some("<BLQ") {
                pstresc[i].clone()
            } else {
                pstresn[i].map(|v| pretty(math::signif(v, 3)))
            },
            adtm,
            nfrlt,
            avisitn: nfrlt.map(|n| avisit(n).0),
            atpt: ptpt[i].clone(),
            atptn: ptptnum[i],
            atptref: reference.map(|d| avisit(d.nfrlt).1),
            dtype: None,
            afrlt: adtm.map(|a| hours(a, first)),
            arrlt: adtm.and_then(|a| {
                prev.map(|p| hours(a, p.adtm))
                    .or_else(|| next.map(|n| hours(a, n.adtm)))
            }),
            nrrlt: nfrlt.map(|n| n - nom_prev.unwrap_or(f64::NAN)),
            alloq: plloq[i],
            dosea: reference.map(|d| d.dose),
            pcrftdtm: reference.map(|d| d.adtm),
            fanldtm: Some(first),
            srcdom: "PC",
            srcseq: pseq[i],
            nxrlt: nfrlt.zip(nom_next).map(|(n, x)| n - x),
            axrlt: adtm.zip(next).map(|(a, n)| hours(a, n.adtm)),
            next_visit: next.map(|d| avisit(d.nfrlt).1),
            next_adtm: next.map(|d| d.adtm),
            next_dose: next.map(|d| d.dose),
            blq,
            ..Rec::default()
        });
    }

    let usubjid = sl.chars("USUBJID")?;
    let mut out: Vec<Rec> = Vec::new();
    for (r, mut recs) in per.into_iter().enumerate() {
        if recs.is_empty() {
            continue;
        }
        let u = usubjid[r].as_deref().unwrap_or("");
        // NRRLT with no previous nominal dose: NFRLT minus the smallest previous-dose NFRLT.
        let min_prev = recs
            .iter()
            .filter_map(|x| match (x.nfrlt, x.nrrlt) {
                (Some(n), Some(rr)) if !rr.is_nan() => Some(n - rr),
                _ => None,
            })
            .fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.min(x))));
        for x in &mut recs {
            if x.nrrlt.is_some_and(f64::is_nan) {
                x.nrrlt = x.nfrlt.zip(min_prev).map(|(n, m)| n - m);
            }
        }
        // Dose records up to the last sample date.
        let maxdate = recs.iter().filter_map(|x| x.adtm.map(DateTime::date)).max();
        let first = recs.iter().filter_map(|x| x.fanldtm).min();
        for d in doses.get(u).into_iter().flatten() {
            if maxdate.is_some_and(|m| d.adtm.date() > m) {
                continue;
            }
            let (n, v) = avisit(d.nfrlt);
            recs.push(Rec {
                sl: r,
                evid: 1,
                paramcd: "DOSE".into(),
                param: format!("{} dose", d.drug),
                aval: Some(d.dose),
                avalu: d.unit.clone(),
                adtm: Some(d.adtm),
                nfrlt: Some(d.nfrlt),
                avisitn: Some(n),
                atpt: Some("Dose".into()),
                atptn: Some(0.0),
                atptref: Some(v),
                afrlt: first.map(|f| d.adtm.hours_since(f)),
                arrlt: Some(0.0),
                nrrlt: Some(0.0),
                dosea: Some(d.dose),
                pcrftdtm: Some(d.adtm),
                fanldtm: first,
                srcdom: "EX",
                srcseq: d.seq,
                nxrlt: Some(0.0),
                ..Rec::default()
            });
        }
        // HALFLLOQ records for BLQ samples.
        let halves: Vec<Rec> = recs
            .iter()
            .filter(|x| x.blq)
            .map(|x| Rec {
                aval: x.alloq.map(|l| l * 0.5),
                dtype: Some("HALFLLOQ"),
                ..x.clone()
            })
            .collect();
        recs.extend(halves);
        // COPY records: samples at the next dose's nominal time become its pre-dose records.
        let copies: Vec<Rec> = recs
            .iter()
            .filter(|x| {
                x.evid == 0
                    && x.nfrlt.is_some_and(|n| n > 0.0)
                    && x.nxrlt == Some(0.0)
                    && x.next_visit.is_some()
            })
            .map(|x| Rec {
                atptref: x.next_visit.clone(),
                arrlt: x.axrlt,
                nrrlt: x.nxrlt,
                pcrftdtm: x.next_adtm,
                dosea: x.next_dose,
                atpt: Some("Pre-dose".into()),
                atptn: Some(-0.5),
                dtype: Some(if x.dtype.is_some() {
                    "COPY/HALFLLOQ"
                } else {
                    "COPY"
                }),
                ..x.clone()
            })
            .collect();
        recs.extend(copies);
        for x in &mut recs {
            x.basetype = x.atptref.as_ref().map(|a| format!("{a} Baseline"));
            // Adaptation: only pre-dose records taken before their reference dose, so that a
            // crossover's later pre-dose samples are baselines of their own period (through
            // their COPY records) and not of the previous one. No effect on repeated dosing.
            x.ablfl = x.atpt.as_deref() == Some("Pre-dose")
                && x.aval.is_some()
                && x.arrlt.is_none_or(|a| a <= 0.0);
        }
        // Baseline per parameter, specimen and baseline type.
        let mut bases: HashMap<(String, Option<String>, Option<String>), f64> = HashMap::new();
        for x in recs.iter().filter(|x| x.ablfl) {
            if let Some(a) = x.aval {
                bases.insert(
                    (x.paramcd.clone(), x.parcat1.clone(), x.basetype.clone()),
                    a,
                );
            }
        }
        for x in &mut recs {
            x.base = bases
                .get(&(x.paramcd.clone(), x.parcat1.clone(), x.basetype.clone()))
                .copied();
        }
        // ASEQ order: ADTM, BASETYPE, EVID, AVISITN, ATPTN, PARCAT1, DTYPE (missing last).
        let key = |x: &Rec| {
            (
                (x.adtm.is_none(), x.adtm),
                (x.basetype.is_none(), x.basetype.clone()),
                x.evid,
                (x.avisitn.is_none(), x.avisitn.map(|v| (v * 1000.0) as i64)),
                (x.atptn.is_none(), x.atptn.map(|v| (v * 1000.0) as i64)),
                (x.parcat1.is_none(), x.parcat1.clone()),
                (x.dtype.is_none(), x.dtype),
            )
        };
        recs.sort_by_key(key);
        out.extend(recs);
    }

    let n = out.len();
    let mut aseq = Vec::with_capacity(n);
    for (i, x) in out.iter().enumerate() {
        let k = if i > 0 && out[i - 1].sl == x.sl {
            aseq.last().copied().flatten().unwrap_or(0) + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }
    let sl_rows: Vec<usize> = out.iter().map(|x| x.sl).collect();
    let copy_sl = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.chars(name)?;
        Ok(Values::Char(
            sl_rows.iter().map(|&r| v[r].clone()).collect(),
        ))
    };
    let y = |b: bool| b.then(|| "Y".to_string());
    let adt: Vec<Option<Date>> = out.iter().map(|x| x.adtm.map(DateTime::date)).collect();
    let mut cols: Vec<(&str, Values)> = vec![
        ("STUDYID", copy_sl("STUDYID")?),
        ("USUBJID", copy_sl("USUBJID")?),
        ("SUBJID", copy_sl("SUBJID")?),
        ("SITEID", copy_sl("SITEID")?),
        ("TRT01P", copy_sl("TRT01P")?),
        ("TRT01A", copy_sl("TRT01A")?),
        ("TRTP", copy_sl("TRT01P")?),
        ("TRTA", copy_sl("TRT01A")?),
        ("ASEQ", Values::Int(aseq)),
        (
            "PARAMCD",
            Values::Char(out.iter().map(|x| Some(x.paramcd.clone())).collect()),
        ),
        (
            "PARAM",
            Values::Char(out.iter().map(|x| Some(x.param.clone())).collect()),
        ),
        (
            "PARCAT1",
            Values::Char(out.iter().map(|x| x.parcat1.clone()).collect()),
        ),
        ("AVAL", Values::Num(out.iter().map(|x| x.aval).collect())),
        (
            "AVALU",
            Values::Char(out.iter().map(|x| x.avalu.clone()).collect()),
        ),
        (
            "AVALCAT1",
            Values::Char(out.iter().map(|x| x.avalcat1.clone()).collect()),
        ),
        (
            "BASETYPE",
            Values::Char(out.iter().map(|x| x.basetype.clone()).collect()),
        ),
        ("BASE", Values::Num(out.iter().map(|x| x.base).collect())),
        (
            "CHG",
            Values::Num(
                out.iter()
                    .map(|x| x.aval.zip(x.base).map(|(a, b)| a - b))
                    .collect(),
            ),
        ),
        (
            "DTYPE",
            Values::Char(out.iter().map(|x| x.dtype.map(str::to_string)).collect()),
        ),
        (
            "ABLFL",
            Values::Char(out.iter().map(|x| y(x.ablfl)).collect()),
        ),
        ("ANL01FL", Values::Char(vec![Some("Y".into()); n])),
        (
            "ANL02FL",
            Values::Char(out.iter().map(|x| y(x.dtype.is_none())).collect()),
        ),
        (
            "ADTM",
            Values::DateTime(out.iter().map(|x| x.adtm).collect()),
        ),
        ("ADT", Values::Date(adt.clone())),
        (
            "ADY",
            Values::Int(
                adt.iter()
                    .zip(&sl_rows)
                    .map(|(a, &r)| rel_day(*a, trtsdt[r]))
                    .collect(),
            ),
        ),
        (
            "AVISIT",
            Values::Char(
                out.iter()
                    .map(|x| x.avisitn.map(|v| format!("Day {v}")))
                    .collect(),
            ),
        ),
        (
            "AVISITN",
            Values::Num(out.iter().map(|x| x.avisitn).collect()),
        ),
        (
            "ATPT",
            Values::Char(out.iter().map(|x| x.atpt.clone()).collect()),
        ),
        ("ATPTN", Values::Num(out.iter().map(|x| x.atptn).collect())),
        (
            "ATPTREF",
            Values::Char(out.iter().map(|x| x.atptref.clone()).collect()),
        ),
        ("NFRLT", Values::Num(out.iter().map(|x| x.nfrlt).collect())),
        ("AFRLT", Values::Num(out.iter().map(|x| x.afrlt).collect())),
        ("NRRLT", Values::Num(out.iter().map(|x| x.nrrlt).collect())),
        ("ARRLT", Values::Num(out.iter().map(|x| x.arrlt).collect())),
        (
            "MRRLT",
            Values::Num(out.iter().map(|x| x.arrlt.map(|a| a.max(0.0))).collect()),
        ),
        ("FRLTU", Values::Char(vec![Some("HOURS".into()); n])),
        ("RRLTU", Values::Char(vec![Some("HOURS".into()); n])),
        ("ALLOQ", Values::Num(out.iter().map(|x| x.alloq).collect())),
        ("DOSEA", Values::Num(out.iter().map(|x| x.dosea).collect())),
        ("DOSEU", Values::Char(vec![Some("mg".into()); n])),
        (
            "PCRFTDTM",
            Values::DateTime(out.iter().map(|x| x.pcrftdtm).collect()),
        ),
        (
            "FANLDTM",
            Values::DateTime(out.iter().map(|x| x.fanldtm).collect()),
        ),
        (
            "SRCDOM",
            Values::Char(out.iter().map(|x| Some(x.srcdom.to_string())).collect()),
        ),
        ("SRCVAR", Values::Char(vec![Some("SEQ".into()); n])),
        (
            "SRCSEQ",
            Values::Int(out.iter().map(|x| x.srcseq).collect()),
        ),
    ];
    let periods = crate::Periods::new(adsl)?;
    let keys: Vec<(usize, Option<Date>)> =
        sl_rows.iter().copied().zip(adt.iter().copied()).collect();
    crate::by_period(&mut cols, &periods, &keys);
    Ok(build("ADPC", cols)?)
}

/// Derives ADPP.
pub fn adpp(pp: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let pp = Frame::new(pp);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;
    let u = pp.chars("USUBJID")?;
    let rows: Vec<(usize, usize)> = (0..pp.n())
        .filter_map(|i| u[i].as_deref().and_then(|x| sl_idx.get(x)).map(|&r| (i, r)))
        .collect();
    let adt: Vec<Option<Date>> = {
        let v = pp.chars("PPRFTDTC")?;
        rows.iter().map(|&(i, _)| dt(&v[i])).collect()
    };
    let ady: Vec<Option<i64>> = adt
        .iter()
        .zip(&rows)
        .map(|(a, &(_, r))| rel_day(*a, trtsdt[r]))
        .collect();
    let adt_keys = adt.clone();
    let from_pp = |name: &str| -> Result<Values, DeriveError> {
        let v = pp.chars(name)?;
        Ok(Values::Char(
            rows.iter().map(|&(i, _)| v[i].clone()).collect(),
        ))
    };
    let from_sl = |name: &str| -> Result<Values, DeriveError> {
        let v = sl.chars(name)?;
        Ok(Values::Char(
            rows.iter().map(|&(_, r)| v[r].clone()).collect(),
        ))
    };
    let stresn = pp.nums("PPSTRESN")?;
    let seq = pp.ints("PPSEQ")?;
    let n = rows.len();
    let mut cols: Vec<(&str, Values)> = vec![
        ("STUDYID", from_sl("STUDYID")?),
        ("USUBJID", from_sl("USUBJID")?),
        ("SUBJID", from_sl("SUBJID")?),
        ("SITEID", from_sl("SITEID")?),
        ("TRT01P", from_sl("TRT01P")?),
        ("TRT01A", from_sl("TRT01A")?),
        ("TRTP", from_sl("TRT01P")?),
        ("TRTA", from_sl("TRT01A")?),
        ("PARAMCD", from_pp("PPTESTCD")?),
        ("PARAM", from_pp("PPTEST")?),
        ("PARCAT1", from_pp("PPCAT")?),
        (
            "AVAL",
            Values::Num(rows.iter().map(|&(i, _)| stresn[i]).collect()),
        ),
        ("AVALU", from_pp("PPSTRESU")?),
        ("ADT", Values::Date(adt)),
        ("ADY", Values::Int(ady.clone())),
        (
            "AVISIT",
            Values::Char(ady.iter().map(|d| d.map(|d| format!("Day {d}"))).collect()),
        ),
        (
            "AVISITN",
            Values::Num(ady.iter().map(|d| d.map(|d| d as f64)).collect()),
        ),
        ("PPSTINT", from_pp("PPSTINT")?),
        ("PPENINT", from_pp("PPENINT")?),
        ("SRCDOM", Values::Char(vec![Some("PP".into()); n])),
        ("SRCVAR", Values::Char(vec![Some("SEQ".into()); n])),
        (
            "SRCSEQ",
            Values::Int(rows.iter().map(|&(i, _)| seq[i]).collect()),
        ),
    ];
    let periods = crate::Periods::new(adsl)?;
    let keys: Vec<(usize, Option<Date>)> = rows.iter().map(|&(_, r)| r).zip(adt_keys).collect();
    crate::by_period(&mut cols, &periods, &keys);
    Ok(build("ADPP", cols)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn time_points_and_numbers_format_as_r() {
        assert_eq!(tpt_hours("Pre-dose"), Some(0.0));
        assert_eq!(tpt_hours("0.25h Post-dose"), Some(0.25));
        assert_eq!(tpt_hours("72h Post-dose"), Some(72.0));
        assert_eq!(pretty(66.4), "66.4");
        assert_eq!(pretty(1072.0), "1072");
        assert_eq!(pretty(0.999), "0.999");
        assert_eq!(pretty(100000.0), "1e+05");
        assert_eq!(pretty(123000.0), "123000");
    }
}
