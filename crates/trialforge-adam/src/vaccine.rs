//! Vaccine analysis datasets, following admiralvaccine's templates (`ad_adis.R`, `ad_adce.R`,
//! `ad_adface.R`), adapted only where noted (the adaptations are repeated in
//! `scripts/admiral_equivalence.R`).
//!
//! - **ADIS** from IS: for each titre the original value and its LOG10, 4FOLD and LOG10 4FOLD
//!   parameters (admiralvaccine `derive_var_aval_adis()`: below the LLOQ, LLOQ / 2 or, for
//!   4FOLD, the LLOQ; above the ULOQ, the ULOQ; rounded to 2 decimals), `DTYPE`, baseline,
//!   `CHG`, the fold rise `R2BASE`, `SERCAT1`, `BASECAT1` and `CRIT1`. Adaptations: the
//!   baseline is the `ISBLFL` record (the template hard-codes visit 10); the log parameters
//!   compare with log10 of the natural-scale thresholds (the template compares log titres with
//!   the LLOQ itself, so every log record is seronegative); no `BASETYPE` copies.
//! - **ADCE** from the reactogenicity records of CE: analysis dates, the vaccination period,
//!   `ASEV`/`ASEVN`, `AOCC01FL` (the most severe reaction per period) and the duration.
//! - **ADFACE** from FACE and the diary temperatures in VS: fever records (38.0 °C or more,
//!   `derive_fever_records()`), severities of redness and swelling from their diameters
//!   (`derive_diam_to_sev_records()`: ≤ 2 cm none, ≤ 5 mild, ≤ 10 moderate, more severe),
//!   maximum severity, diameter and temperature per vaccination (`derive_extreme_records()`,
//!   the latest record on ties), parameters by test and object (`derive_vars_params()`), the
//!   maximum flags `ANL02FL`/`ANL03FL` (`derive_vars_max_flag()`), the event flags
//!   `EVENTFL`/`EVENTDFL` (`derive_vars_event_flag()`, cut-off 2.5) and the blanking of
//!   derived records' `FA` variables (`post_process_reacto()`). Adaptations: the parameter
//!   lookup covers every reaction this archetype records, and fever records have an `ADY`
//!   (the template merges the reference date before adding them).

use std::collections::{BTreeMap, HashMap};

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values, math};

use crate::onco::adsl_columns;
use crate::ophtha::title_case;
use crate::{DeriveError, Frame, Periods, dt, index, rel_day};

/// `--CAT` of reactogenicity records.
const REACTOGENICITY: &str = "REACTOGENICITY";

fn y(b: bool) -> Option<String> {
    b.then(|| "Y".to_string())
}

fn text(rows: usize, f: impl Fn(usize) -> Option<String>) -> Values {
    Values::Char((0..rows).map(f).collect())
}

fn num(rows: usize, f: impl Fn(usize) -> Option<f64>) -> Values {
    Values::Num((0..rows).map(f).collect())
}

/// The ADSL keys and demographics carried onto each record.
const ADSL_VARS: [&str; 8] = [
    "STUDYID", "USUBJID", "SUBJID", "SITEID", "AGE", "SEX", "RACE", "SAFFL",
];

// ---- ADIS ------------------------------------------------------------------------------------

struct IsRec {
    sl: usize,
    i: usize,
    derived: usize,
    aval: Option<f64>,
    adt: Option<Date>,
}

/// Derives ADIS.
pub fn adis(is: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let f = Frame::new(is);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;
    let vax: Vec<Vec<Option<Date>>> = ["VAX01DT", "VAX02DT"]
        .iter()
        .filter(|v| sl.has(v))
        .map(|v| sl.dates(v))
        .collect::<Result<_, _>>()?;
    let periods = Periods::new(adsl)?;
    let (u, seq, testcd, test, cat, orres, stresn, stresu, lloq, uloq, blfl, visitnum, visit, dtc) = (
        f.chars("USUBJID")?,
        f.ints("ISSEQ")?,
        f.chars("ISTESTCD")?,
        f.chars("ISTEST")?,
        f.chars("ISCAT")?,
        f.chars("ISORRES")?,
        f.nums("ISSTRESN")?,
        f.chars("ISSTRESU")?,
        f.nums("ISLLOQ")?,
        f.nums("ISULOQ")?,
        f.chars("ISBLFL")?,
        f.ints("VISITNUM")?,
        f.chars("VISIT")?,
        f.chars("ISDTC")?,
    );
    let is_log = |d: usize| d == 1 || d == 3;
    // derive_var_aval_adis(): the rule by the result's position against the limits.
    let value = |i: usize, d: usize| -> Option<f64> {
        let (lo, hi) = (lloq[i]?, uloq[i]?);
        let lower = if d < 2 { lo / 2.0 } else { lo };
        let has = |c: char| orres[i].as_deref().is_some_and(|s| s.contains(c));
        let natural = match stresn[i] {
            Some(v) if v < lo => lower,
            Some(v) if v < hi => v,
            Some(_) => hi,
            None if has('<') => lower,
            None if has('>') => hi,
            None => return None,
        };
        let v = if is_log(d) {
            math::log10(natural)
        } else {
            natural
        };
        Some(math::round_dp(v, 2))
    };
    let mut recs: Vec<IsRec> = Vec::new();
    let mut rows: Vec<usize> = (0..f.n()).collect();
    rows.sort_by(|&a, &b| (&u[a], seq[a]).cmp(&(&u[b], seq[b])));
    for d in 0..4 {
        for &i in &rows {
            let Some(&r) = u[i].as_deref().and_then(|x| sl_idx.get(x)) else {
                continue;
            };
            recs.push(IsRec {
                sl: r,
                i,
                derived: d,
                aval: value(i, d),
                adt: dt(&dtc[i]),
            });
        }
    }
    recs.sort_by(|a, b| (a.sl, a.derived, seq[a.i]).cmp(&(b.sl, b.derived, seq[b.i])));
    let paramcd = |r: &IsRec| {
        let cd = testcd[r.i].clone().unwrap_or_default();
        match r.derived {
            0 => cd,
            1 => format!("{cd}L"),
            2 => format!("{cd}F"),
            _ => format!("{}LF", &cd[..cd.len().min(6)]),
        }
    };
    let param = |r: &IsRec| {
        let t = test[r.i].clone().unwrap_or_default();
        match r.derived {
            0 => t,
            1 => format!("LOG10 ({t})"),
            2 => format!("4FOLD ({t})"),
            _ => format!("LOG10 4FOLD ({t})"),
        }
    };
    let threshold =
        |r: &IsRec| lloq[r.i].map(|l| if is_log(r.derived) { math::log10(l) } else { l });
    let ablfl = |r: &IsRec| blfl[r.i].as_deref() == Some("Y");
    // Baseline value and visit per (subject, parameter).
    let mut base: HashMap<(usize, usize), (Option<f64>, Option<i64>)> = HashMap::new();
    for r in recs.iter().filter(|r| ablfl(r)) {
        base.entry((r.sl, r.derived))
            .or_insert((r.aval, visitnum[r.i]));
    }
    let b = |r: &IsRec| base.get(&(r.sl, r.derived)).copied();
    let post =
        |r: &IsRec| matches!((visitnum[r.i], b(r)), (Some(v), Some((_, Some(bv)))) if v > bv);
    let n = recs.len();
    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let period = |r: &IsRec| periods.of(r.sl, r.adt);
    let atptref = |r: &IsRec| {
        let k = vax
            .iter()
            .enumerate()
            .filter(|(_, v)| matches!((v[r.sl], r.adt), (Some(d), Some(a)) if d < a))
            .map(|(k, _)| k + 1)
            .max()
            .unwrap_or(1);
        Some(format!("VACCINATION {k}"))
    };
    let sercat = |r: &IsRec| -> Option<&'static str> {
        if !ablfl(r) {
            return None;
        }
        Some(match (r.aval, threshold(r)) {
            (Some(a), Some(t)) if a < t => "S-",
            (Some(_), Some(_)) => "S+",
            _ => "UNKNOWN",
        })
    };
    let dtype = |r: &IsRec| -> Option<String> {
        if r.derived > 1 {
            return None;
        }
        let i = r.i;
        let has = |c: char| orres[i].as_deref().is_some_and(|s| s.contains(c));
        let (lo, hi) = (lloq[i]?, uloq[i]?);
        if stresn[i].is_some_and(|v| v < lo) || has('<') {
            Some("HALFLLOQ".into())
        } else if stresn[i].is_some_and(|v| v > hi) || has('>') {
            Some("ULOQ".into())
        } else {
            None
        }
    };
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_VARS)?;
    cols.extend([
        (
            "TRTP",
            text(n, |k| period(&recs[k]).and_then(|p| p.3.clone())),
        ),
        (
            "TRTA",
            text(n, |k| period(&recs[k]).and_then(|p| p.5.clone())),
        ),
        (
            "APERIOD",
            Values::Int((0..n).map(|k| period(&recs[k]).map(|p| p.0)).collect()),
        ),
        ("PARAMCD", text(n, |k| Some(paramcd(&recs[k])))),
        ("PARAM", text(n, |k| Some(param(&recs[k])))),
        (
            "PARAMN",
            num(n, |k| Some([1.0, 11.0, 21.0, 31.0][recs[k].derived])),
        ),
        ("PARCAT1", text(n, |k| cat[recs[k].i].clone())),
        ("AVAL", num(n, |k| recs[k].aval)),
        ("AVALU", text(n, |k| stresu[recs[k].i].clone())),
        ("DTYPE", text(n, |k| dtype(&recs[k]))),
        ("ABLFL", text(n, |k| y(ablfl(&recs[k])))),
        ("BASE", num(n, |k| b(&recs[k]).and_then(|x| x.0))),
        (
            "BASECAT1",
            text(n, |k| {
                let r = &recs[k];
                let cut = if is_log(r.derived) { 1.0 } else { 10.0 };
                b(r).and_then(|x| x.0).map(|v| {
                    if v < cut {
                        "Titer value < 1:10"
                    } else {
                        "Titer value >= 1:10"
                    }
                    .to_string()
                })
            }),
        ),
        (
            "CHG",
            num(n, |k| {
                let r = &recs[k];
                (post(r)).then(|| Some(r.aval? - b(r)?.0?)).flatten()
            }),
        ),
        (
            "R2BASE",
            num(n, |k| {
                let r = &recs[k];
                (post(r)).then(|| Some(r.aval? / b(r)?.0?)).flatten()
            }),
        ),
        ("SERCAT1", text(n, |k| sercat(&recs[k]).map(str::to_string))),
        (
            "SERCAT1N",
            num(n, |k| {
                sercat(&recs[k]).map(|s| match s {
                    "S-" => 1.0,
                    "S+" => 2.0,
                    _ => 3.0,
                })
            }),
        ),
        (
            "CRIT1",
            text(n, |k| recs[k].aval.map(|_| "Titer >= ISLLOQ".to_string())),
        ),
        (
            "CRIT1FL",
            text(n, |k| {
                let r = &recs[k];
                Some(if r.aval? >= threshold(r)? { "Y" } else { "N" }.to_string())
            }),
        ),
        (
            "CRIT1FN",
            num(n, |k| {
                let r = &recs[k];
                Some(if r.aval? >= threshold(r)? { 1.0 } else { 0.0 })
            }),
        ),
        ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
        (
            "ADY",
            Values::Int(recs.iter().map(|r| rel_day(r.adt, trtsdt[r.sl])).collect()),
        ),
        (
            "AVISIT",
            text(n, |k| visit[recs[k].i].as_deref().map(title_case)),
        ),
        ("AVISITN", num(n, |k| visitnum[recs[k].i].map(|v| v as f64))),
        ("ATPTREF", text(n, |k| atptref(&recs[k]))),
        ("SRCDOM", text(n, |_| Some("IS".into()))),
        ("SRCVAR", text(n, |_| Some("ISSTRESN".into()))),
        (
            "SRCSEQ",
            Values::Int(recs.iter().map(|r| seq[r.i]).collect()),
        ),
    ]);
    Ok(build("ADIS", cols)?)
}

// ---- ADCE ------------------------------------------------------------------------------------

/// Derives ADCE.
pub fn adce(ce: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let f = Frame::new(ce);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;
    let periods = Periods::new(adsl)?;
    let (u, seq, term, decod, cat, scat, occur, loc, lat, sev, rel, st, en, tptref) = (
        f.chars("USUBJID")?,
        f.ints("CESEQ")?,
        f.chars("CETERM")?,
        f.chars("CEDECOD")?,
        f.chars("CECAT")?,
        f.chars("CESCAT")?,
        f.chars("CEOCCUR")?,
        f.chars("CELOC")?,
        f.chars("CELAT")?,
        f.chars("CESEV")?,
        f.chars("CEREL")?,
        f.chars("CESTDTC")?,
        f.chars("CEENDTC")?,
        f.chars("CETPTREF")?,
    );
    struct Rec {
        sl: usize,
        i: usize,
        ast: Option<Date>,
        aen: Option<Date>,
        period: Option<(i64, Date, Date)>,
        asevn: Option<i64>,
        occ: bool,
    }
    let mut recs: Vec<Rec> = Vec::new();
    for i in 0..f.n() {
        if cat[i].as_deref() != Some(REACTOGENICITY) {
            continue;
        }
        let Some(&r) = u[i].as_deref().and_then(|x| sl_idx.get(x)) else {
            continue;
        };
        let ast = dt(&st[i]);
        recs.push(Rec {
            sl: r,
            i,
            ast,
            aen: dt(&en[i]),
            period: periods.of(r, ast).map(|p| (p.0, p.1, p.2)),
            asevn: match sev[i].as_deref() {
                Some("MILD") => Some(1),
                Some("MODERATE") => Some(2),
                Some("SEVERE") => Some(3),
                Some("DEATH THREATENING") => Some(4),
                _ => None,
            },
            occ: false,
        });
    }
    // AOCC01FL: the most severe (then earliest, by term, by sequence) per subject and period.
    let mut best: HashMap<(usize, i64), usize> = HashMap::new();
    for (k, r) in recs.iter().enumerate() {
        let (Some(p), Some(_)) = (r.period, r.asevn) else {
            continue;
        };
        let key = |r: &Rec| {
            (
                -r.asevn.unwrap_or(0),
                rel_day(r.ast, trtsdt[r.sl]),
                decod[r.i].clone(),
                seq[r.i],
            )
        };
        match best.get(&(r.sl, p.0)) {
            Some(&b) if key(&recs[b]) <= key(r) => {}
            _ => {
                best.insert((r.sl, p.0), k);
            }
        }
    }
    for k in best.into_values() {
        recs[k].occ = true;
    }
    // ASEQ: by CEDECOD, CELAT, CETPTREF, APERIOD within subject (missing last).
    let last = |x: &Option<String>| (x.is_none(), x.clone());
    recs.sort_by(|a, b| {
        (
            a.sl,
            &decod[a.i],
            last(&lat[a.i]),
            &tptref[a.i],
            a.period.map(|p| p.0).is_none(),
            a.period.map(|p| p.0),
        )
            .cmp(&(
                b.sl,
                &decod[b.i],
                last(&lat[b.i]),
                &tptref[b.i],
                b.period.map(|p| p.0).is_none(),
                b.period.map(|p| p.0),
            ))
    });
    let mut aseq = Vec::with_capacity(recs.len());
    for (k, r) in recs.iter().enumerate() {
        let v = if k > 0 && recs[k - 1].sl == r.sl {
            aseq[k - 1] + 1
        } else {
            1
        };
        aseq.push(v);
    }
    let n = recs.len();
    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let chr = |v: &[Option<String>]| text(n, |k| v[recs[k].i].clone());
    let dur = |r: &Rec| match (r.ast, r.aen) {
        (Some(a), Some(e)) => Some(f64::from(e.days_since(a) + 1)),
        _ => None,
    };
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_VARS)?;
    cols.extend(adsl_columns(
        &sl,
        &sl_rows,
        &["TRT01P", "TRT01A", "TRTSDT", "TRTEDT"],
    )?);
    cols.extend([
        (
            "CESEQ",
            Values::Int(recs.iter().map(|r| seq[r.i]).collect()),
        ),
        ("CETERM", chr(term)),
        ("CEDECOD", chr(decod)),
        ("CECAT", chr(cat)),
        ("CESCAT", chr(scat)),
        ("CEOCCUR", chr(occur)),
        ("CELOC", chr(loc)),
        ("CELAT", chr(lat)),
        ("CETPTREF", chr(tptref)),
        ("ASEQ", Values::Int(aseq.into_iter().map(Some).collect())),
        ("ASTDT", Values::Date(recs.iter().map(|r| r.ast).collect())),
        ("AENDT", Values::Date(recs.iter().map(|r| r.aen).collect())),
        (
            "ASTDY",
            Values::Int(recs.iter().map(|r| rel_day(r.ast, trtsdt[r.sl])).collect()),
        ),
        (
            "AENDY",
            Values::Int(recs.iter().map(|r| rel_day(r.aen, trtsdt[r.sl])).collect()),
        ),
        (
            "APERIOD",
            Values::Int(recs.iter().map(|r| r.period.map(|p| p.0)).collect()),
        ),
        (
            "APERSDT",
            Values::Date(recs.iter().map(|r| r.period.map(|p| p.1)).collect()),
        ),
        (
            "APEREDT",
            Values::Date(recs.iter().map(|r| r.period.map(|p| p.2)).collect()),
        ),
        (
            "APERSTDY",
            Values::Int(
                recs.iter()
                    .map(|r| Some(i64::from(r.ast?.days_since(r.period?.1)) + 1))
                    .collect(),
            ),
        ),
        ("AREL", chr(rel)),
        ("ASEV", chr(sev)),
        ("ASEVN", Values::Int(recs.iter().map(|r| r.asevn).collect())),
        ("AOCC01FL", text(n, |k| y(recs[k].occ))),
        ("ADURN", num(n, |k| dur(&recs[k]))),
        (
            "ADURU",
            text(n, |k| dur(&recs[k]).map(|_| "days".to_string())),
        ),
    ]);
    Ok(build("ADCE", cols)?)
}

// ---- ADFACE ----------------------------------------------------------------------------------

/// The template's lookup: (FATESTCD, FAOBJ, PARAMCD); PARAMN is the position.
const PARAMS: [(&str, &str, &str); 32] = [
    ("SEV", "REDNESS", "SEVREDN"),
    ("DIAMETER", "REDNESS", "DIARE"),
    ("MAXDIAM", "REDNESS", "MDIRE"),
    ("MAXTEMP", "FEVER", "MAXTEMP"),
    ("OCCUR", "FEVER", "OCFEVER"),
    ("OCCUR", "ERYTHEMA", "OCERYTH"),
    ("MAXSEV", "SWELLING", "MAXSWEL"),
    ("MAXSEV", "REDNESS", "MAXREDN"),
    ("MAXSEV", "FATIGUE", "MAXSFAT"),
    ("MAXSEV", "HEADACHE", "MAXSHEA"),
    ("MAXSEV", "NEW OR WORSENED JOINT PAIN", "MSEVNWJP"),
    ("MAXSEV", "NEW OR WORSENED MUSCLE PAIN", "MSEVNWMP"),
    ("OCCUR", "REDNESS", "OCISR"),
    ("OCCUR", "SWELLING", "OCINS"),
    ("OCCUR", "PAIN AT INJECTION SITE", "OCPIS"),
    ("OCCUR", "FATIGUE", "OCFATIG"),
    ("OCCUR", "HEADACHE", "OCHEAD"),
    ("OCCUR", "CHILLS", "OCCHILLS"),
    ("OCCUR", "DIARRHEA", "OCDIAR"),
    ("OCCUR", "NEW OR WORSENED JOINT PAIN", "OCCNWJP"),
    ("OCCUR", "NEW OR WORSENED MUSCLE PAIN", "OCCNWMP"),
    ("SEV", "SWELLING", "SEVSWEL"),
    ("SEV", "PAIN AT INJECTION SITE", "SEVPIS"),
    ("SEV", "FATIGUE", "SEVFAT"),
    ("SEV", "HEADACHE", "SEVHEAD"),
    ("SEV", "DIARRHEA", "SEVDIAR"),
    ("SEV", "NEW OR WORSENED JOINT PAIN", "SEVNWJP"),
    ("SEV", "NEW OR WORSENED MUSCLE PAIN", "SEVNWMP"),
    ("MAXDIAM", "SWELLING", "MDISW"),
    ("MAXSEV", "PAIN AT INJECTION SITE", "MAXSPIS"),
    ("OCCUR", "VOMITING", "OCCVOM"),
    ("DIAMETER", "SWELLING", "DIASWEL"),
];

/// Lookup additions for reactions the template's lookup lacks (PARAMN 33 on).
const MORE_PARAMS: [(&str, &str, &str); 5] = [
    ("SEV", "CHILLS", "SEVCHILL"),
    ("SEV", "VOMITING", "SEVVOM"),
    ("MAXSEV", "CHILLS", "MAXSCHIL"),
    ("MAXSEV", "DIARRHEA", "MAXSDIAR"),
    ("MAXSEV", "VOMITING", "MAXSVOM"),
];

fn lookup(testcd: &str, obj: &str) -> Option<(String, f64)> {
    PARAMS
        .iter()
        .chain(&MORE_PARAMS)
        .position(|(t, o, _)| *t == testcd && *o == obj)
        .map(|k| {
            let cd = PARAMS.iter().chain(&MORE_PARAMS).nth(k).expect("found").2;
            (cd.to_string(), k as f64 + 1.0)
        })
}

/// stringr `str_to_sentence()` for ASCII text.
fn sentence(s: &str) -> String {
    let lower = s.to_ascii_lowercase();
    let mut c = lower.chars();
    match c.next() {
        Some(f) => f.to_ascii_uppercase().to_string() + c.as_str(),
        None => String::new(),
    }
}

#[derive(Clone, Debug)]
struct FaRec {
    sl: usize,
    testcd: String,
    test: String,
    obj: String,
    group: Option<String>,
    eval: Option<String>,
    orres: Option<String>,
    avalc: Option<String>,
    aval: Option<f64>,
    temp: Option<f64>,
    site: bool,
    loc: Option<String>,
    lat: Option<String>,
    tptref: Option<String>,
    tpt: Option<String>,
    tptnum: Option<f64>,
    adt: Option<Date>,
    extrt: Option<String>,
    exdose: Option<f64>,
    paramcd: Option<String>,
    paramn: Option<f64>,
    anl02: bool,
    anl03: bool,
    eventfl: Option<&'static str>,
    eventdfl: Option<&'static str>,
}

/// AVAL from FASTRESN, or the severity as a number.
fn sev_to_numeric(avalc: Option<&str>, aval: Option<f64>) -> Option<f64> {
    match avalc {
        Some("NONE") => Some(0.0),
        Some("MILD") => Some(1.0),
        Some("MODERATE") => Some(2.0),
        Some("SEVERE") => Some(3.0),
        _ => aval,
    }
}

/// Derives ADFACE.
pub fn adface(
    face: &Dataset,
    ex: &Dataset,
    vs: Option<&Dataset>,
    adsl: &Dataset,
) -> Result<Dataset, DeriveError> {
    let f = Frame::new(face);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let trtsdt = sl.dates("TRTSDT")?;
    let periods = Periods::new(adsl)?;

    // EX variables by (subject, link group, location, laterality).
    let e = Frame::new(ex);
    type ExVars = (Option<String>, Option<String>, Option<String>, Option<f64>);
    let mut exv: HashMap<(String, String), Vec<ExVars>> = HashMap::new();
    if e.has("EXLNKGRP") {
        let (eu, eg, el, ela, et, ed) = (
            e.chars("USUBJID")?,
            e.chars("EXLNKGRP")?,
            e.chars("EXLOC")?,
            e.chars("EXLAT")?,
            e.chars("EXTRT")?,
            e.nums("EXDOSE")?,
        );
        for i in 0..e.n() {
            if let (Some(x), Some(g)) = (eu[i].clone(), eg[i].clone()) {
                exv.entry((x, g)).or_default().push((
                    el[i].clone(),
                    ela[i].clone(),
                    et[i].clone(),
                    ed[i],
                ));
            }
        }
    }

    let (
        u,
        testcd,
        test,
        obj,
        fcat,
        fscat,
        grp,
        eval,
        orres,
        stresc,
        stresn,
        loc,
        lat,
        dtc,
        tpt,
        tptnum,
        tptref,
    ) = (
        f.chars("USUBJID")?,
        f.chars("FATESTCD")?,
        f.chars("FATEST")?,
        f.chars("FAOBJ")?,
        f.chars("FACAT")?,
        f.chars("FASCAT")?,
        f.chars("FALNKGRP")?,
        f.chars("FAEVAL")?,
        f.chars("FAORRES")?,
        f.chars("FASTRESC")?,
        f.nums("FASTRESN")?,
        f.chars("FALOC")?,
        f.chars("FALAT")?,
        f.chars("FADTC")?,
        f.chars("FATPT")?,
        f.nums("FATPTNUM")?,
        f.chars("FATPTREF")?,
    );
    let mut base: Vec<FaRec> = Vec::new();
    // derive_vars_merged_vaccine(): administration-site records, then systemic ones.
    for want_site in [true, false] {
        for i in 0..f.n() {
            let site = fscat[i].as_deref() == Some("ADMINISTRATION SITE");
            let systemic = fscat[i].as_deref() == Some("SYSTEMIC");
            if fcat[i].as_deref() != Some(REACTOGENICITY)
                || !(site || systemic)
                || site != want_site
            {
                continue;
            }
            let Some(&r) = u[i].as_deref().and_then(|x| sl_idx.get(x)) else {
                continue;
            };
            let key = (
                u[i].clone().unwrap_or_default(),
                tptref[i].clone().unwrap_or_default(),
            );
            let ex_match = exv.get(&key).and_then(|v| {
                v.iter()
                    .find(|x| !site || (x.0 == loc[i] && x.1 == lat[i]))
                    .map(|x| (x.2.clone(), x.3))
            });
            let avalc = stresc[i].clone();
            base.push(FaRec {
                sl: r,
                testcd: testcd[i].clone().unwrap_or_default(),
                test: test[i].clone().unwrap_or_default(),
                obj: obj[i].clone().unwrap_or_default().to_ascii_uppercase(),
                group: grp[i].clone(),
                eval: eval[i].clone(),
                orres: orres[i].clone(),
                aval: sev_to_numeric(avalc.as_deref(), stresn[i]),
                avalc,
                temp: None,
                site,
                loc: loc[i].clone(),
                lat: lat[i].clone(),
                tptref: tptref[i].clone(),
                tpt: tpt[i].clone(),
                tptnum: tptnum[i],
                adt: dt(&dtc[i]),
                extrt: ex_match.as_ref().and_then(|x| x.0.clone()),
                exdose: ex_match.and_then(|x| x.1),
                paramcd: None,
                paramn: None,
                anl02: false,
                anl03: false,
                eventfl: None,
                eventdfl: None,
            });
        }
    }
    // derive_fever_records(): an occurrence record per diary temperature.
    if let Some(vs) = vs.map(Frame::new)
        && vs.has("VSCAT")
        && !base.iter().any(|r| r.obj == "FEVER")
    {
        let (vu, vt, vc, vv, vun, veval, vdtc, vtpt, vtptn, vref) = (
            vs.chars("USUBJID")?,
            vs.chars("VSTESTCD")?,
            vs.chars("VSCAT")?,
            vs.nums("VSSTRESN")?,
            vs.chars("VSSTRESU")?,
            vs.chars("VSEVAL")?,
            vs.chars("VSDTC")?,
            vs.chars("VSTPT")?,
            vs.nums("VSTPTNUM")?,
            vs.chars("VSTPTREF")?,
        );
        for i in 0..vs.n() {
            if vc[i].as_deref() != Some(REACTOGENICITY) || vt[i].as_deref() != Some("TEMP") {
                continue;
            }
            let Some(&r) = vu[i].as_deref().and_then(|x| sl_idx.get(x)) else {
                continue;
            };
            let res = match (vun[i].as_deref(), vv[i]) {
                (Some("C"), Some(t)) => Some(if t >= 38.0 { "Y" } else { "N" }.to_string()),
                _ => None,
            };
            base.push(FaRec {
                sl: r,
                testcd: "OCCUR".into(),
                test: "Occurrence Indicator".into(),
                obj: "FEVER".into(),
                group: None,
                eval: veval[i].clone(),
                orres: res.clone(),
                aval: sev_to_numeric(res.as_deref(), None),
                avalc: res,
                temp: vv[i],
                site: false,
                loc: None,
                lat: None,
                tptref: vref[i].clone(),
                tpt: vtpt[i].clone(),
                tptnum: vtptn[i],
                adt: dt(&vdtc[i]),
                extrt: None,
                exdose: None,
                paramcd: None,
                paramn: None,
                anl02: false,
                anl03: false,
                eventfl: None,
                eventdfl: None,
            });
        }
    }
    // derive_diam_to_sev_records(): severities of redness and swelling, placed first.
    let mut recs: Vec<FaRec> = base
        .iter()
        .filter(|r| r.testcd == "DIAMETER" && matches!(r.obj.as_str(), "REDNESS" | "SWELLING"))
        .map(|r| {
            let a = r.aval;
            let avalc = match a {
                Some(v) if (0.0..=2.0).contains(&v) => Some("NONE".to_string()),
                Some(v) if v > 2.0 && v <= 5.0 => Some("MILD".to_string()),
                Some(v) if v > 5.0 && v <= 10.0 => Some("MODERATE".to_string()),
                Some(v) if v > 10.0 => Some("SEVERE".to_string()),
                _ => r.avalc.clone(),
            };
            FaRec {
                testcd: "SEV".into(),
                test: "Severity/Intensity".into(),
                aval: match avalc.as_deref() {
                    Some("NONE") => Some(0.0),
                    Some("MILD") => Some(1.0),
                    Some("MODERATE") => Some(2.0),
                    Some("SEVERE") => Some(3.0),
                    _ => None,
                },
                avalc,
                ..r.clone()
            }
        })
        .collect();
    recs.extend(base);
    // derive_extreme_records(): the latest record with the largest value per group.
    type Group = (usize, String, Option<String>);
    let extreme = |recs: &[FaRec],
                   keep: &dyn Fn(&FaRec) -> bool,
                   group: &dyn Fn(&FaRec) -> Group,
                   value: &dyn Fn(&FaRec) -> Option<f64>,
                   code: &str,
                   name: &str|
     -> Vec<FaRec> {
        let mut best: BTreeMap<(usize, String, Option<String>), usize> = BTreeMap::new();
        for (k, r) in recs.iter().enumerate().filter(|(_, r)| keep(r)) {
            let g = group(r);
            let better = match best.get(&g) {
                None => true,
                // arrange() puts missing values last, so a missing value wins "last".
                Some(&b) => match (value(&recs[b]), value(r)) {
                    (_, None) => true,
                    (None, Some(_)) => false,
                    (Some(x), Some(v)) => v >= x,
                },
            };
            if better {
                best.insert(g, k);
            }
        }
        best.into_values()
            .map(|k| FaRec {
                testcd: code.into(),
                test: name.into(),
                ..recs[k].clone()
            })
            .collect()
    };
    let by_ref = |r: &FaRec| (r.sl, r.obj.clone(), r.tptref.clone());
    let mut extra = extreme(
        &recs,
        &|r| r.testcd == "SEV",
        &by_ref,
        &|r| r.aval,
        "MAXSEV",
        "Maximum Severity",
    );
    extra.extend(extreme(
        &recs,
        &|r| r.testcd == "DIAMETER" && matches!(r.obj.as_str(), "REDNESS" | "SWELLING"),
        &|r| (r.sl, r.obj.clone(), r.group.clone()),
        &|r| r.aval,
        "MAXDIAM",
        "Maximum Diameter",
    ));
    extra.extend(extreme(
        &recs,
        &|r| r.obj == "FEVER",
        &by_ref,
        &|r| r.temp,
        "MAXTEMP",
        "Maximum Temperature",
    ));
    recs.extend(extra);
    // derive_vars_params().
    for r in recs.iter_mut() {
        if let Some((cd, n)) = lookup(&r.testcd, &r.obj) {
            r.paramcd = Some(cd);
            r.paramn = Some(n);
        }
    }
    // derive_vars_max_flag(): the first of the largest positive values per group.
    let flag = |recs: &mut Vec<FaRec>, by_ref: bool, set: fn(&mut FaRec)| {
        type FlagGroup = (usize, String, Option<String>, Option<String>);
        let mut groups: BTreeMap<FlagGroup, Vec<usize>> = BTreeMap::new();
        for (k, r) in recs.iter().enumerate().filter(|(_, r)| r.aval.is_some()) {
            let g = (
                r.sl,
                r.obj.clone(),
                if by_ref { r.tptref.clone() } else { None },
                r.paramcd.clone(),
            );
            groups.entry(g).or_default().push(k);
        }
        for idx in groups.values() {
            let first = idx.iter().copied().min_by(|&a, &b| {
                let (x, y) = (&recs[a], &recs[b]);
                y.aval
                    .partial_cmp(&x.aval)
                    .expect("finite")
                    .then_with(|| x.tpt.cmp(&y.tpt))
                    .then(a.cmp(&b))
            });
            if let Some(k) = first.filter(|&k| recs[k].aval.is_some_and(|v| v > 0.0)) {
                set(&mut recs[k]);
            }
        }
    };
    flag(&mut recs, true, |r| r.anl02 = true);
    flag(&mut recs, false, |r| r.anl03 = true);
    // derive_vars_event_flag(): with R's missing-value logic.
    let event = |r: &FaRec| -> Option<bool> {
        let named = r
            .avalc
            .as_deref()
            .is_some_and(|c| ["Y", "MILD", "MODERATE", "SEVERE"].contains(&c));
        match r.aval {
            Some(v) => Some(v > 2.5 || named),
            None => named.then_some(true),
        }
    };
    let no_event = |r: &FaRec| -> Option<bool> {
        let n = r.avalc.as_deref() == Some("N");
        match r.aval {
            Some(v) => Some(v <= 2.5 || n),
            None => n.then_some(true),
        }
    };
    let mut groups: BTreeMap<(usize, String, Option<String>), Vec<usize>> = BTreeMap::new();
    for (k, r) in recs.iter().enumerate() {
        groups
            .entry((r.sl, r.obj.clone(), r.tptref.clone()))
            .or_default()
            .push(k);
    }
    for idx in groups.values() {
        let any = idx.iter().any(|&k| event(&recs[k]) == Some(true));
        for &k in idx {
            let r = &recs[k];
            let row = |r: &FaRec| {
                if event(r) == Some(true) {
                    Some("Y")
                } else if no_event(r) == Some(true) {
                    Some("N")
                } else {
                    None
                }
            };
            let fl = if any {
                Some("Y")
            } else if no_event(r) == Some(true) {
                Some("N")
            } else {
                None
            };
            let dfl = if r.testcd.contains("MAX") {
                None
            } else {
                row(r)
            };
            recs[k].eventfl = fl;
            recs[k].eventdfl = dfl;
        }
    }
    // Output order: subject, parameter, time point reference, time point.
    recs.sort_by(|a, b| {
        (
            a.sl,
            a.paramn.map(|p| p as i64),
            &a.tptref,
            a.tptnum.map(|t| t as i64),
        )
            .cmp(&(
                b.sl,
                b.paramn.map(|p| p as i64),
                &b.tptref,
                b.tptnum.map(|t| t as i64),
            ))
    });
    // post_process_reacto(): derived records lose their FA variables.
    let blank = |r: &FaRec| {
        matches!(r.testcd.as_str(), "MAXDIAM" | "MAXSEV" | "MAXTEMP")
            || (matches!(r.testcd.as_str(), "OCCUR" | "SEV")
                && matches!(r.obj.as_str(), "FEVER" | "REDNESS" | "SWELLING"))
    };
    let n = recs.len();
    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let fa = |k: usize, v: Option<String>| if blank(&recs[k]) { None } else { v };
    let period = |r: &FaRec| periods.of(r.sl, r.adt);
    let param = |r: &FaRec| {
        let parts: Vec<&str> = [
            Some(r.obj.as_str()),
            Some(r.test.as_str()),
            r.loc.as_deref(),
            r.lat.as_deref(),
        ]
        .into_iter()
        .flatten()
        .collect();
        sentence(&parts.join(" "))
    };
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_VARS)?;
    cols.extend([
        (
            "TRTP",
            text(n, |k| period(&recs[k]).and_then(|p| p.3.clone())),
        ),
        (
            "TRTA",
            text(n, |k| period(&recs[k]).and_then(|p| p.5.clone())),
        ),
        (
            "APERIOD",
            Values::Int((0..n).map(|k| period(&recs[k]).map(|p| p.0)).collect()),
        ),
        (
            "APERSDT",
            Values::Date((0..n).map(|k| period(&recs[k]).map(|p| p.1)).collect()),
        ),
        (
            "APEREDT",
            Values::Date((0..n).map(|k| period(&recs[k]).map(|p| p.2)).collect()),
        ),
        ("FATESTCD", text(n, |k| fa(k, Some(recs[k].testcd.clone())))),
        ("FATEST", text(n, |k| fa(k, Some(recs[k].test.clone())))),
        ("FAOBJ", text(n, |k| fa(k, Some(recs[k].obj.clone())))),
        ("FALNKGRP", text(n, |k| fa(k, recs[k].group.clone()))),
        ("FAEVAL", text(n, |k| fa(k, recs[k].eval.clone()))),
        ("FAORRES", text(n, |k| fa(k, recs[k].orres.clone()))),
        ("PARAMCD", text(n, |k| recs[k].paramcd.clone())),
        ("PARAM", text(n, |k| Some(param(&recs[k])))),
        ("PARAMN", num(n, |k| recs[k].paramn)),
        ("PARCAT1", text(n, |_| Some(REACTOGENICITY.to_string()))),
        (
            "PARCAT2",
            text(n, |k| {
                Some(
                    if recs[k].site {
                        "ADMINISTRATION SITE"
                    } else {
                        "SYSTEMIC"
                    }
                    .to_string(),
                )
            }),
        ),
        ("AVAL", num(n, |k| recs[k].aval)),
        ("AVALC", text(n, |k| recs[k].avalc.clone())),
        ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
        (
            "ADY",
            Values::Int(recs.iter().map(|r| rel_day(r.adt, trtsdt[r.sl])).collect()),
        ),
        ("ATPT", text(n, |k| recs[k].tpt.clone())),
        ("ATPTN", num(n, |k| recs[k].tptnum)),
        ("ATPTREF", text(n, |k| recs[k].tptref.clone())),
        ("EXTRT", text(n, |k| recs[k].extrt.clone())),
        ("EXDOSE", num(n, |k| recs[k].exdose)),
        ("ANL01FL", text(n, |_| Some("Y".into()))),
        ("ANL02FL", text(n, |k| y(recs[k].anl02))),
        ("ANL03FL", text(n, |k| y(recs[k].anl03))),
        ("EVENTFL", text(n, |k| recs[k].eventfl.map(str::to_string))),
        (
            "EVENTDFL",
            text(n, |k| recs[k].eventdfl.map(str::to_string)),
        ),
    ]);
    Ok(build("ADFACE", cols)?)
}
