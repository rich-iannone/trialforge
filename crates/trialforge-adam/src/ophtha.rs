//! ADBCVA (best corrected visual acuity analysis) from OE and ADSL, following admiralophtha's
//! `ad_adbcva.R` template:
//!
//! - `AFEYE` (admiralophtha `derive_var_afeye()`): `Study Eye` when `OELAT` is the subject's
//!   `STUDYEYE`, otherwise `Fellow Eye` (location `EYE` only).
//! - Letters parameters `SBCVA`/`FBCVA` from OE, and logMAR parameters `SBCVALOG`/`FBCVALOG`
//!   computed from them (`-0.02 × letters + 1.7`, admiralophtha `convert_etdrs_to_logmar()`).
//!   Adaptation: logMAR is rounded to 2 decimals, the scale's precision. Unrounded, values such
//!   as `-0.0600000000000001` make `AVALC` (R `as.character()`) depend on the platform's long
//!   double support.
//! - `AVISIT`: `Screening` for screening visits, otherwise the title-cased `VISIT`;
//!   `AVISITN` = `VISITNUM`.
//! - `ONTRTFL` (admiral `derive_var_ontrtfl()`, excluding the Baseline visit): `ADT` after
//!   `TRTSDT` and no more than 112 days after `TRTEDT`. Adaptation: the template uses no window
//!   after the last dose, which would drop every assessment after the last injection; 112 days
//!   is the longest maintenance interval.
//! - `ABLFL`: the last record on or before `TRTSDT` per parameter; `BASE`, `BASEC`, `CHG`, `PCHG`
//!   (admiral `derive_var_base()`, `derive_var_chg()`, `derive_var_pchg()`).
//! - `ANL01FL`: the last record (by `ADT`, `AVAL`) per parameter and analysis visit; `ANL02FL`:
//!   the last per parameter and `ABLFL`; `WORS01FL`: the worst (lowest, then earliest) letters
//!   score per parameter; all among baseline and on-treatment records.
//! - `CRIT1`–`CRIT8` (admiral `derive_vars_crit_flag()`) and `AVALCAT1`/`AVALCA1N` (Snellen
//!   equivalents) on the letters parameters.

use std::collections::{BTreeMap, HashMap};

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::{DeriveError, Frame, dt, index, rel_day};

/// `OETESTCD` of best corrected visual acuity.
const BCVA: &str = "NUMLCOR";
/// Days after `TRTEDT` that still count as on treatment.
const ON_TREATMENT_WINDOW: i32 = 112;

/// A response criterion on `CHG`: (description, condition).
pub(crate) type Criterion = (&'static str, fn(f64) -> bool);

/// The template's response criteria on `CHG`.
pub(crate) const CRITERIA: [Criterion; 8] = [
    ("0 <= CHG <= 5", |c| (0.0..=5.0).contains(&c)),
    ("-5 <= CHG <= -1", |c| (-5.0..=-1.0).contains(&c)),
    ("10 <= CHG <= 15", |c| (10.0..=15.0).contains(&c)),
    ("CHG <= -20", |c| c <= -20.0),
    ("CHG <= 5", |c| c <= 5.0),
    ("CHG <= 10", |c| c <= 10.0),
    ("CHG >= -15", |c| c >= -15.0),
    ("CHG >= 15", |c| c >= 15.0),
];

/// The Snellen equivalent of an ETDRS letters score: (`AVALCA1N`, `AVALCAT1`).
pub(crate) fn snellen(letters: f64) -> Option<(f64, &'static str)> {
    const BANDS: [(f64, f64, f64, &str); 21] = [
        (0.0, 3.0, 1000.0, "< 20/800"),
        (4.0, 8.0, 800.0, "20/800"),
        (9.0, 13.0, 640.0, "20/640"),
        (14.0, 18.0, 500.0, "20/500"),
        (19.0, 23.0, 400.0, "20/400"),
        (24.0, 28.0, 320.0, "20/320"),
        (29.0, 33.0, 250.0, "20/250"),
        (34.0, 38.0, 200.0, "20/200"),
        (39.0, 43.0, 160.0, "20/160"),
        (44.0, 48.0, 125.0, "20/125"),
        (49.0, 53.0, 100.0, "20/100"),
        (54.0, 58.0, 80.0, "20/80"),
        (59.0, 63.0, 63.0, "20/63"),
        (64.0, 68.0, 50.0, "20/50"),
        (69.0, 73.0, 40.0, "20/40"),
        (74.0, 78.0, 32.0, "20/32"),
        (79.0, 83.0, 25.0, "20/25"),
        (84.0, 88.0, 20.0, "20/20"),
        (89.0, 93.0, 16.0, "20/16"),
        (94.0, 97.0, 12.0, "20/12"),
        (98.0, f64::INFINITY, 1.0, "> 20/12"),
    ];
    BANDS
        .iter()
        .find(|(lo, hi, _, _)| letters >= *lo && letters <= *hi)
        .map(|(_, _, n, c)| (*n, *c))
}

/// R's `as.character()` of a double: 15 significant digits, trailing zeros dropped.
pub(crate) fn r_character(x: f64) -> String {
    if x == 0.0 {
        return "0".into();
    }
    let sci = format!("{x:.14e}");
    let (mantissa, exp) = sci.split_once('e').expect("scientific format");
    let exp: i32 = exp.parse().expect("exponent");
    let negative = mantissa.starts_with('-');
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    let digits = digits.trim_end_matches('0');
    let digits = if digits.is_empty() { "0" } else { digits };
    let body = if !(-4..15).contains(&exp) {
        // R switches to scientific notation; not reached for visual acuity values.
        let (first, rest) = digits.split_at(1);
        let rest = if rest.is_empty() {
            String::new()
        } else {
            format!(".{rest}")
        };
        format!(
            "{first}{rest}e{}{:02}",
            if exp < 0 { "-" } else { "+" },
            exp.abs()
        )
    } else if exp < 0 {
        format!("0.{}{digits}", "0".repeat((-exp - 1) as usize))
    } else {
        let int_len = exp as usize + 1;
        if digits.len() <= int_len {
            format!("{digits}{}", "0".repeat(int_len - digits.len()))
        } else {
            format!("{}.{}", &digits[..int_len], &digits[int_len..])
        }
    };
    if negative { format!("-{body}") } else { body }
}

/// `stringr::str_to_title()` for ASCII visit names.
pub(crate) fn title_case(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut start = true;
    for c in s.chars() {
        if c.is_ascii_alphanumeric() {
            out.push(if start {
                c.to_ascii_uppercase()
            } else {
                c.to_ascii_lowercase()
            });
            start = false;
        } else {
            out.push(c);
            start = true;
        }
    }
    out
}

#[derive(Clone, Debug)]
struct Rec {
    sl: usize,
    paramcd: &'static str,
    aval: f64,
    adt: Option<Date>,
    visit: Option<String>,
    visitnum: Option<i64>,
    seq: Option<i64>,
    afeye: &'static str,
    ontrt: bool,
    ablfl: bool,
    anl01: bool,
    anl02: bool,
    worst: bool,
}

fn param(paramcd: &str) -> (&'static str, f64, &'static str) {
    match paramcd {
        "SBCVA" => ("Study Eye Visual Acuity Score (letters)", 1.0, "letters"),
        "FBCVA" => ("Fellow Eye Visual Acuity Score (letters)", 2.0, "letters"),
        "SBCVALOG" => ("Study Eye Visual Acuity LogMAR Score", 3.0, "LogMAR"),
        _ => ("Fellow Eye Visual Acuity LogMAR Score", 4.0, "LogMAR"),
    }
}

/// Derives ADBCVA.
pub fn adbcva(oe: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let oe = Frame::new(oe);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let (trtsdt, trtedt, studyeye) = (
        sl.dates("TRTSDT")?,
        sl.dates("TRTEDT")?,
        sl.chars("STUDYEYE")?,
    );
    let (u, testcd, stresn, dtc, seq, loc, lat, visit, visitnum) = (
        oe.chars("USUBJID")?,
        oe.chars("OETESTCD")?,
        oe.nums("OESTRESN")?,
        oe.chars("OEDTC")?,
        oe.ints("OESEQ")?,
        oe.chars("OELOC")?,
        oe.chars("OELAT")?,
        oe.chars("VISIT")?,
        oe.ints("VISITNUM")?,
    );
    let mut recs: Vec<Rec> = Vec::new();
    for i in 0..oe.n() {
        if testcd[i].as_deref() != Some(BCVA) {
            continue;
        }
        let (Some(x), Some(v)) = (u[i].as_deref(), stresn[i]) else {
            continue;
        };
        let Some(&r) = sl_idx.get(x) else { continue };
        let eyes = ["LEFT", "RIGHT", "BILATERAL"];
        let (Some(l), Some(se)) = (lat[i].as_deref(), studyeye[r].as_deref()) else {
            continue;
        };
        if loc[i].as_deref() != Some("EYE") || !eyes.contains(&l) || !eyes.contains(&se) {
            continue;
        }
        let study = l == se;
        let base = Rec {
            sl: r,
            paramcd: if study { "SBCVA" } else { "FBCVA" },
            aval: v,
            adt: dt(&dtc[i]),
            visit: visit[i].clone(),
            visitnum: visitnum[i],
            seq: seq[i],
            afeye: if study { "Study Eye" } else { "Fellow Eye" },
            ontrt: false,
            ablfl: false,
            anl01: false,
            anl02: false,
            worst: false,
        };
        let log = Rec {
            paramcd: if study { "SBCVALOG" } else { "FBCVALOG" },
            aval: trialforge_core::math::round_dp(-0.02 * v + 1.7, 2),
            seq: None,
            ..base.clone()
        };
        recs.push(base);
        recs.push(log);
    }

    let avisit = |r: &Rec| {
        r.visit.as_deref().map(|v| {
            if v.contains("SCREEN") {
                "Screening".to_string()
            } else {
                title_case(v)
            }
        })
    };
    // ONTRTFL.
    for r in recs.iter_mut() {
        let (Some(start), Some(adt)) = (trtsdt[r.sl], r.adt) else {
            continue;
        };
        let baseline_visit = avisit(r).as_deref() == Some("Baseline");
        let end_ok = trtedt[r.sl].is_none_or(|e| adt <= e.add_days(ON_TREATMENT_WINDOW));
        r.ontrt = (adt == start && !baseline_visit) || (start < adt && end_ok);
    }
    // Groups per subject and parameter.
    let mut groups: BTreeMap<(usize, &str), Vec<usize>> = BTreeMap::new();
    for (k, r) in recs.iter().enumerate() {
        groups.entry((r.sl, r.paramcd)).or_default().push(k);
    }
    for idx in groups.values() {
        let b = idx
            .iter()
            .copied()
            .filter(|&k| {
                recs[k]
                    .adt
                    .is_some_and(|d| trtsdt[recs[k].sl].is_some_and(|s| d <= s))
            })
            .max_by_key(|&k| {
                (
                    recs[k].adt,
                    recs[k].visitnum,
                    recs[k].seq.unwrap_or(i64::MAX),
                )
            });
        if let Some(b) = b {
            recs[b].ablfl = true;
        }
        let eligible: Vec<usize> = idx
            .iter()
            .copied()
            .filter(|&k| recs[k].visitnum.is_some() && (recs[k].ontrt || recs[k].ablfl))
            .collect();
        // ANL01FL: last by (ADT, AVAL) per analysis visit.
        let mut by_visit: BTreeMap<Option<String>, usize> = BTreeMap::new();
        for &k in &eligible {
            let key = |k: usize| (recs[k].adt, recs[k].aval);
            match by_visit.get(&avisit(&recs[k])) {
                Some(&b) if key(b) > key(k) => {}
                _ => {
                    by_visit.insert(avisit(&recs[k]), k);
                }
            }
        }
        for &k in by_visit.values() {
            recs[k].anl01 = true;
        }
        // ANL02FL: last by ADT per ABLFL.
        for flag in [true, false] {
            if let Some(k) = eligible
                .iter()
                .copied()
                .filter(|&k| recs[k].ablfl == flag)
                .max_by(|&a, &b| recs[a].adt.cmp(&recs[b].adt).then(a.cmp(&b)))
            {
                recs[k].anl02 = true;
            }
        }
        // WORS01FL: first by (AVAL, ADT), letters only.
        if let Some(k) = eligible
            .iter()
            .copied()
            .filter(|&k| matches!(recs[k].paramcd, "SBCVA" | "FBCVA"))
            .min_by(|&a, &b| {
                recs[a]
                    .aval
                    .total_cmp(&recs[b].aval)
                    .then(recs[a].adt.cmp(&recs[b].adt))
                    .then(a.cmp(&b))
            })
        {
            recs[k].worst = true;
        }
    }
    let mut bases: HashMap<(usize, &str), f64> = HashMap::new();
    for r in recs.iter().filter(|r| r.ablfl) {
        bases.insert((r.sl, r.paramcd), r.aval);
    }
    // ASEQ order: PARAMCD, ADT, AVISITN, VISITNUM.
    recs.sort_by(|a, b| {
        (a.sl, a.paramcd, a.adt, a.visitnum).cmp(&(b.sl, b.paramcd, b.adt, b.visitnum))
    });
    let mut aseq: Vec<Option<i64>> = Vec::with_capacity(recs.len());
    for (i, r) in recs.iter().enumerate() {
        let k = if i > 0 && recs[i - 1].sl == r.sl {
            aseq[i - 1].unwrap_or(0) + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }

    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let from_sl = |name: &str| -> Result<Values, DeriveError> {
        let col = sl
            .ds
            .column(name)
            .ok_or_else(|| DeriveError::Input(format!("ADSL has no {name}")))?;
        Ok(match &col.values {
            Values::Char(v) => Values::Char(sl_rows.iter().map(|&r| v[r].clone()).collect()),
            Values::Num(v) => Values::Num(sl_rows.iter().map(|&r| v[r]).collect()),
            Values::Int(v) => Values::Int(sl_rows.iter().map(|&r| v[r]).collect()),
            Values::Date(v) => Values::Date(sl_rows.iter().map(|&r| v[r]).collect()),
            Values::DateTime(v) => Values::DateTime(sl_rows.iter().map(|&r| v[r]).collect()),
        })
    };
    let base: Vec<Option<f64>> = recs
        .iter()
        .map(|r| bases.get(&(r.sl, r.paramcd)).copied())
        .collect();
    let chg: Vec<Option<f64>> = recs
        .iter()
        .zip(&base)
        .map(|(r, b)| b.map(|b| r.aval - b))
        .collect();
    let letters = |r: &Rec| matches!(r.paramcd, "SBCVA" | "FBCVA");
    let y = |b: bool| b.then(|| "Y".to_string());
    let text = |f: &dyn Fn(&Rec) -> Option<String>| Values::Char(recs.iter().map(f).collect());
    let mut cols: Vec<(&str, Values)> = vec![
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
        ("TRTSDT", from_sl("TRTSDT")?),
        ("TRTEDT", from_sl("TRTEDT")?),
        ("STUDYEYE", from_sl("STUDYEYE")?),
        ("ASEQ", Values::Int(aseq)),
        ("PARAMCD", text(&|r| Some(r.paramcd.to_string()))),
        ("PARAM", text(&|r| Some(param(r.paramcd).0.to_string()))),
        (
            "PARAMN",
            Values::Num(recs.iter().map(|r| Some(param(r.paramcd).1)).collect()),
        ),
        (
            "AVAL",
            Values::Num(recs.iter().map(|r| Some(r.aval)).collect()),
        ),
        ("AVALC", text(&|r| Some(r_character(r.aval)))),
        ("AVALU", text(&|r| Some(param(r.paramcd).2.to_string()))),
        (
            "AVALCAT1",
            text(&|r| {
                letters(r)
                    .then(|| snellen(r.aval))
                    .flatten()
                    .map(|s| s.1.to_string())
            }),
        ),
        (
            "AVALCA1N",
            Values::Num(
                recs.iter()
                    .map(|r| letters(r).then(|| snellen(r.aval)).flatten().map(|s| s.0))
                    .collect(),
            ),
        ),
        ("BASE", Values::Num(base.clone())),
        (
            "BASEC",
            Values::Char(base.iter().map(|b| b.map(r_character)).collect()),
        ),
        ("BASETYPE", text(&|_| Some("LAST".to_string()))),
        ("CHG", Values::Num(chg.clone())),
        (
            "PCHG",
            Values::Num(
                recs.iter()
                    .zip(&base)
                    .map(|(r, b)| match b {
                        Some(b) if *b != 0.0 => Some((r.aval - b) / b.abs() * 100.0),
                        _ => None,
                    })
                    .collect(),
            ),
        ),
        ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
        (
            "ADY",
            Values::Int(
                recs.iter()
                    .zip(&sl_rows)
                    .map(|(r, &s)| rel_day(r.adt, trtsdt[s]))
                    .collect(),
            ),
        ),
        ("AVISIT", text(&avisit)),
        (
            "AVISITN",
            Values::Num(recs.iter().map(|r| r.visitnum.map(|v| v as f64)).collect()),
        ),
        ("AFEYE", text(&|r| Some(r.afeye.to_string()))),
        ("ONTRTFL", text(&|r| y(r.ontrt))),
        ("ABLFL", text(&|r| y(r.ablfl))),
        ("ANL01FL", text(&|r| y(r.anl01))),
        ("ANL02FL", text(&|r| y(r.anl02))),
        ("WORS01FL", text(&|r| y(r.worst))),
    ];
    const NAMES: [(&str, &str); 8] = [
        ("CRIT1", "CRIT1FL"),
        ("CRIT2", "CRIT2FL"),
        ("CRIT3", "CRIT3FL"),
        ("CRIT4", "CRIT4FL"),
        ("CRIT5", "CRIT5FL"),
        ("CRIT6", "CRIT6FL"),
        ("CRIT7", "CRIT7FL"),
        ("CRIT8", "CRIT8FL"),
    ];
    for ((crit, fl), (desc, cond)) in NAMES.iter().zip(CRITERIA) {
        cols.push((crit, text(&|r| letters(r).then(|| desc.to_string()))));
        cols.push((
            fl,
            Values::Char(
                recs.iter()
                    .zip(&chg)
                    .map(|(r, c)| {
                        if !letters(r) {
                            return None;
                        }
                        c.map(|c| if cond(c) { "Y" } else { "N" }.to_string())
                    })
                    .collect(),
            ),
        ));
    }
    cols.extend([
        ("VISIT", text(&|r| r.visit.clone())),
        (
            "VISITNUM",
            Values::Num(recs.iter().map(|r| r.visitnum.map(|v| v as f64)).collect()),
        ),
        ("SRCDOM", text(&|r| r.seq.map(|_| "OE".to_string()))),
        ("SRCVAR", text(&|r| r.seq.map(|_| "OESTRESN".to_string()))),
        ("SRCSEQ", Values::Int(recs.iter().map(|r| r.seq).collect())),
    ]);
    Ok(build("ADBCVA", cols)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_character_matches_r() {
        assert_eq!(r_character(82.0), "82");
        assert_eq!(r_character(1.7), "1.7");
        assert_eq!(r_character(-0.3), "-0.3");
        assert_eq!(r_character(0.0599999999999998), "0.0599999999999998");
        assert_eq!(r_character(100.0), "100");
    }

    #[test]
    fn snellen_bands_cover_the_scale() {
        assert_eq!(snellen(0.0), Some((1000.0, "< 20/800")));
        assert_eq!(snellen(85.0), Some((20.0, "20/20")));
        assert_eq!(snellen(100.0), Some((1.0, "> 20/12")));
    }
}
