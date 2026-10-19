//! Date defects.

use std::collections::BTreeMap;

use trialforge_core::{Date, Draws, Stream, study_day};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, parse_date, pick, seq_key};

/// The largest step back from AESTDTC, in days.
const MAX_DAYS_BACK: i32 = 14;

/// The largest step past DTHDTC, in days.
const MAX_DAYS_AFTER: i32 = 14;

/// `ae_end_before_start`: AEENDTC moves 1–14 days (subtle: 1 day) before AESTDTC, never before
/// informed consent (so only the order is wrong), and AEENDY is recomputed from it. ADAE keeps
/// the clean dates, so its ASTDT/AENDT (and ADURN) no longer match AE.
pub(crate) fn ae_end_before_start(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let refs = cx.subjects();
    let Some(ae) = cx.sdtm("AE") else {
        return Vec::new();
    };
    let (Some(u), Some(seq), Some(st), Some(en)) = (
        chars(ae, "USUBJID"),
        ints(ae, "AESEQ"),
        chars(ae, "AESTDTC"),
        chars(ae, "AEENDTC"),
    ) else {
        return Vec::new();
    };
    let has_endy = ae.column("AEENDY").is_some();

    // (row, start, days from consent to start).
    let candidates: Vec<(usize, Date, i32)> = (0..ae.n_rows())
        .filter(|&i| !cx.is_used("AE", i))
        .filter_map(|i| {
            let start = parse_date(st[i].as_deref()?)?;
            parse_date(en[i].as_deref()?)?;
            let consent = refs.get(u[i].as_deref()?)?.rfic?;
            let slack = start.days_since(consent);
            (slack >= 1).then_some((i, start, slack))
        })
        .collect();
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );

    let mut edits = Vec::new();
    let mut defects = Vec::new();
    for c in chosen {
        let (i, start, slack) = candidates[c];
        let mut d = Draws::new(stream.label("value").index(i as u64));
        let mut back = d.int_range(1, i64::from(slack.min(MAX_DAYS_BACK))) as i32;
        if cx.subtle {
            back = 1;
        }
        let end = start.add_days(-back);
        let original = en[i].clone().expect("a candidate has an end date");
        // Keep any time of day.
        let injected = format!("{end}{}", &original[10..]);
        let x = u[i].clone().expect("a candidate has a subject");
        edits.push(Edit {
            variable: "AEENDTC".into(),
            row: i,
            value: Cell::Char(Some(injected.clone())),
        });
        if has_endy && let Some(rfst) = refs[x.as_str()].rfst {
            edits.push(Edit {
                variable: "AEENDY".into(),
                row: i,
                value: Cell::Int(Some(i64::from(study_day(end, rfst)))),
            });
        }
        let mut expected = vec![Expected {
            rule: "L1.ORDER",
            dataset: "AE".into(),
            row: Some(i),
        }];
        if let Some(s) = seq[i]
            && let Some(j) = cx.adae_row(&x, s)
        {
            expected.push(Expected {
                rule: "L2.ADAE",
                dataset: "ADAE".into(),
                row: Some(j),
            });
        }
        defects.push(Defect {
            id: "ae_end_before_start",
            dataset: "AE".into(),
            row: i,
            key: seq[i].map(|s| format!("AESEQ={s}")),
            usubjid: Some(x),
            variable: "AEENDTC".into(),
            original: Some(original),
            injected: Some(injected),
            paired_row: None,
            expected,
        });
    }
    cx.edit_sdtm("AE", &edits);
    defects
}

/// DS records whose date anchors a derivation or is the death itself, left alone by
/// `record_after_death`.
const DS_FIXED: &[&str] = &["RANDOMIZED", "INFORMED CONSENT OBTAINED", "DEATH"];

/// `record_after_death`: for a subject who died, AEENDTC or the DSSTDTC of a disposition event
/// other than consent, randomisation or the death moves 1–14 days (subtle: 1 day) past DTHDTC;
/// AEENDY or DSSTDY follows. The record is usually also after RFPENDTC (`L1.PART`), a DS record
/// leaves its EPOCH (`L1.EPOCH`), and ADAE keeps the clean AE dates (`L2.ADAE`).
pub(crate) fn record_after_death(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let subjects = cx.subjects();
    // (dataset, row, date variable).
    let mut candidates: Vec<(&'static str, usize, &'static str)> = Vec::new();
    let died = |u: &Option<String>, d: Option<Date>| -> bool {
        let death = u.as_deref().and_then(|x| subjects.get(x)?.death);
        matches!((d, death), (Some(d), Some(death)) if d <= death)
    };
    let date = |v: &Option<String>| v.as_deref().and_then(parse_date);
    if let Some(ae) = cx.sdtm("AE")
        && let (Some(u), Some(st), Some(en)) = (
            chars(ae, "USUBJID"),
            chars(ae, "AESTDTC"),
            chars(ae, "AEENDTC"),
        )
    {
        for i in (0..ae.n_rows()).filter(|&i| !cx.is_used("AE", i)) {
            if date(&en[i]).is_some() && died(&u[i], date(&st[i])) {
                candidates.push(("AE", i, "AEENDTC"));
            }
        }
    }
    if let Some(ds) = cx.sdtm("DS")
        && let (Some(u), Some(dec), Some(st)) = (
            chars(ds, "USUBJID"),
            chars(ds, "DSDECOD"),
            chars(ds, "DSSTDTC"),
        )
    {
        for i in (0..ds.n_rows()).filter(|&i| !cx.is_used("DS", i)) {
            if dec[i].as_deref().is_some_and(|d| !DS_FIXED.contains(&d))
                && died(&u[i], date(&st[i]))
            {
                candidates.push(("DS", i, "DSSTDTC"));
            }
        }
    }
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );

    let mut defects = Vec::new();
    let mut edits: BTreeMap<&str, Vec<Edit>> = BTreeMap::new();
    for c in chosen {
        let (name, i, var) = candidates[c];
        let ds = cx.sdtm(name).expect("a candidate's dataset exists");
        let x = chars(ds, "USUBJID").expect("candidate")[i]
            .clone()
            .expect("a candidate has a subject");
        let s = subjects[&x];
        let death = s.death.expect("a candidate's subject died");
        let mut d = Draws::new(stream.label("value").label(name).index(i as u64));
        let mut after = d.int_range(1, i64::from(MAX_DAYS_AFTER)) as i32;
        if cx.subtle {
            after = 1;
        }
        let moved = death.add_days(after);
        let original = chars(ds, var).expect("candidate")[i]
            .clone()
            .expect("a candidate has a date");
        let injected = format!("{moved}{}", &original[10..]);
        let e = edits.entry(name).or_default();
        e.push(Edit {
            variable: var.into(),
            row: i,
            value: Cell::Char(Some(injected.clone())),
        });
        let dy = var.replace("DTC", "DY");
        if ints(ds, &dy).is_some()
            && let Some(rfst) = s.rfst
        {
            e.push(Edit {
                variable: dy,
                row: i,
                value: Cell::Int(Some(i64::from(study_day(moved, rfst)))),
            });
        }
        let at = |rule| Expected {
            rule,
            dataset: name.into(),
            row: Some(i),
        };
        let mut expected = vec![at("L1.DEATH")];
        if s.rfpen.is_some_and(|end| moved > end) {
            expected.push(at("L1.PART"));
        }
        if name == "DS" && cx.outside_epoch(name, i, &x, moved) {
            expected.push(at("L1.EPOCH"));
        }
        if name == "AE"
            && let Some(q) = ints(ds, "AESEQ").and_then(|q| q[i])
            && let Some(j) = cx.adae_row(&x, q)
        {
            expected.push(Expected {
                rule: "L2.ADAE",
                dataset: "ADAE".into(),
                row: Some(j),
            });
        }
        defects.push(Defect {
            id: "record_after_death",
            dataset: name.into(),
            row: i,
            usubjid: Some(x),
            key: seq_key(ds, i),
            variable: var.into(),
            original: Some(original),
            injected: Some(injected),
            paired_row: None,
            expected,
        });
    }
    for (name, e) in &edits {
        cx.edit_sdtm(name, e);
    }
    defects
}

/// A study-day variable's name, values and the values of its date.
type DayPair<'a> = (String, &'a [Option<i64>], &'a [Option<String>]);

/// `bad_study_day`: a study day computed as if there were a day 0, one less on or after
/// RFSTDTC and one more before it (so day 1 and day -1 both become day 0, which `L0.DY` also
/// flags). One --DY per record, among those `L1.DY` checks.
pub(crate) fn bad_study_day(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let subjects = cx.subjects();
    // (dataset, row, the study-day variables holding the day of their date).
    let mut candidates: Vec<(&str, usize, Vec<String>)> = Vec::new();
    for ds in cx.sdtm.iter() {
        let p = ds.name.as_str();
        let Some(u) = chars(ds, "USUBJID") else {
            continue;
        };
        let pairs: Vec<DayPair<'_>> = [
            (format!("{p}STDY"), format!("{p}STDTC")),
            (format!("{p}ENDY"), format!("{p}ENDTC")),
            (format!("{p}DY"), format!("{p}DTC")),
        ]
        .into_iter()
        .filter_map(|(dy, dtc)| {
            let (days, dates) = (ints(ds, &dy)?, chars(ds, &dtc)?);
            Some((dy, days, dates))
        })
        .collect();
        if pairs.is_empty() {
            continue;
        }
        for i in 0..ds.n_rows() {
            if cx.used.contains(&(p.to_string(), i)) {
                continue;
            }
            let Some(rfst) = u[i].as_deref().and_then(|x| subjects.get(x)?.rfst) else {
                continue;
            };
            let vars: Vec<String> = pairs
                .iter()
                .filter(|(_, days, dates)| {
                    let date = dates[i].as_deref().and_then(parse_date);
                    days[i].is_some() && days[i] == date.map(|d| i64::from(study_day(d, rfst)))
                })
                .map(|(dy, _, _)| dy.clone())
                .collect();
            if !vars.is_empty() {
                candidates.push((p, i, vars));
            }
        }
    }
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );

    let mut defects = Vec::new();
    let mut edits: BTreeMap<String, Vec<Edit>> = BTreeMap::new();
    for c in chosen {
        let (name, i, vars) = &candidates[c];
        let (name, i) = (*name, *i);
        let ds = cx.sdtm(name).expect("a candidate's dataset exists");
        let mut d = Draws::new(stream.label("value").label(name).index(i as u64));
        let var = &vars[d.below(vars.len() as u64) as usize];
        let day = ints(ds, var).expect("candidate")[i].expect("a candidate has a day");
        let wrong = if day > 0 { day - 1 } else { day + 1 };
        let at = |rule| Expected {
            rule,
            dataset: name.into(),
            row: Some(i),
        };
        let mut expected = vec![at("L1.DY")];
        if wrong == 0 {
            expected.push(at("L0.DY"));
        }
        defects.push(Defect {
            id: "bad_study_day",
            dataset: name.into(),
            row: i,
            usubjid: chars(ds, "USUBJID").and_then(|u| u[i].clone()),
            key: seq_key(ds, i),
            variable: var.clone(),
            original: Some(day.to_string()),
            injected: Some(wrong.to_string()),
            paired_row: None,
            expected,
        });
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: var.clone(),
            row: i,
            value: Cell::Int(Some(wrong)),
        });
    }
    for (name, e) in &edits {
        cx.edit_sdtm(name, e);
    }
    defects
}
