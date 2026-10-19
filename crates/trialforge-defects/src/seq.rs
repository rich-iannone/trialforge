//! Key defects.

use std::collections::{BTreeMap, HashMap, HashSet};

use trialforge_cdisc::spec;
use trialforge_core::{Dataset, Draws, Stream, Values};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, parse_date, pick, seq_key};

/// Domains whose --SEQ no other dataset or check refers to, so a duplicate breaks only the
/// uniqueness rules (an AESEQ duplicate would also break ADAE and SUPPAE traceability).
const DOMAINS: &[&str] = &["CM", "MH", "DS", "SE", "TU"];

/// `duplicate_seq`: a record takes the --SEQ of the same subject's previous record.
pub(crate) fn duplicate_seq(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    // (dataset, row, previous row of the subject).
    let mut candidates: Vec<(&str, usize, usize)> = Vec::new();
    for ds in cx
        .sdtm
        .iter()
        .filter(|d| DOMAINS.contains(&d.name.as_str()))
    {
        let name = ds.name.as_str();
        let (Some(u), Some(seq)) = (chars(ds, "USUBJID"), ints(ds, &format!("{name}SEQ"))) else {
            continue;
        };
        let mut last: HashMap<&str, usize> = HashMap::new();
        for i in 0..ds.n_rows() {
            let (Some(x), Some(_)) = (u[i].as_deref(), seq[i]) else {
                continue;
            };
            if let Some(j) = last.insert(x, i) {
                candidates.push((name, i, j));
            }
        }
    }
    // A record is either the copy or the original of at most one duplicate.
    let mut taken = cx.used.clone();
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |c| {
            let (name, i, j) = candidates[c];
            let (a, b) = ((name.to_string(), i), (name.to_string(), j));
            if taken.contains(&a) || taken.contains(&b) {
                return false;
            }
            taken.insert(a);
            taken.insert(b);
            true
        },
    );

    let mut defects = Vec::new();
    let mut edits: BTreeMap<String, Vec<Edit>> = BTreeMap::new();
    for c in chosen {
        let (name, i, j) = candidates[c];
        let ds = cx.sdtm(name).expect("a candidate's dataset exists");
        let var = format!("{name}SEQ");
        let seq = ints(ds, &var).expect("a candidate's dataset has --SEQ");
        let (original, copied) = (seq[i].expect("candidate"), seq[j].expect("candidate"));
        let mut expected = vec![Expected {
            rule: "L0.SEQ",
            dataset: name.into(),
            row: Some(i),
        }];
        if key_collides(ds, &var, i, j) {
            expected.push(Expected {
                rule: "L0.KEY",
                dataset: name.into(),
                row: Some(i),
            });
        }
        defects.push(Defect {
            id: "duplicate_seq",
            dataset: name.into(),
            row: i,
            usubjid: chars(ds, "USUBJID").and_then(|u| u[i].clone()),
            key: Some(format!("{var}={original}")),
            variable: var.clone(),
            original: Some(original.to_string()),
            injected: Some(copied.to_string()),
            paired_row: Some(j),
            expected,
        });
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: var,
            row: i,
            value: Cell::Int(Some(copied)),
        });
    }
    for (name, e) in &edits {
        cx.edit_sdtm(name, e);
    }
    defects
}

/// Whether row `i`, given row `j`'s `seq_var`, has the same natural key as row `j` (only
/// possible when the key includes --SEQ and the rest of it already agrees).
fn key_collides(ds: &Dataset, seq_var: &str, i: usize, j: usize) -> bool {
    let Some(spec) = spec(&ds.name) else {
        return false;
    };
    let keys: Vec<&str> = spec
        .keys
        .iter()
        .map(String::as_str)
        .filter(|k| ds.column(k).is_some())
        .collect();
    !keys.is_empty()
        && keys.contains(&seq_var)
        && keys.iter().filter(|k| **k != seq_var).all(|k| {
            match &ds.column(k).expect("filtered").values {
                Values::Char(v) => v[i] == v[j],
                Values::Num(v) => v[i] == v[j],
                Values::Int(v) => v[i] == v[j],
                Values::Date(v) => v[i] == v[j],
                Values::DateTime(v) => v[i] == v[j],
            }
        })
}

/// Domains whose records `orphan_usubjid` re-assigns: nothing derived (other than ADAE) is
/// keyed by their subject, so a stranger's record trips only the rules listed in its manifest.
const ORPHAN_DOMAINS: &[&str] = &["AE", "CM", "MH"];

/// Attempts at an unused subject identifier before a record is left alone.
const ORPHAN_TRIES: usize = 20;

/// `orphan_usubjid`: one digit of a record's USUBJID changes, giving a subject not in DM
/// (`L1.SUBJ`). The record also falls outside every SE element (`L1.EPOCH`, for a record with
/// an EPOCH), and an AE leaves its SUPPAE qualifiers (`L1.SUPP`) and its ADAE record
/// (`L2.ADAE`) without a parent. Dose-limiting AEs are left alone.
pub(crate) fn orphan_usubjid(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    // DM's subjects, and every orphan made so far (two orphans never share a subject).
    let mut known: HashSet<String> = cx
        .sdtm("DM")
        .and_then(|dm| chars(dm, "USUBJID"))
        .map(|u| u.iter().flatten().cloned().collect())
        .unwrap_or_default();
    let dlts = dose_limiting(cx);
    let mut candidates: Vec<(&str, usize)> = Vec::new();
    for ds in cx
        .sdtm
        .iter()
        .filter(|d| ORPHAN_DOMAINS.contains(&d.name.as_str()))
    {
        let Some(u) = chars(ds, "USUBJID") else {
            continue;
        };
        let name = ds.name.as_str();
        let seq = ints(ds, &format!("{name}SEQ"));
        for i in 0..ds.n_rows() {
            let dlt = name == "AE"
                && u[i]
                    .as_deref()
                    .zip(seq.and_then(|q| q[i]))
                    .is_some_and(|(x, q)| dlts.contains(&(x.to_string(), q.to_string())));
            if !cx.used.contains(&(name.to_string(), i))
                && !dlt
                && u[i]
                    .as_deref()
                    .is_some_and(|x| x.ends_with(|c: char| c.is_ascii_digit()))
            {
                candidates.push((name, i));
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
        let (name, i) = candidates[c];
        let ds = cx.sdtm(name).expect("a candidate's dataset exists");
        let x = chars(ds, "USUBJID").expect("candidate")[i]
            .clone()
            .expect("a candidate has a subject");
        let mut d = Draws::new(stream.label("value").label(name).index(i as u64));
        let Some(orphan) = mistype(&x, &known, &mut d) else {
            continue;
        };
        known.insert(orphan.clone());
        let at = |rule, dataset: &str, row| Expected {
            rule,
            dataset: dataset.into(),
            row: Some(row),
        };
        let mut expected = vec![at("L1.SUBJ", name, i)];
        let date_var = [format!("{name}STDTC"), format!("{name}DTC")]
            .into_iter()
            .find(|v| ds.column(v).is_some());
        let date = date_var.and_then(|v| chars(ds, &v)?[i].as_deref().and_then(parse_date));
        if let Some(date) = date
            && cx.outside_epoch(name, i, &orphan, date)
        {
            expected.push(at("L1.EPOCH", name, i));
        }
        if name == "AE"
            && let Some(q) = ints(ds, "AESEQ").and_then(|q| q[i])
        {
            if let Some(sup) = cx.sdtm("SUPPAE")
                && let (Some(su), Some(rdom), Some(idvar), Some(val)) = (
                    chars(sup, "USUBJID"),
                    chars(sup, "RDOMAIN"),
                    chars(sup, "IDVAR"),
                    chars(sup, "IDVARVAL"),
                )
            {
                let q = q.to_string();
                for k in 0..sup.n_rows() {
                    if su[k].as_deref() == Some(x.as_str())
                        && rdom[k].as_deref() == Some("AE")
                        && idvar[k].as_deref() == Some("AESEQ")
                        && val[k].as_deref() == Some(q.as_str())
                    {
                        expected.push(at("L1.SUPP", "SUPPAE", k));
                    }
                }
            }
            if let Some(j) = cx.adae_row(&x, q) {
                expected.push(at("L2.ADAE", "ADAE", j));
            }
        }
        defects.push(Defect {
            id: "orphan_usubjid",
            dataset: name.into(),
            row: i,
            usubjid: Some(x.clone()),
            key: seq_key(ds, i),
            variable: "USUBJID".into(),
            original: Some(x),
            injected: Some(orphan.clone()),
            paired_row: None,
            expected,
        });
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: "USUBJID".into(),
            row: i,
            value: Cell::Char(Some(orphan)),
        });
    }
    for (name, e) in &edits {
        cx.edit_sdtm(name, e);
    }
    defects
}

/// (USUBJID, AESEQ) of the AEs SUPPAE flags as dose-limiting (`AEDLT` = 'Y'). These are not
/// orphaned: the subject losing its DLT would change ADSL's DLTEVLFL as well.
fn dose_limiting(cx: &Ctx<'_>) -> HashSet<(String, String)> {
    let Some(sup) = cx.sdtm("SUPPAE") else {
        return HashSet::new();
    };
    let (Some(u), Some(q), Some(v), Some(val)) = (
        chars(sup, "USUBJID"),
        chars(sup, "QNAM"),
        chars(sup, "IDVARVAL"),
        chars(sup, "QVAL"),
    ) else {
        return HashSet::new();
    };
    (0..sup.n_rows())
        .filter(|&k| q[k].as_deref() == Some("AEDLT") && val[k].as_deref() == Some("Y"))
        .filter_map(|k| Some((u[k].clone()?, v[k].clone()?)))
        .collect()
}

/// `usubjid` with one digit of its trailing subject number changed, so that it is not in
/// `known`; `None` if [`ORPHAN_TRIES`] attempts all collide.
fn mistype(usubjid: &str, known: &HashSet<String>, d: &mut Draws) -> Option<String> {
    let digits = usubjid.bytes().rev().take_while(u8::is_ascii_digit).count();
    let start = usubjid.len() - digits;
    for _ in 0..ORPHAN_TRIES {
        let at = start + d.below(digits as u64) as usize;
        let old = usubjid.as_bytes()[at] - b'0';
        let new = (old + 1 + d.below(9) as u8) % 10;
        let mut s = usubjid.as_bytes().to_vec();
        s[at] = b'0' + new;
        let s = String::from_utf8(s).expect("ASCII digits");
        if !known.contains(&s) {
            return Some(s);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mistyped_subjects_differ_in_one_digit() {
        let known: HashSet<String> = ["S-1001", "S-1002"].map(String::from).into();
        let mut d = Draws::new(trialforge_core::Stream::root(3).label("t"));
        for _ in 0..100 {
            let m = mistype("S-1001", &known, &mut d).unwrap();
            assert!(!known.contains(&m));
            assert!(m.starts_with("S-"));
            let diff = m
                .bytes()
                .zip("S-1001".bytes())
                .filter(|(a, b)| a != b)
                .count();
            assert_eq!(diff, 1, "{m}");
        }
    }
}
