//! Controlled-terminology defects.

use std::collections::BTreeMap;

use trialforge_cdisc::{Codelist, ct, spec};
use trialforge_core::{Draws, Stream};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, pick};

/// Controlled variables that no check other than L0.CT reads, so an out-of-codelist value is
/// caught by that rule alone (SEX, EPOCH and the test codes, say, feed derivations).
const VARIABLES: &[(&str, &str)] = &[
    ("DM", "RACE"),
    ("DM", "ETHNIC"),
    ("AE", "AESEV"),
    ("AE", "AEACN"),
    ("AE", "AEREL"),
    ("AE", "AEOUT"),
    ("CM", "CMROUTE"),
    ("CM", "CMDOSFRQ"),
    ("CM", "CMDOSU"),
    ("EX", "EXROUTE"),
    ("EX", "EXDOSFRM"),
    ("MH", "MHENRTPT"),
    ("DS", "DSCAT"),
    ("VS", "VSPOS"),
    ("LB", "LBSPEC"),
    ("TU", "TUMETHOD"),
    ("CE", "CEOUT"),
];

/// (dataset, row, the controlled variables holding a valid term, with their codelists).
type Candidate<'a> = (&'a str, usize, Vec<(&'static str, &'static Codelist)>);

/// `out_of_ct`: a controlled value becomes a case variant ("Mild") or a typo ("MIDL") that
/// is not in its codelist (subtle: a case variant wherever one exists).
pub(crate) fn out_of_ct(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let mut candidates: Vec<Candidate<'_>> = Vec::new();
    for ds in cx.sdtm.iter() {
        let name = ds.name.as_str();
        let Some(dspec) = spec(name) else { continue };
        let vars: Vec<(&str, &[Option<String>], &Codelist)> = VARIABLES
            .iter()
            .filter(|(d, _)| *d == name)
            .filter_map(|&(_, v)| {
                let list = ct().codelist(dspec.variable(v)?.codelist.as_deref()?)?;
                Some((v, chars(ds, v)?, list))
            })
            .collect();
        if vars.is_empty() {
            continue;
        }
        for i in 0..ds.n_rows() {
            if cx.used.contains(&(name.to_string(), i)) {
                continue;
            }
            let ok: Vec<(&str, &Codelist)> = vars
                .iter()
                .filter(|(_, x, list)| x[i].as_deref().is_some_and(|x| list.contains(x)))
                .map(|&(v, _, list)| (v, list))
                .collect();
            if !ok.is_empty() {
                candidates.push((name, i, ok));
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
        let (var, list) = vars[d.below(vars.len() as u64) as usize];
        let original = chars(ds, var).expect("candidate")[i]
            .clone()
            .expect("a candidate holds a term");
        let Some(injected) = variant(&original, list, cx.subtle, &mut d) else {
            continue;
        };
        let seq_var = format!("{name}SEQ");
        defects.push(Defect {
            id: "out_of_ct",
            dataset: name.into(),
            row: i,
            usubjid: chars(ds, "USUBJID").and_then(|u| u[i].clone()),
            key: ints(ds, &seq_var)
                .and_then(|s| s[i])
                .map(|s| format!("{seq_var}={s}")),
            variable: var.into(),
            original: Some(original),
            injected: Some(injected.clone()),
            paired_row: None,
            expected: vec![Expected {
                rule: "L0.CT",
                dataset: name.into(),
                row: Some(i),
            }],
        });
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: var.into(),
            row: i,
            value: Cell::Char(Some(injected)),
        });
    }
    for (name, e) in &edits {
        cx.edit_sdtm(name, e);
    }
    defects
}

/// A case variant or a typo of `term` that is not in `list`, each kind with probability one
/// half or, if `case_only`, a case variant (whichever exists, if only one does).
fn variant(term: &str, list: &Codelist, case_only: bool, d: &mut Draws) -> Option<String> {
    let outside = |v: &String| v != term && !list.contains(v);
    let mut case: Vec<String> = vec![title_case(term), term.to_ascii_lowercase()];
    case.dedup();
    case.retain(outside);
    // Adjacent characters swapped, or one doubled.
    let c: Vec<char> = term.chars().collect();
    let mut typo: Vec<String> = Vec::new();
    for p in 0..c.len() {
        if p + 1 < c.len()
            && c[p] != c[p + 1]
            && c[p].is_alphanumeric()
            && c[p + 1].is_alphanumeric()
        {
            let mut t = c.clone();
            t.swap(p, p + 1);
            typo.push(t.into_iter().collect());
        }
        if c[p].is_alphanumeric() {
            let mut t = c.clone();
            t.insert(p, c[p]);
            typo.push(t.into_iter().collect());
        }
    }
    typo.dedup();
    typo.retain(outside);
    let use_case = d.bernoulli(0.5) || case_only;
    let from = match (use_case, case.is_empty(), typo.is_empty()) {
        (_, true, true) => return None,
        (true, false, _) | (false, false, true) => &case,
        _ => &typo,
    };
    Some(from[d.below(from.len() as u64) as usize].clone())
}

/// `NOT RECOVERED/NOT RESOLVED` → `Not Recovered/not Resolved`.
fn title_case(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut start = true;
    for ch in s.chars() {
        if start {
            out.push(ch.to_ascii_uppercase());
        } else {
            out.push(ch.to_ascii_lowercase());
        }
        start = ch == ' ';
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn variants_leave_the_codelist() {
        let list = ct().codelist("AESEV").expect("AESEV codelist");
        for seed in 0..50 {
            let mut d = Draws::new(Stream::root(seed).label("t"));
            for term in &list.terms {
                let v = variant(term, list, false, &mut d).expect("a variant exists");
                assert_ne!(&v, term);
                assert!(!list.contains(&v), "{v}");
                let c = variant(term, list, true, &mut d).expect("a variant exists");
                assert!(c.eq_ignore_ascii_case(term) && c != *term, "{c}");
            }
        }
        assert_eq!(
            title_case("NOT RECOVERED/NOT RESOLVED"),
            "Not Recovered/not Resolved"
        );
    }
}
