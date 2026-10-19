//! ADaM derivation defects.

use std::collections::BTreeMap;

use trialforge_core::{Draws, Stream};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, nums, pick};

/// BDS datasets whose CHG and ABLFL `L2.BDS` recomputes from SDTM.
const BDS: &[&str] = &["ADLB", "ADVS", "ADRE"];

/// The smallest CHG error injected: well above the check's tolerance.
const MIN_ERROR: f64 = 0.001;

/// `{:.6}` without trailing zeros.
fn number(x: f64) -> String {
    let s = format!("{x:.6}");
    let s = s.trim_end_matches('0').trim_end_matches('.');
    if s == "-0" { "0".into() } else { s.into() }
}

/// (dataset, row) of BDS records that pass `keep`, in dataset order.
fn bds_rows<'a>(cx: &'a Ctx<'_>, keep: impl Fn(&str, usize) -> bool) -> Vec<(&'a str, usize)> {
    let mut out = Vec::new();
    for ds in cx.adam.iter().filter(|d| BDS.contains(&d.name.as_str())) {
        let name = ds.name.as_str();
        for i in 0..ds.n_rows() {
            if !cx.is_used(name, i) && keep(name, i) {
                out.push((name, i));
            }
        }
    }
    out
}

fn defect(
    cx: &Ctx<'_>,
    id: &'static str,
    (name, i): (&str, usize),
    variable: &str,
    original: Option<String>,
    injected: Option<String>,
    rule: &'static str,
) -> Defect {
    let ds = cx.adam(name).expect("a candidate's dataset exists");
    let key = match (chars(ds, "PARAMCD"), ints(ds, "SRCSEQ")) {
        (Some(p), Some(q)) => p[i].as_deref().map(|p| match q[i] {
            Some(q) => format!("PARAMCD={p}, SRCSEQ={q}"),
            None => format!("PARAMCD={p}"),
        }),
        _ => None,
    };
    Defect {
        id,
        dataset: name.into(),
        row: i,
        usubjid: chars(ds, "USUBJID").and_then(|u| u[i].clone()),
        key,
        variable: variable.into(),
        original,
        injected,
        paired_row: None,
        expected: vec![Expected {
            rule,
            dataset: name.into(),
            row: Some(i),
        }],
    }
}

/// `chg_mismatch`: a post-baseline CHG becomes BASE - AVAL (the subtraction reversed) or AVAL
/// itself (the subtraction forgotten), each with probability one half where it differs from
/// the right value; subtle: CHG is off by 1% of BASE (at least 0.01).
pub(crate) fn chg_mismatch(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let values = |name: &str, i: usize| -> Option<(f64, f64, f64)> {
        let ds = cx.adam(name)?;
        Some((
            nums(ds, "AVAL")?[i]?,
            nums(ds, "BASE")?[i]?,
            nums(ds, "CHG")?[i]?,
        ))
    };
    let candidates = bds_rows(cx, |name, i| {
        values(name, i)
            .is_some_and(|(_, base, chg)| chg.abs() >= MIN_ERROR || base.abs() >= MIN_ERROR)
    });
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
        let (aval, base, chg) = values(name, i).expect("candidate");
        let mut d = Draws::new(stream.label("value").label(name).index(i as u64));
        let reverse = d.bernoulli(0.5);
        let up = d.bernoulli(0.5);
        let wrong = if cx.subtle {
            let e = (base.abs() * 0.01).max(0.01);
            chg + if up { e } else { -e }
        } else if (reverse && chg.abs() >= MIN_ERROR) || base.abs() < MIN_ERROR {
            base - aval
        } else {
            aval
        };
        // Six decimals, as reported.
        let wrong = (wrong * 1e6).round() / 1e6;
        defects.push(defect(
            cx,
            "chg_mismatch",
            (name, i),
            "CHG",
            Some(number(chg)),
            Some(number(wrong)),
            "L2.BDS",
        ));
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: "CHG".into(),
            row: i,
            value: Cell::Num(Some(wrong)),
        });
    }
    for (name, e) in &edits {
        cx.edit_adam(name, e);
    }
    defects
}

/// `missing_baseline`: ABLFL is cleared on a baseline record (BASE and CHG are left as derived).
pub(crate) fn missing_baseline(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let candidates = bds_rows(cx, |name, i| {
        cx.adam(name)
            .and_then(|ds| chars(ds, "ABLFL"))
            .is_some_and(|f| f[i].as_deref() == Some("Y"))
    });
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
        defects.push(defect(
            cx,
            "missing_baseline",
            (name, i),
            "ABLFL",
            Some("Y".into()),
            None,
            "L2.BDS",
        ));
        edits.entry(name.to_string()).or_default().push(Edit {
            variable: "ABLFL".into(),
            row: i,
            value: Cell::Char(None),
        });
    }
    for (name, e) in &edits {
        cx.edit_adam(name, e);
    }
    defects
}

/// `trtemfl_wrong`: ADAE TRTEMFL flips, 'Y' to missing or missing to 'Y'.
pub(crate) fn trtemfl_wrong(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let Some(ad) = cx.adam("ADAE") else {
        return Vec::new();
    };
    let (Some(u), Some(seq), Some(fl)) = (
        chars(ad, "USUBJID"),
        ints(ad, "AESEQ"),
        chars(ad, "TRTEMFL"),
    ) else {
        return Vec::new();
    };
    let candidates: Vec<usize> = (0..ad.n_rows())
        .filter(|&i| !cx.is_used("ADAE", i) && u[i].is_some() && seq[i].is_some())
        .collect();
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );
    let mut defects = Vec::new();
    let mut edits = Vec::new();
    for c in chosen {
        let i = candidates[c];
        let flipped = match fl[i].as_deref() {
            Some("Y") => None,
            _ => Some("Y".to_string()),
        };
        let mut d = defect(
            cx,
            "trtemfl_wrong",
            ("ADAE", i),
            "TRTEMFL",
            fl[i].clone(),
            flipped.clone(),
            "L2.ADAE",
        );
        d.key = seq[i].map(|s| format!("AESEQ={s}"));
        defects.push(d);
        edits.push(Edit {
            variable: "TRTEMFL".into(),
            row: i,
            value: Cell::Char(flipped),
        });
    }
    cx.edit_adam("ADAE", &edits);
    defects
}

/// `cnsr_wrong`: ADTTE CNSR flips between 0 (event) and 1 (censored). The time-to-first
/// exacerbation parameter is checked by `L2.ADEXAC`, the others by `L2.ADTTE`.
pub(crate) fn cnsr_wrong(cx: &mut Ctx<'_>, stream: Stream, amount: Amount) -> Vec<Defect> {
    let Some(tte) = cx.adam("ADTTE") else {
        return Vec::new();
    };
    let (Some(pc), Some(cnsr)) = (chars(tte, "PARAMCD"), ints(tte, "CNSR")) else {
        return Vec::new();
    };
    let candidates: Vec<usize> = (0..tte.n_rows())
        .filter(|&i| !cx.is_used("ADTTE", i) && pc[i].is_some() && matches!(cnsr[i], Some(0 | 1)))
        .collect();
    let chosen = pick(
        &mut Draws::new(stream.label("select")),
        candidates.len(),
        amount,
        |_| true,
    );
    let mut defects = Vec::new();
    let mut edits = Vec::new();
    for c in chosen {
        let i = candidates[c];
        let was = cnsr[i].expect("candidate");
        let rule = if pc[i].as_deref() == Some("TTFEXAC") {
            "L2.ADEXAC"
        } else {
            "L2.ADTTE"
        };
        let mut d = defect(
            cx,
            "cnsr_wrong",
            ("ADTTE", i),
            "CNSR",
            Some(was.to_string()),
            Some((1 - was).to_string()),
            rule,
        );
        d.key = pc[i].as_deref().map(|p| format!("PARAMCD={p}"));
        defects.push(d);
        edits.push(Edit {
            variable: "CNSR".into(),
            row: i,
            value: Cell::Int(Some(1 - was)),
        });
    }
    cx.edit_adam("ADTTE", &edits);
    defects
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn numbers_are_trimmed() {
        assert_eq!(number(1.5), "1.5");
        assert_eq!(number(-0.0000001), "0");
        assert_eq!(number(12.0), "12");
        assert_eq!(number(0.1 + 0.2), "0.3");
    }
}
