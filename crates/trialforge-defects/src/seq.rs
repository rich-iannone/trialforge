//! Key defects.

use std::collections::{BTreeMap, HashMap};

use trialforge_cdisc::spec;
use trialforge_core::{Dataset, Draws, Stream, Values};

use crate::{Amount, Cell, Ctx, Defect, Edit, Expected, chars, ints, pick};

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
    for d in &defects {
        cx.used.insert((d.dataset.clone(), d.row));
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
