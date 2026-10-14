//! L0: structural checks against the metadata bundle.

use std::collections::{HashMap, HashSet};

use trialforge_cdisc::{Core, VarType, ct, spec};
use trialforge_core::{Dataset, Values};

use crate::table::{T, iso_ok};
use crate::{Level, Report};

const L: Level = Level::L0;

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    for ds in sdtm {
        dataset(ds, true, r);
    }
    for ds in adam {
        dataset(ds, false, r);
    }
}

fn dataset(ds: &Dataset, is_sdtm: bool, r: &mut Report) {
    let name = ds.name.as_str();
    let Some(spec) = spec(name) else {
        r.add(
            L,
            "L0.SPEC",
            name,
            None,
            None,
            "dataset is not in the metadata bundle".into(),
        );
        return;
    };
    let t = T::new(ds);

    // Presence and population of required/expected variables.
    for v in &spec.variables {
        let col = ds.column(&v.name);
        match (v.core, col) {
            (Core::Req, None) => r.add(
                L,
                "L0.REQ",
                name,
                Some(&v.name),
                None,
                "required variable is missing".into(),
            ),
            (Core::Exp, None) if is_sdtm => r.add(
                L,
                "L0.EXP",
                name,
                Some(&v.name),
                None,
                "expected variable is missing".into(),
            ),
            (Core::Req, Some(c)) => {
                let nulls = match &c.values {
                    Values::Char(x) => x
                        .iter()
                        .filter(|v| v.as_deref().is_none_or(str::is_empty))
                        .count(),
                    Values::Num(x) => x.iter().filter(|v| v.is_none()).count(),
                    Values::Int(x) => x.iter().filter(|v| v.is_none()).count(),
                    Values::Date(x) => x.iter().filter(|v| v.is_none()).count(),
                    Values::DateTime(x) => x.iter().filter(|v| v.is_none()).count(),
                };
                if nulls > 0 {
                    r.add(
                        L,
                        "L0.REQ",
                        name,
                        Some(&v.name),
                        None,
                        format!("required variable has {nulls} null value(s)"),
                    );
                }
            }
            _ => {}
        }
    }

    for c in ds.columns() {
        let Some(v) = spec.variable(&c.name) else {
            r.add(
                L,
                "L0.SPEC",
                name,
                Some(&c.name),
                None,
                "variable is not in the metadata bundle".into(),
            );
            continue;
        };
        // Type.
        let is_char = matches!(c.values, Values::Char(_));
        if is_char != (v.var_type == VarType::Char) {
            r.add(
                L,
                "L0.TYPE",
                name,
                Some(&c.name),
                None,
                format!("expected {:?}", v.var_type),
            );
        }
        if let Values::Char(vals) = &c.values {
            for (i, x) in vals.iter().enumerate() {
                let Some(x) = x else { continue };
                if x.len() > 200 || !x.is_ascii() {
                    r.add(
                        L,
                        "L0.LEN",
                        name,
                        Some(&c.name),
                        Some(i),
                        "value is longer than 200 bytes or not ASCII".into(),
                    );
                }
                // Controlled terminology.
                if let Some(cl) = &v.codelist {
                    if let Some(list) = ct().codelist(cl) {
                        if !list.contains(x) {
                            r.add(
                                L,
                                "L0.CT",
                                name,
                                Some(&c.name),
                                Some(i),
                                format!("{x:?} is not in codelist {cl}"),
                            );
                        }
                    }
                }
                // ISO 8601.
                if c.name.ends_with("DTC") && !iso_ok(x) {
                    r.add(
                        L,
                        "L0.ISO",
                        name,
                        Some(&c.name),
                        Some(i),
                        format!("{x:?} is not ISO 8601"),
                    );
                }
            }
        }
        // No study day 0.
        if c.name.ends_with("DY") && c.name != "VISITDY" {
            if let Values::Int(vals) = &c.values {
                for (i, x) in vals.iter().enumerate() {
                    if *x == Some(0) {
                        r.add(
                            L,
                            "L0.DY",
                            name,
                            Some(&c.name),
                            Some(i),
                            "study day 0 does not exist".into(),
                        );
                    }
                }
            }
        }
    }

    // DOMAIN matches the dataset.
    if let Some(dom) = t.c("DOMAIN") {
        for (i, d) in dom.iter().enumerate() {
            if d.as_deref() != Some(name) {
                r.add(
                    L,
                    "L0.DOMAIN",
                    name,
                    Some("DOMAIN"),
                    Some(i),
                    format!("DOMAIN is {d:?}"),
                );
            }
        }
    }

    // Natural keys unique.
    let keys: Vec<&str> = spec
        .keys
        .iter()
        .map(String::as_str)
        .filter(|k| ds.column(k).is_some())
        .collect();
    if !keys.is_empty() {
        let render = |i: usize| -> String {
            keys.iter()
                .map(|k| match &ds.column(k).expect("filtered").values {
                    Values::Char(v) => v[i].clone().unwrap_or_default(),
                    Values::Num(v) => v[i].map(|x| x.to_string()).unwrap_or_default(),
                    Values::Int(v) => v[i].map(|x| x.to_string()).unwrap_or_default(),
                    Values::Date(v) => v[i].map(|x| x.to_string()).unwrap_or_default(),
                    Values::DateTime(v) => v[i].map(|x| x.to_string()).unwrap_or_default(),
                })
                .collect::<Vec<_>>()
                .join("|")
        };
        let mut seen = HashSet::new();
        for i in 0..ds.n_rows() {
            if !seen.insert(render(i)) {
                r.add(
                    L,
                    "L0.KEY",
                    name,
                    None,
                    Some(i),
                    format!("duplicate natural key ({})", keys.join(", ")),
                );
            }
        }
    }

    // --SEQ unique within subject.
    let seq_var = format!("{name}SEQ");
    if let (Some(u), Some(seq)) = (t.c("USUBJID"), t.i(&seq_var)) {
        let mut seen: HashMap<(&str, i64), usize> = HashMap::new();
        for i in 0..t.n() {
            if let (Some(u), Some(s)) = (u[i].as_deref(), seq[i]) {
                if seen.insert((u, s), i).is_some() {
                    r.add(
                        L,
                        "L0.SEQ",
                        name,
                        Some(&seq_var),
                        Some(i),
                        format!("{seq_var} {s} repeats for {u}"),
                    );
                }
            }
        }
    }
}
