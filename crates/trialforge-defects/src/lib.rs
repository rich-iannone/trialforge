//! Seeded data-quality defects, with a manifest of what was injected and where.
//!
//! A study is generated and verified clean first; [`inject`] then corrupts copies of its SDTM
//! and ADaM datasets. A defect touches one record of one dataset and nothing downstream: an
//! SDTM defect does not propagate to ADaM (the ADaM stays derived from the clean SDTM, so the
//! disagreement is itself what a validator should catch). Variables that depend on the edited
//! one in the same record are kept consistent, so a defect reads as one mistake rather than
//! several: moving `AEENDTC` also moves `AEENDY`.
//!
//! Each [`Defect`] lists the findings it should produce ([`Expected`]): the rule of
//! `trialforge-check` that catches it and the dataset and row. [`score`] matches a
//! validator's findings against them, for precision and recall per defect type.
//!
//! Randomness comes from the stream `study → "defects" → <defect id>`, so adding a defect type
//! to a request does not move the others. Defect types are applied in [`CATALOGUE`] order and a
//! record receives at most one defect.

mod ct;
mod dates;
mod score;
mod seq;

use std::collections::{BTreeMap, HashSet};

use serde::{Deserialize, Deserializer, Serialize};
use serde_json::Value;
use trialforge_core::{Dataset, Date, DateTime, Draws, Stream, Values};

pub use score::{Observed, Score, TypeScore, score};

/// A defect type.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct DefectType {
    /// Identifier, e.g. `ae_end_before_start`.
    pub id: &'static str,
    /// What the defect does.
    pub description: &'static str,
    /// The `trialforge-check` rules that can catch it.
    pub rules: &'static [&'static str],
}

/// Every defect type, in the order they are applied.
pub const CATALOGUE: &[DefectType] = &[
    DefectType {
        id: "ae_end_before_start",
        description: "AEENDTC moved to before AESTDTC (AEENDY follows); ADAE keeps the clean dates",
        rules: &["L1.ORDER", "L2.ADAE"],
    },
    DefectType {
        id: "duplicate_seq",
        description: "a record's --SEQ set to that of the subject's previous record",
        rules: &["L0.SEQ", "L0.KEY"],
    },
    DefectType {
        id: "out_of_ct",
        description: "a controlled term replaced by a case variant or typo outside its codelist",
        rules: &["L0.CT"],
    },
];

/// Named defect profiles (not available yet).
const PROFILES: &[&str] = &["light", "realistic", "adversarial"];

/// The defect type called `id`.
pub fn defect_type(id: &str) -> Option<&'static DefectType> {
    CATALOGUE.iter().find(|t| t.id == id)
}

/// How many defects of one type to inject.
#[derive(Clone, Copy, Debug, PartialEq, Serialize)]
#[serde(untagged)]
pub enum Amount {
    /// This many records (all eligible ones if the study has fewer).
    Count(u64),
    /// Each eligible record independently with this probability.
    Rate(f64),
}

/// Which defects to inject: defect id → [`Amount`]. Empty means none.
///
/// Deserialises from `"none"` or from a map such as `{"out_of_ct": 5, "duplicate_seq": 0.02}`,
/// where an integer is a count and any other number a rate in `[0, 1]`.
#[derive(Clone, Debug, Default, PartialEq, Serialize)]
#[serde(transparent)]
pub struct DefectSpec {
    /// Amount per defect id.
    pub defects: BTreeMap<String, Amount>,
}

impl DefectSpec {
    /// Whether no defect type is requested.
    pub fn is_empty(&self) -> bool {
        self.defects.is_empty()
    }

    /// Parses and validates a spec from JSON (see the type docs).
    pub fn from_value(v: &Value) -> Result<Self, String> {
        let usage = "pass \"none\" or a map of defect id to count or rate";
        match v {
            Value::Null => Ok(Self::default()),
            Value::String(s) if s == "none" => Ok(Self::default()),
            Value::String(s) if PROFILES.contains(&s.as_str()) => Err(format!(
                "defect profile {s:?} is not available yet; {usage}"
            )),
            Value::String(s) => Err(format!("unknown defect profile {s:?}; {usage}")),
            Value::Object(m) => {
                let mut defects = BTreeMap::new();
                for (id, a) in m {
                    if defect_type(id).is_none() {
                        let ids: Vec<&str> = CATALOGUE.iter().map(|t| t.id).collect();
                        return Err(format!(
                            "unknown defect {id:?}; available: {}",
                            ids.join(", ")
                        ));
                    }
                    let amount = match (a.as_u64(), a.as_f64()) {
                        (Some(k), _) => Amount::Count(k),
                        (None, Some(p)) if (0.0..=1.0).contains(&p) => Amount::Rate(p),
                        _ => {
                            return Err(format!(
                                "defect {id:?} needs a count (a non-negative integer) or a rate \
                                 in [0, 1], not {a}"
                            ));
                        }
                    };
                    defects.insert(id.clone(), amount);
                }
                Ok(Self { defects })
            }
            _ => Err(format!("`defects`: {usage}")),
        }
    }
}

impl<'de> Deserialize<'de> for DefectSpec {
    fn deserialize<D: Deserializer<'de>>(d: D) -> Result<Self, D::Error> {
        Self::from_value(&Value::deserialize(d)?).map_err(serde::de::Error::custom)
    }
}

/// A finding a defect should produce.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct Expected {
    /// The `trialforge-check` rule, e.g. `L1.ORDER`.
    pub rule: &'static str,
    /// Dataset.
    pub dataset: String,
    /// Zero-based row, if the finding is about one.
    pub row: Option<usize>,
}

/// One injected defect.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Defect {
    /// Defect type id.
    pub id: &'static str,
    /// Dataset edited.
    pub dataset: String,
    /// Zero-based row edited.
    pub row: usize,
    /// `USUBJID` of the record.
    pub usubjid: Option<String>,
    /// The record's identifying key within the subject as it was before the edit, e.g.
    /// `AESEQ=3` (`None` for a subject-level record).
    pub key: Option<String>,
    /// Variable edited.
    pub variable: String,
    /// Value before.
    pub original: Option<String>,
    /// Value after.
    pub injected: Option<String>,
    /// Zero-based row of another record the defect involves, which a validator may name
    /// instead (the record whose --SEQ a `duplicate_seq` copied).
    pub paired_row: Option<usize>,
    /// Findings this defect should produce.
    pub expected: Vec<Expected>,
}

/// The defects of a study and the request that produced them.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Manifest {
    /// What was requested.
    pub spec: DefectSpec,
    /// What was injected, by type in [`CATALOGUE`] order, then dataset and row.
    pub defects: Vec<Defect>,
}

/// Injects the defects of `spec` into `sdtm` and `adam`, drawing from `stream` (the study's
/// `"defects"` stream). Types not in [`CATALOGUE`] are rejected by [`DefectSpec::from_value`].
pub fn inject(
    stream: Stream,
    spec: &DefectSpec,
    sdtm: &mut [Dataset],
    adam: &mut [Dataset],
) -> Manifest {
    let mut cx = Ctx {
        sdtm,
        adam,
        used: HashSet::new(),
    };
    let mut defects = Vec::new();
    for t in CATALOGUE {
        let Some(&amount) = spec.defects.get(t.id) else {
            continue;
        };
        let s = stream.label(t.id);
        defects.extend(match t.id {
            "ae_end_before_start" => dates::ae_end_before_start(&mut cx, s, amount),
            "duplicate_seq" => seq::duplicate_seq(&mut cx, s, amount),
            "out_of_ct" => ct::out_of_ct(&mut cx, s, amount),
            _ => unreachable!("every catalogue entry is implemented"),
        });
    }
    Manifest {
        spec: spec.clone(),
        defects,
    }
}

/// The datasets being corrupted, and the records already touched.
pub(crate) struct Ctx<'a> {
    pub sdtm: &'a mut [Dataset],
    pub adam: &'a mut [Dataset],
    /// (dataset, row) of every record that already has a defect.
    pub used: HashSet<(String, usize)>,
}

impl Ctx<'_> {
    pub fn is_used(&self, dataset: &str, row: usize) -> bool {
        self.used.contains(&(dataset.to_string(), row))
    }

    pub fn sdtm(&self, name: &str) -> Option<&Dataset> {
        self.sdtm.iter().find(|d| d.name == name)
    }

    pub fn adam(&self, name: &str) -> Option<&Dataset> {
        self.adam.iter().find(|d| d.name == name)
    }

    /// Applies `edits` to the SDTM dataset `name`.
    pub fn edit_sdtm(&mut self, name: &str, edits: &[Edit]) {
        if let Some(ds) = self.sdtm.iter_mut().find(|d| d.name == name) {
            apply(ds, edits);
        }
    }
}

/// A new value for one cell.
pub(crate) enum Cell {
    Char(Option<String>),
    Int(Option<i64>),
}

/// One cell edit.
pub(crate) struct Edit {
    pub variable: String,
    pub row: usize,
    pub value: Cell,
}

fn apply(ds: &mut Dataset, edits: &[Edit]) {
    if edits.is_empty() {
        return;
    }
    let mut columns = ds.columns().to_vec();
    for e in edits {
        let col = columns
            .iter_mut()
            .find(|c| c.name == e.variable)
            .expect("an edited variable exists");
        match (&mut col.values, &e.value) {
            (Values::Char(v), Cell::Char(x)) => v[e.row].clone_from(x),
            (Values::Int(v), Cell::Int(x)) => v[e.row] = *x,
            _ => unreachable!("an edit matches its column's type"),
        }
    }
    *ds =
        Dataset::new(ds.name.clone(), ds.label.clone(), columns).expect("an edit keeps the shape");
}

/// Chooses among `n` candidates, in candidate order: `k` of them, visited in a random order (a
/// partial Fisher–Yates shuffle), or each with probability `p`. `accept` can turn a candidate
/// down, for example when it shares a record with one already chosen.
pub(crate) fn pick(
    d: &mut Draws,
    n: usize,
    amount: Amount,
    mut accept: impl FnMut(usize) -> bool,
) -> Vec<usize> {
    let mut chosen = Vec::new();
    match amount {
        Amount::Count(k) => {
            let mut order: Vec<usize> = (0..n).collect();
            for t in 0..n {
                if chosen.len() as u64 >= k {
                    break;
                }
                let j = t + d.below((n - t) as u64) as usize;
                order.swap(t, j);
                if accept(order[t]) {
                    chosen.push(order[t]);
                }
            }
        }
        Amount::Rate(p) => {
            for c in 0..n {
                // Draw for every candidate, so one rejection does not shift the rest.
                if d.bernoulli(p) && accept(c) {
                    chosen.push(c);
                }
            }
        }
    }
    chosen.sort_unstable();
    chosen
}

/// Character values of `var`.
pub(crate) fn chars<'a>(ds: &'a Dataset, var: &str) -> Option<&'a [Option<String>]> {
    ds.column(var)?.as_char()
}

/// Integer values of `var`.
pub(crate) fn ints<'a>(ds: &'a Dataset, var: &str) -> Option<&'a [Option<i64>]> {
    ds.column(var)?.as_int()
}

/// The date of a complete ISO 8601 date or date-time.
pub(crate) fn parse_date(s: &str) -> Option<Date> {
    Date::parse_iso(s).or_else(|| DateTime::parse_iso(s).map(DateTime::date))
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn spec_parses_counts_rates_and_none() {
        let s = DefectSpec::from_value(&json!({"out_of_ct": 3, "duplicate_seq": 0.25})).unwrap();
        assert_eq!(s.defects["out_of_ct"], Amount::Count(3));
        assert_eq!(s.defects["duplicate_seq"], Amount::Rate(0.25));
        assert!(DefectSpec::from_value(&json!("none")).unwrap().is_empty());
        assert_eq!(
            serde_json::to_value(&s).unwrap(),
            json!({"duplicate_seq": 0.25, "out_of_ct": 3})
        );
    }

    #[test]
    fn spec_rejects_bad_input() {
        for bad in [
            json!({"nope": 1}),
            json!({"out_of_ct": -1}),
            json!({"out_of_ct": 1.5}),
            json!({"out_of_ct": "many"}),
            json!("realistic"),
            json!("bogus"),
            json!([1]),
        ] {
            assert!(DefectSpec::from_value(&bad).is_err(), "{bad}");
        }
        let e = DefectSpec::from_value(&json!("light")).unwrap_err();
        assert!(e.contains("not available yet"), "{e}");
    }

    #[test]
    fn pick_takes_k_or_all() {
        let mut d = Draws::new(Stream::root(1).label("t"));
        let c = pick(&mut d, 10, Amount::Count(4), |_| true);
        assert_eq!(c.len(), 4);
        assert!(c.windows(2).all(|w| w[0] < w[1]));
        assert_eq!(pick(&mut d, 3, Amount::Count(9), |_| true), vec![0, 1, 2]);
        assert_eq!(
            pick(&mut d, 5, Amount::Rate(1.0), |i| i != 2),
            vec![0, 1, 3, 4]
        );
        assert!(pick(&mut d, 5, Amount::Rate(0.0), |_| true).is_empty());
        assert!(pick(&mut d, 0, Amount::Count(2), |_| true).is_empty());
    }
}
