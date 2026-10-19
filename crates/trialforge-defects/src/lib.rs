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
//! record receives at most one defect: once a defect is injected, its record, its paired record
//! and every record it expects a finding at are off limits to later ones.
//!
//! A request is either a map of defect id to [`Amount`] or one of the [`PROFILES`].

mod adam;
mod ct;
mod dates;
mod score;
mod seq;

use std::collections::{BTreeMap, HashMap, HashSet};

use serde::{Deserialize, Deserializer, Serialize, Serializer};
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
    DefectType {
        id: "orphan_usubjid",
        description: "an AE, CM or MH record's USUBJID changed (one digit) to a subject not in DM",
        rules: &["L1.SUBJ", "L1.EPOCH", "L1.SUPP", "L2.ADAE"],
    },
    DefectType {
        id: "record_after_death",
        description: "AEENDTC, or a DS record's DSSTDTC, moved to after DTHDTC (--DY follows)",
        rules: &["L1.DEATH", "L1.PART", "L1.EPOCH", "L2.ADAE"],
    },
    DefectType {
        id: "bad_study_day",
        description: "a --DY counted as if there were a day 0: one less after RFSTDTC, one more before",
        rules: &["L1.DY", "L0.DY"],
    },
    DefectType {
        id: "chg_mismatch",
        description: "a BDS CHG that is not AVAL - BASE (reversed, AVAL itself, or 1% off)",
        rules: &["L2.BDS"],
    },
    DefectType {
        id: "missing_baseline",
        description: "ABLFL cleared on a BDS baseline record",
        rules: &["L2.BDS"],
    },
    DefectType {
        id: "trtemfl_wrong",
        description: "ADAE TRTEMFL flipped",
        rules: &["L2.ADAE"],
    },
    DefectType {
        id: "cnsr_wrong",
        description: "ADTTE CNSR flipped",
        rules: &["L2.ADTTE", "L2.ADEXAC"],
    },
];

/// A named defect request.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Profile {
    /// Name, e.g. `realistic`.
    pub name: &'static str,
    /// What it is for.
    pub description: &'static str,
    /// Amount per defect id.
    pub defects: &'static [(&'static str, Amount)],
    /// Whether each type uses its subtlest variant: an AE end one day before its start, a
    /// record one day after death, a case-only CT difference, a CHG 1% off.
    pub subtle: bool,
}

/// The named profiles.
pub const PROFILES: &[Profile] = &[
    Profile {
        name: "light",
        description: "about five defects of the common kinds",
        defects: &[
            ("ae_end_before_start", Amount::Count(1)),
            ("duplicate_seq", Amount::Count(1)),
            ("out_of_ct", Amount::Count(2)),
            ("bad_study_day", Amount::Count(1)),
        ],
        subtle: false,
    },
    Profile {
        name: "realistic",
        description: "rates per eligible record: tens of defects of most kinds in 300 subjects",
        defects: &[
            ("ae_end_before_start", Amount::Rate(0.02)),
            ("duplicate_seq", Amount::Rate(0.001)),
            ("out_of_ct", Amount::Rate(0.0001)),
            ("orphan_usubjid", Amount::Rate(0.002)),
            ("record_after_death", Amount::Rate(0.02)),
            ("bad_study_day", Amount::Rate(0.0001)),
            ("chg_mismatch", Amount::Rate(0.0001)),
            ("missing_baseline", Amount::Rate(0.0005)),
            ("trtemfl_wrong", Amount::Rate(0.01)),
            ("cnsr_wrong", Amount::Rate(0.005)),
        ],
        subtle: false,
    },
    Profile {
        name: "adversarial",
        description: "three of every kind, each in its subtlest form",
        defects: &[
            ("ae_end_before_start", Amount::Count(3)),
            ("duplicate_seq", Amount::Count(3)),
            ("out_of_ct", Amount::Count(3)),
            ("orphan_usubjid", Amount::Count(3)),
            ("record_after_death", Amount::Count(3)),
            ("bad_study_day", Amount::Count(3)),
            ("chg_mismatch", Amount::Count(3)),
            ("missing_baseline", Amount::Count(3)),
            ("trtemfl_wrong", Amount::Count(3)),
            ("cnsr_wrong", Amount::Count(3)),
        ],
        subtle: true,
    },
];

/// The defect type called `id`.
pub fn defect_type(id: &str) -> Option<&'static DefectType> {
    CATALOGUE.iter().find(|t| t.id == id)
}

/// The profile called `name`.
pub fn profile(name: &str) -> Option<&'static Profile> {
    PROFILES.iter().find(|p| p.name == name)
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

/// Which defects to inject: defect id → [`Amount`], from a map or a [`Profile`]. Empty means
/// none.
///
/// Deserialises from `"none"`, a profile name, or a map such as
/// `{"out_of_ct": 5, "duplicate_seq": 0.02}`, where an integer is a count and any other number
/// a rate in `[0, 1]`. Serialises as the map of amounts.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DefectSpec {
    /// The profile the amounts come from, if any.
    pub profile: Option<&'static str>,
    /// Amount per defect id.
    pub defects: BTreeMap<String, Amount>,
}

impl Serialize for DefectSpec {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        self.defects.serialize(s)
    }
}

impl DefectSpec {
    /// Whether no defect type is requested.
    pub fn is_empty(&self) -> bool {
        self.defects.is_empty()
    }

    /// Whether defects take their subtlest form (see [`Profile::subtle`]).
    pub fn subtle(&self) -> bool {
        self.profile.and_then(profile).is_some_and(|p| p.subtle)
    }

    /// Parses and validates a spec from JSON (see the type docs).
    pub fn from_value(v: &Value) -> Result<Self, String> {
        let names: Vec<&str> = PROFILES.iter().map(|p| p.name).collect();
        let usage = format!(
            "pass \"none\", a profile ({}) or a map of defect id to count or rate",
            names.join(", ")
        );
        match v {
            Value::Null => Ok(Self::default()),
            Value::String(s) if s == "none" => Ok(Self::default()),
            Value::String(s) => match profile(s) {
                Some(p) => Ok(Self {
                    profile: Some(p.name),
                    defects: p.defects.iter().map(|&(id, a)| (id.into(), a)).collect(),
                }),
                None => Err(format!("unknown defect profile {s:?}; {usage}")),
            },
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
                Ok(Self {
                    profile: None,
                    defects,
                })
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
        subtle: spec.subtle(),
    };
    let mut defects = Vec::new();
    for t in CATALOGUE {
        let Some(&amount) = spec.defects.get(t.id) else {
            continue;
        };
        let s = stream.label(t.id);
        let new = match t.id {
            "ae_end_before_start" => dates::ae_end_before_start(&mut cx, s, amount),
            "duplicate_seq" => seq::duplicate_seq(&mut cx, s, amount),
            "out_of_ct" => ct::out_of_ct(&mut cx, s, amount),
            "orphan_usubjid" => seq::orphan_usubjid(&mut cx, s, amount),
            "record_after_death" => dates::record_after_death(&mut cx, s, amount),
            "bad_study_day" => dates::bad_study_day(&mut cx, s, amount),
            "chg_mismatch" => adam::chg_mismatch(&mut cx, s, amount),
            "missing_baseline" => adam::missing_baseline(&mut cx, s, amount),
            "trtemfl_wrong" => adam::trtemfl_wrong(&mut cx, s, amount),
            "cnsr_wrong" => adam::cnsr_wrong(&mut cx, s, amount),
            _ => unreachable!("every catalogue entry is implemented"),
        };
        for d in &new {
            cx.used.insert((d.dataset.clone(), d.row));
            if let Some(j) = d.paired_row {
                cx.used.insert((d.dataset.clone(), j));
            }
            for e in &d.expected {
                if let Some(r) = e.row {
                    cx.used.insert((e.dataset.clone(), r));
                }
            }
        }
        defects.extend(new);
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
    /// (dataset, row) of every record that already has a defect or an expected finding.
    pub used: HashSet<(String, usize)>,
    /// Whether to inject the subtlest variant of each type.
    pub subtle: bool,
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

    /// Applies `edits` to the ADaM dataset `name`.
    pub fn edit_adam(&mut self, name: &str, edits: &[Edit]) {
        if let Some(ds) = self.adam.iter_mut().find(|d| d.name == name) {
            apply(ds, edits);
        }
    }

    /// Reference dates per subject, from DM.
    pub fn subjects(&self) -> HashMap<String, Subject> {
        let Some(dm) = self.sdtm("DM") else {
            return HashMap::new();
        };
        let Some(u) = chars(dm, "USUBJID") else {
            return HashMap::new();
        };
        let col = |v: &str| -> Vec<Option<Date>> {
            match chars(dm, v) {
                Some(x) => x
                    .iter()
                    .map(|x| x.as_deref().and_then(parse_date))
                    .collect(),
                None => vec![None; dm.n_rows()],
            }
        };
        let (rfic, rfst, rfpen, death) = (
            col("RFICDTC"),
            col("RFSTDTC"),
            col("RFPENDTC"),
            col("DTHDTC"),
        );
        (0..dm.n_rows())
            .filter_map(|i| {
                Some((
                    u[i].clone()?,
                    Subject {
                        rfic: rfic[i],
                        rfst: rfst[i],
                        rfpen: rfpen[i],
                        death: death[i],
                    },
                ))
            })
            .collect()
    }

    /// The ADAE row of the AE record (`usubjid`, `aeseq`).
    pub fn adae_row(&self, usubjid: &str, aeseq: i64) -> Option<usize> {
        let ad = self.adam("ADAE")?;
        let (u, s) = (chars(ad, "USUBJID")?, ints(ad, "AESEQ")?);
        (0..ad.n_rows()).find(|&j| u[j].as_deref() == Some(usubjid) && s[j] == Some(aeseq))
    }

    /// Whether row `row` of the SDTM dataset `name`, dated `date` and belonging to `usubjid`,
    /// lies outside every SE element of its EPOCH (what `L1.EPOCH` flags).
    pub fn outside_epoch(&self, name: &str, row: usize, usubjid: &str, date: Date) -> bool {
        let (Some(se), Some(ds)) = (self.sdtm("SE"), self.sdtm(name)) else {
            return false;
        };
        let (Some(su), Some(sep), Some(sst), Some(sen)) = (
            chars(se, "USUBJID"),
            chars(se, "EPOCH"),
            chars(se, "SESTDTC"),
            chars(se, "SEENDTC"),
        ) else {
            return false;
        };
        let Some(epoch) = chars(ds, "EPOCH").and_then(|e| e[row].as_deref()) else {
            return false;
        };
        let date_of = |v: &Option<String>| v.as_deref().and_then(parse_date);
        !(0..se.n_rows()).any(|k| {
            su[k].as_deref() == Some(usubjid)
                && sep[k].as_deref() == Some(epoch)
                && matches!((date_of(&sst[k]), date_of(&sen[k])), (Some(a), Some(b)) if a <= date && date <= b)
        })
    }
}

/// A subject's reference dates.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Subject {
    pub rfic: Option<Date>,
    pub rfst: Option<Date>,
    pub rfpen: Option<Date>,
    pub death: Option<Date>,
}

/// A new value for one cell.
pub(crate) enum Cell {
    Char(Option<String>),
    Int(Option<i64>),
    Num(Option<f64>),
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
            (Values::Num(v), Cell::Num(x)) => v[e.row] = *x,
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

/// Numeric (`Num`) values of `var`.
pub(crate) fn nums<'a>(ds: &'a Dataset, var: &str) -> Option<&'a [Option<f64>]> {
    match &ds.column(var)?.values {
        Values::Num(v) => Some(v),
        _ => None,
    }
}

/// `{var}={value}` for a record's --SEQ, if it has one.
pub(crate) fn seq_key(ds: &Dataset, row: usize) -> Option<String> {
    let var = format!("{}SEQ", ds.name);
    ints(ds, &var)?[row].map(|s| format!("{var}={s}"))
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
            json!("bogus"),
            json!([1]),
        ] {
            assert!(DefectSpec::from_value(&bad).is_err(), "{bad}");
        }
        let e = DefectSpec::from_value(&json!("heavy")).unwrap_err();
        assert!(e.contains("light, realistic, adversarial"), "{e}");
    }

    #[test]
    fn profiles_expand_to_catalogue_types() {
        for p in PROFILES {
            let s = DefectSpec::from_value(&json!(p.name)).unwrap();
            assert_eq!(s.profile, Some(p.name));
            assert_eq!(s.subtle(), p.subtle);
            assert_eq!(s.defects.len(), p.defects.len());
            for (id, a) in p.defects {
                assert!(defect_type(id).is_some(), "{id}");
                assert!(matches!(a, Amount::Count(_) | Amount::Rate(0.0..=1.0)));
            }
            // Serialised as its amounts.
            assert!(serde_json::to_value(&s).unwrap().is_object());
        }
        let all = profile("adversarial").unwrap();
        assert_eq!(all.defects.len(), CATALOGUE.len());
        assert!(
            !DefectSpec::from_value(&json!({"out_of_ct": 1}))
                .unwrap()
                .subtle()
        );
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
