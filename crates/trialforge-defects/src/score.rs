//! Scoring a validator's findings against the injected defects.

use std::collections::{BTreeMap, HashMap};

use serde::Serialize;

use crate::{Defect, defect_type};

/// One finding of the validator being scored.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Observed<'a> {
    /// Level: `L0`, `L1` and `L2` are `trialforge-check` levels and must match the expected
    /// rule; any other level (Pointblank's `L3`, say) is matched by dataset and row only, and
    /// may also name a defect's [`paired_row`](Defect::paired_row).
    pub level: &'a str,
    /// Rule identifier.
    pub rule: &'a str,
    /// Dataset (matched case-insensitively).
    pub dataset: &'a str,
    /// Zero-based row.
    pub row: Option<usize>,
}

/// Precision and recall for one defect type.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct TypeScore {
    /// Defects of this type injected.
    pub injected: usize,
    /// Of those, how many at least one finding matched.
    pub detected: usize,
    /// `detected / injected`.
    pub recall: Option<f64>,
    /// Findings matched to defects of this type.
    pub findings: usize,
    /// Unmatched `trialforge-check` findings from a rule that catches this type.
    pub false_positives: usize,
    /// `findings / (findings + false_positives)`.
    pub precision: Option<f64>,
}

/// How well a validator's findings recover the injected defects.
#[derive(Clone, Debug, PartialEq, Serialize)]
pub struct Score {
    /// Defects injected.
    pub injected: usize,
    /// Defects at least one finding matched.
    pub detected: usize,
    /// `detected / injected`.
    pub recall: Option<f64>,
    /// Findings scored.
    pub findings: usize,
    /// Findings that matched a defect.
    pub matched: usize,
    /// `matched / findings`.
    pub precision: Option<f64>,
    /// Per defect type.
    pub by_type: BTreeMap<&'static str, TypeScore>,
    /// Indices (into the defects) of the defects no finding matched.
    pub missed: Vec<usize>,
    /// Indices (into the findings) of the findings that matched no defect.
    pub unmatched: Vec<usize>,
}

/// (upper-case dataset, row) → the defects expecting a finding there, with the rule (`None`
/// for a paired row, which only other validators' findings match).
type Index<'a> = HashMap<(String, Option<usize>), Vec<(usize, Option<&'a str>)>>;

fn ratio(a: usize, b: usize) -> Option<f64> {
    (b > 0).then(|| a as f64 / b as f64)
}

/// Scores `findings` against `defects`. A finding matches a defect when it names one of the
/// defect's [`Expected`](crate::Expected) datasets and rows (and, for a `trialforge-check`
/// level, its rule), or for another validator its paired row; a defect is detected when any
/// finding matches it.
pub fn score(defects: &[Defect], findings: &[Observed<'_>]) -> Score {
    let mut index = Index::new();
    for (k, d) in defects.iter().enumerate() {
        for e in &d.expected {
            index
                .entry((e.dataset.to_ascii_uppercase(), e.row))
                .or_default()
                .push((k, Some(e.rule)));
        }
        if let Some(row) = d.paired_row {
            index
                .entry((d.dataset.to_ascii_uppercase(), Some(row)))
                .or_default()
                .push((k, None));
        }
    }
    let mut detected = vec![false; defects.len()];
    let mut by_type: BTreeMap<&'static str, TypeScore> = BTreeMap::new();
    for d in defects {
        let t = by_type.entry(d.id).or_insert(TypeScore {
            injected: 0,
            detected: 0,
            recall: None,
            findings: 0,
            false_positives: 0,
            precision: None,
        });
        t.injected += 1;
    }
    let mut unmatched = Vec::new();
    for (fi, f) in findings.iter().enumerate() {
        let internal = matches!(f.level, "L0" | "L1" | "L2");
        let hits: Vec<usize> = index
            .get(&(f.dataset.to_ascii_uppercase(), f.row))
            .into_iter()
            .flatten()
            .filter(|(_, rule)| !internal || *rule == Some(f.rule))
            .map(|(k, _)| *k)
            .collect();
        match hits.first() {
            Some(&first) => {
                for &k in &hits {
                    detected[k] = true;
                }
                by_type
                    .get_mut(defects[first].id)
                    .expect("every defect's type is scored")
                    .findings += 1;
            }
            None => {
                unmatched.push(fi);
                if internal {
                    for (id, t) in &mut by_type {
                        if defect_type(id).is_some_and(|dt| dt.rules.contains(&f.rule)) {
                            t.false_positives += 1;
                        }
                    }
                }
            }
        }
    }
    for (k, d) in defects.iter().enumerate() {
        if detected[k] {
            by_type.get_mut(d.id).expect("scored").detected += 1;
        }
    }
    for t in by_type.values_mut() {
        t.recall = ratio(t.detected, t.injected);
        t.precision = ratio(t.findings, t.findings + t.false_positives);
    }
    let n_detected = detected.iter().filter(|x| **x).count();
    let matched = findings.len() - unmatched.len();
    Score {
        injected: defects.len(),
        detected: n_detected,
        recall: ratio(n_detected, defects.len()),
        findings: findings.len(),
        matched,
        precision: ratio(matched, findings.len()),
        by_type,
        missed: (0..defects.len()).filter(|&k| !detected[k]).collect(),
        unmatched,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Expected;

    fn defect(id: &'static str, dataset: &str, row: usize, rule: &'static str) -> Defect {
        Defect {
            id,
            dataset: dataset.into(),
            row,
            usubjid: None,
            key: None,
            variable: String::new(),
            original: None,
            injected: None,
            paired_row: None,
            expected: vec![Expected {
                rule,
                dataset: dataset.into(),
                row: Some(row),
            }],
        }
    }

    #[test]
    fn matches_by_rule_or_by_row() {
        let defects = [
            defect("out_of_ct", "AE", 3, "L0.CT"),
            defect("out_of_ct", "DM", 0, "L0.CT"),
            Defect {
                paired_row: Some(6),
                ..defect("duplicate_seq", "CM", 7, "L0.SEQ")
            },
        ];
        let f = |level, rule, dataset, row| Observed {
            level,
            rule,
            dataset,
            row: Some(row),
        };
        let findings = [
            f("L0", "L0.CT", "AE", 3),       // matches defect 0
            f("L0", "L0.KEY", "CM", 7),      // right row, wrong rule
            f("L3", "col_vals_in", "cm", 7), // L3: dataset and row suffice
            f("L0", "L0.CT", "LB", 9),       // a false positive for out_of_ct
            f("L3", "SDTM-149", "CM", 6),    // L3 may name the paired row
            f("L0", "L0.SEQ", "CM", 6),      // trialforge's own rules may not
        ];
        let s = score(&defects, &findings);
        assert_eq!((s.injected, s.detected, s.matched), (3, 2, 3));
        assert_eq!(s.missed, vec![1]);
        assert_eq!(s.unmatched, vec![1, 3, 5]);
        assert_eq!(s.precision, Some(0.5));
        let ct = &s.by_type["out_of_ct"];
        assert_eq!(
            (ct.injected, ct.detected, ct.findings, ct.false_positives),
            (2, 1, 1, 1)
        );
        assert_eq!(ct.recall, Some(0.5));
        assert_eq!(ct.precision, Some(0.5));
        let seq = &s.by_type["duplicate_seq"];
        // The L0.KEY and paired-row L0.SEQ findings are duplicate_seq rules that matched
        // nothing.
        assert_eq!((seq.recall, seq.precision), (Some(1.0), Some(0.5)));
        assert_eq!(score(&[], &[]).recall, None);
    }
}
