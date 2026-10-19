//! Defect injection: every injected defect is caught by the rule its manifest names, and
//! nothing else fails (recall and precision 1 at L0–L2).

use serde_json::json;
use trialforge::defects::{DefectSpec, Observed};
use trialforge::{ARCHETYPES, Study, StudyOptions, study};

fn opts(n: usize, seed: u64, defects: serde_json::Value) -> StudyOptions {
    StudyOptions {
        n,
        seed,
        studyid: None,
        controls: json!(null),
        defects: Some(DefectSpec::from_value(&defects).unwrap()),
    }
}

fn observed(r: &trialforge::check::Report) -> Vec<(String, &'static str, String, Option<usize>)> {
    r.findings
        .iter()
        .map(|f| (format!("{:?}", f.level), f.rule, f.dataset.clone(), f.row))
        .collect()
}

/// Asserts every expected finding is reported, and every finding is expected.
fn assert_exact(s: &Study) {
    let report = trialforge::verify(s);
    let found = observed(&report);
    let m = s.defects.as_ref().expect("defects requested");
    for d in &m.defects {
        for e in &d.expected {
            assert!(
                found
                    .iter()
                    .any(|f| f.1 == e.rule && f.2 == e.dataset && f.3 == e.row),
                "{}: {} expected {e:?} for {d:?}\n{}",
                s.archetype,
                d.id,
                report.summary()
            );
        }
    }
    let obs: Vec<Observed<'_>> = found
        .iter()
        .map(|f| Observed {
            level: &f.0,
            rule: f.1,
            dataset: &f.2,
            row: f.3,
        })
        .collect();
    let score = s.score(&obs).unwrap();
    let unmatched: Vec<_> = score
        .unmatched
        .iter()
        .map(|&i| &report.findings[i])
        .collect();
    assert!(
        unmatched.is_empty(),
        "{} (seed {}): findings no defect explains: {unmatched:#?}",
        s.archetype,
        s.seed
    );
    assert_eq!(score.recall, (!m.defects.is_empty()).then_some(1.0));
}

#[test]
fn every_archetype_catches_exactly_its_defects() {
    let spec = json!({"ae_end_before_start": 3, "duplicate_seq": 3, "out_of_ct": 5});
    for (archetype, _) in ARCHETYPES {
        for seed in 1..=3 {
            let s = study(archetype, &opts(40, seed, spec.clone())).unwrap();
            assert!(
                !s.defects.as_ref().unwrap().defects.is_empty(),
                "{archetype}"
            );
            assert_exact(&s);
        }
    }
}

#[test]
fn rates_and_counts() {
    let s = study(
        "parallel_continuous",
        &opts(
            60,
            7,
            json!({"ae_end_before_start": 1.0, "out_of_ct": 0.05}),
        ),
    )
    .unwrap();
    let m = s.defects.as_ref().unwrap();
    let ae = s.dataset("AE").unwrap();
    let n_ae = m
        .defects
        .iter()
        .filter(|d| d.id == "ae_end_before_start")
        .count();
    // Every AE with an end date after the day of consent.
    assert!(n_ae > 0 && n_ae <= ae.n_rows());
    assert!(m.defects.iter().any(|d| d.id == "out_of_ct"));
    assert_exact(&s);

    let s = study("parallel_continuous", &opts(60, 7, json!({"out_of_ct": 4}))).unwrap();
    assert_eq!(s.defects.as_ref().unwrap().defects.len(), 4);
    assert_exact(&s);
}

#[test]
fn defects_are_deterministic_and_independent() {
    let spec = json!({"ae_end_before_start": 4, "out_of_ct": 4});
    let a = study("parallel_continuous", &opts(80, 3, spec.clone())).unwrap();
    let b = study("parallel_continuous", &opts(80, 3, spec)).unwrap();
    assert_eq!(a, b);
    // A type applied earlier does not move when a later one is added.
    let alone = study(
        "parallel_continuous",
        &opts(80, 3, json!({"ae_end_before_start": 4})),
    )
    .unwrap();
    let first = |s: &Study| -> Vec<_> {
        s.defects
            .as_ref()
            .unwrap()
            .defects
            .iter()
            .filter(|d| d.id == "ae_end_before_start")
            .cloned()
            .collect()
    };
    assert_eq!(first(&a), first(&alone));
    // The truth manifest and every untouched dataset are the clean study's.
    let clean = study(
        "parallel_continuous",
        &StudyOptions {
            defects: None,
            ..opts(80, 3, json!("none"))
        },
    )
    .unwrap();
    assert!(clean.defects.is_none());
    assert_eq!(a.truth, clean.truth);
    let touched: Vec<&str> = a
        .defects
        .as_ref()
        .unwrap()
        .defects
        .iter()
        .map(|d| d.dataset.as_str())
        .collect();
    for (x, y) in a
        .sdtm
        .iter()
        .chain(&a.adam)
        .zip(clean.sdtm.iter().chain(&clean.adam))
    {
        assert_eq!(touched.contains(&x.name.as_str()), x != y, "{}", x.name);
    }
    // An empty request is a clean study.
    let none = study("parallel_continuous", &opts(80, 3, json!("none"))).unwrap();
    assert_eq!(none, clean);
}

#[test]
fn manifest_table_and_files() {
    let s = study(
        "parallel_continuous",
        &opts(40, 2, json!({"duplicate_seq": 2, "out_of_ct": 2})),
    )
    .unwrap();
    let t = s.defects_table(0).unwrap();
    assert_eq!(t.n_rows(), 4);
    let rows = |t: &trialforge::Dataset| t.column("row").unwrap().as_int().unwrap().to_vec();
    let one = s.defects_table(1).unwrap();
    assert!(
        rows(&t)
            .iter()
            .zip(rows(&one))
            .all(|(a, b)| a.unwrap() + 1 == b.unwrap())
    );
    let j = s.defects_json().unwrap();
    assert_eq!(j["spec"], json!({"duplicate_seq": 2, "out_of_ct": 2}));
    assert_eq!(j["defects"].as_array().unwrap().len(), 4);
    let dir = std::env::temp_dir().join(format!("tf-defects-{}", std::process::id()));
    let written = s.write(&dir, &["csv"]).unwrap();
    assert!(written.iter().any(|p| p.ends_with("defects.json")));
    std::fs::remove_dir_all(&dir).unwrap();

    let clean = study("parallel_continuous", &opts(40, 2, json!({}))).unwrap();
    assert!(clean.defects_table(0).is_none());
    assert!(clean.score(&[]).is_err());
}
