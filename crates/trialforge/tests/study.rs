//! End-to-end tests. Every successful `study()` has already passed verification levels 0–2
//! (generation fails otherwise), so the property tests only need to ask for random studies.
//!
//! The Phase 1 exit criterion runs the property tests with 1,000 cases:
//! `PROPTEST_CASES=1000 cargo test --release -p trialforge --test study`.

use proptest::prelude::*;
use serde_json::json;
use trialforge::{Dataset, StudyOptions, Values, study};

fn opts(n: usize, seed: u64, controls: serde_json::Value) -> StudyOptions {
    StudyOptions {
        n,
        seed,
        studyid: None,
        controls,
    }
}

fn xpt(s: &trialforge::Study) -> Vec<Vec<u8>> {
    s.sdtm
        .iter()
        .chain(&s.adam)
        .map(|d| trialforge::io::xpt::to_bytes(d, s.cutoff).unwrap())
        .collect()
}

fn chars<'a>(d: &'a Dataset, v: &str) -> &'a [Option<String>] {
    d.column(v)
        .and_then(|c| c.as_char())
        .unwrap_or_else(|| panic!("{}.{v}", d.name))
}

#[test]
fn generation_is_deterministic() {
    for archetype in ["parallel_continuous", "oncology_solid_tumor"] {
        let a = study(archetype, &opts(120, 1, json!(null))).unwrap();
        let b = study(archetype, &opts(120, 1, json!(null))).unwrap();
        let c = study(archetype, &opts(120, 2, json!(null))).unwrap();
        assert_eq!(a, b, "{archetype}");
        assert_ne!(a.sdtm, c.sdtm, "{archetype}");
        assert_eq!(xpt(&a), xpt(&b), "{archetype}");
    }
}

#[test]
fn parallel_study_has_every_phase1_domain() {
    let s = study("parallel_continuous", &opts(60, 3, json!(null))).unwrap();
    let names: Vec<&str> = s
        .sdtm
        .iter()
        .chain(&s.adam)
        .map(|d| d.name.as_str())
        .collect();
    for want in [
        "TA", "TE", "TV", "TI", "TS", "DM", "SE", "SV", "EX", "CM", "AE", "DS", "MH", "LB", "VS",
        "SUPPAE", "ADSL", "ADAE", "ADLB", "ADVS",
    ] {
        assert!(names.contains(&want), "missing {want}: {names:?}");
    }
}

#[test]
fn endpoint_effect_is_recovered_at_scale() {
    for (indication, tol) in [
        ("type_2_diabetes", 0.08),
        ("hypertension", 1.5),
        ("obesity", 1.0),
    ] {
        let s = study(
            "parallel_continuous",
            &opts(1000, 7, json!({"indication": indication})),
        )
        .unwrap();
        let e = &s.truth["endpoint"];
        let (target, realised) = (
            e["target_difference"].as_f64().unwrap(),
            e["realised_difference"].as_f64().unwrap(),
        );
        assert!(
            (target - realised).abs() < tol,
            "{indication}: target {target}, realised {realised}"
        );
    }
}

#[test]
fn active_drug_raises_its_signature_adverse_events() {
    let s = study("parallel_continuous", &opts(1000, 11, json!(null))).unwrap();
    let adae = s.dataset("ADAE").unwrap();
    let (trt, pt, fl) = (
        chars(adae, "TRTA"),
        chars(adae, "AEDECOD"),
        chars(adae, "TRTEMFL"),
    );
    let count = |arm: &str| {
        (0..adae.n_rows())
            .filter(|&i| {
                trt[i].as_deref() == Some(arm)
                    && pt[i].as_deref() == Some("Nausea")
                    && fl[i].as_deref() == Some("Y")
            })
            .count()
    };
    let (active, placebo) = (count("Drug A 10 mg"), count("Placebo"));
    assert!(
        active as f64 > 1.8 * placebo as f64,
        "nausea: active {active}, placebo {placebo}"
    );
}

#[test]
fn baseline_inclusion_criterion_holds() {
    let s = study("parallel_continuous", &opts(300, 5, json!(null))).unwrap();
    let lb = s.dataset("LB").unwrap();
    let (tc, visit) = (chars(lb, "LBTESTCD"), chars(lb, "VISIT"));
    let Values::Num(v) = &lb.column("LBSTRESN").unwrap().values else {
        panic!()
    };
    for i in 0..lb.n_rows() {
        if tc[i].as_deref() == Some("HBA1C") && visit[i].as_deref() == Some("SCREENING") {
            let x = v[i].unwrap();
            assert!((7.0..=10.5).contains(&x), "screening HbA1c {x}");
        }
    }
}

#[test]
fn oncology_truth_reports_overall_survival() {
    let s = study("oncology_solid_tumor", &opts(1000, 4, json!(null))).unwrap();
    let os = &s.truth["overall_survival"];
    let a = os["by_arm"]["Drug A"]["km_median_months"].as_f64().unwrap();
    let c = os["by_arm"]["Control"]["km_median_months"]
        .as_f64()
        .unwrap();
    assert!(a > c, "Drug A should live longer: {a} vs {c}");
    let ci = os["stratified_cox_hr"]["ci95"].as_array().unwrap();
    let (lo, hi) = (ci[0].as_f64().unwrap(), ci[1].as_f64().unwrap());
    assert!(
        lo < 0.75 && 0.75 < hi,
        "OS HR CI [{lo}, {hi}] should cover the target 0.75"
    );
    let orr = &s.truth["confirmed_orr"];
    assert!(
        orr["Drug A"]["proportion"].as_f64().unwrap()
            > orr["Control"]["proportion"].as_f64().unwrap()
    );
}

/// Rewrites column `var` of `dataset` with `f`.
fn edit(sets: &mut [Dataset], dataset: &str, var: &str, f: impl FnOnce(&mut Values)) {
    let d = sets.iter_mut().find(|d| d.name == dataset).unwrap();
    let mut columns = d.columns().to_vec();
    let c = columns.iter_mut().find(|c| c.name == var).unwrap();
    f(&mut c.values);
    *d = Dataset::new(d.name.clone(), d.label.clone(), columns).unwrap();
}

/// Replaces the first `from` in character column `var` with `to`.
fn corrupt(sets: &mut [Dataset], dataset: &str, var: &str, from: &str, to: &str) {
    edit(sets, dataset, var, |v| {
        let Values::Char(v) = v else {
            panic!("{var} is character")
        };
        let x = v.iter_mut().find(|x| x.as_deref() == Some(from)).unwrap();
        *x = Some(to.to_string());
    });
}

#[test]
fn oncology_checks_catch_inconsistent_responses() {
    let base = study("oncology_solid_tumor", &opts(60, 2, json!(null))).unwrap();
    assert!(trialforge::verify(&base).passed());
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };

    // An RS response that RECIST 1.1 does not give.
    let mut s = base.clone();
    corrupt(&mut s.sdtm, "RS", "RSSTRESC", "SD", "PR");
    assert!(rules(&s).contains(&"L1.RECIST"), "{:?}", rules(&s));

    // A best overall response that the responses do not support.
    let mut s = base.clone();
    corrupt(&mut s.adam, "ADRS", "AVALC", "MISSING", "SD");
    assert!(rules(&s).contains(&"L2.ADRS"), "{:?}", rules(&s));

    // A survival record flipped from event to censored.
    let mut s = base.clone();
    edit(&mut s.adam, "ADTTE", "CNSR", |v| {
        let Values::Int(v) = v else {
            panic!("CNSR is integer")
        };
        let x = v.iter_mut().find(|x| **x == Some(0)).unwrap();
        *x = Some(1);
    });
    assert!(rules(&s).contains(&"L2.ADTTE"), "{:?}", rules(&s));
}

#[test]
fn bad_input_is_rejected() {
    use trialforge::Error::Invalid;
    assert!(matches!(
        study("nope", &opts(10, 1, json!(null))),
        Err(Invalid(_))
    ));
    assert!(matches!(
        study("parallel_continuous", &opts(10, 1, json!({"weeks": 1}))),
        Err(Invalid(_))
    ));
    assert!(matches!(
        study(
            "parallel_continuous",
            &opts(10, 1, json!({"colour": "red"}))
        ),
        Err(Invalid(_))
    ));
    assert!(matches!(
        study("parallel_continuous", &opts(0, 1, json!(null))),
        Err(Invalid(_))
    ));
    assert!(matches!(
        study(
            "parallel_continuous",
            &opts(10, trialforge::MAX_SEED + 1, json!(null))
        ),
        Err(Invalid(_))
    ));
    let s = study("parallel_continuous", &opts(5, 1, json!(null))).unwrap();
    assert!(matches!(
        s.write(&std::env::temp_dir().join("tf-bad-format"), &["sas7bdat"]),
        Err(Invalid(_))
    ));
}

#[test]
fn controls_are_documented() {
    for (name, _) in trialforge::ARCHETYPES {
        let docs = trialforge::controls(name).unwrap();
        for d in docs.as_array().unwrap() {
            assert!(!d["description"].as_str().unwrap().is_empty());
            assert!(!d["default"].is_null() || d["name"] == "change");
        }
    }
}

prop_compose! {
    fn parallel_controls()(
        indication in prop::sample::select(vec!["type_2_diabetes", "hypertension", "obesity"]),
        weeks in 4u32..=52,
        change in prop::option::of((-5.0f64..5.0, -15.0f64..5.0)),
        ratio in (1u32..=3, 1u32..=3),
        dose in 0.5f64..500.0,
        dropout in 0.0f64..2.0,
        ae in 0.0f64..6.0,
        missed in 0.0f64..0.5,
    ) -> serde_json::Value {
        json!({
            "indication": indication, "weeks": weeks, "change": change.map(|(a, b)| [a, b]),
            "ratio": [ratio.0, ratio.1], "dose": dose, "dropout_rate": dropout,
            "ae_rate_multiplier": ae, "missed_visit": missed,
        })
    }
}

prop_compose! {
    fn oncology_controls()(
        os_median in 3.0f64..40.0,
        pfs_fraction in 0.1f64..0.9,
        hr in (0.3f64..1.5, 0.3f64..1.5),
        orr in (0.0f64..0.9, 0.0f64..0.9),
        shape in (0.5f64..3.0, 0.5f64..3.0),
        correlation in 0.0f64..0.95,
        withdrawal in 0.0f64..1.0,
        ae in 0.0f64..1.0,
        enrollment in 1.0f64..36.0,
        extra in 1.0f64..36.0,
    ) -> serde_json::Value {
        json!({
            "os_median_months": os_median, "pfs_median_months": os_median * pfs_fraction,
            "os_hr": hr.0, "pfs_hr": hr.1, "orr": [orr.0, orr.1],
            "os_shape": shape.0, "pfs_shape": shape.1, "os_pfs_correlation": correlation,
            "withdrawal_prob": withdrawal, "ae_stop_prob": ae,
            "enrollment_months": enrollment, "cutoff_months": enrollment + extra,
        })
    }
}

proptest! {
    #![proptest_config(ProptestConfig { cases: 48, ..ProptestConfig::default() })]

    #[test]
    fn random_parallel_specs_verify(c in parallel_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("parallel_continuous", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_oncology_specs_verify(c in oncology_controls(), n in 1usize..80, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("oncology_solid_tumor", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }
}
