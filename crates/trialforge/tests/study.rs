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
    for archetype in [
        "parallel_continuous",
        "oncology_solid_tumor",
        "phase1_sad_mad",
    ] {
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
            let nullable = ["change", "ni_margin", "placebo_change", "emax"];
            assert!(!d["default"].is_null() || nullable.iter().any(|n| d["name"] == *n));
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
        active in any::<bool>(),
        margin in prop::option::of(0.01f64..5.0),
    ) -> serde_json::Value {
        json!({
            "comparator": if active { "active" } else { "placebo" }, "ni_margin": margin,
            "indication": indication, "weeks": weeks, "change": change.map(|(a, b)| [a, b]),
            "ratio": [ratio.0, ratio.1], "dose": dose, "dropout_rate": dropout,
            "ae_rate_multiplier": ae, "missed_visit": missed,
        })
    }
}

#[test]
fn phase1_pk_is_verified_against_an_independent_nca() {
    let base = study("phase1_sad_mad", &opts(16, 2, json!({"doses": [10, 100]}))).unwrap();
    assert!(trialforge::verify(&base).passed());
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let pp = base.sdtm.iter().find(|d| d.name == "PP").unwrap();
    let cd = chars(pp, "PPTESTCD");
    assert!(cd.iter().any(|c| c.as_deref() == Some("AUCIFO")));
    // A wrong PP value.
    let mut s = base.clone();
    edit(&mut s.sdtm, "PP", "PPSTRESN", |v| {
        let Values::Num(v) = v else {
            panic!("PPSTRESN is numeric")
        };
        v[0] = v[0].map(|x| x * 1.01);
    });
    assert!(rules(&s).contains(&"L2.PP"), "{:?}", rules(&s));
    // A concentration reported both as BLQ and as a number.
    let mut s = base.clone();
    corrupt(&mut s.sdtm, "PC", "PCSTRESC", "<BLQ", "1.5");
    assert!(rules(&s).contains(&"L1.PK"), "{:?}", rules(&s));
    // The multiple-dose part has steady-state parameters and accumulation ratios.
    let mad = study(
        "phase1_sad_mad",
        &opts(16, 2, json!({"part": "mad", "doses": [50, 100]})),
    )
    .unwrap();
    let pp = mad.sdtm.iter().find(|d| d.name == "PP").unwrap();
    let cd = chars(pp, "PPTESTCD");
    for want in ["AUCTAU", "CLFTAU", "ARAUC"] {
        assert!(cd.iter().any(|c| c.as_deref() == Some(want)), "{want}");
    }
    let slope = &base.truth["pk"]["power_model"];
    assert!(slope.get("CMAX").is_some());
}

#[test]
fn crossover_periods_are_verified_and_be_is_recovered() {
    let base = study("crossover_be", &opts(48, 2, json!({"design": "2x2x4"}))).unwrap();
    assert!(trialforge::verify(&base).passed());
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let ta = base.sdtm.iter().find(|d| d.name == "TA").unwrap();
    assert!(
        chars(ta, "EPOCH")
            .iter()
            .any(|e| e.as_deref() == Some("TREATMENT 4"))
    );
    // A period treatment swapped in ADSL.
    let mut s = base.clone();
    edit(&mut s.adam, "ADSL", "TRT02A", |v| {
        let Values::Char(v) = v else {
            panic!("TRT02A is character")
        };
        let other = v
            .iter()
            .flatten()
            .find(|x| Some(*x) != v[0].as_ref())
            .cloned();
        v[0] = other;
    });
    assert!(rules(&s).contains(&"L2.PERIOD"), "{:?}", rules(&s));
    // A large study recovers the GMR within its CI.
    let big = study("crossover_be", &opts(240, 3, json!({"gmr": 0.9}))).unwrap();
    let auc = &big.truth["bioequivalence"]["AUCIFO"];
    let ci = auc["ci90"].as_array().unwrap();
    let (lo, hi) = (ci[0].as_f64().unwrap(), ci[1].as_f64().unwrap());
    assert!(lo < 0.9 && 0.9 < hi, "GMR CI [{lo}, {hi}] should cover 0.9");
}

#[test]
fn dose_response_and_non_inferiority_are_recovered() {
    let s = study("dose_ranging_ph2", &opts(600, 1, json!(null))).unwrap();
    let fit = &s.truth["dose_response"]["fitted"];
    let ed50 = fit["ed50"].as_f64().unwrap();
    assert!(
        (3.0..20.0).contains(&ed50),
        "ED50 {ed50} should be near 7.5"
    );
    let ni = study(
        "parallel_continuous",
        &opts(300, 1, json!({"comparator": "active", "ni_margin": 0.3})),
    )
    .unwrap();
    assert_eq!(ni.truth["non_inferiority"]["non_inferior"], json!(true));
    let ta = ni.sdtm.iter().find(|d| d.name == "TS").unwrap();
    let vals = chars(ta, "TSVAL");
    assert!(vals.iter().any(|v| v.as_deref() == Some("ACTIVE")));
}

#[test]
fn adas_cog_totals_are_verified() {
    let base = study("cns_alzheimers", &opts(60, 3, json!({"item_missing": 0.1}))).unwrap();
    assert!(trialforge::verify(&base).passed());
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut s = base.clone();
    edit(&mut s.adam, "ADADAS", "AVAL", |v| {
        let Values::Num(v) = v else {
            panic!("AVAL is numeric")
        };
        let x = v.iter_mut().flatten().last().unwrap();
        *x += 1.0;
    });
    assert!(rules(&s).contains(&"L2.ADADAS"), "{:?}", rules(&s));
    let ft = base.sdtm.iter().find(|d| d.name == "FT").unwrap();
    assert!(
        chars(ft, "FTSTAT")
            .iter()
            .any(|x| x.as_deref() == Some("NOT DONE"))
    );
}

#[test]
fn cv_outcomes_are_event_driven_and_verified() {
    let base = study("cv_outcomes", &opts(800, 2, json!({"target_events": 60}))).unwrap();
    assert!(trialforge::verify(&base).passed());
    let events = base.truth["outcomes"]["realised_first_mace_events"]
        .as_u64()
        .unwrap();
    assert!(
        (60..=75).contains(&events),
        "{events} first MACE events for a target of 60"
    );
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
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

prop_compose! {
    fn cv_controls()(
        target in 1usize..200,
        rate in 0.005f64..0.5,
        hr in 0.3f64..2.0,
        components in (0.0f64..1.0, 0.0f64..1.0, 0.01f64..1.0),
        ncv in 0.0f64..0.2,
        frailty in 0.0f64..1.5,
        enrollment in 1.0f64..36.0,
        extra in 1.0f64..60.0,
        stop in 0.0f64..1.0,
        withdraw in 0.0f64..1.0,
    ) -> serde_json::Value {
        json!({
            "target_events": target, "mace_rate": rate, "hr": hr,
            "components": [components.0, components.1, components.2], "non_cv_death_rate": ncv,
            "frailty_sd": frailty, "enrollment_months": enrollment, "max_months": enrollment + extra,
            "discontinuation_rate": stop, "withdrawal_rate": withdraw,
        })
    }
}

prop_compose! {
    fn cns_controls()(
        weeks in prop::sample::select(vec![8u32, 16, 24, 48]),
        doses in prop::collection::vec(1.0f64..100.0, 1..4),
        baseline in 10.0f64..45.0,
        missing in 0.0f64..0.3,
        decline in 0.0f64..10.0,
        dropout in 0.0f64..3.0,
        seed_change in -5.0f64..8.0,
    ) -> serde_json::Value {
        let change: Vec<f64> = (0..=doses.len()).map(|k| seed_change - k as f64).collect();
        json!({
            "weeks": weeks, "doses": doses, "baseline_mean": baseline, "change": change,
            "item_missing": missing, "decline_sd": decline, "dropout_rate": dropout,
        })
    }
}

prop_compose! {
    fn dose_ranging_controls()(
        indication in prop::sample::select(vec!["type_2_diabetes", "hypertension", "obesity"]),
        weeks in 4u32..=52,
        doses in prop::collection::btree_set(1u32..500, 1..6),
        ed50 in 0.5f64..200.0,
        hill in 0.3f64..4.0,
        emax in prop::option::of(-20.0f64..5.0),
        dropout in 0.0f64..2.0,
    ) -> serde_json::Value {
        let mut d: Vec<f64> = vec![0.0];
        d.extend(doses.into_iter().map(f64::from));
        json!({
            "indication": indication, "weeks": weeks, "doses": d, "ed50": ed50, "hill": hill,
            "emax": emax, "dropout_rate": dropout,
        })
    }
}

prop_compose! {
    fn crossover_controls()(
        design in prop::sample::select(vec!["2x2", "2x2x4", "food_effect"]),
        dose in 1.0f64..1000.0,
        gmr in 0.5f64..2.0,
        ka_ratio in 0.3f64..3.0,
        intra in 0.0f64..0.6,
        inter in 0.0f64..0.8,
        days in 5u32..20,
        lloq in 0.01f64..5.0,
        dropout in 0.0f64..5.0,
    ) -> serde_json::Value {
        json!({
            "design": design, "dose": dose, "gmr": gmr, "ka_ratio": ka_ratio,
            "intra_cv": intra, "inter_cv": inter, "days_between_doses": days, "lloq": lloq,
            "dropout_rate": dropout,
        })
    }
}

prop_compose! {
    fn phase1_controls()(
        mad in any::<bool>(),
        doses in prop::collection::vec(1.0f64..1000.0, 1..5),
        placebo in 0.0f64..0.5,
        days in 2u32..10,
        compartments in 1u8..=2,
        cl in 1.0f64..50.0,
        v in 10.0f64..300.0,
        ka in 0.2f64..5.0,
        cv in (0.0f64..0.8, 0.0f64..0.8, 0.0f64..0.8),
        beta in 0.5f64..1.5,
        lloq in 0.01f64..5.0,
    ) -> serde_json::Value {
        json!({
            "part": if mad { "mad" } else { "sad" }, "doses": doses, "placebo_fraction": placebo,
            "mad_days": days, "compartments": compartments, "cl": cl, "v": v, "ka": ka,
            "cv": [cv.0, cv.1, cv.2], "proportionality": beta, "lloq": lloq,
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
    fn random_cv_specs_verify(c in cv_controls(), n in 1usize..300, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("cv_outcomes", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_cns_specs_verify(c in cns_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("cns_alzheimers", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_dose_ranging_specs_verify(c in dose_ranging_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("dose_ranging_ph2", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_crossover_specs_verify(c in crossover_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("crossover_be", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_phase1_specs_verify(c in phase1_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("phase1_sad_mad", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_oncology_specs_verify(c in oncology_controls(), n in 1usize..80, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("oncology_solid_tumor", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }
}
