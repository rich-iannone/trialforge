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

/// A study generated, or rejected only because an interim cutoff precedes every
/// randomisation (possible for small random specs).
fn generated_or_empty_cut(s: &Result<trialforge::Study, trialforge::Error>) -> bool {
    match s {
        Ok(_) => true,
        Err(trialforge::Error::Invalid(e)) => e.contains("randomised by the data cutoff"),
        Err(_) => false,
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
            let nullable = [
                "change",
                "ni_margin",
                "placebo_change",
                "emax",
                "crossover_os_hr",
                "interim_months",
                "interim_events",
            ];
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

#[test]
fn ophthalmology_bcva_is_verified_and_recovered() {
    let s = study(
        "ophthalmology",
        &opts(1500, 4, json!({"bcva_change": [8.0, 5.0]})),
    )
    .unwrap();
    let names: Vec<&str> = s
        .sdtm
        .iter()
        .chain(&s.adam)
        .map(|d| d.name.as_str())
        .collect();
    for want in ["OE", "SC", "ADBCVA"] {
        assert!(names.contains(&want), "missing {want}: {names:?}");
    }
    let ex = s.dataset("EX").unwrap();
    assert!(
        chars(ex, "EXROUTE")
            .iter()
            .all(|r| r.as_deref() == Some("INTRAVITREAL"))
    );
    let b = &s.truth["bcva"]["realised"];
    let (c, d) = (
        b["Comparator 2 mg"]["mean_change"].as_f64().unwrap(),
        b["Drug A 6 mg"]["mean_change"].as_f64().unwrap(),
    );
    assert!((c - 8.0).abs() < 1.0 && (d - 5.0).abs() < 1.0, "{c} {d}");
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    corrupt(&mut bad.adam, "ADBCVA", "AFEYE", "Study Eye", "Fellow Eye");
    edit(&mut bad.adam, "ADBCVA", "CHG", |v| {
        let Values::Num(v) = v else {
            panic!("CHG is numeric")
        };
        let x = v.iter_mut().find(|x| x.is_some()).unwrap();
        *x = Some(x.unwrap() + 1.0);
    });
    assert!(rules(&bad).contains(&"L2.ADBCVA"), "{:?}", rules(&bad));
}

#[test]
fn pediatric_growth_is_verified_and_recovered() {
    let s = study(
        "pediatric_growth",
        &opts(1500, 6, json!({"height_velocity": [10.0, 9.0]})),
    )
    .unwrap();
    let dm = s.dataset("DM").unwrap();
    assert!(chars(dm, "BRTHDTC").iter().all(Option::is_some));
    let advs = s.dataset("ADVS").unwrap();
    let params = chars(advs, "PARAMCD");
    for want in ["HGTSDS", "HGTPCTL", "WGTASDS", "BMIPCTL", "HTVEL"] {
        assert!(
            params.iter().any(|p| p.as_deref() == Some(want)),
            "missing {want}"
        );
    }
    let hv = &s.truth["height_velocity"]["realised"];
    let (c, d) = (
        hv["Somatropin 0.034 mg/kg/day"]["mean"].as_f64().unwrap(),
        hv["Drug A 0.24 mg/kg/week"]["mean"].as_f64().unwrap(),
    );
    assert!((c - 10.0).abs() < 0.4 && (d - 9.0).abs() < 0.4, "{c} {d}");
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    edit(&mut bad.adam, "ADVS", "AVAL", |v| {
        let Values::Num(v) = v else {
            panic!("AVAL is numeric")
        };
        let i = (0..v.len())
            .find(|&i| params[i].as_deref() == Some("HGTSDS"))
            .unwrap();
        v[i] = v[i].map(|x| x + 0.1);
    });
    assert!(rules(&bad).contains(&"L2.GROWTH"), "{:?}", rules(&bad));
    let mut bad = s.clone();
    edit(&mut bad.sdtm, "DM", "AGE", |v| {
        let Values::Int(v) = v else {
            panic!("AGE is integer")
        };
        v[0] = v[0].map(|a| a + 1);
    });
    assert!(rules(&bad).contains(&"L1.AGE"), "{:?}", rules(&bad));
}

#[test]
fn respiratory_exacerbations_are_verified_and_recovered() {
    let s = study(
        "respiratory_exacerbation",
        &opts(2000, 8, json!({"rate_ratio": 0.6, "dispersion": 0.4})),
    )
    .unwrap();
    let names: Vec<&str> = s
        .sdtm
        .iter()
        .chain(&s.adam)
        .map(|d| d.name.as_str())
        .collect();
    for want in ["RE", "CE", "ADRE", "ADEXAC", "ADTTE"] {
        assert!(names.contains(&want), "missing {want}: {names:?}");
    }
    let rr = s.truth["exacerbations"]["rate_ratio"]["estimate"]
        .as_f64()
        .unwrap();
    assert!((rr - 0.6).abs() < 0.08, "rate ratio {rr}");
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    edit(&mut bad.adam, "ADEXAC", "AVAL", |v| {
        let Values::Num(v) = v else {
            panic!("AVAL is numeric")
        };
        v[0] = v[0].map(|x| x + 1.0);
    });
    assert!(rules(&bad).contains(&"L2.ADEXAC"), "{:?}", rules(&bad));
}

#[test]
fn vaccine_immunogenicity_and_reactogenicity_are_verified() {
    let s = study("vaccine", &opts(1000, 9, json!({"gmt": 320.0}))).unwrap();
    let names: Vec<&str> = s
        .sdtm
        .iter()
        .chain(&s.adam)
        .map(|d| d.name.as_str())
        .collect();
    for want in ["IS", "CE", "FACE", "ADIS", "ADCE", "ADFACE"] {
        assert!(names.contains(&want), "missing {want}: {names:?}");
    }
    let gmt = s.truth["vaccine"]["realised"]["Vaccine A"]["gmt"]
        .as_f64()
        .unwrap();
    assert!((gmt / 320.0 - 1.0).abs() < 0.15, "GMT {gmt}");
    let face = s.dataset("FACE").unwrap();
    assert!(
        chars(face, "DOMAIN")
            .iter()
            .all(|d| d.as_deref() == Some("FA"))
    );
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    edit(&mut bad.adam, "ADIS", "AVAL", |v| {
        let Values::Num(v) = v else {
            panic!("AVAL is numeric")
        };
        v[0] = v[0].map(|x| x * 2.0);
    });
    assert!(rules(&bad).contains(&"L2.VACCINE"), "{:?}", rules(&bad));
    let mut bad = s.clone();
    corrupt(&mut bad.adam, "ADCE", "AOCC01FL", "Y", "N");
    assert!(rules(&bad).contains(&"L2.VACCINE"), "{:?}", rules(&bad));
}

#[test]
fn dose_escalation_follows_its_design_and_flags_dlts() {
    for method in ["3+3", "boin", "mtpi2"] {
        let s = study("dose_escalation", &opts(36, 10, json!({"method": method}))).unwrap();
        let e = &s.truth["escalation"];
        let cohorts = e["cohorts"].as_array().unwrap();
        let treated: u64 = cohorts
            .iter()
            .map(|c| c["subjects"].as_u64().unwrap())
            .sum();
        assert_eq!(
            treated as usize,
            s.dataset("DM").unwrap().n_rows(),
            "{method}"
        );
        assert!(treated <= 36, "{method}");
        // Escalation starts at the lowest dose and never skips a level upwards.
        assert_eq!(cohorts[0]["dose"], "Drug A 10 mg", "{method}");
        let sl = s.dataset("ADSL").unwrap();
        assert!(sl.column("DLTEVLFL").is_some(), "{method}");
    }
    let s = study("dose_escalation", &opts(36, 10, json!(null))).unwrap();
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    corrupt(&mut bad.adam, "ADSL", "DLTEVLFL", "Y", "N");
    assert!(rules(&bad).contains(&"L2.DLT"), "{:?}", rules(&bad));
    assert!(matches!(
        study(
            "dose_escalation",
            &opts(36, 1, json!({"dlt_rates": [0.5, 0.1]}))
        ),
        Err(trialforge::Error::Invalid(_))
    ));
}

#[test]
fn lymphoma_follows_lugano_and_recovers_targets() {
    let s = study("heme_lymphoma", &opts(1500, 12, json!(null))).unwrap();
    let rs = s.dataset("RS").unwrap();
    assert!(
        chars(rs, "RSCAT")
            .iter()
            .all(|c| c.as_deref() == Some("LUGANO CLASSIFICATION"))
    );
    let tr = s.dataset("TR").unwrap();
    for want in ["LPERP", "PPD", "SUMPPD", "FDP5PS"] {
        assert!(
            chars(tr, "TRTESTCD")
                .iter()
                .any(|t| t.as_deref() == Some(want)),
            "{want}"
        );
    }
    let orr = s.truth["orr"]["Drug A"]["proportion"].as_f64().unwrap();
    assert!((orr - 0.6).abs() < 0.06, "ORR {orr}");
    assert!(s.dataset("ADTR").is_none());
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    corrupt(&mut bad.sdtm, "RS", "RSSTRESC", "CMR", "PMR");
    assert!(rules(&bad).contains(&"L1.LUGANO"), "{:?}", rules(&bad));
}

#[test]
fn multiple_cuts_are_snapshots_of_the_same_study() {
    let ts_value = |s: &trialforge::Study, parm: &str| {
        let ts = s.dataset("TS").unwrap();
        let i = (0..ts.n_rows())
            .find(|&i| chars(ts, "TSPARMCD")[i].as_deref() == Some(parm))
            .unwrap();
        chars(ts, "TSVAL")[i].clone().unwrap()
    };
    let subjects = |s: &trialforge::Study| -> Vec<(String, String)> {
        let sl = s.dataset("ADSL").unwrap();
        (0..sl.n_rows())
            .map(|i| {
                (
                    chars(sl, "USUBJID")[i].clone().unwrap(),
                    chars(sl, "TRT01P")[i].clone().unwrap(),
                )
            })
            .collect()
    };
    for archetype in ["oncology_solid_tumor", "heme_lymphoma"] {
        let fin = study(archetype, &opts(400, 4, json!(null))).unwrap();
        let mid = study(archetype, &opts(400, 4, json!({"interim_months": 10}))).unwrap();
        assert!(mid.cutoff < fin.cutoff);
        assert_eq!(ts_value(&mid, "DCUTDTC"), mid.cutoff.to_string());
        assert_eq!(ts_value(&mid, "DCUTDESC"), "INTERIM ANALYSIS");
        assert_eq!(ts_value(&mid, "PLANSUB"), "400");
        // Subjects randomised by the interim cutoff, with the same IDs and arms.
        let (m, f) = (subjects(&mid), subjects(&fin));
        assert!(m.len() < f.len(), "{} of {}", m.len(), f.len());
        assert!(m.iter().all(|x| f.contains(x)));
        assert_eq!(
            mid.truth["interim"]["n_randomised"].as_u64(),
            Some(m.len() as u64)
        );
        let deaths = |s: &trialforge::Study| {
            s.truth["overall_survival"]["by_arm"]
                .as_object()
                .unwrap()
                .values()
                .map(|a| a["events"].as_u64().unwrap())
                .sum::<u64>()
        };
        assert!(deaths(&mid) < deaths(&fin));
    }
    // An event-driven interim: the cut falls at the interim event count.
    let mid = study(
        "cv_outcomes",
        &opts(1500, 3, json!({"interim_events": 100})),
    )
    .unwrap();
    assert_eq!(
        mid.truth["outcomes"]["realised_first_mace_events"].as_u64(),
        Some(100)
    );
    assert_eq!(ts_value(&mid, "DCUTDESC"), "INTERIM ANALYSIS");
    assert!(matches!(
        study(
            "oncology_solid_tumor",
            &opts(40, 1, json!({"interim_months": 31}))
        ),
        Err(trialforge::Error::Invalid(_))
    ));
}

#[test]
fn crossover_on_progression_switches_control_subjects() {
    let c = json!({"crossover_prob": 0.8, "os_hr": 0.6});
    let s = study("oncology_solid_tumor", &opts(800, 13, c)).unwrap();
    let x = &s.truth["crossover"];
    let crossed = x["crossed_over"].as_u64().unwrap();
    assert!(crossed > 100, "{crossed} crossed over");
    let ta = s.dataset("TA").unwrap();
    assert!(
        chars(ta, "EPOCH")
            .iter()
            .any(|e| e.as_deref() == Some("CONTINUATION TREATMENT"))
    );
    let ex = s.dataset("EX").unwrap();
    let sl = s.dataset("ADSL").unwrap();
    let control: std::collections::HashSet<&str> = (0..sl.n_rows())
        .filter(|&i| chars(sl, "TRT01P")[i].as_deref() == Some("Control"))
        .filter_map(|i| chars(sl, "USUBJID")[i].as_deref())
        .collect();
    let switched = (0..ex.n_rows()).any(|i| {
        chars(ex, "EXTRT")[i].as_deref() == Some("DRUG A")
            && chars(ex, "USUBJID")[i]
                .as_deref()
                .is_some_and(|u| control.contains(u))
    });
    assert!(switched);
    // The intention-to-treat OS hazard ratio is diluted towards 1.
    let hr = s.truth["overall_survival"]["stratified_cox_hr"]["hr"]
        .as_f64()
        .unwrap();
    assert!(hr > 0.65, "OS HR {hr}");
    let rules = |s: &trialforge::Study| {
        trialforge::verify(s)
            .counts
            .keys()
            .copied()
            .collect::<Vec<_>>()
    };
    let mut bad = s.clone();
    corrupt(&mut bad.adam, "ADSL", "CROSSFL", "Y", "N");
    assert!(rules(&bad).contains(&"L2.CROSSOVER"), "{:?}", rules(&bad));
}

prop_compose! {
    fn escalation_controls()(
        method in prop::sample::select(vec!["3+3", "boin", "mtpi2"]),
        levels in 1usize..7,
        base in 0.0f64..0.3,
        step in 0.0f64..0.2,
        target in 0.1f64..0.4,
        cohort in 1u32..=6,
        weeks in prop::sample::select(vec![4u32, 8, 12, 24]),
        dropout in 0.0f64..3.0,
    ) -> serde_json::Value {
        let doses: Vec<f64> = (0..levels).map(|k| 10.0 * (k as f64 + 1.0)).collect();
        let rates: Vec<f64> = (0..levels).map(|k| (base + step * k as f64).min(1.0)).collect();
        json!({
            "method": method, "doses": doses, "dlt_rates": rates, "target": target,
            "cohort_size": cohort, "weeks": weeks, "dropout_rate": dropout,
        })
    }
}

prop_compose! {
    fn vaccine_controls()(
        doses in 1u32..=2,
        interval in 21u32..=84,
        gmt in 10.0f64..5000.0,
        sero in 0.0f64..1.0,
        frac in 0.0f64..1.0,
        sd in 0.05f64..1.5,
        reacto in 0.0f64..3.0,
        fever in 0.0f64..1.0,
        fu in 56u32..400,
        dropout in 0.0f64..3.0,
        missed in 0.0f64..0.5,
    ) -> serde_json::Value {
        json!({
            "doses": doses, "interval_days": interval, "gmt": gmt, "seropositive": sero,
            "first_dose_fraction": frac, "titer_sd": sd, "reactogenicity": reacto,
            "fever_rate": fever, "followup_days": fu, "dropout_rate": dropout, "missed_visit": missed,
        })
    }
}

prop_compose! {
    fn respiratory_controls()(
        weeks in prop::sample::select(vec![12u32, 24, 52, 104]),
        rate in 0.0f64..5.0,
        rr in 0.2f64..2.0,
        dispersion in 0.0f64..3.0,
        severe in 0.0f64..1.0,
        prior in 0.5f64..4.0,
        fev1 in 0.5f64..3.0,
        change in (-0.5f64..0.5, -0.5f64..0.5),
        dropout in 0.0f64..3.0,
    ) -> serde_json::Value {
        json!({
            "weeks": weeks, "exacerbation_rate": rate, "rate_ratio": rr, "dispersion": dispersion,
            "severe_fraction": severe, "prior_exacerbations_rr": prior, "fev1_baseline": fev1,
            "fev1_change": [change.0, change.1], "dropout_rate": dropout,
        })
    }
}

prop_compose! {
    fn pediatric_controls()(
        weeks in prop::sample::select(vec![26u32, 52, 78, 104]),
        ages in (3i64..=15, 0i64..=6),
        velocity in (1.0f64..20.0, 1.0f64..20.0),
        sd in 0.0f64..5.0,
        sds in -4.5f64..-1.0,
        max_sds in -4.0f64..0.0,
        untreated in 0.5f64..10.0,
        dropout in 0.0f64..3.0,
        missed in 0.0f64..0.5,
    ) -> serde_json::Value {
        json!({
            "weeks": weeks, "age_range": [ages.0, (ages.0 + ages.1).min(15)],
            "height_velocity": [velocity.0, velocity.1], "velocity_sd": sd, "height_sds": sds,
            "max_screening_sds": max_sds, "untreated_velocity": untreated,
            "dropout_rate": dropout, "missed_visit": missed,
        })
    }
}

prop_compose! {
    fn ophthalmology_controls()(
        weeks in prop::sample::select(vec![12u32, 24, 48, 52, 100]),
        loading in 1u32..=3,
        intervals in (prop::sample::select(vec![4u32, 8, 12, 16]), prop::sample::select(vec![4u32, 8, 12, 16])),
        baseline in 30.0f64..75.0,
        change in (-10.0f64..20.0, -10.0f64..20.0),
        response in 0.0f64..20.0,
        cst in (-300.0f64..100.0, -300.0f64..100.0),
        dropout in 0.0f64..3.0,
        missed in 0.0f64..0.5,
    ) -> serde_json::Value {
        json!({
            "weeks": weeks, "loading_doses": loading, "interval_weeks": [intervals.0, intervals.1],
            "bcva_baseline": baseline, "bcva_change": [change.0, change.1], "response_sd": response,
            "cst_change": [cst.0, cst.1], "dropout_rate": dropout, "missed_visit": missed,
        })
    }
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
        interim in prop::option::of(0.0f64..1.0),
    ) -> serde_json::Value {
        let interim = interim
            .filter(|_| target > 1)
            .map(|f| ((target as f64 * f) as usize).clamp(1, target - 1));
        json!({
            "interim_events": interim,
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
        crossover in prop::option::of((0.0f64..1.0, 0.3f64..2.0)),
        interim in prop::option::of(0.0f64..1.0),
    ) -> serde_json::Value {
        let (cross, cross_hr) = crossover.map_or((0.0, None), |(p, h)| (p, Some(h)));
        // An interim cutoff in [1, cutoff_months).
        let interim = interim.map(|f| 1.0 + f * (enrollment + extra - 1.0) * 0.999);
        json!({
            "interim_months": interim,
            "crossover_prob": cross, "crossover_os_hr": cross_hr,
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
        prop_assert!(generated_or_empty_cut(&s), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_ophthalmology_specs_verify(c in ophthalmology_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("ophthalmology", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_pediatric_specs_verify(c in pediatric_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("pediatric_growth", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_respiratory_specs_verify(c in respiratory_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("respiratory_exacerbation", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_vaccine_specs_verify(c in vaccine_controls(), n in 1usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("vaccine", &opts(n, seed, c.clone()));
        prop_assert!(s.is_ok(), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_escalation_specs_verify(c in escalation_controls(), n in 3usize..60, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("dose_escalation", &opts(n, seed, c.clone()));
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
    fn random_lymphoma_specs_verify(c in oncology_controls(), n in 1usize..80, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("heme_lymphoma", &opts(n, seed, c.clone()));
        prop_assert!(generated_or_empty_cut(&s), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }

    #[test]
    fn random_oncology_specs_verify(c in oncology_controls(), n in 1usize..80, seed in 0u64..=trialforge::MAX_SEED) {
        let s = study("oncology_solid_tumor", &opts(n, seed, c.clone()));
        prop_assert!(generated_or_empty_cut(&s), "{c}: {}", s.err().map(|e| e.to_string()).unwrap_or_default());
    }
}
