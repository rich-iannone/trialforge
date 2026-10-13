//! End-to-end checks on the oncology spike: the SDTM/ADaM invariants that verification levels
//! 0–2 will enforce on every generation, run here across many seeds.

use std::collections::{HashMap, HashSet};

use trialforge::{Dataset, Date, StudyOptions, Values, study};

fn opts(n: usize, seed: u64) -> StudyOptions {
    StudyOptions {
        n,
        seed,
        studyid: None,
    }
}

fn chars<'a>(ds: &'a Dataset, name: &str) -> Vec<Option<&'a str>> {
    match &ds
        .column(name)
        .unwrap_or_else(|| panic!("{}.{name}", ds.name))
        .values
    {
        Values::Char(v) => v.iter().map(Option::as_deref).collect(),
        _ => panic!("{}.{name} is not character", ds.name),
    }
}

fn ints(ds: &Dataset, name: &str) -> Vec<Option<i64>> {
    match &ds.column(name).unwrap().values {
        Values::Int(v) => v.clone(),
        _ => panic!("{}.{name} is not integer", ds.name),
    }
}

fn dates(ds: &Dataset, name: &str) -> Vec<Option<Date>> {
    match &ds.column(name).unwrap().values {
        Values::Date(v) => v.clone(),
        _ => panic!("{}.{name} is not a date", ds.name),
    }
}

fn iso(v: Option<&str>) -> Option<Date> {
    v.map(|s| Date::parse_iso(s).unwrap_or_else(|| panic!("not an ISO date: {s}")))
}

/// Per-subject DM facts keyed by USUBJID.
struct Dm {
    rfst: Date,
    rfxen: Date,
    rfic: Date,
    rfpen: Date,
    death: Option<Date>,
}

fn dm_index(dm: &Dataset) -> HashMap<String, Dm> {
    let u = chars(dm, "USUBJID");
    let (rfst, rfxst, rfxen, rfic, rfpen, dth) = (
        chars(dm, "RFSTDTC"),
        chars(dm, "RFXSTDTC"),
        chars(dm, "RFXENDTC"),
        chars(dm, "RFICDTC"),
        chars(dm, "RFPENDTC"),
        chars(dm, "DTHDTC"),
    );
    (0..dm.n_rows())
        .map(|i| {
            assert_eq!(rfst[i], rfxst[i], "RFSTDTC is first dose in this design");
            (
                u[i].unwrap().to_string(),
                Dm {
                    rfst: iso(rfst[i]).unwrap(),
                    rfxen: iso(rfxen[i]).unwrap(),
                    rfic: iso(rfic[i]).unwrap(),
                    rfpen: iso(rfpen[i]).unwrap(),
                    death: iso(dth[i]),
                },
            )
        })
        .collect()
}

fn study_day(d: Date, rf: Date) -> i64 {
    i64::from(trialforge::core::study_day(d, rf))
}

fn check_sequences(ds: &Dataset, seq: &str) {
    let u = chars(ds, "USUBJID");
    let s = ints(ds, seq);
    let mut expected: HashMap<&str, i64> = HashMap::new();
    for i in 0..ds.n_rows() {
        let next = expected.entry(u[i].unwrap()).or_insert(1);
        assert_eq!(
            s[i],
            Some(*next),
            "{}.{seq} runs 1..k within subject, in order",
            ds.name
        );
        *next += 1;
    }
}

fn check_study(n: usize, seed: u64) {
    let st = study("oncology_solid_tumor", &opts(n, seed)).unwrap();
    let (dm, ex, ds) = (
        st.dataset("DM").unwrap(),
        st.dataset("EX").unwrap(),
        st.dataset("DS").unwrap(),
    );
    let adsl = st.dataset("ADSL").unwrap();
    let cutoff = st.cutoff;

    // DM: one row per subject, unique and sorted keys.
    assert_eq!(dm.n_rows(), n);
    let subj = dm_index(dm);
    assert_eq!(subj.len(), n, "USUBJID unique in DM");
    let u = chars(dm, "USUBJID");
    assert!(u.windows(2).all(|w| w[0] < w[1]), "DM sorted by USUBJID");

    // DM internal dates.
    let dmdy = ints(dm, "DMDY");
    for (i, id) in u.iter().enumerate() {
        let s = &subj[id.unwrap()];
        assert!(s.rfic < s.rfst, "consent before first dose");
        assert!(s.rfst <= s.rfxen && s.rfxen <= s.rfpen && s.rfpen <= cutoff);
        if let Some(d) = s.death {
            assert_eq!(d, s.rfpen, "death ends participation");
            assert!(s.rfxen < d, "no dosing on or after death");
        }
        assert_eq!(dmdy[i], Some(study_day(s.rfic, s.rfst)));
        assert_ne!(dmdy[i], Some(0));
    }
    let dthfl = chars(dm, "DTHFL");
    for (i, id) in u.iter().enumerate() {
        assert_eq!(dthfl[i] == Some("Y"), subj[id.unwrap()].death.is_some());
    }

    // EX: in DM, sequenced, within the exposure window, study days right.
    check_sequences(ex, "EXSEQ");
    let (eu, est, edy) = (
        chars(ex, "USUBJID"),
        chars(ex, "EXSTDTC"),
        ints(ex, "EXSTDY"),
    );
    let mut last_dose: HashMap<&str, Date> = HashMap::new();
    for i in 0..ex.n_rows() {
        let s = &subj[eu[i].unwrap()];
        let d = iso(est[i]).unwrap();
        assert!(
            s.rfst <= d && d <= s.rfxen,
            "dose within RFXSTDTC..RFXENDTC"
        );
        assert_eq!(edy[i], Some(study_day(d, s.rfst)));
        let e = last_dose.entry(eu[i].unwrap()).or_insert(d);
        *e = (*e).max(d);
    }
    assert_eq!(last_dose.len(), n, "every subject is dosed");
    for (id, d) in &last_dose {
        assert_eq!(*d, subj[*id].rfxen, "RFXENDTC is the last dose");
    }

    // DS: in DM, sequenced, milestones present once, nothing after death or cutoff.
    check_sequences(ds, "DSSEQ");
    let (du, dd, dt, dep, ddy) = (
        chars(ds, "USUBJID"),
        chars(ds, "DSDECOD"),
        chars(ds, "DSSTDTC"),
        chars(ds, "EPOCH"),
        ints(ds, "DSSTDY"),
    );
    let mut milestones: HashMap<(&str, &str), u32> = HashMap::new();
    let mut seen = HashSet::new();
    for i in 0..ds.n_rows() {
        let id = du[i].unwrap();
        let s = &subj[id];
        let d = iso(dt[i]).unwrap();
        assert!(d <= cutoff);
        if let Some(death) = s.death {
            assert!(d <= death, "no disposition after death");
        }
        assert_eq!(ddy[i], Some(study_day(d, s.rfst)));
        assert!(
            seen.insert((id, dep[i].unwrap(), dd[i].unwrap())),
            "no duplicate DS events"
        );
        if matches!(dd[i], Some("INFORMED CONSENT OBTAINED" | "RANDOMIZED")) {
            *milestones.entry((id, dd[i].unwrap())).or_default() += 1;
            if dd[i] == Some("RANDOMIZED") {
                assert!(
                    s.rfic < d && d <= s.rfst,
                    "consent < randomisation <= first dose"
                );
            }
        }
    }
    assert_eq!(milestones.len(), 2 * n);
    assert!(milestones.values().all(|&k| k == 1));

    // ADSL agrees with SDTM.
    assert_eq!(adsl.n_rows(), n);
    assert_eq!(chars(adsl, "USUBJID"), u);
    let (trtsdt, trtedt, dthdt, lstalv) = (
        dates(adsl, "TRTSDT"),
        dates(adsl, "TRTEDT"),
        dates(adsl, "DTHDT"),
        dates(adsl, "LSTALVDT"),
    );
    let (eosstt, dcsreas, saffl) = (
        chars(adsl, "EOSSTT"),
        chars(adsl, "DCSREAS"),
        chars(adsl, "SAFFL"),
    );
    for (i, id) in u.iter().enumerate() {
        let s = &subj[id.unwrap()];
        assert_eq!(trtsdt[i], Some(s.rfst));
        assert_eq!(trtedt[i], Some(s.rfxen));
        assert_eq!(dthdt[i], s.death);
        assert_eq!(saffl[i], Some("Y"));
        assert!(lstalv[i] >= trtedt[i]);
        if let Some(d) = s.death {
            assert_eq!(dcsreas[i], Some("DEATH"));
            assert_eq!(lstalv[i], Some(d));
        }
        assert_eq!(eosstt[i] == Some("DISCONTINUED"), dcsreas[i].is_some());
    }
}

#[test]
fn invariants_hold_across_seeds() {
    for seed in 0..40 {
        check_study(120, seed);
    }
    check_study(1, 7);
    check_study(1_000, 2026);
}

#[test]
fn generation_is_deterministic() {
    let a = study("oncology_solid_tumor", &opts(300, 1)).unwrap();
    let b = study("oncology_solid_tumor", &opts(300, 1)).unwrap();
    let c = study("oncology_solid_tumor", &opts(300, 2)).unwrap();
    assert_eq!(a, b);
    assert_ne!(a.sdtm, c.sdtm);
    let xpt = |s: &trialforge::Study| {
        s.sdtm
            .iter()
            .chain(&s.adam)
            .map(|d| trialforge::io::xpt::to_bytes(d, s.cutoff).unwrap())
            .collect::<Vec<_>>()
    };
    assert_eq!(xpt(&a), xpt(&b));
}

#[test]
fn arms_are_balanced_and_outcomes_plausible() {
    let st = study("oncology_solid_tumor", &opts(300, 1)).unwrap();
    let adsl = st.dataset("ADSL").unwrap();
    let arm = chars(adsl, "ARMCD");
    let on_a = arm.iter().filter(|a| **a == Some("A")).count();
    assert!((146..=154).contains(&on_a), "1:1 permuted blocks: {on_a}");
    let deaths = chars(adsl, "DTHFL")
        .iter()
        .filter(|f| **f == Some("Y"))
        .count();
    assert!(
        (60..=240).contains(&deaths),
        "deaths by the 30-month cutoff: {deaths}"
    );
    let reasons: HashSet<_> = chars(adsl, "DCTREAS").into_iter().flatten().collect();
    for r in ["PROGRESSIVE DISEASE", "DEATH", "ADVERSE EVENT"] {
        assert!(
            reasons.contains(r),
            "missing treatment discontinuation reason {r}"
        );
    }
}

#[test]
fn unknown_archetypes_and_formats_are_rejected() {
    assert!(matches!(
        study("nope", &opts(10, 1)),
        Err(trialforge::Error::Invalid(_))
    ));
    let st = study("oncology_solid_tumor", &opts(10, 1)).unwrap();
    let dir = std::env::temp_dir().join("trialforge-unknown-format");
    assert!(matches!(
        st.write(&dir, &["sas7bdat"]),
        Err(trialforge::Error::Invalid(_))
    ));
}
