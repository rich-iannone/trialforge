//! L1: referential and temporal checks across domains.

use std::collections::{HashMap, HashSet};

use trialforge_core::{Dataset, Date, study_day};

use crate::table::{T, get};
use crate::{Level, Report};

const L: Level = Level::L1;

/// Per-subject reference dates from DM.
struct Ref {
    rfst: Option<Date>,
    rfxst: Option<Date>,
    rfxen: Option<Date>,
    rfic: Option<Date>,
    rfpen: Option<Date>,
    death: Option<Date>,
}

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let Some(dm) = get(sdtm, "DM") else {
        r.add(L, "L1.SUBJ", "DM", None, None, "SDTM has no DM".into());
        return;
    };
    let u = dm.c("USUBJID").unwrap_or_default();
    let col = |v: &str| dm.d(v).unwrap_or_else(|| vec![None; dm.n()]);
    let (rfst, rfxst, rfxen, rfic, rfpen, dth) = (
        col("RFSTDTC"),
        col("RFXSTDTC"),
        col("RFXENDTC"),
        col("RFICDTC"),
        col("RFPENDTC"),
        col("DTHDTC"),
    );
    let refs: HashMap<&str, Ref> = (0..dm.n())
        .filter_map(|i| {
            Some((
                u[i].as_deref()?,
                Ref {
                    rfst: rfst[i],
                    rfxst: rfxst[i],
                    rfxen: rfxen[i],
                    rfic: rfic[i],
                    rfpen: rfpen[i],
                    death: dth[i],
                },
            ))
        })
        .collect();

    // DTHFL iff DTHDTC.
    if let (Some(fl), Some(dtc)) = (dm.c("DTHFL"), dm.c("DTHDTC")) {
        for i in 0..dm.n() {
            if (fl[i].as_deref() == Some("Y")) != dtc[i].is_some() {
                r.add(
                    L,
                    "L1.DTHFL",
                    "DM",
                    Some("DTHFL"),
                    Some(i),
                    "DTHFL must be Y exactly when DTHDTC is present".into(),
                );
            }
        }
    }

    // Every subject record belongs to a DM subject.
    for ds in sdtm.iter().chain(adam).filter(|d| d.name != "DM") {
        let t = T::new(ds);
        if let Some(us) = t.c("USUBJID") {
            for (i, x) in us.iter().enumerate() {
                if x.as_deref().is_some_and(|x| !refs.contains_key(x)) {
                    r.add(
                        L,
                        "L1.SUBJ",
                        t.name(),
                        Some("USUBJID"),
                        Some(i),
                        format!("{x:?} is not in DM"),
                    );
                }
            }
        }
    }

    for ds in sdtm {
        let t = T::new(ds);
        study_days(&t, &refs, r);
        start_end(&t, r);
        within_participation(&t, &refs, r);
    }
    exposure(sdtm, &refs, r);
    visits(sdtm, r);
    epochs(sdtm, r);
    arms(sdtm, r);
    supp(sdtm, r);
    baseline_flags(sdtm, &refs, r);
    disposition(sdtm, r);
    seriousness(sdtm, r);
    crate::onco::l1(sdtm, r);
}

/// Pairs of (study day, date) variables in a dataset.
fn day_pairs(t: &T<'_>) -> Vec<(String, String)> {
    let p = t.name().to_string();
    let mut out = Vec::new();
    for (dy, dtc) in [
        (format!("{p}STDY"), format!("{p}STDTC")),
        (format!("{p}ENDY"), format!("{p}ENDTC")),
        (format!("{p}DY"), format!("{p}DTC")),
    ] {
        if t.has(&dy) && t.has(&dtc) {
            out.push((dy, dtc));
        }
    }
    out
}

fn study_days(t: &T<'_>, refs: &HashMap<&str, Ref>, r: &mut Report) {
    let Some(u) = t.c("USUBJID") else { return };
    for (dy, dtc) in day_pairs(t) {
        let (days, dates) = (t.i(&dy).unwrap_or_default(), t.d(&dtc).unwrap_or_default());
        for (i, subject) in u.iter().enumerate() {
            let rf = subject
                .as_deref()
                .and_then(|x| refs.get(x))
                .and_then(|x| x.rfst);
            let expected = match (dates.get(i).copied().flatten(), rf) {
                (Some(d), Some(rf)) => Some(i64::from(study_day(d, rf))),
                _ => None,
            };
            if expected.is_some() && days.get(i).copied().flatten() != expected {
                r.add(
                    L,
                    "L1.DY",
                    t.name(),
                    Some(&dy),
                    Some(i),
                    format!("{dy} should be {expected:?}"),
                );
            }
        }
    }
}

fn start_end(t: &T<'_>, r: &mut Report) {
    let p = t.name();
    let (s, e) = (format!("{p}STDTC"), format!("{p}ENDTC"));
    let (Some(st), Some(en)) = (t.d(&s), t.d(&e)) else {
        return;
    };
    for i in 0..t.n() {
        if let (Some(a), Some(b)) = (st[i], en[i]) {
            if b < a {
                r.add(
                    L,
                    "L1.ORDER",
                    p,
                    Some(&e),
                    Some(i),
                    format!("{e} is before {s}"),
                );
            }
        }
    }
}

fn within_participation(t: &T<'_>, refs: &HashMap<&str, Ref>, r: &mut Report) {
    let p = t.name();
    if matches!(
        p,
        "DM" | "CM" | "MH" | "TA" | "TE" | "TV" | "TI" | "TS" | "SUPPAE"
    ) {
        return;
    }
    let Some(u) = t.c("USUBJID") else { return };
    let date_vars: Vec<String> = [format!("{p}STDTC"), format!("{p}DTC"), format!("{p}ENDTC")]
        .into_iter()
        .filter(|v| t.has(v))
        .collect();
    for v in &date_vars {
        let dates = t.d(v).unwrap_or_default();
        for i in 0..t.n() {
            let (Some(d), Some(rf)) = (dates[i], u[i].as_deref().and_then(|x| refs.get(x))) else {
                continue;
            };
            if rf.rfic.is_some_and(|c| d < c) {
                r.add(
                    L,
                    "L1.PART",
                    p,
                    Some(v),
                    Some(i),
                    "date is before informed consent".into(),
                );
            }
            if rf.rfpen.is_some_and(|e| d > e) {
                r.add(
                    L,
                    "L1.PART",
                    p,
                    Some(v),
                    Some(i),
                    "date is after the end of participation".into(),
                );
            }
            if rf.death.is_some_and(|dd| d > dd) {
                r.add(
                    L,
                    "L1.DEATH",
                    p,
                    Some(v),
                    Some(i),
                    "date is after death".into(),
                );
            }
        }
    }
}

fn exposure(sdtm: &[Dataset], refs: &HashMap<&str, Ref>, r: &mut Report) {
    let Some(ex) = get(sdtm, "EX") else { return };
    let (Some(u), Some(st), Some(en)) = (ex.c("USUBJID"), ex.d("EXSTDTC"), ex.d("EXENDTC")) else {
        return;
    };
    let mut span: HashMap<&str, (Date, Date)> = HashMap::new();
    for i in 0..ex.n() {
        let (Some(x), Some(a), Some(b)) = (u[i].as_deref(), st[i], en[i]) else {
            continue;
        };
        let e = span.entry(x).or_insert((a, b));
        e.0 = e.0.min(a);
        e.1 = e.1.max(b);
        if let Some(rf) = refs.get(x) {
            if rf.rfxst.is_some_and(|s| a < s) || rf.rfxen.is_some_and(|f| b > f) {
                r.add(
                    L,
                    "L1.EXPO",
                    "EX",
                    None,
                    Some(i),
                    "exposure outside RFXSTDTC..RFXENDTC".into(),
                );
            }
        }
    }
    for (x, rf) in refs {
        match span.get(x) {
            Some((a, b)) => {
                if rf.rfxst != Some(*a) || rf.rfxen != Some(*b) {
                    r.add(
                        L,
                        "L1.EXPO",
                        "DM",
                        None,
                        None,
                        format!("{x}: RFXSTDTC/RFXENDTC do not match EX"),
                    );
                }
            }
            None if rf.rfxst.is_some() => r.add(
                L,
                "L1.EXPO",
                "DM",
                None,
                None,
                format!("{x}: RFXSTDTC without EX records"),
            ),
            None => {}
        }
    }
}

fn visits(sdtm: &[Dataset], r: &mut Report) {
    let Some(tv) = get(sdtm, "TV") else { return };
    let (Some(num), Some(name), Some(day)) = (tv.i("VISITNUM"), tv.c("VISIT"), tv.i("VISITDY"))
    else {
        return;
    };
    let planned: HashMap<i64, (Option<&str>, Option<i64>)> = (0..tv.n())
        .filter_map(|i| Some((num[i]?, (name[i].as_deref(), day[i]))))
        .collect();
    // SV visit dates, for matching findings.
    let mut sv_dates: HashMap<(&str, i64), Date> = HashMap::new();
    if let Some(sv) = get(sdtm, "SV") {
        if let (Some(u), Some(n), Some(d)) = (sv.c("USUBJID"), sv.i("VISITNUM"), sv.d("SVSTDTC")) {
            for i in 0..sv.n() {
                if let (Some(u), Some(n), Some(d)) = (u[i].as_deref(), n[i], d[i]) {
                    sv_dates.insert((u, n), d);
                }
            }
        }
    }
    for ds in sdtm
        .iter()
        .filter(|d| matches!(d.name.as_str(), "SV" | "LB" | "VS" | "TU" | "TR" | "RS"))
    {
        let t = T::new(ds);
        let (Some(u), Some(n), Some(v)) = (t.c("USUBJID"), t.i("VISITNUM"), t.c("VISIT")) else {
            continue;
        };
        let vdy = t.i("VISITDY");
        let date_var = if t.name() == "SV" {
            "SVSTDTC".to_string()
        } else {
            format!("{}DTC", t.name())
        };
        let dates = t.d(&date_var).unwrap_or_default();
        for i in 0..t.n() {
            let Some(vn) = n[i] else { continue };
            match planned.get(&vn) {
                None => r.add(
                    L,
                    "L1.VISIT",
                    t.name(),
                    Some("VISITNUM"),
                    Some(i),
                    format!("VISITNUM {vn} is not in TV"),
                ),
                Some((pname, pday)) => {
                    if v[i].as_deref() != *pname || vdy.as_ref().is_some_and(|d| d[i] != *pday) {
                        r.add(
                            L,
                            "L1.VISIT",
                            t.name(),
                            Some("VISIT"),
                            Some(i),
                            "VISIT/VISITDY do not match TV".into(),
                        );
                    }
                }
            }
            if t.name() != "SV" {
                let sv = u[i].as_deref().and_then(|u| sv_dates.get(&(u, vn)));
                if sv.copied() != dates.get(i).copied().flatten() {
                    r.add(
                        L,
                        "L1.SV",
                        t.name(),
                        Some(&date_var),
                        Some(i),
                        "date does not match the SV visit date".into(),
                    );
                }
            }
        }
    }
}

fn epochs(sdtm: &[Dataset], r: &mut Report) {
    let Some(se) = get(sdtm, "SE") else { return };
    let (Some(u), Some(ep), Some(st), Some(en)) = (
        se.c("USUBJID"),
        se.c("EPOCH"),
        se.d("SESTDTC"),
        se.d("SEENDTC"),
    ) else {
        return;
    };
    let mut spans: HashMap<&str, Vec<(&str, Date, Date)>> = HashMap::new();
    for i in 0..se.n() {
        if let (Some(u), Some(e), Some(a), Some(b)) =
            (u[i].as_deref(), ep[i].as_deref(), st[i], en[i])
        {
            spans.entry(u).or_default().push((e, a, b));
        }
    }
    for ds in sdtm
        .iter()
        .filter(|d| !matches!(d.name.as_str(), "SE" | "TA" | "TE" | "TV" | "TI" | "TS"))
    {
        let t = T::new(ds);
        let (Some(us), Some(eps)) = (t.c("USUBJID"), t.c("EPOCH")) else {
            continue;
        };
        let p = t.name();
        let var = [
            format!("{p}STDTC"),
            format!("{p}DTC"),
            "SVSTDTC".to_string(),
        ]
        .into_iter()
        .find(|v| t.has(v));
        let Some(var) = var else { continue };
        let dates = t.d(&var).unwrap_or_default();
        for i in 0..t.n() {
            let (Some(x), Some(e), Some(d)) = (us[i].as_deref(), eps[i].as_deref(), dates[i])
            else {
                continue;
            };
            let ok = spans
                .get(x)
                .is_some_and(|s| s.iter().any(|(ee, a, b)| *ee == e && *a <= d && d <= *b));
            if !ok {
                r.add(
                    L,
                    "L1.EPOCH",
                    p,
                    Some("EPOCH"),
                    Some(i),
                    format!("{e} does not cover {d}"),
                );
            }
        }
    }
}

fn arms(sdtm: &[Dataset], r: &mut Report) {
    let (Some(dm), Some(ta)) = (get(sdtm, "DM"), get(sdtm, "TA")) else {
        return;
    };
    let ta_arms: HashSet<&str> = ta
        .c("ARMCD")
        .unwrap_or_default()
        .iter()
        .flatten()
        .map(String::as_str)
        .collect();
    for (i, a) in dm.c("ARMCD").unwrap_or_default().iter().enumerate() {
        if a.as_deref().is_some_and(|a| !ta_arms.contains(a)) {
            r.add(
                L,
                "L1.ARM",
                "DM",
                Some("ARMCD"),
                Some(i),
                format!("{a:?} is not in TA"),
            );
        }
    }
}

fn supp(sdtm: &[Dataset], r: &mut Report) {
    let (Some(ae), Some(sup)) = (get(sdtm, "AE"), get(sdtm, "SUPPAE")) else {
        return;
    };
    let keys: HashSet<(&str, String)> = ae
        .c("USUBJID")
        .unwrap_or_default()
        .iter()
        .zip(ae.i("AESEQ").unwrap_or_default())
        .filter_map(|(u, s)| Some((u.as_deref()?, s?.to_string())))
        .collect();
    let (Some(u), Some(idvar), Some(val), Some(rdom)) = (
        sup.c("USUBJID"),
        sup.c("IDVAR"),
        sup.c("IDVARVAL"),
        sup.c("RDOMAIN"),
    ) else {
        return;
    };
    for i in 0..sup.n() {
        let ok = rdom[i].as_deref() == Some("AE")
            && idvar[i].as_deref() == Some("AESEQ")
            && matches!((u[i].as_deref(), &val[i]), (Some(u), Some(v)) if keys.contains(&(u, v.clone())));
        if !ok {
            r.add(
                L,
                "L1.SUPP",
                "SUPPAE",
                Some("IDVARVAL"),
                Some(i),
                "does not resolve to an AE record".into(),
            );
        }
    }
}

fn baseline_flags(sdtm: &[Dataset], refs: &HashMap<&str, Ref>, r: &mut Report) {
    for name in ["LB", "VS"] {
        let Some(t) = get(sdtm, name) else { continue };
        let fl = format!("{name}LOBXFL");
        let (Some(u), Some(tc), Some(f), Some(d)) = (
            t.c("USUBJID"),
            t.c(&format!("{name}TESTCD")),
            t.c(&fl),
            t.d(&format!("{name}DTC")),
        ) else {
            continue;
        };
        let mut seen = HashSet::new();
        for i in 0..t.n() {
            if f[i].as_deref() != Some("Y") {
                continue;
            }
            let (Some(x), Some(c)) = (u[i].as_deref(), tc[i].as_deref()) else {
                continue;
            };
            if !seen.insert((x, c)) {
                r.add(
                    L,
                    "L1.LOBXFL",
                    name,
                    Some(&fl),
                    Some(i),
                    "more than one baseline record per test".into(),
                );
            }
            if let (Some(dd), Some(rf)) = (d[i], refs.get(x).and_then(|r| r.rfxst)) {
                if dd > rf {
                    r.add(
                        L,
                        "L1.LOBXFL",
                        name,
                        Some(&fl),
                        Some(i),
                        "baseline record is after first exposure".into(),
                    );
                }
            }
        }
    }
}

fn disposition(sdtm: &[Dataset], r: &mut Report) {
    let Some(ds) = get(sdtm, "DS") else { return };
    let (Some(u), Some(dec), Some(scat)) = (ds.c("USUBJID"), ds.c("DSDECOD"), ds.c("DSSCAT"))
    else {
        return;
    };
    let mut counts: HashMap<(&str, &str), usize> = HashMap::new();
    for i in 0..ds.n() {
        let Some(x) = u[i].as_deref() else { continue };
        let what = match (dec[i].as_deref(), scat[i].as_deref()) {
            (Some("RANDOMIZED"), _) => "RANDOMIZED",
            (Some("INFORMED CONSENT OBTAINED"), _) => "CONSENT",
            (_, Some(s)) => s,
            _ => continue,
        };
        *counts.entry((x, what)).or_default() += 1;
    }
    let subjects: HashSet<&str> = u.iter().flatten().map(String::as_str).collect();
    for x in subjects {
        for (what, exactly_one) in [
            ("RANDOMIZED", true),
            ("CONSENT", true),
            ("END OF TREATMENT", false),
            ("END OF STUDY", false),
        ] {
            let k = counts.get(&(x, what)).copied().unwrap_or(0);
            if k > 1 || (exactly_one && k != 1) {
                r.add(
                    L,
                    "L1.DS",
                    "DS",
                    None,
                    None,
                    format!("{x}: {k} {what} record(s)"),
                );
            }
        }
    }
}

fn seriousness(sdtm: &[Dataset], r: &mut Report) {
    let Some(ae) = get(sdtm, "AE") else { return };
    let Some(ser) = ae.c("AESER") else { return };
    let crit: Vec<&[Option<String>]> = [
        "AESCONG", "AESDISAB", "AESDTH", "AESHOSP", "AESLIFE", "AESMIE",
    ]
    .iter()
    .filter_map(|v| ae.c(v))
    .collect();
    for i in 0..ae.n() {
        let any_y = crit.iter().any(|c| c[i].as_deref() == Some("Y"));
        let any_set = crit.iter().any(|c| c[i].is_some());
        match ser[i].as_deref() {
            Some("Y") if !any_y => r.add(
                L,
                "L1.SER",
                "AE",
                Some("AESER"),
                Some(i),
                "serious AE without a seriousness criterion".into(),
            ),
            Some("N") if any_set => r.add(
                L,
                "L1.SER",
                "AE",
                Some("AESER"),
                Some(i),
                "non-serious AE with seriousness criteria".into(),
            ),
            _ => {}
        }
    }
}
