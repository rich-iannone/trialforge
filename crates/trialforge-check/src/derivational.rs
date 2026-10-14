//! L2: ADaM recomputed from SDTM (independently of `trialforge-adam`) and compared.

use std::collections::HashMap;

use trialforge_cdisc::TEAE_WINDOW_DAYS;
use trialforge_core::{Dataset, Date};

use crate::table::{T, get};
use crate::{Level, Report};

const L: Level = Level::L2;
const TOL: f64 = 1e-6;

/// (start, end) of a source record.
type Span = (Option<Date>, Option<Date>);
/// (value, date, reference low, reference high, test code) of a source finding.
type SourceFinding<'a> = (
    Option<f64>,
    Option<Date>,
    Option<f64>,
    Option<f64>,
    Option<&'a str>,
);

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    adsl(sdtm, adam, r);
    adae(sdtm, adam, r);
    for (bds, src) in [("ADLB", "LB"), ("ADVS", "VS"), ("ADRE", "RE")] {
        if let (Some(b), Some(s)) = (get(adam, bds), get(sdtm, src)) {
            basic(&b, &s, src, adam, r);
        }
    }
    crate::onco::l2(sdtm, adam, r);
    crate::pk::l2(sdtm, adam, r);
    crate::periods::run(sdtm, adam, r);
    crate::adas::run(sdtm, adam, r);
    crate::outcomes::run(sdtm, adam, r);
    crate::ophtha::run(sdtm, adam, r);
    crate::growth::run(sdtm, adam, r);
    crate::exacerbation::run(sdtm, adam, r);
    crate::vaccine::run(sdtm, adam, r);
    crate::escalation::run(sdtm, adam, r);
    crate::switching::run(sdtm, adam, r);
}

fn adsl(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(sl), Some(ex), Some(ds), Some(dm)) = (
        get(adam, "ADSL"),
        get(sdtm, "EX"),
        get(sdtm, "DS"),
        get(sdtm, "DM"),
    ) else {
        return;
    };
    // Exposure span from EX.
    let mut span: HashMap<&str, (Date, Date)> = HashMap::new();
    if let (Some(u), Some(a), Some(b)) = (ex.c("USUBJID"), ex.d("EXSTDTC"), ex.d("EXENDTC")) {
        for i in 0..ex.n() {
            if let (Some(u), Some(a), Some(b)) = (u[i].as_deref(), a[i], b[i]) {
                let e = span.entry(u).or_insert((a, b));
                e.0 = e.0.min(a);
                e.1 = e.1.max(b);
            }
        }
    }
    // End-of-study and randomisation from DS.
    let mut eos: HashMap<&str, &str> = HashMap::new();
    let mut rand: HashMap<&str, Date> = HashMap::new();
    if let (Some(u), Some(d), Some(s), Some(dt)) = (
        ds.c("USUBJID"),
        ds.c("DSDECOD"),
        ds.c("DSSCAT"),
        ds.d("DSSTDTC"),
    ) {
        for i in 0..ds.n() {
            let Some(x) = u[i].as_deref() else { continue };
            if d[i].as_deref() == Some("RANDOMIZED") {
                if let Some(dd) = dt[i] {
                    rand.insert(x, dd);
                }
            }
            if s[i].as_deref() == Some("END OF STUDY") {
                if let Some(dd) = d[i].as_deref() {
                    eos.insert(x, dd);
                }
            }
        }
    }
    let dm_dth: HashMap<&str, Option<Date>> = match (dm.c("USUBJID"), dm.d("DTHDTC")) {
        (Some(u), Some(d)) => u
            .iter()
            .zip(d)
            .filter_map(|(u, d)| Some((u.as_deref()?, d)))
            .collect(),
        _ => HashMap::new(),
    };
    let u = sl.c("USUBJID").unwrap_or_default();
    let (sdt, edt, rdt, dth) = (
        sl.d("TRTSDT"),
        sl.d("TRTEDT"),
        sl.d("RANDDT"),
        sl.d("DTHDT"),
    );
    let (saffl, eosstt, dcsreas) = (sl.c("SAFFL"), sl.c("EOSSTT"), sl.c("DCSREAS"));
    let dur = sl.f("TRTDURD");
    for i in 0..sl.n() {
        let Some(x) = u[i].as_deref() else { continue };
        let sp = span.get(x);
        if sdt.as_ref().map(|v| v[i]) != Some(sp.map(|s| s.0))
            || edt.as_ref().map(|v| v[i]) != Some(sp.map(|s| s.1))
        {
            r.add(
                L,
                "L2.ADSL",
                "ADSL",
                Some("TRTSDT"),
                Some(i),
                "TRTSDT/TRTEDT do not match EX".into(),
            );
        }
        if let Some(d) = &dur {
            let expected = sp.map(|s| f64::from(s.1.days_since(s.0) + 1));
            if d[i] != expected {
                r.add(
                    L,
                    "L2.ADSL",
                    "ADSL",
                    Some("TRTDURD"),
                    Some(i),
                    "TRTDURD does not match EX".into(),
                );
            }
        }
        if let Some(f) = saffl {
            if f[i].as_deref() != Some(if sp.is_some() { "Y" } else { "N" }) {
                r.add(
                    L,
                    "L2.ADSL",
                    "ADSL",
                    Some("SAFFL"),
                    Some(i),
                    "SAFFL does not match EX".into(),
                );
            }
        }
        if rdt.as_ref().map(|v| v[i]) != Some(rand.get(x).copied()) {
            r.add(
                L,
                "L2.ADSL",
                "ADSL",
                Some("RANDDT"),
                Some(i),
                "RANDDT does not match DS".into(),
            );
        }
        if dth.as_ref().map(|v| v[i]) != Some(dm_dth.get(x).copied().flatten()) {
            r.add(
                L,
                "L2.ADSL",
                "ADSL",
                Some("DTHDT"),
                Some(i),
                "DTHDT does not match DM".into(),
            );
        }
        if let (Some(st), Some(re)) = (eosstt, dcsreas) {
            let expected = match eos.get(x) {
                None => ("ONGOING", None),
                Some(&"COMPLETED") => ("COMPLETED", None),
                Some(d) => ("DISCONTINUED", Some(*d)),
            };
            if st[i].as_deref() != Some(expected.0) || re[i].as_deref() != expected.1 {
                r.add(
                    L,
                    "L2.ADSL",
                    "ADSL",
                    Some("EOSSTT"),
                    Some(i),
                    "EOSSTT/DCSREAS do not match DS".into(),
                );
            }
        }
    }
}

fn adae(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ad), Some(ae), Some(sl)) = (get(adam, "ADAE"), get(sdtm, "AE"), get(adam, "ADSL"))
    else {
        return;
    };
    if ad.n() != ae.n() {
        r.add(
            L,
            "L2.ADAE",
            "ADAE",
            None,
            None,
            format!("{} ADAE records for {} AE records", ad.n(), ae.n()),
        );
        return;
    }
    let sl_idx = sl.by_subject();
    let (sdt, edt) = (
        sl.d("TRTSDT").unwrap_or_default(),
        sl.d("TRTEDT").unwrap_or_default(),
    );
    // AE by (USUBJID, AESEQ).
    let mut src: HashMap<(&str, i64), Span> = HashMap::new();
    if let (Some(u), Some(s), Some(a), Some(b)) = (
        ae.c("USUBJID"),
        ae.i("AESEQ"),
        ae.d("AESTDTC"),
        ae.d("AEENDTC"),
    ) {
        for i in 0..ae.n() {
            if let (Some(u), Some(s)) = (u[i].as_deref(), s[i]) {
                src.insert((u, s), (a[i], b[i]));
            }
        }
    }
    let (u, seq) = (
        ad.c("USUBJID").unwrap_or_default(),
        ad.i("AESEQ").unwrap_or_default(),
    );
    let (ast, aen, fl, durn) = (
        ad.d("ASTDT").unwrap_or_default(),
        ad.d("AENDT").unwrap_or_default(),
        ad.c("TRTEMFL").unwrap_or_default(),
        ad.f("ADURN").unwrap_or_default(),
    );
    for i in 0..ad.n() {
        let (Some(x), Some(s)) = (u[i].as_deref(), seq[i]) else {
            continue;
        };
        let Some((a, b)) = src.get(&(x, s)) else {
            r.add(
                L,
                "L2.ADAE",
                "ADAE",
                Some("AESEQ"),
                Some(i),
                "no matching AE record".into(),
            );
            continue;
        };
        if ast[i] != *a || aen[i] != *b {
            r.add(
                L,
                "L2.ADAE",
                "ADAE",
                Some("ASTDT"),
                Some(i),
                "ASTDT/AENDT do not match AE".into(),
            );
        }
        let slr = sl_idx.get(x).copied();
        let em = match (a, slr.and_then(|k| sdt[k]), slr.and_then(|k| edt[k])) {
            (Some(a), Some(s), Some(e)) => *a >= s && *a <= e.add_days(TEAE_WINDOW_DAYS),
            _ => false,
        };
        if (fl[i].as_deref() == Some("Y")) != em {
            r.add(
                L,
                "L2.ADAE",
                "ADAE",
                Some("TRTEMFL"),
                Some(i),
                "TRTEMFL does not match the dates".into(),
            );
        }
        let expected = match (a, b) {
            (Some(a), Some(b)) => Some(f64::from(b.days_since(*a) + 1)),
            _ => None,
        };
        if durn[i] != expected {
            r.add(
                L,
                "L2.ADAE",
                "ADAE",
                Some("ADURN"),
                Some(i),
                "ADURN does not match the dates".into(),
            );
        }
    }
}

fn basic(b: &T<'_>, s: &T<'_>, src: &str, adam: &[Dataset], r: &mut Report) {
    let name = b.name().to_string();
    let Some(sl) = get(adam, "ADSL") else { return };
    let sl_idx = sl.by_subject();
    let sdt = sl.d("TRTSDT").unwrap_or_default();

    // Source records by (USUBJID, --SEQ).
    let mut source: HashMap<(&str, i64), SourceFinding<'_>> = HashMap::new();
    let stest = s.c(&format!("{src}TESTCD")).unwrap_or_default();
    let (su, sseq, sval, sdtc) = (
        s.c("USUBJID").unwrap_or_default(),
        s.i(&format!("{src}SEQ")).unwrap_or_default(),
        s.f(&format!("{src}STRESN")).unwrap_or_default(),
        s.d(&format!("{src}DTC")).unwrap_or_default(),
    );
    let lo = s
        .f(&format!("{src}STNRLO"))
        .unwrap_or_else(|| vec![None; s.n()]);
    let hi = s
        .f(&format!("{src}STNRHI"))
        .unwrap_or_else(|| vec![None; s.n()]);
    for i in 0..s.n() {
        if let (Some(u), Some(q)) = (su[i].as_deref(), sseq[i]) {
            source.insert(
                (u, q),
                (sval[i], sdtc[i], lo[i], hi[i], stest[i].as_deref()),
            );
        }
    }

    let u = b.c("USUBJID").unwrap_or_default();
    let pc = b.c("PARAMCD").unwrap_or_default();
    let srcseq = b.i("SRCSEQ").unwrap_or_default();
    let (aval, base, chg) = (
        b.f("AVAL").unwrap_or_default(),
        b.f("BASE").unwrap_or_default(),
        b.f("CHG").unwrap_or_default(),
    );
    let adt = b.d("ADT").unwrap_or_default();
    let (ablfl, anl, avn) = (
        b.c("ABLFL").unwrap_or_default(),
        b.c("ANL01FL").unwrap_or_default(),
        b.f("AVISITN").unwrap_or_default(),
    );
    let anrind = b.c("ANRIND");

    // Traceability, AVAL, ADT, ANRIND.
    for i in 0..b.n() {
        let (Some(x), Some(q)) = (u[i].as_deref(), srcseq[i]) else {
            r.add(
                L,
                "L2.TRACE",
                &name,
                Some("SRCSEQ"),
                Some(i),
                "missing traceability".into(),
            );
            continue;
        };
        let Some((v, d, lo, hi, testcd)) = source.get(&(x, q)) else {
            r.add(
                L,
                "L2.TRACE",
                &name,
                Some("SRCSEQ"),
                Some(i),
                format!("{src}SEQ {q} does not exist for {x}"),
            );
            continue;
        };
        // Parameters derived from the record (growth indicators) are checked by `growth`.
        let derived = pc[i].as_deref() != *testcd;
        if (aval[i] != *v && !derived) || adt[i] != *d {
            r.add(
                L,
                "L2.BDS",
                &name,
                Some("AVAL"),
                Some(i),
                "AVAL/ADT do not match the source".into(),
            );
        }
        if let Some(ind) = anrind {
            let expected = match (v, lo, hi) {
                (Some(v), Some(lo), Some(hi)) => Some(if v < lo {
                    "LOW"
                } else if v > hi {
                    "HIGH"
                } else {
                    "NORMAL"
                }),
                _ => None,
            };
            if ind[i].as_deref() != expected {
                r.add(
                    L,
                    "L2.BDS",
                    &name,
                    Some("ANRIND"),
                    Some(i),
                    "ANRIND does not match the reference range".into(),
                );
            }
        }
    }

    // Baseline: the last record on or before TRTSDT per subject and parameter.
    let mut groups: HashMap<(&str, &str), Vec<usize>> = HashMap::new();
    for i in 0..b.n() {
        if let (Some(x), Some(p)) = (u[i].as_deref(), pc[i].as_deref()) {
            groups.entry((x, p)).or_default().push(i);
        }
    }
    let mut anl_seen: HashMap<(&str, &str, i64), usize> = HashMap::new();
    for ((x, p), rows) in &groups {
        let trt = sl_idx.get(x).and_then(|&k| sdt[k]);
        let base_row = rows
            .iter()
            .copied()
            .filter(|&i| matches!((adt[i], trt), (Some(a), Some(t)) if a <= t) && aval[i].is_some())
            .max_by_key(|&i| (adt[i], srcseq[i]));
        for &i in rows {
            if (ablfl[i].as_deref() == Some("Y")) != (Some(i) == base_row) {
                r.add(
                    L,
                    "L2.BDS",
                    &name,
                    Some("ABLFL"),
                    Some(i),
                    "ABLFL is not the last pre-treatment record".into(),
                );
            }
            let bv = base_row.and_then(|k| aval[k]);
            if base[i] != bv {
                r.add(
                    L,
                    "L2.BDS",
                    &name,
                    Some("BASE"),
                    Some(i),
                    "BASE does not equal the baseline AVAL".into(),
                );
            }
            let post = matches!((adt[i], trt), (Some(a), Some(t)) if a > t);
            let expected = match (post, aval[i], bv) {
                (true, Some(a), Some(bv)) => Some(a - bv),
                _ => None,
            };
            let ok = match (chg[i], expected) {
                (Some(c), Some(e)) => (c - e).abs() < TOL,
                (None, None) => true,
                _ => false,
            };
            if !ok {
                r.add(
                    L,
                    "L2.BDS",
                    &name,
                    Some("CHG"),
                    Some(i),
                    "CHG does not equal AVAL - BASE".into(),
                );
            }
            if anl[i].as_deref() == Some("Y") {
                if let Some(v) = avn[i] {
                    if let Some(prev) = anl_seen.insert((x, p, v as i64), i) {
                        r.add(
                            L,
                            "L2.BDS",
                            &name,
                            Some("ANL01FL"),
                            Some(i),
                            format!("second ANL01FL record (first at row {})", prev + 1),
                        );
                    }
                }
            }
        }
    }
}
