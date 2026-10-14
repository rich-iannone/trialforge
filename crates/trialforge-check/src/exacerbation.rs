//! L2 `L2.ADEXAC`: for a study with exacerbations (CE `CECAT = EXACERBATION`), the counts
//! (all and severe) from first dose to `LSTALVDT`, the observation years and the annualized
//! rate in ADEXAC are recomputed from CE and ADSL, and ADTTE `TTFEXAC` is the first
//! exacerbation after randomisation, else censored at `LSTALVDT`.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L2;
const TOL: f64 = 1e-9;

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ce), Some(sl)) = (get(sdtm, "CE"), get(adam, "ADSL")) else {
        return;
    };
    let (Some(u), Some(cat), Some(sev), Some(dtc)) = (
        ce.c("USUBJID"),
        ce.c("CECAT"),
        ce.c("CESEV"),
        ce.d("CESTDTC"),
    ) else {
        return;
    };
    let mut onsets: HashMap<&str, Vec<(Date, bool)>> = HashMap::new();
    for i in 0..ce.n() {
        if let (Some(x), Some("EXACERBATION"), Some(d)) =
            (u[i].as_deref(), cat[i].as_deref(), dtc[i])
        {
            onsets
                .entry(x)
                .or_default()
                .push((d, sev[i].as_deref() == Some("SEVERE")));
        }
    }
    if onsets.is_empty() && get(adam, "ADEXAC").is_none() {
        return;
    }
    let (Some(su), Some(start), Some(end), Some(rand)) = (
        sl.c("USUBJID"),
        sl.d("TRTSDT"),
        sl.d("LSTALVDT"),
        sl.d("RANDDT"),
    ) else {
        return;
    };
    // Expected per subject and parameter.
    let mut expected: HashMap<(&str, &str), f64> = HashMap::new();
    let mut tte: HashMap<&str, (Date, i64)> = HashMap::new();
    for i in 0..sl.n() {
        let Some(x) = su[i].as_deref() else { continue };
        let events = onsets.get(x).map(Vec::as_slice).unwrap_or_default();
        if let (Some(s), Some(e)) = (start[i], end[i]) {
            let within: Vec<&(Date, bool)> =
                events.iter().filter(|o| o.0 >= s && o.0 <= e).collect();
            let n = within.len() as f64;
            let years = f64::from(e.days_since(s) + 1) / 365.25;
            expected.insert((x, "EXACN"), n);
            expected.insert(
                (x, "SEVEXACN"),
                within.iter().filter(|o| o.1).count() as f64,
            );
            expected.insert((x, "FUPYR"), years);
            expected.insert((x, "EXACRATE"), n / years);
        }
        if let Some(rd) = rand[i] {
            let first = events.iter().map(|o| o.0).min();
            let (date, cnsr) = match first {
                Some(d) => (d.max(rd), 0),
                None => (end[i].map_or(rd, |e| e.max(rd)), 1),
            };
            tte.insert(x, (date, cnsr));
        }
    }
    let Some(ad) = get(adam, "ADEXAC") else {
        r.add(
            L,
            "L2.ADEXAC",
            "ADEXAC",
            None,
            None,
            "exacerbations in CE but no ADEXAC".into(),
        );
        return;
    };
    let (Some(au), Some(pc), Some(aval)) = (ad.c("USUBJID"), ad.c("PARAMCD"), ad.f("AVAL")) else {
        return;
    };
    let mut seen = 0;
    for i in 0..ad.n() {
        let (Some(x), Some(p)) = (au[i].as_deref(), pc[i].as_deref()) else {
            continue;
        };
        let want = expected.get(&(x, p)).copied();
        seen += usize::from(want.is_some());
        if !matches!((want, aval[i]), (Some(w), Some(a)) if (w - a).abs() < TOL) {
            r.add(
                L,
                "L2.ADEXAC",
                "ADEXAC",
                Some("AVAL"),
                Some(i),
                format!("{x} {p}: expected {want:?}, found {:?}", aval[i]),
            );
        }
    }
    if seen != expected.len() {
        r.add(
            L,
            "L2.ADEXAC",
            "ADEXAC",
            Some("PARAMCD"),
            None,
            format!("{seen} ADEXAC records for {} expected", expected.len()),
        );
    }
    let Some(t) = get(adam, "ADTTE") else { return };
    let (Some(tu), Some(tp), Some(adt), Some(cnsr)) =
        (t.c("USUBJID"), t.c("PARAMCD"), t.d("ADT"), t.i("CNSR"))
    else {
        return;
    };
    for i in 0..t.n() {
        let (Some(x), Some("TTFEXAC")) = (tu[i].as_deref(), tp[i].as_deref()) else {
            continue;
        };
        let want = tte.get(x).copied();
        if want != adt[i].zip(cnsr[i]) {
            r.add(
                L,
                "L2.ADEXAC",
                "ADTTE",
                Some("ADT"),
                Some(i),
                format!("{x} TTFEXAC: expected {want:?}"),
            );
        }
    }
}
