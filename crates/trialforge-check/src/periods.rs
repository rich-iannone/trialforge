//! L2 `L2.PERIOD`: crossover period variables. ADSL `TRTkkP`/`TRTkkA` are the TA elements of
//! the subject's arm in `TREATMENT k`, `APkkSDT`/`APkkEDT` the subject's SE dates for it; every
//! analysis record's `APERIOD` is the period whose dates contain the record's date, and its
//! `TRTA` the actual treatment of the latest period started on or before that date.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date};

use crate::table::{T, get};
use crate::{Level, Report};

const L: Level = Level::L2;

/// (start, end) of a subject's element.
type Span = (Option<Date>, Option<Date>);
/// A subject's period: (period, start, end, actual treatment).
type SubjectPeriod = (i64, Date, Date, Option<String>);

fn period(epoch: Option<&str>) -> Option<i64> {
    epoch?.strip_prefix("TREATMENT ")?.parse().ok()
}

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ta), Some(se), Some(dm), Some(sl)) = (
        get(sdtm, "TA"),
        get(sdtm, "SE"),
        get(sdtm, "DM"),
        get(adam, "ADSL"),
    ) else {
        return;
    };
    let (Some(ta_arm), Some(ta_ep), Some(ta_el)) = (ta.c("ARMCD"), ta.c("EPOCH"), ta.c("ELEMENT"))
    else {
        return;
    };
    let mut element: HashMap<(&str, i64), &str> = HashMap::new();
    for i in 0..ta.n() {
        if let (Some(a), Some(k), Some(e)) = (
            ta_arm[i].as_deref(),
            period(ta_ep[i].as_deref()),
            ta_el[i].as_deref(),
        ) {
            element.insert((a, k), e);
        }
    }
    let n_periods = element.keys().map(|k| k.1).max().unwrap_or(0);
    if n_periods == 0 {
        return;
    }
    let mut spans: HashMap<(&str, i64), Span> = HashMap::new();
    if let (Some(u), Some(ep), Some(st), Some(en)) = (
        se.c("USUBJID"),
        se.c("EPOCH"),
        se.d("SESTDTC"),
        se.d("SEENDTC"),
    ) {
        for i in 0..se.n() {
            if let (Some(u), Some(k)) = (u[i].as_deref(), period(ep[i].as_deref())) {
                spans.insert((u, k), (st[i], en[i]));
            }
        }
    }
    let dm_idx = dm.by_subject();
    let (Some(armcd), Some(actarmcd)) = (dm.c("ARMCD"), dm.c("ACTARMCD")) else {
        return;
    };
    let Some(u) = sl.c("USUBJID") else { return };
    // Per ADSL subject: (k, start, end, TRTkkA).
    let mut subject: HashMap<&str, Vec<SubjectPeriod>> = HashMap::new();
    for k in 1..=n_periods {
        let (p, a) = (sl.c(&format!("TRT{k:02}P")), sl.c(&format!("TRT{k:02}A")));
        let (sd, ed) = (sl.d(&format!("AP{k:02}SDT")), sl.d(&format!("AP{k:02}EDT")));
        let (Some(p), Some(a), Some(sd), Some(ed)) = (p, a, sd, ed) else {
            r.add(
                L,
                "L2.PERIOD",
                "ADSL",
                None,
                None,
                format!("period {k} variables are missing"),
            );
            continue;
        };
        for i in 0..sl.n() {
            let Some(x) = u[i].as_deref() else { continue };
            let Some(&d) = dm_idx.get(x) else { continue };
            let want =
                |cd: &Option<String>| cd.as_deref().and_then(|c| element.get(&(c, k))).copied();
            // Actual treatment only for periods the subject entered (SE).
            let actual = want(&actarmcd[d]).filter(|_| spans.contains_key(&(x, k)));
            if p[i].as_deref() != want(&armcd[d]) || a[i].as_deref() != actual {
                r.add(
                    L,
                    "L2.PERIOD",
                    "ADSL",
                    Some("TRT01P"),
                    Some(i),
                    format!("period {k} treatments do not match TA"),
                );
            }
            let span = spans.get(&(x, k)).copied().unwrap_or((None, None));
            if (sd[i], ed[i]) != span {
                r.add(
                    L,
                    "L2.PERIOD",
                    "ADSL",
                    Some("AP01SDT"),
                    Some(i),
                    format!("period {k} dates do not match SE"),
                );
            }
            if let (Some(s), Some(e)) = (sd[i], ed[i]) {
                subject.entry(x).or_default().push((k, s, e, a[i].clone()));
            }
        }
    }
    for ds in adam.iter().filter(|d| d.name != "ADSL") {
        let t = T::new(ds);
        let Some(ap) = t.i("APERIOD") else { continue };
        let date_var = if t.has("ASTDT") { "ASTDT" } else { "ADT" };
        let (Some(us), Some(dates)) = (t.c("USUBJID"), t.d(date_var)) else {
            continue;
        };
        let trta = t.c("TRTA");
        for i in 0..t.n() {
            let Some(ps) = us[i].as_deref().and_then(|x| subject.get(x)) else {
                continue;
            };
            let d = dates[i];
            let want = d
                .and_then(|d| ps.iter().find(|p| p.1 <= d && d <= p.2))
                .map(|p| p.0);
            if ap[i] != want {
                r.add(
                    L,
                    "L2.PERIOD",
                    t.name(),
                    Some("APERIOD"),
                    Some(i),
                    format!("APERIOD should be {want:?}"),
                );
            }
            if let Some(trta) = trta {
                let attributed = d
                    .and_then(|d| ps.iter().rev().find(|p| p.1 <= d))
                    .or(ps.first())
                    .and_then(|p| p.3.as_deref());
                if trta[i].as_deref() != attributed {
                    r.add(
                        L,
                        "L2.PERIOD",
                        t.name(),
                        Some("TRTA"),
                        Some(i),
                        "TRTA is not the attributed period's treatment".into(),
                    );
                }
            }
        }
    }
}
