//! Dose-finding checks.
//!
//! - L1 `L1.DLT`: an AE flagged as a DLT (SUPPAE `AEDLT` = 'Y') starts within the 28-day DLT
//!   window from the first dose, and only AEs in the window carry the qualifier.
//! - L2 `L2.DLT`: ADAE `DLTFL` is the AE's `AEDLT`, and ADSL `DLTEVLFL` is 'Y' exactly for
//!   subjects with a DLT or with participation through the window.

use std::collections::{HashMap, HashSet};

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ae), Some(sup), Some(dm)) = (get(sdtm, "AE"), get(sdtm, "SUPPAE"), get(sdtm, "DM"))
    else {
        return;
    };
    let (Some(su), Some(q), Some(idv), Some(val)) = (
        sup.c("USUBJID"),
        sup.c("QNAM"),
        sup.c("IDVARVAL"),
        sup.c("QVAL"),
    ) else {
        return;
    };
    let mut dlt: HashMap<(&str, i64), &str> = HashMap::new();
    for i in 0..sup.n() {
        if let (Some(x), Some("AEDLT"), Some(v), Some(f)) = (
            su[i].as_deref(),
            q[i].as_deref(),
            idv[i].as_deref(),
            val[i].as_deref(),
        ) {
            if let Ok(seq) = v.parse() {
                dlt.insert((x, seq), f);
            }
        }
    }
    if dlt.is_empty() {
        return;
    }
    let first: HashMap<&str, (Option<Date>, Option<Date>)> =
        match (dm.c("USUBJID"), dm.d("RFXSTDTC"), dm.d("RFPENDTC")) {
            (Some(u), Some(a), Some(b)) => (0..dm.n())
                .filter_map(|i| Some((u[i].as_deref()?, (a[i], b[i]))))
                .collect(),
            _ => return,
        };
    let (Some(au), Some(aseq), Some(ast)) = (ae.c("USUBJID"), ae.i("AESEQ"), ae.d("AESTDTC"))
    else {
        return;
    };
    let mut with_dlt: HashSet<&str> = HashSet::new();
    for i in 0..ae.n() {
        let (Some(x), Some(s), Some(d)) = (au[i].as_deref(), aseq[i], ast[i]) else {
            continue;
        };
        let Some(&(Some(f), _)) = first.get(x) else {
            continue;
        };
        let in_window = d >= f && d.days_since(f) < 28;
        match dlt.get(&(x, s)) {
            Some(flag) => {
                if !in_window {
                    r.add(
                        Level::L1,
                        "L1.DLT",
                        "SUPPAE",
                        Some("QVAL"),
                        None,
                        format!("{x} AESEQ {s}: AEDLT outside the DLT window"),
                    );
                }
                if *flag == "Y" {
                    with_dlt.insert(x);
                }
            }
            None if in_window => {
                r.add(
                    Level::L1,
                    "L1.DLT",
                    "SUPPAE",
                    Some("QNAM"),
                    None,
                    format!("{x} AESEQ {s}: AE in the DLT window without AEDLT"),
                );
            }
            None => {}
        }
    }
    if let Some(ad) = get(adam, "ADAE") {
        if let (Some(u), Some(s), Some(fl)) = (ad.c("USUBJID"), ad.i("AESEQ"), ad.c("DLTFL")) {
            for i in 0..ad.n() {
                let want = u[i].as_deref().zip(s[i]).and_then(|k| dlt.get(&k).copied());
                if fl[i].as_deref() != want {
                    r.add(
                        Level::L2,
                        "L2.DLT",
                        "ADAE",
                        Some("DLTFL"),
                        Some(i),
                        format!("DLTFL should be {want:?}"),
                    );
                }
            }
        }
    }
    if let Some(sl) = get(adam, "ADSL") {
        let (Some(u), Some(fl)) = (sl.c("USUBJID"), sl.c("DLTEVLFL")) else {
            r.add(
                Level::L2,
                "L2.DLT",
                "ADSL",
                Some("DLTEVLFL"),
                None,
                "dose-finding ADSL lacks DLTEVLFL".into(),
            );
            return;
        };
        for i in 0..sl.n() {
            let Some(x) = u[i].as_deref() else { continue };
            let full = first.get(x).is_some_and(
                |(a, b)| matches!((a, b), (Some(a), Some(b)) if b.days_since(*a) >= 27),
            );
            let want = if full || with_dlt.contains(x) {
                "Y"
            } else {
                "N"
            };
            if fl[i].as_deref() != Some(want) {
                r.add(
                    Level::L2,
                    "L2.DLT",
                    "ADSL",
                    Some("DLTEVLFL"),
                    Some(i),
                    format!("{x}: DLTEVLFL should be {want}"),
                );
            }
        }
    }
}
