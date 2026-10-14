//! L2 `L2.ADTTE`: for an outcomes study, every endpoint in ADTTE is recomputed from CE and
//! ADSL: the first qualifying CE event (or `DTHDT` for all-cause death) is the event; otherwise
//! the subject is censored at the later of `LSTALVDT` and `RANDDT`.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L2;

/// (RANDDT, LSTALVDT, DTHDT) of a subject.
type Dates = (Option<Date>, Option<Date>, Option<Date>);

/// PARAMCD → the CE Preferred Terms that are its events (empty: all-cause death).
const ENDPOINTS: [(&str, &[&str]); 5] = [
    (
        "MACE",
        &["Cardiac death", "Myocardial infarction", "Ischaemic stroke"],
    ),
    ("CVDEATH", &["Cardiac death"]),
    ("MI", &["Myocardial infarction"]),
    ("STROKE", &["Ischaemic stroke"]),
    ("DEATH", &[]),
];

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ce), Some(sl), Some(tte)) = (get(sdtm, "CE"), get(adam, "ADSL"), get(adam, "ADTTE"))
    else {
        return;
    };
    let (Some(u), Some(decod), Some(dtc), Some(cat)) = (
        ce.c("USUBJID"),
        ce.c("CEDECOD"),
        ce.d("CESTDTC"),
        ce.c("CECAT"),
    ) else {
        return;
    };
    let adjudicated = |i: usize| cat[i].as_deref() == Some("ADJUDICATED ENDPOINT EVENT");
    if !(0..ce.n()).any(adjudicated) {
        return;
    }
    let mut first: HashMap<(&str, &str), Date> = HashMap::new();
    for i in (0..ce.n()).filter(|&i| adjudicated(i)) {
        if let (Some(x), Some(pt), Some(d)) = (u[i].as_deref(), decod[i].as_deref(), dtc[i]) {
            let e = first.entry((x, pt)).or_insert(d);
            *e = (*e).min(d);
        }
    }
    let (Some(su), Some(rand), Some(alive), Some(dth)) = (
        sl.c("USUBJID"),
        sl.d("RANDDT"),
        sl.d("LSTALVDT"),
        sl.d("DTHDT"),
    ) else {
        return;
    };
    let subj: HashMap<&str, Dates> = (0..sl.n())
        .filter_map(|i| Some((su[i].as_deref()?, (rand[i], alive[i], dth[i]))))
        .collect();
    let (Some(tu), Some(pc), Some(adt), Some(cnsr), Some(aval)) = (
        tte.c("USUBJID"),
        tte.c("PARAMCD"),
        tte.d("ADT"),
        tte.i("CNSR"),
        tte.f("AVAL"),
    ) else {
        return;
    };
    let mut seen = 0;
    for i in 0..tte.n() {
        let (Some(x), Some(p)) = (tu[i].as_deref(), pc[i].as_deref()) else {
            continue;
        };
        let Some((_, pts)) = ENDPOINTS.iter().find(|(cd, _)| *cd == p) else {
            continue;
        };
        let Some(&(Some(rd), alive, death)) = subj.get(x) else {
            continue;
        };
        seen += 1;
        let event = if pts.is_empty() {
            death
        } else {
            pts.iter()
                .filter_map(|pt| first.get(&(x, *pt)).copied())
                .min()
        };
        let (date, censored) = match event {
            Some(d) => (d, 0),
            None => (alive.map_or(rd, |a| a.max(rd)), 1),
        };
        let date = date.max(rd);
        let want_aval = f64::from(date.days_since(rd) + 1);
        if adt[i] != Some(date) || cnsr[i] != Some(censored) || aval[i] != Some(want_aval) {
            r.add(
                L,
                "L2.ADTTE",
                "ADTTE",
                Some("ADT"),
                Some(i),
                format!(
                    "{p} for {x}: expected ADT {date}, CNSR {censored}, found {:?}, {:?}",
                    adt[i], cnsr[i]
                ),
            );
        }
    }
    if seen != ENDPOINTS.len() * subj.values().filter(|s| s.0.is_some()).count() {
        r.add(
            L,
            "L2.ADTTE",
            "ADTTE",
            Some("PARAMCD"),
            None,
            "an endpoint record is missing".into(),
        );
    }
}
