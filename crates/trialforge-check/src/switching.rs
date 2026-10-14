//! L2 `L2.CROSSOVER`: for crossover on progression, ADSL `CROSSDT` is the first EX
//! administration of a treatment other than the subject's first (the switch), and `CROSSFL`
//! says whether there is one.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ex), Some(sl)) = (get(sdtm, "EX"), get(adam, "ADSL")) else {
        return;
    };
    let (Some(su), Some(fl), Some(dt)) = (sl.c("USUBJID"), sl.c("CROSSFL"), sl.d("CROSSDT")) else {
        return;
    };
    let (Some(u), Some(trt), Some(st)) = (ex.c("USUBJID"), ex.c("EXTRT"), ex.d("EXSTDTC")) else {
        return;
    };
    let mut first: HashMap<&str, (Date, &str)> = HashMap::new();
    let mut rows: Vec<(&str, Date, &str)> = (0..ex.n())
        .filter_map(|i| Some((u[i].as_deref()?, st[i]?, trt[i].as_deref()?)))
        .collect();
    rows.sort();
    let mut switch: HashMap<&str, Date> = HashMap::new();
    for (x, d, t) in rows {
        let f = first.entry(x).or_insert((d, t));
        if f.1 != t {
            switch.entry(x).or_insert(d);
        }
    }
    for i in 0..sl.n() {
        let Some(x) = su[i].as_deref() else { continue };
        let want = switch.get(x).copied();
        let want_fl = if want.is_some() { "Y" } else { "N" };
        if dt[i] != want || fl[i].as_deref() != Some(want_fl) {
            r.add(
                Level::L2,
                "L2.CROSSOVER",
                "ADSL",
                Some("CROSSDT"),
                Some(i),
                format!("{x}: CROSSDT should be {want:?}"),
            );
        }
    }
}
