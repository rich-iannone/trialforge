//! L2 `L2.ADADAS`: the ADAS-Cog total (`ACTOT`) of each visit is recomputed from the FT items
//! (prorated sum × 11 / n with 8 to 10 items, missing with fewer), every item record matches
//! FT, `CHG` is `AVAL - BASE`, and a `LOCF` record repeats the previous analysis value.

use std::collections::{BTreeMap, HashMap};

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L2;
const TOL: f64 = 1e-9;

/// FT items by (USUBJID, FTSEQ): (FTSTRESN, FTTESTCD).
type Items<'a> = HashMap<(&'a str, i64), (Option<f64>, Option<&'a str>)>;
/// Analysis records by (USUBJID, PARAMCD): (AVISITN, AVAL, LOCF).
type Analysis<'a> = BTreeMap<(&'a str, &'a str), Vec<(f64, f64, bool)>>;

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ft), Some(ad)) = (get(sdtm, "FT"), get(adam, "ADADAS")) else {
        return;
    };
    let (Some(fu), Some(fv), Some(fd), Some(fseq), Some(fcd)) = (
        ft.c("USUBJID"),
        ft.f("FTSTRESN"),
        ft.d("FTDTC"),
        ft.i("FTSEQ"),
        ft.c("FTTESTCD"),
    ) else {
        return;
    };
    let mut items: Items<'_> = HashMap::new();
    let mut totals: HashMap<(&str, Date), (f64, usize)> = HashMap::new();
    for i in 0..ft.n() {
        let Some(x) = fu[i].as_deref() else { continue };
        if let Some(q) = fseq[i] {
            items.insert((x, q), (fv[i], fcd[i].as_deref()));
        }
        if let (Some(v), Some(d)) = (fv[i], fd[i]) {
            let e = totals.entry((x, d)).or_insert((0.0, 0));
            e.0 += v;
            e.1 += 1;
        }
    }
    let (
        Some(u),
        Some(pc),
        Some(aval),
        Some(adt),
        Some(seq),
        Some(base),
        Some(chg),
        Some(avn),
        Some(dtype),
        Some(anl),
    ) = (
        ad.c("USUBJID"),
        ad.c("PARAMCD"),
        ad.f("AVAL"),
        ad.d("ADT"),
        ad.i("SRCSEQ"),
        ad.f("BASE"),
        ad.f("CHG"),
        ad.f("AVISITN"),
        ad.c("DTYPE"),
        ad.c("ANL01FL"),
    )
    else {
        return;
    };
    // Analysis values per subject and parameter by visit, for LOCF.
    let mut analysis: Analysis<'_> = BTreeMap::new();
    for i in 0..ad.n() {
        let (Some(x), Some(p)) = (u[i].as_deref(), pc[i].as_deref()) else {
            continue;
        };
        let locf = dtype[i].as_deref() == Some("LOCF");
        if p == "ACTOT" {
            if !locf {
                let want = adt[i]
                    .and_then(|d| totals.get(&(x, d)))
                    .and_then(|&(s, n)| match n {
                        11 => Some(s),
                        8..=10 => Some(s * 11.0 / n as f64),
                        _ => None,
                    });
                let ok = matches!((aval[i], want), (Some(a), Some(w)) if (a - w).abs() < TOL);
                if !ok {
                    r.add(
                        L,
                        "L2.ADADAS",
                        "ADADAS",
                        Some("AVAL"),
                        Some(i),
                        format!("ACTOT should be {want:?}"),
                    );
                }
            }
        } else {
            let src = seq[i].and_then(|q| items.get(&(x, q)));
            if src.is_none_or(|(v, cd)| *v != aval[i] || *cd != Some(p)) {
                r.add(
                    L,
                    "L2.ADADAS",
                    "ADADAS",
                    Some("AVAL"),
                    Some(i),
                    "item does not match FT".into(),
                );
            }
        }
        let want_chg = match (avn[i], aval[i], base[i]) {
            (Some(w), Some(a), Some(b)) if w > 0.0 => Some(a - b),
            _ => None,
        };
        let ok = match (chg[i], want_chg) {
            (Some(c), Some(w)) => (c - w).abs() < TOL,
            (None, None) => true,
            _ => false,
        };
        if !ok {
            r.add(
                L,
                "L2.ADADAS",
                "ADADAS",
                Some("CHG"),
                Some(i),
                "CHG is not AVAL - BASE".into(),
            );
        }
        if anl[i].as_deref() == Some("Y") {
            if let (Some(w), Some(a)) = (avn[i], aval[i]) {
                analysis.entry((x, p)).or_default().push((w, a, locf));
            }
        }
    }
    for ((x, p), mut v) in analysis {
        v.sort_by(|a, b| a.0.total_cmp(&b.0));
        for k in 1..v.len() {
            if v[k].2 && (v[k].1 - v[k - 1].1).abs() > TOL {
                r.add(
                    L,
                    "L2.ADADAS",
                    "ADADAS",
                    Some("DTYPE"),
                    None,
                    format!(
                        "{x} {p}: LOCF at {} does not carry the previous value",
                        v[k].0
                    ),
                );
            }
            if v[k].0 == v[k - 1].0 {
                r.add(
                    L,
                    "L2.ADADAS",
                    "ADADAS",
                    Some("ANL01FL"),
                    None,
                    format!("{x} {p}: two analysis records at {}", v[k].0),
                );
            }
        }
    }
}
