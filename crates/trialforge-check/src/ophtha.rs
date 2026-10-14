//! L2 `L2.ADBCVA`: every visual-acuity record in ADBCVA is recomputed from OE and the study
//! eye in SC: the parameter and `AFEYE` (study or fellow eye, letters or logMAR), `AVAL`, the
//! baseline (last record on or before the first dose), `CHG`, the analysis flag (baseline and
//! on-treatment records) and the Snellen category and criterion flags of the letters score.
//! Every OE visual-acuity record must appear in both scales.

use std::collections::HashMap;

use trialforge_core::{Dataset, Date};

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L2;
const TOL: f64 = 1e-9;
/// Assessments up to this many days after the last injection are on treatment.
const WINDOW: i32 = 112;

/// The Snellen denominator of an ETDRS letters score.
fn snellen(letters: f64) -> f64 {
    let lines = [
        1000.0, 800.0, 640.0, 500.0, 400.0, 320.0, 250.0, 200.0, 160.0, 125.0, 100.0, 80.0, 63.0,
        50.0, 40.0, 32.0, 25.0, 20.0, 16.0, 12.0,
    ];
    if letters >= 98.0 {
        1.0
    } else if letters < 4.0 {
        1000.0
    } else {
        // 4–8 letters is 20/800, then one line per 5 letters (94–97 is 20/12).
        lines[(((letters - 4.0) / 5.0).floor() as usize + 1).min(19)]
    }
}

/// The criterion flags on a letters change, in CRIT1–CRIT8 order.
fn criteria(c: f64) -> [bool; 8] {
    [
        (0.0..=5.0).contains(&c),
        (-5.0..=-1.0).contains(&c),
        (10.0..=15.0).contains(&c),
        c <= -20.0,
        c <= 5.0,
        c <= 10.0,
        c >= -15.0,
        c >= 15.0,
    ]
}

/// A source assessment: (letters, date, visit, visit number).
type Assessment<'a> = (f64, Option<Date>, Option<&'a str>, i64);

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(oe), Some(sc), Some(sl), Some(ad)) = (
        get(sdtm, "OE"),
        get(sdtm, "SC"),
        get(adam, "ADSL"),
        get(adam, "ADBCVA"),
    ) else {
        return;
    };
    // Study eye from SC.
    let mut eye: HashMap<&str, &str> = HashMap::new();
    if let (Some(u), Some(cd), Some(res)) = (sc.c("USUBJID"), sc.c("SCTESTCD"), sc.c("SCSTRESC")) {
        for i in 0..sc.n() {
            if let (Some(x), Some("FOCID"), Some(f)) =
                (u[i].as_deref(), cd[i].as_deref(), res[i].as_deref())
            {
                eye.insert(x, if f == "OS" { "LEFT" } else { "RIGHT" });
            }
        }
    }
    let (Some(su), Some(sdt), Some(edt)) = (sl.c("USUBJID"), sl.d("TRTSDT"), sl.d("TRTEDT")) else {
        return;
    };
    let dates: HashMap<&str, (Option<Date>, Option<Date>)> = (0..sl.n())
        .filter_map(|i| Some((su[i].as_deref()?, (sdt[i], edt[i]))))
        .collect();
    // Source assessments by (subject, parameter): letters, date, visit, visit number.
    let mut source: HashMap<(&str, &str), Vec<Assessment<'_>>> = HashMap::new();
    let (Some(u), Some(cd), Some(v), Some(dtc), Some(lat), Some(visit), Some(vn)) = (
        oe.c("USUBJID"),
        oe.c("OETESTCD"),
        oe.f("OESTRESN"),
        oe.d("OEDTC"),
        oe.c("OELAT"),
        oe.c("VISIT"),
        oe.i("VISITNUM"),
    ) else {
        return;
    };
    for i in 0..oe.n() {
        let (Some(x), Some("NUMLCOR"), Some(val), Some(l), Some(num)) = (
            u[i].as_deref(),
            cd[i].as_deref(),
            v[i],
            lat[i].as_deref(),
            vn[i],
        ) else {
            continue;
        };
        let Some(&se) = eye.get(x) else { continue };
        let study = l == se;
        let a = (val, dtc[i], visit[i].as_deref(), num);
        source
            .entry((x, if study { "SBCVA" } else { "FBCVA" }))
            .or_default()
            .push(a);
        let log = (
            trialforge_core::math::round_dp(1.7 - 0.02 * val, 2),
            dtc[i],
            visit[i].as_deref(),
            num,
        );
        source
            .entry((x, if study { "SBCVALOG" } else { "FBCVALOG" }))
            .or_default()
            .push(log);
    }
    // Expected values per (subject, parameter, visit number): (AVAL, ABLFL, BASE, ANL01FL).
    type Expected = (f64, bool, Option<f64>, bool);
    let mut expected: HashMap<(&str, &str, i64), Expected> = HashMap::new();
    for (&(x, p), recs) in &source {
        let Some(&(Some(start), end)) = dates.get(x) else {
            continue;
        };
        let base = recs
            .iter()
            .filter(|a| a.1.is_some_and(|d| d <= start))
            .max_by_key(|a| (a.1, a.3));
        for a in recs {
            let is_base = base.is_some_and(|b| b.3 == a.3);
            let on = a.1.is_some_and(|d| {
                (d == start && a.2 != Some("BASELINE"))
                    || (d > start && end.is_none_or(|e| d <= e.add_days(WINDOW)))
            });
            expected.insert(
                (x, p, a.3),
                (a.0, is_base, base.map(|b| b.0), is_base || on),
            );
        }
    }
    let (
        Some(au),
        Some(pc),
        Some(aval),
        Some(base),
        Some(chg),
        Some(abl),
        Some(anl),
        Some(avn),
        Some(cat),
    ) = (
        ad.c("USUBJID"),
        ad.c("PARAMCD"),
        ad.f("AVAL"),
        ad.f("BASE"),
        ad.f("CHG"),
        ad.c("ABLFL"),
        ad.c("ANL01FL"),
        ad.i("VISITNUM"),
        ad.f("AVALCA1N"),
    )
    else {
        return;
    };
    let afeye = ad.c("AFEYE");
    let crit: Vec<Option<&[Option<String>]>> =
        (1..=8).map(|k| ad.c(&format!("CRIT{k}FL"))).collect();
    let mut seen = 0;
    let close = |a: Option<f64>, b: Option<f64>| match (a, b) {
        (Some(a), Some(b)) => (a - b).abs() < TOL,
        (None, None) => true,
        _ => false,
    };
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(num)) = (au[i].as_deref(), pc[i].as_deref(), avn[i]) else {
            continue;
        };
        let Some(&(want, is_base, want_base, want_anl)) = expected.get(&(x, p, num)) else {
            r.add(
                L,
                "L2.ADBCVA",
                "ADBCVA",
                Some("PARAMCD"),
                Some(i),
                format!("{x} {p} at visit {num} has no OE source"),
            );
            continue;
        };
        seen += 1;
        let want_chg = want_base.map(|b| want - b);
        let flag = |b: bool| b.then_some("Y");
        let mut bad: Vec<&str> = Vec::new();
        if !close(aval[i], Some(want)) {
            bad.push("AVAL");
        }
        if !close(base[i], want_base) {
            bad.push("BASE");
        }
        if !close(chg[i], want_chg) {
            bad.push("CHG");
        }
        if abl[i].as_deref() != flag(is_base) {
            bad.push("ABLFL");
        }
        if anl[i].as_deref() != flag(want_anl) {
            bad.push("ANL01FL");
        }
        let want_eye = if p.starts_with('S') {
            "Study Eye"
        } else {
            "Fellow Eye"
        };
        if afeye.and_then(|a| a[i].as_deref()) != Some(want_eye) {
            bad.push("AFEYE");
        }
        let letters = matches!(p, "SBCVA" | "FBCVA");
        if !close(cat[i], letters.then(|| snellen(want))) {
            bad.push("AVALCA1N");
        }
        for (k, col) in crit.iter().enumerate() {
            let got = col.and_then(|c| c[i].as_deref());
            let exp = if letters {
                want_chg.map(|c| if criteria(c)[k] { "Y" } else { "N" })
            } else {
                None
            };
            if got != exp {
                bad.push("CRITyFL");
                break;
            }
        }
        for v in bad {
            r.add(
                L,
                "L2.ADBCVA",
                "ADBCVA",
                Some(v),
                Some(i),
                format!("{x} {p} at visit {num}: {v} does not follow from OE"),
            );
        }
    }
    let total: usize = source.values().map(Vec::len).sum();
    if seen != total {
        r.add(
            L,
            "L2.ADBCVA",
            "ADBCVA",
            Some("PARAMCD"),
            None,
            format!("{seen} ADBCVA records for {total} OE visual-acuity results in both scales"),
        );
    }
}
