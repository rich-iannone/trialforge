//! L2 `L2.VACCINE`: the vaccine analysis datasets recomputed from SDTM.
//!
//! - ADSL `VAXkkDT` are the subject's vaccination dates (EX, in visit order).
//! - ADIS: every titre's original and LOG10 values (half the LLOQ below it, the ULOQ above),
//!   the baseline, `R2BASE` and the LLOQ criterion.
//! - ADCE: `ASEV` is `CESEV`, `APERIOD` is the vaccination period of the start date, and
//!   `AOCC01FL` marks one most severe reaction per subject and period.
//! - ADFACE: the maximum severity of each reaction after each vaccination (graded reactions
//!   from FACE; redness and swelling from their largest diameter), the fever occurrence of each
//!   diary temperature (VS), and `EVENTFL` agreeing within each reaction and vaccination.

use std::collections::{BTreeMap, HashMap};

use trialforge_core::{Dataset, Date, math};

use crate::table::get;
use crate::{Level, Report};

const L: Level = Level::L2;
const TOL: f64 = 1e-9;
const RULE: &str = "L2.VACCINE";

fn sev_rank(s: &str) -> Option<f64> {
    match s {
        "NONE" => Some(0.0),
        "MILD" => Some(1.0),
        "MODERATE" => Some(2.0),
        "SEVERE" => Some(3.0),
        _ => None,
    }
}

fn diameter_rank(cm: f64) -> f64 {
    if cm <= 2.0 {
        0.0
    } else if cm <= 5.0 {
        1.0
    } else if cm <= 10.0 {
        2.0
    } else {
        3.0
    }
}

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(ex), Some(sl)) = (get(sdtm, "EX"), get(adam, "ADSL")) else {
        return;
    };
    if !sl.has("VAX01DT") {
        return;
    }
    // Vaccination dates per subject, in visit order.
    let mut vax: HashMap<&str, Vec<(Option<i64>, Date)>> = HashMap::new();
    if let (Some(u), Some(vn), Some(st)) = (ex.c("USUBJID"), ex.i("VISITNUM"), ex.d("EXSTDTC")) {
        for i in 0..ex.n() {
            if let (Some(x), Some(d)) = (u[i].as_deref(), st[i]) {
                vax.entry(x).or_default().push((vn[i], d));
            }
        }
    }
    for v in vax.values_mut() {
        v.sort();
    }
    let su = sl.c("USUBJID").unwrap_or_default();
    for (k, var) in ["VAX01DT", "VAX02DT"].iter().enumerate() {
        let Some(d) = sl.d(var) else { continue };
        for i in 0..sl.n() {
            let want = su[i]
                .as_deref()
                .and_then(|x| vax.get(x))
                .and_then(|v| v.get(k))
                .map(|v| v.1);
            if d[i] != want {
                r.add(
                    L,
                    RULE,
                    "ADSL",
                    Some(var),
                    Some(i),
                    format!("{var} should be {want:?}"),
                );
            }
        }
    }
    adis(sdtm, adam, r);
    adce(sdtm, adam, r);
    adface(sdtm, adam, r);
}

fn adis(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(is), Some(ad)) = (get(sdtm, "IS"), get(adam, "ADIS")) else {
        return;
    };
    let (Some(u), Some(seq), Some(orres), Some(stresn), Some(lo), Some(hi), Some(bl)) = (
        is.c("USUBJID"),
        is.i("ISSEQ"),
        is.c("ISORRES"),
        is.f("ISSTRESN"),
        is.f("ISLLOQ"),
        is.f("ISULOQ"),
        is.c("ISBLFL"),
    ) else {
        return;
    };
    // (subject, ISSEQ) → (natural value, LLOQ, baseline?).
    let mut src: HashMap<(&str, i64), (f64, f64, bool)> = HashMap::new();
    for i in 0..is.n() {
        let (Some(x), Some(q), Some(l), Some(h)) = (u[i].as_deref(), seq[i], lo[i], hi[i]) else {
            continue;
        };
        let text = orres[i].as_deref().unwrap_or_default();
        let v = match stresn[i] {
            Some(v) => v.clamp(l / 2.0, h),
            None if text.starts_with('<') => l / 2.0,
            None if text.starts_with('>') => h,
            None => continue,
        };
        src.insert((x, q), (v, l, bl[i].as_deref() == Some("Y")));
    }
    let (Some(au), Some(pc), Some(pn), Some(aval), Some(sseq), Some(base), Some(ratio), Some(crit)) = (
        ad.c("USUBJID"),
        ad.c("PARAMCD"),
        ad.f("PARAMN"),
        ad.f("AVAL"),
        ad.i("SRCSEQ"),
        ad.f("BASE"),
        ad.f("R2BASE"),
        ad.c("CRIT1FL"),
    ) else {
        return;
    };
    let mut baselines: HashMap<(&str, &str), f64> = HashMap::new();
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(q)) = (au[i].as_deref(), pc[i].as_deref(), sseq[i]) else {
            continue;
        };
        if src.get(&(x, q)).is_some_and(|s| s.2) {
            if let Some(a) = aval[i] {
                baselines.insert((x, p), a);
            }
        }
    }
    let mut seen = 0;
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(q)) = (au[i].as_deref(), pc[i].as_deref(), sseq[i]) else {
            continue;
        };
        let Some(&(v, l, is_base)) = src.get(&(x, q)) else {
            r.add(
                L,
                RULE,
                "ADIS",
                Some("SRCSEQ"),
                Some(i),
                format!("{x} ISSEQ {q} is not in IS"),
            );
            continue;
        };
        // PARAMN 1 original, 11 LOG10, 21 4FOLD, 31 LOG10 4FOLD; 4FOLD uses the LLOQ itself.
        let kind = pn[i].unwrap_or(1.0) as i64 / 10;
        let (log, fold) = (kind % 2 == 1, kind >= 2);
        let natural = if fold && v < l { l } else { v };
        let want = math::round_dp(if log { math::log10(natural) } else { natural }, 2);
        seen += 1;
        let close = |a: Option<f64>, b: Option<f64>| match (a, b) {
            (Some(a), Some(b)) => (a - b).abs() < TOL,
            (None, None) => true,
            _ => false,
        };
        if !close(aval[i], Some(want)) {
            r.add(
                L,
                RULE,
                "ADIS",
                Some("AVAL"),
                Some(i),
                format!("{x} {p}: AVAL should be {want}"),
            );
        }
        let b = baselines.get(&(x, p)).copied();
        if !close(base[i], b) {
            r.add(
                L,
                RULE,
                "ADIS",
                Some("BASE"),
                Some(i),
                format!("{x} {p}: BASE should be {b:?}"),
            );
        }
        if !is_base && ratio[i].is_some() && !close(ratio[i], b.map(|b| want / b)) {
            r.add(
                L,
                RULE,
                "ADIS",
                Some("R2BASE"),
                Some(i),
                format!("{x} {p}: R2BASE is not AVAL / BASE"),
            );
        }
        let thr = if log { math::log10(l) } else { l };
        let want_fl = if want >= thr { "Y" } else { "N" };
        if crit[i].as_deref() != Some(want_fl) {
            r.add(
                L,
                RULE,
                "ADIS",
                Some("CRIT1FL"),
                Some(i),
                format!("{x} {p}: CRIT1FL should be {want_fl}"),
            );
        }
    }
    if seen != 4 * src.len() {
        r.add(
            L,
            RULE,
            "ADIS",
            Some("PARAMCD"),
            None,
            format!("{seen} ADIS records for {} IS results", src.len()),
        );
    }
}

fn adce(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(sl), Some(ad)) = (get(adam, "ADSL"), get(adam, "ADCE")) else {
        return;
    };
    let su = sl.c("USUBJID").unwrap_or_default();
    let mut periods: HashMap<&str, Vec<(i64, Date, Date)>> = HashMap::new();
    for k in 1..=2 {
        let (Some(s), Some(e)) = (sl.d(&format!("AP{k:02}SDT")), sl.d(&format!("AP{k:02}EDT")))
        else {
            continue;
        };
        for i in 0..sl.n() {
            if let (Some(x), Some(a), Some(b)) = (su[i].as_deref(), s[i], e[i]) {
                periods.entry(x).or_default().push((k, a, b));
            }
        }
    }
    let n_ce = get(sdtm, "CE").map_or(0, |ce| {
        ce.c("CECAT").map_or(0, |c| {
            c.iter()
                .filter(|x| x.as_deref() == Some("REACTOGENICITY"))
                .count()
        })
    });
    if ad.n() != n_ce {
        r.add(
            L,
            RULE,
            "ADCE",
            None,
            None,
            format!(
                "{} ADCE records for {n_ce} reactogenicity CE records",
                ad.n()
            ),
        );
    }
    let (Some(u), Some(st), Some(per), Some(asevn), Some(occ)) = (
        ad.c("USUBJID"),
        ad.d("ASTDT"),
        ad.i("APERIOD"),
        ad.i("ASEVN"),
        ad.c("AOCC01FL"),
    ) else {
        return;
    };
    let mut most: BTreeMap<(&str, i64), (i64, usize)> = BTreeMap::new();
    for i in 0..ad.n() {
        let Some(x) = u[i].as_deref() else { continue };
        let want = st[i].and_then(|d| {
            periods
                .get(x)
                .and_then(|p| p.iter().find(|p| p.1 <= d && d <= p.2))
                .map(|p| p.0)
        });
        if per[i] != want {
            r.add(
                L,
                RULE,
                "ADCE",
                Some("APERIOD"),
                Some(i),
                format!("{x}: APERIOD should be {want:?}"),
            );
        }
        if let (Some(p), Some(s)) = (per[i], asevn[i]) {
            let e = most.entry((x, p)).or_insert((s, 0));
            e.0 = e.0.max(s);
            if occ[i].as_deref() == Some("Y") {
                e.1 += 1;
            }
        }
    }
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(s)) = (u[i].as_deref(), per[i], asevn[i]) else {
            continue;
        };
        if occ[i].as_deref() == Some("Y") && most.get(&(x, p)).is_some_and(|m| s != m.0) {
            r.add(
                L,
                RULE,
                "ADCE",
                Some("AOCC01FL"),
                Some(i),
                format!("{x}: AOCC01FL is not on the most severe reaction"),
            );
        }
    }
    for ((x, p), (_, flags)) in most {
        if flags != 1 {
            r.add(
                L,
                RULE,
                "ADCE",
                Some("AOCC01FL"),
                None,
                format!("{x} period {p}: {flags} AOCC01FL records"),
            );
        }
    }
}

fn adface(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let (Some(face), Some(ad)) = (get(sdtm, "FACE"), get(adam, "ADFACE")) else {
        return;
    };
    // Maximum severity per (subject, object, vaccination) from FACE.
    let mut max_sev: HashMap<(&str, &str, &str), f64> = HashMap::new();
    if let (Some(u), Some(t), Some(o), Some(res), Some(num), Some(rf)) = (
        face.c("USUBJID"),
        face.c("FATESTCD"),
        face.c("FAOBJ"),
        face.c("FASTRESC"),
        face.f("FASTRESN"),
        face.c("FATPTREF"),
    ) {
        for i in 0..face.n() {
            let (Some(x), Some(t), Some(o), Some(rf)) = (
                u[i].as_deref(),
                t[i].as_deref(),
                o[i].as_deref(),
                rf[i].as_deref(),
            ) else {
                continue;
            };
            let v = match t {
                "SEV" => res[i].as_deref().and_then(sev_rank),
                "DIAMETER" => num[i].map(diameter_rank),
                _ => None,
            };
            if let Some(v) = v {
                let e = max_sev.entry((x, o, rf)).or_insert(v);
                *e = e.max(v);
            }
        }
    }
    // Fever per diary temperature: (subject, vaccination, day) → occurred.
    let mut fever: HashMap<(&str, &str, i64), bool> = HashMap::new();
    if let Some(vs) = get(sdtm, "VS") {
        if let (Some(u), Some(c), Some(v), Some(rf), Some(tp)) = (
            vs.c("USUBJID"),
            vs.c("VSCAT"),
            vs.f("VSSTRESN"),
            vs.c("VSTPTREF"),
            vs.i("VSTPTNUM"),
        ) {
            for i in 0..vs.n() {
                if let (Some(x), Some("REACTOGENICITY"), Some(t), Some(rf), Some(tp)) = (
                    u[i].as_deref(),
                    c[i].as_deref(),
                    v[i],
                    rf[i].as_deref(),
                    tp[i],
                ) {
                    fever.insert((x, rf, tp), t >= 38.0);
                }
            }
        }
    }
    let (Some(u), Some(pc), Some(param), Some(aval), Some(avalc), Some(rf), Some(tp), Some(ev)) = (
        ad.c("USUBJID"),
        ad.c("PARAMCD"),
        ad.c("PARAM"),
        ad.f("AVAL"),
        ad.c("AVALC"),
        ad.c("ATPTREF"),
        ad.i("ATPTN"),
        ad.c("EVENTFL"),
    ) else {
        return;
    };
    let mut maxima = 0;
    let mut event: HashMap<(&str, String, &str), (bool, bool)> = HashMap::new();
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(name), Some(rf)) = (
            u[i].as_deref(),
            pc[i].as_deref(),
            param[i].as_deref(),
            rf[i].as_deref(),
        ) else {
            continue;
        };
        // The object is the start of PARAM, before the test name.
        let obj = name.to_ascii_uppercase();
        let obj = [
            "OCCURRENCE INDICATOR",
            "SEVERITY/INTENSITY",
            "DIAMETER",
            "MAXIMUM SEVERITY",
            "MAXIMUM DIAMETER",
            "MAXIMUM TEMPERATURE",
        ]
        .iter()
        .find_map(|t| obj.find(&format!(" {t}")).map(|k| obj[..k].to_string()))
        .unwrap_or_default();
        let e = event.entry((x, obj.clone(), rf)).or_insert((false, false));
        if ev[i].as_deref() == Some("Y") {
            e.0 = true;
        } else {
            e.1 = true;
        }
        if name.to_ascii_uppercase().contains(" MAXIMUM SEVERITY") {
            maxima += 1;
            let want = max_sev.get(&(x, obj.as_str(), rf)).copied();
            if aval[i] != want {
                r.add(
                    L,
                    RULE,
                    "ADFACE",
                    Some("AVAL"),
                    Some(i),
                    format!("{x} {p} {rf}: maximum severity should be {want:?}"),
                );
            }
        }
        if p == "OCFEVER" {
            let want = tp[i]
                .and_then(|t| fever.get(&(x, rf, t)))
                .map(|f| if *f { "Y" } else { "N" });
            if avalc[i].as_deref() != want {
                r.add(
                    L,
                    RULE,
                    "ADFACE",
                    Some("AVALC"),
                    Some(i),
                    format!("{x} fever {rf}: AVALC should be {want:?}"),
                );
            }
        }
    }
    for ((x, obj, rf), (yes, other)) in event {
        if yes && other {
            r.add(
                L,
                RULE,
                "ADFACE",
                Some("EVENTFL"),
                None,
                format!("{x} {obj} {rf}: EVENTFL is not the same for the whole reaction"),
            );
        }
    }
    if maxima != max_sev.len() {
        r.add(
            L,
            RULE,
            "ADFACE",
            Some("PARAMCD"),
            None,
            format!(
                "{maxima} maximum severity records for {} reactions",
                max_sev.len()
            ),
        );
    }
}
