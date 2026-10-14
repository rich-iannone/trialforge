//! Paediatric growth checks.
//!
//! - L1 `L1.AGE`: `DM.AGE` is the completed years from `BRTHDTC` to `RFICDTC`.
//! - L2 `L2.GROWTH`: every height, weight and BMI in VS has its CDC 2000 z-score and percentile
//!   in ADVS, recomputed here from the LMS parameters at the child's age in days
//!   (`ADT − BRTHDT + 1`, also `AAGECUR`), with the extended percentiles for BMI at or above the
//!   95th; and every post-baseline height has its annualized velocity, `(height − baseline
//!   height) / (days between them / 365.25)`.

use std::collections::HashMap;

use trialforge_cdisc::growth::{Measure, lms};
use trialforge_core::{Dataset, Date, math};

use crate::table::get;
use crate::{Level, Report};

const TOL: f64 = 1e-6;

/// (VS test, measure, z-score PARAMCD, percentile PARAMCD).
const INDICATORS: [(&str, Measure, &str, &str); 3] = [
    ("WEIGHT", Measure::Weight, "WGTASDS", "WGTAPCTL"),
    ("HEIGHT", Measure::Height, "HGTSDS", "HGTPCTL"),
    ("BMI", Measure::Bmi, "BMISDS", "BMIPCTL"),
];

fn completed_years(birth: Date, on: Date) -> i64 {
    let ((by, bm, bd), (y, m, d)) = (birth.ymd(), on.ymd());
    i64::from(y - by) - i64::from((m, d) < (bm, bd))
}

/// (z-score, percentile) of `x`.
fn indicator(measure: Measure, sex: &str, age: i64, x: f64) -> Option<(f64, f64)> {
    let p = lms(measure, sex, age)?;
    let z = (math::pow(x / p.m, p.l) - 1.0) / (p.l * p.s);
    match (measure, p.p95, p.sigma) {
        (Measure::Bmi, Some(p95), Some(sigma)) if x >= p95 => {
            let pctl = 90.0 + 10.0 * math::norm_cdf((x - p95) / sigma);
            let z = if pctl < 100.0 {
                math::norm_inv(pctl / 100.0)
            } else {
                8.21
            };
            Some((z, pctl))
        }
        _ => Some((z, 100.0 * math::norm_cdf(z))),
    }
}

pub(crate) fn run(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let Some(dm) = get(sdtm, "DM") else { return };
    let (Some(du), Some(birth), Some(consent), Some(age)) = (
        dm.c("USUBJID"),
        dm.d("BRTHDTC"),
        dm.d("RFICDTC"),
        dm.i("AGE"),
    ) else {
        return;
    };
    for i in 0..dm.n() {
        if let (Some(b), Some(c), Some(a)) = (birth[i], consent[i], age[i]) {
            if completed_years(b, c) != a {
                r.add(
                    Level::L1,
                    "L1.AGE",
                    "DM",
                    Some("AGE"),
                    Some(i),
                    format!("{:?}: AGE {a} is not the age at RFICDTC", du[i]),
                );
            }
        }
    }
    let (Some(vs), Some(sl), Some(ad)) = (get(sdtm, "VS"), get(adam, "ADSL"), get(adam, "ADVS"))
    else {
        return;
    };
    let (Some(su), Some(sex), Some(brth), Some(start)) =
        (sl.c("USUBJID"), sl.c("SEX"), sl.d("BRTHDT"), sl.d("TRTSDT"))
    else {
        return;
    };
    let subject: HashMap<&str, (&str, Date, Option<Date>)> = (0..sl.n())
        .filter_map(|i| Some((su[i].as_deref()?, (sex[i].as_deref()?, brth[i]?, start[i]))))
        .collect();
    // Expected values by (subject, PARAMCD, VSSEQ), and AAGECUR by (subject, VSSEQ).
    let mut expected: HashMap<(&str, &str, i64), f64> = HashMap::new();
    let mut ages: HashMap<(&str, i64), i64> = HashMap::new();
    let (Some(vu), Some(vt), Some(vv), Some(vd), Some(vq)) = (
        vs.c("USUBJID"),
        vs.c("VSTESTCD"),
        vs.f("VSSTRESN"),
        vs.d("VSDTC"),
        vs.i("VSSEQ"),
    ) else {
        return;
    };
    let mut heights: HashMap<&str, Vec<(Date, i64, f64)>> = HashMap::new();
    for i in 0..vs.n() {
        let (Some(x), Some(t), Some(v), Some(d), Some(q)) =
            (vu[i].as_deref(), vt[i].as_deref(), vv[i], vd[i], vq[i])
        else {
            continue;
        };
        let Some(&(s, b, _)) = subject.get(x) else {
            continue;
        };
        let a = i64::from(d.days_since(b)) + 1;
        ages.insert((x, q), a);
        if t == "HEIGHT" {
            heights.entry(x).or_default().push((d, q, v));
        }
        let Some(&(_, m, zcd, pcd)) = INDICATORS.iter().find(|ind| ind.0 == t) else {
            continue;
        };
        if let Some((z, p)) = indicator(m, s, a, v) {
            expected.insert((x, zcd, q), z);
            expected.insert((x, pcd, q), p);
        }
    }
    for (x, hs) in &heights {
        let Some(&(_, _, Some(st))) = subject.get(x) else {
            continue;
        };
        let Some(&(bd, _, bh)) = hs.iter().filter(|h| h.0 <= st).max_by_key(|h| (h.0, h.1)) else {
            continue;
        };
        for &(d, q, h) in hs.iter().filter(|h| h.0 > st) {
            expected.insert(
                (x, "HTVEL", q),
                (h - bh) / (f64::from(d.days_since(bd)) / 365.25),
            );
        }
    }
    let (Some(au), Some(pc), Some(aval), Some(seq), Some(aage)) = (
        ad.c("USUBJID"),
        ad.c("PARAMCD"),
        ad.f("AVAL"),
        ad.i("SRCSEQ"),
        ad.f("AAGECUR"),
    ) else {
        if !expected.is_empty() {
            r.add(
                Level::L2,
                "L2.GROWTH",
                "ADVS",
                Some("AAGECUR"),
                None,
                "paediatric ADVS lacks AAGECUR or growth parameters".into(),
            );
        }
        return;
    };
    let mut seen = 0;
    for i in 0..ad.n() {
        let (Some(x), Some(p), Some(q)) = (au[i].as_deref(), pc[i].as_deref(), seq[i]) else {
            continue;
        };
        if ages.get(&(x, q)).map(|a| *a as f64) != aage[i] {
            r.add(
                Level::L2,
                "L2.GROWTH",
                "ADVS",
                Some("AAGECUR"),
                Some(i),
                format!("{x}: AAGECUR is not ADT - BRTHDT + 1"),
            );
        }
        let derived = INDICATORS.iter().any(|ind| ind.2 == p || ind.3 == p) || p == "HTVEL";
        if !derived {
            continue;
        }
        seen += 1;
        let ok = matches!((expected.get(&(x, p, q)), aval[i]), (Some(e), Some(a)) if (e - a).abs() < TOL);
        if !ok {
            r.add(
                Level::L2,
                "L2.GROWTH",
                "ADVS",
                Some("AVAL"),
                Some(i),
                format!(
                    "{x} {p} (VSSEQ {q}): expected {:?}, found {:?}",
                    expected.get(&(x, p, q)),
                    aval[i]
                ),
            );
        }
    }
    if seen != expected.len() {
        r.add(
            Level::L2,
            "L2.GROWTH",
            "ADVS",
            Some("PARAMCD"),
            None,
            format!(
                "{seen} growth records in ADVS for {} expected",
                expected.len()
            ),
        );
    }
}
