//! Pharmacokinetic checks.
//!
//! - L1 `L1.PK`: PC results are consistent (BLQ iff no numeric result, numeric results at or
//!   above the LLOQ), every reference time is an EX administration, and post-dose samples
//!   follow their reference dose.
//! - L2 `L2.PP`: every PP parameter is recomputed from PC and EX by this crate's own NCA (the
//!   same documented rules as PKNCA's defaults, implemented independently of `trialforge-core`).
//! - L2 `L2.ADPC`, `L2.ADPP`: analysis values trace back to PC, EX and PP.

use std::collections::{BTreeMap, HashMap, HashSet};

use trialforge_core::{Dataset, DateTime};

use crate::table::{T, get};
use crate::{Level, Report};

fn dtm(s: &Option<String>) -> Option<DateTime> {
    s.as_deref().and_then(DateTime::parse_iso)
}

/// L1 checks on PC.
pub(crate) fn l1(sdtm: &[Dataset], r: &mut Report) {
    const L: Level = Level::L1;
    let (Some(pc), Some(ex)) = (get(sdtm, "PC"), get(sdtm, "EX")) else {
        return;
    };
    let mut admin: HashSet<(&str, DateTime)> = HashSet::new();
    if let (Some(u), Some(st)) = (ex.c("USUBJID"), ex.c("EXSTDTC")) {
        for i in 0..ex.n() {
            if let (Some(u), Some(t)) = (u[i].as_deref(), dtm(&st[i])) {
                admin.insert((u, t));
            }
        }
    }
    let (Some(u), Some(stresc), Some(stresn), Some(lloq), Some(dtc), Some(rft), Some(tptnum)) = (
        pc.c("USUBJID"),
        pc.c("PCSTRESC"),
        pc.f("PCSTRESN"),
        pc.f("PCLLOQ"),
        pc.c("PCDTC"),
        pc.c("PCRFTDTC"),
        pc.f("PCTPTNUM"),
    ) else {
        return;
    };
    for i in 0..pc.n() {
        let blq = stresc[i].as_deref() == Some("<BLQ");
        if blq == stresn[i].is_some() {
            r.add(
                L,
                "L1.PK",
                "PC",
                Some("PCSTRESN"),
                Some(i),
                "PCSTRESN must be missing exactly when the result is <BLQ".into(),
            );
        }
        if let (Some(v), Some(l)) = (stresn[i], lloq[i]) {
            if v < l {
                r.add(
                    L,
                    "L1.PK",
                    "PC",
                    Some("PCSTRESN"),
                    Some(i),
                    format!("{v} is below the LLOQ {l}"),
                );
            }
        }
        let (Some(x), Some(reference), Some(t)) = (u[i].as_deref(), dtm(&rft[i]), dtm(&dtc[i]))
        else {
            r.add(
                L,
                "L1.PK",
                "PC",
                Some("PCRFTDTC"),
                Some(i),
                "sample or reference time is missing".into(),
            );
            continue;
        };
        if !admin.contains(&(x, reference)) {
            r.add(
                L,
                "L1.PK",
                "PC",
                Some("PCRFTDTC"),
                Some(i),
                "PCRFTDTC is not an EX administration".into(),
            );
        }
        let pre = tptnum[i].is_some_and(|n| n < 0.0);
        if pre == (t > reference) {
            r.add(
                L,
                "L1.PK",
                "PC",
                Some("PCDTC"),
                Some(i),
                "pre-dose samples must precede, and others follow, their reference dose".into(),
            );
        }
    }
}

/// PC samples per (subject, reference dose): (hours after dose, concentration, PCTPTNUM).
type Samples = HashMap<(String, DateTime), Vec<(f64, f64, f64)>>;
/// A source record's fields by (USUBJID, --SEQ).
type BySeq<'a, V> = HashMap<(&'a str, i64), V>;

/// A profile point: (hours after dose, concentration with BLQ as 0).
type Pt = (f64, f64);

/// An independent NCA on one profile (PKNCA defaults).
mod nca {
    use super::Pt;

    /// Drops BLQ zeros strictly between the first and last measurable concentrations.
    pub fn clean(pts: &[Pt]) -> Vec<Pt> {
        let pos: Vec<usize> = (0..pts.len()).filter(|&i| pts[i].1 > 0.0).collect();
        let (Some(&a), Some(&b)) = (pos.first(), pos.last()) else {
            return pts.to_vec();
        };
        pts.iter()
            .enumerate()
            .filter(|(i, p)| !(p.1 == 0.0 && *i > a && *i < b))
            .map(|(_, p)| *p)
            .collect()
    }

    pub fn cmax(pts: &[Pt]) -> Option<(f64, f64)> {
        let mut best: Option<(f64, f64)> = None;
        for &(t, c) in pts {
            if c > 0.0 && best.is_none_or(|(_, b)| c > b) {
                best = Some((t, c));
            }
        }
        best
    }

    pub fn last(pts: &[Pt]) -> Option<(f64, f64)> {
        pts.iter().rev().find(|p| p.1 > 0.0).copied()
    }

    /// (λz, adjusted R², points).
    pub fn lambda_z(pts: &[Pt]) -> Option<(f64, f64, usize)> {
        let (tmax, _) = cmax(pts)?;
        let tail: Vec<(f64, f64)> = pts
            .iter()
            .filter(|p| p.0 > tmax && p.1 > 0.0)
            .map(|p| (p.0, trialforge_core::math::ln(p.1)))
            .collect();
        let mut fits = Vec::new();
        for n in 3..=tail.len() {
            let s = &tail[tail.len() - n..];
            let k = n as f64;
            let (sx, sy) = s.iter().fold((0.0, 0.0), |a, p| (a.0 + p.0, a.1 + p.1));
            let (mx, my) = (sx / k, sy / k);
            let (mut sxx, mut sxy, mut syy) = (0.0, 0.0, 0.0);
            for p in s {
                sxx += (p.0 - mx) * (p.0 - mx);
                sxy += (p.0 - mx) * (p.1 - my);
                syy += (p.1 - my) * (p.1 - my);
            }
            let b = sxy / sxx;
            let ssr: f64 = s
                .iter()
                .map(|p| {
                    let e = p.1 - (my - b * mx) - b * p.0;
                    e * e
                })
                .sum();
            let r2 = 1.0 - ssr / syy;
            fits.push((-b, 1.0 - (1.0 - r2) * (k - 1.0) / (k - 2.0), n));
        }
        let top = fits
            .iter()
            .map(|f| f.1)
            .filter(|x| !x.is_nan())
            .fold(f64::NEG_INFINITY, f64::max);
        fits.into_iter()
            .filter(|f| f.0 > 0.0 && f.1 > top - 1e-4)
            .max_by_key(|f| f.2)
    }

    fn segment(a: Pt, b: Pt) -> f64 {
        if b.1 < a.1 && b.1 > 0.0 {
            (a.1 - b.1) / trialforge_core::math::ln(a.1 / b.1) * (b.0 - a.0)
        } else {
            (a.1 + b.1) / 2.0 * (b.0 - a.0)
        }
    }

    /// AUC from 0 to the last measurable concentration, plus (for `all`) the linear tail to
    /// the next sample.
    pub fn auc(pts: &[Pt], all: bool) -> Option<f64> {
        if pts.len() < 2 || pts[0].0 > 0.0 {
            return None;
        }
        let Some((tl, _)) = last(pts) else {
            return Some(0.0);
        };
        let mut total = 0.0;
        for w in pts.windows(2) {
            if w[1].0 <= tl {
                if !(w[0].1 == 0.0 && w[1].1 == 0.0) {
                    total += segment(w[0], w[1]);
                }
            } else if all && w[0].0 == tl {
                total += (w[0].1 + w[1].1) / 2.0 * (w[1].0 - w[0].0);
            }
        }
        Some(total)
    }

    /// AUC over [0, end] with the concentration at `end` interpolated (lin up/log down) when
    /// `end` falls before the last measurable concentration, and 0 after it. With nothing
    /// measurable the AUC is 0, as in PKNCA (`aucint.last` of an all-BLQ interval).
    pub fn auc_to(pts: &[Pt], end: f64) -> Option<f64> {
        let tl = last(pts).map(|(t, _)| t);
        let mut kept: Vec<Pt> = pts.iter().copied().filter(|p| p.0 < end).collect();
        let at = if let Some(p) = pts.iter().find(|p| p.0 == end) {
            p.1
        } else if tl.is_some_and(|tl| end < tl) {
            let i = pts.iter().rposition(|p| p.0 < end)?;
            let (a, b) = (pts[i], pts[i + 1]);
            if b.1 < a.1 && b.1 > 0.0 {
                let ln = trialforge_core::math::ln;
                trialforge_core::math::exp(
                    ln(a.1) + (end - a.0) / (b.0 - a.0) * (ln(b.1) - ln(a.1)),
                )
            } else {
                a.1 + (end - a.0) / (b.0 - a.0) * (b.1 - a.1)
            }
        } else {
            0.0
        };
        kept.push((end, at));
        auc(&kept, false)
    }
}

/// Recomputes PP from PC and EX.
fn expected_pp(sdtm: &[Dataset]) -> HashMap<(String, DateTime, String), f64> {
    let mut out = HashMap::new();
    let (Some(pc), Some(ex)) = (get(sdtm, "PC"), get(sdtm, "EX")) else {
        return out;
    };
    // Administrations with a dose, per subject, in time order.
    let mut doses: BTreeMap<String, Vec<(DateTime, f64)>> = BTreeMap::new();
    // Subjects with a repeated-dose regimen (any EXDOSFRQ other than ONCE).
    let mut repeated: HashSet<String> = HashSet::new();
    if let (Some(u), Some(st), Some(dose), Some(frq)) = (
        ex.c("USUBJID"),
        ex.c("EXSTDTC"),
        ex.f("EXDOSE"),
        ex.c("EXDOSFRQ"),
    ) {
        for i in 0..ex.n() {
            if let (Some(u), Some(t), Some(d)) = (u[i].clone(), dtm(&st[i]), dose[i]) {
                if frq[i].as_deref() != Some("ONCE") {
                    repeated.insert(u.clone());
                }
                doses.entry(u).or_default().push((t, d));
            }
        }
    }
    let (Some(u), Some(stresn), Some(dtc), Some(rft), Some(tptnum)) = (
        pc.c("USUBJID"),
        pc.f("PCSTRESN"),
        pc.c("PCDTC"),
        pc.c("PCRFTDTC"),
        pc.f("PCTPTNUM"),
    ) else {
        return out;
    };
    let mut samples: Samples = HashMap::new();
    for i in 0..pc.n() {
        let (Some(x), Some(t), Some(rf)) = (u[i].clone(), dtm(&dtc[i]), dtm(&rft[i])) else {
            continue;
        };
        samples.entry((x, rf)).or_default().push((
            t.hours_since(rf).max(0.0),
            stresn[i].unwrap_or(0.0),
            tptnum[i].unwrap_or(0.0),
        ));
    }
    for (x, ds) in &mut doses {
        ds.sort_by_key(|d| d.0);
        if ds.iter().all(|d| d.1 <= 0.0) {
            continue;
        }
        let n = ds.len();
        let single = n == 1 || !repeated.contains(x);
        let tau = 24.0;
        let mut first: Option<(f64, f64)> = None;
        for (k, (rf, dose)) in ds.iter().enumerate() {
            let Some(s) = samples.get(&(x.clone(), *rf)) else {
                continue;
            };
            if s.len() < 3 {
                continue;
            }
            let mut pts: Vec<Pt> = s.iter().map(|p| (p.0, p.1)).collect();
            if k > 0 && !single {
                let trough = samples
                    .get(&(x.clone(), ds[k - 1].0))
                    .and_then(|prev| prev.iter().find(|p| p.2 == tau));
                if let Some(tr) = trough {
                    pts.insert(0, (0.0, tr.1));
                }
            }
            pts.sort_by(|a, b| a.0.total_cmp(&b.0));
            let pts = nca::clean(&pts);
            let mut put = |cd: &str, v: Option<f64>| {
                if let Some(v) = v.filter(|v| v.is_finite()) {
                    out.insert((x.clone(), *rf, cd.to_string()), v);
                }
            };
            let lz = nca::lambda_z(&pts);
            if single {
                let (tmax, cmax) = nca::cmax(&pts).unzip();
                let cmax = cmax.or((!pts.is_empty()).then_some(0.0));
                let (tlast, clast) = nca::last(&pts).unzip();
                let auclast = nca::auc(&pts, false);
                let inf = lz.zip(auclast).zip(clast).map(|((l, a), c)| a + c / l.0);
                put("CMAX", cmax);
                put("TMAX", tmax);
                put("CLST", clast);
                put("TLST", tlast);
                put("AUCLST", auclast);
                put("AUCALL", nca::auc(&pts, true));
                put("LAMZ", lz.map(|l| l.0));
                put("LAMZHL", lz.map(|l| std::f64::consts::LN_2 / l.0));
                put("LAMZNPT", lz.map(|l| l.2 as f64));
                put("R2ADJ", lz.map(|l| l.1));
                put("AUCIFO", inf);
                put("AUCPEO", inf.zip(auclast).map(|(i, a)| (i - a) / i * 100.0));
                put("CLFO", inf.map(|i| dose * 1000.0 / i));
                put("VZFO", inf.zip(lz).map(|(i, l)| dose * 1000.0 / (i * l.0)));
                put("CMAXD", cmax.map(|c| c / dose));
                put("AUCIFOD", inf.map(|i| i / dose));
            } else if k == 0 || k + 1 == n {
                let window: Vec<Pt> = pts.iter().copied().filter(|p| p.0 <= tau).collect();
                // An all-BLQ window has Cmax 0 and no Tmax (PKNCA).
                let (tmax, cmax) = nca::cmax(&window).unzip();
                let cmax = cmax.or((!window.is_empty()).then_some(0.0));
                let auctau = nca::auc_to(&pts, tau);
                put("CMAX", cmax);
                put("TMAX", tmax);
                put("AUCTAU", auctau);
                if k == 0 {
                    first = auctau.zip(cmax);
                } else {
                    put(
                        "CLFTAU",
                        auctau.filter(|a| *a > 0.0).map(|a| dose * 1000.0 / a),
                    );
                    put("LAMZ", lz.map(|l| l.0));
                    put("LAMZHL", lz.map(|l| std::f64::consts::LN_2 / l.0));
                    put("LAMZNPT", lz.map(|l| l.2 as f64));
                    put("R2ADJ", lz.map(|l| l.1));
                    if let (Some((a1, c1)), Some(a), Some(c)) = (first, auctau, cmax) {
                        put("ARAUC", (a1 > 0.0).then(|| a / a1));
                        put("ARCMAX", (c1 > 0.0).then(|| c / c1));
                    }
                }
            }
        }
    }
    out
}

/// L2 checks.
pub(crate) fn l2(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    const L: Level = Level::L2;
    if let Some(pp) = get(sdtm, "PP") {
        let want = expected_pp(sdtm);
        let mut seen = 0;
        if let (Some(u), Some(cd), Some(v), Some(rf)) = (
            pp.c("USUBJID"),
            pp.c("PPTESTCD"),
            pp.f("PPSTRESN"),
            pp.c("PPRFTDTC"),
        ) {
            for i in 0..pp.n() {
                let key = (
                    u[i].clone().unwrap_or_default(),
                    dtm(&rf[i]).unwrap_or(DateTime::from_seconds(0)),
                    cd[i].clone().unwrap_or_default(),
                );
                match (want.get(&key), v[i]) {
                    (Some(w), Some(x)) if (w - x).abs() <= 1e-9 * w.abs().max(1.0) => seen += 1,
                    (w, x) => r.add(
                        L,
                        "L2.PP",
                        "PP",
                        Some("PPSTRESN"),
                        Some(i),
                        format!("{}: NCA gives {w:?}, PP has {x:?}", key.2),
                    ),
                }
            }
        }
        if seen != want.len() && seen == pp.n() {
            r.add(
                L,
                "L2.PP",
                "PP",
                None,
                None,
                format!(
                    "{} parameters computable from PC are missing",
                    want.len() - seen
                ),
            );
        }
    }
    if let (Some(adpc), Some(pc)) = (get(adam, "ADPC"), get(sdtm, "PC")) {
        check_adpc(&adpc, &pc, r);
    }
    if let (Some(adpp), Some(pp)) = (get(adam, "ADPP"), get(sdtm, "PP")) {
        let src: BySeq<'_, (Option<&str>, Option<f64>)> = match (
            pp.c("USUBJID"),
            pp.i("PPSEQ"),
            pp.c("PPTESTCD"),
            pp.f("PPSTRESN"),
        ) {
            (Some(u), Some(q), Some(c), Some(v)) => (0..pp.n())
                .filter_map(|i| Some(((u[i].as_deref()?, q[i]?), (c[i].as_deref(), v[i]))))
                .collect(),
            _ => HashMap::new(),
        };
        if let (Some(u), Some(q), Some(c), Some(v)) = (
            adpp.c("USUBJID"),
            adpp.i("SRCSEQ"),
            adpp.c("PARAMCD"),
            adpp.f("AVAL"),
        ) {
            for i in 0..adpp.n() {
                let ok = u[i]
                    .as_deref()
                    .zip(q[i])
                    .and_then(|k| src.get(&k))
                    .is_some_and(|(cd, x)| *cd == c[i].as_deref() && *x == v[i]);
                if !ok {
                    r.add(
                        L,
                        "L2.ADPP",
                        "ADPP",
                        Some("AVAL"),
                        Some(i),
                        "PARAMCD/AVAL do not match PP".into(),
                    );
                }
            }
        }
    }
}

fn check_adpc(adpc: &T<'_>, pc: &T<'_>, r: &mut Report) {
    const L: Level = Level::L2;
    let src: BySeq<'_, (Option<f64>, Option<f64>, Option<DateTime>)> = match (
        pc.c("USUBJID"),
        pc.i("PCSEQ"),
        pc.f("PCSTRESN"),
        pc.f("PCLLOQ"),
        pc.c("PCDTC"),
    ) {
        (Some(u), Some(q), Some(v), Some(l), Some(d)) => (0..pc.n())
            .filter_map(|i| Some(((u[i].as_deref()?, q[i]?), (v[i], l[i], dtm(&d[i])))))
            .collect(),
        _ => HashMap::new(),
    };
    let (Some(u), Some(dom), Some(q), Some(aval), Some(dtype), Some(base), Some(chg), Some(ablfl)) = (
        adpc.c("USUBJID"),
        adpc.c("SRCDOM"),
        adpc.i("SRCSEQ"),
        adpc.f("AVAL"),
        adpc.c("DTYPE"),
        adpc.f("BASE"),
        adpc.f("CHG"),
        adpc.c("ABLFL"),
    ) else {
        return;
    };
    for i in 0..adpc.n() {
        if dom[i].as_deref() != Some("PC") {
            continue;
        }
        let Some((v, lloq, _)) = u[i].as_deref().zip(q[i]).and_then(|k| src.get(&k)) else {
            r.add(
                L,
                "L2.ADPC",
                "ADPC",
                Some("SRCSEQ"),
                Some(i),
                "no PC record".into(),
            );
            continue;
        };
        let want = if dtype[i].as_deref().is_some_and(|d| d.ends_with("HALFLLOQ")) {
            lloq.map(|l| l * 0.5)
        } else {
            *v
        };
        if aval[i] != want {
            r.add(
                L,
                "L2.ADPC",
                "ADPC",
                Some("AVAL"),
                Some(i),
                format!("AVAL should be {want:?}"),
            );
        }
        if ablfl[i].is_some() && base[i] != aval[i] {
            r.add(
                L,
                "L2.ADPC",
                "ADPC",
                Some("BASE"),
                Some(i),
                "baseline record's BASE is not its AVAL".into(),
            );
        }
        let ok = match (aval[i], base[i], chg[i]) {
            (Some(a), Some(b), Some(c)) => (a - b - c).abs() < 1e-9,
            (_, _, None) => aval[i].is_none() || base[i].is_none(),
            _ => false,
        };
        if !ok {
            r.add(
                L,
                "L2.ADPC",
                "ADPC",
                Some("CHG"),
                Some(i),
                "CHG is not AVAL - BASE".into(),
            );
        }
    }
}
