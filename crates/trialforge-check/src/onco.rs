//! Oncology checks.
//!
//! - L1 `L1.RECIST`: every RS response and TR sum of diameters is recomputed from the TR
//!   measurements and lesion states under RECIST 1.1.
//! - L2 `L2.ADRS`, `L2.ADTTE`, `L2.ADTR`: response parameters, time-to-event endpoints and
//!   sums of diameters are recomputed from SDTM (and ADSL) and compared.

use std::collections::{BTreeMap, HashMap, HashSet};

use trialforge_core::{Dataset, Date};

use crate::table::{T, get};
use crate::{Level, Report};

const TOL: f64 = 1e-6;

/// One subject's tumour data at one visit, from TR.
#[derive(Default)]
struct Visit {
    date: Option<Date>,
    /// Target lesions: (link ID, measurement).
    target: Vec<(String, Option<f64>)>,
    /// (row, value) of SUMDIAM.
    sum: Option<(usize, Option<f64>)>,
    non_target: Vec<String>,
    new_lesion: bool,
}

/// TR grouped by subject and visit number.
fn tr_visits<'a>(tr: &T<'a>) -> HashMap<&'a str, BTreeMap<i64, Visit>> {
    let mut out: HashMap<&str, BTreeMap<i64, Visit>> = HashMap::new();
    let (Some(u), Some(grp), Some(test), Some(lnk), Some(orres), Some(val), Some(vn), Some(dtc)) = (
        tr.c("USUBJID"),
        tr.c("TRGRPID"),
        tr.c("TRTESTCD"),
        tr.c("TRLNKID"),
        tr.c("TRORRES"),
        tr.f("TRSTRESN"),
        tr.i("VISITNUM"),
        tr.d("TRDTC"),
    ) else {
        return out;
    };
    for i in 0..tr.n() {
        let (Some(x), Some(v)) = (u[i].as_deref(), vn[i]) else {
            continue;
        };
        let visit = out.entry(x).or_default().entry(v).or_default();
        visit.date = visit.date.or(dtc[i]);
        match (grp[i].as_deref(), test[i].as_deref()) {
            (_, Some("SUMDIAM")) => visit.sum = Some((i, val[i])),
            (Some("TARGET"), Some("LDIAM" | "SAXIS")) => {
                visit
                    .target
                    .push((lnk[i].clone().unwrap_or_default(), val[i]));
            }
            (Some("NON-TARGET"), Some("TUMSTATE")) => {
                visit.non_target.push(orres[i].clone().unwrap_or_default());
            }
            (Some("NEW"), Some("TUMSTATE")) => {
                visit.new_lesion |= orres[i].as_deref() == Some("PRESENT");
            }
            _ => {}
        }
    }
    out
}

/// RECIST 1.1 target, non-target and overall response.
fn recist(
    target: &[(String, Option<f64>)],
    nodes: &HashSet<&str>,
    baseline: f64,
    nadir: f64,
    non_target: &[String],
    new_lesion: bool,
) -> (&'static str, Option<&'static str>, &'static str) {
    let measured: Vec<f64> = target.iter().filter_map(|(_, v)| *v).collect();
    let complete = measured.len() == target.len();
    let sum: f64 = measured.iter().sum();
    let increase = sum - nadir;
    let gone = target.iter().all(|(id, v)| match v {
        Some(v) if nodes.contains(id.as_str()) => *v < 10.0,
        Some(v) => *v == 0.0,
        None => false,
    });
    let tl = if increase >= 5.0 && sum >= nadir * 1.2 {
        "PD"
    } else if !complete {
        "NE"
    } else if gone {
        "CR"
    } else if sum <= baseline * 0.7 {
        "PR"
    } else {
        "SD"
    };
    let nt = if non_target.is_empty() {
        None
    } else if non_target.iter().any(|s| s == "UNEQUIVOCAL PROGRESSION") {
        Some("PD")
    } else if non_target.iter().all(|s| s == "ABSENT") {
        Some("CR")
    } else {
        Some("NON-CR/NON-PD")
    };
    let overall = match (tl, nt) {
        _ if new_lesion => "PD",
        ("PD", _) | (_, Some("PD")) => "PD",
        ("NE", _) => "NE",
        ("CR", None | Some("CR")) => "CR",
        ("CR", _) => "PR",
        (t, _) => t,
    };
    (tl, nt, overall)
}

/// L1: TR sums and RS responses agree with RECIST 1.1 applied to TR.
pub(crate) fn l1(sdtm: &[Dataset], r: &mut Report) {
    const L: Level = Level::L1;
    let (Some(tr), Some(rs)) = (get(sdtm, "TR"), get(sdtm, "RS")) else {
        return;
    };
    // Nodal lesions, from TU.
    let mut nodes: HashMap<&str, HashSet<&str>> = HashMap::new();
    if let Some(tu) = get(sdtm, "TU")
        && let (Some(u), Some(l), Some(loc)) = (tu.c("USUBJID"), tu.c("TULNKID"), tu.c("TULOC"))
    {
        for i in 0..tu.n() {
            if loc[i].as_deref() == Some("LYMPH NODE")
                && let (Some(u), Some(l)) = (u[i].as_deref(), l[i].as_deref())
            {
                nodes.entry(u).or_default().insert(l);
            }
        }
    }
    let mut responses: HashMap<(&str, i64, &str), (usize, &str)> = HashMap::new();
    if let (Some(u), Some(vn), Some(test), Some(res)) = (
        rs.c("USUBJID"),
        rs.i("VISITNUM"),
        rs.c("RSTESTCD"),
        rs.c("RSSTRESC"),
    ) {
        for i in 0..rs.n() {
            if let (Some(x), Some(v), Some(t), Some(s)) = (
                u[i].as_deref(),
                vn[i],
                test[i].as_deref(),
                res[i].as_deref(),
            ) {
                responses.insert((x, v, t), (i, s));
            }
        }
    }
    let empty = HashSet::new();
    let mut matched = 0;
    for (x, visits) in tr_visits(&tr) {
        let nodes = nodes.get(x).unwrap_or(&empty);
        let mut baseline = None;
        let mut nadir = f64::INFINITY;
        for (&vn, v) in &visits {
            let complete = v.target.iter().all(|(_, m)| m.is_some());
            let sum: f64 = v.target.iter().filter_map(|(_, m)| *m).sum();
            match v.sum {
                Some((row, value)) => {
                    if !complete || value.is_none_or(|s| (s - sum).abs() > TOL) {
                        r.add(
                            L,
                            "L1.RECIST",
                            "TR",
                            Some("TRSTRESN"),
                            Some(row),
                            "SUMDIAM is not the sum of the target lesion diameters".into(),
                        );
                    }
                }
                None if complete && !v.target.is_empty() => r.add(
                    L,
                    "L1.RECIST",
                    "TR",
                    Some("TRTESTCD"),
                    None,
                    format!("{x} visit {vn}: no SUMDIAM although every target lesion was measured"),
                ),
                None => {}
            }
            let Some(base) = baseline else {
                // The first assessment is the baseline: no responses.
                baseline = Some(sum);
                nadir = sum;
                if let Some((row, _)) = responses.get(&(x, vn, "OVRLRESP")) {
                    r.add(
                        L,
                        "L1.RECIST",
                        "RS",
                        Some("VISITNUM"),
                        Some(*row),
                        "response at baseline".into(),
                    );
                }
                continue;
            };
            let (tl, nt, overall) =
                recist(&v.target, nodes, base, nadir, &v.non_target, v.new_lesion);
            if complete {
                nadir = nadir.min(sum);
            }
            for (test, expected) in [
                ("TRGRESP", Some(tl)),
                ("NTRGRESP", nt),
                ("OVRLRESP", Some(overall)),
            ] {
                let got = responses.get(&(x, vn, test));
                if got.is_some() {
                    matched += 1;
                }
                if got.map(|g| g.1) != expected {
                    r.add(
                        L,
                        "L1.RECIST",
                        "RS",
                        Some("RSSTRESC"),
                        got.map(|g| g.0),
                        format!(
                            "{x} visit {vn} {test}: RECIST 1.1 gives {expected:?}, RS has {:?}",
                            got.map(|g| g.1)
                        ),
                    );
                }
            }
        }
    }
    if matched != responses.len() {
        r.add(
            L,
            "L1.RECIST",
            "RS",
            None,
            None,
            format!(
                "{} RS records have no TR assessment",
                responses.len() - matched
            ),
        );
    }
}

/// A subject's investigator overall responses on or after randomisation, in date order.
fn overall_responses<'a>(
    rs: &T<'a>,
    randdt: &HashMap<&str, Option<Date>>,
) -> HashMap<&'a str, Vec<(Date, &'a str)>> {
    let mut out: HashMap<&str, Vec<(Date, &str)>> = HashMap::new();
    let (Some(u), Some(test), Some(eval), Some(res), Some(dtc)) = (
        rs.c("USUBJID"),
        rs.c("RSTESTCD"),
        rs.c("RSEVAL"),
        rs.c("RSSTRESC"),
        rs.d("RSDTC"),
    ) else {
        return out;
    };
    for i in 0..rs.n() {
        if test[i].as_deref() != Some("OVRLRESP") || eval[i].as_deref() != Some("INVESTIGATOR") {
            continue;
        }
        let (Some(x), Some(v), Some(d)) = (u[i].as_deref(), res[i].as_deref(), dtc[i]) else {
            continue;
        };
        if randdt.get(x).copied().flatten().is_some_and(|rd| d >= rd) {
            out.entry(x).or_default().push((d, v));
        }
    }
    for v in out.values_mut() {
        v.sort();
    }
    out
}

/// Responses up to and including the first PD (admiralonco's `ANL02FL`).
fn up_to_pd<'a>(v: &[(Date, &'a str)]) -> Vec<(Date, &'a str)> {
    let end = v
        .iter()
        .position(|(_, r)| *r == "PD")
        .map_or(v.len(), |p| p + 1);
    v[..end].to_vec()
}

/// The first confirmed CR or PR (date), per admiralonco: confirmed by a CR (for a CR) or a
/// CR/PR (for a PR) at least 28 days later; in between only CR/PR/NE (CR/NE for a CR), at
/// most one NE, and for a PR no PR after a CR.
fn first_confirmed(v: &[(Date, &str)], only: Option<&str>) -> Option<Date> {
    'outer: for (i, &(d0, r0)) in v.iter().enumerate() {
        if !matches!(r0, "CR" | "PR") || only.is_some_and(|o| o != r0) {
            continue;
        }
        let mut ne = 0;
        let mut first_cr: Option<Date> = None;
        let mut last_pr: Option<Date> = None;
        for &(d, r) in &v[i + 1..] {
            let ok = if r0 == "CR" {
                matches!(r, "CR" | "NE")
            } else {
                matches!(r, "CR" | "PR" | "NE")
            };
            ne += usize::from(r == "NE");
            if r == "CR" && first_cr.is_none() {
                first_cr = Some(d);
            }
            if r == "PR" {
                last_pr = Some(d);
            }
            let confirms = d.days_since(d0) >= 28
                && if r0 == "CR" {
                    r == "CR"
                } else {
                    matches!(r, "CR" | "PR")
                };
            if confirms {
                let pr_after_cr = matches!((first_cr, last_pr), (Some(c), Some(p)) if c <= p);
                if ok && ne <= 1 && !(r0 == "PR" && pr_after_cr) {
                    return Some(d0);
                }
                continue 'outer;
            }
            if !ok {
                // A disallowed record before any confirmation: this start cannot confirm, but a
                // later confirmation window would include it too.
                continue 'outer;
            }
        }
    }
    None
}

/// (AVALC, ADT) of best overall response from the response list.
fn best_response(v: &[(Date, &str)], randdt: Date, confirmed: bool) -> (String, Option<Date>) {
    let late = |d: Date| d.days_since(randdt) >= 42;
    let first =
        |pred: &dyn Fn(Date, &str) -> bool| v.iter().find(|(d, r)| pred(*d, r)).map(|(d, _)| *d);
    let candidates: [(&str, Option<Date>); 2] = if confirmed {
        [
            ("CR", first_confirmed(v, Some("CR"))),
            ("PR", first_confirmed(v, Some("PR"))),
        ]
    } else {
        [
            ("CR", first(&|_, r| r == "CR")),
            ("PR", first(&|_, r| r == "PR")),
        ]
    };
    if let Some((c, Some(d))) = candidates.iter().find(|(_, d)| d.is_some()) {
        return (c.to_string(), Some(*d));
    }
    let rest: [(&str, Option<Date>); 4] = [
        (
            "SD",
            first(&|d, r| matches!(r, "CR" | "PR" | "SD") && late(d)),
        ),
        (
            "NON-CR/NON-PD",
            first(&|d, r| r == "NON-CR/NON-PD" && late(d)),
        ),
        ("PD", first(&|_, r| r == "PD")),
        (
            "NE",
            first(&|_, r| matches!(r, "CR" | "PR" | "SD" | "NON-CR/NON-PD" | "NE")),
        ),
    ];
    match rest.iter().find(|(_, d)| d.is_some()) {
        Some((c, d)) => (c.to_string(), *d),
        None => ("MISSING".into(), None),
    }
}

/// Expected subject-level oncology facts, recomputed from SDTM and ADSL.
struct Expected {
    randdt: Option<Date>,
    death: Option<Date>,
    last_alive: Option<Date>,
    pd: Option<Date>,
    rsp: Option<Date>,
    crsp: Option<Date>,
    last_assessment: Option<Date>,
    bor: (String, Option<Date>),
    cbor: (String, Option<Date>),
    measurable: bool,
}

/// L2: ADRS, ADTTE and ADTR against an independent recomputation.
pub(crate) fn l2(sdtm: &[Dataset], adam: &[Dataset], r: &mut Report) {
    let Some(sl) = get(adam, "ADSL") else { return };
    let (Some(slu), Some(rand), Some(dth), Some(alive)) = (
        sl.c("USUBJID"),
        sl.d("RANDDT"),
        sl.d("DTHDT"),
        sl.d("LSTALVDT"),
    ) else {
        return;
    };
    let randdt: HashMap<&str, Option<Date>> = (0..sl.n())
        .filter_map(|i| Some((slu[i].as_deref()?, rand[i])))
        .collect();

    if let Some(rs) = get(sdtm, "RS") {
        let ovr = overall_responses(&rs, &randdt);
        let mut measurable: HashSet<&str> = HashSet::new();
        if let Some(tu) = get(sdtm, "TU")
            && let (Some(u), Some(res), Some(visit)) =
                (tu.c("USUBJID"), tu.c("TUSTRESC"), tu.c("VISIT"))
        {
            for i in 0..tu.n() {
                if res[i].as_deref() == Some("TARGET")
                    && visit[i].as_deref() == Some("SCREENING")
                    && let Some(x) = u[i].as_deref()
                {
                    measurable.insert(x);
                }
            }
        }
        let none = Vec::new();
        let expected: HashMap<&str, Expected> = (0..sl.n())
            .filter_map(|i| {
                let x = slu[i].as_deref()?;
                let all = ovr.get(x).unwrap_or(&none);
                let v = &up_to_pd(all);
                let rd = rand[i]?;
                Some((
                    x,
                    Expected {
                        randdt: rand[i],
                        death: dth[i],
                        last_alive: alive[i],
                        pd: all.iter().find(|(_, r)| *r == "PD").map(|(d, _)| *d),
                        rsp: v
                            .iter()
                            .find(|(_, r)| matches!(*r, "CR" | "PR"))
                            .map(|(d, _)| *d),
                        crsp: first_confirmed(v, None),
                        last_assessment: all.last().map(|(d, _)| *d),
                        bor: best_response(v, rd, false),
                        cbor: best_response(v, rd, true),
                        measurable: measurable.contains(x),
                    },
                ))
            })
            .collect();
        if let Some(adrs) = get(adam, "ADRS") {
            check_adrs(&adrs, &expected, r);
        }
        if let Some(adtte) = get(adam, "ADTTE") {
            check_adtte(&adtte, &expected, r);
        }
    }
    if let (Some(adtr), Some(tr)) = (get(adam, "ADTR"), get(sdtm, "TR")) {
        check_adtr(&adtr, &tr, r);
    }
}

fn check_adrs(adrs: &T<'_>, expected: &HashMap<&str, Expected>, r: &mut Report) {
    const L: Level = Level::L2;
    let (Some(u), Some(pc), Some(avalc), Some(adt)) = (
        adrs.c("USUBJID"),
        adrs.c("PARAMCD"),
        adrs.c("AVALC"),
        adrs.d("ADT"),
    ) else {
        return;
    };
    let yn = |d: Option<Date>| (if d.is_some() { "Y" } else { "N" }.to_string(), d);
    let mut seen: HashSet<(&str, &str)> = HashSet::new();
    for i in 0..adrs.n() {
        let (Some(x), Some(p)) = (u[i].as_deref(), pc[i].as_deref()) else {
            continue;
        };
        let Some(e) = expected.get(x) else { continue };
        let want: Option<(String, Option<Date>)> = match p {
            "PD" => Some(yn(e.pd)),
            "RSP" => Some(yn(e.rsp)),
            "CRSP" => Some(yn(e.crsp)),
            "BOR" => Some(e.bor.clone()),
            "CBOR" => Some(e.cbor.clone()),
            "BCP" => Some(yn(e
                .bor
                .1
                .filter(|_| matches!(e.bor.0.as_str(), "CR" | "PR")))),
            "CBCP" => Some(yn(e
                .cbor
                .1
                .filter(|_| matches!(e.cbor.0.as_str(), "CR" | "PR")))),
            "DEATH" => Some(yn(e.death)),
            "MDIS" => Some(((if e.measurable { "Y" } else { "N" }).to_string(), None)),
            _ => None,
        };
        if !matches!(p, "OVR") && !seen.insert((x, p)) {
            r.add(
                L,
                "L2.ADRS",
                "ADRS",
                Some("PARAMCD"),
                Some(i),
                format!("second {p} record for {x}"),
            );
        }
        if p == "LSTA" && adt[i] != e.last_assessment {
            r.add(
                L,
                "L2.ADRS",
                "ADRS",
                Some("ADT"),
                Some(i),
                "LSTA is not the last assessment".into(),
            );
        }
        if let Some((value, date)) = want
            && (avalc[i].as_deref() != Some(value.as_str()) || adt[i] != date)
        {
            r.add(
                L,
                "L2.ADRS",
                "ADRS",
                Some("AVALC"),
                Some(i),
                format!(
                    "{p} for {x}: expected {value} on {date:?}, found {:?} on {:?}",
                    avalc[i], adt[i]
                ),
            );
        }
    }
    for x in expected.keys() {
        for p in ["PD", "RSP", "BOR", "CRSP", "CBOR", "DEATH", "LSTA", "MDIS"] {
            if !seen.contains(&(*x, p)) {
                r.add(
                    L,
                    "L2.ADRS",
                    "ADRS",
                    Some("PARAMCD"),
                    None,
                    format!("{x} has no {p} record"),
                );
            }
        }
    }
}

fn check_adtte(adtte: &T<'_>, expected: &HashMap<&str, Expected>, r: &mut Report) {
    const L: Level = Level::L2;
    let (Some(u), Some(pc), Some(aval), Some(adt), Some(cnsr), Some(start)) = (
        adtte.c("USUBJID"),
        adtte.c("PARAMCD"),
        adtte.f("AVAL"),
        adtte.d("ADT"),
        adtte.i("CNSR"),
        adtte.d("STARTDT"),
    ) else {
        return;
    };
    // (start, end, censored) per subject and parameter.
    let endpoint = |e: &Expected, p: &str| -> Option<(Date, Date, i64)> {
        let rd = e.randdt?;
        let (start, event, censor) = match p {
            "OS" => (rd, e.death, Some(e.last_alive.map_or(rd, |a| a.max(rd)))),
            "PFS" => (
                rd,
                [e.pd, e.death].into_iter().flatten().min(),
                Some(e.last_assessment.map_or(rd, |a| a.max(rd))),
            ),
            _ => (
                e.rsp?,
                [e.pd, e.death].into_iter().flatten().min(),
                e.last_assessment,
            ),
        };
        let (end, c) = match (event, censor) {
            (Some(d), _) => (d, 0),
            (None, Some(d)) => (d, 1),
            (None, None) => return None,
        };
        Some((start, end.max(start), c))
    };
    let mut seen: HashSet<(&str, &str)> = HashSet::new();
    for i in 0..adtte.n() {
        let (Some(x), Some(p)) = (u[i].as_deref(), pc[i].as_deref()) else {
            continue;
        };
        seen.insert((x, p));
        let want = expected.get(x).and_then(|e| endpoint(e, p));
        let got = match (start[i], adt[i], cnsr[i]) {
            (Some(s), Some(d), Some(c)) => Some((s, d, c)),
            _ => None,
        };
        if want != got {
            r.add(
                L,
                "L2.ADTTE",
                "ADTTE",
                Some("ADT"),
                Some(i),
                format!("{p} for {x}: expected (STARTDT, ADT, CNSR) {want:?}, found {got:?}"),
            );
            continue;
        }
        if let Some((s, d, _)) = want
            && aval[i] != Some(f64::from(d.days_since(s) + 1))
        {
            r.add(
                L,
                "L2.ADTTE",
                "ADTTE",
                Some("AVAL"),
                Some(i),
                "AVAL is not ADT - STARTDT + 1".into(),
            );
        }
    }
    for (x, e) in expected {
        for p in ["OS", "PFS", "RSD"] {
            if endpoint(e, p).is_some() && !seen.contains(&(*x, p)) {
                r.add(
                    L,
                    "L2.ADTTE",
                    "ADTTE",
                    Some("PARAMCD"),
                    None,
                    format!("{x} has no {p} record"),
                );
            }
        }
    }
}

fn check_adtr(adtr: &T<'_>, tr: &T<'_>, r: &mut Report) {
    const L: Level = Level::L2;
    let (
        Some(u),
        Some(pc),
        Some(aval),
        Some(base),
        Some(chg),
        Some(pchg),
        Some(avisitn),
        Some(anl),
    ) = (
        adtr.c("USUBJID"),
        adtr.c("PARAMCD"),
        adtr.f("AVAL"),
        adtr.f("BASE"),
        adtr.f("CHG"),
        adtr.f("PCHG"),
        adtr.f("AVISITN"),
        adtr.c("ANL01FL"),
    )
    else {
        return;
    };
    let visits = tr_visits(tr);
    let close = |a: Option<f64>, b: Option<f64>| match (a, b) {
        (Some(a), Some(b)) => (a - b).abs() < TOL,
        (None, None) => true,
        _ => false,
    };
    for i in 0..adtr.n() {
        if pc[i].as_deref() != Some("SDIAM") {
            continue;
        }
        let (Some(x), Some(av)) = (u[i].as_deref(), avisitn[i]) else {
            continue;
        };
        let Some(vs) = visits.get(x) else {
            r.add(
                L,
                "L2.ADTR",
                "ADTR",
                Some("USUBJID"),
                Some(i),
                "no TR records".into(),
            );
            continue;
        };
        let Some((&first, bl)) = vs.iter().next() else {
            continue;
        };
        let vn = if av == 0.0 { first } else { av as i64 };
        let Some(v) = vs.get(&vn) else {
            r.add(
                L,
                "L2.ADTR",
                "ADTR",
                Some("AVISITN"),
                Some(i),
                format!("no TR visit {vn}"),
            );
            continue;
        };
        let sum = |v: &Visit| v.target.iter().filter_map(|(_, m)| *m).sum::<f64>();
        let complete = v.target.iter().all(|(_, m)| m.is_some());
        let b = sum(bl);
        let a = sum(v);
        let checks = [
            ("AVAL", close(aval[i], Some(a))),
            ("BASE", close(base[i], Some(b))),
            ("CHG", close(chg[i], Some(a - b))),
            (
                "PCHG",
                close(pchg[i], (b != 0.0).then(|| (a - b) / b.abs() * 100.0)),
            ),
            ("ANL01FL", (anl[i].as_deref() == Some("Y")) == complete),
        ];
        for (var, ok) in checks {
            if !ok {
                r.add(
                    L,
                    "L2.ADTR",
                    "ADTR",
                    Some(var),
                    Some(i),
                    format!("SDIAM {var} does not match TR"),
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn confirmation() {
        let d = Date::from_days;
        assert_eq!(
            first_confirmed(&[(d(0), "PR"), (d(56), "PR")], None),
            Some(d(0))
        );
        assert_eq!(
            first_confirmed(
                &[(d(0), "PR"), (d(20), "NE"), (d(30), "NE"), (d(56), "PR")],
                None
            ),
            None
        );
        assert_eq!(
            first_confirmed(&[(d(0), "PR"), (d(10), "CR"), (d(40), "PR")], None),
            None
        );
        assert_eq!(
            first_confirmed(&[(d(0), "SD"), (d(56), "PR"), (d(112), "CR")], None),
            Some(d(56))
        );
        assert_eq!(first_confirmed(&[(d(0), "CR"), (d(56), "PR")], None), None);
    }
}
