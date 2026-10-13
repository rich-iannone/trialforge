//! Oncology BDS datasets, following admiralonco's templates (`ad_adrs.R`, `ad_adtte.R`,
//! `ad_adtr.R`, admiralonco 1.5):
//! - ADRS from RS (investigator `OVRLRESP`), TU and ADSL: OVR, PD, RSP, CB, BOR, BCP, CRSP, CCB,
//!   CBOR, CBCP, DEATH, LSTA and MDIS.
//! - ADTTE from ADRS and ADSL: OS, PFS and RSD (duration of response).
//! - ADTR from TR, TU, RS and ADSL: one parameter per target lesion plus SDIAM.
//!
//! One adaptation: nodal target lesions are measured by short axis (`TRTESTCD = SAXIS`, RECIST
//! 1.1) and become `NLDIAMnn`, where the template reads `LPERP`.

use std::collections::HashMap;

use trialforge_cdisc::build;
use trialforge_core::{Dataset, Date, Values};

use crate::{DeriveError, Frame, dt, index, rel_day};

/// admiralonco's `aval_resp()`.
fn aval_resp(c: &str) -> Option<f64> {
    Some(match c {
        "CR" => 1.0,
        "PR" => 2.0,
        "SD" => 3.0,
        "NON-CR/NON-PD" => 4.0,
        "PD" => 5.0,
        "NE" => 6.0,
        "MISSING" => 7.0,
        _ => return None,
    })
}

/// The template's `worst_resp()`, used to pick one response per date.
fn worst_resp(c: &str) -> u8 {
    match c {
        "NE" => 1,
        "CR" => 2,
        "PR" => 3,
        "SD" => 4,
        "NON-CR/NON-PD" => 5,
        "PD" => 6,
        _ => 0,
    }
}

/// Sort key with missing values last, as `dplyr::arrange()`.
fn last<T: Ord + Copy>(v: Option<T>) -> (bool, Option<T>) {
    (v.is_none(), v)
}

/// An overall-response (OVR) record.
#[derive(Clone, Debug)]
struct Ovr {
    avalc: String,
    adt: Option<Date>,
    visit: Option<String>,
    visitnum: Option<i64>,
    rsseq: Option<i64>,
    anl01: bool,
    anl02: bool,
}

/// A best-response event: the value it sets and which OVR records satisfy it.
type Event<'a> = (&'a str, &'a dyn Fn(usize, &Ovr) -> bool);

/// Per subject, `PARAMCD` → (`AVALC`, `ADT`, `ASEQ`) of the ADRS `ANL01FL` records.
type ParamFacts<'a> = HashMap<&'a str, (Option<&'a str>, Option<Date>, Option<i64>)>;

/// One ADRS record.
#[derive(Clone, Debug)]
struct RsRec {
    paramcd: &'static str,
    avalc: Option<String>,
    aval: Option<f64>,
    adt: Option<Date>,
    avisit: Option<String>,
    avisitn: Option<i64>,
    srcseq: Option<i64>,
    from_rs: bool,
    anl01: bool,
    anl02: bool,
}

impl RsRec {
    /// A parameter record copied from an OVR record (`derive_extreme_records` /
    /// `derive_extreme_event` keep the selected record's variables).
    fn from_ovr(paramcd: &'static str, o: &Ovr, avalc: &str, aval: Option<f64>) -> Self {
        Self {
            paramcd,
            avalc: Some(avalc.to_string()),
            aval,
            adt: o.adt,
            avisit: o.visit.clone(),
            avisitn: o.visitnum,
            srcseq: o.rsseq,
            from_rs: true,
            anl01: true,
            anl02: o.anl02,
        }
    }

    /// A parameter record with no source record (the subject-level fallback event).
    fn bare(paramcd: &'static str, avalc: Option<&str>, aval: Option<f64>) -> Self {
        Self {
            paramcd,
            avalc: avalc.map(str::to_string),
            aval,
            adt: None,
            avisit: None,
            avisitn: None,
            srcseq: None,
            from_rs: false,
            anl01: true,
            anl02: false,
        }
    }
}

fn yn(b: bool) -> (&'static str, Option<f64>) {
    if b {
        ("Y", Some(1.0))
    } else {
        ("N", Some(0.0))
    }
}

/// `PARAM`, `PARCAT1`, `PARCAT2`, `PARCAT3` of an ADRS parameter.
fn adrs_param(
    cd: &str,
) -> (
    &'static str,
    Option<&'static str>,
    Option<&'static str>,
    Option<&'static str>,
) {
    let resp = (
        Some("Tumor Response"),
        Some("Investigator"),
        Some("RECIST 1.1"),
    );
    let (param, cats) = match cd {
        "OVR" => ("Overall Response by Investigator", resp),
        "PD" => ("Disease Progression by Investigator", resp),
        "RSP" => ("Response by Investigator (confirmation not required)", resp),
        "CB" => (
            "Clinical Benefit by Investigator (confirmation for response not required)",
            resp,
        ),
        "BOR" => (
            "Best Overall Response by Investigator (confirmation not required)",
            resp,
        ),
        "BCP" => (
            "Best Overall Response of CR/PR by Investigator (confirmation not required)",
            resp,
        ),
        "CRSP" => ("Confirmed Response by Investigator", resp),
        "CCB" => ("Confirmed Clinical Benefit by Investigator", resp),
        "CBOR" => ("Best Confirmed Overall Response by Investigator", resp),
        "CBCP" => (
            "Best Confirmed Overall Response of CR/PR by Investigator",
            resp,
        ),
        "DEATH" => ("Death", (Some("Reference Event"), None, None)),
        "LSTA" => ("Last Disease Assessment by Investigator", resp),
        _ => (
            "Measurable Disease at Baseline by Investigator",
            (None, Some("Investigator"), Some("RECIST 1.1")),
        ),
    };
    (param, cats.0, cats.1, cats.2)
}

/// Whether `ovr[i]` is a confirmed response (admiralonco's `crsp_y_cr` / `crsp_y_pr`, equally
/// `cbor_cr` / `cbor_pr`): a later record confirms it at least 28 days on (CR for a CR; CR or PR
/// for a PR), and the records up to that one are all CR/NE (CR) or CR/PR/NE (PR), with at most
/// one NE and, for a PR, no PR after a CR. `ovr` is in date order.
fn confirmed(ovr: &[&Ovr], i: usize) -> bool {
    let (ri, Some(di)) = (ovr[i].avalc.as_str(), ovr[i].adt) else {
        return false;
    };
    if ri != "CR" && ri != "PR" {
        return false;
    }
    let confirms = |o: &&Ovr| {
        let ok = if ri == "CR" {
            o.avalc == "CR"
        } else {
            o.avalc == "CR" || o.avalc == "PR"
        };
        ok && o.adt.is_some_and(|d| d >= di.add_days(28))
    };
    let Some(j) = ovr[i + 1..].iter().position(confirms) else {
        return false;
    };
    let joined = &ovr[i + 1..=i + 1 + j];
    let allowed = |r: &str| {
        if ri == "CR" {
            matches!(r, "CR" | "NE")
        } else {
            matches!(r, "CR" | "PR" | "NE")
        }
    };
    let ne = joined.iter().filter(|o| o.avalc == "NE").count();
    if !joined.iter().all(|o| allowed(&o.avalc)) || ne > 1 {
        return false;
    }
    if ri == "PR" {
        let first_cr = joined
            .iter()
            .filter(|o| o.avalc == "CR")
            .filter_map(|o| o.adt)
            .min();
        let last_pr = joined
            .iter()
            .filter(|o| o.avalc == "PR")
            .filter_map(|o| o.adt)
            .max();
        if let (Some(c), Some(p)) = (first_cr, last_pr) {
            return c > p;
        }
    }
    true
}

/// One subject's ADRS records (before sequencing).
fn adrs_subject(
    mut ovr: Vec<Ovr>,
    randdt: Option<Date>,
    dthdt: Option<Date>,
    mdis: bool,
) -> Vec<RsRec> {
    ovr.sort_by_key(|o| (last(o.adt), last(o.rsseq)));
    // ANL01FL: the worst response per date, among records on or after randomisation.
    let eligible = |o: &Ovr| {
        aval_resp(&o.avalc).is_some() && matches!((o.adt, randdt), (Some(a), Some(r)) if a >= r)
    };
    let mut by_date: HashMap<Date, usize> = HashMap::new();
    for (k, o) in ovr.iter().enumerate() {
        if !eligible(o) {
            continue;
        }
        let adt = o.adt.expect("eligible records have a date");
        let key = |o: &Ovr| (worst_resp(&o.avalc), o.rsseq);
        match by_date.get(&adt) {
            Some(&b) if key(&ovr[b]) > key(o) => {}
            _ => {
                by_date.insert(adt, k);
            }
        }
    }
    for &k in by_date.values() {
        ovr[k].anl01 = true;
    }
    // ANL02FL: up to and including the first PD (all records without a PD).
    let first_pd = ovr.iter().position(|o| o.avalc == "PD");
    for (k, o) in ovr.iter_mut().enumerate() {
        o.anl02 = first_pd.is_none_or(|p| k <= p);
    }

    let mut out: Vec<RsRec> = ovr
        .iter()
        .map(|o| RsRec {
            anl01: o.anl01,
            ..RsRec::from_ovr("OVR", o, &o.avalc, aval_resp(&o.avalc))
        })
        .collect();
    let sel: Vec<&Ovr> = ovr.iter().filter(|o| o.anl01 && o.anl02).collect();
    let after42 = |o: &Ovr| matches!((o.adt, randdt), (Some(a), Some(r)) if a >= r.add_days(42));
    let first = |pred: &dyn Fn(usize, &Ovr) -> bool| {
        sel.iter()
            .enumerate()
            .find(|(k, o)| pred(*k, o))
            .map(|(_, o)| *o)
    };
    let flag = |cd: &'static str, o: Option<&Ovr>| match o {
        Some(o) => RsRec::from_ovr(cd, o, "Y", Some(1.0)),
        None => RsRec::bare(cd, Some("N"), Some(0.0)),
    };

    // PD: the first progression among ANL01FL records.
    let pd = ovr.iter().find(|o| o.anl01 && o.avalc == "PD");
    out.push(flag("PD", pd));
    let rsp_y = |_: usize, o: &Ovr| o.avalc == "CR" || o.avalc == "PR";
    let cb_y = |_: usize, o: &Ovr| {
        matches!(o.avalc.as_str(), "CR" | "PR" | "SD" | "NON-CR/NON-PD") && after42(o)
    };
    let crsp_y = |k: usize, _: &Ovr| confirmed(&sel, k);
    out.push(flag("RSP", first(&rsp_y)));
    out.push(flag("CB", first(&|k, o| rsp_y(k, o) || cb_y(k, o))));

    // Best overall response: the first event (in priority order) that any record satisfies.
    let best = |cd: &'static str, events: &[Event<'_>]| {
        events
            .iter()
            .find_map(|(value, pred)| {
                first(pred).map(|o| RsRec::from_ovr(cd, o, value, aval_resp(value)))
            })
            .unwrap_or_else(|| RsRec::bare(cd, Some("MISSING"), aval_resp("MISSING")))
    };
    let bor_cr = |_: usize, o: &Ovr| o.avalc == "CR";
    let bor_pr = |_: usize, o: &Ovr| o.avalc == "PR";
    let bor_sd = |_: usize, o: &Ovr| matches!(o.avalc.as_str(), "CR" | "PR" | "SD") && after42(o);
    let bor_non = |_: usize, o: &Ovr| o.avalc == "NON-CR/NON-PD" && after42(o);
    let bor_pd = |_: usize, o: &Ovr| o.avalc == "PD";
    let bor_ne = |_: usize, o: &Ovr| {
        matches!(
            o.avalc.as_str(),
            "CR" | "PR" | "SD" | "NON-CR/NON-PD" | "NE"
        )
    };
    let cbor_cr = |k: usize, o: &Ovr| o.avalc == "CR" && confirmed(&sel, k);
    let cbor_pr = |k: usize, o: &Ovr| o.avalc == "PR" && confirmed(&sel, k);
    let tail: [Event<'_>; 4] = [
        ("SD", &bor_sd),
        ("NON-CR/NON-PD", &bor_non),
        ("PD", &bor_pd),
        ("NE", &bor_ne),
    ];
    let bor_events: Vec<Event<'_>> = [
        ("CR", &bor_cr as &dyn Fn(usize, &Ovr) -> bool),
        ("PR", &bor_pr),
    ]
    .into_iter()
    .chain(tail)
    .collect();
    let cbor_events: Vec<Event<'_>> = [
        ("CR", &cbor_cr as &dyn Fn(usize, &Ovr) -> bool),
        ("PR", &cbor_pr),
    ]
    .into_iter()
    .chain(tail)
    .collect();

    let bor = best("BOR", &bor_events);
    let bcp = best_of_cr_pr("BCP", &bor);
    out.push(bor);
    out.push(bcp);
    out.push(flag("CRSP", first(&crsp_y)));
    out.push(flag("CCB", first(&|k, o| crsp_y(k, o) || cb_y(k, o))));
    let cbor = best("CBOR", &cbor_events);
    let cbcp = best_of_cr_pr("CBCP", &cbor);
    out.push(cbor);
    out.push(cbcp);

    let (v, a) = yn(dthdt.is_some());
    out.push(RsRec {
        adt: dthdt,
        ..RsRec::bare("DEATH", Some(v), a)
    });
    // LSTA: the last ANL01FL OVR record, as it is.
    out.push(match ovr.iter().rev().find(|o| o.anl01) {
        Some(o) => RsRec::from_ovr("LSTA", o, &o.avalc, aval_resp(&o.avalc)),
        None => RsRec::bare("LSTA", None, None),
    });
    let (v, a) = yn(mdis);
    out.push(RsRec::bare("MDIS", Some(v), a));
    out
}

/// `BCP` / `CBCP`: whether the best response is CR or PR (the best-response record if so).
fn best_of_cr_pr(cd: &'static str, best: &RsRec) -> RsRec {
    if matches!(best.avalc.as_deref(), Some("CR" | "PR")) {
        RsRec {
            paramcd: cd,
            avalc: Some("Y".into()),
            aval: Some(1.0),
            ..best.clone()
        }
    } else {
        RsRec::bare(cd, Some("N"), Some(0.0))
    }
}

/// Columns copied from ADSL onto every record.
fn adsl_columns<'a>(
    sl: &Frame<'a>,
    rows: &[usize],
    names: &[&'a str],
) -> Result<Vec<(&'a str, Values)>, DeriveError> {
    let mut out = Vec::new();
    for &name in names {
        if !sl.has(name) {
            continue;
        }
        let v = match &sl.ds.column(name).expect("checked").values {
            Values::Char(v) => Values::Char(rows.iter().map(|&r| v[r].clone()).collect()),
            Values::Num(v) => Values::Num(rows.iter().map(|&r| v[r]).collect()),
            Values::Int(v) => Values::Int(rows.iter().map(|&r| v[r]).collect()),
            Values::Date(v) => Values::Date(rows.iter().map(|&r| v[r]).collect()),
        };
        out.push((name, v));
    }
    Ok(out)
}

const ADSL_KEYS: [&str; 10] = [
    "STUDYID", "USUBJID", "SUBJID", "SITEID", "ARMCD", "ARM", "TRT01P", "TRT01A", "STRAT1R",
    "ITTFL",
];

/// Derives ADRS.
pub fn adrs(rs: &Dataset, tu: Option<&Dataset>, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let rs = Frame::new(rs);
    let sl = Frame::new(adsl);
    let sl_usubjid = sl.chars("USUBJID")?;
    let (randdt, dthdt) = (sl.dates("RANDDT")?, sl.dates("DTHDT")?);

    let mut ovr: HashMap<&str, Vec<Ovr>> = HashMap::new();
    let (u, testcd, eval, stresc, dtc) = (
        rs.chars("USUBJID")?,
        rs.chars("RSTESTCD")?,
        rs.chars("RSEVAL")?,
        rs.chars("RSSTRESC")?,
        rs.chars("RSDTC")?,
    );
    let (visit, visitnum, seq) = (rs.chars("VISIT")?, rs.ints("VISITNUM")?, rs.ints("RSSEQ")?);
    for i in 0..rs.n() {
        if testcd[i].as_deref() != Some("OVRLRESP") || eval[i].as_deref() != Some("INVESTIGATOR") {
            continue;
        }
        let Some(u) = u[i].as_deref() else { continue };
        ovr.entry(u).or_default().push(Ovr {
            avalc: stresc[i].clone().unwrap_or_default(),
            adt: dt(&dtc[i]),
            visit: visit[i].clone(),
            visitnum: visitnum[i],
            rsseq: seq[i],
            anl01: false,
            anl02: false,
        });
    }
    // MDIS: an investigator target lesion at screening.
    let mut measurable: std::collections::HashSet<&str> = std::collections::HashSet::new();
    if let Some(tu) = tu.map(Frame::new) {
        let (u, eval, stresc, visit) = (
            tu.chars("USUBJID")?,
            tu.chars("TUEVAL")?,
            tu.chars("TUSTRESC")?,
            tu.chars("VISIT")?,
        );
        for i in 0..tu.n() {
            if eval[i].as_deref() == Some("INVESTIGATOR")
                && stresc[i].as_deref() == Some("TARGET")
                && visit[i].as_deref() == Some("SCREENING")
                && let Some(u) = u[i].as_deref()
            {
                measurable.insert(u);
            }
        }
    }

    let mut rows: Vec<(usize, RsRec)> = Vec::new();
    for (r, u) in sl_usubjid.iter().enumerate() {
        let Some(u) = u.as_deref() else { continue };
        let mut recs = adrs_subject(
            ovr.remove(u).unwrap_or_default(),
            randdt[r],
            dthdt[r],
            measurable.contains(u),
        );
        recs.sort_by(|a, b| {
            (a.paramcd, last(a.adt), last(a.avisitn), last(a.srcseq)).cmp(&(
                b.paramcd,
                last(b.adt),
                last(b.avisitn),
                last(b.srcseq),
            ))
        });
        rows.extend(recs.into_iter().map(|rec| (r, rec)));
    }

    let sl_rows: Vec<usize> = rows.iter().map(|(r, _)| *r).collect();
    let mut aseq = Vec::with_capacity(rows.len());
    let mut k = 0;
    for (i, (r, _)) in rows.iter().enumerate() {
        k = if i > 0 && rows[i - 1].0 == *r {
            k + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }
    let recs: Vec<&RsRec> = rows.iter().map(|(_, rec)| rec).collect();
    let param = |f: fn(&str) -> Option<&'static str>| {
        Values::Char(
            recs.iter()
                .map(|r| f(r.paramcd).map(str::to_string))
                .collect(),
        )
    };
    let y = |b: bool| b.then(|| "Y".to_string());
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_KEYS)?;
    cols.extend(adsl_columns(&sl, &sl_rows, &["RANDDT"])?);
    cols.extend(vec![
        ("ASEQ", Values::Int(aseq)),
        (
            "PARAMCD",
            Values::Char(recs.iter().map(|r| Some(r.paramcd.to_string())).collect()),
        ),
        ("PARAM", param(|cd| Some(adrs_param(cd).0))),
        ("PARCAT1", param(|cd| adrs_param(cd).1)),
        ("PARCAT2", param(|cd| adrs_param(cd).2)),
        ("PARCAT3", param(|cd| adrs_param(cd).3)),
        (
            "AVALC",
            Values::Char(recs.iter().map(|r| r.avalc.clone()).collect()),
        ),
        ("AVAL", Values::Num(recs.iter().map(|r| r.aval).collect())),
        ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
        (
            "AVISIT",
            Values::Char(recs.iter().map(|r| r.avisit.clone()).collect()),
        ),
        (
            "AVISITN",
            Values::Num(recs.iter().map(|r| r.avisitn.map(|v| v as f64)).collect()),
        ),
        (
            "ANL01FL",
            Values::Char(recs.iter().map(|r| y(r.anl01)).collect()),
        ),
        (
            "ANL02FL",
            Values::Char(recs.iter().map(|r| y(r.anl02)).collect()),
        ),
        (
            "SRCDOM",
            Values::Char(
                recs.iter()
                    .map(|r| r.from_rs.then(|| "RS".to_string()))
                    .collect(),
            ),
        ),
        (
            "SRCVAR",
            Values::Char(
                recs.iter()
                    .map(|r| r.from_rs.then(|| "RSSTRESC".to_string()))
                    .collect(),
            ),
        ),
        (
            "SRCSEQ",
            Values::Int(recs.iter().map(|r| r.srcseq).collect()),
        ),
    ]);
    Ok(build("ADRS", cols)?)
}

/// A dated record that can end or censor a time-to-event parameter.
#[derive(Clone, Copy)]
struct Source {
    date: Date,
    event: bool,
    desc: &'static str,
    cnsdesc: Option<&'static str>,
    srcdom: &'static str,
    srcvar: &'static str,
    srcseq: Option<i64>,
}

/// admiral's `derive_param_tte()`: the earliest event (ties: the first source listed), else the
/// latest censoring date (ties: the last source listed); never before the start date.
fn tte(events: &[Option<Source>], censors: &[Option<Source>]) -> Option<Source> {
    let mut event: Option<Source> = None;
    for e in events.iter().flatten() {
        if event.is_none_or(|b| e.date < b.date) {
            event = Some(*e);
        }
    }
    let mut censor: Option<Source> = None;
    for c in censors.iter().flatten() {
        if censor.is_none_or(|b| c.date >= b.date) {
            censor = Some(*c);
        }
    }
    event.or(censor)
}

/// Derives ADTTE.
pub fn adtte(adrs: &Dataset, adsl: &Dataset) -> Result<Dataset, DeriveError> {
    let rs = Frame::new(adrs);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let (randdt, lstalvdt) = (sl.dates("RANDDT")?, sl.dates("LSTALVDT")?);

    // Per subject: (PARAMCD → (AVALC, ADT, ASEQ)) for ANL01FL records.
    let (u, paramcd, avalc, anl01, adt, aseq) = (
        rs.chars("USUBJID")?,
        rs.chars("PARAMCD")?,
        rs.chars("AVALC")?,
        rs.chars("ANL01FL")?,
        rs.dates("ADT")?,
        rs.ints("ASEQ")?,
    );
    let mut params: Vec<ParamFacts<'_>> = vec![HashMap::new(); sl.n()];
    for i in 0..rs.n() {
        if anl01[i].as_deref() != Some("Y") {
            continue;
        }
        let (Some(u), Some(cd)) = (u[i].as_deref(), paramcd[i].as_deref()) else {
            continue;
        };
        if let Some(&r) = sl_idx.get(u) {
            params[r].insert(cd, (avalc[i].as_deref(), adt[i], aseq[i]));
        }
    }

    struct Rec {
        sl: usize,
        paramcd: &'static str,
        start: Date,
        s: Source,
    }
    let mut recs = Vec::new();
    for r in 0..sl.n() {
        let p = &params[r];
        let adrs_src = |cd: &str, event: bool, desc: &'static str, yes: bool| {
            p.get(cd).and_then(|&(v, d, seq)| {
                let ok = !yes || v == Some("Y");
                d.filter(|_| ok).map(|date| Source {
                    date,
                    event,
                    desc,
                    cnsdesc: (!event).then_some(desc),
                    srcdom: "ADRS",
                    srcvar: "ADT",
                    srcseq: seq,
                })
            })
        };
        let death = adrs_src("DEATH", true, "Death", true);
        let pd = adrs_src("PD", true, "Disease Progression", true);
        let lasta = adrs_src("LSTA", false, "Last Tumor Assessment", false);
        let adsl_src =
            |d: Option<Date>, desc: &'static str, cns: &'static str, var: &'static str| {
                d.map(|date| Source {
                    date,
                    event: false,
                    desc,
                    cnsdesc: Some(cns),
                    srcdom: "ADSL",
                    srcvar: var,
                    srcseq: None,
                })
            };
        let alive = adsl_src(lstalvdt[r], "Alive", "Alive During Study", "LSTALVDT");
        let rand = adsl_src(randdt[r], "Randomization", "Randomization", "RANDDT");
        let Some(start) = randdt[r] else { continue };
        if let Some(s) = tte(&[death], &[alive, rand]) {
            recs.push(Rec {
                sl: r,
                paramcd: "OS",
                start,
                s,
            });
        }
        if let Some(s) = tte(&[pd, death], &[lasta, rand]) {
            recs.push(Rec {
                sl: r,
                paramcd: "PFS",
                start,
                s,
            });
        }
        let resp = p.get("RSP").filter(|v| v.0 == Some("Y")).and_then(|v| v.1);
        if let Some(start) = resp
            && let Some(s) = tte(&[pd, death], &[lasta])
        {
            recs.push(Rec {
                sl: r,
                paramcd: "RSD",
                start,
                s,
            });
        }
    }
    // Already in USUBJID, PARAMCD order.
    let sl_rows: Vec<usize> = recs.iter().map(|r| r.sl).collect();
    let mut aseq = Vec::with_capacity(recs.len());
    for (i, r) in recs.iter().enumerate() {
        let k = if i > 0 && recs[i - 1].sl == r.sl {
            aseq.last().copied().flatten().unwrap_or(0) + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }
    let adt: Vec<Date> = recs.iter().map(|r| r.s.date.max(r.start)).collect();
    let param = |cd: &str| match cd {
        "OS" => "Overall Survival",
        "PFS" => "Progression Free Survival",
        _ => "Duration of Response",
    };
    let mut cols = adsl_columns(&sl, &sl_rows, &ADSL_KEYS)?;
    cols.extend(adsl_columns(&sl, &sl_rows, &["AGE", "SEX", "RANDDT"])?);
    cols.extend(vec![
        ("ASEQ", Values::Int(aseq)),
        (
            "PARAMCD",
            Values::Char(recs.iter().map(|r| Some(r.paramcd.to_string())).collect()),
        ),
        (
            "PARAM",
            Values::Char(
                recs.iter()
                    .map(|r| Some(param(r.paramcd).to_string()))
                    .collect(),
            ),
        ),
        (
            "AVAL",
            Values::Num(
                recs.iter()
                    .zip(&adt)
                    .map(|(r, d)| Some(f64::from(d.days_since(r.start) + 1)))
                    .collect(),
            ),
        ),
        (
            "AVALU",
            Values::Char(vec![Some("DAYS".to_string()); recs.len()]),
        ),
        (
            "STARTDT",
            Values::Date(recs.iter().map(|r| Some(r.start)).collect()),
        ),
        ("ADT", Values::Date(adt.iter().map(|d| Some(*d)).collect())),
        (
            "CNSR",
            Values::Int(recs.iter().map(|r| Some(i64::from(!r.s.event))).collect()),
        ),
        (
            "EVNTDESC",
            Values::Char(recs.iter().map(|r| Some(r.s.desc.to_string())).collect()),
        ),
        (
            "CNSDTDSC",
            Values::Char(
                recs.iter()
                    .map(|r| r.s.cnsdesc.map(str::to_string))
                    .collect(),
            ),
        ),
        (
            "SRCDOM",
            Values::Char(recs.iter().map(|r| Some(r.s.srcdom.to_string())).collect()),
        ),
        (
            "SRCVAR",
            Values::Char(recs.iter().map(|r| Some(r.s.srcvar.to_string())).collect()),
        ),
        (
            "SRCSEQ",
            Values::Int(recs.iter().map(|r| r.s.srcseq).collect()),
        ),
    ]);
    Ok(build("ADTTE", cols)?)
}

/// One ADTR record.
#[derive(Clone, Debug, Default)]
struct TrRec {
    paramcd: String,
    param: String,
    tuloc: Option<String>,
    nodal: Option<bool>,
    lsexp: Option<String>,
    lsass: Option<String>,
    aval: Option<f64>,
    base: Option<f64>,
    chg: Option<f64>,
    pchg: Option<f64>,
    nadir: Option<f64>,
    adt: Option<Date>,
    ady: Option<i64>,
    avisit: Option<String>,
    avisitn: Option<i64>,
    ablfl: bool,
    pdfl: bool,
    anl01: bool,
    anl02: bool,
    anl03: bool,
    anl04: bool,
    srcseq: Option<i64>,
    sum: bool,
}

/// Derives ADTR.
pub fn adtr(
    tr: &Dataset,
    tu: &Dataset,
    rs: &Dataset,
    adsl: &Dataset,
) -> Result<Dataset, DeriveError> {
    let tr = Frame::new(tr);
    let tu = Frame::new(tu);
    let rs = Frame::new(rs);
    let sl = Frame::new(adsl);
    let sl_idx = index(sl.chars("USUBJID")?);
    let randdt = sl.dates("RANDDT")?;

    // TU: lesion location by (USUBJID, TULNKID).
    let mut loc: HashMap<(&str, &str), Option<String>> = HashMap::new();
    let (tuu, tulnk, tuloc, tueval) = (
        tu.chars("USUBJID")?,
        tu.chars("TULNKID")?,
        tu.chars("TULOC")?,
        tu.chars("TUEVAL")?,
    );
    for i in 0..tu.n() {
        if tueval[i].as_deref() != Some("INVESTIGATOR") {
            continue;
        }
        if let (Some(u), Some(l)) = (tuu[i].as_deref(), tulnk[i].as_deref()) {
            loc.insert((u, l), tuloc[i].clone());
        }
    }
    // RS: dates with an investigator overall response of PD.
    let mut pd_dates: std::collections::HashSet<(&str, Date)> = std::collections::HashSet::new();
    let (rsu, rstest, rseval, rsres, rsdtc) = (
        rs.chars("USUBJID")?,
        rs.chars("RSTESTCD")?,
        rs.chars("RSEVAL")?,
        rs.chars("RSSTRESC")?,
        rs.chars("RSDTC")?,
    );
    for i in 0..rs.n() {
        if rstest[i].as_deref() == Some("OVRLRESP")
            && rseval[i].as_deref() == Some("INVESTIGATOR")
            && rsres[i].as_deref() == Some("PD")
            && let (Some(u), Some(d)) = (rsu[i].as_deref(), dt(&rsdtc[i]))
        {
            pd_dates.insert((u, d));
        }
    }

    // Per-lesion records from investigator target-lesion measurements.
    let (u, eval, grp, testcd, lnk, stresn, dtc, visit, visitnum, seq) = (
        tr.chars("USUBJID")?,
        tr.chars("TREVAL")?,
        tr.chars("TRGRPID")?,
        tr.chars("TRTESTCD")?,
        tr.chars("TRLNKID")?,
        tr.nums("TRSTRESN")?,
        tr.chars("TRDTC")?,
        tr.chars("VISIT")?,
        tr.ints("VISITNUM")?,
        tr.ints("TRSEQ")?,
    );
    let mut per: Vec<Vec<TrRec>> = vec![Vec::new(); sl.n()];
    for i in 0..tr.n() {
        let nodal_test = match testcd[i].as_deref() {
            Some("LDIAM") => false,
            Some("SAXIS") => true,
            _ => continue,
        };
        if eval[i].as_deref() != Some("INVESTIGATOR") || grp[i].as_deref() != Some("TARGET") {
            continue;
        }
        let (Some(uu), Some(l)) = (u[i].as_deref(), lnk[i].as_deref()) else {
            continue;
        };
        let Some(&r) = sl_idx.get(uu) else { continue };
        let nr = l
            .get(1..)
            .unwrap_or("")
            .trim_start_matches(|c: char| !c.is_ascii_digit());
        let tuloc = loc.get(&(uu, l)).cloned().flatten();
        let adt = dt(&dtc[i]);
        let baseline = visit[i].as_deref() == Some("SCREENING");
        let (paramcd, param) = if nodal_test {
            (
                format!("NLDIAM{nr}"),
                format!("Target Lesion {nr} Analysis Short Axis"),
            )
        } else {
            (
                format!("LDIAM{nr}"),
                format!("Target Lesion {nr} Analysis Diameter"),
            )
        };
        per[r].push(TrRec {
            paramcd,
            param,
            nodal: tuloc.as_ref().map(|l| l == "LYMPH NODE"),
            tuloc,
            lsexp: Some(l.to_string()),
            lsass: stresn[i].map(|_| l.to_string()),
            aval: stresn[i],
            adt,
            ady: rel_day(adt, randdt[r]),
            avisit: if baseline {
                Some("BASELINE".into())
            } else {
                visit[i].clone()
            },
            avisitn: if baseline { Some(0) } else { visitnum[i] },
            anl01: stresn[i].is_some(),
            srcseq: seq[i],
            ..TrRec::default()
        });
    }

    let mut out: Vec<(usize, TrRec)> = Vec::new();
    for (r, lesions) in per.into_iter().enumerate() {
        let usubjid = sl.chars("USUBJID")?[r].as_deref().unwrap_or("");
        // Lesions entering the sum: LDIAM of non-nodal and NLDIAM of nodal lesions.
        let in_sum = |l: &TrRec| {
            (l.paramcd.starts_with("LDIAM") && l.nodal == Some(false))
                || (l.paramcd.starts_with("NLDIAM") && l.nodal == Some(true))
        };
        let mut expected: Vec<&str> = lesions
            .iter()
            .filter(|l| l.avisit.as_deref() == Some("BASELINE") && in_sum(l))
            .filter_map(|l| l.lsexp.as_deref())
            .collect();
        expected.sort_unstable();
        let lsexp = (!expected.is_empty()).then(|| expected.join(", "));
        // One SDIAM per analysis visit, in (AVISIT, AVISITN) order as dplyr's summarise.
        let mut visits: Vec<(Option<String>, Option<i64>)> = lesions
            .iter()
            .filter(|l| in_sum(l))
            .map(|l| (l.avisit.clone(), l.avisitn))
            .collect();
        visits.sort();
        visits.dedup();
        let mut sums: Vec<TrRec> = visits
            .into_iter()
            .map(|(avisit, avisitn)| {
                let at: Vec<&TrRec> = lesions
                    .iter()
                    .filter(|l| in_sum(l) && l.avisit == avisit && l.avisitn == avisitn)
                    .collect();
                let mut assessed: Vec<&str> = at
                    .iter()
                    .filter(|l| l.anl01 && l.avisit == avisit)
                    .filter_map(|l| l.lsexp.as_deref())
                    .collect();
                assessed.sort_unstable();
                let lsass = (!assessed.is_empty()).then(|| assessed.join(", "));
                let adt = at.iter().filter_map(|l| l.adt).min();
                TrRec {
                    paramcd: "SDIAM".into(),
                    param: "Target Lesions Sum of Diameters by Investigator".into(),
                    aval: Some(at.iter().filter_map(|l| l.aval).sum()),
                    adt,
                    ady: at.iter().filter_map(|l| l.ady).min(),
                    anl01: lsexp.is_some() && lsass == lsexp,
                    lsexp: lsexp.clone(),
                    lsass,
                    avisit,
                    avisitn,
                    pdfl: adt.is_some_and(|d| pd_dates.contains(&(usubjid, d))),
                    sum: true,
                    ..TrRec::default()
                }
            })
            .collect();
        // Baseline: the last record on or before day 1.
        let base_pos = sums
            .iter()
            .enumerate()
            .filter(|(_, s)| s.ady.is_some_and(|d| d <= 1))
            .max_by_key(|(k, s)| (s.ady, std::cmp::Reverse(*k)))
            .map(|(k, _)| k);
        let base = base_pos.and_then(|k| sums[k].aval);
        if let Some(k) = base_pos {
            sums[k].ablfl = true;
        }
        let snapshot: Vec<(Option<i64>, Option<f64>, bool)> =
            sums.iter().map(|s| (s.ady, s.aval, s.anl01)).collect();
        for s in &mut sums {
            s.base = base;
            if let (Some(a), Some(b)) = (s.aval, base) {
                s.chg = Some(a - b);
                s.pchg = (b != 0.0).then(|| (a - b) / b.abs() * 100.0);
            }
            s.nadir = snapshot
                .iter()
                .filter(|(ady, _, anl)| {
                    *anl && matches!((ady, s.ady), (Some(j), Some(i)) if j < &i)
                })
                .filter_map(|(_, v, _)| *v)
                .min_by(f64::total_cmp);
            s.anl04 = s.anl01 || s.pdfl;
        }
        // ANL02FL: the best (smallest) PCHG after randomisation; ties go to the earliest, as the
        // template's preceding derive_var_relative_flag() leaves the records in date order.
        let post = |s: &TrRec| matches!((s.adt, randdt[r]), (Some(a), Some(rd)) if a > rd);
        let best = sums
            .iter()
            .enumerate()
            .filter(|(_, s)| s.anl01 && post(s))
            .min_by(|(i, a), (j, b)| {
                let key = |s: &TrRec| (s.pchg.is_none(), s.pchg.unwrap_or(0.0));
                let (ka, kb) = (key(a), key(b));
                ka.0.cmp(&kb.0)
                    .then(ka.1.total_cmp(&kb.1))
                    .then(last(a.adt).cmp(&last(b.adt)))
                    .then(i.cmp(j))
            })
            .map(|(k, _)| k);
        if let Some(k) = best {
            sums[k].anl02 = true;
        }
        // ANL03FL: records before the first PD among ANL01FL / PDFL records.
        let mut cand: Vec<usize> = (0..sums.len()).filter(|&k| sums[k].anl04).collect();
        cand.sort_by_key(|&k| (last(sums[k].ady), k));
        let first_pd = cand.iter().position(|&k| sums[k].pdfl);
        for (pos, &k) in cand.iter().enumerate() {
            sums[k].anl03 = first_pd.is_none_or(|p| pos < p);
        }

        let mut recs: Vec<TrRec> = lesions.into_iter().chain(sums).collect();
        recs.sort_by(|a, b| {
            (&a.paramcd, last(a.avisitn), last(a.srcseq)).cmp(&(
                &b.paramcd,
                last(b.avisitn),
                last(b.srcseq),
            ))
        });
        out.extend(recs.into_iter().map(|rec| (r, rec)));
    }

    let sl_rows: Vec<usize> = out.iter().map(|(r, _)| *r).collect();
    let mut aseq: Vec<Option<i64>> = Vec::with_capacity(out.len());
    for (i, (r, _)) in out.iter().enumerate() {
        let k = if i > 0 && out[i - 1].0 == *r {
            aseq[i - 1].unwrap_or(0) + 1
        } else {
            1
        };
        aseq.push(Some(k));
    }
    let recs: Vec<&TrRec> = out.iter().map(|(_, rec)| rec).collect();
    let y = |b: bool| b.then(|| "Y".to_string());
    let sum_only = |f: fn(&TrRec) -> bool| {
        Values::Char(
            recs.iter()
                .map(|r| if r.sum { y(f(r)) } else { None })
                .collect(),
        )
    };
    let cat = |v: &str| Values::Char(vec![Some(v.to_string()); recs.len()]);
    let mut cols = adsl_columns(
        &sl,
        &sl_rows,
        &[
            "STUDYID", "USUBJID", "SUBJID", "SITEID", "TRT01P", "TRT01A", "STRAT1R", "RANDDT",
        ],
    )?;
    cols.extend(vec![
        ("ASEQ", Values::Int(aseq)),
        (
            "PARAMCD",
            Values::Char(recs.iter().map(|r| Some(r.paramcd.clone())).collect()),
        ),
        (
            "PARAM",
            Values::Char(recs.iter().map(|r| Some(r.param.clone())).collect()),
        ),
        ("PARCAT1", cat("Target Lesion(s)")),
        ("PARCAT2", cat("Investigator")),
        ("PARCAT3", cat("RECIST 1.1")),
        (
            "TULOC",
            Values::Char(recs.iter().map(|r| r.tuloc.clone()).collect()),
        ),
        (
            "TULOCGR1",
            Values::Char(
                recs.iter()
                    .map(|r| {
                        r.nodal
                            .map(|n| if n { "NODAL" } else { "NON-NODAL" }.to_string())
                    })
                    .collect(),
            ),
        ),
        (
            "LSEXP",
            Values::Char(recs.iter().map(|r| r.lsexp.clone()).collect()),
        ),
        (
            "LSASS",
            Values::Char(recs.iter().map(|r| r.lsass.clone()).collect()),
        ),
        ("AVAL", Values::Num(recs.iter().map(|r| r.aval).collect())),
        ("BASE", Values::Num(recs.iter().map(|r| r.base).collect())),
        ("CHG", Values::Num(recs.iter().map(|r| r.chg).collect())),
        ("PCHG", Values::Num(recs.iter().map(|r| r.pchg).collect())),
        ("NADIR", Values::Num(recs.iter().map(|r| r.nadir).collect())),
        (
            "CHGNAD",
            Values::Num(
                recs.iter()
                    .map(|r| r.aval.zip(r.nadir).map(|(a, n)| a - n))
                    .collect(),
            ),
        ),
        (
            "PCHGNAD",
            Values::Num(
                recs.iter()
                    .map(|r| match (r.aval, r.nadir) {
                        (Some(a), Some(n)) if n != 0.0 => Some(100.0 * (a - n) / n),
                        _ => None,
                    })
                    .collect(),
            ),
        ),
        ("ADT", Values::Date(recs.iter().map(|r| r.adt).collect())),
        ("ADY", Values::Int(recs.iter().map(|r| r.ady).collect())),
        (
            "AVISIT",
            Values::Char(recs.iter().map(|r| r.avisit.clone()).collect()),
        ),
        (
            "AVISITN",
            Values::Num(recs.iter().map(|r| r.avisitn.map(|v| v as f64)).collect()),
        ),
        ("ABLFL", sum_only(|r| r.ablfl)),
        ("PDFL", sum_only(|r| r.pdfl)),
        (
            "ANL01FL",
            Values::Char(recs.iter().map(|r| y(r.anl01)).collect()),
        ),
        ("ANL02FL", sum_only(|r| r.anl02)),
        ("ANL03FL", sum_only(|r| r.anl03)),
        ("ANL04FL", sum_only(|r| r.anl04)),
        (
            "SRCDOM",
            Values::Char(
                recs.iter()
                    .map(|r| (!r.sum).then(|| "TR".to_string()))
                    .collect(),
            ),
        ),
        (
            "SRCVAR",
            Values::Char(
                recs.iter()
                    .map(|r| (!r.sum).then(|| "TRSTRESN".to_string()))
                    .collect(),
            ),
        ),
        (
            "SRCSEQ",
            Values::Int(recs.iter().map(|r| r.srcseq).collect()),
        ),
    ]);
    Ok(build("ADTR", cols)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn o(day: i32, r: &str) -> Ovr {
        Ovr {
            avalc: r.into(),
            adt: Some(Date::from_days(day)),
            visit: None,
            visitnum: None,
            rsseq: Some(i64::from(day)),
            anl01: true,
            anl02: true,
        }
    }

    fn conf(v: &[(i32, &str)]) -> bool {
        let recs: Vec<Ovr> = v.iter().map(|&(d, r)| o(d, r)).collect();
        let refs: Vec<&Ovr> = recs.iter().collect();
        (0..refs.len()).any(|i| confirmed(&refs, i))
    }

    #[test]
    fn confirmation_matches_admiralonco() {
        assert!(conf(&[(0, "PR"), (56, "PR")]));
        assert!(!conf(&[(0, "PR"), (20, "PR")]));
        assert!(conf(&[(0, "PR"), (20, "NE"), (56, "PR")]));
        assert!(!conf(&[(0, "PR"), (20, "NE"), (30, "NE"), (56, "PR")]));
        assert!(!conf(&[(0, "PR"), (56, "SD")]));
        assert!(conf(&[(0, "CR"), (56, "CR")]));
        assert!(!conf(&[(0, "CR"), (56, "PR")]));
        assert!(conf(&[(0, "PR"), (56, "CR")]));
        assert!(!conf(&[(0, "PR"), (10, "CR"), (40, "PR")]));
    }

    #[test]
    fn best_response_needs_42_days_for_sd() {
        let rand = Date::from_days(0);
        let bor = |v: &[(i32, &str)]| {
            let recs = adrs_subject(
                v.iter()
                    .map(|&(d, r)| Ovr {
                        anl01: false,
                        anl02: false,
                        ..o(d, r)
                    })
                    .collect(),
                Some(rand),
                None,
                true,
            );
            let b = recs
                .iter()
                .find(|r| r.paramcd == "BOR")
                .expect("BOR")
                .avalc
                .clone();
            let c = recs
                .iter()
                .find(|r| r.paramcd == "CBOR")
                .expect("CBOR")
                .avalc
                .clone();
            (b.expect("value"), c.expect("value"))
        };
        assert_eq!(bor(&[(30, "SD"), (60, "PD")]), ("PD".into(), "PD".into()));
        assert_eq!(bor(&[(30, "SD")]), ("NE".into(), "NE".into()));
        assert_eq!(bor(&[(56, "PR"), (112, "SD")]), ("PR".into(), "SD".into()));
        assert_eq!(bor(&[]), ("MISSING".into(), "MISSING".into()));
    }

    #[test]
    fn tte_prefers_events_and_latest_censor() {
        let s = |d: i32, event: bool, desc: &'static str| Source {
            date: Date::from_days(d),
            event,
            desc,
            cnsdesc: None,
            srcdom: "X",
            srcvar: "X",
            srcseq: None,
        };
        let r = tte(
            &[Some(s(50, true, "PD")), Some(s(50, true, "Death"))],
            &[Some(s(90, false, "LSTA"))],
        );
        assert_eq!(r.map(|r| r.desc), Some("PD"));
        let r = tte(
            &[None],
            &[
                Some(s(10, false, "Alive")),
                Some(s(10, false, "Randomization")),
            ],
        );
        assert_eq!(r.map(|r| r.desc), Some("Randomization"));
    }
}
