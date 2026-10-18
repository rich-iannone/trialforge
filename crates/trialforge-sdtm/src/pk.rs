//! Pharmacokinetic domains: PC (concentrations) and PP (NCA parameters).
//!
//! PP is computed from the reported PC values (BLQ as 0), one profile per reference dose, with
//! [`trialforge_core::nca`] (PKNCA's defaults):
//! - a single dose: Cmax, Tmax, Clast, Tlast, AUClast, AUCall, λz (+ half-life, points,
//!   adjusted R²), AUCinf.obs, %extrapolated, CL/F, Vz/F and the dose-normalised Cmax and AUC;
//! - the first of several doses: Cmax, Tmax and AUCτ;
//! - the last of several doses: the same plus CL/F over τ, λz and the accumulation ratios.
//!
//! A profile with every sample BLQ has Cmax 0 and AUCs of 0, as in PKNCA. Ratios with a zero
//! denominator (CL/F over τ, the accumulation ratios) are left out: PKNCA gives `Inf`, which
//! PP can't hold.
//!
//! Doses in between have trough samples only. A profile after the first dose takes the
//! previous dose's τ sample (drawn just before this dose) as its pre-dose value at time 0.

use std::collections::BTreeMap;

use trialforge_cdisc::{MetadataError, build};
use trialforge_core::nca::{AucType, Profile};
use trialforge_core::{Dataset, DateTime, math};
use trialforge_model::Epoch;
use trialforge_processes::pk::fmt_hours;
use trialforge_processes::timeline::{PkSample, Timeline};

use crate::{Context, chars, constant, dy, ints, nums};

/// The specimen.
pub const SPECIMEN: &str = "PLASMA";
/// The bioanalytical laboratory.
const LAB: &str = "TRIALFORGE BIOANALYTICAL";
/// Concentration unit.
pub const CONC_UNIT: &str = "ng/mL";
/// How BLQ results are reported.
pub const BLQ: &str = "<BLQ";

pub(crate) fn emit(ctx: &Context<'_>) -> Result<Vec<Dataset>, MetadataError> {
    Ok(vec![pc(ctx)?, pp(ctx)?])
}

/// `PCTESTCD` for an analyte: its letters and digits, at most 8.
pub fn testcd(analyte: &str) -> String {
    analyte
        .chars()
        .filter(char::is_ascii_alphanumeric)
        .take(8)
        .collect::<String>()
        .to_ascii_uppercase()
}

/// A concentration as reported (3 significant digits, no trailing zeros).
fn fmt_conc(x: f64) -> String {
    let s = format!("{}", math::signif(x, 3));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

fn pc(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        s: &'a PkSample,
        visit: String,
        visitdy: Option<i64>,
        epoch: Epoch,
        rf: trialforge_core::Date,
        first_dose: DateTime,
    }
    let schedule = ctx.design.schedule.as_ref();
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let Some(first_dose) = t.exposure.first().and_then(|e| e.time) else {
            continue;
        };
        for (k, x) in t.pk.iter().enumerate() {
            let pv = schedule
                .and_then(|sch| sch.visit(x.visitnum))
                .expect("PK samples are taken at planned visits");
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                s: x,
                visit: pv.name.clone(),
                visitdy: pv.day.map(i64::from),
                epoch: t.epoch_of(x.time.date()),
                rf: t.first_dose(),
                first_dose,
            });
        }
    }
    let result = |r: &Row<'_>| Some(r.s.conc.map_or_else(|| BLQ.to_string(), fmt_conc));
    let tptref = |r: &Row<'_>| {
        let day = r.s.reference.date().days_since(r.first_dose.date()) + 1;
        format!("DAY {day} DOSE")
    };
    build(
        "PC",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "PC")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("PCSEQ", ints(&rows, |r| Some(r.seq))),
            ("PCTESTCD", chars(&rows, |r| Some(testcd(&r.s.analyte)))),
            ("PCTEST", chars(&rows, |r| Some(r.s.analyte.clone()))),
            ("PCORRES", chars(&rows, result)),
            ("PCORRESU", constant(&rows, CONC_UNIT)),
            ("PCSTRESC", chars(&rows, result)),
            ("PCSTRESN", nums(&rows, |r| r.s.conc)),
            ("PCSTRESU", constant(&rows, CONC_UNIT)),
            ("PCNAM", constant(&rows, LAB)),
            ("PCSPEC", constant(&rows, SPECIMEN)),
            ("PCLLOQ", nums(&rows, |r| Some(r.s.lloq))),
            ("VISITNUM", ints(&rows, |r| Some(r.s.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("PCDTC", chars(&rows, |r| Some(r.s.time.to_string()))),
            ("PCDY", ints(&rows, |r| dy(r.s.time.date(), r.rf))),
            ("PCTPT", chars(&rows, |r| Some(r.s.tpt.clone()))),
            ("PCTPTNUM", nums(&rows, |r| Some(r.s.tptnum))),
            (
                "PCELTM",
                chars(&rows, |r| {
                    Some(format!("PT{}H", fmt_hours(r.s.tptnum.max(0.0))))
                }),
            ),
            ("PCTPTREF", chars(&rows, |r| Some(tptref(r)))),
            (
                "PCRFTDTC",
                chars(&rows, |r| Some(r.s.reference.to_string())),
            ),
        ],
    )
}

/// One NCA result.
pub struct Param {
    /// `PPTESTCD`.
    pub testcd: &'static str,
    /// The value.
    pub value: f64,
    /// `PPSTRESU`, if any.
    pub unit: Option<&'static str>,
    /// Planned interval (start, end) in hours, for interval parameters.
    pub interval: Option<(f64, f64)>,
}

/// NCA of one subject's profiles: (reference dose time, parameters).
pub fn nca(t: &Timeline, tau: f64, single_doses: bool) -> Vec<(DateTime, Vec<Param>)> {
    let doses: Vec<(DateTime, f64)> = t
        .exposure
        .iter()
        .filter_map(|e| e.time.map(|time| (time, e.dose)))
        .collect();
    if doses.iter().all(|(_, d)| *d <= 0.0) {
        return Vec::new();
    }
    let mut by_ref: BTreeMap<DateTime, Vec<&PkSample>> = BTreeMap::new();
    for x in &t.pk {
        by_ref.entry(x.reference).or_default().push(x);
    }
    let n_doses = doses.len();
    let mut out = Vec::new();
    let mut first_tau: Option<(f64, f64)> = None; // (AUCτ, Cmax) of the first dose
    for (k, (dose_time, dose)) in doses.iter().enumerate() {
        let Some(samples) = by_ref.get(dose_time) else {
            continue;
        };
        if samples.len() < 3 {
            continue;
        }
        let mut time: Vec<f64> = Vec::new();
        let mut conc: Vec<f64> = Vec::new();
        // After the first dose, the previous dose's τ sample (taken just before this dose) is
        // this profile's pre-dose concentration at time 0.
        if k > 0
            && !single_doses
            && let Some(trough) = by_ref
                .get(&doses[k - 1].0)
                .and_then(|prev| prev.iter().find(|x| x.tptnum == tau))
        {
            time.push(0.0);
            conc.push(trough.conc.unwrap_or(0.0));
        }
        time.extend(
            samples
                .iter()
                .map(|x| x.time.hours_since(*dose_time).max(0.0)),
        );
        conc.extend(samples.iter().map(|x| x.conc.unwrap_or(0.0)));
        let p = Profile::new(&time, &conc);
        let mut params: Vec<Param> = Vec::new();
        let mut push = |testcd, value: Option<f64>, unit, interval| {
            if let Some(value) = value.filter(|v| v.is_finite()) {
                params.push(Param {
                    testcd,
                    value,
                    unit,
                    interval,
                });
            }
        };
        // Each profile is a single dose (single-dose and crossover studies), or one of
        // several repeated doses.
        let single = single_doses || n_doses == 1;
        let last = k + 1 == n_doses && !single;
        if single || last {
            let hl = p.half_life();
            if single {
                let inf = hl.and_then(|h| p.auc(0.0, f64::INFINITY, AucType::Inf(h.lambda_z)));
                let auclast = p.auc(0.0, f64::INFINITY, AucType::Last);
                push("CMAX", p.cmax(), Some(CONC_UNIT), None);
                push("TMAX", p.tmax(), Some("h"), None);
                push(
                    "CLST",
                    (p.clast() > 0.0).then(|| p.clast()),
                    Some(CONC_UNIT),
                    None,
                );
                push("TLST", p.tlast(), Some("h"), None);
                push("AUCLST", auclast, Some("h*ng/mL"), None);
                push(
                    "AUCALL",
                    p.auc(0.0, f64::INFINITY, AucType::All),
                    Some("h*ng/mL"),
                    None,
                );
                push("LAMZ", hl.map(|h| h.lambda_z), Some("/h"), None);
                push("LAMZHL", hl.map(|h| h.half_life), Some("h"), None);
                push("LAMZNPT", hl.map(|h| h.n_points as f64), None, None);
                push("R2ADJ", hl.map(|h| h.adj_r_squared), None, None);
                push("AUCIFO", inf, Some("h*ng/mL"), None);
                push(
                    "AUCPEO",
                    inf.zip(auclast).map(|(i, l)| (i - l) / i * 100.0),
                    Some("%"),
                    None,
                );
                // Dose in mg, AUC in h·ng/mL: CL (L/h) = dose · 1000 / AUC.
                push("CLFO", inf.map(|i| dose * 1000.0 / i), Some("L/h"), None);
                push(
                    "VZFO",
                    inf.zip(hl).map(|(i, h)| dose * 1000.0 / (i * h.lambda_z)),
                    Some("L"),
                    None,
                );
                push("CMAXD", p.cmax().map(|c| c / dose), Some("ng/mL/mg"), None);
                push("AUCIFOD", inf.map(|i| i / dose), Some("h*ng/mL/mg"), None);
            } else {
                let window = Profile::new(
                    &time
                        .iter()
                        .copied()
                        .filter(|x| *x <= tau)
                        .collect::<Vec<_>>(),
                    &time
                        .iter()
                        .zip(&conc)
                        .filter(|(x, _)| **x <= tau)
                        .map(|(_, c)| *c)
                        .collect::<Vec<_>>(),
                );
                let auctau = p.auc(0.0, tau, AucType::Last);
                let cmax = window.cmax();
                push("CMAX", cmax, Some(CONC_UNIT), Some((0.0, tau)));
                push("TMAX", window.tmax(), Some("h"), Some((0.0, tau)));
                push("AUCTAU", auctau, Some("h*ng/mL"), Some((0.0, tau)));
                push(
                    "CLFTAU",
                    auctau.filter(|a| *a > 0.0).map(|a| dose * 1000.0 / a),
                    Some("L/h"),
                    Some((0.0, tau)),
                );
                push("LAMZ", hl.map(|h| h.lambda_z), Some("/h"), None);
                push("LAMZHL", hl.map(|h| h.half_life), Some("h"), None);
                push("LAMZNPT", hl.map(|h| h.n_points as f64), None, None);
                push("R2ADJ", hl.map(|h| h.adj_r_squared), None, None);
                if let (Some((a1, c1)), Some(a), Some(c)) = (first_tau, auctau, cmax) {
                    let ratio = |x: f64, first: f64| (first > 0.0).then(|| x / first);
                    push("ARAUC", ratio(a, a1), None, Some((0.0, tau)));
                    push("ARCMAX", ratio(c, c1), None, Some((0.0, tau)));
                }
            }
        } else if k == 0 {
            let window_t: Vec<f64> = time.iter().copied().filter(|x| *x <= tau).collect();
            let window_c: Vec<f64> = time
                .iter()
                .zip(&conc)
                .filter(|(x, _)| **x <= tau)
                .map(|(_, c)| *c)
                .collect();
            let window = Profile::new(&window_t, &window_c);
            let auctau = p.auc(0.0, tau, AucType::Last);
            push("CMAX", window.cmax(), Some(CONC_UNIT), Some((0.0, tau)));
            push("TMAX", window.tmax(), Some("h"), Some((0.0, tau)));
            push("AUCTAU", auctau, Some("h*ng/mL"), Some((0.0, tau)));
            if let (Some(a), Some(c)) = (auctau, window.cmax()) {
                first_tau = Some((a, c));
            }
        }
        out.push((*dose_time, params));
    }
    out
}

/// A PP result as reported (4 significant digits).
fn fmt_param(x: f64) -> String {
    let s = format!("{}", math::signif(x, 4));
    s.strip_suffix(".0").map_or(s.clone(), str::to_string)
}

fn pp(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        analyte: &'a str,
        testcd: &'static str,
        value: f64,
        unit: Option<&'static str>,
        interval: Option<(f64, f64)>,
        reference: DateTime,
    }
    let ct = trialforge_cdisc::ct();
    let name = |cd: &str| {
        let pos = ct
            .codelist("PKPARMCD")
            .and_then(|c| c.terms.iter().position(|t| t == cd))
            .expect("PKPARMCD in CT");
        ct.codelist("PKPARM").expect("PKPARM in CT").terms[pos].clone()
    };
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let Some(analyte) = t.pk.first().map(|x| x.analyte.as_str()) else {
            continue;
        };
        let arm = &ctx.design.arms[s.arm];
        let tau = arm.treatment.cycle_days as f64 * 24.0;
        let once = |t: &trialforge_model::Treatment| t.frequency == "ONCE";
        let single = once(&arm.treatment) && arm.periods.iter().all(|p| once(&p.treatment));
        let mut k = 0;
        for (reference, params) in nca(t, tau, single) {
            for Param {
                testcd,
                value,
                unit,
                interval,
            } in params
            {
                k += 1;
                rows.push(Row {
                    usubjid: &s.usubjid,
                    seq: k,
                    analyte,
                    testcd,
                    value,
                    unit,
                    interval,
                    reference,
                });
            }
        }
    }
    let dur = |h: f64| format!("PT{}H", fmt_hours(h));
    build(
        "PP",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "PP")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("PPSEQ", ints(&rows, |r| Some(r.seq))),
            ("PPTESTCD", chars(&rows, |r| Some(r.testcd.to_string()))),
            ("PPTEST", chars(&rows, |r| Some(name(r.testcd)))),
            ("PPCAT", chars(&rows, |r| Some(r.analyte.to_string()))),
            ("PPORRES", chars(&rows, |r| Some(fmt_param(r.value)))),
            ("PPORRESU", chars(&rows, |r| r.unit.map(str::to_string))),
            ("PPSTRESC", chars(&rows, |r| Some(fmt_param(r.value)))),
            ("PPSTRESN", nums(&rows, |r| Some(r.value))),
            ("PPSTRESU", chars(&rows, |r| r.unit.map(str::to_string))),
            ("PPSPEC", constant(&rows, SPECIMEN)),
            ("PPRFTDTC", chars(&rows, |r| Some(r.reference.to_string()))),
            ("PPSTINT", chars(&rows, |r| r.interval.map(|i| dur(i.0)))),
            ("PPENINT", chars(&rows, |r| r.interval.map(|i| dur(i.1)))),
        ],
    )
}
