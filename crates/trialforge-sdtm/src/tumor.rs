//! Oncology domains: TU, TR, RS (investigator assessment), under RECIST 1.1 or, for lymphoma,
//! the Lugano 2014 classification (perpendicular diameters, SPD, the FDG-PET Deauville score,
//! and CT- and PET-CT-based overall responses).

use trialforge_cdisc::{MetadataError, build};
use trialforge_core::{Dataset, Date};
use trialforge_model::Epoch;
use trialforge_processes::timeline::{Lesion, Timeline};

use crate::{Context, chars, constant, dy, fmt_dp, ints, nums};

const METHOD: &str = "CT SCAN";
const EVALUATOR: &str = "INVESTIGATOR";

pub(crate) fn emit(ctx: &Context<'_>) -> Result<Vec<Dataset>, MetadataError> {
    Ok(vec![tu(ctx)?, tr(ctx)?, rs(ctx)?])
}

fn visit_name(ctx: &Context<'_>, visitnum: i64) -> String {
    ctx.design
        .schedule
        .as_ref()
        .and_then(|sch| sch.visit(visitnum))
        .expect("tumour assessments are taken at planned visits")
        .name
        .clone()
}

fn tu(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        l: &'a Lesion,
        visit: String,
        epoch: Epoch,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, l) in t.lesions.iter().enumerate() {
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                l,
                visit: visit_name(ctx, l.visitnum),
                epoch: t.epoch_of(l.date),
                rf: t.first_dose(),
            });
        }
    }
    build(
        "TU",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "TU")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("TUSEQ", ints(&rows, |r| Some(r.seq))),
            ("TULNKID", chars(&rows, |r| Some(r.l.id.clone()))),
            ("TUTESTCD", constant(&rows, "TUMIDENT")),
            ("TUTEST", constant(&rows, "Tumor Identification")),
            (
                "TUORRES",
                chars(&rows, |r| Some(r.l.kind.as_str().to_string())),
            ),
            (
                "TUSTRESC",
                chars(&rows, |r| Some(r.l.kind.as_str().to_string())),
            ),
            ("TULOC", chars(&rows, |r| Some(r.l.location.clone()))),
            ("TULAT", chars(&rows, |r| r.l.laterality.clone())),
            ("TUMETHOD", constant(&rows, METHOD)),
            ("TUEVAL", constant(&rows, EVALUATOR)),
            ("VISITNUM", ints(&rows, |r| Some(r.l.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("TUDTC", chars(&rows, |r| Some(r.l.date.to_string()))),
            ("TUDY", ints(&rows, |r| dy(r.l.date, r.rf))),
        ],
    )
}

/// One TR record before sequencing.
struct TrRow<'a> {
    usubjid: &'a str,
    seq: i64,
    group: &'static str,
    link: Option<String>,
    testcd: &'static str,
    /// Measurement, or a lesion state.
    result: TrResult,
    /// Unit of a measurement.
    unit: Option<&'static str>,
    method: &'static str,
    visitnum: i64,
    visit: String,
    epoch: Epoch,
    date: Date,
    rf: Date,
}

enum TrResult {
    Mm(Option<f64>),
    State(String),
}

fn tr_rows<'a>(ctx: &Context<'_>, usubjid: &'a str, t: &'a Timeline) -> Vec<TrRow<'a>> {
    let node = |id: &str| t.lesions.iter().any(|l| l.id == id && l.node);
    let mut rows = Vec::new();
    for a in &t.tumor {
        let visit = visit_name(ctx, a.visitnum);
        let lugano = !a.perpendicular.is_empty();
        let mut push = |group, link: Option<&str>, testcd, result: TrResult| {
            let unit = match (&result, testcd) {
                (TrResult::Mm(_), "PPD" | "SUMPPD") => Some("mm2"),
                (TrResult::Mm(_), "FDP5PS") => None,
                (TrResult::Mm(_), _) => Some("mm"),
                _ => None,
            };
            rows.push(TrRow {
                usubjid,
                seq: 0,
                group,
                link: link.map(str::to_string),
                testcd,
                result,
                unit,
                method: if testcd == "FDP5PS" {
                    "PET/CT SCAN"
                } else {
                    METHOD
                },
                visitnum: a.visitnum,
                visit: visit.clone(),
                epoch: t.epoch_of(a.date),
                date: a.date,
                rf: t.first_dose(),
            });
        };
        if lugano {
            for ((id, l), (_, w)) in a.target.iter().zip(&a.perpendicular) {
                push("TARGET", Some(id), "LDIAM", TrResult::Mm(*l));
                push("TARGET", Some(id), "LPERP", TrResult::Mm(*w));
                push(
                    "TARGET",
                    Some(id),
                    "PPD",
                    TrResult::Mm(l.zip(*w).map(|(l, w)| l * w)),
                );
            }
            if let Some(sum) = a.sum {
                push("TARGET", None, "SUMPPD", TrResult::Mm(Some(sum)));
            }
            if let Some(score) = a.deauville {
                push("PET", None, "FDP5PS", TrResult::Mm(Some(f64::from(score))));
            }
        } else {
            for (id, v) in &a.target {
                let testcd = if node(id) { "SAXIS" } else { "LDIAM" };
                push("TARGET", Some(id), testcd, TrResult::Mm(*v));
            }
            if let Some(sum) = a.sum {
                push("TARGET", None, "SUMDIAM", TrResult::Mm(Some(sum)));
            }
        }
        for (id, state) in &a.non_target {
            push(
                "NON-TARGET",
                Some(id),
                "TUMSTATE",
                TrResult::State(state.clone()),
            );
        }
        for id in &a.new_lesions {
            push(
                "NEW",
                Some(id),
                "TUMSTATE",
                TrResult::State("PRESENT".into()),
            );
        }
    }
    for (k, r) in rows.iter_mut().enumerate() {
        r.seq = k as i64 + 1;
    }
    rows
}

fn tr(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        rows.extend(tr_rows(ctx, &s.usubjid, t));
    }
    let test = |cd: &str| match cd {
        "LDIAM" => "Longest Diameter",
        "SAXIS" => "Short Axis",
        "SUMDIAM" => "Sum of Diameter",
        "LPERP" => "Longest Perpendicular",
        "PPD" => "Product of Perpendicular Diameters",
        "SUMPPD" => "Sum of Products of Perpendicular Diam",
        "FDP5PS" => "FDG PET 5PS Score",
        _ => "Tumor State",
    };
    let result = |r: &TrRow<'_>| match &r.result {
        TrResult::Mm(v) => v.map(|v| fmt_dp(v, 0)),
        TrResult::State(s) => Some(s.clone()),
    };
    let mm = |r: &TrRow<'_>| matches!(r.result, TrResult::Mm(Some(_)));
    let not_done = |r: &TrRow<'_>| matches!(r.result, TrResult::Mm(None));
    build(
        "TR",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "TR")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("TRSEQ", ints(&rows, |r| Some(r.seq))),
            ("TRGRPID", chars(&rows, |r| Some(r.group.to_string()))),
            (
                "TRLNKGRP",
                chars(&rows, |r| Some(format!("A{}", r.visitnum))),
            ),
            ("TRLNKID", chars(&rows, |r| r.link.clone())),
            ("TRTESTCD", chars(&rows, |r| Some(r.testcd.to_string()))),
            ("TRTEST", chars(&rows, |r| Some(test(r.testcd).to_string()))),
            ("TRORRES", chars(&rows, result)),
            (
                "TRORRESU",
                chars(&rows, |r| {
                    mm(r).then(|| r.unit.map(str::to_string)).flatten()
                }),
            ),
            ("TRSTRESC", chars(&rows, result)),
            (
                "TRSTRESN",
                nums(&rows, |r| match r.result {
                    TrResult::Mm(v) => v,
                    TrResult::State(_) => None,
                }),
            ),
            (
                "TRSTRESU",
                chars(&rows, |r| {
                    mm(r).then(|| r.unit.map(str::to_string)).flatten()
                }),
            ),
            (
                "TRSTAT",
                chars(&rows, |r| not_done(r).then(|| "NOT DONE".to_string())),
            ),
            (
                "TRREASND",
                chars(&rows, |r| {
                    not_done(r).then(|| "IMAGE QUALITY INADEQUATE".to_string())
                }),
            ),
            ("TRMETHOD", chars(&rows, |r| Some(r.method.to_string()))),
            ("TREVAL", constant(&rows, EVALUATOR)),
            ("VISITNUM", ints(&rows, |r| Some(r.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("TRDTC", chars(&rows, |r| Some(r.date.to_string()))),
            ("TRDY", ints(&rows, |r| dy(r.date, r.rf))),
        ],
    )
}

fn rs(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        testcd: &'static str,
        scat: Option<&'static str>,
        result: String,
        visitnum: i64,
        visit: String,
        epoch: Epoch,
        date: Date,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut k = 0;
        for a in &t.tumor {
            let Some((tl, nt, ovr)) = &a.response else {
                continue;
            };
            let visit = visit_name(ctx, a.visitnum);
            let tests = if a.perpendicular.is_empty() {
                [
                    ("TRGRESP", None, Some(tl)),
                    ("NTRGRESP", None, nt.as_ref()),
                    ("OVRLRESP", None, Some(ovr)),
                ]
            } else {
                // Lugano: the CT-based and the PET-CT-based overall response.
                [
                    ("OVRLRESP", Some("CT"), Some(tl)),
                    ("OVRLRESP", Some("PET-CT"), Some(ovr)),
                    ("OVRLRESP", None, None),
                ]
            };
            for (testcd, scat, result) in tests {
                let Some(result) = result else { continue };
                k += 1;
                rows.push(Row {
                    usubjid: &s.usubjid,
                    seq: k,
                    testcd,
                    scat,
                    result: result.clone(),
                    visitnum: a.visitnum,
                    visit: visit.clone(),
                    epoch: t.epoch_of(a.date),
                    date: a.date,
                    rf: t.first_dose(),
                });
            }
        }
    }
    let test = |cd: &str| match cd {
        "TRGRESP" => "Target Response",
        "NTRGRESP" => "Non-target Response",
        _ => "Overall Response",
    };
    let lugano = rows.iter().any(|r| r.scat.is_some());
    let mut cols = vec![
        ("STUDYID", constant(&rows, &ctx.design.studyid)),
        ("DOMAIN", constant(&rows, "RS")),
        ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
        ("RSSEQ", ints(&rows, |r| Some(r.seq))),
        (
            "RSLNKGRP",
            chars(&rows, |r| Some(format!("A{}", r.visitnum))),
        ),
        ("RSTESTCD", chars(&rows, |r| Some(r.testcd.to_string()))),
        ("RSTEST", chars(&rows, |r| Some(test(r.testcd).to_string()))),
        (
            "RSCAT",
            constant(
                &rows,
                if lugano {
                    "LUGANO CLASSIFICATION"
                } else {
                    "RECIST 1.1"
                },
            ),
        ),
        ("RSORRES", chars(&rows, |r| Some(r.result.clone()))),
        ("RSSTRESC", chars(&rows, |r| Some(r.result.clone()))),
        ("RSEVAL", constant(&rows, EVALUATOR)),
        ("VISITNUM", ints(&rows, |r| Some(r.visitnum))),
        ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
        (
            "EPOCH",
            chars(&rows, |r| Some(r.epoch.as_str().to_string())),
        ),
        ("RSDTC", chars(&rows, |r| Some(r.date.to_string()))),
        ("RSDY", ints(&rows, |r| dy(r.date, r.rf))),
    ];
    if lugano {
        cols.push(("RSSCAT", chars(&rows, |r| r.scat.map(str::to_string))));
    }
    build("RS", cols)
}
