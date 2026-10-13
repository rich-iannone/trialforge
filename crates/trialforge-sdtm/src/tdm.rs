//! Trial design model: TA, TE, TV, TI, TS.

use std::collections::BTreeSet;

use trialforge_cdisc::{MetadataError, build};
use trialforge_core::Dataset;
use trialforge_model::{Epoch, VisitKind};

use crate::{Context, chars, constant, ints};

/// `ETCD` and `ELEMENT` of an epoch's element for arm `arm`.
pub fn element(ctx: &Context<'_>, epoch: Epoch, arm: usize) -> (String, String) {
    match epoch {
        Epoch::Screening => ("SCRN".into(), "Screening".into()),
        Epoch::Treatment => {
            let a = &ctx.design.arms[arm];
            let etcd: String = a.armcd.chars().take(8).collect();
            let t = &a.treatment;
            let dose = trialforge_core::math::round_dp(t.dose, 4);
            (etcd, format!("{} {} {}", t.name, dose, t.unit))
        }
        Epoch::FollowUp => ("FU".into(), "Follow-up".into()),
    }
}

/// The order of an epoch's element within every arm.
pub fn taetord(epoch: Epoch) -> i64 {
    match epoch {
        Epoch::Screening => 1,
        Epoch::Treatment => 2,
        Epoch::FollowUp => 3,
    }
}

const EPOCHS: [Epoch; 3] = [Epoch::Screening, Epoch::Treatment, Epoch::FollowUp];

pub(crate) fn emit(ctx: &Context<'_>) -> Result<Vec<Dataset>, MetadataError> {
    Ok(vec![ta(ctx)?, te(ctx)?, tv(ctx)?, ti(ctx)?, ts(ctx)?])
}

fn studyid<'a>(ctx: &Context<'a>) -> &'a str {
    &ctx.design.studyid
}

fn ta(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row {
        armcd: String,
        arm: String,
        order: i64,
        etcd: String,
        element: String,
        branch: Option<String>,
        epoch: Epoch,
    }
    let mut rows = Vec::new();
    for (i, a) in ctx.design.arms.iter().enumerate() {
        for epoch in EPOCHS {
            let (etcd, element) = element(ctx, epoch, i);
            rows.push(Row {
                armcd: a.armcd.clone(),
                arm: a.arm.clone(),
                order: taetord(epoch),
                etcd,
                element,
                branch: (epoch == Epoch::Screening).then(|| format!("Randomized to {}", a.arm)),
                epoch,
            });
        }
    }
    build(
        "TA",
        vec![
            ("STUDYID", constant(&rows, studyid(ctx))),
            ("DOMAIN", constant(&rows, "TA")),
            ("ARMCD", chars(&rows, |r| Some(r.armcd.clone()))),
            ("ARM", chars(&rows, |r| Some(r.arm.clone()))),
            ("TAETORD", ints(&rows, |r| Some(r.order))),
            ("ETCD", chars(&rows, |r| Some(r.etcd.clone()))),
            ("ELEMENT", chars(&rows, |r| Some(r.element.clone()))),
            ("TABRANCH", chars(&rows, |r| r.branch.clone())),
            ("TATRANS", chars(&rows, |_| None)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
        ],
    )
}

fn te(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let mut rows: Vec<(String, String, &str, &str)> = vec![(
        "SCRN".into(),
        "Screening".into(),
        "Informed consent",
        "First dose of study treatment",
    )];
    for i in 0..ctx.design.arms.len() {
        let (etcd, element) = element(ctx, Epoch::Treatment, i);
        rows.push((
            etcd,
            element,
            "First dose of study treatment",
            "End of study treatment",
        ));
    }
    rows.push((
        "FU".into(),
        "Follow-up".into(),
        "End of study treatment",
        "End of study participation",
    ));
    build(
        "TE",
        vec![
            ("STUDYID", constant(&rows, studyid(ctx))),
            ("DOMAIN", constant(&rows, "TE")),
            ("ETCD", chars(&rows, |r| Some(r.0.clone()))),
            ("ELEMENT", chars(&rows, |r| Some(r.1.clone()))),
            ("TESTRL", chars(&rows, |r| Some(r.2.to_string()))),
            ("TEENRL", chars(&rows, |r| Some(r.3.to_string()))),
        ],
    )
}

fn tv(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let visits = ctx
        .design
        .schedule
        .as_ref()
        .map(|s| s.visits.as_slice())
        .unwrap_or_default();
    let rule = |v: &trialforge_model::PlannedVisit| match (v.kind, v.day) {
        (VisitKind::Screening, _) => "Informed consent obtained".to_string(),
        (VisitKind::Baseline, _) => "Day 1, before first dose".to_string(),
        (VisitKind::EarlyTermination, _) => "At early termination".to_string(),
        (_, Some(d)) => format!("Day {d} (window -{} to +{} days)", v.window.0, v.window.1),
        (_, None) => "Unscheduled".to_string(),
    };
    // One record per planned visit per arm (every arm shares the schedule).
    let rows: Vec<(&trialforge_model::Arm, &trialforge_model::PlannedVisit)> = ctx
        .design
        .arms
        .iter()
        .flat_map(|a| visits.iter().map(move |v| (a, v)))
        .collect();
    build(
        "TV",
        vec![
            ("STUDYID", constant(&rows, studyid(ctx))),
            ("DOMAIN", constant(&rows, "TV")),
            ("VISITNUM", ints(&rows, |(_, v)| Some(v.visitnum))),
            ("VISIT", chars(&rows, |(_, v)| Some(v.name.clone()))),
            ("VISITDY", ints(&rows, |(_, v)| v.day.map(i64::from))),
            ("ARMCD", chars(&rows, |(a, _)| Some(a.armcd.clone()))),
            ("ARM", chars(&rows, |(a, _)| Some(a.arm.clone()))),
            ("TVSTRL", chars(&rows, |(_, v)| Some(rule(v)))),
        ],
    )
}

fn ti(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let rows = &ctx.design.info.criteria;
    build(
        "TI",
        vec![
            ("STUDYID", constant(rows, studyid(ctx))),
            ("DOMAIN", constant(rows, "TI")),
            ("IETESTCD", chars(rows, |r| Some(r.0.clone()))),
            ("IETEST", chars(rows, |r| Some(r.2.clone()))),
            ("IECAT", chars(rows, |r| Some(r.1.clone()))),
            ("TIVERS", constant(rows, "1")),
        ],
    )
}

fn ts(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let d = ctx.design;
    let info = &d.info;
    let tsparm = |cd: &str| -> String {
        let ct = trialforge_cdisc::ct();
        let pos = ct
            .codelist("TSPARMCD")
            .and_then(|c| c.terms.iter().position(|t| t == cd))
            .expect("TSPARMCD in CT");
        ct.codelist("TSPARM").expect("TSPARM in CT").terms[pos].clone()
    };
    let n = ctx.subjects.len();
    let first_consent = ctx.subjects.iter().map(|s| s.consent).min();
    let last_contact = ctx.timelines.iter().map(|t| t.last_contact).max();
    let countries: BTreeSet<&str> = d.sites.iter().map(|s| s.country.as_str()).collect();
    let investigational = d
        .arms
        .iter()
        .find(|a| a.treatment.dose > 0.0)
        .unwrap_or(&d.arms[0]);
    let yes_no = |b: bool| if b { "Y" } else { "N" }.to_string();

    // (TSPARMCD, value, reference terminology)
    let mut params: Vec<(&str, String, Option<&str>)> = vec![
        ("TITLE", info.title.clone(), None),
        ("SPONSOR", info.sponsor.clone(), None),
        ("TPHASE", info.phase.clone(), Some("CDISC CT")),
        ("TTYPE", "EFFICACY".into(), Some("CDISC CT")),
        ("STYPE", "INTERVENTIONAL".into(), Some("CDISC CT")),
        ("INTMODEL", "PARALLEL".into(), Some("CDISC CT")),
        ("INTTYPE", "DRUG".into(), Some("CDISC CT")),
        ("RANDOM", "Y".into(), Some("CDISC CT")),
        ("TBLIND", info.blinding.clone(), Some("CDISC CT")),
        ("TCNTRL", info.control.clone(), Some("CDISC CT")),
        ("ADDON", yes_no(info.add_on), Some("CDISC CT")),
        ("ADAPT", "N".into(), Some("CDISC CT")),
        ("HLTSUBJI", "N".into(), Some("CDISC CT")),
        ("TDIGRP", info.indication.clone(), None),
        ("THERAREA", info.therapeutic_area.clone(), None),
        ("OBJPRIM", info.primary_objective.clone(), None),
        ("TRT", investigational.treatment.name.clone(), None),
        (
            "AGEMIN",
            format!("P{}Y", d.population.age_min),
            Some("ISO 8601"),
        ),
        (
            "AGEMAX",
            format!("P{}Y", d.population.age_max),
            Some("ISO 8601"),
        ),
        ("SEXPOP", "BOTH".into(), Some("CDISC CT")),
        ("PLANSUB", n.to_string(), None),
        ("ACTSUB", n.to_string(), None),
        ("NARMS", d.arms.len().to_string(), None),
        ("SDTIGVER", "3.4".into(), None),
        ("SDTMVER", "2.0".into(), None),
    ];
    if let Some(sch) = &d.schedule {
        params.push((
            "LENGTH",
            format!("P{}D", sch.treatment_days + sch.followup_days),
            Some("ISO 8601"),
        ));
    }
    if let Some(dc) = first_consent {
        params.push(("SSTDTC", dc.to_string(), Some("ISO 8601")));
    }
    if d.schedule.is_some() {
        if let Some(lc) = last_contact {
            params.push(("SENDTC", lc.to_string(), Some("ISO 8601")));
        }
    } else {
        params.push(("DCUTDTC", d.cutoff.to_string(), Some("ISO 8601")));
        params.push(("DCUTDESC", "PRIMARY ANALYSIS".into(), None));
    }
    for c in &countries {
        params.push(("FCNTRY", (*c).to_string(), Some("ISO 3166-1 alpha-3")));
    }

    struct Row {
        seq: i64,
        cd: String,
        parm: String,
        val: String,
        reference: Option<String>,
    }
    let mut rows: Vec<Row> = Vec::new();
    for (cd, val, reference) in params {
        let seq = rows.iter().filter(|r| r.cd == cd).count() as i64 + 1;
        rows.push(Row {
            seq,
            cd: cd.into(),
            parm: tsparm(cd),
            val,
            reference: reference.map(str::to_string),
        });
    }
    let ct_version = trialforge_cdisc::ct().version;
    build(
        "TS",
        vec![
            ("STUDYID", constant(&rows, studyid(ctx))),
            ("DOMAIN", constant(&rows, "TS")),
            ("TSSEQ", ints(&rows, |r| Some(r.seq))),
            ("TSPARMCD", chars(&rows, |r| Some(r.cd.clone()))),
            ("TSPARM", chars(&rows, |r| Some(r.parm.clone()))),
            ("TSVAL", chars(&rows, |r| Some(r.val.clone()))),
            ("TSVALNF", chars(&rows, |_| None)),
            ("TSVALCD", chars(&rows, |_| None)),
            ("TSVCDREF", chars(&rows, |r| r.reference.clone())),
            (
                "TSVCDVER",
                chars(&rows, |r| match r.reference.as_deref() {
                    Some("CDISC CT") => Some(ct_version.to_string()),
                    _ => None,
                }),
            ),
        ],
    )
}
