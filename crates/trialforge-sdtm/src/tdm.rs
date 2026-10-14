//! Trial design model: TA, TE, TV, TI, TS.

use std::collections::BTreeSet;

use trialforge_cdisc::{MetadataError, build};
use trialforge_core::Dataset;
use trialforge_model::{Epoch, VisitKind};

use crate::{Context, chars, constant, ints};

/// The element name of a treatment, e.g. `DRUG A 10 mg`.
fn treatment_element(t: &trialforge_model::Treatment) -> String {
    let dose = trialforge_core::math::round_dp(t.dose, 4);
    format!("{} {} {}", t.name, dose, t.unit)
}

/// The elements of arm `arm`, in order: (epoch, `ETCD`, `ELEMENT`).
pub fn arm_elements(ctx: &Context<'_>, arm: usize) -> Vec<(Epoch, String, String)> {
    let a = &ctx.design.arms[arm];
    let mut out = vec![(
        Epoch::Screening,
        "SCRN".to_string(),
        "Screening".to_string(),
    )];
    if a.periods.is_empty() {
        let etcd: String = a.armcd.chars().take(8).collect();
        out.push((Epoch::Treatment, etcd, treatment_element(&a.treatment)));
        if let Some((_, element)) = &a.switch_to {
            out.push((Epoch::Continuation, "XOVER".into(), element.clone()));
        }
    } else {
        let n = a.periods.len();
        for (k, p) in a.periods.iter().enumerate() {
            out.push((
                Epoch::Period(k as u8 + 1),
                p.etcd.clone(),
                p.element.clone(),
            ));
            if k + 1 < n {
                let w = if n == 2 { 0 } else { k as u8 + 1 };
                out.push((Epoch::Washout(w), "WASHOUT".into(), "Washout".into()));
            }
        }
    }
    out.push((Epoch::FollowUp, "FU".into(), "Follow-up".into()));
    out
}

/// `ETCD` and `ELEMENT` of an epoch's element for arm `arm`.
pub fn element(ctx: &Context<'_>, epoch: Epoch, arm: usize) -> (String, String) {
    arm_elements(ctx, arm)
        .into_iter()
        .find(|(e, _, _)| *e == epoch)
        .map(|(_, c, n)| (c, n))
        .expect("every timeline epoch is an element of the arm")
}

/// The order of an epoch's element within arm `arm` (`TAETORD`).
pub fn taetord(ctx: &Context<'_>, epoch: Epoch, arm: usize) -> i64 {
    arm_elements(ctx, arm)
        .iter()
        .position(|(e, _, _)| *e == epoch)
        .map_or(0, |p| p as i64 + 1)
}

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
        for (k, (epoch, etcd, element)) in arm_elements(ctx, i).into_iter().enumerate() {
            rows.push(Row {
                armcd: a.armcd.clone(),
                arm: a.arm.clone(),
                order: k as i64 + 1,
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
    let crossover = ctx.design.arms.iter().any(|a| !a.periods.is_empty());
    let mut rows: Vec<(String, String, &str, &str)> = vec![(
        "SCRN".into(),
        "Screening".into(),
        "Informed consent",
        "First dose of study treatment",
    )];
    let mut washout = false;
    for i in 0..ctx.design.arms.len() {
        for (epoch, etcd, element) in arm_elements(ctx, i) {
            if rows.iter().any(|r| r.0 == etcd) {
                continue;
            }
            match epoch {
                Epoch::Treatment => rows.push((
                    etcd,
                    element,
                    "First dose of study treatment",
                    "End of study treatment",
                )),
                Epoch::Period(_) => rows.push((
                    etcd,
                    element,
                    "Dose of the period's treatment",
                    "End of the treatment period",
                )),
                Epoch::Continuation => rows.push((
                    etcd,
                    element,
                    "First dose after switching on progression",
                    "End of the switch treatment",
                )),
                Epoch::Washout(_) => washout = true,
                _ => {}
            }
        }
    }
    if washout {
        rows.push((
            "WASHOUT".into(),
            "Washout".into(),
            "End of the treatment period",
            "Dose of the next period's treatment",
        ));
    }
    rows.push((
        "FU".into(),
        "Follow-up".into(),
        if crossover {
            "End of the last treatment period"
        } else {
            "End of study treatment"
        },
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
    ];
    for t in &info.trial_types {
        params.push(("TTYPE", t.clone(), Some("CDISC CT")));
    }
    params.extend(vec![
        ("STYPE", "INTERVENTIONAL".into(), Some("CDISC CT")),
        (
            "INTMODEL",
            info.intervention_model.clone(),
            Some("CDISC CT"),
        ),
        ("INTTYPE", info.intervention_type.clone(), Some("CDISC CT")),
        ("RANDOM", "Y".into(), Some("CDISC CT")),
        ("TBLIND", info.blinding.clone(), Some("CDISC CT")),
        ("TCNTRL", info.control.clone(), Some("CDISC CT")),
        ("ADDON", yes_no(info.add_on), Some("CDISC CT")),
        ("ADAPT", "N".into(), Some("CDISC CT")),
        ("HLTSUBJI", yes_no(info.healthy_subjects), Some("CDISC CT")),
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
        ("PLANSUB", d.n.to_string(), None),
        ("ACTSUB", n.to_string(), None),
        ("NARMS", d.arms.len().to_string(), None),
        ("SDTIGVER", "3.4".into(), None),
        ("SDTMVER", "2.0".into(), None),
    ]);
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
    if d.schedule.is_some()
        && let Some(lc) = last_contact
    {
        params.push(("SENDTC", lc.to_string(), Some("ISO 8601")));
    }
    // Event-driven studies have a data cutoff; an interim snapshot reports its own.
    if d.schedule.is_none() || d.interim {
        params.push(("DCUTDTC", d.cutoff.to_string(), Some("ISO 8601")));
        let desc = if d.interim {
            "INTERIM ANALYSIS"
        } else {
            "PRIMARY ANALYSIS"
        };
        params.push(("DCUTDESC", desc.into(), None));
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
