//! Subject-level domains: DM, SE, SV, EX, CM, AE, SUPPAE, DS, MH, LB, VS, FT, OE, SC, CE.

use trialforge_cdisc::{MetadataError, TEAE_WINDOW_DAYS, build, dictionaries};
use trialforge_core::{Dataset, Date};
use trialforge_model::{Epoch, Subject};
use trialforge_processes::reference;
use trialforge_processes::timeline::{DispositionKind, Finding, FindingDomain, Timeline};

use crate::tdm::{element, taetord};
use crate::{Context, chars, constant, dy, fmt_dp, ints, nums, yn};

pub(crate) fn emit(ctx: &Context<'_>) -> Result<Vec<Dataset>, MetadataError> {
    let (ae, suppae) = ae(ctx)?;
    Ok(vec![
        dm(ctx)?,
        se(ctx)?,
        sv(ctx)?,
        ex(ctx)?,
        cm(ctx)?,
        ae,
        ds(ctx)?,
        mh(ctx)?,
        lb(ctx)?,
        vs(ctx)?,
        ft(ctx)?,
        oe(ctx)?,
        re(ctx)?,
        crate::vaccine::is(ctx)?,
        sc(ctx)?,
        if ctx.timelines.iter().any(|t| !t.diaries.is_empty()) {
            crate::vaccine::ce(ctx)?
        } else {
            ce(ctx)?
        },
        crate::vaccine::face(ctx)?,
        suppae,
    ])
}

fn iso(d: Date) -> Option<String> {
    Some(d.to_string())
}

fn dm(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let d = ctx.design;
    let pop = &d.population;
    let rows: Vec<(&Subject, &Timeline)> = ctx.sorted().collect();
    let arm = |s: &Subject| &d.arms[s.arm];
    let mut cols = vec![
        ("STUDYID", constant(&rows, &d.studyid)),
        ("DOMAIN", constant(&rows, "DM")),
        ("USUBJID", chars(&rows, |(s, _)| Some(s.usubjid.clone()))),
        ("SUBJID", chars(&rows, |(s, _)| Some(s.subjid.clone()))),
        ("RFSTDTC", chars(&rows, |(_, t)| iso(t.first_dose()))),
        ("RFENDTC", chars(&rows, |(_, t)| iso(t.last_contact))),
        ("RFXSTDTC", chars(&rows, |(_, t)| iso(t.first_dose()))),
        ("RFXENDTC", chars(&rows, |(_, t)| iso(t.last_dose()))),
        ("RFICDTC", chars(&rows, |(s, _)| iso(s.consent))),
        ("RFPENDTC", chars(&rows, |(_, t)| iso(t.last_contact))),
        ("DTHDTC", chars(&rows, |(_, t)| t.death.and_then(iso))),
        (
            "DTHFL",
            chars(&rows, |(_, t)| t.death.map(|_| "Y".to_string())),
        ),
        (
            "SITEID",
            chars(&rows, |(s, _)| Some(d.sites[s.site].siteid.clone())),
        ),
        ("AGE", ints(&rows, |(s, _)| Some(s.age))),
        ("AGEU", constant(&rows, "YEARS")),
        (
            "SEX",
            chars(&rows, |(s, _)| Some(pop.sex.levels[s.sex].0.clone())),
        ),
        (
            "RACE",
            chars(&rows, |(s, _)| Some(pop.race.levels[s.race].0.clone())),
        ),
        (
            "ETHNIC",
            chars(&rows, |(s, _)| Some(pop.ethnic.levels[s.ethnic].0.clone())),
        ),
        ("ARMCD", chars(&rows, |(s, _)| Some(arm(s).armcd.clone()))),
        ("ARM", chars(&rows, |(s, _)| Some(arm(s).arm.clone()))),
        (
            "ACTARMCD",
            chars(&rows, |(s, _)| Some(arm(s).armcd.clone())),
        ),
        ("ACTARM", chars(&rows, |(s, _)| Some(arm(s).arm.clone()))),
        (
            "COUNTRY",
            chars(&rows, |(s, _)| Some(d.sites[s.site].country.clone())),
        ),
        ("DMDTC", chars(&rows, |(s, _)| iso(s.consent))),
        ("DMDY", ints(&rows, |(s, t)| dy(s.consent, t.first_dose()))),
    ];
    if rows.iter().any(|(_, t)| t.birth_date.is_some()) {
        cols.push(("BRTHDTC", chars(&rows, |(_, t)| t.birth_date.and_then(iso))));
    }
    build("DM", cols)
}

fn se(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row {
        usubjid: String,
        seq: i64,
        etcd: String,
        element: String,
        order: i64,
        epoch: Epoch,
        start: Date,
        end: Date,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, e) in t.elements.iter().enumerate() {
            let (etcd, element) = element(ctx, e.epoch, s.arm);
            let order = taetord(ctx, e.epoch, s.arm);
            rows.push(Row {
                order,
                usubjid: s.usubjid.clone(),
                seq: k as i64 + 1,
                etcd,
                element,
                epoch: e.epoch,
                start: e.start,
                end: e.end,
                rf: t.first_dose(),
            });
        }
    }
    build(
        "SE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "SE")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.clone()))),
            ("SESEQ", ints(&rows, |r| Some(r.seq))),
            ("ETCD", chars(&rows, |r| Some(r.etcd.clone()))),
            ("ELEMENT", chars(&rows, |r| Some(r.element.clone()))),
            ("TAETORD", ints(&rows, |r| Some(r.order))),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("SESTDTC", chars(&rows, |r| iso(r.start))),
            ("SEENDTC", chars(&rows, |r| iso(r.end))),
            ("SESTDY", ints(&rows, |r| dy(r.start, r.rf))),
            ("SEENDY", ints(&rows, |r| dy(r.end, r.rf))),
        ],
    )
}

fn sv(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        visitnum: i64,
        visit: String,
        visitdy: Option<i64>,
        epoch: Epoch,
        date: Date,
        end: Date,
        rf: Date,
    }
    let schedule = ctx.design.schedule.as_ref();
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for v in &t.visits {
            let pv = schedule
                .and_then(|sch| sch.visit(v.visitnum))
                .expect("visits come from the schedule");
            rows.push(Row {
                usubjid: &s.usubjid,
                visitnum: v.visitnum,
                visit: pv.name.clone(),
                visitdy: pv.day.map(i64::from),
                epoch: t.epoch_of(v.date),
                date: v.date,
                end: v.end.unwrap_or(v.date),
                rf: t.first_dose(),
            });
        }
    }
    rows.sort_by(|a, b| (a.usubjid, a.visitnum).cmp(&(b.usubjid, b.visitnum)));
    build(
        "SV",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "SV")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("VISITNUM", ints(&rows, |r| Some(r.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("SVSTDTC", chars(&rows, |r| iso(r.date))),
            ("SVENDTC", chars(&rows, |r| iso(r.end))),
            ("SVSTDY", ints(&rows, |r| dy(r.date, r.rf))),
            ("SVENDY", ints(&rows, |r| dy(r.end, r.rf))),
        ],
    )
}

fn ex(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        s: &'a Subject,
        seq: i64,
        start: Date,
        end: Date,
        time: Option<trialforge_core::DateTime>,
        dose: f64,
        epoch: Epoch,
        rf: Date,
        visit: Option<&'a trialforge_model::PlannedVisit>,
        e: &'a trialforge_processes::timeline::Exposure,
    }
    let schedule = ctx.design.schedule.as_ref();
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, e) in t.exposure.iter().enumerate() {
            let study_day = e.start.days_since(t.first_dose()) + 1;
            rows.push(Row {
                s,
                seq: k as i64 + 1,
                start: e.start,
                end: e.end,
                time: e.time,
                dose: e.dose,
                epoch: t.epoch_of(e.start),
                rf: t.first_dose(),
                // Timed administrations belong to their recorded visit, else the dosing day's.
                visit: match e.visitnum {
                    Some(n) => schedule.and_then(|sch| sch.visit(n)),
                    None => e
                        .time
                        .and_then(|_| schedule?.visits.iter().find(|v| v.day == Some(study_day))),
                },
                e,
            });
        }
    }
    let timed = rows.iter().any(|r| r.time.is_some());
    let dtc = |r: &Row<'_>, d: Date| match r.time {
        Some(t) => Some(t.to_string()),
        None => iso(d),
    };
    let trt = |r: &Row<'_>| {
        let arm = &ctx.design.arms[r.s.arm];
        match (&arm.switch_to, r.e.switched) {
            (Some((t, _)), true) => t,
            _ => &arm.treatment,
        }
    };
    let mut cols = vec![
        ("STUDYID", constant(&rows, &ctx.design.studyid)),
        ("DOMAIN", constant(&rows, "EX")),
        ("USUBJID", chars(&rows, |r| Some(r.s.usubjid.clone()))),
        ("EXSEQ", ints(&rows, |r| Some(r.seq))),
        ("EXTRT", chars(&rows, |r| Some(trt(r).name.clone()))),
        ("EXDOSE", nums(&rows, |r| Some(r.dose))),
        ("EXDOSU", chars(&rows, |r| Some(trt(r).unit.clone()))),
        ("EXDOSFRM", chars(&rows, |r| Some(trt(r).form.clone()))),
        ("EXDOSFRQ", chars(&rows, |r| Some(trt(r).frequency.clone()))),
        ("EXROUTE", chars(&rows, |r| Some(trt(r).route.clone()))),
        (
            "EPOCH",
            chars(&rows, |r| Some(r.epoch.as_str().to_string())),
        ),
        ("EXSTDTC", chars(&rows, |r| dtc(r, r.start))),
        ("EXENDTC", chars(&rows, |r| dtc(r, r.end))),
        ("EXSTDY", ints(&rows, |r| dy(r.start, r.rf))),
        ("EXENDY", ints(&rows, |r| dy(r.end, r.rf))),
    ];
    if rows.iter().any(|r| r.e.location.is_some()) {
        cols.push((
            "EXLOC",
            chars(&rows, |r| r.e.location.map(|l| l.0.to_string())),
        ));
        cols.push((
            "EXLAT",
            chars(&rows, |r| r.e.location.map(|l| l.1.to_string())),
        ));
    }
    if rows.iter().any(|r| r.e.group.is_some()) {
        cols.push(("EXLNKGRP", chars(&rows, |r| r.e.group.clone())));
        cols.push((
            "EXLNKID",
            chars(&rows, |r| {
                let g = r.e.group.as_ref()?;
                Some(match r.e.location {
                    Some((loc, lat)) => format!("{g}-{loc}-{lat}"),
                    None => g.clone(),
                })
            }),
        ));
    }
    if timed {
        cols.push(("VISITNUM", ints(&rows, |r| r.visit.map(|v| v.visitnum))));
        cols.push(("VISIT", chars(&rows, |r| r.visit.map(|v| v.name.clone()))));
        cols.push((
            "VISITDY",
            ints(&rows, |r| r.visit.and_then(|v| v.day.map(i64::from))),
        ));
    }
    build("EX", cols)
}

fn cm(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        m: &'a trialforge_processes::timeline::Medication,
        start: Option<Date>,
        t: &'a Timeline,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut meds: Vec<_> = t.medications.iter().collect();
        meds.sort_by(|a, b| (&a.start, &a.ingredient).cmp(&(&b.start, &b.ingredient)));
        for (k, m) in meds.into_iter().enumerate() {
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                m,
                start: Date::parse_iso(&m.start),
                t,
            });
        }
    }
    let ing = |r: &Row<'_>| {
        dictionaries::ingredient(&r.m.ingredient).expect("medications come from the dictionary")
    };
    build(
        "CM",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "CM")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("CMSEQ", ints(&rows, |r| Some(r.seq))),
            ("CMTRT", chars(&rows, |r| Some(r.m.ingredient.clone()))),
            ("CMDECOD", chars(&rows, |r| Some(r.m.ingredient.clone()))),
            ("CMCLAS", chars(&rows, |r| Some(ing(r).atc_class.clone()))),
            (
                "CMINDC",
                chars(&rows, |r| Some(r.m.indication.to_uppercase())),
            ),
            ("CMDOSE", nums(&rows, |r| Some(ing(r).dose))),
            ("CMDOSU", chars(&rows, |r| Some(ing(r).unit.clone()))),
            ("CMDOSFRQ", chars(&rows, |r| Some(ing(r).frequency.clone()))),
            ("CMROUTE", chars(&rows, |r| Some(ing(r).route.clone()))),
            (
                "EPOCH",
                chars(&rows, |r| {
                    r.start.map(|d| r.t.epoch_of(d).as_str().to_string())
                }),
            ),
            ("CMSTDTC", chars(&rows, |r| Some(r.m.start.clone()))),
            ("CMENDTC", chars(&rows, |r| r.m.end.and_then(iso))),
            (
                "CMSTDY",
                ints(&rows, |r| r.start.and_then(|d| dy(d, r.t.first_dose()))),
            ),
            (
                "CMENDY",
                ints(&rows, |r| r.m.end.and_then(|d| dy(d, r.t.first_dose()))),
            ),
            (
                "CMENRTPT",
                chars(&rows, |r| r.m.end.is_none().then(|| "ONGOING".to_string())),
            ),
            (
                "CMENTPT",
                chars(&rows, |r| {
                    r.m.end.is_none().then(|| "END OF STUDY".to_string())
                }),
            ),
        ],
    )
}

/// Whether an AE starting on `start` is treatment-emergent for this subject.
pub fn treatment_emergent(t: &Timeline, start: Date) -> bool {
    start >= t.first_dose() && start <= t.last_dose().add_days(TEAE_WINDOW_DAYS)
}

fn ae(ctx: &Context<'_>) -> Result<(Dataset, Dataset), MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        a: &'a trialforge_processes::timeline::AdverseEvent,
        soc: &'static str,
        t: &'a Timeline,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut events: Vec<_> = t.adverse_events.iter().collect();
        events.sort_by(|a, b| (a.start, &a.pt).cmp(&(b.start, &b.pt)));
        for (k, a) in events.into_iter().enumerate() {
            let soc = &dictionaries::ae_term(&a.pt)
                .expect("AE terms come from the dictionary")
                .soc;
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                a,
                soc,
                t,
            });
        }
    }
    let crit = |r: &Row<'_>, f: fn(&trialforge_processes::timeline::SeriousCriteria) -> bool| {
        r.a.serious.as_ref().and_then(|c| yn(f(c)))
    };
    let ae = build(
        "AE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "AE")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("AESEQ", ints(&rows, |r| Some(r.seq))),
            ("AETERM", chars(&rows, |r| Some(r.a.pt.to_uppercase()))),
            ("AEDECOD", chars(&rows, |r| Some(r.a.pt.clone()))),
            ("AEBODSYS", chars(&rows, |r| Some(r.soc.to_string()))),
            ("AESOC", chars(&rows, |r| Some(r.soc.to_string()))),
            (
                "AESEV",
                chars(&rows, |r| Some(r.a.severity.as_str().to_string())),
            ),
            ("AESER", chars(&rows, |r| yn(r.a.serious.is_some()))),
            ("AEACN", chars(&rows, |r| Some(r.a.action.to_string()))),
            ("AEREL", chars(&rows, |r| Some(r.a.causality.to_string()))),
            ("AEOUT", chars(&rows, |r| Some(r.a.outcome.to_string()))),
            ("AESCONG", chars(&rows, |r| crit(r, |_| false))),
            ("AESDISAB", chars(&rows, |r| crit(r, |c| c.disability))),
            ("AESDTH", chars(&rows, |r| crit(r, |c| c.death))),
            ("AESHOSP", chars(&rows, |r| crit(r, |c| c.hospitalization))),
            ("AESLIFE", chars(&rows, |r| crit(r, |c| c.life_threatening))),
            (
                "AESMIE",
                chars(&rows, |r| crit(r, |c| c.medically_important)),
            ),
            ("AECONTRT", chars(&rows, |r| yn(r.a.treated))),
            (
                "EPOCH",
                chars(&rows, |r| {
                    Some(r.t.epoch_of(r.a.start).as_str().to_string())
                }),
            ),
            ("AESTDTC", chars(&rows, |r| iso(r.a.start))),
            ("AEENDTC", chars(&rows, |r| r.a.end.and_then(iso))),
            ("AESTDY", ints(&rows, |r| dy(r.a.start, r.t.first_dose()))),
            (
                "AEENDY",
                ints(&rows, |r| r.a.end.and_then(|d| dy(d, r.t.first_dose()))),
            ),
            (
                "AEENRTPT",
                chars(&rows, |r| r.a.end.is_none().then(|| "ONGOING".to_string())),
            ),
            (
                "AEENTPT",
                chars(&rows, |r| {
                    r.a.end.is_none().then(|| "END OF STUDY".to_string())
                }),
            ),
        ],
    )?;
    // One AETRTEM qualifier per AE, and AEDLT for AEs in a DLT observation window.
    let supp: Vec<(&Row<'_>, bool)> = rows
        .iter()
        .flat_map(|r| std::iter::once((r, false)).chain(r.a.dlt.map(|_| (r, true))))
        .collect();
    let suppae = build(
        "SUPPAE",
        vec![
            ("STUDYID", constant(&supp, &ctx.design.studyid)),
            ("RDOMAIN", constant(&supp, "AE")),
            (
                "USUBJID",
                chars(&supp, |(r, _)| Some(r.usubjid.to_string())),
            ),
            ("IDVAR", constant(&supp, "AESEQ")),
            ("IDVARVAL", chars(&supp, |(r, _)| Some(r.seq.to_string()))),
            (
                "QNAM",
                chars(&supp, |(_, d)| {
                    Some(if *d { "AEDLT" } else { "AETRTEM" }.to_string())
                }),
            ),
            (
                "QLABEL",
                chars(&supp, |(_, d)| {
                    Some(
                        if *d {
                            "Dose Limiting Toxicity"
                        } else {
                            "Treatment Emergent Flag"
                        }
                        .to_string(),
                    )
                }),
            ),
            (
                "QVAL",
                chars(&supp, |(r, d)| {
                    if *d {
                        r.a.dlt.and_then(yn)
                    } else {
                        yn(treatment_emergent(r.t, r.a.start))
                    }
                }),
            ),
            (
                "QORIG",
                chars(&supp, |(_, d)| {
                    Some(if *d { "Assigned" } else { "Derived" }.to_string())
                }),
            ),
            (
                "QEVAL",
                chars(&supp, |(_, d)| d.then(|| "INVESTIGATOR".to_string())),
            ),
        ],
    )?;
    Ok((ae, suppae))
}

fn ds(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        decod: String,
        kind: DispositionKind,
        epoch: Epoch,
        date: Date,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut events: Vec<_> = t.disposition.iter().collect();
        events.sort_by_key(|e| {
            (
                e.date,
                e.kind != DispositionKind::Milestone,
                e.kind == DispositionKind::EndOfStudy,
            )
        });
        for (k, e) in events.into_iter().enumerate() {
            // Milestones on or before first dose belong to screening.
            let epoch = if e.kind == DispositionKind::Milestone && e.date <= t.first_dose() {
                Epoch::Screening
            } else {
                t.epoch_of(e.date)
            };
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                decod: e.decod.clone(),
                kind: e.kind,
                epoch,
                date: e.date,
                rf: t.first_dose(),
            });
        }
    }
    build(
        "DS",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "DS")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("DSSEQ", ints(&rows, |r| Some(r.seq))),
            ("DSTERM", chars(&rows, |r| Some(r.decod.clone()))),
            ("DSDECOD", chars(&rows, |r| Some(r.decod.clone()))),
            (
                "DSCAT",
                chars(&rows, |r| {
                    Some(
                        if r.kind == DispositionKind::Milestone {
                            "PROTOCOL MILESTONE"
                        } else {
                            "DISPOSITION EVENT"
                        }
                        .into(),
                    )
                }),
            ),
            (
                "DSSCAT",
                chars(&rows, |r| match r.kind {
                    DispositionKind::Milestone => None,
                    DispositionKind::EndOfTreatment => Some("END OF TREATMENT".into()),
                    DispositionKind::EndOfStudy => Some("END OF STUDY".into()),
                }),
            ),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("DSSTDTC", chars(&rows, |r| iso(r.date))),
            ("DSSTDY", ints(&rows, |r| dy(r.date, r.rf))),
        ],
    )
}

fn mh(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        h: &'a trialforge_processes::timeline::HistoryItem,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut items: Vec<_> = t.history.iter().collect();
        items.sort_by(|a, b| (&a.category, &a.start, &a.pt).cmp(&(&b.category, &b.start, &b.pt)));
        for (k, h) in items.into_iter().enumerate() {
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                h,
            });
        }
    }
    let soc = |r: &Row<'_>| {
        dictionaries::ae_term(&r.h.pt)
            .expect("MH terms come from the dictionary")
            .soc
            .clone()
    };
    build(
        "MH",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "MH")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("MHSEQ", ints(&rows, |r| Some(r.seq))),
            ("MHTERM", chars(&rows, |r| Some(r.h.pt.to_uppercase()))),
            ("MHDECOD", chars(&rows, |r| Some(r.h.pt.clone()))),
            ("MHBODSYS", chars(&rows, |r| Some(soc(r)))),
            ("MHCAT", chars(&rows, |r| Some(r.h.category.clone()))),
            ("MHSTDTC", chars(&rows, |r| Some(r.h.start.clone()))),
            (
                "MHENRTPT",
                chars(&rows, |r| {
                    Some(if r.h.ongoing { "ONGOING" } else { "BEFORE" }.to_string())
                }),
            ),
            ("MHENTPT", constant(&rows, "SCREENING")),
        ],
    )
}

/// One reported finding, with everything an LB or VS row needs.
pub(crate) struct FindingRow<'a> {
    pub usubjid: &'a str,
    pub seq: i64,
    pub f: &'a Finding,
    pub visit: String,
    pub visitdy: Option<i64>,
    pub epoch: Epoch,
    pub rf: Date,
    pub baseline: bool,
    pub female: bool,
}

pub(crate) fn finding_rows<'a>(
    ctx: &'a Context<'_>,
    domain: FindingDomain,
    order: impl Fn(&str) -> usize,
) -> Vec<FindingRow<'a>> {
    let schedule = ctx.design.schedule.as_ref();
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut fs: Vec<&Finding> = t.findings.iter().filter(|f| f.domain == domain).collect();
        fs.sort_by(|a, b| {
            (order(&a.testcd), a.laterality, a.visitnum).cmp(&(
                order(&b.testcd),
                b.laterality,
                b.visitnum,
            ))
        });
        // Last observation on or before first dose, per test (and eye).
        let baseline: Vec<Option<usize>> = fs
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let last = fs
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| {
                        g.testcd == f.testcd
                            && g.laterality == f.laterality
                            && g.date <= t.first_dose()
                            && !g.value.is_nan()
                    })
                    .map(|(j, _)| j)
                    .next_back();
                (last == Some(i)).then_some(i)
            })
            .collect();
        let female = ctx.design.population.sex.levels[s.sex].0 == "F";
        for (k, f) in fs.iter().enumerate() {
            let pv = schedule
                .and_then(|sch| sch.visit(f.visitnum))
                .expect("findings are taken at planned visits");
            rows.push(FindingRow {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                f,
                visit: pv.name.clone(),
                visitdy: pv.day.map(i64::from),
                epoch: t.epoch_of(f.date),
                rf: t.first_dose(),
                baseline: baseline[k].is_some(),
                female,
            });
        }
    }
    rows
}

fn lb(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let panel = reference::lab_panel();
    let pos = |cd: &str| {
        panel
            .iter()
            .position(|a| a.testcd == cd)
            .unwrap_or(usize::MAX)
    };
    let rows = finding_rows(ctx, FindingDomain::Lb, pos);
    let a = |r: &FindingRow<'_>| &panel[pos(&r.f.testcd)];
    // Reported in conventional units; the standard value derives from the rounded report.
    let conv = |r: &FindingRow<'_>| {
        let an = a(r);
        trialforge_core::math::round_dp((r.f.value - an.offset) / an.factor, an.conv_dp)
    };
    let si = |r: &FindingRow<'_>| {
        let an = a(r);
        trialforge_core::math::round_dp(conv(r) * an.factor + an.offset, an.si_dp)
    };
    let range = |r: &FindingRow<'_>| if r.female { a(r).range_f } else { a(r).range_m };
    let conv_range =
        |r: &FindingRow<'_>, x: f64| fmt_dp((x - a(r).offset) / a(r).factor, a(r).conv_dp);
    let nrind = |r: &FindingRow<'_>| {
        let (lo, hi) = range(r);
        let v = si(r);
        Some(
            if v < lo {
                "LOW"
            } else if v > hi {
                "HIGH"
            } else {
                "NORMAL"
            }
            .to_string(),
        )
    };
    build(
        "LB",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "LB")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("LBSEQ", ints(&rows, |r| Some(r.seq))),
            ("LBTESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            ("LBTEST", chars(&rows, |r| Some(a(r).test.clone()))),
            ("LBCAT", chars(&rows, |r| Some(a(r).category.clone()))),
            (
                "LBORRES",
                chars(&rows, |r| Some(fmt_dp(conv(r), a(r).conv_dp))),
            ),
            ("LBORRESU", chars(&rows, |r| Some(a(r).conv_unit.clone()))),
            (
                "LBORNRLO",
                chars(&rows, |r| Some(conv_range(r, range(r).0))),
            ),
            (
                "LBORNRHI",
                chars(&rows, |r| Some(conv_range(r, range(r).1))),
            ),
            (
                "LBSTRESC",
                chars(&rows, |r| Some(fmt_dp(si(r), a(r).si_dp))),
            ),
            ("LBSTRESN", nums(&rows, |r| Some(si(r)))),
            ("LBSTRESU", chars(&rows, |r| Some(a(r).si_unit.clone()))),
            ("LBSTNRLO", nums(&rows, |r| Some(range(r).0))),
            ("LBSTNRHI", nums(&rows, |r| Some(range(r).1))),
            ("LBNRIND", chars(&rows, nrind)),
            ("LBSPEC", chars(&rows, |r| Some(a(r).specimen.clone()))),
            (
                "LBLOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            (
                "LBBLFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("LBDTC", chars(&rows, |r| iso(r.f.date))),
            ("LBDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}

fn vs(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let vitals: Vec<&reference::VitalSign> = reference::vital_signs()
        .iter()
        .chain(reference::extra_vital_signs())
        .collect();
    let pos = |cd: &str| {
        vitals
            .iter()
            .position(|v| v.testcd == cd)
            .unwrap_or(usize::MAX)
    };
    let found = finding_rows(ctx, FindingDomain::Vs, pos);
    let v = |testcd: &str| vitals[pos(testcd)];
    /// One VS record: a scheduled measurement or a reactogenicity diary temperature.
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        testcd: String,
        value: f64,
        dp: u32,
        baseline: bool,
        visit: Option<(i64, String, Option<i64>)>,
        epoch: Epoch,
        date: Date,
        rf: Date,
        diary: Option<(u8, usize, String)>,
    }
    let mut rows: Vec<Row<'_>> = Vec::new();
    let mut k = 0;
    for (s, t) in ctx.sorted() {
        let mut seq = 0;
        while k < found.len() && found[k].usubjid == s.usubjid {
            let r = &found[k];
            seq = r.seq;
            rows.push(Row {
                usubjid: r.usubjid,
                seq: r.seq,
                testcd: r.f.testcd.clone(),
                value: r.f.value,
                dp: r.f.dp.unwrap_or(v(&r.f.testcd).dp),
                baseline: r.baseline,
                visit: Some((r.f.visitnum, r.visit.clone(), r.visitdy)),
                epoch: r.epoch,
                date: r.f.date,
                rf: r.rf,
                diary: None,
            });
            k += 1;
        }
        for temp in crate::vaccine::temperatures(t) {
            seq += 1;
            rows.push(Row {
                usubjid: &s.usubjid,
                seq,
                testcd: "TEMP".into(),
                value: temp.value,
                dp: 1,
                baseline: false,
                visit: None,
                epoch: temp.t.epoch_of(temp.date),
                date: temp.date,
                rf: t.first_dose(),
                diary: Some((
                    temp.diary.vaccination,
                    temp.day,
                    temp.diary.reference.to_string(),
                )),
            });
        }
    }
    let val = |r: &Row<'_>| trialforge_core::math::round_dp(r.value, r.dp);
    let mut cols = vec![
        ("STUDYID", constant(&rows, &ctx.design.studyid)),
        ("DOMAIN", constant(&rows, "VS")),
        ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
        ("VSSEQ", ints(&rows, |r| Some(r.seq))),
        ("VSTESTCD", chars(&rows, |r| Some(r.testcd.clone()))),
        ("VSTEST", chars(&rows, |r| Some(v(&r.testcd).test.clone()))),
        (
            "VSPOS",
            chars(&rows, |r| {
                r.diary
                    .is_none()
                    .then(|| v(&r.testcd).position.clone())
                    .flatten()
            }),
        ),
        ("VSORRES", chars(&rows, |r| Some(fmt_dp(val(r), r.dp)))),
        (
            "VSORRESU",
            chars(&rows, |r| Some(v(&r.testcd).unit.clone())),
        ),
        ("VSSTRESC", chars(&rows, |r| Some(fmt_dp(val(r), r.dp)))),
        ("VSSTRESN", nums(&rows, |r| Some(val(r)))),
        (
            "VSSTRESU",
            chars(&rows, |r| Some(v(&r.testcd).unit.clone())),
        ),
        (
            "VSLOBXFL",
            chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
        ),
        (
            "VSBLFL",
            chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
        ),
        ("VISITNUM", ints(&rows, |r| r.visit.as_ref().map(|x| x.0))),
        (
            "VISIT",
            chars(&rows, |r| r.visit.as_ref().map(|x| x.1.clone())),
        ),
        (
            "VISITDY",
            ints(&rows, |r| r.visit.as_ref().and_then(|x| x.2)),
        ),
        (
            "EPOCH",
            chars(&rows, |r| Some(r.epoch.as_str().to_string())),
        ),
        ("VSDTC", chars(&rows, |r| iso(r.date))),
        ("VSDY", ints(&rows, |r| dy(r.date, r.rf))),
    ];
    if rows.iter().any(|r| r.diary.is_some()) {
        fn diary<'r>(r: &'r Row<'_>) -> Option<&'r (u8, usize, String)> {
            r.diary.as_ref()
        }
        cols.extend([
            (
                "VSCAT",
                chars(&rows, |r| {
                    diary(r).map(|_| crate::vaccine::REACTOGENICITY.to_string())
                }),
            ),
            (
                "VSSCAT",
                chars(&rows, |r| diary(r).map(|_| "SYSTEMIC".to_string())),
            ),
            (
                "VSEVAL",
                chars(&rows, |r| diary(r).map(|_| "STUDY SUBJECT".to_string())),
            ),
            (
                "VSTPT",
                chars(&rows, |r| diary(r).map(|d| format!("DAY {}", d.1))),
            ),
            ("VSTPTNUM", ints(&rows, |r| diary(r).map(|d| d.1 as i64))),
            (
                "VSTPTREF",
                chars(&rows, |r| diary(r).map(|d| format!("VACCINATION {}", d.0))),
            ),
            ("VSRFTDTC", chars(&rows, |r| diary(r).map(|d| d.2.clone()))),
        ]);
    }
    build("VS", cols)
}

fn ft(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let items = reference::adas_items();
    let pos = |cd: &str| {
        items
            .iter()
            .position(|v| v.testcd == cd)
            .unwrap_or(usize::MAX)
    };
    let rows = finding_rows(ctx, FindingDomain::Ft, pos);
    let item = |r: &FindingRow<'_>| &items[pos(&r.f.testcd)];
    let done = |r: &FindingRow<'_>| !r.f.value.is_nan();
    let result = |r: &FindingRow<'_>| done(r).then(|| fmt_dp(r.f.value, 0));
    build(
        "FT",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "FT")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("FTSEQ", ints(&rows, |r| Some(r.seq))),
            ("FTTESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            ("FTTEST", chars(&rows, |r| Some(item(r).test.clone()))),
            ("FTCAT", constant(&rows, "ADAS-COG")),
            ("FTORRES", chars(&rows, result)),
            ("FTSTRESC", chars(&rows, result)),
            ("FTSTRESN", nums(&rows, |r| done(r).then_some(r.f.value))),
            (
                "FTSTAT",
                chars(&rows, |r| (!done(r)).then(|| "NOT DONE".to_string())),
            ),
            (
                "FTREASND",
                chars(&rows, |r| {
                    (!done(r)).then(|| "UNABLE TO COMPLETE".to_string())
                }),
            ),
            (
                "FTLOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("FTDTC", chars(&rows, |r| iso(r.f.date))),
            ("FTDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}

/// `FOCID` (CDISC OEFOCUS) of an eye.
fn focus(laterality: &str) -> String {
    if laterality == "LEFT" { "OS" } else { "OD" }.to_string()
}

fn oe(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let tests = reference::oe_tests();
    let pos = |cd: &str| {
        tests
            .iter()
            .position(|v| v.testcd == cd)
            .unwrap_or(usize::MAX)
    };
    let rows = finding_rows(ctx, FindingDomain::Oe, pos);
    let test = |r: &FindingRow<'_>| &tests[pos(&r.f.testcd)];
    let result = |r: &FindingRow<'_>| Some(fmt_dp(r.f.value, 0));
    let unit = |r: &FindingRow<'_>| test(r).unit.map(str::to_string);
    let lat = |r: &FindingRow<'_>| r.f.laterality.map(str::to_string);
    build(
        "OE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "OE")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("FOCID", chars(&rows, |r| r.f.laterality.map(focus))),
            ("OESEQ", ints(&rows, |r| Some(r.seq))),
            ("OETESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            ("OETEST", chars(&rows, |r| Some(test(r).test.to_string()))),
            (
                "OETSTDTL",
                chars(&rows, |r| test(r).detail.map(str::to_string)),
            ),
            (
                "OECAT",
                chars(&rows, |r| Some(test(r).category.to_string())),
            ),
            ("OEORRES", chars(&rows, result)),
            ("OEORRESU", chars(&rows, unit)),
            ("OESTRESC", chars(&rows, result)),
            ("OESTRESN", nums(&rows, |r| Some(r.f.value))),
            ("OESTRESU", chars(&rows, unit)),
            (
                "OELOC",
                chars(&rows, |r| Some(test(r).location.to_string())),
            ),
            ("OELAT", chars(&rows, lat)),
            (
                "OEMETHOD",
                chars(&rows, |r| Some(test(r).method.to_string())),
            ),
            (
                "OELOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("OEDTC", chars(&rows, |r| iso(r.f.date))),
            ("OEDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}

/// RE: respiratory function (FEV1 by spirometry).
fn re(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let rows = finding_rows(ctx, FindingDomain::Re, |_| 0);
    let result = |r: &FindingRow<'_>| Some(fmt_dp(r.f.value, r.f.dp.unwrap_or(2)));
    build(
        "RE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "RE")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("RESEQ", ints(&rows, |r| Some(r.seq))),
            ("RETESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            (
                "RETEST",
                constant(&rows, "Forced Expiratory Volume in 1 Second"),
            ),
            ("REORRES", chars(&rows, result)),
            ("REORRESU", constant(&rows, "L")),
            ("RESTRESC", chars(&rows, result)),
            ("RESTRESN", nums(&rows, |r| Some(r.f.value))),
            ("RESTRESU", constant(&rows, "L")),
            ("REMETHOD", constant(&rows, "SPIROMETRY")),
            (
                "RELOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("REDTC", chars(&rows, |r| iso(r.f.date))),
            ("REDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}

/// SC: the study eye selected at screening (pharmaverse/admiralophtha convention:
/// `SCTESTCD = FOCID`, `SCSTRESC` in OEFOCUS).
fn sc(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let rows: Vec<(&Subject, &Timeline, &str)> = ctx
        .sorted()
        .filter_map(|(s, t)| Some((s, t, t.study_eye?)))
        .collect();
    build(
        "SC",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "SC")),
            ("USUBJID", chars(&rows, |(s, _, _)| Some(s.usubjid.clone()))),
            ("SCSEQ", ints(&rows, |_| Some(1))),
            ("SCTESTCD", constant(&rows, "FOCID")),
            (
                "SCTEST",
                constant(&rows, "Focus of Study-Specific Interest"),
            ),
            ("SCCAT", constant(&rows, "STUDY EYE SELECTION")),
            (
                "SCORRES",
                chars(&rows, |(_, _, e)| {
                    Some(
                        if *e == "LEFT" {
                            "Left Eye"
                        } else {
                            "Right Eye"
                        }
                        .to_string(),
                    )
                }),
            ),
            ("SCSTRESC", chars(&rows, |(_, _, e)| Some(focus(e)))),
            ("EPOCH", constant(&rows, Epoch::Screening.as_str())),
            ("SCDTC", chars(&rows, |(s, _, _)| iso(s.screening))),
            (
                "SCDY",
                ints(&rows, |(s, t, _)| dy(s.screening, t.first_dose())),
            ),
        ],
    )
}

fn ce(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        usubjid: &'a str,
        seq: i64,
        e: &'a trialforge_processes::timeline::ClinicalEvent,
        epoch: Epoch,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, e) in t.clinical_events.iter().enumerate() {
            rows.push(Row {
                usubjid: &s.usubjid,
                seq: k as i64 + 1,
                e,
                epoch: t.epoch_of(e.date),
                rf: t.first_dose(),
            });
        }
    }
    let soc = |r: &Row<'_>| {
        dictionaries::ae_term(&r.e.pt)
            .expect("clinical events come from the dictionary")
            .soc
            .clone()
    };
    let mut cols = vec![
        ("STUDYID", constant(&rows, &ctx.design.studyid)),
        ("DOMAIN", constant(&rows, "CE")),
        ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
        ("CESEQ", ints(&rows, |r| Some(r.seq))),
        ("CETERM", chars(&rows, |r| Some(r.e.term.clone()))),
        ("CEDECOD", chars(&rows, |r| Some(r.e.pt.clone()))),
        ("CECAT", chars(&rows, |r| Some(r.e.category.to_string()))),
        ("CEBODSYS", chars(&rows, |r| Some(soc(r)))),
        (
            "EPOCH",
            chars(&rows, |r| Some(r.epoch.as_str().to_string())),
        ),
        ("CESTDTC", chars(&rows, |r| iso(r.e.date))),
        ("CESTDY", ints(&rows, |r| dy(r.e.date, r.rf))),
    ];
    // Graded events with a resolution date (exacerbations) carry severity and end.
    if rows.iter().any(|r| r.e.severity.is_some()) {
        cols.push((
            "CESEV",
            chars(&rows, |r| r.e.severity.map(|s| s.as_str().to_string())),
        ));
        cols.push(("CEENDTC", chars(&rows, |r| r.e.end.and_then(iso))));
        cols.push(("CEENDY", ints(&rows, |r| r.e.end.and_then(|d| dy(d, r.rf)))));
    }
    build("CE", cols)
}
