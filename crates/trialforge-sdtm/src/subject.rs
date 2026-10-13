//! Subject-level domains: DM, SE, SV, EX, CM, AE, SUPPAE, DS, MH, LB, VS.

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
    build(
        "DM",
        vec![
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
        ],
    )
}

fn se(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row {
        usubjid: String,
        seq: i64,
        etcd: String,
        element: String,
        epoch: Epoch,
        start: Date,
        end: Date,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, e) in t.elements.iter().enumerate() {
            let (etcd, element) = element(ctx, e.epoch, s.arm);
            rows.push(Row {
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
            ("TAETORD", ints(&rows, |r| Some(taetord(r.epoch)))),
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
            ("SVENDTC", chars(&rows, |r| iso(r.date))),
            ("SVSTDY", ints(&rows, |r| dy(r.date, r.rf))),
            ("SVENDY", ints(&rows, |r| dy(r.date, r.rf))),
        ],
    )
}

fn ex(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    struct Row<'a> {
        s: &'a Subject,
        seq: i64,
        start: Date,
        end: Date,
        dose: f64,
        epoch: Epoch,
        rf: Date,
    }
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        for (k, e) in t.exposure.iter().enumerate() {
            rows.push(Row {
                s,
                seq: k as i64 + 1,
                start: e.start,
                end: e.end,
                dose: e.dose,
                epoch: t.epoch_of(e.start),
                rf: t.first_dose(),
            });
        }
    }
    let trt = |r: &Row<'_>| &ctx.design.arms[r.s.arm].treatment;
    build(
        "EX",
        vec![
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
            ("EXSTDTC", chars(&rows, |r| iso(r.start))),
            ("EXENDTC", chars(&rows, |r| iso(r.end))),
            ("EXSTDY", ints(&rows, |r| dy(r.start, r.rf))),
            ("EXENDY", ints(&rows, |r| dy(r.end, r.rf))),
        ],
    )
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
    let suppae = build(
        "SUPPAE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("RDOMAIN", constant(&rows, "AE")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("IDVAR", constant(&rows, "AESEQ")),
            ("IDVARVAL", chars(&rows, |r| Some(r.seq.to_string()))),
            ("QNAM", constant(&rows, "AETRTEM")),
            ("QLABEL", constant(&rows, "Treatment Emergent Flag")),
            (
                "QVAL",
                chars(&rows, |r| yn(treatment_emergent(r.t, r.a.start))),
            ),
            ("QORIG", constant(&rows, "Derived")),
            ("QEVAL", chars(&rows, |_| None)),
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
struct FindingRow<'a> {
    usubjid: &'a str,
    seq: i64,
    f: &'a Finding,
    visit: String,
    visitdy: Option<i64>,
    epoch: Epoch,
    rf: Date,
    baseline: bool,
    female: bool,
}

fn finding_rows<'a>(
    ctx: &'a Context<'_>,
    domain: FindingDomain,
    order: impl Fn(&str) -> usize,
) -> Vec<FindingRow<'a>> {
    let schedule = ctx.design.schedule.as_ref();
    let mut rows = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut fs: Vec<&Finding> = t.findings.iter().filter(|f| f.domain == domain).collect();
        fs.sort_by(|a, b| (order(&a.testcd), a.visitnum).cmp(&(order(&b.testcd), b.visitnum)));
        // Last observation on or before first dose, per test.
        let baseline: Vec<Option<usize>> = fs
            .iter()
            .enumerate()
            .map(|(i, f)| {
                let last = fs
                    .iter()
                    .enumerate()
                    .filter(|(_, g)| g.testcd == f.testcd && g.date <= t.first_dose())
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
    let vitals = reference::vital_signs();
    let pos = |cd: &str| {
        vitals
            .iter()
            .position(|v| v.testcd == cd)
            .unwrap_or(usize::MAX)
    };
    let rows = finding_rows(ctx, FindingDomain::Vs, pos);
    let v = |r: &FindingRow<'_>| &vitals[pos(&r.f.testcd)];
    let val = |r: &FindingRow<'_>| trialforge_core::math::round_dp(r.f.value, v(r).dp);
    build(
        "VS",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "VS")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("VSSEQ", ints(&rows, |r| Some(r.seq))),
            ("VSTESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            ("VSTEST", chars(&rows, |r| Some(v(r).test.clone()))),
            ("VSPOS", chars(&rows, |r| v(r).position.clone())),
            ("VSORRES", chars(&rows, |r| Some(fmt_dp(val(r), v(r).dp)))),
            ("VSORRESU", chars(&rows, |r| Some(v(r).unit.clone()))),
            ("VSSTRESC", chars(&rows, |r| Some(fmt_dp(val(r), v(r).dp)))),
            ("VSSTRESN", nums(&rows, |r| Some(val(r)))),
            ("VSSTRESU", chars(&rows, |r| Some(v(r).unit.clone()))),
            (
                "VSLOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            (
                "VSBLFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("VSDTC", chars(&rows, |r| iso(r.f.date))),
            ("VSDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}
