//! Vaccine domains: IS (antibody titres), FACE (the reactogenicity diary) and the
//! reactogenicity records of CE (one per solicited reaction per vaccination) and VS (daily
//! temperatures), following pharmaverse's vaccine SDTM.

use trialforge_cdisc::{MetadataError, build, dictionaries};
use trialforge_core::{Dataset, Date};
use trialforge_model::Subject;
use trialforge_processes::timeline::{Diary, FindingDomain, Reaction, Severity, Timeline};

use crate::subject::finding_rows;
use crate::{Context, chars, constant, dy, fmt_dp, ints, nums};

/// `--CAT` of reactogenicity records.
pub const REACTOGENICITY: &str = "REACTOGENICITY";

fn iso(d: Date) -> Option<String> {
    Some(d.to_string())
}

fn scat(site: bool) -> &'static str {
    if site {
        "ADMINISTRATION SITE"
    } else {
        "SYSTEMIC"
    }
}

/// The CE severity of a measured reaction: by its largest diameter (cm).
fn diameter_severity(cm: f64) -> Severity {
    if cm > 10.0 {
        Severity::Severe
    } else if cm > 5.0 {
        Severity::Moderate
    } else {
        Severity::Mild
    }
}

/// IS: neutralising antibody titres.
pub(crate) fn is(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let rows = finding_rows(ctx, FindingDomain::Is, |_| 0);
    let (lloq, uloq) = (10.0, 20480.0);
    let result = |r: &crate::subject::FindingRow<'_>| {
        Some(if r.f.value <= 0.0 {
            format!("<{lloq}")
        } else if r.f.value.is_infinite() {
            format!(">{uloq}")
        } else {
            fmt_dp(r.f.value, 0)
        })
    };
    let quantified = |r: &crate::subject::FindingRow<'_>| r.f.value > 0.0 && r.f.value.is_finite();
    build(
        "IS",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "IS")),
            ("USUBJID", chars(&rows, |r| Some(r.usubjid.to_string()))),
            ("ISSEQ", ints(&rows, |r| Some(r.seq))),
            ("ISTESTCD", chars(&rows, |r| Some(r.f.testcd.clone()))),
            (
                "ISTEST",
                constant(&rows, "Neutralizing Microbial-induced Antibody"),
            ),
            ("ISCAT", constant(&rows, "SEROLOGY")),
            ("ISORRES", chars(&rows, result)),
            ("ISORRESU", constant(&rows, "titer")),
            ("ISSTRESC", chars(&rows, result)),
            (
                "ISSTRESN",
                nums(&rows, |r| quantified(r).then_some(r.f.value)),
            ),
            ("ISSTRESU", constant(&rows, "titer")),
            ("ISBDAGNT", constant(&rows, "VIRUS X")),
            ("ISSPEC", constant(&rows, "SERUM")),
            ("ISMETHOD", constant(&rows, "MICRONEUTRALIZATION ASSAY")),
            (
                "ISBLFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            (
                "ISLOBXFL",
                chars(&rows, |r| r.baseline.then(|| "Y".to_string())),
            ),
            ("ISLLOQ", nums(&rows, |_| Some(lloq))),
            ("ISULOQ", nums(&rows, |_| Some(uloq))),
            ("VISITNUM", ints(&rows, |r| Some(r.f.visitnum))),
            ("VISIT", chars(&rows, |r| Some(r.visit.clone()))),
            ("VISITDY", ints(&rows, |r| r.visitdy)),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.epoch.as_str().to_string())),
            ),
            ("ISDTC", chars(&rows, |r| iso(r.f.date))),
            ("ISDY", ints(&rows, |r| dy(r.f.date, r.rf))),
        ],
    )
}

/// One CE row: a solicited reaction after one vaccination.
struct CeRow<'a> {
    s: &'a Subject,
    t: &'a Timeline,
    seq: i64,
    diary: &'a Diary,
    r: Reaction,
}

/// Fever (38.0 °C or more) as a reaction, graded by the maximum temperature (FDA toxicity
/// grading: 38.0–38.4 mild, 38.5–38.9 moderate, 39.0 or more severe).
fn fever(diary: &Diary) -> Reaction {
    let days = diary
        .temperatures
        .iter()
        .map(|&(date, c)| trialforge_processes::timeline::DiaryDay {
            date,
            occurred: c >= 38.0,
            severity: (c >= 38.0).then_some(if c >= 39.0 {
                Severity::Severe
            } else if c >= 38.5 {
                Severity::Moderate
            } else {
                Severity::Mild
            }),
            diameter: None,
        })
        .collect();
    Reaction {
        term: "FEVER",
        pt: "Pyrexia",
        site: false,
        measured: false,
        days,
    }
}

impl CeRow<'_> {
    fn on(&self) -> Vec<&trialforge_processes::timeline::DiaryDay> {
        self.r.days.iter().filter(|d| d.occurred).collect()
    }
    fn occurred(&self) -> bool {
        !self.on().is_empty()
    }
    fn severity(&self) -> Option<Severity> {
        let on = self.on();
        if self.r.measured {
            on.iter()
                .filter_map(|d| d.diameter)
                .reduce(f64::max)
                .map(diameter_severity)
        } else {
            on.iter().filter_map(|d| d.severity).max()
        }
    }
}

/// CE for a vaccine study: every solicited reaction (fever from the temperatures) after every
/// vaccination, with whether it occurred and, if so, its dates and maximum severity.
pub(crate) fn ce(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let mut rows: Vec<CeRow<'_>> = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut seq = 0;
        for diary in &t.diaries {
            for r in diary.reactions.iter().cloned().chain([fever(diary)]) {
                seq += 1;
                rows.push(CeRow {
                    s,
                    t,
                    seq,
                    diary,
                    r,
                });
            }
        }
    }
    let soc = |r: &CeRow<'_>| {
        dictionaries::ae_term(r.r.pt)
            .expect("solicited reactions come from the dictionary")
            .soc
            .clone()
    };
    let first = |r: &CeRow<'_>| r.on().first().map(|d| d.date);
    let last = |r: &CeRow<'_>| r.on().last().map(|d| d.date);
    let reference = |r: &CeRow<'_>| r.diary.reference.date();
    let site = |r: &CeRow<'_>| r.r.site.then_some(r.diary.location);
    build(
        "CE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "CE")),
            ("USUBJID", chars(&rows, |r| Some(r.s.usubjid.clone()))),
            ("CESEQ", ints(&rows, |r| Some(r.seq))),
            (
                "CELNKGRP",
                chars(&rows, |r| {
                    Some(format!("VACCINATION {}-{}", r.diary.vaccination, r.r.term))
                }),
            ),
            ("CETERM", chars(&rows, |r| Some(r.r.term.to_string()))),
            ("CEDECOD", chars(&rows, |r| Some(r.r.pt.to_string()))),
            ("CECAT", constant(&rows, REACTOGENICITY)),
            ("CESCAT", chars(&rows, |r| Some(scat(r.r.site).to_string()))),
            ("CEPRESP", constant(&rows, "Y")),
            (
                "CEOCCUR",
                chars(&rows, |r| {
                    Some(if r.occurred() { "Y" } else { "N" }.to_string())
                }),
            ),
            ("CEBODSYS", chars(&rows, |r| Some(soc(r)))),
            ("CELOC", chars(&rows, |r| site(r).map(|l| l.0.to_string()))),
            ("CELAT", chars(&rows, |r| site(r).map(|l| l.1.to_string()))),
            (
                "CESEV",
                chars(&rows, |r| r.severity().map(|s| s.as_str().to_string())),
            ),
            (
                "CEREL",
                chars(&rows, |r| r.occurred().then(|| "RELATED".to_string())),
            ),
            (
                "CEOUT",
                chars(&rows, |r| {
                    r.occurred().then(|| "RECOVERED/RESOLVED".to_string())
                }),
            ),
            (
                "EPOCH",
                chars(&rows, |r| {
                    Some(
                        r.t.epoch_of(first(r).unwrap_or(reference(r)))
                            .as_str()
                            .to_string(),
                    )
                }),
            ),
            ("CESTDTC", chars(&rows, |r| first(r).and_then(iso))),
            ("CEENDTC", chars(&rows, |r| last(r).and_then(iso))),
            (
                "CESTDY",
                ints(&rows, |r| first(r).and_then(|d| dy(d, r.t.first_dose()))),
            ),
            (
                "CEENDY",
                ints(&rows, |r| last(r).and_then(|d| dy(d, r.t.first_dose()))),
            ),
            ("CETPT", constant(&rows, "DAY 7")),
            ("CETPTNUM", ints(&rows, |_| Some(7))),
            (
                "CETPTREF",
                chars(&rows, |r| {
                    Some(format!("VACCINATION {}", r.diary.vaccination))
                }),
            ),
            (
                "CERFTDTC",
                chars(&rows, |r| Some(r.diary.reference.to_string())),
            ),
            ("CEEVINTX", constant(&rows, "SINCE VACCINATION")),
        ],
    )
}

/// One FACE row: one finding about a solicited reaction on one diary day.
struct FaRow<'a> {
    s: &'a Subject,
    t: &'a Timeline,
    seq: i64,
    diary: &'a Diary,
    r: &'a Reaction,
    day: usize,
    test: (&'static str, &'static str),
    result: String,
    numeric: Option<f64>,
    unit: Option<&'static str>,
}

/// FACE: the diary, one record per reaction, test and day.
pub(crate) fn face(ctx: &Context<'_>) -> Result<Dataset, MetadataError> {
    let mut rows: Vec<FaRow<'_>> = Vec::new();
    for (s, t) in ctx.sorted() {
        let mut seq = 0;
        for diary in &t.diaries {
            for r in &diary.reactions {
                let tests: [(&str, &str); 3] = [
                    ("OCCUR", "Occurrence Indicator"),
                    ("SEV", "Severity/Intensity"),
                    ("DIAMETER", "Diameter"),
                ];
                for test in tests {
                    for (k, d) in r.days.iter().enumerate() {
                        let (result, numeric, unit) = match test.0 {
                            "OCCUR" => (if d.occurred { "Y" } else { "N" }.to_string(), None, None),
                            "SEV" => match d.severity {
                                Some(sev) => (sev.as_str().to_string(), None, None),
                                None => continue,
                            },
                            _ => match d.diameter {
                                Some(cm) => (fmt_dp(cm, 1), Some(cm), Some("cm")),
                                None => continue,
                            },
                        };
                        seq += 1;
                        rows.push(FaRow {
                            s,
                            t,
                            seq,
                            diary,
                            r,
                            day: k + 1,
                            test,
                            result,
                            numeric,
                            unit,
                        });
                    }
                }
            }
        }
    }
    let date = |r: &FaRow<'_>| r.r.days[r.day - 1].date;
    let site = |r: &FaRow<'_>| r.r.site.then_some(r.diary.location);
    let group = |r: &FaRow<'_>| format!("VACCINATION {}-{}", r.diary.vaccination, r.r.term);
    build(
        "FACE",
        vec![
            ("STUDYID", constant(&rows, &ctx.design.studyid)),
            ("DOMAIN", constant(&rows, "FA")),
            ("USUBJID", chars(&rows, |r| Some(r.s.usubjid.clone()))),
            ("FASEQ", ints(&rows, |r| Some(r.seq))),
            ("FALNKGRP", chars(&rows, |r| Some(group(r)))),
            (
                "FALNKID",
                chars(&rows, |r| Some(format!("{}-DAY {}", group(r), r.day))),
            ),
            ("FATESTCD", chars(&rows, |r| Some(r.test.0.to_string()))),
            ("FATEST", chars(&rows, |r| Some(r.test.1.to_string()))),
            ("FAOBJ", chars(&rows, |r| Some(r.r.term.to_string()))),
            ("FACAT", constant(&rows, REACTOGENICITY)),
            ("FASCAT", chars(&rows, |r| Some(scat(r.r.site).to_string()))),
            ("FAORRES", chars(&rows, |r| Some(r.result.clone()))),
            ("FAORRESU", chars(&rows, |r| r.unit.map(str::to_string))),
            ("FASTRESC", chars(&rows, |r| Some(r.result.clone()))),
            ("FASTRESN", nums(&rows, |r| r.numeric)),
            ("FASTRESU", chars(&rows, |r| r.unit.map(str::to_string))),
            ("FALOC", chars(&rows, |r| site(r).map(|l| l.0.to_string()))),
            ("FALAT", chars(&rows, |r| site(r).map(|l| l.1.to_string()))),
            ("FAEVAL", constant(&rows, "STUDY SUBJECT")),
            (
                "EPOCH",
                chars(&rows, |r| Some(r.t.epoch_of(date(r)).as_str().to_string())),
            ),
            ("FADTC", chars(&rows, |r| iso(date(r)))),
            ("FADY", ints(&rows, |r| dy(date(r), r.t.first_dose()))),
            ("FATPT", chars(&rows, |r| Some(format!("DAY {}", r.day)))),
            ("FATPTNUM", ints(&rows, |r| Some(r.day as i64))),
            (
                "FATPTREF",
                chars(&rows, |r| {
                    Some(format!("VACCINATION {}", r.diary.vaccination))
                }),
            ),
            (
                "FARFTDTC",
                chars(&rows, |r| Some(r.diary.reference.to_string())),
            ),
        ],
    )
}

/// A reactogenicity temperature for VS: (subject, timeline, vaccination, day, date, °C).
pub(crate) struct Temperature<'a> {
    pub t: &'a Timeline,
    pub diary: &'a Diary,
    pub day: usize,
    pub date: Date,
    pub value: f64,
}

/// The daily diary temperatures of a subject, in vaccination and day order.
pub(crate) fn temperatures(t: &Timeline) -> Vec<Temperature<'_>> {
    t.diaries
        .iter()
        .flat_map(|diary| {
            diary
                .temperatures
                .iter()
                .enumerate()
                .map(move |(k, &(date, value))| Temperature {
                    t,
                    diary,
                    day: k + 1,
                    date,
                    value,
                })
        })
        .collect()
}
