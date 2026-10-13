//! Observation layer: turns subject timelines into SDTM records, including SUPP--, RELREC and
//! the trial design model datasets.
//!
//! The spike emits DM, EX and DS from the oncology course. Records are sorted by the SDTM
//! natural keys (`USUBJID`, then `--SEQ`), and `--DY` is relative to `RFSTDTC` (first dose).

use trialforge_cdisc::{DM, DS, EX, MetadataError};
use trialforge_core::{Dataset, Date, Values, study_day};
use trialforge_model::{Design, Subject};
use trialforge_processes::oncology::{Course, TreatmentStop};

/// The SDTM datasets of the oncology spike.
#[derive(Clone, Debug, PartialEq)]
pub struct Sdtm {
    /// Demographics.
    pub dm: Dataset,
    /// Exposure.
    pub ex: Dataset,
    /// Disposition.
    pub ds: Dataset,
}

/// Emits DM, EX and DS. `subjects` and `courses` are parallel.
pub fn emit(
    design: &Design,
    subjects: &[Subject],
    courses: &[Course],
) -> Result<Sdtm, MetadataError> {
    assert_eq!(subjects.len(), courses.len(), "one course per subject");
    let mut order: Vec<usize> = (0..subjects.len()).collect();
    order.sort_by(|&a, &b| subjects[a].usubjid.cmp(&subjects[b].usubjid));
    Ok(Sdtm {
        dm: dm(design, subjects, courses, &order)?,
        ex: ex(design, subjects, courses, &order)?,
        ds: ds(design, subjects, courses, &order)?,
    })
}

fn chars<I: IntoIterator<Item = S>, S: Into<String>>(it: I) -> Values {
    Values::Char(it.into_iter().map(|s| Some(s.into())).collect())
}

fn opt_chars<I: IntoIterator<Item = Option<String>>>(it: I) -> Values {
    Values::Char(it.into_iter().collect())
}

fn iso(d: Date) -> String {
    d.to_string()
}

fn dm(
    design: &Design,
    subjects: &[Subject],
    courses: &[Course],
    order: &[usize],
) -> Result<Dataset, MetadataError> {
    let pop = &design.population;
    let rows: Vec<(&Subject, &Course)> =
        order.iter().map(|&i| (&subjects[i], &courses[i])).collect();
    let n = rows.len();
    let arm = |s: &Subject| &design.arms[s.arm];
    DM.build(vec![
        (
            "STUDYID",
            chars(std::iter::repeat_n(design.studyid.as_str(), n)),
        ),
        ("DOMAIN", chars(std::iter::repeat_n("DM", n))),
        (
            "USUBJID",
            chars(rows.iter().map(|(s, _)| s.usubjid.as_str())),
        ),
        ("SUBJID", chars(rows.iter().map(|(s, _)| s.subjid.as_str()))),
        (
            "RFSTDTC",
            chars(rows.iter().map(|(s, _)| iso(s.treatment_start))),
        ),
        (
            "RFENDTC",
            chars(rows.iter().map(|(_, c)| iso(c.last_contact))),
        ),
        (
            "RFXSTDTC",
            chars(rows.iter().map(|(s, _)| iso(s.treatment_start))),
        ),
        (
            "RFXENDTC",
            chars(
                rows.iter()
                    .map(|(_, c)| iso(*c.doses.last().expect("at least one dose"))),
            ),
        ),
        ("RFICDTC", chars(rows.iter().map(|(s, _)| iso(s.consent)))),
        (
            "RFPENDTC",
            chars(rows.iter().map(|(_, c)| iso(c.last_contact))),
        ),
        (
            "DTHDTC",
            opt_chars(rows.iter().map(|(_, c)| c.death.map(iso))),
        ),
        (
            "DTHFL",
            opt_chars(rows.iter().map(|(_, c)| c.death.map(|_| "Y".to_string()))),
        ),
        (
            "SITEID",
            chars(
                rows.iter()
                    .map(|(s, _)| design.sites[s.site].siteid.as_str()),
            ),
        ),
        (
            "AGE",
            Values::Int(rows.iter().map(|(s, _)| Some(s.age)).collect()),
        ),
        ("AGEU", chars(std::iter::repeat_n("YEARS", n))),
        (
            "SEX",
            chars(rows.iter().map(|(s, _)| pop.sex.levels[s.sex].0.as_str())),
        ),
        (
            "RACE",
            chars(rows.iter().map(|(s, _)| pop.race.levels[s.race].0.as_str())),
        ),
        (
            "ETHNIC",
            chars(
                rows.iter()
                    .map(|(s, _)| pop.ethnic.levels[s.ethnic].0.as_str()),
            ),
        ),
        (
            "ARMCD",
            chars(rows.iter().map(|(s, _)| arm(s).armcd.as_str())),
        ),
        ("ARM", chars(rows.iter().map(|(s, _)| arm(s).arm.as_str()))),
        (
            "ACTARMCD",
            chars(rows.iter().map(|(s, _)| arm(s).armcd.as_str())),
        ),
        (
            "ACTARM",
            chars(rows.iter().map(|(s, _)| arm(s).arm.as_str())),
        ),
        (
            "COUNTRY",
            chars(
                rows.iter()
                    .map(|(s, _)| design.sites[s.site].country.as_str()),
            ),
        ),
        ("DMDTC", chars(rows.iter().map(|(s, _)| iso(s.consent)))),
        (
            "DMDY",
            Values::Int(
                rows.iter()
                    .map(|(s, _)| Some(i64::from(study_day(s.consent, s.treatment_start))))
                    .collect(),
            ),
        ),
    ])
}

fn ex(
    design: &Design,
    subjects: &[Subject],
    courses: &[Course],
    order: &[usize],
) -> Result<Dataset, MetadataError> {
    let (mut usubjid, mut seq, mut extrt, mut dose, mut unit, mut form, mut freq, mut route) = (
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    let (mut dtc, mut dy) = (Vec::new(), Vec::new());
    for &i in order {
        let (s, c) = (&subjects[i], &courses[i]);
        let t = &design.arms[s.arm].treatment;
        for (k, &d) in c.doses.iter().enumerate() {
            usubjid.push(Some(s.usubjid.clone()));
            seq.push(Some(k as i64 + 1));
            extrt.push(Some(t.name.clone()));
            dose.push(Some(t.dose));
            unit.push(Some(t.unit.clone()));
            form.push(Some(t.form.clone()));
            freq.push(Some(t.frequency.clone()));
            route.push(Some(t.route.clone()));
            dtc.push(Some(iso(d)));
            dy.push(Some(i64::from(study_day(d, s.treatment_start))));
        }
    }
    let n = usubjid.len();
    EX.build(vec![
        (
            "STUDYID",
            chars(std::iter::repeat_n(design.studyid.as_str(), n)),
        ),
        ("DOMAIN", chars(std::iter::repeat_n("EX", n))),
        ("USUBJID", Values::Char(usubjid)),
        ("EXSEQ", Values::Int(seq)),
        ("EXTRT", Values::Char(extrt)),
        ("EXDOSE", Values::Num(dose)),
        ("EXDOSU", Values::Char(unit)),
        ("EXDOSFRM", Values::Char(form)),
        ("EXDOSFRQ", Values::Char(freq)),
        ("EXROUTE", Values::Char(route)),
        ("EPOCH", chars(std::iter::repeat_n("TREATMENT", n))),
        ("EXSTDTC", Values::Char(dtc.clone())),
        ("EXENDTC", Values::Char(dtc)),
        ("EXSTDY", Values::Int(dy.clone())),
        ("EXENDY", Values::Int(dy)),
    ])
}

/// One DS record before sequencing.
struct DsRecord {
    term: String,
    category: &'static str,
    epoch: &'static str,
    date: Date,
}

fn ds_records(s: &Subject, c: &Course) -> Vec<DsRecord> {
    let milestone = |term: &str, date| DsRecord {
        term: term.into(),
        category: "PROTOCOL MILESTONE",
        epoch: "SCREENING",
        date,
    };
    let mut recs = vec![
        milestone("INFORMED CONSENT OBTAINED", s.consent),
        milestone("RANDOMIZED", s.randomization),
    ];
    let treatment_term = match &c.treatment_stop_reason {
        TreatmentStop::Ongoing => None,
        TreatmentStop::Progression => Some("PROGRESSIVE DISEASE".to_string()),
        TreatmentStop::Death => Some("DEATH".to_string()),
        TreatmentStop::Withdrawal(reason) => Some(reason.clone()),
        TreatmentStop::AdverseEvent => Some("ADVERSE EVENT".to_string()),
    };
    if let Some(term) = treatment_term {
        recs.push(DsRecord {
            term,
            category: "DISPOSITION EVENT",
            epoch: "TREATMENT",
            date: c.treatment_stop,
        });
    }
    let study_end = match (c.death, c.withdrawal) {
        (Some(d), _) => Some(("DEATH".to_string(), d)),
        (None, Some(w)) => c.withdrawal_reason.clone().map(|r| (r, w)),
        (None, None) => None,
    };
    if let Some((term, date)) = study_end {
        recs.push(DsRecord {
            term,
            category: "DISPOSITION EVENT",
            epoch: "FOLLOW-UP",
            date,
        });
    }
    recs
}

fn ds(
    design: &Design,
    subjects: &[Subject],
    courses: &[Course],
    order: &[usize],
) -> Result<Dataset, MetadataError> {
    let (mut usubjid, mut seq, mut term, mut cat, mut epoch, mut dtc, mut dy) = (
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
        Vec::new(),
    );
    for &i in order {
        let s = &subjects[i];
        for (k, r) in ds_records(s, &courses[i]).into_iter().enumerate() {
            usubjid.push(Some(s.usubjid.clone()));
            seq.push(Some(k as i64 + 1));
            term.push(Some(r.term));
            cat.push(Some(r.category.to_string()));
            epoch.push(Some(r.epoch.to_string()));
            dtc.push(Some(iso(r.date)));
            dy.push(Some(i64::from(study_day(r.date, s.treatment_start))));
        }
    }
    let n = usubjid.len();
    DS.build(vec![
        (
            "STUDYID",
            chars(std::iter::repeat_n(design.studyid.as_str(), n)),
        ),
        ("DOMAIN", chars(std::iter::repeat_n("DS", n))),
        ("USUBJID", Values::Char(usubjid)),
        ("DSSEQ", Values::Int(seq)),
        ("DSTERM", Values::Char(term.clone())),
        ("DSDECOD", Values::Char(term)),
        ("DSCAT", Values::Char(cat)),
        ("EPOCH", Values::Char(epoch)),
        ("DSSTDTC", Values::Char(dtc)),
        ("DSSTDY", Values::Int(dy)),
    ])
}
