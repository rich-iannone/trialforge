//! A subject's observed course, independent of any SDTM layout.
//!
//! Every process engine produces a [`Timeline`] per subject; the SDTM emitters read only
//! timelines (plus the design), so all archetypes share one observation layer.

use trialforge_core::{Date, DateTime};
use trialforge_model::Epoch;

/// A constant-dose interval of study treatment (one EX record).
#[derive(Clone, Debug, PartialEq)]
pub struct Exposure {
    /// First day of the interval.
    pub start: Date,
    /// Last day of the interval.
    pub end: Date,
    /// Dose per administration.
    pub dose: f64,
    /// Time of a single administration (`start == end`), when the study records dosing times.
    pub time: Option<DateTime>,
    /// The crossover period (1-based) the administration belongs to.
    pub period: Option<u8>,
    /// Site of administration (`EXLOC`, `EXLAT`), e.g. an injection into one eye or arm.
    pub location: Option<(&'static str, &'static str)>,
    /// Link group (`EXLNKGRP`), e.g. `VACCINATION 1`.
    pub group: Option<String>,
    /// The visit of the administration (`VISITNUM`), when recorded.
    pub visitnum: Option<i64>,
    /// Whether this is the arm's switch treatment (crossover on progression) rather than its
    /// own.
    pub switched: bool,
}

/// A span of one trial element (one SE record).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ElementSpan {
    /// The epoch the element belongs to.
    pub epoch: Epoch,
    /// Start date.
    pub start: Date,
    /// End date (the start of the next element, or the end of participation).
    pub end: Date,
}

/// An attended visit (one SV record).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct VisitOccurrence {
    /// `VISITNUM` of the planned visit.
    pub visitnum: i64,
    /// Visit date.
    pub date: Date,
    /// Last date of a visit that spans days (e.g. an in-house PK profile); `None` = `date`.
    pub end: Option<Date>,
}

/// Which disposition a DS record describes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DispositionKind {
    /// A protocol milestone (consent, randomisation).
    Milestone,
    /// End of study treatment.
    EndOfTreatment,
    /// End of study participation.
    EndOfStudy,
}

/// A disposition event (one DS record).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DispositionEvent {
    /// `DSDECOD`.
    pub decod: String,
    /// What it describes.
    pub kind: DispositionKind,
    /// Date.
    pub date: Date,
}

/// AE severity.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum Severity {
    /// `MILD`.
    Mild,
    /// `MODERATE`.
    Moderate,
    /// `SEVERE`.
    Severe,
}

impl Severity {
    /// The CT submission value.
    pub fn as_str(self) -> &'static str {
        match self {
            Severity::Mild => "MILD",
            Severity::Moderate => "MODERATE",
            Severity::Severe => "SEVERE",
        }
    }
}

/// Why an AE is serious.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct SeriousCriteria {
    /// `AESDTH`.
    pub death: bool,
    /// `AESLIFE`.
    pub life_threatening: bool,
    /// `AESHOSP`.
    pub hospitalization: bool,
    /// `AESDISAB`.
    pub disability: bool,
    /// `AESMIE`.
    pub medically_important: bool,
}

/// An adverse event (one AE record).
#[derive(Clone, Debug, PartialEq)]
pub struct AdverseEvent {
    /// Preferred Term (coded); also used as the verbatim term.
    pub pt: String,
    /// Start date.
    pub start: Date,
    /// End date; `None` if ongoing at the end of participation.
    pub end: Option<Date>,
    /// Severity.
    pub severity: Severity,
    /// Seriousness criteria; `None` if not serious.
    pub serious: Option<SeriousCriteria>,
    /// `AEREL`.
    pub causality: &'static str,
    /// `AEACN`.
    pub action: &'static str,
    /// `AEOUT`.
    pub outcome: &'static str,
    /// Whether a concomitant medication was started for it.
    pub treated: bool,
    /// In a dose-escalation study, whether the AE is a dose-limiting toxicity (`None` outside
    /// the DLT observation window, or in other studies).
    pub dlt: Option<bool>,
}

/// A medication start or end that may be a partial ISO 8601 date (e.g. `2019-04`).
pub type PartialDate = String;

/// A concomitant or prior medication (one CM record).
#[derive(Clone, Debug, PartialEq)]
pub struct Medication {
    /// RxNorm ingredient (`CMDECOD`; also the verbatim `CMTRT`).
    pub ingredient: String,
    /// `CMINDC`.
    pub indication: String,
    /// Start (complete or partial ISO 8601).
    pub start: PartialDate,
    /// End; `None` if ongoing.
    pub end: Option<Date>,
}

/// A medical-history condition (one MH record).
#[derive(Clone, Debug, PartialEq)]
pub struct HistoryItem {
    /// Preferred Term.
    pub pt: String,
    /// `MHCAT`.
    pub category: String,
    /// Onset (partial ISO 8601, e.g. `2015` or `2015-06`).
    pub start: PartialDate,
    /// Whether ongoing at screening.
    pub ongoing: bool,
}

/// The findings domain a measurement belongs to.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FindingDomain {
    /// Laboratory.
    Lb,
    /// Vital signs.
    Vs,
    /// Functional tests (e.g. ADAS-Cog items); a NaN value is an item not done.
    Ft,
    /// Ophthalmic examinations (e.g. visual acuity), one per eye.
    Oe,
    /// Respiratory system findings (e.g. FEV1).
    Re,
    /// Immunogenicity specimen assessments (e.g. antibody titres); the value is the measured
    /// titre, reported as below or above the limits of quantitation outside them.
    Is,
}

/// One measurement (one LB, VS, FT or OE record).
#[derive(Clone, Debug, PartialEq)]
pub struct Finding {
    /// Domain.
    pub domain: FindingDomain,
    /// `--TESTCD`.
    pub testcd: String,
    /// The visit it was taken at.
    pub visitnum: i64,
    /// Date.
    pub date: Date,
    /// Value in standard units, unrounded.
    pub value: f64,
    /// Laterality (`--LAT`, e.g. `LEFT`) of a finding on a paired organ.
    pub laterality: Option<&'static str>,
    /// Decimals reported, overriding the test's default.
    pub dp: Option<u32>,
}

/// What a lesion is, in RECIST 1.1 terms (`TUSTRESC`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum LesionKind {
    /// Measured target lesion.
    Target,
    /// Followed but not measured.
    NonTarget,
    /// New lesion (progression).
    New,
}

impl LesionKind {
    /// The `TUSTRESC` value.
    pub fn as_str(self) -> &'static str {
        match self {
            LesionKind::Target => "TARGET",
            LesionKind::NonTarget => "NON-TARGET",
            LesionKind::New => "NEW",
        }
    }
}

/// An identified lesion (one TU record).
#[derive(Clone, Debug, PartialEq)]
pub struct Lesion {
    /// Link ID, e.g. `T01`, `NT02`, `NEW01` (`TULNKID`).
    pub id: String,
    /// Kind.
    pub kind: LesionKind,
    /// Anatomical location (CT `LOC`).
    pub location: String,
    /// Laterality (CT `LAT`) for paired organs.
    pub laterality: Option<String>,
    /// Whether it is a lymph node (measured by short axis).
    pub node: bool,
    /// Visit it was identified at.
    pub visitnum: i64,
    /// Date it was identified.
    pub date: Date,
}

/// One tumour assessment (TR records, and RS records unless it is the baseline).
#[derive(Clone, Debug, PartialEq)]
pub struct TumorAssessment {
    /// Visit.
    pub visitnum: i64,
    /// Date (scan date).
    pub date: Date,
    /// Target lesion measurements in mm: (lesion ID, value; `None` = not evaluable).
    pub target: Vec<(String, Option<f64>)>,
    /// Sum of target diameters, if every target lesion was measured.
    pub sum: Option<f64>,
    /// Non-target lesion states: (lesion ID, `PRESENT` / `ABSENT` / `UNEQUIVOCAL PROGRESSION`).
    pub non_target: Vec<(String, String)>,
    /// New lesions first seen at this assessment.
    pub new_lesions: Vec<String>,
    /// Responses; `None` at baseline. RECIST 1.1: (target, non-target, overall). Lugano 2014:
    /// (CT-based, `None`, PET-CT-based).
    pub response: Option<(String, Option<String>, String)>,
    /// Lugano: target lesion perpendicular diameters in mm (`LPERP`), as `target` holds the
    /// longest diameters; empty under RECIST.
    pub perpendicular: Vec<(String, Option<f64>)>,
    /// Lugano: the FDG-PET Deauville five-point score (`None` if not done or under RECIST).
    pub deauville: Option<u8>,
}

/// A clinical event (one CE record): an adjudicated endpoint event or a disease exacerbation.
#[derive(Clone, Debug, PartialEq)]
pub struct ClinicalEvent {
    /// Preferred Term (`CEDECOD`).
    pub pt: String,
    /// Reported term (`CETERM`).
    pub term: String,
    /// `CECAT`.
    pub category: &'static str,
    /// Date of onset.
    pub date: Date,
    /// Date of resolution; `None` if not recorded or ongoing.
    pub end: Option<Date>,
    /// Severity (`CESEV`), if graded.
    pub severity: Option<Severity>,
}

/// One day of a solicited reaction in a post-vaccination diary.
#[derive(Clone, Debug, PartialEq)]
pub struct DiaryDay {
    /// Date.
    pub date: Date,
    /// Whether the reaction was present.
    pub occurred: bool,
    /// Graded severity (graded reactions, when present).
    pub severity: Option<Severity>,
    /// Largest diameter in cm (redness and swelling, when present).
    pub diameter: Option<f64>,
}

/// A solicited reaction followed in a post-vaccination diary (FACE; summarised in CE).
#[derive(Clone, Debug, PartialEq)]
pub struct Reaction {
    /// Reported term (`CETERM`, `FAOBJ`), e.g. `REDNESS`.
    pub term: &'static str,
    /// Preferred Term (`CEDECOD`).
    pub pt: &'static str,
    /// Administration-site reaction (else systemic).
    pub site: bool,
    /// Measured by diameter (redness, swelling) rather than graded.
    pub measured: bool,
    /// One entry per diary day.
    pub days: Vec<DiaryDay>,
}

/// The reactogenicity diary after one vaccination.
#[derive(Clone, Debug, PartialEq)]
pub struct Diary {
    /// The vaccination (1-based).
    pub vaccination: u8,
    /// Vaccination time (`--RFTDTC`).
    pub reference: DateTime,
    /// Injection site (`--LOC`, `--LAT`).
    pub location: (&'static str, &'static str),
    /// Solicited reactions, fever excluded.
    pub reactions: Vec<Reaction>,
    /// Daily body temperature (°C), for fever.
    pub temperatures: Vec<(Date, f64)>,
}

/// A pharmacokinetic sample (one PC record).
#[derive(Clone, Debug, PartialEq)]
pub struct PkSample {
    /// Analyte name (`PCTEST`; also the treatment name, `EXTRT`).
    pub analyte: String,
    /// The visit of the reference dose (`VISITNUM`).
    pub visitnum: i64,
    /// Planned time point (`PCTPT`), e.g. `Pre-dose` or `2h Post-dose`.
    pub tpt: String,
    /// Planned hours after the reference dose (`PCTPTNUM`); -0.5 before the first dose.
    pub tptnum: f64,
    /// Reference dose time (`PCRFTDTC`).
    pub reference: DateTime,
    /// Actual sampling time (`PCDTC`).
    pub time: DateTime,
    /// Measured concentration (standard units); `None` below the limit of quantitation.
    pub conc: Option<f64>,
    /// Lower limit of quantitation (standard units).
    pub lloq: f64,
}

/// A subject's whole observed course.
#[derive(Clone, Debug, PartialEq)]
pub struct Timeline {
    /// Study-treatment intervals, in order. Never empty.
    pub exposure: Vec<Exposure>,
    /// Elements, in order.
    pub elements: Vec<ElementSpan>,
    /// Attended visits, in date order.
    pub visits: Vec<VisitOccurrence>,
    /// Disposition events, in date order.
    pub disposition: Vec<DispositionEvent>,
    /// Date of death, if observed.
    pub death: Option<Date>,
    /// End of participation (last contact).
    pub last_contact: Date,
    /// The last day of the treatment epoch.
    pub treatment_end: Date,
    /// Adverse events.
    pub adverse_events: Vec<AdverseEvent>,
    /// Medications.
    pub medications: Vec<Medication>,
    /// Medical history.
    pub history: Vec<HistoryItem>,
    /// Lab and vital-sign measurements.
    pub findings: Vec<Finding>,
    /// Identified lesions (oncology).
    pub lesions: Vec<Lesion>,
    /// Tumour assessments in date order, baseline first (oncology).
    pub tumor: Vec<TumorAssessment>,
    /// Pharmacokinetic samples in time order.
    pub pk: Vec<PkSample>,
    /// Adjudicated clinical events in date order.
    pub clinical_events: Vec<ClinicalEvent>,
    /// The study eye (`LEFT` or `RIGHT`) of an ophthalmology study: the treated eye.
    pub study_eye: Option<&'static str>,
    /// Date of birth (`BRTHDTC`), recorded in paediatric studies.
    pub birth_date: Option<Date>,
    /// Reactogenicity diaries, one per vaccination.
    pub diaries: Vec<Diary>,
}

impl Timeline {
    /// First dose date.
    pub fn first_dose(&self) -> Date {
        self.exposure.first().expect("at least one exposure").start
    }

    /// Last dose date.
    pub fn last_dose(&self) -> Date {
        self.exposure.last().expect("at least one exposure").end
    }

    /// The epoch a date falls in: before first dose → screening; up to the end of treatment →
    /// treatment (for a crossover, the last period or washout element starting on or before
    /// the date); after → follow-up.
    pub fn epoch_of(&self, date: Date) -> Epoch {
        // After a switch: the element containing the date (the treatment element includes its
        // last day; later elements start on their first).
        if self.elements.iter().any(|e| e.epoch == Epoch::Continuation) {
            if date < self.first_dose() {
                return Epoch::Screening;
            }
            let n = self.elements.len();
            return self
                .elements
                .iter()
                .enumerate()
                .filter(|(_, e)| e.epoch != Epoch::Screening)
                .find(|(k, e)| {
                    e.start <= date
                        && (date < e.end
                            || (date == e.end && (e.epoch == Epoch::Treatment || k + 1 == n)))
                })
                .map_or(Epoch::FollowUp, |(_, e)| e.epoch);
        }
        if self
            .elements
            .iter()
            .any(|e| matches!(e.epoch, Epoch::Period(_)))
        {
            if date < self.first_dose() {
                return Epoch::Screening;
            }
            if date > self.treatment_end {
                return Epoch::FollowUp;
            }
            return self
                .elements
                .iter()
                .filter(|e| matches!(e.epoch, Epoch::Period(_) | Epoch::Washout(_)))
                .rev()
                .find(|e| e.start <= date)
                .map_or(Epoch::Screening, |e| e.epoch);
        }
        if date < self.first_dose() {
            Epoch::Screening
        } else if date <= self.treatment_end {
            Epoch::Treatment
        } else {
            Epoch::FollowUp
        }
    }
}
