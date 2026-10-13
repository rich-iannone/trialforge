//! A subject's observed course, independent of any SDTM layout.
//!
//! Every process engine produces a [`Timeline`] per subject; the SDTM emitters read only
//! timelines (plus the design), so all archetypes share one observation layer.

use trialforge_core::Date;
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
}

/// One measurement (one LB or VS record).
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
    /// treatment; after → follow-up.
    pub fn epoch_of(&self, date: Date) -> Epoch {
        if date < self.first_dose() {
            Epoch::Screening
        } else if date <= self.treatment_end {
            Epoch::Treatment
        } else {
            Epoch::FollowUp
        }
    }
}
