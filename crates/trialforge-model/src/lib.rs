//! The study model: design (arms, sites, randomisation, calendar), enrollment, and subjects
//! with their baseline covariates.
//!
//! Every per-subject draw comes from that subject's own stream (`study → "subject" → i`), so a
//! subject's baseline does not depend on other subjects. Enrollment order and randomisation
//! are study-level sequential processes: changing `n` reorders enrollment, and permuted blocks
//! assign arms in that order.

use trialforge_core::{Date, Draws, Stream};

/// A study treatment as administered (feeds EX).
#[derive(Clone, Debug, PartialEq)]
pub struct Treatment {
    /// `EXTRT`, e.g. `DRUG A`.
    pub name: String,
    /// `EXDOSE`.
    pub dose: f64,
    /// `EXDOSU`, e.g. `mg`.
    pub unit: String,
    /// `EXDOSFRM`, e.g. `SOLUTION`.
    pub form: String,
    /// `EXROUTE`, e.g. `INTRAVENOUS`.
    pub route: String,
    /// `EXDOSFRQ`, e.g. `Q3W`.
    pub frequency: String,
    /// Days between administrations.
    pub cycle_days: i32,
}

/// A treatment period of a crossover sequence.
#[derive(Clone, Debug, PartialEq)]
pub struct Period {
    /// `ETCD` of the period's element (one per treatment, shared across sequences).
    pub etcd: String,
    /// `ELEMENT` of the period's element, e.g. `TFX-101 100 mg tablet (test)`.
    pub element: String,
    /// The treatment given in the period.
    pub treatment: Treatment,
}

/// A study arm.
#[derive(Clone, Debug, PartialEq)]
pub struct Arm {
    /// `ARMCD`.
    pub armcd: String,
    /// `ARM`.
    pub arm: String,
    /// Allocation ratio weight (e.g. 1 for 1:1).
    pub ratio: u32,
    /// The treatment given in this arm (the first period's, for a crossover sequence).
    pub treatment: Treatment,
    /// Treatment periods, in order, when the arm is a crossover sequence; empty otherwise.
    pub periods: Vec<Period>,
    /// A treatment subjects may switch to after progression (crossover on progression), and
    /// its element name.
    pub switch_to: Option<(Treatment, String)>,
}

/// A site with its country (ISO 3166-1 alpha-3).
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    /// `SITEID`.
    pub siteid: String,
    /// `COUNTRY`.
    pub country: String,
}

/// A categorical distribution over labelled values.
#[derive(Clone, Debug, PartialEq)]
pub struct Categorical {
    /// (value, weight) pairs. Weights need not sum to 1.
    pub levels: Vec<(String, f64)>,
}

impl Categorical {
    /// A distribution from (value, weight) pairs.
    pub fn new<S: Into<String>>(levels: impl IntoIterator<Item = (S, f64)>) -> Self {
        Self {
            levels: levels.into_iter().map(|(v, w)| (v.into(), w)).collect(),
        }
    }

    fn draw(&self, d: &mut Draws) -> usize {
        let weights: Vec<f64> = self.levels.iter().map(|(_, w)| *w).collect();
        d.weighted(&weights)
    }
}

/// Baseline population distributions.
#[derive(Clone, Debug, PartialEq)]
pub struct Population {
    /// Mean age in years.
    pub age_mean: f64,
    /// Standard deviation of age.
    pub age_sd: f64,
    /// Minimum age (inclusive).
    pub age_min: i64,
    /// Maximum age (inclusive).
    pub age_max: i64,
    /// `SEX` (CT C66731).
    pub sex: Categorical,
    /// `RACE` (CT C74457).
    pub race: Categorical,
    /// `ETHNIC` (CT C66790).
    pub ethnic: Categorical,
    /// The randomisation stratification factor; level names become `STRAT1R`.
    pub strata: Categorical,
}

/// A trial epoch (CT `EPOCH`).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Epoch {
    /// Screening.
    Screening,
    /// Treatment.
    Treatment,
    /// Follow-up.
    FollowUp,
    /// Treatment period `n` (1-based) of a crossover.
    Period(u8),
    /// The washout after period `n` of a crossover (`n = 0` when there is a single washout).
    Washout(u8),
    /// Treatment after switching, e.g. crossover to the experimental drug on progression.
    Continuation,
}

impl Epoch {
    /// The CT submission value.
    pub fn as_str(self) -> &'static str {
        match self {
            Epoch::Screening => "SCREENING",
            Epoch::Treatment => "TREATMENT",
            Epoch::FollowUp => "FOLLOW-UP",
            Epoch::Period(n) => match n {
                1 => "TREATMENT 1",
                2 => "TREATMENT 2",
                3 => "TREATMENT 3",
                4 => "TREATMENT 4",
                _ => "TREATMENT",
            },
            Epoch::Continuation => "CONTINUATION TREATMENT",
            Epoch::Washout(n) => match n {
                1 => "WASHOUT 1",
                2 => "WASHOUT 2",
                3 => "WASHOUT 3",
                _ => "WASHOUT",
            },
        }
    }
}

/// What a planned visit is for.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VisitKind {
    /// Screening (at consent).
    Screening,
    /// Baseline (Day 1, first dose).
    Baseline,
    /// A scheduled on-treatment visit.
    Treatment,
    /// The last scheduled on-treatment visit.
    EndOfTreatment,
    /// Post-treatment follow-up.
    FollowUp,
    /// Early termination (unscheduled; happens when a subject leaves early).
    EarlyTermination,
}

/// A planned visit (one TV record).
#[derive(Clone, Debug, PartialEq)]
pub struct PlannedVisit {
    /// `VISITNUM`.
    pub visitnum: i64,
    /// `VISIT`.
    pub name: String,
    /// Planned study day (`VISITDY`); `None` for early termination.
    pub day: Option<i32>,
    /// Allowed deviation from the planned day: (days early, days late).
    pub window: (i32, i32),
    /// The epoch the visit falls in.
    pub epoch: Epoch,
    /// What the visit is for.
    pub kind: VisitKind,
    /// Whether the safety lab panel is drawn.
    pub labs: bool,
    /// Whether vital signs are taken.
    pub vitals: bool,
    /// Whether the efficacy endpoint is measured.
    pub endpoint: bool,
}

/// The visit schedule of a fixed-duration study.
#[derive(Clone, Debug, PartialEq)]
pub struct Schedule {
    /// Planned visits in `VISITNUM` order, including early termination.
    pub visits: Vec<PlannedVisit>,
    /// Planned treatment duration in days (Day 1 to the end-of-treatment visit).
    pub treatment_days: i32,
    /// Planned follow-up duration in days after the last dose.
    pub followup_days: i32,
}

impl Schedule {
    /// The visit with the given `VISITNUM`.
    pub fn visit(&self, visitnum: i64) -> Option<&PlannedVisit> {
        self.visits.iter().find(|v| v.visitnum == visitnum)
    }
}

/// Trial-level descriptors (feed TS and TI).
#[derive(Clone, Debug, PartialEq)]
pub struct TrialInfo {
    /// `TITLE`.
    pub title: String,
    /// `TPHASE` (CT value, e.g. `PHASE III TRIAL`).
    pub phase: String,
    /// `SPONSOR`.
    pub sponsor: String,
    /// `TDIGRP`: the condition studied.
    pub indication: String,
    /// `THERAREA`.
    pub therapeutic_area: String,
    /// `OBJPRIM`.
    pub primary_objective: String,
    /// `TBLIND` (e.g. `DOUBLE BLIND`, `OPEN LABEL`).
    pub blinding: String,
    /// `TCNTRL` (e.g. `PLACEBO`, `ACTIVE`).
    pub control: String,
    /// `ADDON`: added on to existing treatment.
    pub add_on: bool,
    /// `TTYPE` values (CT), e.g. `EFFICACY`, or `SAFETY` and `PHARMACOKINETIC`.
    pub trial_types: Vec<String>,
    /// `INTMODEL` (CT), e.g. `PARALLEL` or `SEQUENTIAL`.
    pub intervention_model: String,
    /// `HLTSUBJI`: whether the subjects are healthy volunteers.
    pub healthy_subjects: bool,
    /// `INTTYPE` (CT), e.g. `DRUG` or `VACCINE`.
    pub intervention_type: String,
    /// Inclusion and exclusion criteria: (`IETESTCD`, `IECAT`, `IETEST`).
    pub criteria: Vec<(String, String, String)>,
}

/// A sequential enrollment cohort (e.g. a dose level of a phase 1 study): its subjects are
/// randomised within the cohort to fixed numbers per arm.
#[derive(Clone, Debug, PartialEq)]
pub struct Cohort {
    /// (arm index, number of subjects) in the cohort.
    pub arms: Vec<(usize, u32)>,
    /// First possible randomisation, in days after [`Design::first_randomization`].
    pub start_day: i32,
    /// Length of the cohort's enrollment window in days.
    pub enrollment_days: i32,
}

impl Cohort {
    /// Subjects in the cohort.
    pub fn size(&self) -> usize {
        self.arms.iter().map(|(_, k)| *k as usize).sum()
    }
}

/// Day offsets around randomisation, each drawn uniformly from an inclusive range.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Timing {
    /// Days from the screening assessments to randomisation.
    pub screening_lead: (i64, i64),
    /// Extra days from informed consent to screening.
    pub consent_lead: (i64, i64),
    /// Days from randomisation to first dose.
    pub dose_lag: (i64, i64),
}

/// Study-level design.
#[derive(Clone, Debug, PartialEq)]
pub struct Design {
    /// `STUDYID`.
    pub studyid: String,
    /// Root seed.
    pub seed: u64,
    /// Number of randomised subjects.
    pub n: usize,
    /// Arms, in display order.
    pub arms: Vec<Arm>,
    /// Sites; subjects are spread across them uniformly.
    pub sites: Vec<Site>,
    /// Date of the first possible randomisation.
    pub first_randomization: Date,
    /// Length of the enrollment window in days.
    pub enrollment_days: i32,
    /// Data cutoff date.
    pub cutoff: Date,
    /// The cutoff is an interim snapshot rather than the primary analysis.
    pub interim: bool,
    /// Permuted block size (a multiple of the summed arm ratios).
    pub block_size: u32,
    /// Baseline distributions.
    pub population: Population,
    /// Screening, consent and dosing offsets.
    pub timing: Timing,
    /// Visit schedule; `None` for event-driven designs without fixed visits.
    pub schedule: Option<Schedule>,
    /// Sequential cohorts; when non-empty they replace permuted-block randomisation and their
    /// sizes add up to `n`.
    pub cohorts: Vec<Cohort>,
    /// Trial-level descriptors.
    pub info: TrialInfo,
}

impl Design {
    /// The root stream for this design.
    pub fn stream(&self) -> Stream {
        Stream::root(self.seed)
    }

    /// Subject `i`'s stream.
    pub fn subject_stream(&self, i: usize) -> Stream {
        self.stream().label("subject").index(i as u64)
    }

    /// Checks the design for internal consistency.
    pub fn validate(&self) -> Result<(), DesignError> {
        let err = |m: String| Err(DesignError(m));
        if self.n == 0 {
            return err("`n` must be at least 1".into());
        }
        if self.arms.is_empty() || self.arms.iter().any(|a| a.ratio == 0) {
            return err("every arm needs a positive ratio".into());
        }
        let ratio_sum: u32 = self.arms.iter().map(|a| a.ratio).sum();
        if !self.cohorts.is_empty() {
            let total: usize = self.cohorts.iter().map(Cohort::size).sum();
            if total != self.n {
                return err(format!(
                    "cohort sizes add up to {total}, not n = {}",
                    self.n
                ));
            }
            if self.cohorts.iter().any(|c| {
                c.enrollment_days < 1
                    || c.start_day < 0
                    || c.arms.iter().any(|(a, _)| *a >= self.arms.len())
            }) {
                return err(
                    "every cohort needs valid arms and a positive enrollment window".into(),
                );
            }
        } else if self.block_size == 0 || self.block_size % ratio_sum != 0 {
            return err(format!(
                "block size {} is not a multiple of the summed arm ratios ({ratio_sum})",
                self.block_size
            ));
        }
        if self.sites.is_empty() {
            return err("the design needs at least one site".into());
        }
        if self.enrollment_days < 1 {
            return err("the enrollment window must be at least one day".into());
        }
        if self.cutoff <= self.first_randomization.add_days(self.enrollment_days) {
            return err("the data cutoff must fall after the end of enrollment".into());
        }
        let t = &self.timing;
        for (name, (lo, hi)) in [
            ("screening_lead", t.screening_lead),
            ("consent_lead", t.consent_lead),
            ("dose_lag", t.dose_lag),
        ] {
            if lo < 0 || lo > hi {
                return err(format!(
                    "timing `{name}` must be a non-negative range lo..=hi"
                ));
            }
        }
        if t.screening_lead.0 < 1 {
            return err("screening must happen at least one day before randomisation".into());
        }
        if let Some(sch) = &self.schedule {
            if !sch.visits.windows(2).all(|w| w[0].visitnum < w[1].visitnum) {
                return err("schedule visits must be in increasing VISITNUM order".into());
            }
            let scheduled: Vec<i32> = sch.visits.iter().filter_map(|v| v.day).collect();
            if !scheduled.windows(2).all(|w| w[0] < w[1]) {
                return err("scheduled visit days must increase".into());
            }
            if sch
                .visits
                .iter()
                .filter(|v| v.kind == VisitKind::Screening)
                .count()
                != 1
            {
                return err("the schedule needs exactly one screening visit".into());
            }
        }
        if self.studyid.is_empty() || self.studyid.len() > 20 {
            return err("`studyid` must be 1–20 characters".into());
        }
        Ok(())
    }
}

/// An invalid design.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DesignError(pub String);

impl core::fmt::Display for DesignError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(&self.0)
    }
}

impl std::error::Error for DesignError {}

/// An enrolled, randomised subject and their baseline.
#[derive(Clone, Debug, PartialEq)]
pub struct Subject {
    /// Index of the subject's stream (stable for a given design).
    pub index: usize,
    /// `SUBJID`: site ID plus a three-digit sequence within the site.
    pub subjid: String,
    /// `USUBJID`: `STUDYID-SUBJID`.
    pub usubjid: String,
    /// Index into [`Design::sites`].
    pub site: usize,
    /// `AGE` in years.
    pub age: i64,
    /// Index into [`Population::sex`].
    pub sex: usize,
    /// Index into [`Population::race`].
    pub race: usize,
    /// Index into [`Population::ethnic`].
    pub ethnic: usize,
    /// Index into [`Population::strata`].
    pub stratum: usize,
    /// Index into [`Design::arms`].
    pub arm: usize,
    /// Informed consent date.
    pub consent: Date,
    /// Screening assessment date (the screening visit).
    pub screening: Date,
    /// Randomisation date.
    pub randomization: Date,
    /// First dose date.
    pub treatment_start: Date,
}

/// Enrolls and randomises `design.n` subjects. The result is in enrollment order.
pub fn enroll(design: &Design) -> Result<Vec<Subject>, DesignError> {
    design.validate()?;
    let pop = &design.population;
    let mut subjects: Vec<Subject> = (0..design.n)
        .map(|i| {
            let mut d = Draws::new(design.subject_stream(i).label("baseline"));
            // Draw order is part of the output contract; append new draws at the end.
            let site = d.below(design.sites.len() as u64) as usize;
            let age =
                (d.normal(pop.age_mean, pop.age_sd).round() as i64).clamp(pop.age_min, pop.age_max);
            let sex = pop.sex.draw(&mut d);
            let race = pop.race.draw(&mut d);
            let ethnic = pop.ethnic.draw(&mut d);
            let stratum = pop.strata.draw(&mut d);
            let rand_offset = d.below(design.enrollment_days as u64) as i32;
            // In a cohort design the subject's cohort sets the window (an independent draw on
            // the same stream position would shift every later draw, so reuse the uniform).
            let rand_offset = cohort_of(design, i).map_or(rand_offset, |c| {
                c.start_day + (rand_offset as i64 % i64::from(c.enrollment_days)) as i32
            });
            let t = &design.timing;
            let scan_lead = d.int_range(t.screening_lead.0, t.screening_lead.1) as i32;
            let consent_lead = d.int_range(t.consent_lead.0, t.consent_lead.1) as i32;
            let dose_lag = d.int_range(t.dose_lag.0, t.dose_lag.1) as i32;
            let randomization = design.first_randomization.add_days(rand_offset);
            Subject {
                index: i,
                subjid: String::new(),
                usubjid: String::new(),
                site,
                age,
                sex,
                race,
                ethnic,
                stratum,
                arm: usize::MAX,
                consent: randomization.add_days(-(scan_lead + consent_lead)),
                screening: randomization.add_days(-scan_lead),
                randomization,
                treatment_start: randomization.add_days(dose_lag),
            }
        })
        .collect();

    // Enrollment order: by randomisation date, ties broken by stream index.
    subjects.sort_by_key(|s| (s.randomization, s.index));

    // Subject IDs: site ID + sequence within the site, in enrollment order.
    let mut per_site = vec![0u32; design.sites.len()];
    for s in &mut subjects {
        per_site[s.site] += 1;
        s.subjid = format!("{}{:03}", design.sites[s.site].siteid, per_site[s.site]);
        s.usubjid = format!("{}-{}", design.studyid, s.subjid);
    }

    randomize(design, &mut subjects);
    Ok(subjects)
}

/// The cohort of subject stream `i` (cohorts take consecutive stream indices).
fn cohort_of(design: &Design, i: usize) -> Option<&Cohort> {
    let mut start = 0;
    for c in &design.cohorts {
        if i < start + c.size() {
            return Some(c);
        }
        start += c.size();
    }
    None
}

/// Stratified permuted-block randomisation, in enrollment order; or, for cohort designs, a
/// shuffle of each cohort's fixed allocation.
fn randomize(design: &Design, subjects: &mut [Subject]) {
    if !design.cohorts.is_empty() {
        let stream = design.stream().label("randomization").label("cohort");
        let mut start = 0;
        for (k, c) in design.cohorts.iter().enumerate() {
            let mut slots: Vec<usize> = c
                .arms
                .iter()
                .flat_map(|(a, n)| std::iter::repeat_n(*a, *n as usize))
                .collect();
            Draws::new(stream.index(k as u64)).shuffle(&mut slots);
            let mut members: Vec<&mut Subject> = subjects
                .iter_mut()
                .filter(|s| s.index >= start && s.index < start + c.size())
                .collect();
            members.sort_by_key(|s| (s.randomization, s.index));
            for (s, arm) in members.into_iter().zip(slots) {
                s.arm = arm;
            }
            start += c.size();
        }
        return;
    }
    let ratio_sum: u32 = design.arms.iter().map(|a| a.ratio).sum();
    let reps = design.block_size / ratio_sum;
    let template: Vec<usize> = design
        .arms
        .iter()
        .enumerate()
        .flat_map(|(i, a)| std::iter::repeat_n(i, (a.ratio * reps) as usize))
        .collect();
    let stream = design.stream().label("randomization");
    for (k, (stratum, _)) in design.population.strata.levels.iter().enumerate() {
        let mut block_no = 0u64;
        let mut block: Vec<usize> = Vec::new();
        for s in subjects.iter_mut().filter(|s| s.stratum == k) {
            if block.is_empty() {
                block = template.clone();
                Draws::new(stream.label(stratum).index(block_no)).shuffle(&mut block);
                block.reverse(); // pop() then deals from the front of the shuffled block
                block_no += 1;
            }
            s.arm = block.pop().expect("block refilled above");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn design(n: usize) -> Design {
        let trt = |name: &str| Treatment {
            name: name.into(),
            dose: 1.0,
            unit: "mg".into(),
            form: "TABLET".into(),
            route: "ORAL".into(),
            frequency: "QD".into(),
            cycle_days: 1,
        };
        Design {
            studyid: "TEST01".into(),
            seed: 42,
            n,
            arms: vec![
                Arm {
                    armcd: "A".into(),
                    arm: "Drug A".into(),
                    ratio: 2,
                    treatment: trt("A"),
                    periods: Vec::new(),
                    switch_to: None,
                },
                Arm {
                    armcd: "P".into(),
                    arm: "Placebo".into(),
                    ratio: 1,
                    treatment: trt("P"),
                    periods: Vec::new(),
                    switch_to: None,
                },
            ],
            sites: vec![
                Site {
                    siteid: "101".into(),
                    country: "USA".into(),
                },
                Site {
                    siteid: "102".into(),
                    country: "CAN".into(),
                },
            ],
            first_randomization: Date::from_ymd(2024, 1, 1).unwrap(),
            enrollment_days: 365,
            cutoff: Date::from_ymd(2026, 1, 1).unwrap(),
            interim: false,
            block_size: 6,
            population: Population {
                age_mean: 60.0,
                age_sd: 10.0,
                age_min: 18,
                age_max: 85,
                sex: Categorical::new([("F", 1.0), ("M", 1.0)]),
                race: Categorical::new([("WHITE", 1.0)]),
                ethnic: Categorical::new([("NOT HISPANIC OR LATINO", 1.0)]),
                strata: Categorical::new([("LOW", 1.0), ("HIGH", 1.0)]),
            },
            timing: Timing {
                screening_lead: (1, 14),
                consent_lead: (0, 14),
                dose_lag: (0, 3),
            },
            schedule: None,
            cohorts: Vec::new(),
            info: TrialInfo {
                title: "Test".into(),
                phase: "PHASE III TRIAL".into(),
                sponsor: "Test".into(),
                indication: "Test".into(),
                therapeutic_area: "Test".into(),
                primary_objective: "Test".into(),
                blinding: "DOUBLE BLIND".into(),
                control: "PLACEBO".into(),
                add_on: false,
                trial_types: vec!["EFFICACY".into()],
                intervention_model: "PARALLEL".into(),
                healthy_subjects: false,
                intervention_type: "DRUG".into(),
                criteria: vec![],
            },
        }
    }

    #[test]
    fn enrollment_is_deterministic_and_ordered() {
        let a = enroll(&design(120)).unwrap();
        assert_eq!(a, enroll(&design(120)).unwrap());
        assert!(
            a.windows(2)
                .all(|w| w[0].randomization <= w[1].randomization)
        );
        assert!(a.iter().all(|s| s.consent < s.randomization));
        assert!(a.iter().all(|s| s.treatment_start >= s.randomization));
    }

    #[test]
    fn subject_ids_are_unique() {
        let subjects = enroll(&design(300)).unwrap();
        let mut ids: Vec<&str> = subjects.iter().map(|s| s.usubjid.as_str()).collect();
        ids.sort_unstable();
        ids.dedup();
        assert_eq!(ids.len(), 300);
    }

    #[test]
    fn blocks_are_balanced_within_strata() {
        let subjects = enroll(&design(300)).unwrap();
        for k in 0..2 {
            let in_stratum: Vec<_> = subjects.iter().filter(|s| s.stratum == k).collect();
            for block in in_stratum.chunks(6).filter(|b| b.len() == 6) {
                let on_a = block.iter().filter(|s| s.arm == 0).count();
                assert_eq!(on_a, 4, "2:1 allocation in every complete block of 6");
            }
        }
    }

    #[test]
    fn baselines_are_keyed_by_subject_stream() {
        let small = enroll(&design(50)).unwrap();
        let large = enroll(&design(51)).unwrap();
        for s in &small {
            let t = large.iter().find(|t| t.index == s.index).unwrap();
            assert_eq!(
                (s.age, s.sex, s.stratum, s.randomization),
                (t.age, t.sex, t.stratum, t.randomization)
            );
        }
    }

    #[test]
    fn cohorts_enrol_in_sequence_with_fixed_allocation() {
        let mut d = design(16);
        d.cohorts = vec![
            Cohort {
                arms: vec![(0, 6), (1, 2)],
                start_day: 0,
                enrollment_days: 7,
            },
            Cohort {
                arms: vec![(0, 6), (1, 2)],
                start_day: 28,
                enrollment_days: 7,
            },
        ];
        let s = enroll(&d).unwrap();
        for (k, c) in s.chunks(8).enumerate() {
            assert_eq!(c.iter().filter(|x| x.arm == 0).count(), 6);
            assert!(c.iter().all(|x| x.index / 8 == k));
            let lo = d.first_randomization.add_days(28 * k as i32);
            assert!(
                c.iter()
                    .all(|x| x.randomization >= lo && x.randomization < lo.add_days(7))
            );
        }
    }

    #[test]
    fn invalid_designs_are_rejected() {
        let mut d = design(10);
        d.block_size = 4;
        assert!(enroll(&d).is_err());
    }
}
