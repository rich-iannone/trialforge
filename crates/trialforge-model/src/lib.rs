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

/// A study arm.
#[derive(Clone, Debug, PartialEq)]
pub struct Arm {
    /// `ARMCD`.
    pub armcd: String,
    /// `ARM`.
    pub arm: String,
    /// Allocation ratio weight (e.g. 1 for 1:1).
    pub ratio: u32,
    /// The treatment given in this arm.
    pub treatment: Treatment,
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
    /// Permuted block size (a multiple of the summed arm ratios).
    pub block_size: u32,
    /// Baseline distributions.
    pub population: Population,
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
        if self.block_size == 0 || self.block_size % ratio_sum != 0 {
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
            let scan_lead = d.int_range(1, 14) as i32;
            let consent_lead = d.int_range(0, 14) as i32;
            let dose_lag = d.int_range(0, 3) as i32;
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

/// Stratified permuted-block randomisation, in enrollment order.
fn randomize(design: &Design, subjects: &mut [Subject]) {
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
                },
                Arm {
                    armcd: "P".into(),
                    arm: "Placebo".into(),
                    ratio: 1,
                    treatment: trt("P"),
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
    fn invalid_designs_are_rejected() {
        let mut d = design(10);
        d.block_size = 4;
        assert!(enroll(&d).is_err());
    }
}
