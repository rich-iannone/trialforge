//! Vaccinations, reactogenicity diaries and antibody titres, layered on a fixed-duration
//! course.
//!
//! - **Vaccinations**: an intramuscular injection into the deltoid at each attended visit on a
//!   planned vaccination day (up to the course's last dose), with a clock time.
//! - **Reactogenicity**: for 7 days from each vaccination (day 1 is the vaccination day) the
//!   subject records each solicited reaction. A reaction occurs after a vaccination with its
//!   arm's probability (systemic reactions to the active vaccine are more frequent after the
//!   second dose); it starts on day 1–3 and lasts a few days. Graded reactions take a peak
//!   severity that eases over the days; redness and swelling are measured as a diameter
//!   (0.5 cm steps). Body temperature is recorded daily; fever is a temperature of 38.0 °C or
//!   more.
//! - **Antibody titres** at the assay visits: a baseline titre (seropositive or below the
//!   LLOQ), a log-normal response after the last dose with the arm's geometric mean (part of
//!   it after the first dose), waning afterwards, and assay noise. Titres are reported as the
//!   reciprocal of the last positive two-fold dilution, below the LLOQ or above the ULOQ.

use trialforge_core::{Date, DateTime, Draws, math};
use trialforge_model::{Design, Subject};

use crate::timeline::{
    Diary, DiaryDay, Exposure, Finding, FindingDomain, Reaction, Severity, Timeline,
};

/// Days recorded in a reactogenicity diary.
pub const DIARY_DAYS: i32 = 7;

/// A solicited reaction and its probability per vaccination, per arm.
#[derive(Clone, Debug, PartialEq)]
pub struct Solicited {
    /// `CETERM`/`FAOBJ`.
    pub term: &'static str,
    /// `CEDECOD`.
    pub pt: &'static str,
    /// Administration-site reaction.
    pub site: bool,
    /// Measured by diameter.
    pub measured: bool,
    /// Probability per vaccination, per arm.
    pub rate: Vec<f64>,
}

/// Vaccine parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct VaccineParams {
    /// Planned study days of the vaccinations.
    pub vaccination_days: Vec<i32>,
    /// Planned study days of the antibody assays.
    pub assay_days: Vec<i32>,
    /// Dose (mL per injection).
    pub dose: f64,
    /// Solicited reactions (fever excluded).
    pub solicited: Vec<Solicited>,
    /// Probability of fever per vaccination, per arm.
    pub fever_rate: Vec<f64>,
    /// Multiplier on systemic reaction and fever probabilities after a later dose, per arm.
    pub later_dose_multiplier: Vec<f64>,
    /// Share of subjects seropositive at baseline.
    pub seropositive: f64,
    /// Geometric mean titre 28 days after the last dose, per arm (placebo: baseline).
    pub gmt_by_arm: Vec<f64>,
    /// Whether each arm's vaccine induces a response.
    pub responds: Vec<bool>,
    /// Share of the log-titre rise reached after the first of several doses.
    pub first_dose_fraction: f64,
    /// SD of individual log10 titres after vaccination.
    pub log10_sd: f64,
    /// Log10 titre lost per year after the last dose.
    pub waning_per_year: f64,
    /// Lower and upper limits of quantitation (titre).
    pub limits: (f64, f64),
}

/// The reported titre: the last two-fold dilution of `lloq` reached, 0 below the LLOQ,
/// infinite above the ULOQ.
fn reported(titre: f64, limits: (f64, f64)) -> f64 {
    if titre < limits.0 {
        0.0
    } else if titre >= limits.1 {
        f64::INFINITY
    } else {
        let steps = math::ln(titre / limits.0) / core::f64::consts::LN_2;
        limits.0 * math::pow(2.0, (steps + 1e-9).floor())
    }
}

/// Simulates vaccinations, diaries and titres onto `timelines` (in subject order).
pub fn simulate(
    design: &Design,
    p: &VaccineParams,
    subjects: &[Subject],
    timelines: &mut [Timeline],
) {
    let schedule = design
        .schedule
        .as_ref()
        .expect("a vaccine study has a schedule");
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let stream = design.subject_stream(s.index).label("vaccine");
        let mut d = Draws::new(stream.label("latent"));
        // Draw order is part of the output contract; append new draws at the end.
        let arm_side = if d.uniform() < 0.85 { "LEFT" } else { "RIGHT" };
        let seropositive = d.uniform() < p.seropositive;
        let baseline_log = if seropositive {
            d.normal(math::ln(4.0 * p.limits.0) / core::f64::consts::LN_10, 0.35)
        } else {
            d.normal(math::ln(p.limits.0 / 2.0) / core::f64::consts::LN_10, 0.1)
        };
        // Reported titres are floored to a dilution step: centre the latent half a step up.
        let gmt = p.gmt_by_arm.get(s.arm).copied().unwrap_or(p.limits.0);
        let peak_log = d.normal(
            math::ln(gmt * math::sqrt(2.0)) / core::f64::consts::LN_10,
            p.log10_sd,
        );
        let responds = p.responds.get(s.arm).copied().unwrap_or(false);

        // Vaccinations at attended visits on planned vaccination days.
        let last = t.last_dose();
        let mut vaccinations: Vec<(u8, i64, Date)> = Vec::new();
        for v in &t.visits {
            let Some(day) = schedule.visit(v.visitnum).and_then(|pv| pv.day) else {
                continue;
            };
            if let Some(k) = p.vaccination_days.iter().position(|&x| x == day) {
                if v.date <= last {
                    vaccinations.push((k as u8 + 1, v.visitnum, v.date));
                }
            }
        }
        if vaccinations.is_empty() {
            vaccinations.push((1, 0, t.first_dose()));
        }
        let mut exposure = Vec::new();
        let mut diaries = Vec::new();
        let n_doses = p.vaccination_days.len();
        for &(k, visitnum, date) in &vaccinations {
            let mut v = Draws::new(stream.label("vaccination").index(u64::from(k)));
            let minutes = v.int_range(8 * 60, 16 * 60) as i32;
            let time = DateTime::new(date, (minutes / 60) as u32, (minutes % 60) as u32);
            exposure.push(Exposure {
                start: date,
                end: date,
                dose: p.dose,
                time: Some(time),
                period: None,
                location: Some(("DELTOID MUSCLE", arm_side)),
                group: Some(format!("VACCINATION {k}")),
                visitnum: (visitnum > 0).then_some(visitnum),
                switched: false,
            });
            let later = if k > 1 {
                p.later_dose_multiplier.get(s.arm).copied().unwrap_or(1.0)
            } else {
                1.0
            };
            let days: Vec<Date> = (0..DIARY_DAYS)
                .map(|i| date.add_days(i))
                .filter(|d| *d <= t.last_contact)
                .collect();
            let mut reactions = Vec::new();
            for sol in &p.solicited {
                let mut r =
                    Draws::new(stream.label("reaction").index(u64::from(k)).label(sol.term));
                let rate = sol.rate.get(s.arm).copied().unwrap_or(0.0)
                    * if sol.site { 1.0 } else { later };
                let occurs = r.uniform() < rate.min(0.95);
                let onset = [0, 0, 0, 1, 1, 2][r.below(6) as usize];
                let duration = 1 + r.below(3) as i32 + i32::from(r.uniform() < 0.3);
                let peak = [Severity::Mild, Severity::Moderate, Severity::Severe]
                    [r.weighted(&[0.6, 0.32, 0.08])];
                let diameter = math::round_dp(2.0 * r.lognormal(3.5, 0.7), 0) / 2.0;
                let reaction_days = days
                    .iter()
                    .enumerate()
                    .map(|(i, &date)| {
                        let i = i as i32;
                        let on = occurs && i >= onset && i < onset + duration;
                        let since = i - onset;
                        let severity = (on && !sol.measured).then_some(match (peak, since / 2) {
                            (Severity::Severe, 0) => Severity::Severe,
                            (Severity::Severe, 1) | (Severity::Moderate, 0) => Severity::Moderate,
                            _ => Severity::Mild,
                        });
                        let diam = (on && sol.measured).then(|| {
                            (math::round_dp(2.0 * diameter * math::pow(0.7, f64::from(since)), 0)
                                / 2.0)
                                .max(0.5)
                        });
                        DiaryDay {
                            date,
                            occurred: on,
                            severity,
                            diameter: diam,
                        }
                    })
                    .collect();
                reactions.push(Reaction {
                    term: sol.term,
                    pt: sol.pt,
                    site: sol.site,
                    measured: sol.measured,
                    days: reaction_days,
                });
            }
            let mut f = Draws::new(stream.label("fever").index(u64::from(k)));
            let fever =
                f.uniform() < (p.fever_rate.get(s.arm).copied().unwrap_or(0.0) * later).min(0.95);
            let onset = f.below(2) as usize;
            let length = 1 + f.below(2) as usize;
            let temperatures = days
                .iter()
                .enumerate()
                .map(|(i, &date)| {
                    let normal = f.normal(36.6, 0.25).min(37.7);
                    let high = (38.0 + f.exponential(2.0)).min(40.5);
                    let febrile = fever && i >= onset && i < onset + length;
                    (date, math::round_dp(if febrile { high } else { normal }, 1))
                })
                .collect();
            diaries.push(Diary {
                vaccination: k,
                reference: time,
                location: ("DELTOID MUSCLE", arm_side),
                reactions,
                temperatures,
            });
        }

        // Antibody titres at the assay visits.
        let first = vaccinations[0].2;
        let last_vax = vaccinations.last().map_or(first, |v| v.2);
        let assays: Vec<(i64, Date)> = t
            .visits
            .iter()
            .filter(|v| {
                schedule
                    .visit(v.visitnum)
                    .and_then(|pv| pv.day)
                    .is_some_and(|day| p.assay_days.contains(&day))
            })
            .map(|v| (v.visitnum, v.date))
            .collect();
        for (visitnum, date) in assays {
            let mut a = Draws::new(stream.label("assay").index(visitnum as u64));
            let noise = a.normal(0.0, 0.1);
            let log = if !responds || date <= first {
                baseline_log
            } else {
                // Doses received before this date (a dose on the assay day comes after the draw).
                let received = vaccinations.iter().filter(|v| v.2 < date).count();
                let fraction = if received >= n_doses {
                    1.0
                } else {
                    p.first_dose_fraction
                };
                let rise = (peak_log - baseline_log).max(0.0) * fraction;
                let waning = if received == n_doses {
                    p.waning_per_year * (f64::from(date.days_since(last_vax)) - 28.0).max(0.0)
                        / 365.25
                } else {
                    0.0
                };
                baseline_log + rise - waning
            };
            let titre = math::pow(10.0, log + noise);
            t.findings.push(Finding {
                domain: FindingDomain::Is,
                testcd: "MBNAB".into(),
                visitnum,
                date,
                value: reported(titre, p.limits),
                laterality: None,
                dp: None,
            });
        }
        t.exposure = exposure;
        t.diaries = diaries;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn titres_are_reported_as_dilutions() {
        let lim = (10.0, 20480.0);
        assert_eq!(reported(9.9, lim), 0.0);
        assert_eq!(reported(10.0, lim), 10.0);
        assert_eq!(reported(39.9, lim), 20.0);
        assert_eq!(reported(40.0, lim), 40.0);
        assert!(reported(30000.0, lim).is_infinite());
    }
}
