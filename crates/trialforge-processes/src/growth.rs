//! Paediatric growth, layered on a fixed-duration course: a date of birth for each child and
//! height, weight and BMI at every visit with vital signs, replacing the course's adult
//! anthropometry.
//!
//! - **Birth date**: uniform over the year that makes the child's completed age at consent
//!   equal to `AGE`.
//! - **Height**: the child's baseline height-for-age z-score z₀ (CDC 2000 charts) is drawn and
//!   re-drawn until the screening height meets the inclusion maximum; height at first dose is
//!   the chart's height for z₀ at that age. On treatment the child grows at an annualized
//!   velocity v ~ N(v_arm, σ_v) (at least 1 cm/yr); before treatment at the untreated
//!   velocity. Measurements add stadiometer noise and are reported to 0.1 cm.
//! - **Weight**: the chart's BMI for the child's own BMI-for-age z-score, times height², with
//!   noise; BMI is reported from the reported weight and height.

use trialforge_cdisc::growth::{Measure, lms};
use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject, VisitKind};

use crate::timeline::{Finding, FindingDomain, Timeline};

/// Growth parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct GrowthParams {
    /// Baseline height-for-age z-score: population (mean, SD).
    pub height_sds: (f64, f64),
    /// Inclusion: the screening height-for-age z-score must be at most this.
    pub max_screening_sds: f64,
    /// Mean annualized height velocity on treatment, per arm (cm/yr).
    pub velocity_by_arm: Vec<f64>,
    /// SD of children's height velocity around their arm's mean (cm/yr).
    pub velocity_sd: f64,
    /// Height velocity before treatment (cm/yr).
    pub untreated_velocity: f64,
}

/// Days of age (admiral `derive_vars_duration()`: the end day counts).
pub fn age_days(birth: Date, date: Date) -> i64 {
    i64::from(date.days_since(birth)) + 1
}

/// The chart's value at z-score `z`: M (1 + L S z)^(1/L).
fn at_z(measure: Measure, sex: &str, age: i64, z: f64) -> f64 {
    let p = lms(measure, sex, age).expect("ages are within the CDC charts");
    p.m * math::pow(1.0 + p.l * p.s * z, 1.0 / p.l)
}

/// The z-score of `x`: ((x / M)^L − 1) / (L S).
fn z_of(measure: Measure, sex: &str, age: i64, x: f64) -> f64 {
    let p = lms(measure, sex, age).expect("ages are within the CDC charts");
    (math::pow(x / p.m, p.l) - 1.0) / (p.l * p.s)
}

/// A birth date that makes the completed age at `consent` equal to `age` years: `offset`
/// (0–364) days before the latest such date.
fn birth_date(consent: Date, age: i64, offset: i64) -> Date {
    let (y, m, d) = consent.ymd();
    let year = y - age as i32;
    // The latest birth date: the same day `age` years earlier (28 February for 29 February).
    let latest = Date::from_ymd(year, m, d)
        .or_else(|| Date::from_ymd(year, m, d - 1))
        .expect("valid date");
    let earliest = Date::from_ymd(year - 1, m, d)
        .or_else(|| Date::from_ymd(year - 1, m, d - 1))
        .expect("valid date")
        .add_days(1);
    let span = i64::from(latest.days_since(earliest));
    latest.add_days(-(offset.min(span) as i32))
}

/// Adds birth dates and growth measurements to `timelines` (in subject order).
pub fn simulate(
    design: &Design,
    p: &GrowthParams,
    subjects: &[Subject],
    timelines: &mut [Timeline],
) {
    let schedule = design
        .schedule
        .as_ref()
        .expect("a growth study has a schedule");
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let sex = design.population.sex.levels[s.sex].0.as_str();
        let stream = design.subject_stream(s.index).label("growth");
        let mut d = Draws::new(stream.label("latent"));
        // Draw order is part of the output contract; append new draws at the end.
        let birth = birth_date(s.consent, s.age, d.int_range(0, 364));
        let mut z0 = d.normal(p.height_sds.0, p.height_sds.1);
        let velocity = d
            .normal(
                p.velocity_by_arm.get(s.arm).copied().unwrap_or(0.0),
                p.velocity_sd,
            )
            .max(1.0);
        let bmi_z = d.normal(0.0, 0.8).clamp(-2.5, 2.5);

        let first = t.first_dose();
        let years = |date: Date| f64::from(date.days_since(first)) / 365.25;
        let noise = |visitnum: i64| {
            let mut v = Draws::new(stream.label("visit").index(visitnum as u64));
            (v.normal(0.0, 0.2), v.normal(0.0, 0.015))
        };
        let height_at = |z0: f64, date: Date, e: f64| {
            let h0 = at_z(Measure::Height, sex, age_days(birth, first), z0);
            let grown = if date >= first {
                velocity * years(date)
            } else {
                p.untreated_velocity * years(date)
            };
            math::round_dp(h0 + grown + e, 1)
        };
        // Inclusion on the screening height.
        let screening = t.visits.iter().find(|v| {
            schedule
                .visit(v.visitnum)
                .is_some_and(|pv| pv.kind == VisitKind::Screening)
        });
        if let Some(scr) = screening {
            let e = noise(scr.visitnum).0;
            let mut redraw = Draws::new(stream.label("inclusion"));
            for _ in 0..200 {
                let h = height_at(z0, scr.date, e);
                if z_of(Measure::Height, sex, age_days(birth, scr.date), h) <= p.max_screening_sds {
                    break;
                }
                z0 = redraw.normal(p.height_sds.0, p.height_sds.1);
            }
        }

        t.findings.retain(|f| {
            !(f.domain == FindingDomain::Vs && matches!(f.testcd.as_str(), "HEIGHT" | "WEIGHT"))
        });
        let measured: Vec<(i64, Date)> = t
            .visits
            .iter()
            .filter(|v| schedule.visit(v.visitnum).is_some_and(|pv| pv.vitals))
            .map(|v| (v.visitnum, v.date))
            .collect();
        for (visitnum, date) in measured {
            let (eh, ew) = noise(visitnum);
            let height = height_at(z0, date, eh);
            let bmi = at_z(Measure::Bmi, sex, age_days(birth, date), bmi_z);
            let weight = math::round_dp(bmi * (height / 100.0) * (height / 100.0) * (1.0 + ew), 1);
            let reported_bmi = math::round_dp(weight / ((height / 100.0) * (height / 100.0)), 1);
            for (testcd, value) in [
                ("HEIGHT", height),
                ("WEIGHT", weight),
                ("BMI", reported_bmi),
            ] {
                t.findings.push(Finding {
                    domain: FindingDomain::Vs,
                    testcd: testcd.into(),
                    visitnum,
                    date,
                    value,
                    laterality: None,
                    dp: Some(1),
                });
            }
        }
        t.birth_date = Some(birth);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn completed_years(birth: Date, on: Date) -> i64 {
        let ((by, bm, bd), (y, m, d)) = (birth.ymd(), on.ymd());
        i64::from(y - by) - i64::from((m, d) < (bm, bd))
    }

    #[test]
    fn birth_dates_give_the_recorded_age() {
        for consent in ["2024-02-29", "2024-03-01", "2023-12-31", "2024-01-01"] {
            let c = Date::parse_iso(consent).unwrap();
            for offset in [0, 1, 180, 364] {
                let b = birth_date(c, 7, offset);
                assert_eq!(completed_years(b, c), 7, "{consent} {offset} {b}");
            }
        }
    }

    #[test]
    fn chart_values_invert() {
        let h = at_z(Measure::Height, "M", 2500, -2.5);
        assert!((z_of(Measure::Height, "M", 2500, h) + 2.5).abs() < 1e-9);
        assert!(h > 95.0 && h < 115.0, "{h}");
    }
}
