//! Ophthalmic examinations of both eyes and intravitreal injections into the study eye,
//! layered on a fixed-duration course (OE findings at every attended visit).
//!
//! Each subject has a study eye (left or right, equally likely) and a fellow eye. Best
//! corrected visual acuity (ETDRS letters, 0–100) follows
//!
//! - study eye: L + (Δ_arm + u)·g(t) + e, with the subject level L ~ N(μ, σ) re-drawn until the
//!   screening value meets the inclusion range, a subject response deviation u ~ N(0, σ_u), the
//!   gain building up as g(t) = (1 − e^(−t/τ)) / (1 − e^(−T/τ)) over the treatment period T,
//!   and test–retest noise e;
//! - fellow eye: a level with a slow decline and the same noise.
//!
//! Central subfield thickness (µm, SD-OCT) has the same structure with its own effect per arm.
//! Values are reported as integers (letters, µm), clamped to their scales.
//!
//! Dosing replaces the course's daily exposure: one injection at each attended visit on a
//! planned injection day of the subject's arm, up to the course's last dose.

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject, VisitKind};

use crate::timeline::{Exposure, Finding, FindingDomain, Timeline};

/// `OETESTCD` of best corrected visual acuity (ETDRS letters correct).
pub const BCVA: &str = "NUMLCOR";
/// `OETESTCD` of central subfield thickness.
pub const CST: &str = "THICK";

/// Ophthalmology parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct OphthaParams {
    /// Study-eye BCVA at baseline: population (mean, SD) in letters.
    pub bcva_baseline: (f64, f64),
    /// Inclusion range for the study eye's screening BCVA (letters).
    pub bcva_inclusion: (f64, f64),
    /// Mean change in study-eye BCVA at the end of treatment, per arm (letters).
    pub bcva_change_by_arm: Vec<f64>,
    /// SD of the subject deviation from the arm's mean BCVA change (letters).
    pub response_sd: f64,
    /// Test–retest SD of a BCVA measurement (letters).
    pub bcva_noise_sd: f64,
    /// Mean change in study-eye central subfield thickness at the end of treatment, per arm
    /// (µm).
    pub cst_change_by_arm: Vec<f64>,
    /// Onset time constant of the treatment effect, in days (τ).
    pub onset_days: f64,
    /// Treatment period length in days (T).
    pub treatment_days: i32,
    /// Planned injection study days, per arm.
    pub injection_days: Vec<Vec<i32>>,
}

/// Simulates ophthalmic examinations and injections onto `timelines` (in subject order).
pub fn simulate(
    design: &Design,
    p: &OphthaParams,
    subjects: &[Subject],
    timelines: &mut [Timeline],
) {
    let schedule = design
        .schedule
        .as_ref()
        .expect("an ophthalmology study has a schedule");
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let stream = design.subject_stream(s.index).label("ophtha");
        let mut d = Draws::new(stream.label("latent"));
        // Draw order is part of the output contract; append new draws at the end.
        let study_eye = if d.uniform() < 0.5 { "LEFT" } else { "RIGHT" };
        let fellow_eye = if study_eye == "LEFT" { "RIGHT" } else { "LEFT" };
        let mut bcva_level = d.normal(p.bcva_baseline.0, p.bcva_baseline.1);
        let fellow_level = d.normal(72.0, 12.0);
        let response = d.normal(0.0, p.response_sd);
        let cst_level = d.normal(420.0, 90.0).max(250.0);
        let fellow_cst = d.normal(275.0, 25.0).max(200.0);
        let cst_response = d.normal(0.0, 40.0);

        // Measurement noise per visit: study BCVA, fellow BCVA, study CST, fellow CST.
        let noise = |visitnum: i64| {
            let mut v = Draws::new(stream.label("visit").index(visitnum as u64));
            [
                v.normal(0.0, p.bcva_noise_sd),
                v.normal(0.0, p.bcva_noise_sd),
                v.normal(0.0, 20.0),
                v.normal(0.0, 10.0),
            ]
        };
        // Inclusion: re-draw the study-eye level until the screening value qualifies.
        let screening = t.visits.iter().find(|v| {
            schedule
                .visit(v.visitnum)
                .is_some_and(|pv| pv.kind == VisitKind::Screening)
        });
        if let Some(scr) = screening {
            let e = noise(scr.visitnum)[0];
            let (lo, hi) = p.bcva_inclusion;
            let mut redraw = Draws::new(stream.label("inclusion"));
            for _ in 0..200 {
                if (lo..=hi).contains(&math::round_dp(bcva_level + e, 0)) {
                    break;
                }
                bcva_level = redraw.normal(p.bcva_baseline.0, p.bcva_baseline.1);
            }
            bcva_level = bcva_level.clamp(lo - e, hi - e);
        }

        let first = t.first_dose();
        let full = 1.0 - math::exp(-f64::from(p.treatment_days) / p.onset_days);
        let ramp = |date: Date| {
            let days = f64::from(date.days_since(first)).max(0.0);
            ((1.0 - math::exp(-days / p.onset_days)) / full).min(1.0)
        };
        let years = |date: Date| f64::from(date.days_since(first)).max(0.0) / 365.25;
        let gain = p.bcva_change_by_arm.get(s.arm).copied().unwrap_or(0.0) + response;
        let thinning = p.cst_change_by_arm.get(s.arm).copied().unwrap_or(0.0) + cst_response;
        let letters = |x: f64| math::round_dp(x, 0).clamp(0.0, 100.0);
        let microns = |x: f64| math::round_dp(x, 0).max(150.0);
        let visits: Vec<(i64, Date)> = t.visits.iter().map(|v| (v.visitnum, v.date)).collect();
        for (visitnum, date) in visits {
            let e = noise(visitnum);
            let g = ramp(date);
            let values = [
                (BCVA, study_eye, letters(bcva_level + gain * g + e[0])),
                (BCVA, fellow_eye, letters(fellow_level - years(date) + e[1])),
                (CST, study_eye, microns(cst_level + thinning * g + e[2])),
                (CST, fellow_eye, microns(fellow_cst + e[3])),
            ];
            for (testcd, eye, value) in values {
                t.findings.push(Finding {
                    domain: FindingDomain::Oe,
                    testcd: testcd.to_string(),
                    visitnum,
                    date,
                    value,
                    laterality: Some(eye),
                    dp: None,
                });
            }
        }

        // One injection per attended visit on a planned injection day, up to the last dose.
        let last = t.last_dose();
        let dose = design.arms[s.arm].treatment.dose;
        let planned = p.injection_days.get(s.arm).cloned().unwrap_or_default();
        let injections: Vec<Exposure> = t
            .visits
            .iter()
            .filter(|v| {
                v.date <= last
                    && schedule
                        .visit(v.visitnum)
                        .and_then(|pv| pv.day)
                        .is_some_and(|day| planned.contains(&day))
            })
            .map(|v| Exposure {
                start: v.date,
                end: v.date,
                dose,
                time: None,
                period: None,
                location: Some(("EYE", study_eye)),
                group: None,
                visitnum: None,
                switched: false,
            })
            .collect();
        if !injections.is_empty() {
            t.exposure = injections;
        }
        t.study_eye = Some(study_eye);
    }
}
