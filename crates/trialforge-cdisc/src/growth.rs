//! Growth references for paediatric anthropometry: the CDC 2000 growth charts' LMS parameters
//! for height, weight and BMI for age (2 to 20 years; public domain), as packaged in
//! admiralpeds.
//!
//! [`lms`] returns the parameters for an age in days exactly as admiralpeds' ADVS template
//! prepares them: monthly ages converted to days with `round(AGE * 30.4375)` (R rounds half to
//! even), L, M, S (and the extended-BMI `Sigma` and `P95`) linearly interpolated to every day
//! (`derive_interp_records()`, R `approx()`), restricted to 731–7305 days, and an age matched to
//! its day (ages outside the range take the nearest end).

use std::sync::OnceLock;

use crate::read_tsv;

/// An anthropometric measure with a growth reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Measure {
    /// Height for age (cm).
    Height,
    /// Weight for age (kg).
    Weight,
    /// BMI for age (kg/m²).
    Bmi,
}

/// LMS parameters at one age.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Lms {
    /// Box-Cox power.
    pub l: f64,
    /// Median.
    pub m: f64,
    /// Coefficient of variation.
    pub s: f64,
    /// Extended BMI: scale of the half-normal above the 95th percentile.
    pub sigma: Option<f64>,
    /// Extended BMI: the 95th percentile.
    pub p95: Option<f64>,
}

/// First and last day of the per-day reference.
pub const AGE_DAYS: (i64, i64) = (731, 7305);

/// R's `approx()` (linear, rule 1) at `v`.
fn approx(x: &[f64], y: &[f64], v: f64) -> f64 {
    let (mut i, mut j) = (0, x.len() - 1);
    while i + 1 < j {
        let ij = (i + j) / 2;
        if v < x[ij] {
            j = ij;
        } else {
            i = ij;
        }
    }
    if v == x[j] {
        return y[j];
    }
    if v == x[i] {
        return y[i];
    }
    y[i] + (y[j] - y[i]) * ((v - x[i]) / (x[j] - x[i]))
}

/// Per-day LMS for one measure: index `[sex][day - 731]`, sex 0 = male, 1 = female.
fn table(measure: Measure) -> &'static [Vec<Lms>; 2] {
    static H: OnceLock<[Vec<Lms>; 2]> = OnceLock::new();
    static W: OnceLock<[Vec<Lms>; 2]> = OnceLock::new();
    static B: OnceLock<[Vec<Lms>; 2]> = OnceLock::new();
    let (cell, text) = match measure {
        Measure::Height => (&H, include_str!("../assets/growth/cdc_htage.tsv")),
        Measure::Weight => (&W, include_str!("../assets/growth/cdc_wtage.tsv")),
        Measure::Bmi => (&B, include_str!("../assets/growth/cdc_bmiage.tsv")),
    };
    cell.get_or_init(|| {
        let rows = read_tsv(text);
        let num = |s: &str| s.parse::<f64>().expect("numeric LMS cell");
        let per_sex = |sex: &str| {
            let mut r: Vec<_> = rows.iter().filter(|r| r["SEX"] == sex).collect();
            r.sort_by(|a, b| num(a["AGE"]).total_cmp(&num(b["AGE"])));
            let days: Vec<f64> = r
                .iter()
                .map(|x| (num(x["AGE"]) * 30.4375).round_ties_even())
                .collect();
            let col = |name: &str| -> Option<Vec<f64>> {
                r.first()?.get(name)?;
                Some(r.iter().map(|x| num(x[name])).collect())
            };
            let (l, m, s) = (
                col("L").expect("L"),
                col("M").expect("M"),
                col("S").expect("S"),
            );
            let (sigma, p95) = (col("Sigma"), col("P95"));
            (AGE_DAYS.0..=AGE_DAYS.1)
                .map(|d| {
                    let v = d as f64;
                    Lms {
                        l: approx(&days, &l, v),
                        m: approx(&days, &m, v),
                        s: approx(&days, &s, v),
                        sigma: sigma.as_ref().map(|y| approx(&days, y, v)),
                        p95: p95.as_ref().map(|y| approx(&days, y, v)),
                    }
                })
                .collect::<Vec<Lms>>()
        };
        [per_sex("1"), per_sex("2")]
    })
}

/// The LMS parameters of `measure` for a child of `sex` (`M` or `F`) aged `age_days`.
pub fn lms(measure: Measure, sex: &str, age_days: i64) -> Option<Lms> {
    let k = match sex {
        "M" => 0,
        "F" => 1,
        _ => return None,
    };
    let d = age_days.clamp(AGE_DAYS.0, AGE_DAYS.1);
    table(measure)[k].get((d - AGE_DAYS.0) as usize).copied()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reference_ages_are_interpolated_by_day() {
        // 24.5 months is day round(745.71875) = 746; day 746 is that row exactly.
        let a = lms(Measure::Height, "M", 746).unwrap();
        assert_eq!(a.m, 86.86160934);
        // Between rows the parameters are linear in the day.
        let b = lms(Measure::Height, "M", 760).unwrap();
        let c = lms(Measure::Height, "M", 776).unwrap();
        assert!(a.m < b.m && b.m < c.m);
        assert!(lms(Measure::Bmi, "F", 4000).unwrap().p95.is_some());
        assert!(lms(Measure::Weight, "F", 4000).unwrap().p95.is_none());
        assert_eq!(
            lms(Measure::Height, "M", 100),
            lms(Measure::Height, "M", 731)
        );
    }
}
