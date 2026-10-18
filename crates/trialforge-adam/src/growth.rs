//! Paediatric anthropometric indicators for ADVS, following admiralpeds'
//! `derive_params_growth_age()` with the CDC 2000 growth charts (children aged 2 years or more):
//!
//! - z-score `((X / M)^L - 1) / (L × S)` and percentile `Φ(z) × 100`, with L, M, S for the
//!   child's sex and current age in days (`AAGECUR`, admiral `derive_vars_duration()`: the days
//!   from birth to `ADT`, plus one);
//! - BMI at or above the 95th percentile uses the CDC extended percentiles:
//!   `90 + 10 × Φ((X - P95) / Sigma)`, and the z-score is its normal quantile (8.21 if infinite).
//!
//! Adaptation: only the CDC reference is used (the archetype enrols children from 3 years), so
//! the WHO charts for children under 2 are not needed.

use trialforge_cdisc::growth::{Measure, lms};
use trialforge_core::math;

/// (VS test, z-score PARAMCD and PARAM, percentile PARAMCD and PARAM).
pub(crate) type Indicator = (&'static str, Measure, [(&'static str, &'static str); 2]);

pub(crate) const INDICATORS: [Indicator; 3] = [
    (
        "WEIGHT",
        Measure::Weight,
        [
            ("WGTASDS", "Weight-for-age z-score"),
            ("WGTAPCTL", "Weight-for-age percentile"),
        ],
    ),
    (
        "HEIGHT",
        Measure::Height,
        [
            ("HGTSDS", "Height-for-age z-score"),
            ("HGTPCTL", "Height-for-age percentile"),
        ],
    ),
    (
        "BMI",
        Measure::Bmi,
        [
            ("BMISDS", "BMI-for-age z-score"),
            ("BMIPCTL", "BMI-for-age percentile"),
        ],
    ),
];

/// Annualized height velocity: PARAMCD and PARAM.
pub(crate) const HTVEL: (&str, &str) = ("HTVEL", "Annualized Height Velocity (cm/yr)");

/// (z-score, percentile) of `x` for a child of `sex` aged `age_days`.
pub(crate) fn indicator(measure: Measure, sex: &str, age_days: i64, x: f64) -> Option<(f64, f64)> {
    let p = lms(measure, sex, age_days)?;
    let z = (math::pow(x / p.m, p.l) - 1.0) / (p.l * p.s);
    let mut pctl = math::norm_cdf(z) * 100.0;
    let mut sds = z;
    if measure == Measure::Bmi
        && let (Some(p95), Some(sigma)) = (p.p95, p.sigma)
        && x >= p95
    {
        pctl = 90.0 + 10.0 * math::norm_cdf((x - p95) / sigma);
        sds = if pctl >= 100.0 {
            8.21
        } else {
            math::norm_inv(pctl / 100.0)
        };
    }
    Some((sds, pctl))
}
