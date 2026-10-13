//! Calibration: maps clinical targets (medians, hazard ratios, confirmed response rates) to the
//! latent parameters of the process engines.
//!
//! - [`median_scale`]: the Weibull scale that gives a target median under a mixture of strata
//!   with proportional hazard multipliers (solved exactly by bisection on the mixture survival).
//! - [`oncology`]: PFS-median and confirmed-ORR calibration for the oncology course, by Monte
//!   Carlo with common random numbers over a large synthetic cohort.

pub mod oncology;

/// Bisection for an increasing function: the `x` in `[lo, hi]` with `f(x) ≈ target`.
pub fn bisect(
    f: impl Fn(f64) -> f64,
    mut lo: f64,
    mut hi: f64,
    target: f64,
    iterations: u32,
) -> f64 {
    for _ in 0..iterations {
        let mid = 0.5 * (lo + hi);
        if f(mid) < target {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// The Weibull scale (same units as `median`) such that the population median is `median`, when
/// a share `weights[s]` of the population has hazard multiplier `hrs[s]`:
/// `Σ_s w_s exp(-hrs[s] (t / scale)^shape) = 0.5` at `t = median`.
pub fn median_scale(median: f64, shape: f64, weights: &[f64], hrs: &[f64]) -> f64 {
    let total: f64 = weights.iter().sum();
    let surv = |u: f64| -> f64 {
        weights
            .iter()
            .zip(hrs)
            .map(|(w, h)| w / total * trialforge_core::math::exp(-h * u))
            .sum()
    };
    // surv is decreasing in u = (median / scale)^shape; find surv(u) = 0.5.
    let u = bisect(|u| 1.0 - surv(u), 0.0, 100.0, 0.5, 200);
    median / trialforge_core::math::pow(u, 1.0 / shape)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn median_scale_single_stratum_is_closed_form() {
        // S(m) = exp(-(m / scale)^k) = 0.5 → scale = m / ln(2)^(1/k).
        let s = median_scale(10.0, 1.5, &[1.0], &[1.0]);
        assert!((s - 10.0 / 2f64.ln().powf(1.0 / 1.5)).abs() < 1e-9);
    }

    #[test]
    fn median_scale_mixture_hits_the_median() {
        let (w, h, k) = ([0.45, 0.55], [0.8, 1.25], 1.2);
        let s = median_scale(15.0, k, &w, &h);
        let surv: f64 = w
            .iter()
            .zip(h)
            .map(|(w, h)| w * (-h * (15.0f64 / s).powf(k)).exp())
            .sum();
        assert!((surv - 0.5).abs() < 1e-9);
    }
}
