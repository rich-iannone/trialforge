//! Small survival statistics used by the truth manifest and checks: Kaplan–Meier medians and a
//! (stratified) Cox proportional-hazards model for one binary covariate.

use crate::math;

/// One subject's survival observation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Survival {
    /// Follow-up time.
    pub time: f64,
    /// Whether the event happened (otherwise censored).
    pub event: bool,
}

/// The Kaplan–Meier median (the first time survival drops to 0.5 or below), if reached.
pub fn km_median(data: &[Survival]) -> Option<f64> {
    let mut d: Vec<Survival> = data.to_vec();
    d.sort_by(|a, b| a.time.total_cmp(&b.time).then(b.event.cmp(&a.event)));
    let mut at_risk = d.len() as f64;
    let mut surv = 1.0;
    let mut i = 0;
    while i < d.len() {
        let t = d[i].time;
        let (mut events, mut total) = (0.0, 0.0);
        while i < d.len() && d[i].time == t {
            events += f64::from(u8::from(d[i].event));
            total += 1.0;
            i += 1;
        }
        if events > 0.0 {
            surv *= 1.0 - events / at_risk;
            if surv <= 0.5 {
                return Some(t);
            }
        }
        at_risk -= total;
    }
    None
}

/// A Cox model fit for one binary covariate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoxFit {
    /// Log hazard ratio.
    pub beta: f64,
    /// Standard error of `beta`.
    pub se: f64,
}

impl CoxFit {
    /// Hazard ratio.
    pub fn hr(&self) -> f64 {
        math::exp(self.beta)
    }

    /// Wald confidence interval for the hazard ratio at the given level (e.g. 0.95).
    pub fn hr_ci(&self, level: f64) -> (f64, f64) {
        let z = math::norm_inv(0.5 + level / 2.0);
        (
            math::exp(self.beta - z * self.se),
            math::exp(self.beta + z * self.se),
        )
    }
}

/// Fits `h(t) = h_s(t) exp(beta * x)` by partial likelihood with Breslow ties, stratified by
/// `strata` (one baseline hazard per stratum). `x` is 0/1. Returns `None` if the model is not
/// estimable (no events, or no variation in `x` among those at risk at events).
pub fn cox_binary(data: &[Survival], x: &[f64], strata: &[usize]) -> Option<CoxFit> {
    assert!(
        data.len() == x.len() && x.len() == strata.len(),
        "parallel inputs"
    );
    let n_strata = strata.iter().copied().max().map_or(0, |m| m + 1);
    // Per stratum, indices sorted by descending time (events before censorings at ties).
    let mut groups: Vec<Vec<usize>> = vec![Vec::new(); n_strata];
    for (i, &s) in strata.iter().enumerate() {
        groups[s].push(i);
    }
    for g in &mut groups {
        g.sort_by(|&a, &b| {
            data[b]
                .time
                .total_cmp(&data[a].time)
                .then(data[a].event.cmp(&data[b].event))
        });
    }
    let mut beta = 0.0;
    for _ in 0..50 {
        let (mut score, mut info) = (0.0, 0.0);
        for g in &groups {
            // Walk from the longest time down, accumulating risk-set sums; handle tied times.
            let (mut s0, mut s1) = (0.0, 0.0);
            let mut k = 0;
            while k < g.len() {
                let t = data[g[k]].time;
                let mut j = k;
                while j < g.len() && data[g[j]].time == t {
                    let w = math::exp(beta * x[g[j]]);
                    s0 += w;
                    s1 += w * x[g[j]];
                    j += 1;
                }
                let d: Vec<usize> = g[k..j].iter().copied().filter(|&i| data[i].event).collect();
                if !d.is_empty() {
                    let mean = s1 / s0;
                    // x is binary, so x² = x and the risk-set second moment equals the first.
                    let var = mean - mean * mean;
                    for &i in &d {
                        score += x[i] - mean;
                        info += var;
                    }
                }
                k = j;
            }
        }
        if info.is_nan() || info <= 1e-12 {
            return None;
        }
        let step = score / info;
        beta += step.clamp(-5.0, 5.0);
        if step.abs() < 1e-10 {
            return Some(CoxFit {
                beta,
                se: 1.0 / math::sqrt(info),
            });
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(time: f64, event: bool) -> Survival {
        Survival { time, event }
    }

    #[test]
    fn km_median_simple() {
        assert_eq!(
            km_median(&[s(1.0, true), s(2.0, true), s(3.0, false), s(4.0, true)]),
            Some(2.0)
        );
        assert_eq!(km_median(&[s(1.0, false), s(2.0, false)]), None);
    }

    /// Values checked against R: survival::coxph(Surv(t, e) ~ x, ties = "breslow").
    #[test]
    fn cox_matches_r() {
        let t = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0];
        let e = [
            true, true, false, true, true, true, false, true, true, false,
        ];
        let x = [1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0];
        let data: Vec<Survival> = t.iter().zip(e).map(|(&t, e)| s(t, e)).collect();
        let fit = cox_binary(&data, &x, &[0; 10]).unwrap();
        assert!((fit.beta - EXPECTED_BETA).abs() < 1e-6, "beta {}", fit.beta);
        assert!((fit.se - EXPECTED_SE).abs() < 1e-6, "se {}", fit.se);
    }

    const EXPECTED_BETA: f64 = 0.621445757130;
    const EXPECTED_SE: f64 = 0.770453077363;

    #[test]
    fn cox_recovers_a_known_hazard_ratio() {
        use crate::{Draws, Stream};
        let mut d = Draws::new(Stream::root(3).label("cox"));
        let (mut data, mut x, mut strata) = (Vec::new(), Vec::new(), Vec::new());
        for i in 0..20_000 {
            let arm = (i % 2) as f64;
            let stratum = (i / 2) % 2;
            let rate =
                0.1 * if arm == 1.0 { 0.6 } else { 1.0 } * if stratum == 1 { 2.0 } else { 1.0 };
            let t = d.exponential(rate);
            let c = d.uniform_range(0.0, 30.0);
            data.push(s(t.min(c), t <= c));
            x.push(arm);
            strata.push(stratum);
        }
        let fit = cox_binary(&data, &x, &strata).unwrap();
        assert!((fit.hr() - 0.6).abs() < 0.03, "hr {}", fit.hr());
        let (lo, hi) = fit.hr_ci(0.95);
        assert!(lo < 0.6 && 0.6 < hi);
    }
}
