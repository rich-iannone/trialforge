//! Dose escalation: the sequence of cohorts a phase 1 dose-finding design treats, and the
//! dose-limiting toxicities (DLTs) layered on each subject's fixed-duration course.
//!
//! Each subject has a latent tolerance u ~ U(0, 1) on their own stream; at dose level d they
//! have a DLT in the observation window when u < p_d, the true DLT probability of the level.
//! The design decides the next cohort's level from the DLTs seen so far:
//!
//! - **3+3**: 0/3 escalates; 1/3 adds 3 at the same level, then ≤ 1/6 escalates; ≥ 2 DLTs
//!   de-escalate, expanding the lower level to 6 if it had 3; the MTD is the highest level
//!   with ≤ 1/6 (the top level after ≤ 1/6 there).
//! - **BOIN** (Liu & Yuan 2015): escalate when the level's DLT rate is ≤ λe, de-escalate when
//!   ≥ λd (φ1 = 0.6φ, φ2 = 1.4φ); levels with P(p > φ | data) > 0.95 (Beta(1, 1) prior, at
//!   least 3 treated) are eliminated with all higher ones.
//! - **mTPI-2** (Guo et al. 2017): the decision of the interval with the largest unit
//!   probability mass, with the equivalence interval [φ − 0.05, φ + 0.05] and the rest of
//!   [0, 1] cut into intervals of the same width; same elimination rule.
//!
//! BOIN and mTPI-2 treat cohorts until the maximum sample size and select the MTD by
//! isotonic regression of (x + 0.05) / (n + 0.1) over the untreated-free, non-eliminated
//! levels: the level closest to φ (on ties, the higher if below φ, else the lower).

use trialforge_core::{Date, Draws, math};
use trialforge_model::{Design, Subject};

use crate::timeline::{AdverseEvent, SeriousCriteria, Severity, Timeline};

/// Days of the DLT observation window (cycle 1).
pub const DLT_WINDOW_DAYS: i32 = 28;

/// The dose-finding design.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Method {
    /// The traditional 3+3.
    ThreeThree,
    /// Bayesian optimal interval.
    Boin,
    /// Modified toxicity probability interval-2.
    Mtpi2,
}

/// Escalation parameters.
#[derive(Clone, Debug, PartialEq)]
pub struct EscalationParams {
    /// The design.
    pub method: Method,
    /// True DLT probability per dose level.
    pub dlt_rates: Vec<f64>,
    /// Target DLT rate φ (BOIN, mTPI-2).
    pub target: f64,
    /// Subjects per cohort (BOIN, mTPI-2; the 3+3 always treats 3).
    pub cohort_size: usize,
    /// Maximum number of subjects.
    pub max_n: usize,
}

/// What the escalation did.
#[derive(Clone, Debug, PartialEq)]
pub struct Trace {
    /// Treated cohorts in order: (dose level, subjects).
    pub cohorts: Vec<(usize, usize)>,
    /// The selected MTD level; `None` if every level was too toxic.
    pub mtd: Option<usize>,
}

/// P(X ≥ k) for X ~ Binomial(n, p).
fn binom_tail(n: u32, k: u32, p: f64) -> f64 {
    let mut coef = 1.0;
    let mut total = 0.0;
    for j in 0..=n {
        if j > 0 {
            coef = coef * f64::from(n - j + 1) / f64::from(j);
        }
        if j >= k {
            total += coef * math::pow(p, f64::from(j)) * math::pow(1.0 - p, f64::from(n - j));
        }
    }
    total
}

/// The regularized incomplete beta I_t(a, b) for integer a, b ≥ 1: P(Bin(a + b − 1, t) ≥ a).
fn beta_cdf(t: f64, a: u32, b: u32) -> f64 {
    if t <= 0.0 {
        0.0
    } else if t >= 1.0 {
        1.0
    } else {
        binom_tail(a + b - 1, a, t)
    }
}

/// P(p > φ | x DLTs in n) under a Beta(1, 1) prior.
pub fn overdose_probability(x: u32, n: u32, phi: f64) -> f64 {
    1.0 - beta_cdf(phi, 1 + x, 1 + n - x)
}

/// BOIN's escalation and de-escalation boundaries (λe, λd).
pub fn boin_boundaries(phi: f64) -> (f64, f64) {
    let (p1, p2) = (0.6 * phi, 1.4 * phi);
    let e = math::ln((1.0 - p1) / (1.0 - phi)) / math::ln(phi * (1.0 - p1) / (p1 * (1.0 - phi)));
    let d = math::ln((1.0 - phi) / (1.0 - p2)) / math::ln(p2 * (1.0 - phi) / (phi * (1.0 - p2)));
    (e, d)
}

/// The mTPI-2 decision: -1 de-escalate, 0 stay, 1 escalate.
fn mtpi2(x: u32, n: u32, phi: f64) -> i32 {
    let eps = 0.05;
    let (lo, hi) = (phi - eps, phi + eps);
    let w = 2.0 * eps;
    let mass = |a: f64, b: f64| beta_cdf(b, 1 + x, 1 + n - x) - beta_cdf(a, 1 + x, 1 + n - x);
    let mut best = (mass(lo, hi) / w, 0);
    let mut b = lo;
    while b > 0.0 {
        let a = (b - w).max(0.0);
        let upm = mass(a, b) / (b - a);
        if upm > best.0 {
            best = (upm, 1);
        }
        b = a;
    }
    let mut a = hi;
    while a < 1.0 {
        let b = (a + w).min(1.0);
        let upm = mass(a, b) / (b - a);
        if upm > best.0 {
            best = (upm, -1);
        }
        a = b;
    }
    best.1
}

/// Isotonic (pool-adjacent-violators) estimates of `rates` with `weights`.
fn isotonic(rates: &[f64], weights: &[f64]) -> Vec<f64> {
    let mut blocks: Vec<(f64, f64, usize)> = Vec::new(); // (value, weight, count)
    for (&r, &w) in rates.iter().zip(weights) {
        blocks.push((r, w, 1));
        while blocks.len() > 1 && blocks[blocks.len() - 2].0 > blocks[blocks.len() - 1].0 {
            let (v2, w2, c2) = blocks.pop().expect("two blocks");
            let (v1, w1, c1) = blocks.pop().expect("two blocks");
            blocks.push(((v1 * w1 + v2 * w2) / (w1 + w2), w1 + w2, c1 + c2));
        }
    }
    blocks
        .into_iter()
        .flat_map(|(v, _, c)| std::iter::repeat_n(v, c))
        .collect()
}

/// Model-based MTD selection over treated, non-eliminated levels.
fn select(n: &[u32], x: &[u32], allowed: usize, phi: f64) -> Option<usize> {
    let levels: Vec<usize> = (0..allowed).filter(|&l| n[l] > 0).collect();
    if levels.is_empty() {
        return None;
    }
    let rates: Vec<f64> = levels
        .iter()
        .map(|&l| (f64::from(x[l]) + 0.05) / (f64::from(n[l]) + 0.1))
        .collect();
    let weights: Vec<f64> = levels.iter().map(|&l| f64::from(n[l]) + 0.1).collect();
    let est = isotonic(&rates, &weights);
    let mut best = 0;
    for k in 1..levels.len() {
        let (d, bd) = ((est[k] - phi).abs(), (est[best] - phi).abs());
        // Ties: the higher level when the estimate is below the target, else the lower.
        if d < bd - 1e-12 || ((d - bd).abs() <= 1e-12 && est[k] < phi) {
            best = k;
        }
    }
    Some(levels[best])
}

/// Escalation state: subjects and DLTs per level, cohorts treated so far.
struct State {
    n: Vec<u32>,
    x: Vec<u32>,
    cohorts: Vec<(usize, usize)>,
    treated: usize,
}

impl State {
    fn treat(&mut self, level: usize, size: usize, dlt: &impl Fn(usize, usize) -> bool) {
        let dlts = (0..size).filter(|k| dlt(self.treated + k, level)).count() as u32;
        self.treated += size;
        self.n[level] += size as u32;
        self.x[level] += dlts;
        self.cohorts.push((level, size));
    }
}

/// Runs the design. `dlt(i)` says whether the i-th subject treated (0-based) at `level` has a
/// DLT.
pub fn escalate(p: &EscalationParams, dlt: impl Fn(usize, usize) -> bool) -> Trace {
    let top = p.dlt_rates.len();
    let mut st = State {
        n: vec![0; top],
        x: vec![0; top],
        cohorts: Vec::new(),
        treated: 0,
    };
    let mtd = match p.method {
        Method::ThreeThree => {
            let mut level = 0;
            let mut mtd;
            loop {
                if st.treated + 3 > p.max_n {
                    mtd = (0..top).rev().find(|&l| st.n[l] >= 6 && st.x[l] <= 1);
                    break;
                }
                st.treat(level, 3, &dlt);
                let tolerated = (st.n[level] == 3 && st.x[level] == 0)
                    || (st.n[level] == 6 && st.x[level] <= 1);
                if st.n[level] == 3 && st.x[level] == 1
                    || (level + 1 == top && st.n[level] == 3 && st.x[level] == 0)
                {
                    continue; // expand to 6 at this level (the top level is always confirmed on 6)
                }
                if tolerated {
                    if level + 1 == top {
                        mtd = Some(level);
                        break;
                    }
                    level += 1;
                    continue;
                }
                // Too toxic: de-escalate, expanding a lower level to 6 when it has 3.
                let mut l = level;
                loop {
                    if l == 0 {
                        mtd = None;
                        break;
                    }
                    l -= 1;
                    if st.n[l] >= 6 {
                        mtd = (st.x[l] <= 1).then_some(l);
                        if mtd.is_some() {
                            break;
                        }
                        continue;
                    }
                    if st.treated + 3 > p.max_n {
                        mtd = None;
                        break;
                    }
                    st.treat(l, 3, &dlt);
                    if st.x[l] <= 1 {
                        mtd = Some(l);
                        break;
                    }
                }
                break;
            }
            mtd
        }
        Method::Boin | Method::Mtpi2 => {
            let (le, ld) = boin_boundaries(p.target);
            let mut allowed = top;
            let mut level = 0;
            while st.treated < p.max_n && allowed > 0 {
                let size = p.cohort_size.min(p.max_n - st.treated);
                st.treat(level, size, &dlt);
                if st.n[level] >= 3
                    && overdose_probability(st.x[level], st.n[level], p.target) > 0.95
                {
                    allowed = level;
                    if level == 0 {
                        break;
                    }
                    level -= 1;
                    continue;
                }
                let decision = match p.method {
                    Method::Boin => {
                        let rate = f64::from(st.x[level]) / f64::from(st.n[level]);
                        if rate <= le {
                            1
                        } else if rate >= ld {
                            -1
                        } else {
                            0
                        }
                    }
                    _ => mtpi2(st.x[level], st.n[level], p.target),
                };
                if decision == 1 && level + 1 < allowed {
                    level += 1;
                } else if decision == -1 && level > 0 {
                    level -= 1;
                }
            }
            select(&st.n, &st.x, allowed, p.target)
        }
    };
    Trace {
        cohorts: st.cohorts,
        mtd,
    }
}

/// Probabilities of selecting each level as the MTD (last entry: no MTD) and the mean sample
/// size over `trials` simulated trials with fresh DLT draws from `draws`.
pub fn operating_characteristics(
    p: &EscalationParams,
    trials: usize,
    draws: &mut Draws,
) -> (Vec<f64>, f64) {
    let mut picks = vec![0usize; p.dlt_rates.len() + 1];
    let mut total = 0usize;
    for _ in 0..trials {
        let u: Vec<f64> = (0..p.max_n.max(3 * p.dlt_rates.len() * 2))
            .map(|_| draws.uniform())
            .collect();
        let t = escalate(p, |i, level| u[i] < p.dlt_rates[level]);
        picks[t.mtd.unwrap_or(p.dlt_rates.len())] += 1;
        total += t.cohorts.iter().map(|c| c.1).sum::<usize>();
    }
    (
        picks.iter().map(|&k| k as f64 / trials as f64).collect(),
        total as f64 / trials as f64,
    )
}

/// DLT PTs and weights.
const DLT_TERMS: [(&str, f64); 7] = [
    ("Neutropenia", 0.30),
    ("Thrombocytopenia", 0.20),
    ("Alanine aminotransferase increased", 0.15),
    ("Fatigue", 0.10),
    ("Stomatitis", 0.10),
    ("Diarrhoea", 0.10),
    ("Rash", 0.05),
];

/// Adds each DLT as a severe AE in the subject's DLT window and flags every AE starting in the
/// window (`dlt`). `level` gives each subject's dose level (by arm).
pub fn simulate(
    design: &Design,
    p: &EscalationParams,
    subjects: &[Subject],
    timelines: &mut [Timeline],
) {
    for (s, t) in subjects.iter().zip(timelines.iter_mut()) {
        let level = s.arm;
        let has_dlt = tolerance(design, s.index) < p.dlt_rates[level];
        let first = t.first_dose();
        let window_end = first.add_days(DLT_WINDOW_DAYS - 1);
        for a in &mut t.adverse_events {
            if a.start >= first && a.start <= window_end {
                a.dlt = Some(false);
            }
        }
        if has_dlt {
            let mut d = Draws::new(design.subject_stream(s.index).label("dlt-event"));
            let w: Vec<f64> = DLT_TERMS.iter().map(|x| x.1).collect();
            let pt = DLT_TERMS[d.weighted(&w)].0;
            let day = d.int_range(2, i64::from(DLT_WINDOW_DAYS)) as i32;
            let serious = d.uniform() < 0.3;
            let duration = d.int_range(4, 14) as i32;
            let start: Date = first.add_days(day - 1).min(t.last_contact).max(first);
            let end = start.add_days(duration - 1);
            t.adverse_events.push(AdverseEvent {
                pt: pt.to_string(),
                start,
                end: (end <= t.last_contact).then_some(end),
                severity: Severity::Severe,
                serious: serious.then_some(SeriousCriteria {
                    hospitalization: true,
                    ..SeriousCriteria::default()
                }),
                causality: "RELATED",
                action: if start < t.last_dose() {
                    "DRUG INTERRUPTED"
                } else {
                    "DOSE NOT CHANGED"
                },
                outcome: if end <= t.last_contact {
                    "RECOVERED/RESOLVED"
                } else {
                    "NOT RECOVERED/NOT RESOLVED"
                },
                treated: false,
                dlt: Some(true),
            });
            t.adverse_events
                .sort_by(|a, b| (a.start, &a.pt).cmp(&(b.start, &b.pt)));
        }
    }
}

/// Subject `index`'s latent DLT tolerance.
pub fn tolerance(design: &Design, index: usize) -> f64 {
    Draws::new(design.subject_stream(index).label("dlt")).uniform()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn boin_boundaries_match_the_published_table() {
        let (e, d) = boin_boundaries(0.3);
        assert!(
            (e - 0.236).abs() < 1e-3 && (d - 0.358).abs() < 1e-3,
            "{e} {d}"
        );
        let (e, d) = boin_boundaries(0.25);
        assert!(
            (e - 0.197).abs() < 1e-3 && (d - 0.298).abs() < 1e-3,
            "{e} {d}"
        );
    }

    #[test]
    fn three_plus_three_follows_its_rules() {
        let p = EscalationParams {
            method: Method::ThreeThree,
            dlt_rates: vec![0.0, 0.0, 1.0],
            target: 0.25,
            cohort_size: 3,
            max_n: 60,
        };
        let t = escalate(&p, |_, level| level == 2);
        assert_eq!(t.cohorts, vec![(0, 3), (1, 3), (2, 3), (1, 3)]);
        assert_eq!(t.mtd, Some(1));
    }

    #[test]
    fn interval_designs_find_the_target_dose() {
        for method in [Method::Boin, Method::Mtpi2] {
            let p = EscalationParams {
                method,
                dlt_rates: vec![0.05, 0.1, 0.25, 0.45, 0.6],
                target: 0.25,
                cohort_size: 3,
                max_n: 36,
            };
            let mut d = Draws::new(trialforge_core::Stream::root(1).label("oc"));
            let (pick, _) = operating_characteristics(&p, 400, &mut d);
            assert!(pick[2] > 0.4, "{method:?}: {pick:?}");
        }
    }

    #[test]
    fn overdose_probability_is_the_beta_tail() {
        // Beta(1 + 3, 1 + 0) = Beta(4, 1): P(p > 0.3) = 1 − 0.3^4.
        assert!((overdose_probability(3, 3, 0.3) - (1.0 - 0.0081)).abs() < 1e-12);
    }
}
