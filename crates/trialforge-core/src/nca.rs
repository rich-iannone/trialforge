//! Non-compartmental analysis of one concentration–time profile, following PKNCA's defaults
//! (PKNCA 0.12) so that results can be reproduced there:
//!
//! - BLQ values are zeros; leading and trailing zeros are kept and zeros between the first and
//!   last measurable concentration are dropped (`conc.blq`).
//! - AUCs use the "lin up/log down" trapezoidal rule; AUClast stops at the last measurable
//!   concentration, AUCall continues linearly to the next (zero) sample, and AUCinf.obs adds
//!   Clast.obs / λz.
//! - λz: log-linear fits to the last 3, 4, … points after Tmax (positive concentrations only);
//!   among fits with λz > 0 and adjusted R² within 1e-4 of the best, the one with the most
//!   points.
//! - Tmax is the first time of Cmax.
//!
//! Times are hours after the reference dose, increasing; concentrations are ≥ 0.

use crate::math;

/// The terminal phase fit.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct HalfLife {
    /// λz (1/h).
    pub lambda_z: f64,
    /// R² of the fit.
    pub r_squared: f64,
    /// Adjusted R².
    pub adj_r_squared: f64,
    /// Number of points.
    pub n_points: usize,
    /// First time in the fit.
    pub time_first: f64,
    /// Last time in the fit.
    pub time_last: f64,
    /// ln 2 / λz.
    pub half_life: f64,
}

/// Which AUC.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AucType {
    /// To the last measurable concentration.
    Last,
    /// To the last time, including trailing BLQ.
    All,
    /// To infinity, extrapolating Clast.obs / λz.
    Inf(f64),
}

/// A cleaned profile: (time, conc) after the BLQ rules.
#[derive(Clone, Debug, PartialEq)]
pub struct Profile {
    /// Times (h).
    pub time: Vec<f64>,
    /// Concentrations.
    pub conc: Vec<f64>,
}

impl Profile {
    /// Applies PKNCA's default BLQ rules to raw (time, conc) pairs (BLQ given as 0).
    pub fn new(time: &[f64], conc: &[f64]) -> Self {
        let mut p = Profile {
            time: time.to_vec(),
            conc: conc.to_vec(),
        };
        if p.time.is_empty() {
            return p;
        }
        let (tfirst, tlast) = if p.conc.iter().all(|c| *c == 0.0) {
            let t = p.time.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            (t, t + 1.0)
        } else {
            (
                tfirst(&p.time, &p.conc).expect("measurable"),
                tlast(&p.time, &p.conc).expect("measurable"),
            )
        };
        let keep: Vec<bool> = p
            .time
            .iter()
            .zip(&p.conc)
            .map(|(t, c)| !(*c == 0.0 && tfirst < *t && *t < tlast))
            .collect();
        let mut k = keep.iter();
        p.time.retain(|_| *k.next().expect("same length"));
        let mut k = keep.iter();
        p.conc.retain(|_| *k.next().expect("same length"));
        p
    }

    /// Cmax (`None` if every concentration is zero).
    pub fn cmax(&self) -> Option<f64> {
        self.conc
            .iter()
            .copied()
            .filter(|c| *c > 0.0)
            .fold(None, |m, c| Some(m.map_or(c, |m: f64| m.max(c))))
    }

    /// First time of Cmax.
    pub fn tmax(&self) -> Option<f64> {
        let m = self.cmax()?;
        self.time
            .iter()
            .zip(&self.conc)
            .find(|(_, c)| **c == m)
            .map(|(t, _)| *t)
    }

    /// Time of the last measurable concentration.
    pub fn tlast(&self) -> Option<f64> {
        tlast(&self.time, &self.conc)
    }

    /// Last measurable concentration (0 if none).
    pub fn clast(&self) -> f64 {
        self.tlast()
            .and_then(|t| self.time.iter().position(|x| *x == t))
            .map_or(0.0, |i| self.conc[i])
    }

    /// The λz fit, if at least 3 positive points follow Tmax and a fit has λz > 0.
    pub fn half_life(&self) -> Option<HalfLife> {
        let tmax = self.tmax()?;
        let pts: Vec<(f64, f64)> = self
            .time
            .iter()
            .zip(&self.conc)
            .filter(|(t, c)| **c > 0.0 && **t > tmax)
            .map(|(t, c)| (*t, math::ln(*c)))
            .collect();
        if pts.len() < 3 {
            return None;
        }
        // Candidates use the last i points (i = 3..=n).
        let fits: Vec<HalfLife> = (3..=pts.len())
            .map(|i| fit(&pts[pts.len() - i..]))
            .collect();
        let best = fits
            .iter()
            .map(|f| f.adj_r_squared)
            .filter(|x| !x.is_nan())
            .fold(f64::NEG_INFINITY, f64::max);
        let ok: Vec<&HalfLife> = fits
            .iter()
            .filter(|f| f.lambda_z > 0.0 && f.adj_r_squared > best - 1e-4)
            .collect();
        ok.into_iter().max_by_key(|f| f.n_points).copied()
    }

    /// AUC over `[start, end]` (`end` may be infinite), PKNCA's `pk.calc.auc`.
    pub fn auc(&self, start: f64, end: f64, kind: AucType) -> Option<f64> {
        let n = self.time.len();
        if n < 2 {
            return None;
        }
        if self.conc.iter().all(|c| *c == 0.0) {
            return Some(0.0);
        }
        if start < self.time[0] {
            return None;
        }
        let mut time = Vec::with_capacity(n + 2);
        let mut conc = Vec::with_capacity(n + 2);
        time.push(start);
        conc.push(self.conc_at(start, kind)?);
        for (t, c) in self.time.iter().zip(&self.conc) {
            if *t != start {
                time.push(*t);
                conc.push(*c);
            }
        }
        if end.is_finite() {
            if let Some(c_end) = self.conc_at(end, kind) {
                let keep: Vec<bool> = time.iter().map(|t| *t != end).collect();
                let (mut t2, mut c2) = (Vec::new(), Vec::new());
                for i in 0..time.len() {
                    if keep[i] {
                        t2.push(time[i]);
                        c2.push(conc[i]);
                    }
                }
                t2.push(end);
                c2.push(c_end);
                time = t2;
                conc = c2;
            }
        }
        let (time, conc): (Vec<f64>, Vec<f64>) = time
            .into_iter()
            .zip(conc)
            .filter(|(t, _)| start <= *t && *t <= end)
            .unzip();
        let Some(tl) = tlast(&time, &conc) else {
            return Some(0.0);
        };
        let methods = interval_methods(&conc, &time, tl, kind);
        let mut total = 0.0;
        for i in 0..time.len() - 1 {
            let (c1, c2, t1, t2) = (conc[i], conc[i + 1], time[i], time[i + 1]);
            total += match methods[i] {
                Method::Linear => (c1 + c2) / 2.0 * (t2 - t1),
                Method::Log => (c1 - c2) / math::ln(c1 / c2) * (t2 - t1),
                Method::Zero => 0.0,
            };
        }
        if let AucType::Inf(lambda_z) = kind {
            total += self.clast() / lambda_z;
        }
        Some(total)
    }

    /// PKNCA's `interp.extrap.conc` at `t` (lin up/log down).
    fn conc_at(&self, t: f64, kind: AucType) -> Option<f64> {
        let tl = self.tlast()?;
        if t <= tl {
            if t < self.time[0] {
                return Some(0.0);
            }
            if let Some(i) = self.time.iter().position(|x| *x == t) {
                return Some(self.conc[i]);
            }
            let methods = interval_methods(&self.conc, &self.time, tl, AucType::Last);
            let before = self.time.iter().rposition(|x| *x <= t)?;
            let after = self.time.iter().position(|x| t <= *x)?;
            let (t1, t2, c1, c2) = (
                self.time[before],
                self.time[after],
                self.conc[before],
                self.conc[after],
            );
            Some(match methods[before] {
                Method::Linear => c1 + (t - t1) / (t2 - t1) * (c2 - c1),
                Method::Log => {
                    math::exp(math::ln(c1) + (t - t1) / (t2 - t1) * (math::ln(c2) - math::ln(c1)))
                }
                Method::Zero => 0.0,
            })
        } else {
            match kind {
                AucType::Inf(lambda_z) => Some(self.clast() * math::exp(-lambda_z * (t - tl))),
                AucType::Last => Some(0.0),
                AucType::All => {
                    let tmax_all = *self.time.last().expect("non-empty");
                    if tl == tmax_all {
                        return Some(0.0);
                    }
                    let prev = self.time.iter().rposition(|x| *x <= t)?;
                    if self.conc[prev] == 0.0 {
                        return Some(0.0);
                    }
                    let next = self.time.iter().position(|x| *x >= t)?;
                    let (t1, t2, c1) = (self.time[prev], self.time[next], self.conc[prev]);
                    Some(c1 + (t - t1) / (t2 - t1) * (0.0 - c1))
                }
            }
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
enum Method {
    Linear,
    Log,
    Zero,
}

/// PKNCA's `choose_interval_method` for "lin up/log down" (one entry per interval).
fn interval_methods(conc: &[f64], time: &[f64], tlast: f64, kind: AucType) -> Vec<Method> {
    let n = conc.len();
    let idx_tlast = time
        .iter()
        .position(|t| *t == tlast)
        .expect("tlast is a sample time");
    let mut m: Vec<Method> = (0..n - 1)
        .map(|i| {
            let zero = conc[i] == 0.0 && conc[i + 1] == 0.0;
            let down = conc[i + 1] < conc[i] && conc[i + 1] != 0.0;
            if zero {
                Method::Zero
            } else if down {
                Method::Log
            } else {
                Method::Linear
            }
        })
        .collect();
    for x in m.iter_mut().skip(idx_tlast + 1) {
        *x = Method::Zero;
    }
    if idx_tlast < n - 1 {
        m[idx_tlast] = if matches!(kind, AucType::All) {
            Method::Linear
        } else {
            Method::Zero
        };
    }
    m
}

fn tfirst(time: &[f64], conc: &[f64]) -> Option<f64> {
    time.iter()
        .zip(conc)
        .find(|(_, c)| **c > 0.0)
        .map(|(t, _)| *t)
}

fn tlast(time: &[f64], conc: &[f64]) -> Option<f64> {
    time.iter()
        .zip(conc)
        .filter(|(_, c)| **c > 0.0)
        .map(|(t, _)| *t)
        .next_back()
}

/// Least-squares fit of log concentration on time.
fn fit(pts: &[(f64, f64)]) -> HalfLife {
    let n = pts.len() as f64;
    let mt = pts.iter().map(|p| p.0).sum::<f64>() / n;
    let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
    let sxx: f64 = pts.iter().map(|p| (p.0 - mt) * (p.0 - mt)).sum();
    let sxy: f64 = pts.iter().map(|p| (p.0 - mt) * (p.1 - my)).sum();
    let syy: f64 = pts.iter().map(|p| (p.1 - my) * (p.1 - my)).sum();
    let slope = sxy / sxx;
    let intercept = my - slope * mt;
    let ss_res: f64 = pts
        .iter()
        .map(|p| {
            let r = p.1 - intercept - slope * p.0;
            r * r
        })
        .sum();
    let r_squared = 1.0 - ss_res / syy;
    let adj_r_squared = 1.0 - (1.0 - r_squared) * (n - 1.0) / (n - 2.0);
    HalfLife {
        lambda_z: -slope,
        r_squared,
        adj_r_squared,
        n_points: pts.len(),
        time_first: pts[0].0,
        time_last: pts[pts.len() - 1].0,
        half_life: core::f64::consts::LN_2 / -slope,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Reference values from PKNCA 0.12.1 (pk.calc.* defaults) on the same data.
    #[test]
    fn matches_pknca_on_a_single_dose_profile() {
        let time = [0.0, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 12.0, 24.0];
        let conc = [0.0, 2.1, 3.4, 4.0, 3.1, 2.3, 1.7, 0.95, 0.0];
        let p = Profile::new(&time, &conc);
        assert_eq!(p.cmax(), Some(4.0));
        assert_eq!(p.tmax(), Some(2.0));
        assert_eq!(p.tlast(), Some(12.0));
        let hl = p.half_life().unwrap();
        let close = |a: f64, b: f64| assert!((a - b).abs() < 1e-9, "{a} vs {b}");
        close(hl.lambda_z, 0.147850764966697);
        close(hl.adj_r_squared, 0.999883736218633);
        assert_eq!((hl.n_points, hl.time_first), (4, 4.0));
        let inf = AucType::Inf(hl.lambda_z);
        close(
            p.auc(0.0, f64::INFINITY, AucType::Last).unwrap(),
            27.1472193794398,
        );
        close(
            p.auc(0.0, f64::INFINITY, AucType::All).unwrap(),
            32.8472193794398,
        );
        close(p.auc(0.0, f64::INFINITY, inf).unwrap(), 33.5726173150810);
        close(p.auc(0.0, 10.0, AucType::Last).unwrap(), 24.9419289568902);
        close(p.auc(0.0, 10.0, AucType::All).unwrap(), 24.9419289568902);
        close(p.auc(0.0, 16.0, AucType::Last).unwrap(), 27.1472193794398);
    }

    #[test]
    fn middle_blq_is_dropped_and_all_blq_gives_zero() {
        let p = Profile::new(&[0.0, 1.0, 2.0, 3.0, 4.0], &[0.0, 5.0, 0.0, 2.0, 0.0]);
        assert_eq!(p.time, vec![0.0, 1.0, 3.0, 4.0]);
        let z = Profile::new(&[0.0, 1.0, 2.0], &[0.0, 0.0, 0.0]);
        assert_eq!(z.auc(0.0, f64::INFINITY, AucType::Last), Some(0.0));
        assert_eq!(z.half_life(), None);
    }
}
