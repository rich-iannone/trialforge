//! Portable random draws.
//!
//! Every transcendental function goes through the pure-Rust `libm` crate rather than the
//! platform maths library, so a given stream produces bit-identical draws on every OS. Use
//! [`Draws`] for all randomness in an output path; never call `f64::ln`, `f64::exp` and
//! friends on random values directly.

use rand_chacha::ChaCha8Rng;
use rand_core::Rng;

use crate::rng::Stream;

/// A source of random draws for one stream.
pub struct Draws {
    rng: ChaCha8Rng,
}

impl Draws {
    /// Draws from `stream`.
    pub fn new(stream: Stream) -> Self {
        Self { rng: stream.rng() }
    }

    /// A raw 64-bit draw.
    pub fn next_u64(&mut self) -> u64 {
        self.rng.next_u64()
    }

    /// Uniform on `[0, 1)`, with 53 bits of precision.
    pub fn uniform(&mut self) -> f64 {
        (self.rng.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    /// Uniform on `(0, 1]`; safe to take the log of.
    pub fn uniform_pos(&mut self) -> f64 {
        1.0 - self.uniform()
    }

    /// Uniform on `[lo, hi)`.
    pub fn uniform_range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Uniform integer on `[0, n)`, without modulo bias (Lemire's method). `n` must be > 0.
    pub fn below(&mut self, n: u64) -> u64 {
        assert!(n > 0, "`below` needs n > 0");
        let threshold = n.wrapping_neg() % n;
        loop {
            let m = u128::from(self.rng.next_u64()) * u128::from(n);
            if (m as u64) >= threshold {
                return (m >> 64) as u64;
            }
        }
    }

    /// Uniform integer on `[lo, hi]` (inclusive).
    pub fn int_range(&mut self, lo: i64, hi: i64) -> i64 {
        assert!(lo <= hi, "`int_range` needs lo <= hi");
        lo + self.below((hi - lo) as u64 + 1) as i64
    }

    /// `true` with probability `p`.
    pub fn bernoulli(&mut self, p: f64) -> bool {
        self.uniform() < p
    }

    /// Exponential with the given rate.
    pub fn exponential(&mut self, rate: f64) -> f64 {
        -libm::log(self.uniform_pos()) / rate
    }

    /// Normal with the given mean and standard deviation (Box–Muller; one draw per call).
    pub fn normal(&mut self, mean: f64, sd: f64) -> f64 {
        let r = libm::sqrt(-2.0 * libm::log(self.uniform_pos()));
        let theta = 2.0 * core::f64::consts::PI * self.uniform();
        mean + sd * r * libm::cos(theta)
    }

    /// An index into `weights`, chosen with probability proportional to its weight.
    pub fn weighted(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        assert!(total > 0.0, "`weighted` needs a positive total weight");
        let mut target = self.uniform() * total;
        for (i, &w) in weights.iter().enumerate() {
            if target < w {
                return i;
            }
            target -= w;
        }
        // Floating-point slack: fall back to the last positive weight.
        weights.iter().rposition(|&w| w > 0.0).unwrap_or(0)
    }

    /// Shuffle `items` in place (Fisher–Yates).
    pub fn shuffle<T>(&mut self, items: &mut [T]) {
        for i in (1..items.len()).rev() {
            let j = self.below(i as u64 + 1) as usize;
            items.swap(i, j);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn draws() -> Draws {
        Draws::new(Stream::root(7).label("test"))
    }

    #[test]
    fn uniform_in_range() {
        let mut d = draws();
        for _ in 0..10_000 {
            let u = d.uniform();
            assert!((0.0..1.0).contains(&u));
        }
    }

    #[test]
    fn below_and_int_range_cover_bounds() {
        let mut d = draws();
        let mut seen = [false; 5];
        for _ in 0..1_000 {
            seen[d.below(5) as usize] = true;
            let x = d.int_range(-2, 2);
            assert!((-2..=2).contains(&x));
        }
        assert!(seen.iter().all(|&s| s));
    }

    #[test]
    fn moments_are_plausible() {
        let mut d = draws();
        let n = 50_000;
        let mean_exp = (0..n).map(|_| d.exponential(2.0)).sum::<f64>() / n as f64;
        assert!((mean_exp - 0.5).abs() < 0.02, "{mean_exp}");
        let xs: Vec<f64> = (0..n).map(|_| d.normal(10.0, 3.0)).collect();
        let m = xs.iter().sum::<f64>() / n as f64;
        let v = xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64;
        assert!((m - 10.0).abs() < 0.1, "{m}");
        assert!((v.sqrt() - 3.0).abs() < 0.1, "{v}");
    }

    #[test]
    fn weighted_respects_zero_weights() {
        let mut d = draws();
        for _ in 0..1_000 {
            assert_ne!(d.weighted(&[1.0, 0.0, 2.0]), 1);
        }
    }

    #[test]
    fn shuffle_is_a_permutation() {
        let mut d = draws();
        let mut v: Vec<u32> = (0..20).collect();
        d.shuffle(&mut v);
        let mut sorted = v.clone();
        sorted.sort_unstable();
        assert_eq!(sorted, (0..20).collect::<Vec<_>>());
    }
}
