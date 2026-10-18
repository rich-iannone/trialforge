//! Hierarchical, deterministic RNG streams.
//!
//! Every random draw in trialforge comes from a [`Stream`] addressed by a path from the study
//! seed, for example `study seed → "subject" → 17 → "ae"`. A stream's seed depends only on its
//! path, so:
//!
//! - a subject's own draws do not depend on how many other subjects there are (enrollment
//!   order and block randomisation, being sequential, can still shift when `n` changes);
//! - turning a domain on or off (e.g. EG) does not change any other domain;
//! - user hooks get their own derived seeds and stay deterministic.
//!
//! The mixing functions and constants below are part of the output contract: changing them
//! changes every generated dataset. The golden-value tests at the bottom of this file pin them.

use rand_chacha::ChaCha8Rng;
use rand_core::SeedableRng;

/// Tag mixed in before an integer path component, so `index(n)` and `label(..)` never collide.
const INDEX_TAG: u64 = 0x5452_4941_4c46_4f52; // "TRIALFOR"
/// Tag mixed in before a label path component.
const LABEL_TAG: u64 = 0x4745_5354_5245_414d; // "GESTREAM"

/// A node in the stream tree. Cheap to copy; derive children with [`Stream::index`] and
/// [`Stream::label`], then draw with [`Stream::rng`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Stream {
    key: u64,
}

impl Stream {
    /// The root stream for a study seed.
    pub fn root(seed: u64) -> Self {
        Self {
            key: splitmix64(seed),
        }
    }

    /// The child stream for an integer path component (a subject index, a replicate, …).
    #[must_use]
    pub fn index(self, i: u64) -> Self {
        self.child(splitmix64(INDEX_TAG ^ i))
    }

    /// The child stream for a named path component (a process, a domain, a hook, …).
    #[must_use]
    pub fn label(self, name: &str) -> Self {
        self.child(splitmix64(LABEL_TAG ^ fnv1a64(name.as_bytes())))
    }

    /// The 64-bit seed of this stream, for handing to user hooks.
    pub fn seed(self) -> u64 {
        self.key
    }

    /// A ChaCha8 generator seeded from this stream.
    pub fn rng(self) -> ChaCha8Rng {
        let mut state = self.key;
        let mut seed = [0u8; 32];
        for chunk in seed.as_chunks_mut::<8>().0 {
            state = state.wrapping_add(GOLDEN_GAMMA);
            chunk.copy_from_slice(&mix64(state).to_le_bytes());
        }
        ChaCha8Rng::from_seed(seed)
    }

    fn child(self, component: u64) -> Self {
        Self {
            key: splitmix64(self.key ^ component),
        }
    }
}

const GOLDEN_GAMMA: u64 = 0x9e37_79b9_7f4a_7c15;

/// SplitMix64: one step from `x`.
fn splitmix64(x: u64) -> u64 {
    mix64(x.wrapping_add(GOLDEN_GAMMA))
}

/// The SplitMix64 output finaliser (Stafford variant 13).
fn mix64(mut z: u64) -> u64 {
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

/// 64-bit FNV-1a. Used instead of `std::hash` because that is not stable across Rust releases.
fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in bytes {
        h ^= u64::from(b);
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand_core::Rng;

    #[test]
    fn same_path_same_stream() {
        let a = Stream::root(1).label("subject").index(17).label("ae");
        let b = Stream::root(1).label("subject").index(17).label("ae");
        assert_eq!(a, b);
        assert_eq!(a.rng().next_u64(), b.rng().next_u64());
    }

    #[test]
    fn siblings_and_seeds_differ() {
        let subj = Stream::root(1).label("subject");
        assert_ne!(subj.index(0), subj.index(1));
        assert_ne!(subj.label("ae"), subj.label("lb"));
        assert_ne!(Stream::root(1), Stream::root(2));
        // An index and a label never alias, even for the "same" value.
        assert_ne!(subj.index(fnv1a64(b"ae")), subj.label("ae"));
    }

    #[test]
    fn path_order_matters() {
        let r = Stream::root(1);
        assert_ne!(r.label("a").label("b"), r.label("b").label("a"));
    }

    /// Golden values: if this fails, every generated dataset has changed. Only update these
    /// deliberately, together with the cross-language conformance hashes in `spec/conformance/`.
    #[test]
    fn golden_values() {
        let s = Stream::root(20261011).label("subject").index(0).label("ae");
        assert_eq!(s.seed(), GOLDEN_SEED);
        assert_eq!(s.rng().next_u64(), GOLDEN_FIRST_DRAW);
    }

    const GOLDEN_SEED: u64 = 12342615029409701309;
    const GOLDEN_FIRST_DRAW: u64 = 18063926280474641268;
}
