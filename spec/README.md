# spec/

Language-neutral contracts shared by the Rust core and both bindings. If behaviour is defined
here, the Python and R packages must agree on it exactly.

- `study-spec.schema.json`: the study spec. This is v0, a draft.
- `conformance/`: golden cases that both test suites replay.
  - `rng_streams.json`: RNG stream seeds (`crates/trialforge-core/src/rng.rs`). Seeds are
    decimal strings because R has no unsigned 64-bit integer.
- `designs/`: per-archetype semantics, i.e. controls, defaults, derived truths and expected
  conformance warnings. Empty so far.
- `derivations/`: ADaM derivation rules, each with the admiral function it mirrors. These rules
  are also the source of the define.xml Methods. Empty so far.

Changing a golden value changes generated data. Do it on purpose, update every affected case,
and say so in the commit message.
