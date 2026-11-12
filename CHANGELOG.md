# Changelog

The Python package, the R package and the Rust crates share version numbers. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## 0.1.0 (unreleased)

First release.

- `study()` generates a complete mock study from one of thirteen archetypes:
  `parallel_continuous`, `dose_ranging_ph2`, `cv_outcomes`, `cns_alzheimers`,
  `oncology_solid_tumor`, `heme_lymphoma`, `phase1_sad_mad`, `crossover_be`,
  `dose_escalation`, `pediatric_growth`, `respiratory_exacerbation`, `vaccine` and
  `ophthalmology`. Each takes typed controls (`controls()` lists them), and `study_from_spec()`
  reads the same request from a JSON or YAML study spec.
- Trial design datasets (TA, TE, TV, TI, TS), SDTMIG 3.4 domains and ADaMIG 1.3 datasets,
  written as XPT v5, CSV, Dataset-JSON 1.1 or Parquet (Python only), with a Define-XML 2.1
  `define.xml` per standard (`write()`, `write_submission()`).
- Every study passes trialforge's structural, referential and derivational checks before it is
  returned; `verify()` re-runs them, and `verify(level=3)` (Python) adds Pointblank's SDTMIG
  conformance rules. ADaM matches admiral and its extensions, and PP matches PKNCA.
- Clinical targets (medians, hazard ratios, response rates, effect sizes) are calibrated, and
  the truth manifest (`truth`) reports target and realized values.
- Seeded data-quality defects (`defects=`), with a manifest of every change and `score()` for
  a validator's recall and precision.
- The same seed and controls give byte-identical files from Python and R, on every platform.
- A `trialforge` command-line tool (`generate`, `summary`, `archetypes`, `controls`).
- Documentation at https://rich-iannone.github.io/trialforge/: a user guide and the Python API
  reference, with the R reference under `r/`.
