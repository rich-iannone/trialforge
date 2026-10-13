# trialforge

Synthetic CDISC clinical-trial data from declarative study specs: trial design datasets, SDTM,
ADaM and define.xml for a complete mock study, with interpretable controls for many study
types and output that is verified to be workable.

trialforge has a Rust core with Python and R bindings. The same spec and seed produce identical
data in both languages.

> **Status: pre-alpha.** Two archetypes generate a complete, verified study: trial design
> (TA, TE, TV, TI, TS), SDTM (DM, SE, SV, EX, CM, AE + SUPPAE, DS, MH, LB, VS; TU, TR, RS for
> oncology) and ADaM (ADSL, ADAE, ADLB, ADVS; ADRS, ADTTE, ADTR for oncology), written as XPT
> v5, CSV, Dataset-JSON 1.1 or Parquet, with a Define-XML 2.1 `define.xml` for each standard.

```python
import trialforge as tf

trial = tf.study("parallel_continuous", n=300, seed=1, indication="hypertension", weeks=12)
trial.adam["ADVS"]                       # polars by default
trial.truth["endpoint"]                  # target vs realised treatment effect
trial.write_submission("out/")           # sdtm/ and adam/: XPT + define.xml
trial.verify(level=3)                    # + Pointblank's SDTMIG conformance rules
tf.controls("parallel_continuous")       # every control, with its default
```

```r
library(trialforge)

trial <- study("parallel_continuous", n = 300, seed = 1, indication = "hypertension", weeks = 12)
trial$adam$ADVS
trial$truth$endpoint
write_submission(trial, "out/")
```

Or describe the study in a spec file and use the command-line tool:

```bash
trialforge generate examples/hypertension.yaml -o out/
```

The same `seed` and controls give byte-identical files from Python and R. Every study passes
trialforge's structural, referential and derivational checks before it is returned
(`verify()` re-runs them after you modify the data), and its ADaM matches admiral: admiral and
admiralonco re-derive it from trialforge's SDTM with no differences.

| Archetype | What it simulates |
|---|---|
| `parallel_continuous` | Randomised, double-blind, placebo-controlled study with a continuous endpoint: `type_2_diabetes` (HbA1c), `hypertension` (SBP) or `obesity` (weight) |
| `oncology_solid_tumor` | Randomised solid-tumour study: RECIST 1.1 tumour assessments, confirmed response, progression and death, calibrated to target medians, hazard ratios and response rates |

Oncology controls are clinical targets (`os_median_months`, `os_hr`, `pfs_median_months`,
`pfs_hr`, `orr`, …); trialforge calibrates its latent model to them, and `trial.truth` reports
what the generated study realised (Kaplan-Meier medians, stratified Cox hazard ratios with 95%
CIs, confirmed response rates) next to the targets.

```python
trial = tf.study("oncology_solid_tumor", n=400, seed=1, os_hr=0.7, orr=[0.4, 0.2])
trial.truth["overall_survival"]["stratified_cox_hr"]   # {'hr': ..., 'ci95': [..., ...]}
```

Coding uses open dictionaries: MedDRA System Organ Class names with common Preferred Terms (no
MedDRA codes), and RxNorm ingredients with ATC classes. Bring-your-own MedDRA and WHODrug are
planned.

## Layout

| Path | What |
|---|---|
| `crates/` | Rust workspace. `trialforge` is the facade crate; `trialforge-*` are the engine crates. |
| `bindings/python/` | Python package (pyo3 + maturin). |
| `bindings/r/` | R package (extendr). |
| `spec/` | Language-neutral contracts: the study-spec JSON Schema and golden conformance cases that both bindings replay. |

## Development

Requirements: Rust ≥ 1.85, [uv](https://docs.astral.sh/uv/), and R ≥ 4.2 with `devtools` and
`rextendr`.

```bash
make py-install   # build the extension into .venv (re-run after Rust changes)
make r-install    # build and install the R package from the source tree
make check        # Rust, Python and R gates
make help         # all targets
```

## License

MIT
