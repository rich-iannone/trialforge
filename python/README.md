# trialforge

Synthetic CDISC clinical-trial data from declarative study specs: trial design datasets, SDTM,
ADaM and define.xml for a complete mock study, with interpretable controls for many study
types and output that is verified to be workable.

trialforge's engine is written in Rust. This is the Python package; an R package built on the
same engine gives byte-identical data from the same spec and seed. Source, issues and the R
package are at [github.com/rich-iannone/trialforge](https://github.com/rich-iannone/trialforge).

> **Status: pre-alpha.** The API may change between releases.

## Installation

```bash
pip install "trialforge[pl]"
```

Datasets come back as DataFrames, so install a backend with an extra: `pl` for Polars (the
default output), `pd` for pandas or `pa` for PyArrow. Writing files needs no backend.
`trialforge[verify]` adds Pointblank for `verify(level=3)`.

Wheels cover Python ≥ 3.10 on Linux (x86_64, aarch64), macOS (Intel, Apple silicon) and Windows
(x86_64). Elsewhere pip builds from the source distribution, which needs a Rust toolchain
(≥ 1.88).

## Usage

```python
import polars as pl
import trialforge as tf

trial = tf.study("parallel_continuous", n=300, seed=1, indication="hypertension", weeks=12)

trial.sdtm["DM"].select("USUBJID", "AGE", "SEX", "RACE", "ARM", "RFSTDTC").head(5)
```

```text
shape: (5, 6)
┌─────────────────┬─────┬─────┬───────┬──────────────┬────────────┐
│ USUBJID         ┆ AGE ┆ SEX ┆ RACE  ┆ ARM          ┆ RFSTDTC    │
│ ---             ┆ --- ┆ --- ┆ ---   ┆ ---          ┆ ---        │
│ str             ┆ i64 ┆ str ┆ str   ┆ str          ┆ str        │
╞═════════════════╪═════╪═════╪═══════╪══════════════╪════════════╡
│ SIMPAR01-201001 ┆ 63  ┆ F   ┆ WHITE ┆ Drug A 10 mg ┆ 2024-03-01 │
│ SIMPAR01-201002 ┆ 56  ┆ M   ┆ WHITE ┆ Drug A 10 mg ┆ 2024-03-13 │
│ SIMPAR01-201003 ┆ 65  ┆ F   ┆ WHITE ┆ Placebo      ┆ 2024-03-20 │
│ SIMPAR01-201004 ┆ 48  ┆ F   ┆ WHITE ┆ Drug A 10 mg ┆ 2024-03-21 │
│ SIMPAR01-201005 ┆ 64  ┆ M   ┆ WHITE ┆ Placebo      ┆ 2024-03-31 │
└─────────────────┴─────┴─────┴───────┴──────────────┴────────────┘
```

`trial.sdtm` and `trial.adam` map dataset names to DataFrames (`output="pandas"` or
`output="pyarrow"` in `tf.study()` returns those types instead). ADaM is derived from the SDTM,
and the truth manifest records what the study was built to show and what it realised:

```python
e = trial.truth["endpoint"]
{k: e[k] for k in ("target_difference", "realised_difference")}
```

```text
{'target_difference': -8.0, 'realised_difference': -7.2283}
```

```python
trial.write_submission("out/")  # sdtm/ and adam/: XPT + define.xml
trial.write("out/", formats=["csv"])  # or "xpt", "dataset-json", "parquet"
trial.verify(level=3)  # + Pointblank's SDTMIG conformance rules
tf.archetypes()  # the study types
tf.controls("parallel_continuous")  # every control, with its default
tf.study_from_spec("study.yaml")  # the same request as a JSON or YAML spec
```

Every study passes trialforge's structural, referential and derivational checks before it is
returned (`verify()` re-runs them after you modify the data), and its ADaM matches admiral:
admiral and its extensions re-derive it from trialforge's SDTM with no differences. PK
parameters match PKNCA.

## Archetypes

| Archetype | What it simulates |
|---|---|
| `parallel_continuous` | Randomised study with a continuous endpoint (HbA1c, SBP or weight), placebo- or active-controlled |
| `dose_ranging_ph2` | Phase 2 dose-ranging study with a sigmoid Emax dose-response |
| `cv_outcomes` | Event-driven cardiovascular outcomes trial: MACE and its components, recurrent events |
| `cns_alzheimers` | Alzheimer's disease study with ADAS-Cog (11) items from an IRT model |
| `oncology_solid_tumor` | Solid-tumour study: RECIST 1.1, confirmed response, PFS and OS calibrated to targets |
| `heme_lymphoma` | Relapsed or refractory DLBCL: Lugano 2014 responses, PFS and OS |
| `phase1_sad_mad` | First-in-human single or multiple ascending dose study with population PK and NCA |
| `crossover_be` | Crossover bioequivalence or food-effect study with per-period PK |
| `dose_escalation` | Phase 1 oncology dose escalation (3+3, BOIN or mTPI-2) with DLTs |
| `pediatric_growth` | Growth hormone study with CDC 2000 growth z-scores and height velocity |
| `respiratory_exacerbation` | COPD study: exacerbations (negative binomial) and FEV1 |
| `vaccine` | One- or two-dose vaccine study: antibody titres and reactogenicity diaries |
| `ophthalmology` | Intravitreal treatment of neovascular AMD: visual acuity and retinal thickness per eye |

## Testing a validator with seeded defects

With `defects=`, trialforge corrupts copies of a verified study's datasets and records every
change: where it is, what it was, and which findings a validator should report. `score()`
matches a validator's findings against that record:

```python
dirty = tf.study("parallel_continuous", n=300, seed=1, defects="realistic")
dirty.defects  # one row per injected defect
score = dirty.score(dirty.verify())
score.recall, score.precision
```

```text
(1.0, 1.0)
```

Profiles are `"light"`, `"realistic"` and `"adversarial"`, or ask for counts and rates per
defect type, e.g. `defects={"out_of_ct": 5, "ae_end_before_start": 0.02}`.

## Dictionaries

Coding uses open dictionaries: MedDRA System Organ Class names with common Preferred Terms (no
MedDRA codes), and RxNorm ingredients with ATC classes.

## License

MIT
