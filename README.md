# trialforge

Synthetic CDISC clinical-trial data from declarative study specs: trial design datasets, SDTM,
ADaM and define.xml for a complete mock study, with interpretable controls for many study
types and output that is verified to be workable.

trialforge has a Rust core with Python and R bindings. The same spec and seed produce identical
data in both languages.

> **Status: pre-alpha.** Eight archetypes generate a complete, verified study: trial design
> (TA, TE, TV, TI, TS), SDTM (DM, SE, SV, EX, CM, AE + SUPPAE, DS, MH, LB, VS; TU, TR, RS for
> oncology; PC, PP for pharmacokinetics; FT for functional tests; CE for clinical events) and
> ADaM (ADSL, ADAE, ADLB, ADVS; ADRS, ADTTE, ADTR; ADPC, ADPP; ADADAS), written as XPT v5, CSV, Dataset-JSON 1.1 or
> Parquet, with a Define-XML 2.1 `define.xml` for each standard.

In Python, datasets come back as Polars DataFrames (`output="pandas"` or `output="pyarrow"` in
`tf.study()` returns those types instead):

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

ADaM datasets are derived from the SDTM. Here is one subject's systolic blood pressure in ADVS:

```python
advs = trial.adam["ADVS"]
advs.filter(pl.col("PARAMCD") == "SYSBP", pl.col("USUBJID") == "SIMPAR01-201001").select(
    "USUBJID", "TRTP", "AVISIT", "AVAL", "BASE", "CHG"
)
```

```text
shape: (7, 6)
┌─────────────────┬──────────────┬───────────┬───────┬───────┬───────┐
│ USUBJID         ┆ TRTP         ┆ AVISIT    ┆ AVAL  ┆ BASE  ┆ CHG   │
│ ---             ┆ ---          ┆ ---       ┆ ---   ┆ ---   ┆ ---   │
│ str             ┆ str          ┆ str       ┆ f64   ┆ f64   ┆ f64   │
╞═════════════════╪══════════════╪═══════════╪═══════╪═══════╪═══════╡
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ SCREENING ┆ 145.0 ┆ 145.0 ┆ null  │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ BASELINE  ┆ 145.0 ┆ 145.0 ┆ null  │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ WEEK 2    ┆ 138.0 ┆ 145.0 ┆ -7.0  │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ WEEK 4    ┆ 122.0 ┆ 145.0 ┆ -23.0 │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ WEEK 8    ┆ 116.0 ┆ 145.0 ┆ -29.0 │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ WEEK 12   ┆ 137.0 ┆ 145.0 ┆ -8.0  │
│ SIMPAR01-201001 ┆ Drug A 10 mg ┆ FOLLOW-UP ┆ 130.0 ┆ 145.0 ┆ -15.0 │
└─────────────────┴──────────────┴───────────┴───────┴───────┴───────┘
```

The truth manifest records what the study was built to show and what it realised:

```python
e = trial.truth["endpoint"]
{k: e[k] for k in ("target_difference", "realised_difference")}
```

```text
{'target_difference': -8.0, 'realised_difference': -7.2283}
```

```python
trial.write_submission("out/")           # sdtm/ and adam/: XPT + define.xml
trial.verify(level=3)                    # + Pointblank's SDTMIG conformance rules
tf.controls("parallel_continuous")       # every control, with its default
```

In R, datasets are data frames with variable labels, so they convert directly to tibbles. The
same spec and seed give the same subjects:

```r
library(trialforge)
library(dplyr)

trial <- study("parallel_continuous", n = 300, seed = 1, indication = "hypertension", weeks = 12)

as_tibble(trial$sdtm$DM) |>
  select(USUBJID, AGE, SEX, RACE, ARM, RFSTDTC) |>
  print(n = 5)
```

```text
# A tibble: 300 × 6
  USUBJID           AGE SEX   RACE  ARM          RFSTDTC
  <chr>           <int> <chr> <chr> <chr>        <chr>
1 SIMPAR01-201001    63 F     WHITE Drug A 10 mg 2024-03-01
2 SIMPAR01-201002    56 M     WHITE Drug A 10 mg 2024-03-13
3 SIMPAR01-201003    65 F     WHITE Placebo      2024-03-20
4 SIMPAR01-201004    48 F     WHITE Drug A 10 mg 2024-03-21
5 SIMPAR01-201005    64 M     WHITE Placebo      2024-03-31
# ℹ 295 more rows
```

Mean change from baseline in systolic blood pressure at week 12, by treatment:

```r
as_tibble(trial$adam$ADVS) |>
  filter(PARAMCD == "SYSBP", AVISIT == "WEEK 12", ANL01FL == "Y") |>
  summarise(n = n(), mean_chg = round(mean(CHG), 1), .by = TRTP)
```

```text
# A tibble: 2 × 3
  TRTP             n mean_chg
  <chr>        <int>    <dbl>
1 Drug A 10 mg   147    -13
2 Placebo        145     -5.7
```

```r
write_submission(trial, "out/")
```

Or describe the study in a spec file and use the command-line tool:

```bash
trialforge generate examples/hypertension.yaml -o out/
```

The same `seed` and controls give byte-identical files from Python and R. Every study passes
trialforge's structural, referential and derivational checks before it is returned
(`verify()` re-runs them after you modify the data), and its ADaM matches admiral: admiral and
admiralonco re-derive it from trialforge's SDTM with no differences. PK parameters (PP) match
R's PKNCA, recomputed from trialforge's PC.

| Archetype | What it simulates |
|---|---|
| `parallel_continuous` | Randomised, double-blind study with a continuous endpoint: `type_2_diabetes` (HbA1c), `hypertension` (SBP) or `obesity` (weight); placebo-controlled, or active-controlled with a non-inferiority margin |
| `dose_ranging_ph2` | Phase 2 dose-ranging study (placebo and several doses) with a sigmoid Emax dose-response; the truth manifest refits the Emax model |
| `cv_outcomes` | Event-driven cardiovascular outcomes trial: recurrent MI and stroke, CV death and competing non-CV death in CE; time to first MACE and its components in ADTTE; the study ends at a target number of events |
| `cns_alzheimers` | Placebo-controlled study in mild to moderate Alzheimer's disease: ADAS-Cog (11) items (FT) from an IRT model of cognitive decline, prorated totals and LOCF analysis |
| `oncology_solid_tumor` | Randomised solid-tumour study: RECIST 1.1 tumour assessments, confirmed response, progression and death, calibrated to target medians, hazard ratios and response rates |
| `phase1_sad_mad` | First-in-human single or multiple ascending dose study in healthy volunteers: dose cohorts with placebo, intensive PK (1- or 2-compartment population model), NCA, dose proportionality |
| `crossover_be` | Crossover bioequivalence (2x2 or 2x2x4 replicate) or food-effect study: periods and washouts, per-period PK and NCA, the BE analysis of AUC and Cmax |
| `dose_escalation` | Phase 1 dose escalation in advanced solid tumours with a 3+3, BOIN or mTPI-2 design: cohorts chosen from DLTs, DLT flags in SUPPAE and ADAE, the DLT-evaluable population, the selected MTD and the design's operating characteristics |
| `heme_lymphoma` | Randomised study in relapsed or refractory DLBCL: Lugano 2014 lesions (perpendicular diameters, SPD, FDG-PET Deauville score), CT- and PET-CT-based responses, ADRS (admiralonco), PFS and OS calibrated to targets |
| `pediatric_growth` | Once-weekly vs daily growth hormone in children: height, weight and BMI on the CDC 2000 charts, growth z-scores and percentiles in ADVS (admiralpeds), annualized height velocity |
| `respiratory_exacerbation` | Placebo-controlled COPD study: moderate and severe exacerbations in CE (negative binomial counts), FEV1 in RE; ADEXAC, time to first exacerbation, ADRE |
| `vaccine` | Placebo-controlled one- or two-dose vaccine study: neutralising antibody titres (IS), 7-day reactogenicity diaries (FACE, CE, VS); ADIS, ADCE and ADFACE following admiralvaccine |
| `ophthalmology` | Active-controlled non-inferiority study of intravitreal Drug A in neovascular AMD: study and fellow eye visual acuity (ETDRS letters) and retinal thickness in OE, the study eye in SC, ADBCVA following admiralophtha |

`oncology_solid_tumor` and `heme_lymphoma` take a crossover-on-progression modifier
(`crossover_prob`, `crossover_os_hr`): control subjects who progress may switch to Drug A, in a
`CONTINUATION TREATMENT` element, and their survival follows the Drug A hazard from then on.

Event-driven studies can also be cut early (multiple cuts): `interim_months` on the oncology
archetypes and `interim_events` on `cv_outcomes` give an interim snapshot of the same study,
with the same subjects and courses, cut at the earlier date. Subjects randomised after the cut
are absent, TS records the interim cutoff, and `trial.truth["interim"]` reports it.

Oncology controls are clinical targets (`os_median_months`, `os_hr`, `pfs_median_months`,
`pfs_hr`, `orr`, …); trialforge calibrates its latent model to them, and `trial.truth` reports
what the generated study realised (Kaplan-Meier medians, stratified Cox hazard ratios with 95%
CIs, confirmed response rates) next to the targets.

```python
trial = tf.study("oncology_solid_tumor", n=400, seed=1, os_hr=0.7, orr=[0.4, 0.2])
trial.truth["overall_survival"]["stratified_cox_hr"]
```

```text
{'ci95': [0.4684, 0.8295], 'comparison': 'Drug A vs Control', 'hr': 0.6233}
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
