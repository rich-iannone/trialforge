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
