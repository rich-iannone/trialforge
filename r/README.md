# trialforge

trialforge generates synthetic CDISC clinical-trial data: trial design datasets, SDTM, ADaM and
define.xml for a complete mock study, from a study type (an *archetype*) and a few interpretable
controls. Every study is verified before it is returned, and the same request gives
byte-identical files in R and Python.

trialforge is pre-alpha. It currently provides:

- thirteen archetypes, from parallel-group studies with a continuous endpoint to oncology
  (RECIST 1.1 and Lugano 2014), first-in-human PK, crossover bioequivalence, dose escalation,
  cardiovascular outcomes, Alzheimer's disease, pediatric growth, COPD exacerbations, vaccines
  and ophthalmology; `archetypes()` lists them and `controls()` what each accepts
- clinical targets as controls (median survival, hazard ratios, response rates), with the latent
  model calibrated to them and a truth manifest of target vs. realized values
- ADaM derived from the SDTM following admiral and its extensions, and PK parameters matching
  PKNCA
- XPT v5, CSV and Dataset-JSON 1.1 files, with a Define-XML 2.1 `define.xml` for each standard
- seeded data-quality defects, with a manifest of each one and `score()` to measure how many of
  them a validator finds

## Example

```r
library(trialforge)

trial <- study("parallel_continuous", n = 300, seed = 1, indication = "hypertension", weeks = 12)
trial
#> <trialforge_study> SIMPAR01 (parallel_continuous)
#>   seed: 1   cutoff: 2025-05-09
#>   SDTM: TA (6), TE (4), TV (16), TI (4), TS (35), DM (300), SE (891),
#>         SV (2041), EX (1154), CM (390), AE (148), DS (1200), MH (1097),
#>         LB (29380), VS (12546), SUPPAE (148)
#>   ADaM: ADSL (300), ADAE (148), ADLB (29380), ADVS (12546)

head(trial$sdtm$DM[, c("USUBJID", "AGE", "SEX", "ARM", "RFSTDTC")], 3)
#>           USUBJID AGE SEX          ARM    RFSTDTC
#> 1 SIMPAR01-201001  63   F Drug A 10 mg 2024-03-01
#> 2 SIMPAR01-201002  56   M Drug A 10 mg 2024-03-13
#> 3 SIMPAR01-201003  65   F      Placebo 2024-03-20

trial$truth$endpoint$realised_difference
#> [1] -7.2283

write_submission(trial, "out/")
```

Seeded defects test a validator against a known answer:

```r
dirty <- study("parallel_continuous", n = 300, seed = 1, defects = "realistic")
score(dirty, verify(dirty))
#> <trialforge_score> recall 100.0% (36 of 36 defects), precision 100.0% (53 of 53 findings)
#> ...
```

## Installation

trialforge is not on CRAN yet. Install it from GitHub; this compiles the Rust engine from source,
so it needs a Rust toolchain (1.88 or later; [rustup](https://rustup.rs) is the easiest way to get
one):

```r
# install.packages("remotes")
remotes::install_github("rich-iannone/trialforge", subdir = "r", build = FALSE)
```

`build = FALSE` matters: the package compiles engine crates from elsewhere in the repository, and
building a source tarball first would leave them behind.

## Development

trialforge is developed alongside a Python package of the same name, in the same repository, and
the [Python documentation](https://rich-iannone.github.io/trialforge/) shares this site. Its user
guide covers the archetypes, the truth manifest, verification, seeded defects and study specs;
the R functions have the same names and behavior.
