# Notices

trialforge's own code is licensed under the MIT License (see `LICENSE`). The packages also contain
material from the third parties listed here, which remains under its owners' terms. None of these
organizations is affiliated with trialforge or endorses it.

<!--
  DRAFT for the pre-release licensing review (PLAN.md §11 and §16.2). Comments like this one mark
  what to confirm before the first release; remove them once each point is settled. Keep the
  copy in python/NOTICE.md identical (it ships in the wheels and the sdist).
-->

## Clinical data standards

### CDISC standards

trialforge generates data to the CDISC Study Data Tabulation Model (SDTM) and its implementation
guide SDTMIG 3.4, the Analysis Data Model (ADaM) and ADaMIG 1.3, Define-XML 2.1 and
Dataset-JSON 1.1. Its metadata bundle (`crates/trialforge-cdisc/assets/sdtmig-3.4/` and
`adamig-1.3/`) records, for the datasets and variables trialforge emits, their names, labels,
types, core designations, structures and keys as the implementation guides define them.

CDISC, SDTM, ADaM, Define-XML and Dataset-JSON are standards of the Clinical Data Interchange
Standards Consortium (CDISC), https://www.cdisc.org. The implementation guides are © CDISC.

<!--
  Confirm: that redistributing variable names, labels and core designations (not the guides'
  text) is consistent with CDISC's terms of use for its standards, and whether CDISC asks for a
  specific attribution or trademark statement.
-->

### CDISC Controlled Terminology

`crates/trialforge-cdisc/assets/ct/2024-09-27.tsv` is a subset (859 terms in 56 codelists) of
CDISC Controlled Terminology, release 2024-09-27, which CDISC develops and the U.S. National Cancer
Institute's Enterprise Vocabulary Services (NCI EVS) publishes as part of the NCI Thesaurus,
https://evs.nci.nih.gov/ftp1/CDISC/. The file's `source` column records where each term was
taken from.

<!--
  Confirm: the NCI Thesaurus license (believed to be CC BY 4.0, which requires this attribution)
  and any CDISC terms specific to CT.
-->

## Coding dictionaries

### MedDRA System Organ Class and Preferred Term names

`crates/trialforge-cdisc/assets/dictionaries/meddra_soc.tsv` and `ae_terms.tsv` contain the 27
MedDRA System Organ Class names and 184 common Preferred Term names, each with its primary System
Organ Class. They contain no MedDRA codes and none of the Lowest Level Term, High Level Term or
High Level Group Term levels. This is not MedDRA, and it is not suitable for coding real data: use
a licensed copy of MedDRA for that.

MedDRA® the Medical Dictionary for Regulatory Activities terminology is the international medical
terminology developed under the auspices of the International Council for Harmonisation of
Technical Requirements for Pharmaceuticals for Human Use (ICH). MedDRA® is a registered trademark
of ICH. https://www.meddra.org

<!--
  Confirm with the MedDRA MSSO (mssohelp@meddra.org) that an open-source package may redistribute
  these SOC and PT names without codes, and the attribution they require. If not, replace the
  terms with invented MedDRA-shaped ones. This is the review's highest-risk item.
-->

### RxNorm ingredient names

`crates/trialforge-cdisc/assets/dictionaries/rxnorm_ingredients.tsv` contains 56 RxNorm
ingredient names.

This product uses publicly available data from the U.S. National Library of Medicine (NLM),
National Institutes of Health, Department of Health and Human Services; NLM is not responsible
for the product and does not endorse or recommend this or any other product.
https://www.nlm.nih.gov/research/umls/rxnorm/

<!--
  Confirm: that all 56 names come from RxNorm's own (SAB=RXNORM) content, which NLM releases
  without restriction, rather than from a restricted source vocabulary.
-->

### ATC classification

The same file gives each ingredient one of 32 therapeutic subgroup (level 2) names from the
Anatomical Therapeutic Chemical (ATC) classification system of the WHO Collaborating Centre for
Drug Statistics Methodology, Oslo, https://atcddd.fhi.no. It contains no ATC codes or DDDs.

<!--
  Confirm: the WHOCC's terms for reproducing ATC group names (its copyright notice restricts
  copying the ATC/DDD Index for commercial distribution).
-->

WHODrug is not included. Use a licensed copy of WHODrug Global for coding real data.

## Reference data

### CDC growth charts

`crates/trialforge-cdisc/assets/growth/` holds the LMS parameters of the CDC 2000 growth charts
(height, weight and BMI for age), a work of the U.S. Government in the public domain, from the
National Center for Health Statistics, https://www.cdc.gov/growthcharts/. The files are taken
from the R package admiralpeds 0.4.0 (Apache License 2.0; © Cytel Inc., F. Hoffmann-La Roche AG,
Pfizer Inc. and Atorus Research), which packages them.

Kuczmarski RJ, Ogden CL, Guo SS, et al. 2000 CDC growth charts for the United States: methods and
development. Vital Health Stat 11. 2002;(246):1–190.

### ADAS-Cog

`crates/trialforge-processes/assets/adas_items.tsv` names the items of the 11-item Alzheimer's
Disease Assessment Scale – Cognitive Subscale (ADAS-Cog), with their CDISC Controlled Terminology
codes, and their maximum scores. It contains no instrument text, instructions or scoring manual.

Rosen WG, Mohs RC, Davis KL. A new rating scale for Alzheimer's disease. Am J Psychiatry.
1984;141(11):1356–1364.

<!--
  Confirm: whether item names and maximum scores need permission from the ADAS-Cog's rights
  holders.
-->

### Clinical reference data

The other files in `crates/trialforge-processes/assets/` (the lab panel, vital signs, background
adverse event rates, medical history prevalence and medication indications) were written for
trialforge as plausible defaults for an adult trial population. Their test names and codes are
CDISC Controlled Terminology and their condition names are MedDRA Preferred Term names (see
above).

<!--
  Confirm: that no reference range, rate or prevalence was copied from a copyrighted table.
-->

## Software

### Pointblank

The controlled terminology terms marked `pointblank` (507 of the 859) were taken from the CDISC
Controlled Terminology subset bundled with Pointblank, https://github.com/posit-dev/pointblank.
The SDTM and ADaM metadata were written for trialforge and are only cross-checked against
Pointblank's templates. trialforge is an independent project, not a product of Posit Software,
PBC, which holds the copyright to Pointblank. Pointblank's license:

    MIT License

    Copyright (c) 2024-2026 Posit Software, PBC

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

### Rust crates in the compiled extension

The Python extension module includes the following crates, compiled in. Where a crate offers a
choice that includes the MIT License, trialforge uses it under MIT; ryu is used under Apache-2.0.

<!--
  This table reflects Cargo.lock at 0.1.0. MIT and Apache-2.0 both require shipping each crate's
  license text (with its copyright lines) alongside the binary, and the table alone doesn't do
  that. Before the release, generate the full texts (for example with cargo-about) into a file
  the wheels include, and regenerate it whenever Cargo.lock changes.
-->

| Crate | Version | License |
|---|---|---|
| ahash | 0.8.12 | MIT OR Apache-2.0 |
| arrow-array | 59.3.0 | Apache-2.0 AND MIT |
| arrow-buffer, arrow-data, arrow-ipc, arrow-schema, arrow-select | 59.3.0 | Apache-2.0 |
| base64 | 0.23.1 | MIT OR Apache-2.0 |
| bitflags | 2.13.2 | MIT OR Apache-2.0 |
| bytes | 1.12.1 | MIT |
| cfg-if | 1.0.5 | MIT OR Apache-2.0 |
| chrono | 0.4.45 | MIT OR Apache-2.0 |
| core-foundation-sys (macOS) | 0.8.7 | MIT OR Apache-2.0 |
| equivalent | 1.0.2 | Apache-2.0 OR MIT |
| flatbuffers | 25.12.19 | Apache-2.0 |
| getrandom | 0.3.4 | MIT OR Apache-2.0 |
| half | 2.7.1 | MIT OR Apache-2.0 |
| hashbrown | 0.17.1 | MIT OR Apache-2.0 |
| iana-time-zone | 0.1.65 | MIT OR Apache-2.0 |
| indexmap | 2.14.2 | Apache-2.0 OR MIT |
| itoa | 1.0.18 | MIT OR Apache-2.0 |
| libc | 0.2.190 | MIT OR Apache-2.0 |
| libm | 0.2.16 | MIT |
| memchr | 2.8.3 | Unlicense OR MIT |
| num-bigint | 0.5.1 | MIT OR Apache-2.0 |
| num-complex | 0.4.6 | MIT OR Apache-2.0 |
| num-integer | 0.1.47 | MIT OR Apache-2.0 |
| num-traits | 0.2.19 | MIT OR Apache-2.0 |
| once_cell | 1.21.4 | MIT OR Apache-2.0 |
| parquet | 59.3.0 | Apache-2.0 |
| ppv-lite86 | 0.2.21 | MIT OR Apache-2.0 |
| pyo3, pyo3-ffi | 0.29.3 | MIT OR Apache-2.0 |
| rand_chacha | 0.10.0 | MIT OR Apache-2.0 |
| rand_core | 0.10.1 | MIT OR Apache-2.0 |
| ryu | 1.0.23 | Apache-2.0 OR BSL-1.0 |
| serde, serde_core | 1.0.229 | MIT OR Apache-2.0 |
| serde_json | 1.0.154 | MIT OR Apache-2.0 |
| serde_norway | 0.9.42 | MIT OR Apache-2.0 |
| twox-hash | 2.1.5 | MIT |
| unsafe-libyaml-norway | 0.2.15 | MIT |
| windows-link (Windows) | 0.2.1 | MIT OR Apache-2.0 |
| zerocopy | 0.8.62 | BSD-2-Clause OR Apache-2.0 OR MIT |
| zmij | 1.0.23 | MIT |

The Apache Arrow crates (arrow-\*, parquet) carry this notice:

    Apache Arrow
    Copyright 2016-2026 The Apache Software Foundation

    This product includes software developed at
    The Apache Software Foundation (http://www.apache.org/).

    This product includes software from the chronoutil crate (MIT)
     * Copyright (c) 2020-2022 Oliver Margetts
     * https://github.com/olliemath/chronoutil

    This product includes software from the compact-thrift project (Apache 2.0)
     * Copyright Jörn Horstmann
     * https://github.com/jhorstmann/compact-thrift

## Methods

trialforge reimplements published methods; it contains no code from these sources. Its ADaM
derivations follow the admiral family of R packages (admiral, admiralonco, admiralophtha,
admiralpeds, admiralvaccine; Apache License 2.0, https://pharmaverse.org), and its
non-compartmental analysis follows PKNCA (AGPL-3, https://github.com/humanpred/pknca). Both are
used only to check trialforge's output, in scripts that aren't distributed with the packages.

- RECIST 1.1: Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in
  solid tumours: revised RECIST guideline (version 1.1). Eur J Cancer. 2009;45(2):228–247.
- Lugano 2014: Cheson BD, Fisher RI, Barrington SF, et al. Recommendations for initial evaluation,
  staging, and response assessment of Hodgkin and non-Hodgkin lymphoma: the Lugano
  classification. J Clin Oncol. 2014;32(27):3059–3068.
- BOIN: Liu S, Yuan Y. Bayesian optimal interval designs for phase I clinical trials. J R Stat
  Soc Ser C Appl Stat. 2015;64(3):507–523.
- mTPI-2: Guo W, Wang SJ, Yang S, Lynn H, Ji Y. A Bayesian interval dose-finding design addressing
  Ockham's razor: mTPI-2. Contemp Clin Trials. 2017;58:23–33.
