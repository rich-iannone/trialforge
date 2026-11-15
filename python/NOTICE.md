# Notices

trialforge's own code is licensed under the MIT License (see `LICENSE`). The packages also contain
material from the third parties listed here, which remains under its owners' terms. None of these
organizations is affiliated with trialforge or endorses it.

If you hold rights in any of this material and have a concern about how trialforge uses it, please
open an issue at https://github.com/rich-iannone/trialforge/issues or write to riannone@me.com, and
the material will be changed or removed promptly.

<!-- Keep the copy in python/NOTICE.md identical (it ships in the wheels and the sdist). -->

## Clinical data standards

### CDISC standards

trialforge generates data to the CDISC Study Data Tabulation Model (SDTM) and its implementation
guide SDTMIG 3.4, the Analysis Data Model (ADaM) and ADaMIG 1.3, Define-XML 2.1 and
Dataset-JSON 1.1. Its metadata bundle (`crates/trialforge-cdisc/assets/sdtmig-3.4/` and
`adamig-1.3/`) records, for the datasets and variables trialforge emits, their names, labels,
types, core designations, structures and keys as the implementation guides define them.

CDISC, SDTM, ADaM, Define-XML and Dataset-JSON are standards of the Clinical Data Interchange
Standards Consortium (CDISC), https://www.cdisc.org. The standards and implementation guides are
© CDISC. trialforge contains none of their text, only the metadata that datasets conforming to
them carry.

### CDISC Controlled Terminology

`crates/trialforge-cdisc/assets/ct/2024-09-27.tsv` is a subset (859 terms in 56 codelists) of
CDISC Controlled Terminology, release 2024-09-27, which CDISC develops and the U.S. National Cancer
Institute's Enterprise Vocabulary Services (NCI EVS) publishes as part of the NCI Thesaurus,
https://evs.nci.nih.gov/ftp1/CDISC/. The NCI Thesaurus, including CDISC Controlled Terminology, is
licensed under the Creative Commons Attribution 4.0 International License (CC BY 4.0),
https://creativecommons.org/licenses/by/4.0/.

trialforge includes a selection of the terms. The file's `source` column records where each term
was taken from; the few marked `sponsor extension` are not CDISC terms but extend extensible
codelists.

## Coding dictionaries

### MedDRA System Organ Class and Preferred Term names

`crates/trialforge-cdisc/assets/dictionaries/meddra_soc.tsv` and `ae_terms.tsv` contain the 27
MedDRA System Organ Class names and 184 common Preferred Term names, each with its primary System
Organ Class. They contain no MedDRA codes and none of the Lowest Level Term, High Level Term or
High Level Group Term levels. This is not MedDRA, and it is not suitable for coding real data: use
a licensed copy of MedDRA for that.

The names are included only so that the synthetic adverse event and medical history records read
realistically. At 211 names, the file is well within the fewer than 1,000 terms that the MedDRA
MSSO allows in published examples. The names remain ICH's and are not covered by trialforge's MIT
License.

MedDRA® the Medical Dictionary for Regulatory Activities terminology is the international medical
terminology developed under the auspices of the International Council for Harmonisation of
Technical Requirements for Pharmaceuticals for Human Use (ICH). MedDRA® is a registered trademark
of ICH. https://www.meddra.org

### RxNorm ingredient names

`crates/trialforge-cdisc/assets/dictionaries/rxnorm_ingredients.tsv` contains 56 RxNorm
ingredient names, in upper case. All are RxNorm's own normalized names (SAB=RXNORM), which NLM
creates and which are in the public domain. They were checked against RxNorm on 2026-11-15 and
are not kept up to date with its releases; use RxNorm itself for current data,
https://www.nlm.nih.gov/research/umls/rxnorm/.

This product uses publicly available data courtesy of the U.S. National Library of Medicine
(NLM), National Institutes of Health, Department of Health and Human Services; NLM is not
responsible for the product and does not endorse or recommend this or any other product.

Nelson SJ, Zeng K, Kilbourne J, Powell T, Moore R. Normalized names for clinical drugs: RxNorm at
6 years. J Am Med Inform Assoc. 2011;18(4):441–448.

### ATC classification

The same file gives each ingredient one of 32 therapeutic subgroup (level 2) names from the
Anatomical Therapeutic Chemical (ATC) classification system, used as the drug class (`CMCLAS`) of
the synthetic concomitant medication records. It contains no ATC codes, DDDs or other part of
the ATC/DDD Index. The names remain under the WHO Collaborating Centre's terms and are not covered
by trialforge's MIT License.

Source: WHO Collaborating Centre for Drug Statistics Methodology, Oslo, Norway,
https://atcddd.fhi.no. © Copyright WHO Collaborating Centre for Drug Statistics Methodology.

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
Disease Assessment Scale – Cognitive Subscale (ADAS-Cog), and gives their maximum scores. The item
names and codes are CDISC Controlled Terminology (see above). The ADAS-Cog itself is copyrighted,
and CDISC developed its terminology with the rights holder's permission. trialforge contains none
of the instrument's text, instructions or scoring manual; administering the scale needs the
rights holder's permission.

Rosen WG, Mohs RC, Davis KL. A new rating scale for Alzheimer's disease. Am J Psychiatry.
1984;141(11):1356–1364.

### Clinical reference data

The other files in `crates/trialforge-processes/assets/` (the lab panel, vital signs, background
adverse event rates, medical history prevalence and medication indications) were written for
trialforge as plausible defaults for an adult trial population, not taken from any published
table. Their test names and codes are CDISC Controlled Terminology and their condition names are
MedDRA Preferred Term names (see above).

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

The Python extension module compiles in third-party Rust crates. Where a crate offers a choice of
licenses that includes MIT, trialforge uses it under MIT, and otherwise under the license it
gives. `THIRD_PARTY_LICENSES.txt`, which ships with this file in the wheels and the sdist
(`python/` in the source repository), lists each crate with its version and license, and
reproduces the license texts, copyright lines and NOTICE files the crates ship.

The Apache Arrow crates (arrow-\*, parquet), used under the Apache License 2.0, carry this notice:

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
