# admiral equivalence: re-derive trialforge's ADaM variables from trialforge's SDTM with
# admiral, and compare value by value.
#
# Usage (from the repo root, via ir so nothing touches the user library):
#   ir run --isolated --r-version 4.6 --exclude-newer 2026-09-25 \
#     --with admiral --with haven --with dplyr scripts/admiral_equivalence.R <submission-dir>
#
# <submission-dir> is the output of `write_submission()`: sdtm/*.xpt and adam/*.xpt. Exits
# non-zero if any compared variable differs.

suppressPackageStartupMessages({
  library(admiral)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: admiral_equivalence.R <submission-dir>")
dir <- args[[1]]

read <- function(kind, name) {
  path <- file.path(dir, kind, paste0(name, ".xpt"))
  if (!file.exists(path)) return(NULL)
  df <- haven::read_xpt(path)
  # XPT has no character missing value: blanks are missing.
  df |> mutate(across(where(is.character), ~ na_if(.x, "")))
}

dm <- read("sdtm", "dm")
ex <- read("sdtm", "ex")
ds <- read("sdtm", "ds")
ae <- read("sdtm", "ae")
lb <- read("sdtm", "lb")
vs <- read("sdtm", "vs")
adsl_tf <- read("adam", "adsl")
adae_tf <- read("adam", "adae")
adlb_tf <- read("adam", "adlb")
advs_tf <- read("adam", "advs")

results <- list()
compare <- function(dataset, ours, theirs, keys, vars, tol = 1e-6) {
  joined <- inner_join(ours, theirs, by = keys, suffix = c(".tf", ".admiral"))
  if (nrow(joined) != nrow(ours) || nrow(joined) != nrow(theirs)) {
    stop(sprintf("%s: %d trialforge rows, %d admiral rows, %d joined", dataset, nrow(ours), nrow(theirs), nrow(joined)))
  }
  for (v in vars) {
    a <- joined[[paste0(v, ".tf")]]
    b <- joined[[paste0(v, ".admiral")]]
    if (inherits(a, "Date") || inherits(b, "Date")) {
      a <- as.numeric(as.Date(a))
      b <- as.numeric(as.Date(b))
    }
    same <- if (is.numeric(a) && is.numeric(b)) {
      (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & abs(a - b) < tol)
    } else {
      (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & as.character(a) == as.character(b))
    }
    n_diff <- sum(!same)
    results[[length(results) + 1]] <<- data.frame(dataset = dataset, variable = v, rows = length(a), differences = n_diff)
    if (n_diff > 0) {
      i <- which(!same)[1:min(3, n_diff)]
      message(sprintf("  %s.%s differs on %d rows, e.g.:", dataset, v, n_diff))
      print(data.frame(joined[i, keys, drop = FALSE], trialforge = a[i], admiral = b[i]))
    }
  }
}

# ---- ADSL -----------------------------------------------------------------------------------
ex_dt <- ex |>
  derive_vars_dt(new_vars_prefix = "EXST", dtc = EXSTDTC) |>
  derive_vars_dt(new_vars_prefix = "EXEN", dtc = EXENDTC)

adsl <- dm |>
  select(STUDYID, USUBJID, DTHDTC, RFPENDTC) |>
  derive_vars_merged(
    dataset_add = ex_dt, by_vars = exprs(STUDYID, USUBJID),
    order = exprs(EXSTDT, EXSEQ), mode = "first", new_vars = exprs(TRTSDT = EXSTDT)
  ) |>
  derive_vars_merged(
    dataset_add = ex_dt, by_vars = exprs(STUDYID, USUBJID),
    order = exprs(EXENDT, EXSEQ), mode = "last", new_vars = exprs(TRTEDT = EXENDT)
  ) |>
  derive_var_trtdurd() |>
  derive_var_merged_exist_flag(
    dataset_add = ex, by_vars = exprs(STUDYID, USUBJID), new_var = SAFFL,
    condition = !is.na(EXSTDTC), false_value = "N", missing_value = "N"
  ) |>
  derive_vars_merged(
    dataset_add = ds, by_vars = exprs(STUDYID, USUBJID), filter_add = DSDECOD == "RANDOMIZED",
    new_vars = exprs(RANDDT = convert_dtc_to_dt(DSSTDTC))
  ) |>
  derive_vars_merged(
    dataset_add = ds, by_vars = exprs(STUDYID, USUBJID), filter_add = DSSCAT == "END OF STUDY",
    new_vars = exprs(EOSDT = convert_dtc_to_dt(DSSTDTC), EOSDECOD = DSDECOD)
  ) |>
  mutate(
    EOSSTT = case_when(is.na(EOSDECOD) ~ "ONGOING", EOSDECOD == "COMPLETED" ~ "COMPLETED", TRUE ~ "DISCONTINUED"),
    DCSREAS = if_else(EOSDECOD != "COMPLETED", EOSDECOD, NA_character_)
  ) |>
  derive_vars_dt(new_vars_prefix = "DTH", dtc = DTHDTC) |>
  derive_vars_dt(new_vars_prefix = "RFPEN", dtc = RFPENDTC) |>
  mutate(LSTALVDT = pmax(TRTEDT, RFPENDT, na.rm = TRUE))

compare(
  "ADSL", adsl_tf, adsl, keys = c("STUDYID", "USUBJID"),
  vars = c("TRTSDT", "TRTEDT", "TRTDURD", "SAFFL", "RANDDT", "EOSDT", "EOSSTT", "DCSREAS", "DTHDT", "LSTALVDT")
)

# ---- ADAE -----------------------------------------------------------------------------------
if (!is.null(ae) && !is.null(adae_tf)) {
  adae <- ae |>
    select(STUDYID, USUBJID, AESEQ, AESTDTC, AEENDTC, AEDECOD, AEBODSYS) |>
    derive_vars_merged(dataset_add = adsl, by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(TRTSDT, TRTEDT)) |>
    derive_vars_dt(new_vars_prefix = "AST", dtc = AESTDTC, highest_imputation = "n") |>
    derive_vars_dt(new_vars_prefix = "AEN", dtc = AEENDTC, highest_imputation = "n") |>
    derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ASTDT, AENDT)) |>
    derive_vars_duration(new_var = ADURN, new_var_unit = ADURU, start_date = ASTDT, end_date = AENDT) |>
    derive_var_trtemfl(
      start_date = ASTDT, end_date = AENDT, trt_start_date = TRTSDT, trt_end_date = TRTEDT, end_window = 30
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID), order = exprs(ASTDT, AESEQ), new_var = AOCCFL, mode = "first"),
      filter = TRTEMFL == "Y"
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID, AEBODSYS), order = exprs(ASTDT, AESEQ), new_var = AOCCSFL, mode = "first"),
      filter = TRTEMFL == "Y"
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID, AEDECOD), order = exprs(ASTDT, AESEQ), new_var = AOCCPFL, mode = "first"),
      filter = TRTEMFL == "Y"
    ) |>
    mutate(ADURU = toupper(ADURU))
  compare(
    "ADAE", adae_tf, adae, keys = c("STUDYID", "USUBJID", "AESEQ"),
    vars = c("ASTDT", "AENDT", "ASTDY", "AENDY", "ADURN", "ADURU", "TRTEMFL", "AOCCFL", "AOCCSFL", "AOCCPFL")
  )
}

# ---- BDS (ADLB, ADVS) -------------------------------------------------------------------------
bds <- function(name, src, prefix, ours) {
  if (is.null(src) || is.null(ours)) return(invisible())
  col <- function(v) rlang::sym(paste0(prefix, v))
  has_ranges <- paste0(prefix, "STNRLO") %in% names(src)
  base <- src |>
    transmute(
      STUDYID, USUBJID, PARAMCD = !!col("TESTCD"), AVAL = !!col("STRESN"),
      SRCSEQ = !!col("SEQ"), DTC = !!col("DTC"),
      ANRLO = if (has_ranges) !!col("STNRLO") else NA_real_,
      ANRHI = if (has_ranges) !!col("STNRHI") else NA_real_
    ) |>
    filter(!is.na(AVAL)) |>
    derive_vars_merged(dataset_add = adsl, by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(TRTSDT, TRTEDT)) |>
    derive_vars_dt(new_vars_prefix = "A", dtc = DTC) |>
    derive_vars_dy(reference_date = TRTSDT, source_vars = exprs(ADT)) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID, PARAMCD), order = exprs(ADT, SRCSEQ), new_var = ABLFL, mode = "last"),
      filter = ADT <= TRTSDT & !is.na(AVAL)
    ) |>
    derive_var_base(by_vars = exprs(STUDYID, USUBJID, PARAMCD), source_var = AVAL, new_var = BASE) |>
    restrict_derivation(derivation = derive_var_chg, filter = ADT > TRTSDT) |>
    restrict_derivation(derivation = derive_var_pchg, filter = ADT > TRTSDT)
  vars <- c("ADT", "ADY", "ABLFL", "BASE", "CHG", "PCHG")
  if (has_ranges) {
    base <- base |>
      derive_var_anrind() |>
      derive_var_base(by_vars = exprs(STUDYID, USUBJID, PARAMCD), source_var = ANRIND, new_var = BNRIND)
    vars <- c(vars, "ANRIND", "BNRIND")
  }
  compare(name, ours, base, keys = c("STUDYID", "USUBJID", "PARAMCD", "SRCSEQ"), vars = vars)
}
bds("ADLB", lb, "LB", adlb_tf)
bds("ADVS", vs, "VS", advs_tf)

summary <- do.call(rbind, results)
cat(sprintf(
  "admiral %s: compared %d variables over %d rows; %d variable(s) differ\n",
  as.character(packageVersion("admiral")), nrow(summary), sum(summary$rows), sum(summary$differences > 0)
))
if (any(summary$differences > 0)) {
  print(summary[summary$differences > 0, ], row.names = FALSE)
  quit(status = 1)
}
