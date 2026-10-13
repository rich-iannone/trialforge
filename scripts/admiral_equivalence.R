# admiral equivalence: re-derive trialforge's ADaM variables from trialforge's SDTM with
# admiral, and compare value by value.
#
# Usage (from the repo root, via ir so nothing touches the user library):
#   ir run --isolated --r-version 4.6 --exclude-newer 2026-09-25 \
#     --with admiral --with admiralonco --with haven --with dplyr --with stringr \
#     scripts/admiral_equivalence.R <submission-dir>
#
# Oncology datasets (ADRS, ADTTE, ADTR) are re-derived with admiralonco's templates
# (ad_adrs.R, ad_adtte.R, ad_adtr.R), adapted only where noted.
#
# <submission-dir> is the output of `write_submission()`: sdtm/*.xpt and adam/*.xpt. Exits
# non-zero if any compared variable differs.

suppressPackageStartupMessages({
  library(admiral)
  library(admiralonco)
  library(dplyr)
  library(stringr)
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
tu <- read("sdtm", "tu")
tr <- read("sdtm", "tr")
rs <- read("sdtm", "rs")
adsl_tf <- read("adam", "adsl")
adae_tf <- read("adam", "adae")
adlb_tf <- read("adam", "adlb")
advs_tf <- read("adam", "advs")
adrs_tf <- read("adam", "adrs")
adtte_tf <- read("adam", "adtte")
adtr_tf <- read("adam", "adtr")

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

# ---- Oncology (admiralonco templates) ---------------------------------------------------------
# The templates start from ADSL; use trialforge's, whose RANDDT, DTHDT and LSTALVDT are compared
# with admiral above.
if (!is.null(rs) && !is.null(adrs_tf)) {
  adsl_onc <- adsl_tf |> select(STUDYID, USUBJID, RANDDT, DTHDT, LSTALVDT)
  adsl_vars <- exprs(RANDDT)

  # ad_adrs.R
  adrs <- rs |>
    derive_vars_merged(dataset_add = adsl_onc, new_vars = adsl_vars, by_vars = exprs(STUDYID, USUBJID)) |>
    filter(RSEVAL == "INVESTIGATOR" & RSTESTCD == "OVRLRESP") |>
    mutate(PARAMCD = "OVR") |>
    derive_vars_dt(dtc = RSDTC, new_vars_prefix = "A", highest_imputation = "D", date_imputation = "last") |>
    mutate(AVISIT = VISIT, AVALC = RSSTRESC, AVAL = aval_resp(AVALC))
  worst_resp <- function(arg) {
    case_when(arg == "NE" ~ 1, arg == "CR" ~ 2, arg == "PR" ~ 3, arg == "SD" ~ 4,
              arg == "NON-CR/NON-PD" ~ 5, arg == "PD" ~ 6, TRUE ~ 0)
  }
  adrs <- adrs |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID, ADT), order = exprs(worst_resp(AVALC), RSSEQ), new_var = ANL01FL, mode = "last"),
      filter = !is.na(AVAL) & ADT >= RANDDT
    ) |>
    derive_var_relative_flag(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(ADT, RSSEQ), new_var = ANL02FL,
      condition = AVALC == "PD", mode = "first", selection = "before", inclusive = TRUE
    )
  ovr <- filter(adrs, PARAMCD == "OVR" & ANL01FL == "Y" & ANL02FL == "Y")
  no_data_n <- event(dataset_name = "adsl", condition = TRUE, set_values_to = exprs(AVALC = "N"), keep_source_vars = adsl_vars)
  no_data_missing <- event(dataset_name = "adsl", condition = TRUE, set_values_to = exprs(AVALC = "MISSING"), keep_source_vars = adsl_vars)
  sources <- list(ovr = ovr, adsl = adsl_onc)
  adrs <- adrs |>
    derive_extreme_records(
      dataset_ref = adsl_onc, dataset_add = adrs, by_vars = exprs(STUDYID, USUBJID),
      filter_add = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y", order = exprs(ADT, RSSEQ),
      mode = "first", exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "PD", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", events = list(rsp_y, no_data_n), source_datasets = sources,
      set_values_to = exprs(PARAMCD = "RSP", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(desc(AVALC), ADT, event_nr), tmp_event_nr_var = event_nr,
      mode = "first", events = list(rsp_y, cb_y, no_data_n), source_datasets = sources,
      set_values_to = exprs(PARAMCD = "CB", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y"), check_type = "none"
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", source_datasets = sources,
      events = list(bor_cr, bor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
      set_values_to = exprs(PARAMCD = "BOR", AVAL = aval_resp(AVALC), ANL01FL = "Y")
    )
  # BCP and CBCP read BOR and CBOR from `adrs`, so they need the updated dataset: a chained
  # call would see the previous one (the template chains CBCP, so its CBCP is always "N").
  adrs <- adrs |>
    derive_extreme_records(
      dataset_ref = adsl_onc, dataset_add = adrs, by_vars = exprs(STUDYID, USUBJID),
      filter_add = PARAMCD == "BOR" & AVALC %in% c("CR", "PR"), exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "BCP", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(desc(AVALC), ADT, event_nr), tmp_event_nr_var = event_nr,
      mode = "first", source_datasets = sources, events = list(crsp_y_cr, crsp_y_pr, no_data_n),
      set_values_to = exprs(PARAMCD = "CRSP", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(desc(AVALC), ADT, event_nr), tmp_event_nr_var = event_nr,
      mode = "first", events = list(crsp_y_cr, crsp_y_pr, cb_y, no_data_n), source_datasets = sources,
      set_values_to = exprs(PARAMCD = "CCB", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_extreme_event(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(event_nr, ADT), tmp_event_nr_var = event_nr,
      mode = "first", source_datasets = sources,
      events = list(cbor_cr, cbor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
      set_values_to = exprs(PARAMCD = "CBOR", AVAL = aval_resp(AVALC), ANL01FL = "Y")
    )
  adrs <- adrs |>
    derive_extreme_records(
      dataset_ref = adsl_onc, dataset_add = adrs, by_vars = exprs(STUDYID, USUBJID),
      filter_add = PARAMCD == "CBOR" & AVALC %in% c("CR", "PR"), exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "CBCP", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    )
  adsldth <- adsl_onc |> select(STUDYID, USUBJID, DTHDT, RANDDT)
  adrs <- adrs |>
    derive_extreme_records(
      dataset_ref = adsldth, dataset_add = adsldth, by_vars = exprs(STUDYID, USUBJID),
      filter_add = !is.na(DTHDT), exist_flag = AVALC, false_value = "N",
      set_values_to = exprs(PARAMCD = "DEATH", ANL01FL = "Y", AVAL = yn_to_numeric(AVALC), ADT = DTHDT)
    ) |>
    select(-DTHDT) |>
    derive_extreme_records(
      dataset_ref = adsl_onc, dataset_add = adrs, by_vars = exprs(STUDYID, USUBJID),
      filter_add = PARAMCD == "OVR" & ANL01FL == "Y", order = exprs(ADT, RSSEQ), mode = "last",
      set_values_to = exprs(PARAMCD = "LSTA", ANL01FL = "Y")
    ) |>
    derive_param_exist_flag(
      dataset_ref = adsl_onc |> select(STUDYID, USUBJID, RANDDT), dataset_add = tu,
      condition = TUEVAL == "INVESTIGATOR" & TUSTRESC == "TARGET" & VISIT == "SCREENING",
      false_value = "N", missing_value = "N",
      set_values_to = exprs(PARAMCD = "MDIS", AVAL = yn_to_numeric(AVALC), ANL01FL = "Y")
    ) |>
    derive_var_obs_number(by_vars = exprs(STUDYID, USUBJID), order = exprs(PARAMCD, ADT, VISITNUM, RSSEQ), check_type = "error")
  compare(
    "ADRS", adrs_tf, adrs, keys = c("STUDYID", "USUBJID", "ASEQ"),
    vars = c("PARAMCD", "AVALC", "AVAL", "ADT", "ANL01FL", "ANL02FL")
  )

  # ad_adtte.R (with admiral's ADRS from above)
  adsl_tte <- adsl_onc |>
    derive_vars_merged(
      dataset_add = adrs, filter_add = PARAMCD == "RSP" & AVALC == "Y" & ANL01FL == "Y",
      by_vars = exprs(STUDYID, USUBJID), new_vars = exprs(TEMP_RESPDT = ADT)
    )
  srcs <- list(adsl = adsl_tte, adrs = adrs)
  adtte <- derive_param_tte(
    dataset_adsl = adsl_tte, start_date = RANDDT, event_conditions = list(death_event),
    censor_conditions = list(lastalive_censor, rand_censor), source_datasets = srcs,
    set_values_to = exprs(PARAMCD = "OS")
  ) |>
    derive_param_tte(
      dataset_adsl = adsl_tte, start_date = RANDDT, event_conditions = list(pd_event, death_event),
      censor_conditions = list(lasta_censor, rand_censor), source_datasets = srcs,
      set_values_to = exprs(PARAMCD = "PFS")
    ) |>
    derive_param_tte(
      dataset_adsl = filter(adsl_tte, !is.na(TEMP_RESPDT)), start_date = TEMP_RESPDT,
      event_conditions = list(pd_event, death_event), censor_conditions = list(lasta_censor),
      source_datasets = srcs, set_values_to = exprs(PARAMCD = "RSD")
    ) |>
    derive_vars_duration(new_var = AVAL, start_date = STARTDT, end_date = ADT) |>
    derive_var_obs_number(by_vars = exprs(STUDYID, USUBJID), order = exprs(PARAMCD), check_type = "error")
  compare(
    "ADTTE", adtte_tf, adtte, keys = c("STUDYID", "USUBJID", "PARAMCD"),
    vars = c("ASEQ", "STARTDT", "ADT", "AVAL", "CNSR", "EVNTDESC", "CNSDTDSC", "SRCDOM", "SRCVAR", "SRCSEQ")
  )
}

if (!is.null(tr) && !is.null(adtr_tf)) {
  # ad_adtr.R. Adaptation: nodal target lesions are measured by short axis (TRTESTCD SAXIS,
  # RECIST 1.1), which takes the template's LPERP role (NLDIAMnn); lesion numbers are the
  # digits after the "T" of TRLNKID.
  tu_inv <- filter(tu, TUEVAL == "INVESTIGATOR")
  trl <- tr |>
    filter(TREVAL == "INVESTIGATOR" & TRGRPID == "TARGET" & TRTESTCD %in% c("LDIAM", "SAXIS")) |>
    derive_vars_merged(dataset_add = adsl_onc, new_vars = exprs(RANDDT), by_vars = exprs(STUDYID, USUBJID)) |>
    derive_vars_merged(dataset_add = tu_inv, new_vars = exprs(TULOC), by_vars = exprs(STUDYID, USUBJID, TRLNKID = TULNKID)) |>
    mutate(
      TULOCGR1 = if_else(TULOC == "LYMPH NODE", "NODAL", "NON-NODAL"),
      LSEXP = TRLNKID, LSASS = if_else(!is.na(TRSTRESN), TRLNKID, NA_character_)
    ) |>
    derive_vars_dt(dtc = TRDTC, new_vars_prefix = "A", highest_imputation = "D", date_imputation = "first") |>
    derive_vars_dy(reference_date = RANDDT, source_vars = exprs(ADT)) |>
    mutate(
      AVISIT = if_else(VISIT == "SCREENING", "BASELINE", VISIT),
      AVISITN = if_else(AVISIT == "BASELINE", 0, VISITNUM),
      tmp_lesion_nr = str_sub(TRLNKID, 2)
    )
  adtr <- bind_rows(
    trl |> filter(TRTESTCD == "LDIAM") |> mutate(PARAMCD = paste0("LDIAM", tmp_lesion_nr)),
    trl |> filter(TRTESTCD == "SAXIS") |> mutate(PARAMCD = paste0("NLDIAM", tmp_lesion_nr))
  ) |>
    mutate(AVAL = TRSTRESN, ANL01FL = if_else(!is.na(AVAL), "Y", NA_character_)) |>
    select(-tmp_lesion_nr)
  in_sum <- quote((str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") | (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL"))
  adtr_sum <- derive_summary_records(
    dataset_add = adtr, by_vars = exprs(STUDYID, USUBJID, RANDDT, AVISIT, AVISITN),
    filter_add = !!in_sum,
    set_values_to = exprs(AVAL = sum(AVAL, na.rm = TRUE), ADY = min(ADY, na.rm = TRUE), ADT = min(ADT, na.rm = TRUE), PARAMCD = "SDIAM")
  ) |>
    derive_vars_merged_summary(
      dataset_add = adtr, by_vars = exprs(STUDYID, USUBJID),
      filter_add = AVISIT == "BASELINE" & !!in_sum,
      new_vars = exprs(LSEXP = paste(sort(TRLNKID), collapse = ", "))
    ) |>
    derive_vars_merged_summary(
      dataset_add = adtr, by_vars = exprs(STUDYID, USUBJID, AVISIT),
      filter_add = !!in_sum & ANL01FL == "Y",
      new_vars = exprs(LSASS = paste(sort(TRLNKID), collapse = ", "))
    ) |>
    mutate(ANL01FL = if_else(LSEXP == LSASS, "Y", NA_character_)) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID), order = exprs(ADY), new_var = ABLFL, mode = "last"),
      filter = ADY <= 1
    ) |>
    derive_var_base(by_vars = exprs(STUDYID, USUBJID))
  adtr_sum <- adtr_sum |>
    derive_vars_joined(
      dataset_add = adtr_sum, by_vars = exprs(STUDYID, USUBJID), order = exprs(AVAL), new_vars = exprs(NADIR = AVAL),
      join_vars = exprs(ADY), join_type = "all", filter_add = ANL01FL == "Y", filter_join = ADY.join < ADY,
      mode = "first", check_type = "none"
    ) |>
    derive_var_chg() |>
    derive_var_pchg() |>
    mutate(CHGNAD = AVAL - NADIR, PCHGNAD = if_else(NADIR == 0, NA_real_, 100 * CHGNAD / NADIR)) |>
    derive_var_merged_exist_flag(
      dataset_add = derive_vars_dt(rs, dtc = RSDTC, new_vars_prefix = "A", highest_imputation = "D", flag_imputation = "none"),
      filter_add = RSTESTCD == "OVRLRESP" & RSEVAL == "INVESTIGATOR", by_vars = exprs(STUDYID, USUBJID, ADT),
      new_var = PDFL, condition = RSSTRESC == "PD"
    ) |>
    derive_var_relative_flag(
      by_vars = exprs(STUDYID, USUBJID), order = exprs(ADT), new_var = POSTRNDFL, condition = ADT > RANDDT,
      mode = "first", selection = "after", inclusive = TRUE, flag_no_ref_groups = FALSE
    ) |>
    restrict_derivation(
      derivation = derive_var_extreme_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID), new_var = ANL02FL, order = exprs(PCHG), mode = "first", check_type = "none"),
      filter = ANL01FL == "Y" & POSTRNDFL == "Y"
    ) |>
    select(-POSTRNDFL) |>
    restrict_derivation(
      derivation = derive_var_relative_flag,
      args = params(by_vars = exprs(STUDYID, USUBJID), new_var = ANL03FL, condition = PDFL == "Y", order = exprs(ADY),
                    mode = "first", selection = "before", inclusive = FALSE),
      filter = ANL01FL == "Y" | PDFL == "Y"
    ) |>
    mutate(ANL04FL = if_else(ANL01FL == "Y" | PDFL == "Y", "Y", NA_character_))
  adtr <- bind_rows(adtr, adtr_sum) |>
    derive_var_obs_number(by_vars = exprs(STUDYID, USUBJID), order = exprs(PARAMCD, AVISITN, TRSEQ), check_type = "error")
  compare(
    "ADTR", adtr_tf, adtr, keys = c("STUDYID", "USUBJID", "ASEQ"),
    vars = c("PARAMCD", "AVISITN", "AVAL", "ADT", "ADY", "LSEXP", "LSASS", "BASE", "CHG", "PCHG", "NADIR", "CHGNAD",
             "PCHGNAD", "ABLFL", "PDFL", "ANL01FL", "ANL02FL", "ANL03FL", "ANL04FL")
  )
}

summary <- do.call(rbind, results)
cat(sprintf(
  "admiral %s: compared %d variables over %d rows; %d variable(s) differ\n",
  as.character(packageVersion("admiral")), nrow(summary), sum(summary$rows), sum(summary$differences > 0)
))
if (any(summary$differences > 0)) {
  print(summary[summary$differences > 0, ], row.names = FALSE)
  quit(status = 1)
}
