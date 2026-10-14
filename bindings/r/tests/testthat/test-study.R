par <- study("parallel_continuous", n = 80, seed = 3)
onc <- study("oncology_solid_tumor", n = 60, seed = 3)

test_that("archetypes(), controls() and formats() describe the build", {
  expect_setequal(names(archetypes()), c("parallel_continuous", "oncology_solid_tumor", "phase1_sad_mad", "crossover_be", "dose_ranging_ph2", "cns_alzheimers", "cv_outcomes"))
  ctl <- controls("parallel_continuous")
  expect_true(all(c("indication", "weeks") %in% ctl$name))
  expect_identical(ctl$default[[which(ctl$name == "weeks")]], 24L)
  expect_true(all(c("xpt", "csv") %in% formats()))
})

test_that("a parallel study has every Phase 1 dataset", {
  expect_s3_class(par, "trialforge_study")
  expect_identical(par$studyid, "SIMPAR01")
  expect_identical(par$seed, 3)
  expect_named(par$sdtm, c("TA", "TE", "TV", "TI", "TS", "DM", "SE", "SV", "EX", "CM", "AE", "DS", "MH", "LB", "VS", "SUPPAE"))
  expect_named(par$adam, c("ADSL", "ADAE", "ADLB", "ADVS"))
  expect_identical(nrow(par$sdtm$DM), 80L)
  expect_s3_class(par$adam$ADSL$TRTSDT, "Date")
  expect_type(par$sdtm$DM$AGE, "integer")
  expect_identical(attr(par$sdtm$DM$USUBJID, "label"), "Unique Subject Identifier")
  expect_identical(attr(par$adam$ADSL, "label"), "Subject-Level Analysis Dataset")
  expect_identical(par$truth$endpoint$paramcd, "HBA1C")
  expect_true(verify(par)$passed)
  expect_output(print(par), "SIMPAR01")
})

test_that("an oncology study keeps its strata and cutoff", {
  expect_identical(onc$cutoff, as.Date("2025-08-30"))
  expect_true("STRAT1R" %in% names(onc$adam$ADSL))
  expect_true("overall_survival" %in% names(onc$truth))
  expect_true(all(c("TU", "TR", "RS") %in% names(onc$sdtm)))
  expect_named(onc$adam, c("ADSL", "ADRS", "ADTTE", "ADTR"))
  expect_true(all(c("OS", "PFS") %in% onc$adam$ADTTE$PARAMCD))
  expect_identical(onc$truth$targets$os_hr, 0.75)
  expect_true("stratified_cox_hr" %in% names(onc$truth$overall_survival))
})

test_that("oncology controls are applied and validated", {
  s <- study("oncology_solid_tumor", n = 30, seed = 2, orr = c(0.5, 0.1), os_median_months = 20)
  expect_identical(s$controls$orr, c(0.5, 0.1))
  expect_error(study("oncology_solid_tumor", n = 5, seed = 1, pfs_median_months = 30), "pfs_median_months")
})

test_that("controls are applied and validated", {
  s <- study("parallel_continuous", n = 20, seed = 1, indication = "hypertension", weeks = 12)
  expect_identical(s$controls$weeks, 12L)
  expect_identical(s$truth$endpoint$paramcd, "SYSBP")
  expect_error(study("parallel_continuous", n = 5, seed = 1, wekks = 3), "unknown control")
  expect_error(study("parallel_continuous", n = 5, seed = 1, weeks = "twelve"), "invalid controls")
  expect_error(study("parallel_continuous", 5, 1, NULL, 12), "must be named")
  expect_error(study("parallel_continuous", n = 5, seed = 1, 12), "`studyid` must be a single string")
})

test_that("a missing seed is drawn and recorded", {
  a <- study("parallel_continuous", n = 10)
  expect_true(a$seed >= 0 && a$seed <= 2^53)
  b <- study("parallel_continuous", n = 10, seed = a$seed)
  expect_identical(a$sdtm, b$sdtm)
  expect_no_error(study("parallel_continuous", n = 1, seed = 2^53))
  expect_error(study("parallel_continuous", n = 1, seed = 2^53 + 2), "2\\^53")
})

test_that("bad input is rejected", {
  expect_error(study("nope", n = 10, seed = 1), "unknown archetype")
  expect_error(study("parallel_continuous", n = 0, seed = 1), "at least 1")
  expect_error(write_study(par, tempfile(), formats = "sas7bdat"), "unknown format")
  expect_error(write_study(list(), tempfile()), "trialforge_study")
  expect_error(verify(list()), "trialforge_study")
})

for (trial in list(par, onc)) {
  test_that(paste("XPT files read back with haven match the in-memory data:", trial$archetype), {
    skip_if_not_installed("haven")
    paths <- write_study(trial, tempfile())
    xpts <- paths[grepl("\\.xpt$", paths)]
    frames <- c(trial$sdtm, trial$adam)
    expect_setequal(toupper(tools::file_path_sans_ext(basename(xpts))), names(frames))
    for (path in xpts) {
      name <- toupper(tools::file_path_sans_ext(basename(path)))
      back <- haven::read_xpt(path)
      mem <- frames[[name]]
      expect_identical(names(back), names(mem))
      expect_identical(nrow(back), nrow(mem))
      for (col in names(mem)) {
        expected <- mem[[col]]
        got <- back[[col]]
        # XPT has no character missing value: missing reads back as "".
        if (is.character(expected)) expected[is.na(expected)] <- ""
        if (is.integer(expected)) got <- as.integer(got)
        expect_identical(as.vector(got), as.vector(expected), info = paste0(name, ".", col))
        expect_identical(attr(back[[col]], "label"), attr(mem[[col]], "label"))
      }
    }
    expect_true(file.exists(file.path(dirname(path), "truth.json")))
  })
}

test_that("write_submission() writes SDTM and ADaM folders with define.xml", {
  dir <- tempfile()
  paths <- write_submission(par, dir)
  expect_true(file.exists(file.path(dir, "sdtm", "define.xml")))
  expect_true(file.exists(file.path(dir, "adam", "define.xml")))
  expect_setequal(toupper(tools::file_path_sans_ext(list.files(file.path(dir, "sdtm"), "\\.xpt$"))), names(par$sdtm))
  expect_setequal(toupper(tools::file_path_sans_ext(list.files(file.path(dir, "adam"), "\\.xpt$"))), names(par$adam))
  expect_match(define_xml(par, "adam"), "^<\\?xml version=\"1.0\" encoding=\"UTF-8\"\\?>")
  expect_error(define_xml(par, "send"), "should be one of")
})

test_that("study_from_spec() matches study() for the same settings", {
  path <- test_path("..", "..", "..", "..", "examples", "hypertension.yaml")
  skip_if_not(file.exists(path), "examples/ is not available")
  from_file <- study_from_spec(path)
  direct <- study("parallel_continuous", n = 240, seed = 20261013, studyid = "HTN001",
                  indication = "hypertension", weeks = 12, change = c(-4, -11), dropout_rate = 0.12)
  expect_identical(from_file$sdtm$LB, direct$sdtm$LB)
  from_list <- study_from_spec(list(study = list(archetype = "parallel_continuous", seed = 1), population = list(n = 10)))
  expect_identical(nrow(from_list$sdtm$DM), 10L)
  expect_error(study_from_spec(list(study = list(archetype = "parallel_continuous", seed = 1), dropout = list(a = 1))), "not supported yet")
})

test_that("verify() levels", {
  expect_true(verify(par, level = 2)$passed)
  expect_error(verify(par, level = 3), "only in the Python")
  expect_error(verify(par, level = 1), "must be 2")
})

test_that("a phase 1 study has PK data with date-times", {
  pk <- study("phase1_sad_mad", n = 16, seed = 3, part = "mad", doses = c(50, 200), mad_days = 4)
  expect_true(all(c("PC", "PP") %in% names(pk$sdtm)))
  expect_true(all(c("ADPC", "ADPP") %in% names(pk$adam)))
  expect_s3_class(pk$adam$ADPC$ADTM, "POSIXct")
  expect_identical(attr(pk$adam$ADPC$ADTM, "tzone"), "UTC")
  expect_identical(format(pk$adam$ADPC$ADTM[1], "%H:%M"), substr(pk$sdtm$PC$PCDTC[1], 12, 16))
  expect_true(all(c("AUCTAU", "ARAUC") %in% pk$sdtm$PP$PPTESTCD))
  expect_true(verify(pk)$passed)
})

test_that("a crossover study has periods", {
  be <- study("crossover_be", n = 16, seed = 4, design = "food_effect", gmr = 1.3)
  expect_true(all(c("TRT01P", "TRT02P", "AP02SDT", "TRTSEQP") %in% names(be$adam$ADSL)))
  expect_true("APERIOD" %in% names(be$adam$ADPP))
  expect_identical(be$truth$bioequivalence$target$gmr_auc, 1.3)
  expect_true(verify(be)$passed)
})

test_that("dose-ranging and non-inferiority studies report their analyses", {
  dr <- study("dose_ranging_ph2", n = 60, seed = 2, doses = c(0, 5, 20))
  expect_length(unique(dr$adam$ADSL$TRT01P), 3)
  expect_identical(dr$truth$dose_response$target$ed50, 7.5)
  ni <- study("parallel_continuous", n = 40, seed = 2, comparator = "active", ni_margin = 0.3)
  expect_identical(ni$truth$non_inferiority$margin, 0.3)
})

test_that("an Alzheimer's study has ADAS-Cog items and totals", {
  cns <- study("cns_alzheimers", n = 30, seed = 3)
  expect_true(all(cns$sdtm$FT$FTCAT == "ADAS-COG"))
  expect_true("ACTOT" %in% cns$adam$ADADAS$PARAMCD)
  expect_true(verify(cns)$passed)
})

test_that("a cardiovascular outcomes study is event-driven", {
  cv <- study("cv_outcomes", n = 400, seed = 3, target_events = 30)
  expect_setequal(unique(cv$adam$ADTTE$PARAMCD), c("MACE", "CVDEATH", "MI", "STROKE", "DEATH"))
  expect_gte(cv$truth$outcomes$realised_first_mace_events, 30)
  expect_true(verify(cv)$passed)
})
