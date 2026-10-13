onc <- study("oncology_solid_tumor", n = 120, seed = 3)

test_that("archetypes() lists the oncology archetype", {
  expect_true("oncology_solid_tumor" %in% names(archetypes()))
})

test_that("study() returns the expected structure", {
  expect_s3_class(onc, "trialforge_study")
  expect_identical(onc$studyid, "SIMONC02")
  expect_identical(onc$seed, 3)
  expect_identical(onc$cutoff, as.Date("2025-08-29"))
  expect_named(onc$sdtm, c("DM", "EX", "DS"))
  expect_named(onc$adam, "ADSL")
  expect_identical(nrow(onc$sdtm$DM), 120L)
  expect_s3_class(onc$adam$ADSL$TRTSDT, "Date")
  expect_type(onc$sdtm$DM$AGE, "integer")
  expect_identical(attr(onc$sdtm$DM$USUBJID, "label"), "Unique Subject Identifier")
  expect_identical(attr(onc$adam$ADSL, "label"), "Subject-Level Analysis Dataset")
  expect_output(print(onc), "SIMONC02")
})

test_that("a missing seed is drawn and recorded", {
  a <- study("oncology_solid_tumor", n = 10)
  expect_true(a$seed >= 0 && a$seed <= 2^53)
  b <- study("oncology_solid_tumor", n = 10, seed = a$seed)
  expect_identical(a$sdtm, b$sdtm)
  expect_no_error(study("oncology_solid_tumor", n = 1, seed = 2^53))
  expect_error(study("oncology_solid_tumor", n = 1, seed = 2^53 + 2), "2\\^53")
})

test_that("bad input is rejected", {
  expect_error(study("nope", n = 10, seed = 1), "unknown archetype")
  expect_error(study("oncology_solid_tumor", n = 0, seed = 1), "at least 1")
  expect_error(write_study(onc, tempfile(), formats = "sas7bdat"), "unknown format")
  expect_error(write_study(list(), tempfile()), "trialforge_study")
})

test_that("XPT files read back with haven match the in-memory data", {
  skip_if_not_installed("haven")
  dir <- tempfile()
  paths <- write_study(onc, dir)
  expect_setequal(basename(paths), c("dm.xpt", "ex.xpt", "ds.xpt", "adsl.xpt"))
  frames <- c(onc$sdtm, onc$adam)
  for (path in paths) {
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
})
