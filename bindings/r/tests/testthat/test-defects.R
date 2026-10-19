spec <- list(ae_end_before_start = 3, duplicate_seq = 3, out_of_ct = 5)
dirty <- study("parallel_continuous", n = 60, seed = 4, defects = spec)

test_that("defects are injected and caught exactly", {
  d <- dirty$defects
  expect_s3_class(d, "data.frame")
  expect_named(d, c("id", "dataset", "row", "paired_row", "usubjid", "key", "variable", "original", "injected", "expected"))
  expect_identical(as.vector(table(d$id)[names(spec)]), c(3L, 3L, 5L))
  report <- verify(dirty)
  expect_false(report$passed)
  s <- score(dirty, report)
  expect_s3_class(s, "trialforge_score")
  expect_identical(c(s$injected, s$detected), c(11L, 11L))
  expect_identical(c(s$recall, s$precision), c(1, 1))
  expect_identical(s$missed, integer(0))
  expect_identical(nrow(s$unmatched), 0L)
  expect_setequal(s$by_type$id, names(spec))
  expect_output(print(s), "recall 100.0%")
  expect_output(print(dirty), "defects: 11")
})

test_that("defect rows count from 1, like verify()", {
  s <- study("parallel_continuous", n = 60, seed = 4, defects = list(out_of_ct = 4))
  for (k in seq_len(nrow(s$defects))) {
    d <- s$defects[k, ]
    expect_identical(s$sdtm[[d$dataset]][[d$variable]][d$row], d$injected)
    expect_identical(d$expected, sprintf("L0.CT %s %d", d$dataset, d$row))
  }
  findings <- verify(s)$findings
  expect_setequal(findings$row, s$defects$row)
  # Dropping a finding leaves its defect missed.
  partial <- list(passed = FALSE, findings = findings[-1, ])
  missed <- score(s, partial)$missed
  expect_identical(s$defects$row[missed], findings$row[1])
})

test_that("a clean study has no defects; defects.json is zero-based", {
  clean <- study("parallel_continuous", n = 40, seed = 2)
  expect_null(clean$defects)
  expect_error(score(clean, verify(clean)), "no injected defects")
  expect_null(study("parallel_continuous", n = 40, seed = 2, defects = "none")$defects)
  s <- study("parallel_continuous", n = 40, seed = 2, defects = c(out_of_ct = 2))
  expect_identical(s$truth, clean$truth)
  dir <- tempfile()
  paths <- write_study(s, dir, "csv")
  expect_true(any(basename(paths) == "defects.json"))
  manifest <- jsonlite::fromJSON(file.path(dir, "defects.json"))
  expect_identical(manifest$rows, "zero-based")
  expect_equal(manifest$defects$row + 1L, s$defects$row, ignore_attr = TRUE)
})

test_that("rates, spec files and bad requests", {
  s <- study("parallel_continuous", n = 60, seed = 7, defects = list(out_of_ct = 0.05))
  expect_gt(nrow(s$defects), 0)
  from_spec <- study_from_spec(list(
    study = list(archetype = "parallel_continuous", seed = 7),
    population = list(n = 60),
    data_quality = list(defects = list(out_of_ct = 0.05))
  ))
  expect_identical(from_spec$defects, s$defects)
  expect_error(study("parallel_continuous", n = 10, seed = 1, defects = list(nope = 1)), "unknown defect")
  expect_error(study("parallel_continuous", n = 10, seed = 1, defects = "realistic"), "not available yet")
})
