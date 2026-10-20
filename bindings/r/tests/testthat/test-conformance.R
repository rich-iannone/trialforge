# Replay the shared conformance cases in spec/conformance/ (also replayed by the Python suite).
# The spec directory is outside the package, so these tests skip when the package is checked
# from a built tarball.

# MD5 of a string, via a temporary file so only base R is needed.
digest_md5 <- function(text) {
  f <- tempfile()
  on.exit(unlink(f))
  writeBin(charToRaw(enc2utf8(text)), f)
  unname(tools::md5sum(f))
}

spec_path <- function(name) test_path("..", "..", "..", "..", "spec", "conformance", name)

test_that("RNG stream seeds match the shared golden cases", {
  path <- spec_path("rng_streams.json")
  skip_if_not(file.exists(path), "spec/conformance is not available")
  cases <- jsonlite::read_json(path)$cases
  for (case in cases) {
    expect_identical(stream_seed(case$seed, case$path), case$expected)
  }
})

test_that("written XPT files match the shared golden hashes", {
  path <- spec_path("study_xpt.json")
  skip_if_not(file.exists(path), "spec/conformance is not available")
  for (case in jsonlite::read_json(path)$cases) {
    trial <- do.call(study, c(list(case$archetype, n = case$n, seed = case$seed), case$controls))
    rows <- vapply(c(trial$sdtm, trial$adam), nrow, integer(1))
    expect_identical(as.list(rows), lapply(case$rows, as.integer))
    paths <- write_study(trial, tempfile())
    paths <- paths[grepl("\\.xpt$", paths)]
    md5 <- as.list(stats::setNames(unname(tools::md5sum(paths)), basename(paths)))
    expect_identical(md5[names(case$xpt_md5)], case$xpt_md5)
    for (std in c("sdtm", "adam")) {
      xml_md5 <- digest_md5(define_xml(trial, std))
      expect_identical(xml_md5, case$define_md5[[std]], info = paste(case$archetype, std))
    }
  }
})

test_that("injected defects match the shared golden manifests", {
  path <- spec_path("defects.json")
  skip_if_not(file.exists(path), "spec/conformance is not available")
  # One field of every manifest defect, NULL as NA.
  pluck <- function(defects, name, type) {
    vapply(defects, function(x) if (is.null(x[[name]])) type[NA] else x[[name]], type)
  }
  for (case in jsonlite::read_json(path)$cases) {
    trial <- do.call(study, c(
      list(case$archetype, n = case$n, seed = case$seed, defects = case$defects),
      case$controls
    ))
    dir <- tempfile()
    paths <- write_study(trial, dir)
    defects <- case$manifest$defects
    expect_identical(jsonlite::read_json(file.path(dir, "defects.json")), case$manifest)
    md5 <- as.list(stats::setNames(unname(tools::md5sum(paths)), basename(paths)))
    expect_identical(md5[names(case$xpt_md5)], case$xpt_md5, info = case$archetype)
    # The defect table holds the manifest's fields, with rows counted from 1.
    table <- lapply(trial$defects, as.vector)
    for (name in c("id", "dataset", "usubjid", "key", "variable", "original", "injected")) {
      expect_identical(table[[name]], pluck(defects, name, character(1)), info = name)
    }
    for (name in c("row", "paired_row")) {
      expect_identical(table[[name]], pluck(defects, name, integer(1)) + 1L, info = name)
    }
    expected <- vapply(defects, function(x) {
      at <- vapply(x$expected, function(e) sprintf("%s %s %d", e$rule, e$dataset, e$row + 1L), "")
      paste(at, collapse = "; ")
    }, "")
    expect_identical(table$expected, expected)
    s <- score(trial, verify(trial))
    expect_identical(lapply(s[names(case$score)], as.integer), lapply(case$score, as.integer))
    expect_identical(c(s$recall, s$precision), c(1, 1))
  }
})
