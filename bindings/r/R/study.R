#' Generate a mock study
#'
#' Generates a verified mock clinical study: trial design, SDTM and ADaM datasets. Every study
#' passes trialforge's structural, referential and derivational checks before it is returned.
#'
#' @param archetype The study type, e.g. `"parallel_continuous"`. See [archetypes()].
#' @param n Number of randomized subjects.
#' @param seed Root seed, a whole number in `[0, 2^53]`. With `NULL`, a fresh seed is drawn
#'   and recorded as `$seed` on the result, so the study can always be reproduced.
#' @param studyid `STUDYID`. Each archetype has a default.
#' @param ... Archetype controls, e.g. `indication = "hypertension", weeks = 12`. See
#'   [controls()]. Unknown names are an error.
#' @param defects Data-quality defects to inject after the clean study is generated and
#'   verified: `NULL` or `"none"` for none, or a named list of defect id to a count (a whole
#'   number) or a rate per eligible record (a number strictly between 0 and 1), e.g.
#'   `list(ae_end_before_start = 3, out_of_ct = 0.02)`. Whole numbers are counts, so ask for a
#'   large count to corrupt every eligible record. See [score()].
#'
#' @return A `trialforge_study` object: a list with `studyid`, `archetype`, `seed`, `cutoff`
#'   (a `Date`), `controls`, `truth` (targets vs realised values), `sdtm` and `adam` (named
#'   lists of data frames whose columns carry haven-style `label` attributes), and `defects`:
#'   `NULL`, or a data frame with one row per injected defect (`id`, `dataset`, `row`,
#'   `paired_row` (another record the defect involves, which a validator may name instead),
#'   `usubjid`, `key` (the record's key before the edit, e.g. `AESEQ=3`), `variable`,
#'   `original`, `injected` and `expected` (the findings it should produce, as
#'   `RULE DATASET row`, separated by `; `)). Rows count from 1, as in [verify()].
#' @seealso [write_study()], [verify()], [score()]
#' @export
#' @examples
#' trial <- study("parallel_continuous", n = 50, seed = 1)
#' trial
#' head(trial$adam$ADSL)
#' trial$truth$endpoint$realised_difference
#'
#' dirty <- study("parallel_continuous", n = 50, seed = 1, defects = list(out_of_ct = 3))
#' dirty$defects[, c("dataset", "row", "variable", "original", "injected")]
study <- function(archetype, n = 300, seed = NULL, studyid = NULL, ..., defects = NULL) {
  if (!is.null(studyid) && !(is.character(studyid) && length(studyid) == 1)) {
    stop("`studyid` must be a single string", call. = FALSE)
  }
  controls <- list(...)
  if (length(controls) > 0 && (is.null(names(controls)) || any(names(controls) == ""))) {
    stop("controls must be named, e.g. `weeks = 12`", call. = FALSE)
  }
  payload <- if (length(controls) == 0) {
    NA_character_
  } else {
    as.character(jsonlite::toJSON(controls, auto_unbox = TRUE, null = "null", digits = NA))
  }
  requested <- if (is.null(defects)) {
    NA_character_
  } else {
    if (is.atomic(defects) && !is.null(names(defects))) defects <- as.list(defects)
    as.character(jsonlite::toJSON(defects, auto_unbox = TRUE, digits = NA))
  }
  core <- study_impl(
    archetype,
    n,
    if (is.null(seed)) NA_real_ else as.double(seed),
    if (is.null(studyid)) NA_character_ else studyid,
    payload,
    requested
  )
  new_study(core)
}

#' Generate a study from a study spec
#'
#' @param spec A path to a JSON or YAML study spec, the spec text itself, or a list with the
#'   spec's structure. The spec's `outputs` section is for the command-line tool; write files
#'   with [write_study()] or [write_submission()].
#'
#' @return A `trialforge_study`, as from [study()].
#' @export
#' @examples
#' spec <- list(
#'   study = list(archetype = "parallel_continuous", seed = 1),
#'   population = list(n = 30),
#'   controls = list(indication = "obesity")
#' )
#' study_from_spec(spec)
study_from_spec <- function(spec) {
  text <- if (is.list(spec)) {
    as.character(jsonlite::toJSON(spec, auto_unbox = TRUE, null = "null", digits = NA))
  } else if (length(spec) == 1 && !grepl("\n", spec) && file.exists(spec)) {
    paste(readLines(spec, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  } else {
    paste(spec, collapse = "\n")
  }
  new_study(study_from_spec_impl(text))
}

new_study <- function(core) {
  tables <- function(names) stats::setNames(lapply(names, core$table), names)
  structure(
    list(
      studyid = core$studyid(),
      archetype = core$archetype(),
      seed = core$seed(),
      cutoff = as.Date(core$cutoff_days(), origin = "1970-01-01"),
      controls = jsonlite::fromJSON(core$controls_json(), simplifyVector = TRUE),
      truth = jsonlite::fromJSON(core$truth_json(), simplifyVector = TRUE),
      sdtm = tables(core$sdtm_names()),
      adam = tables(core$adam_names()),
      defects = core$defects_table(),
      .core = core
    ),
    class = "trialforge_study"
  )
}

#' Write a study's datasets to files
#'
#' Writes every dataset of a study to the directory `path` (created if needed), plus
#' `truth.json` (and `defects.json` for a study with defects). Files are named by lower-case dataset name, e.g. `dm.xpt`. This is
#' `Study.write()` in Python; the name differs because `write()` would mask `base::write()`.
#'
#' @param study A `trialforge_study` from [study()].
#' @param path Output directory.
#' @param formats Output formats; see [formats()].
#'
#' @return The paths written, invisibly.
#' @export
#' @examples
#' trial <- study("parallel_continuous", n = 20, seed = 1)
#' write_study(trial, tempfile(), formats = c("xpt", "csv"))
write_study <- function(study, path, formats = "xpt") {
  check_study(study)
  invisible(study$.core$write(path.expand(path), as.character(formats)))
}

#' Write a submission-style package
#'
#' Writes `sdtm/` and `adam/` folders to `path`, each with one XPT file per dataset and a
#' Define-XML 2.1 `define.xml`, plus `truth.json` (and `defects.json`).
#'
#' @param study A `trialforge_study` from [study()].
#' @param path Output directory.
#'
#' @return The paths written, invisibly.
#' @export
#' @examples
#' trial <- study("parallel_continuous", n = 20, seed = 1)
#' write_submission(trial, tempfile())
write_submission <- function(study, path) {
  check_study(study)
  invisible(study$.core$write_submission(path.expand(path)))
}

#' Define-XML for a study
#'
#' @param study A `trialforge_study` from [study()].
#' @param standard `"sdtm"` or `"adam"`.
#'
#' @return The Define-XML 2.1 document as a single string.
#' @export
#' @examples
#' xml <- define_xml(study("parallel_continuous", n = 20, seed = 1), "adam")
#' substr(xml, 1, 200)
define_xml <- function(study, standard = c("sdtm", "adam")) {
  check_study(study)
  study$.core$define_xml(match.arg(standard))
}

#' Verify a study
#'
#' Re-runs trialforge's verification levels 0–2 (structural, referential and derivational
#' checks). Generation already guarantees these pass; this is for after the data has been
#' modified.
#'
#' Level 3 (CDISC conformance rules) runs Pointblank's SDTMIG engine, which exists only in the
#' Python package, so it is not available from R: use `Study.verify(level=3)` in Python.
#'
#' @param study A `trialforge_study` from [study()].
#' @param level Verification level; only `2` is available in R.
#'
#' @return A list with `passed` (logical) and `findings` (a data frame of failed checks).
#' @export
#' @examples
#' verify(study("parallel_continuous", n = 20, seed = 1))$passed
verify <- function(study, level = 2) {
  check_study(study)
  if (identical(level, 3) || identical(level, 3L)) {
    stop(
      "Level 3 runs Pointblank's CDISC conformance engine, which is only in the Python ",
      "package. Use `Study.verify(level=3)` from Python.",
      call. = FALSE
    )
  }
  if (!(identical(level, 2) || identical(level, 2L))) {
    stop("`level` must be 2 (or 3, from Python)", call. = FALSE)
  }
  study$.core$verify()
}

#' Score findings against a study's injected defects
#'
#' Matches the findings of a report to the defects injected by `study(defects = ...)`, for
#' recall and precision per defect type. A finding matches a defect when it names the dataset
#' and row the defect expects a finding at and, for trialforge's own levels (L0–L2), the
#' expected rule; findings at any other level (another validator's) match on dataset and row
#' only, or on the defect's `paired_row`. This is `Study.score()` in Python.
#'
#' @param study A `trialforge_study` from [study()] with defects.
#' @param report A report from [verify()]: a list with a `findings` data frame with `level`,
#'   `rule`, `dataset` and `row` (1-based, `NA` for none) columns.
#'
#' @return A `trialforge_score`: a list with `injected`, `detected`, `recall`
#'   (`detected / injected`), `findings`, `matched`, `precision` (`matched / findings`; `NA`
#'   when undefined), `by_type` (a data frame with the same per defect `id`, where
#'   `false_positives` counts unmatched L0–L2 findings from a rule that catches that type),
#'   `missed` (rows of `study$defects` no finding matched) and `unmatched` (the findings no
#'   defect explains).
#' @export
#' @examples
#' trial <- study("parallel_continuous", n = 50, seed = 1, defects = list(out_of_ct = 3))
#' score(trial, verify(trial))
score <- function(study, report) {
  check_study(study)
  f <- report$findings
  if (!is.data.frame(f) || !all(c("level", "rule", "dataset", "row") %in% names(f))) {
    stop("`report` must come from `verify()`", call. = FALSE)
  }
  raw <- jsonlite::fromJSON(
    study$.core$score(
      as.character(f$level),
      as.character(f$rule),
      as.character(f$dataset),
      as.integer(f$row)
    ),
    simplifyVector = FALSE
  )
  num <- function(x) if (is.null(x)) NA_real_ else x
  types <- raw$by_type
  field <- function(name, na = NA_real_) {
    vapply(types, function(t) if (is.null(t[[name]])) na else as.double(t[[name]]), double(1))
  }
  by_type <- data.frame(
    id = as.character(names(types)),
    injected = as.integer(field("injected")),
    detected = as.integer(field("detected")),
    recall = field("recall"),
    findings = as.integer(field("findings")),
    false_positives = as.integer(field("false_positives")),
    precision = field("precision"),
    stringsAsFactors = FALSE
  )
  rownames(by_type) <- NULL
  unmatched <- f[as.integer(unlist(raw$unmatched)) + 1L, , drop = FALSE]
  rownames(unmatched) <- NULL
  structure(
    list(
      injected = raw$injected,
      detected = raw$detected,
      recall = num(raw$recall),
      findings = raw$findings,
      matched = raw$matched,
      precision = num(raw$precision),
      by_type = by_type,
      missed = as.integer(unlist(raw$missed)) + 1L,
      unmatched = unmatched
    ),
    class = "trialforge_score"
  )
}

#' @export
print.trialforge_score <- function(x, ...) {
  pct <- function(p) if (is.na(p)) "-" else sprintf("%.1f%%", 100 * p)
  cat(sprintf(
    "<trialforge_score> recall %s (%d of %d defects), precision %s (%d of %d findings)\n",
    pct(x$recall), x$detected, x$injected, pct(x$precision), x$matched, x$findings
  ))
  print(x$by_type, row.names = FALSE)
  invisible(x)
}

#' Available study archetypes
#'
#' @return A named character vector: archetype names with one-line descriptions.
#' @export
#' @examples
#' archetypes()
archetypes <- function() {
  a <- archetypes_impl()
  stats::setNames(a$description, a$name)
}

#' Controls of an archetype
#'
#' @param archetype An archetype name; see [archetypes()].
#'
#' @return A data frame with `name`, `description` and `default` (a list column).
#' @export
#' @examples
#' controls("parallel_continuous")
controls <- function(archetype) {
  docs <- jsonlite::fromJSON(controls_impl(archetype), simplifyVector = FALSE)
  data.frame(
    name = vapply(docs, `[[`, character(1), "name"),
    description = vapply(docs, `[[`, character(1), "description"),
    default = I(lapply(docs, function(d) if (is.null(d$default)) NULL else unlist(d$default))),
    stringsAsFactors = FALSE
  )
}

#' Output formats
#'
#' @return The formats [write_study()] supports in this build.
#' @export
#' @examples
#' formats()
formats <- function() {
  formats_impl()
}

check_study <- function(study) {
  if (!inherits(study, "trialforge_study")) {
    stop("`study` must be a trialforge_study from `study()`", call. = FALSE)
  }
}

#' @export
print.trialforge_study <- function(x, ...) {
  rows <- vapply(c(x$sdtm, x$adam), nrow, integer(1))
  fmt <- function(names) paste(sprintf("%s (%d)", names, rows[names]), collapse = ", ")
  cat(sprintf("<trialforge_study> %s (%s)\n", x$studyid, x$archetype))
  cat(sprintf("  seed: %.0f   cutoff: %s\n", x$seed, format(x$cutoff)))
  cat(strwrap(paste("SDTM:", fmt(names(x$sdtm))), indent = 2, exdent = 8), sep = "\n")
  cat(strwrap(paste("ADaM:", fmt(names(x$adam))), indent = 2, exdent = 8), sep = "\n")
  if (!is.null(x$defects)) {
    counts <- table(x$defects$id)
    cat(sprintf(
      "  defects: %d (%s)\n",
      nrow(x$defects), paste(names(counts), counts, collapse = ", ")
    ))
  }
  invisible(x)
}
