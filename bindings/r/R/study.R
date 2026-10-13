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
#'
#' @return A `trialforge_study` object: a list with `studyid`, `archetype`, `seed`, `cutoff`
#'   (a `Date`), `controls`, `truth` (targets vs realised values), and `sdtm` and `adam` (named
#'   lists of data frames whose columns carry haven-style `label` attributes).
#' @seealso [write_study()], [verify()]
#' @export
#' @examples
#' trial <- study("parallel_continuous", n = 50, seed = 1)
#' trial
#' head(trial$adam$ADSL)
#' trial$truth$endpoint$realised_difference
study <- function(archetype, n = 300, seed = NULL, studyid = NULL, ...) {
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
  core <- study_impl(
    archetype,
    n,
    if (is.null(seed)) NA_real_ else as.double(seed),
    if (is.null(studyid)) NA_character_ else studyid,
    payload
  )
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
      .core = core
    ),
    class = "trialforge_study"
  )
}

#' Write a study's datasets to files
#'
#' Writes every dataset of a study to the directory `path` (created if needed), plus
#' `truth.json`. Files are named by lower-case dataset name, e.g. `dm.xpt`. This is
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

#' Verify a study
#'
#' Re-runs trialforge's verification levels 0–2 (structural, referential and derivational
#' checks). Generation already guarantees these pass; this is for after the data has been
#' modified.
#'
#' @param study A `trialforge_study` from [study()].
#'
#' @return A list with `passed` (logical) and `findings` (a data frame of failed checks).
#' @export
#' @examples
#' verify(study("parallel_continuous", n = 20, seed = 1))$passed
verify <- function(study) {
  check_study(study)
  study$.core$verify()
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
  invisible(x)
}
