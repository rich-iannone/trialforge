#' Generate a mock study
#'
#' Generates the SDTM and ADaM datasets of a mock clinical study. This preview supports the
#' `"oncology_solid_tumor"` archetype, which produces DM, EX, DS and ADSL.
#'
#' @param archetype The study type. See [archetypes()].
#' @param n Number of randomized subjects.
#' @param seed Root seed, a whole number in `[0, 2^53]`. With `NULL`, a fresh seed is drawn
#'   and recorded as `$seed` on the result, so the study can always be reproduced.
#' @param studyid `STUDYID`. Each archetype has a default.
#'
#' @return A `trialforge_study` object: a list with `studyid`, `archetype`, `seed`, `cutoff`
#'   (a `Date`), and `sdtm` and `adam` (named lists of data frames whose columns carry
#'   haven-style `label` attributes).
#' @seealso [write_study()]
#' @export
#' @examples
#' trial <- study("oncology_solid_tumor", n = 50, seed = 1)
#' trial
#' head(trial$adam$ADSL)
study <- function(archetype, n = 300, seed = NULL, studyid = NULL) {
  core <- study_impl(
    archetype,
    n,
    if (is.null(seed)) NA_real_ else as.double(seed),
    if (is.null(studyid)) NA_character_ else studyid
  )
  tables <- function(names) stats::setNames(lapply(names, core$table), names)
  structure(
    list(
      studyid = core$studyid(),
      archetype = core$archetype(),
      seed = core$seed(),
      cutoff = as.Date(core$cutoff_days(), origin = "1970-01-01"),
      sdtm = tables(core$sdtm_names()),
      adam = tables(core$adam_names()),
      .core = core
    ),
    class = "trialforge_study"
  )
}

#' Write a study's datasets to files
#'
#' Writes every dataset of a study to the directory `path` (created if needed). Files are named
#' by lower-case dataset name, e.g. `dm.xpt`. This is `Study.write()` in Python; the name
#' differs because `write()` would mask `base::write()`.
#'
#' @param study A `trialforge_study` from [study()].
#' @param path Output directory.
#' @param formats Output formats. Currently only `"xpt"` (SAS transport v5).
#'
#' @return The paths written, invisibly.
#' @export
#' @examples
#' trial <- study("oncology_solid_tumor", n = 20, seed = 1)
#' write_study(trial, tempfile())
write_study <- function(study, path, formats = "xpt") {
  if (!inherits(study, "trialforge_study")) {
    stop("`study` must be a trialforge_study from `study()`", call. = FALSE)
  }
  invisible(study$.core$write(path.expand(path), as.character(formats)))
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

#' @export
print.trialforge_study <- function(x, ...) {
  rows <- vapply(c(x$sdtm, x$adam), nrow, integer(1))
  cat(sprintf("<trialforge_study> %s (%s)\n", x$studyid, x$archetype))
  cat(sprintf("  seed: %.0f   cutoff: %s\n", x$seed, format(x$cutoff)))
  cat(sprintf("  SDTM: %s\n", paste(sprintf("%s (%d)", names(x$sdtm), rows[names(x$sdtm)]), collapse = ", ")))
  cat(sprintf("  ADaM: %s\n", paste(sprintf("%s (%d)", names(x$adam), rows[names(x$adam)]), collapse = ", ")))
  invisible(x)
}
