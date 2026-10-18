# PKNCA equivalence: recompute every PP parameter from trialforge's PC and EX with PKNCA's
# pk.nca() (default options) and compare value by value.
#
# Usage (from the repo root, via ir so nothing touches the user library):
#   ir run --isolated --r-version 4.6 --exclude-newer 2026-09-25 \
#     --with PKNCA --with haven --with dplyr scripts/pknca_equivalence.R <submission-dir>
#
# Profiles are built as trialforge documents them: one per reference dose (PCRFTDTC), times in
# hours after the dose (pre-dose as 0), BLQ as 0, and after the first dose the previous dose's
# 24 h sample as the concentration at time 0. AUCTAU is PKNCA's aucint.last over [0, tau]
# (interpolated at tau from all the data). Exits non-zero if any parameter differs.

suppressPackageStartupMessages({
  library(PKNCA)
  library(dplyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("usage: pknca_equivalence.R <submission-dir>")
read <- function(name) {
  df <- haven::read_xpt(file.path(args[[1]], "sdtm", paste0(name, ".xpt")))
  df |> mutate(across(where(is.character), ~ na_if(.x, "")))
}
pc <- read("pc")
ex <- read("ex")
pp <- read("pp")
tau <- 24
dtm <- function(x) as.POSIXct(x, format = if (nchar(x[1]) > 16) "%Y-%m-%dT%H:%M:%S" else "%Y-%m-%dT%H:%M", tz = "UTC")

# Every profile is a single dose when the regimen is ONCE (single-dose and crossover studies);
# otherwise the first and last of several repeated doses are profiled.
doses <- ex |>
  filter(EXDOSE > 0) |>
  mutate(REF = dtm(EXSTDTC)) |>
  group_by(USUBJID) |>
  arrange(REF, .by_group = TRUE) |>
  mutate(K = row_number(), N = n(), SINGLE = all(EXDOSFRQ == "ONCE") | n() == 1) |>
  ungroup() |>
  select(USUBJID, REF, K, N, SINGLE, DOSE = EXDOSE)
samples <- pc |>
  mutate(REF = dtm(PCRFTDTC), TIME = pmax(as.numeric(difftime(dtm(PCDTC), REF, units = "hours")), 0),
         CONC = coalesce(PCSTRESN, 0)) |>
  inner_join(doses, by = c("USUBJID", "REF"))
troughs <- samples |>
  filter(PCTPTNUM == tau, !SINGLE) |>
  transmute(USUBJID, K = K + 1, TIME = 0, CONC)
profiles <- samples |>
  group_by(USUBJID, K) |>
  filter(n() >= 3) |>
  ungroup() |>
  select(USUBJID, K, N, SINGLE, DOSE, REF, TIME, CONC)
profiles <- bind_rows(
  profiles,
  troughs |> inner_join(distinct(profiles, USUBJID, K, N, SINGLE, DOSE, REF), by = c("USUBJID", "K"))
) |>
  arrange(USUBJID, K, TIME) |>
  mutate(PROFILE = paste(USUBJID, K))

conc <- PKNCAconc(profiles, CONC ~ TIME | PROFILE)
dose <- PKNCAdose(distinct(profiles, PROFILE, DOSE) |> mutate(TIME = 0), DOSE ~ TIME | PROFILE)
single <- c(cmax = TRUE, tmax = TRUE, clast.obs = TRUE, tlast = TRUE, auclast = TRUE, aucall = TRUE,
            half.life = TRUE, lambda.z = TRUE, lambda.z.n.points = TRUE, adj.r.squared = TRUE,
            aucinf.obs = TRUE, aucpext.obs = TRUE, cl.obs = TRUE, vz.obs = TRUE, cmax.dn = TRUE,
            aucinf.obs.dn = TRUE)
kinds <- distinct(profiles, PROFILE, K, N, SINGLE)
iv <- function(start, end, params, which) {
  ids <- kinds$PROFILE[which]
  if (!length(ids)) return(NULL)
  row <- as.data.frame(as.list(params))
  cbind(data.frame(PROFILE = ids, start = start, end = end), row[rep(1, length(ids)), , drop = FALSE])
}
intervals <- bind_rows(
  iv(0, Inf, single, kinds$SINGLE),
  # AUC over the dosing interval interpolates the concentration at tau from all the data.
  iv(0, tau, c(cmax = TRUE, tmax = TRUE, aucint.last = TRUE), !kinds$SINGLE),
  iv(0, Inf, c(half.life = TRUE, lambda.z = TRUE, lambda.z.n.points = TRUE, adj.r.squared = TRUE),
     !kinds$SINGLE & kinds$K == kinds$N)
)
intervals[is.na(intervals)] <- FALSE
res <- as.data.frame(pk.nca(PKNCAdata(conc, dose, intervals = intervals)))

# Map PKNCA results to PPTESTCD (doses in mg, AUC in h*ng/mL, so CL = 1000 * dose / AUC).
code <- c(cmax = "CMAX", tmax = "TMAX", clast.obs = "CLST", tlast = "TLST", auclast = "AUCLST",
          aucall = "AUCALL", lambda.z = "LAMZ", half.life = "LAMZHL", lambda.z.n.points = "LAMZNPT",
          adj.r.squared = "R2ADJ", aucinf.obs = "AUCIFO", aucpext.obs = "AUCPEO", cl.obs = "CLFO",
          vz.obs = "VZFO", cmax.dn = "CMAXD", aucinf.obs.dn = "AUCIFOD")
res <- res |>
  inner_join(kinds, by = "PROFILE") |>
  inner_join(distinct(profiles, PROFILE, USUBJID, DOSE, REF), by = "PROFILE") |>
  mutate(
    multiple = !SINGLE,
    PPTESTCD = case_when(PPTESTCD == "aucint.last" ~ "AUCTAU", TRUE ~ unname(code[PPTESTCD])),
    PPSTRESN = case_when(PPTESTCD %in% c("CLFO", "VZFO") ~ PPORRES * 1000, TRUE ~ PPORRES)
  ) |>
  filter(!is.na(PPTESTCD), !is.na(PPSTRESN)) |>
  # Only the parameters requested per interval (PKNCA adds dependencies such as tlast).
  filter(!multiple | (is.finite(end) & PPTESTCD %in% c("CMAX", "TMAX", "AUCTAU")) |
           (!is.finite(end) & PPTESTCD %in% c("LAMZ", "LAMZHL", "LAMZNPT", "R2ADJ")))
# Derived from PKNCA's results: CL/F over tau and the accumulation ratios.
tau_rows <- res |> filter(multiple, PPTESTCD %in% c("AUCTAU", "CMAX"))
extra <- bind_rows(
  tau_rows |> filter(K == N, PPTESTCD == "AUCTAU") |> mutate(PPTESTCD = "CLFTAU", PPSTRESN = 1000 * DOSE / PPSTRESN),
  tau_rows |>
    group_by(USUBJID, PPTESTCD) |>
    filter(any(K == 1), any(K == N)) |>
    summarise(REF = REF[K == N], PPSTRESN = PPSTRESN[K == N] / PPSTRESN[K == 1], .groups = "drop") |>
    mutate(PPTESTCD = c(AUCTAU = "ARAUC", CMAX = "ARCMAX")[PPTESTCD])
)
# Adaptation (also in crates/trialforge-sdtm/src/pk.rs): a ratio with a zero denominator (an
# all-BLQ interval) is Inf here and left out of PP, which can't hold it.
extra <- extra |> filter(is.finite(PPSTRESN))
pknca <- bind_rows(res |> select(USUBJID, REF, PPTESTCD, PPSTRESN), extra |> select(USUBJID, REF, PPTESTCD, PPSTRESN))

ours <- pp |> mutate(REF = dtm(PPRFTDTC)) |> select(USUBJID, REF, PPTESTCD, OURS = PPSTRESN)
joined <- full_join(ours, rename(pknca, PKNCA = PPSTRESN), by = c("USUBJID", "REF", "PPTESTCD"))
joined$diff <- abs(joined$OURS - joined$PKNCA) / pmax(1, abs(joined$PKNCA))
bad <- joined |> filter(is.na(OURS) | is.na(PKNCA) | diff > 1e-9)
cat(sprintf("PKNCA %s: compared %d parameters over %d profiles; %d differ\n",
            as.character(packageVersion("PKNCA")), nrow(joined), n_distinct(profiles$PROFILE), nrow(bad)))
if (nrow(bad) > 0) {
  print(utils::head(as.data.frame(bad), 20))
  quit(status = 1)
}
