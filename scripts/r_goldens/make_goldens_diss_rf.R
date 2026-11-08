# Live-R goldens for dissrf / seqrf (RP4 T13b + the star seqrf).
#
# Relative-frequency groups (grp.meth="prop") on the phase2 15x15 distance
# matrix. The sort key is R's cmdscale 1-D embedding, dumped as a golden input so
# the Python side uses the identical sortv (avoiding MDS eigenvector-sign drift).
# The seqrf goldens (sequence inputs) are appended in the seqrf commit.
#
# Requires R 4.6.2 + TraMineR 2.2-14.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # phase2_inputs(), write_phase2_csv()
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

inp <- phase2_inputs()
D <- inp$D
w <- inp$wreal

sortv <- as.numeric(cmdscale(D, k = 1))
write_phase2_csv("dissrf_input_sortv", matrix(sortv, nrow = 1), "# dissrf sort key (cmdscale k=1)\n")

gen <- function(tag, weights) {
  rf <- suppressMessages(dissrf(D, k = 3, sortv = sortv, weights = weights, grp.meth = "prop"))
  write_phase2_csv(
    paste0("dissrf_", tag, "_medoids"), matrix(as.vector(rf$medoids) - 1L, nrow = 1),
    "# dissrf group medoids (0-based)\n"
  )
  write_phase2_csv(
    paste0("dissrf_", tag, "_stats"),
    matrix(c(rf$R2, rf$Fstat, rf$pvalue), nrow = 1), "# dissrf R2, Fstat, pvalue\n"
  )
  write_phase2_csv(
    paste0("dissrf_", tag, "_heights"), matrix(rf$heights, nrow = 1), "# dissrf heights\n"
  )
}

gen("uw", NULL)
gen("w", w)

# TraMineR 2.2-14 defines k=1 directly as the overall medoid after normalizing
# supplied weights to n cases. Pin the new release behavior explicitly.
rf1 <- suppressMessages(dissrf(D, k = 1, sortv = sortv, weights = w, grp.meth = "prop"))
write_phase2_csv(
  "dissrf_w_k1_medoid", matrix(as.vector(rf1$medoids) - 1L, nrow = 1),
  "# dissrf weighted k=1 overall medoid (0-based)\n"
)
write_phase2_csv(
  "dissrf_w_k1_sizes", matrix(c(rf1$sizes[["wsum"]], rf1$sizes[["gsize"]]), nrow = 1),
  "# dissrf weighted k=1 normalized wsum and gsize\n"
)

cat("DISS_RF_GOLDENS_DONE\n")

## seqrf goldens (the star item): seqrf wraps dissrf on a sequence diss matrix ----
set.seed(515151)
nrf <- 15
Lrf <- 6
Mrf <- matrix(sample(1:4, nrf * Lrf, replace = TRUE), nrow = nrf, ncol = Lrf)
wrf <- 0.5 + 0.1 * (1:nrf)
sdrf <- suppressMessages(seqdef(Mrf, alphabet = 1:4, weights = wrf))
Drf <- suppressMessages(seqdist(sdrf, method = "LCS", norm = "none"))
sortvrf <- as.numeric(cmdscale(Drf, k = 1))

write_phase2_csv("seqrf_input_codes", Mrf - 1L, "# seqrf input: 0-based codes (alphabet 1:4)\n")
write_phase2_csv("seqrf_input_diss", Drf, "# seqrf input: LCS norm=none distance matrix\n")
write_phase2_csv("seqrf_input_weights", matrix(wrf, nrow = 1), "# seqrf input weights\n")
write_phase2_csv("seqrf_input_sortv", matrix(sortvrf, nrow = 1), "# seqrf sort key (cmdscale k=1)\n")

genrf <- function(tag, weighted) {
  srf <- suppressMessages(
    seqrf(sdrf, diss = Drf, k = 3, sortv = sortvrf, weighted = weighted)
  )
  write_phase2_csv(
    paste0("seqrf_", tag, "_medoids"), matrix(as.vector(srf$rf$medoids) - 1L, nrow = 1),
    "# seqrf group medoids (0-based)\n"
  )
  write_phase2_csv(
    paste0("seqrf_", tag, "_heights"), matrix(srf$rf$heights, nrow = 1), "# seqrf heights\n"
  )
  write_phase2_csv(
    paste0("seqrf_", tag, "_repcodes"), TraMineR:::seqasnum(srf$seqtoplot),
    "# seqrf representative sequence codes (0-based)\n"
  )
}

genrf("uw", FALSE)
genrf("w", TRUE)

## seqrf right-panel box stats (the star plot.seqrf "diss.to.med" boxplot): R draws
## wtd.boxplot.tmr(dist.list, weights.list) -- a WEIGHTED Tukey-hinge box. Pin the
## 5-number summary per group [whislo, q1, med, q3, whishi] so seqrfplot's right box
## is structure-tested against the box R actually draws (not matplotlib's unweighted
## np.percentile). Audit 2026-06-21.
genrf_box <- function(tag, weighted) {
  srf <- suppressMessages(seqrf(sdrf, diss = Drf, k = 3, sortv = sortvrf, weighted = weighted))
  dl <- srf$rf$dist.list
  wl <- srf$rf$weights.list
  bx <- t(vapply(
    seq_along(dl),
    function(i) TraMineR:::wtd.boxplot.stats.tmr(dl[[i]], weights = wl[[i]])$stats,
    numeric(5)
  ))
  write_phase2_csv(
    paste0("seqrf_", tag, "_box"), bx,
    "# plot.seqrf right-panel weighted box stats per group: whislo,q1,med,q3,whishi\n"
  )
}
genrf_box("uw", FALSE)

## seqrf with string sortv (from.start / from.end): deterministic (no MDS), so the
## sortv resolution `order(as.data.frame(seqdata)[,end:beg]); order(...)` is goldened
## directly (RP3/RP4 audit F5 -- previously implemented but unpinned).
genrf_sortv <- function(sv) {
  srf <- suppressMessages(seqrf(sdrf, diss = Drf, k = 3, sortv = sv, weighted = TRUE))
  tag <- gsub("\\.", "", sv)  # "from.start" -> "fromstart"
  write_phase2_csv(
    paste0("seqrf_", tag, "_medoids"), matrix(as.vector(srf$rf$medoids) - 1L, nrow = 1),
    sprintf("# seqrf(sortv=%s, weighted=TRUE) group medoids (0-based)\n", sv)
  )
  write_phase2_csv(
    paste0("seqrf_", tag, "_heights"), matrix(srf$rf$heights, nrow = 1),
    sprintf("# seqrf(sortv=%s) heights\n", sv)
  )
  write_phase2_csv(
    paste0("seqrf_", tag, "_repcodes"), TraMineR:::seqasnum(srf$seqtoplot),
    sprintf("# seqrf(sortv=%s) representative codes (0-based)\n", sv)
  )
}
genrf_sortv("from.start")
genrf_sortv("from.end")

## Weighted-default-MDS counterexample.  The ordinary phase-2 fixture happens
## to have the same frequency groups under cmdscale() and wcmdscale(), so keep
## a small Euclidean case where TraMineR's weighted branch is observably
## load-bearing.  Values are fixed (rather than sampled) to make the oracle
## independent of R's random-number implementation.
mds_points <- matrix(c(
  -0.626453810742332418,  1.511781168450847979,  0.918977371608218241,
   0.183643324222082238,  0.389843236411431093,  0.782136300731067102,
  -0.835628612410047156, -0.621240580541803755,  0.074564983365190601,
   1.595280802137791554, -2.214699887177499882, -1.989351695863372793,
   0.329507771815360506,  1.124930918143108194,  0.619825747894710233,
  -0.820468384118015259, -0.044933609015230851, -0.056128739529000785,
   0.487429052428485277, -0.016190263098946087, -0.155795506705329295,
   0.738324705129217329,  0.943836210685299215, -1.470752383899274429,
   0.575781351653492313,  0.821221195098088552, -0.478150055108620353,
  -0.305388387156356023,  0.593901321217508826,  0.417941560199702411
), ncol = 3L, byrow = TRUE)
mds_weights <- c(
  1.75583379810765039, 0.59076552929498682, 0.90205808905519602,
  0.66374553509308376, 1.26290989297486611, 0.52381470274230102,
  0.93870588320886195, 1.48009329276445523, 0.19689280564047473,
  1.68518047111812352
)
D_mds <- as.matrix(dist(mds_points))
write_phase2_csv(
  "dissrf_weighted_mds_input_D", D_mds,
  "# weighted-default-MDS counterexample distance matrix\n"
)
write_phase2_csv(
  "dissrf_weighted_mds_input_weights", matrix(mds_weights, nrow = 1L),
  "# weighted-default-MDS counterexample case weights\n"
)

cat("SEQRF_GOLDENS_DONE\n")

## R2: complete grp.meth / squared / pow contract -----------------------------
## One long manifest freezes every result component that participates in the
## grouping contract.  Random cases have an explicit R seed so regenerating the
## oracle is byte deterministic; Python tests intentionally compare their
## invariants rather than R's RNG stream.
r2_rows <- list()
r2_add <- function(case, fun, section, row, col, value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
    case = case,
    function_name = fun,
    section = section,
    row = row,
    col = col,
    value = as.numeric(value),
    text = as.character(text),
    stringsAsFactors = FALSE
  )
}

r2_record <- function(case, fun, out) {
  rf <- if (fun == "seqrf") out$rf else out
  ## TraMineR returns the k=1 distance/index/weight vectors bare, whereas the
  ## Python result consistently re-expresses them as one group list.
  dist.list <- if (is.list(rf$dist.list)) rf$dist.list else list(rf$dist.list)
  index.list <- if (is.list(rf$index.list)) rf$index.list else list(rf$index.list)
  weights.list <- if (is.list(rf$weights.list)) rf$weights.list else list(rf$weights.list)
  for (i in seq_along(rf$medoids))
    r2_add(case, fun, "medoid", 1L, i, rf$medoids[[i]] - 1L)
  for (i in seq_along(rf$heights))
    r2_add(case, fun, "height", 1L, i, rf$heights[[i]])
  for (i in seq_along(rf$at))
    r2_add(case, fun, "at", 1L, i, rf$at[[i]])
  for (i in seq_along(index.list)) {
    for (j in seq_along(index.list[[i]])) {
      r2_add(case, fun, "index", i, j, index.list[[i]][[j]] - 1L)
      r2_add(case, fun, "weight", i, j, weights.list[[i]][[j]])
      r2_add(case, fun, "distance", i, j, dist.list[[i]][[j]])
    }
  }
  stats <- c(R2 = rf$R2, Fstat = rf$Fstat, pvalue = rf$pvalue)
  for (i in seq_along(stats))
    r2_add(case, fun, "statistic", 1L, i, stats[[i]], names(stats)[[i]])
  for (i in seq_along(rf$sizes))
    r2_add(case, fun, "size", 1L, i, rf$sizes[[i]], names(rf$sizes)[[i]])
  r2_add(case, fun, "attribute", 1L, 1L, text = paste0("grp.meth=", rf$grp.meth))

  if (!is.null(rf$mdsk)) {
    for (i in seq_along(rf$mdsk)) {
      r2_add(case, fun, "membership", 1L, i, rf$mdsk[[i]] - 1L)
      r2_add(case, fun, "case_medoid", 1L, i, rf$kmedoid.index[[i]] - 1L)
      r2_add(case, fun, "case_distance", 1L, i, rf$kmedoid.dist[[i]])
    }
  }
  if (fun == "seqrf") {
    codes <- TraMineR:::seqasnum(out$seqtoplot)
    for (i in seq_len(nrow(codes)))
      for (j in seq_len(ncol(codes)))
        r2_add(case, fun, "representative", i, j, codes[i, j])
  }
}

r2_diss_case <- function(case, grp.meth, sortv, k = 4L, squared = FALSE,
                         pow = NULL, weights = NULL, seed = NULL, diss = D) {
  if (!is.null(seed)) set.seed(seed)
  out <- suppressMessages(dissrf(
    diss, k = k, sortv = sortv, weights = weights, grp.meth = grp.meth,
    squared = squared, pow = pow
  ))
  r2_record(case, "dissrf", out)
}

r2_seq_case <- function(case, grp.meth, sortv, k = 4L, squared = FALSE,
                        pow = NULL, weights = NULL, weighted = FALSE,
                        seed = NULL, seqdata = sdrf, diss = Drf) {
  if (!is.null(seed)) set.seed(seed)
  out <- suppressMessages(seqrf(
    seqdata, diss = diss, k = k, sortv = sortv, weights = weights,
    weighted = weighted, grp.meth = grp.meth, squared = squared, pow = pow
  ))
  r2_record(case, "seqrf", out)
}

d_ties <- rep(c(0, 0, 1, 1, 2), length.out = nrow(D))
d_integer_weights <- rep(c(1, 2, 3), length.out = nrow(D))
r2_diss_case("d_prop_unweighted_numeric", "prop", sortv, k = 3L)
r2_diss_case("d_prop_integer_numeric", "prop", sortv, k = 3L, weights = d_integer_weights)
r2_diss_case("d_prop_fractional_numeric", "prop", sortv, k = 3L, weights = w)
r2_diss_case(
  "d_prop_fractional_default_mds", "prop", "mds", k = 3L,
  weights = mds_weights, diss = D_mds
)
r2_diss_case("d_first_numeric_sq0_default", "first", sortv, squared = FALSE)
r2_diss_case("d_first_ties_sq1_default", "first", d_ties, squared = TRUE)
r2_diss_case("d_first_null_sq0_pow2", "first", NULL, squared = FALSE, pow = 2L)
r2_diss_case("d_first_mds_sq1_pow1", "first", "mds", squared = TRUE, pow = 1L)
r2_diss_case("d_random_numeric_sq0_default", "random", sortv, squared = FALSE, seed = 6101L)
r2_diss_case("d_random_ties_sq1_default", "random", d_ties, squared = TRUE, seed = 6102L)
r2_diss_case("d_random_null_sq0_pow2", "random", NULL, squared = FALSE, pow = 2L, seed = 6103L)
r2_diss_case("d_random_mds_sq1_pow1", "random", "mds", squared = TRUE, pow = 1L, seed = 6104L)
r2_diss_case("d_k1_forces_prop", "first", sortv, k = 1L)

s_ties <- rep(c(0, 0, 1, 1, 2), length.out = nrow(Drf))
s_integer_weights <- rep(c(3, 2, 1), length.out = nrow(Drf))
r2_seq_case("s_prop_unweighted_numeric", "prop", sortvrf, k = 3L)
r2_seq_case("s_prop_integer_numeric", "prop", sortvrf, k = 3L, weights = s_integer_weights)
r2_seq_case("s_prop_fractional_numeric", "prop", sortvrf, k = 3L, weighted = TRUE)
r2_seq_case(
  "s_prop_fractional_default_mds", "prop", "mds", k = 3L,
  weights = mds_weights, weighted = TRUE,
  seqdata = sdrf[seq_len(nrow(D_mds)), ], diss = D_mds
)
r2_seq_case("s_first_numeric_sq0_default", "first", sortvrf, squared = FALSE)
r2_seq_case("s_first_ties_sq1_default", "first", s_ties, squared = TRUE)
r2_seq_case("s_first_null_sq0_pow2", "first", NULL, squared = FALSE, pow = 2L)
r2_seq_case("s_first_mds_sq1_pow1", "first", "mds", squared = TRUE, pow = 1L)
r2_seq_case("s_random_numeric_sq0_default", "random", sortvrf, squared = FALSE, seed = 7101L)
r2_seq_case("s_random_ties_sq1_default", "random", s_ties, squared = TRUE, seed = 7102L)
r2_seq_case("s_random_null_sq0_pow2", "random", NULL, squared = FALSE, pow = 2L, seed = 7103L)
r2_seq_case("s_random_mds_sq1_pow1", "random", "mds", squared = TRUE, pow = 1L, seed = 7104L)
r2_seq_case("s_k1_forces_prop", "random", sortvrf, k = 1L)

utils::write.csv(
  do.call(rbind, r2_rows),
  file.path("tests", "r_goldens", "diss_rf_r2_cases.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

r2_error <- function(expr) {
  tryCatch(
    {
      suppressMessages(force(expr))
      "VALUE"
    },
    error = function(e) paste0(
      "ERROR:", sub("[[:space:]]+$", "", conditionMessage(e))
    )
  )
}
writeLines(
  c(
    paste0("dissrf_weighted_first=", r2_error(dissrf(
      D, k = 3L, sortv = sortv, weights = w, grp.meth = "first"
    ))),
    paste0("dissrf_weighted_random=", r2_error(dissrf(
      D, k = 3L, sortv = sortv, weights = w, grp.meth = "random"
    ))),
    paste0("seqrf_weighted_first=", r2_error(seqrf(
      sdrf, Drf, k = 3L, sortv = sortvrf, weighted = TRUE, grp.meth = "first"
    ))),
    paste0("seqrf_explicit_weighted_random=", r2_error(seqrf(
      sdrf, Drf, k = 3L, sortv = sortvrf, weights = wrf,
      weighted = TRUE, grp.meth = "random"
    ))),
    paste0("dissrf_invalid_method=", r2_error(dissrf(
      D, k = 3L, sortv = sortv, grp.meth = "bogus"
    )))
  ),
  file.path("tests", "r_goldens", "diss_rf_r2_errors.txt")
)

cat("DISS_RF_R2_GOLDENS_DONE\n")
