# Golden generation: dissmfacw (weighted multi-factor analysis of a distance matrix).
#
# Goldens the deterministic decomposition at %.15g / atol=1e-10.  TraMineR's
# mfac table exposes PseudoF and PseudoR2; the explained SS column is derived as
# PseudoR2 * SCtot using the same public gower_matrix call and weighted trace as
# dissmfacw.R.  Permutation p-values are stochastic and are checked separately
# in Python.  R=1 makes these deterministic fixtures inexpensive and explicit.
#
# Requires live R + TraMineR 2.2-14; see DUPLICATES.md.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_mfac_golden <- function(name, result, diss, weights, provenance) {
  if (is.null(weights)) {
    weights <- rep(1, attr(as.dist(diss), "Size"))
  }
  gmat <- gower_matrix(as.matrix(diss), squared = FALSE, weights = weights)
  sctot <- sum(weights * diag(gmat))
  tab <- cbind(
    PseudoF = result$mfac$PseudoF,
    PseudoR2 = result$mfac$PseudoR2,
    SS = result$mfac$PseudoR2 * sctot
  )

  con <- file(file.path(out_dir, paste0(name, ".csv")), "w")
  cat("# R 4.6.2; TraMineR 2.2-14\n", file = con)
  cat(paste0("# ", provenance, "\n"), file = con)
  cat("# columns: PseudoF, PseudoR2, explained SS (SS = PseudoR2 * SCtot)\n", file = con)
  for (row in seq_len(nrow(tab))) {
    cat(paste(sprintf("%.15g", tab[row, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

inp <- phase2_inputs()
D <- inp$D
n <- inp$n
i <- seq_len(n)
g <- factor(inp$grp) # 3 levels: "1", "2", "3"
h <- factor(ifelse(i %% 2 == 0, "even", "odd"))
x <- seq(-1, 1, length.out = n)

# Existing categorical-only regression guard (weighted).
categorical_data <- data.frame(g = g, h = h)
categorical <- dissmfacw(
  D ~ g + h, data = categorical_data, weights = inp$wreal, R = 1, squared = FALSE
)
write_mfac_golden(
  "dissmfacw_det", categorical, D, inp$wreal,
  "TraMineR::dissmfacw(D ~ g + h, categorical_data, weights=wreal, R=1, squared=FALSE); rows g,h,Total"
)

# Numeric-only design.  The unused column works around TraMineR 2.2-14 dropping
# a one-column data.frame during its complete-case subsetting; it is not in the formula.
numeric_data <- data.frame(x = x, unused = i)
numeric <- dissmfacw(D ~ x, data = numeric_data, R = 1, squared = FALSE)
write_mfac_golden(
  "dissmfacw_numeric_unweighted", numeric, D, NULL,
  "TraMineR::dissmfacw(D ~ x, numeric_data, R=1, squared=FALSE); rows x,Total"
)

# Mixed predictors in formula order, both unweighted and weighted.
mixed_data <- data.frame(g = g, x = x)
mixed_uw <- dissmfacw(D ~ g + x, data = mixed_data, R = 1, squared = FALSE)
write_mfac_golden(
  "dissmfacw_mixed_unweighted", mixed_uw, D, NULL,
  "TraMineR::dissmfacw(D ~ g + x, mixed_data, R=1, squared=FALSE); rows g,x,Total"
)
mixed_w <- dissmfacw(D ~ g + x, data = mixed_data, weights = inp$wreal, R = 1, squared = FALSE)
write_mfac_golden(
  "dissmfacw_mixed_weighted", mixed_w, D, inp$wreal,
  "TraMineR::dissmfacw(D ~ g + x, mixed_data, weights=wreal, R=1, squared=FALSE); rows g,x,Total"
)

# Exact mixed design consumed by pyweightedcluster clustassoc: two clusters and
# a continuous covariate over a six-case block distance matrix.
toy_D <- matrix(1, 6, 6)
diag(toy_D) <- 0
toy_D[1:3, 1:3] <- 0.1
toy_D[4:6, 4:6] <- 0.1
diag(toy_D) <- 0
clustf <- factor(c(1, 1, 1, 2, 2, 2))
covar <- c(0, 0.1, 0.2, 2, 2.1, 2.2)
toy_data <- data.frame(clustf = clustf, covar = covar)
toy <- dissmfacw(toy_D ~ clustf + covar, data = toy_data, R = 1, squared = FALSE)
write_mfac_golden(
  "dissmfacw_clustassoc_mixed", toy, toy_D, NULL,
  "TraMineR::dissmfacw(toy_D ~ clustf + covar, toy_data, R=1, squared=FALSE); rows clustf,covar,Total"
)

# Predictor codes so Python rebuilds the categorical regression fixture.
write_phase2_csv(
  "phase2_input_grp", matrix(inp$grp, nrow = 1),
  "# Phase 2 shared input: latent group labels rep(1:3, each=5)\n"
)

cat("SEQDISSMFACW_GOLDENS_DONE\n")
print(categorical$mfac)
print(numeric$mfac)
print(mixed_uw$mfac)
print(mixed_w$mfac)
print(toy$mfac)
