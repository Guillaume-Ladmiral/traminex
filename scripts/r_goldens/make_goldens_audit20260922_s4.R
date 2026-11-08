# Golden generation: 2026-09-22 whole-codebase audit, package s4 (discrepancy
# analysis, dissimilarity trees, representatives and clustering helpers).
#
# Each section below reproduces one audit finding against TraMineR 2.2-14 and
# writes one fixture under tests/r_goldens/; the consuming test is named in the
# section header.  Values are written with sprintf("%.17g"); R's NA is written
# as "NA".
#
# Permutation tests draw their permutations with sample(n) inside
# TraMineR.permutation (TraMineR-permutation.R:16).  R's and NumPy's random
# streams differ, so a section that needs permuted statistics replaces
# base::sample, for the duration of one call, with a function that replays a
# fixed permutation script.  The script is deterministic arithmetic (no RNG) and
# is written into the fixture, so the Python test replays exactly the same
# permutations.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

num <- function(x) paste(sprintf("%.17g", as.vector(x)), collapse = " ")
txt <- function(x) paste(ifelse(is.na(x), "NA", as.character(x)), collapse = " ")
write_fixture <- function(name, header, lines) {
  writeLines(c(paste("#", header), lines), file.path(out_dir, name))
}

scripted <- new.env()
scripted_sample <- function(x, size, replace = FALSE, prob = NULL) {
  scripted$k <- scripted$k + 1L
  p <- as.integer(scripted$perms[scripted$k, ])
  stopifnot(length(x) == 1L, length(p) == x, missing(size), !replace, is.null(prob))
  p
}
with_scripted_sample <- function(perms, expr) {
  real_sample <- base::sample
  unlockBinding("sample", baseenv())
  assign("sample", scripted_sample, envir = baseenv())
  on.exit({
    assign("sample", real_sample, envir = baseenv())
    lockBinding("sample", baseenv())
  })
  scripted$perms <- perms
  scripted$k <- 0L
  result <- force(expr)
  stopifnot(scripted$k == nrow(perms))
  result
}
# Affine permutations of 1..n for prime n: case i goes to ((a*(i-1) + r) mod n) + 1.
affine_perms <- function(n, count) {
  t(vapply(seq_len(count), function(r) {
    as.integer(((((r %% 11) + 2) * (0:(n - 1)) + r) %% n) + 1)
  }, integer(n)))
}

## --------------------------------------------------------------------------
## S4-01: under weight.permutation = "diss" the Levene residuals are taken at
## the permuted positions (dissassocweighted.R:119-125: wind <- sort(ind[...]);
## WeightedResidualsZ(dissc[wind], ww[wind])).
## Consumed by tests/test_dissassoc_permutation_mechanics.py.
## --------------------------------------------------------------------------
x <- c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3, 4.0, 4.6, 5.2, 6.0, 1.1, 2.9)
y <- c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2, 1.4, 0.1, 2.6, 0.9, 3.3, 0.3)
D <- as.matrix(dist(cbind(x, y)))
g <- c("a", "a", "b", "a", "b", "c", "b", "c", "a", NA, "b", "c", "a", "c")
w <- c(1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 0.4, 1.8, 0.6, 1.2)
R <- 30
perms <- affine_perms(sum(!is.na(g)), R - 1)
res <- with_scripted_sample(perms, dissassoc(D, g, weights = w, R = R, weight.permutation = "diss"))
lines <- c(
  paste("D", num(t(D)), sep = "\t"),
  paste("group", txt(g), sep = "\t"),
  paste("weights", num(w), sep = "\t"),
  paste("R", R, sep = "\t"),
  paste("perm", apply(perms, 1, paste, collapse = " "), sep = "\t"),
  paste("t", apply(res$perms$t, 1, num), sep = "\t"),
  paste("pval", num(res$perms$pval), sep = "\t")
)
write_fixture(
  "dissassoc_diss_scripted_permutations.txt",
  "dissassoc with weights, R = 30 and weight.permutation = diss, replaying the perm rows",
  lines
)

## --------------------------------------------------------------------------
## S4-02: disstree keeps a factor's level order, splits an ordered factor over
## its levels (not over numeric thresholds) and treats a logical predictor as
## an unordered factor (disstree.R:404 factor(pred, ordered = is.ordered(pred)
## || is.numeric(pred)); breaks only when is.numeric(pred), disstree.R:543).
## Consumed by tests/test_disstree_predictor_types.py.
## --------------------------------------------------------------------------
tree_lines <- function(case, tr) {
  out <- character()
  walk <- function(node) {
    if (!is.null(node$split)) {
      s <- node$split
      out <<- c(out, paste(case, "split", node$id,
        paste0("varindex=", s$varindex),
        paste0("breaks=", if (is.null(s$breaks)) "NA" else sprintf("%.17g", s$breaks)),
        paste0("labels=", if (is.null(s$labels)) "NA" else paste(s$labels, collapse = "|")),
        paste0("index=", paste(s$index, collapse = ",")),
        paste0("nagroup=", if (is.null(s$naGroup)) "NA" else s$naGroup),
        paste0("lpop=", sprintf("%.17g", s$info$lpop)),
        paste0("rpop=", sprintf("%.17g", s$info$rpop)),
        paste0("R2=", sprintf("%.17g", s$info$R2)),
        sep = "\t"))
      walk(node$kids[[1]])
      walk(node$kids[[2]])
    }
  }
  walk(tr$root)
  c(out,
    paste(case, "fitted", paste(tr$fitted[, 1], collapse = " "), sep = "\t"),
    paste(case, "rules", paste(disstree.get.rules(tr), collapse = " || "), sep = "\t"))
}
pos <- c(0.0, 0.1, 0.2, 0.5, 0.6, 0.7, 5.0, 5.1, 5.2, 5.5, 5.6, 5.7)
Dpos <- as.matrix(dist(pos))
predictor_cases <- list(
  ordered_string = factor(c("low", "low", "low", "high", "high", "high",
    "mid", "mid", "mid", "mid", "mid", "mid"), levels = c("low", "mid", "high"), ordered = TRUE),
  ordered_integer = factor(c(3, 3, 3, 2, 2, 2, 1, 1, 1, 1, 1, 1), levels = c(3, 1, 2), ordered = TRUE),
  logical = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  factor_level_order = factor(c("a", "a", "a", "b", "b", "b", "c", "c", "c", "d", "d", "d"),
    levels = c("d", "c", "b", "a")),
  ordered_missing_unused = factor(c("x", "x", "x", NA, "y", "y", "z", "z", "z", "z", "z", "z"),
    levels = c("w", "z", "y", "x"), ordered = TRUE)
)
lines <- character()
for (case in names(predictor_cases)) {
  v <- predictor_cases[[case]]
  lines <- c(lines, paste(case, "values", txt(as.character(v)), sep = "\t"))
  tr <- disstree(Dpos ~ v, data = data.frame(v = v), R = 1, min.size = 0, max.depth = 3)
  assigned <- tryCatch(
    paste(disstree.assign(disstree.get.rules(tr), data.frame(v = v)), collapse = " "),
    error = function(e) paste("error", gsub("\\s+", " ", trimws(conditionMessage(e))))
  )
  lines <- c(lines, tree_lines(case, tr), paste(case, "assign", assigned, sep = "\t"))
}
write_fixture(
  "disstree_predictor_types.txt",
  "disstree of D on v with R = 1, min.size = 0 and max.depth = 3 for five predictor types; D is the distance between the pos values",
  c(paste("pos", num(pos), sep = "\t"), lines)
)

## --------------------------------------------------------------------------
## S4-03 and S4-07: a permutation p-value is NA when its column holds an NA,
## observed or permuted (TraMineR-permutation.R:21 sum(t[, i] >= t0[i]) / R),
## and every p-value is NA when R <= 1 (TraMineR-permutation.R:25-26);
## dissmfacw uses the same TraMineR.permutation (dissmfacw.R:124).
## Consumed by tests/test_permutation_na_pvalues.py.
## --------------------------------------------------------------------------
lines <- character()
# (a) A singleton group makes the unweighted Bartlett statistic NaN (1/(ni-1)).
xs <- c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3, 4.0, 4.6, 5.2, 6.0, 2.7)
ys <- c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2, 1.4, 0.1, 2.6, 0.9, 1.9)
Ds <- as.matrix(dist(cbind(xs, ys)))
gs <- c("a", "a", "b", "a", "b", "c", "b", "c", "a", "c", "b", "c", "s")
perms_s <- affine_perms(length(gs), 19)
res <- with_scripted_sample(perms_s, dissassoc(Ds, gs, R = 20))
lines <- c(lines,
  paste("singleton", "D", num(t(Ds)), sep = "\t"),
  paste("singleton", "group", txt(gs), sep = "\t"),
  paste("singleton", "perm", apply(perms_s, 1, paste, collapse = " "), sep = "\t"),
  paste("singleton", "t", apply(res$perms$t, 1, num), sep = "\t"),
  paste("singleton", "pval", num(res$perms$pval), sep = "\t"))
# (b) Duplicated points: the observed grouping is finite, but the second
# scripted permutation groups each duplicate pair, so every within-group
# variance is 0 and that permuted Bartlett statistic is NaN.
Dd <- as.matrix(dist(c(0, 0, 5, 5, 9, 9)))
gd <- c("a", "b", "c", "a", "b", "c")
perms_d <- rbind(c(2, 1, 3, 4, 6, 5), c(1, 3, 5, 2, 4, 6), c(6, 5, 4, 3, 2, 1), c(3, 1, 2, 6, 4, 5))
res <- with_scripted_sample(perms_d, dissassoc(Dd, gd, R = 5))
lines <- c(lines,
  paste("permuted_nan", "D", num(t(Dd)), sep = "\t"),
  paste("permuted_nan", "group", txt(gd), sep = "\t"),
  paste("permuted_nan", "perm", apply(perms_d, 1, paste, collapse = " "), sep = "\t"),
  paste("permuted_nan", "t", apply(res$perms$t, 1, num), sep = "\t"),
  paste("permuted_nan", "pval", num(res$perms$pval), sep = "\t"))
# (c) weight.permutation = "random-sampling" with unit weights and the singleton:
# the observed Bartlett statistic is NaN, so its p-value is NA whatever the
# draws; the other p-values depend on R's draws and are not recorded.
set.seed(20260922)
res <- dissassoc(Ds, gs, weights = rep(1, length(gs)), R = 20, weight.permutation = "random-sampling")
lines <- c(lines,
  paste("random_sampling", "t0", num(res$perms$t0), sep = "\t"),
  paste("random_sampling", "bartlett_pval", num(res$perms$pval[4]), sep = "\t"))
# (d) dissmfacw with R = 1 has no permutation, so every p-value is NA.
mf <- dissmfacw(Ds ~ g + h, data = data.frame(g = factor(gs), h = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1)),
  R = 1, squared = FALSE)
lines <- c(lines,
  paste("dissmfacw_R1", "PseudoF", num(mf$mfac$PseudoF), sep = "\t"),
  paste("dissmfacw_R1", "PseudoR2", num(mf$mfac$PseudoR2), sep = "\t"),
  paste("dissmfacw_R1", "p_value", num(mf$mfac$p_value), sep = "\t"))
# (e) An all-zero dissimilarity matrix makes every pseudo F 0/0 = NaN, so the
# p-values are NA even with permutations, whatever the draws.
set.seed(20260922)
mf <- dissmfacw(matrix(0, 13, 13) ~ g + h, data = data.frame(g = factor(gs), h = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1)),
  R = 5, squared = FALSE)
lines <- c(lines,
  paste("dissmfacw_zero", "PseudoF", num(mf$mfac$PseudoF), sep = "\t"),
  paste("dissmfacw_zero", "p_value", num(mf$mfac$p_value), sep = "\t"))
write_fixture(
  "permutation_na_pvalues.txt",
  "NA permutation p-values: dissassoc with scripted permutations or random sampling; dissmfacw with R = 1",
  lines
)

## --------------------------------------------------------------------------
## S4-04: when a level subset and its complement have mathematically equal
## residual inertia, rounding decides which becomes the left child. R's rounding
## comes from DTNGroupFactorBinary (disstree.R:403-523): sum() accumulates in
## long double in column-major order, and the C kernels tmrWeightedInertiaDist /
## tmrWeightedInterInertia (weightedinertia.cpp) accumulate doubles in member
## order. Every split's SCres, sizes and orientation, and the fitted leaves.
## Consumed by tests/test_disstree_split_arithmetic.py.
## --------------------------------------------------------------------------
split_lines <- function(case, tr) {
  out <- character()
  walk <- function(node) {
    if (!is.null(node$split)) {
      s <- node$split
      out <<- c(out, paste(case, "split", node$id,
        paste0("labels=", paste(s$labels, collapse = "|")),
        paste0("index=", paste(s$index, collapse = ",")),
        paste0("lpop=", sprintf("%.17g", s$info$lpop)),
        paste0("rpop=", sprintf("%.17g", s$info$rpop)),
        paste0("lvar=", sprintf("%.17g", s$info$lvar)),
        paste0("rvar=", sprintf("%.17g", s$info$rvar)),
        paste0("SCres=", sprintf("%.17g", s$info$SCres)),
        sep = "\t"))
      walk(node$kids[[1]])
      walk(node$kids[[2]])
    }
  }
  walk(tr$root)
  c(out, paste(case, "fitted", paste(tr$fitted[, 1], collapse = " "), sep = "\t"))
}
tie_cases <- list()
tie_cases[["audit_p17"]] <- list(
  x = c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3, 4.0, 4.6, 5.2, 6.0, 1.1, 2.9, 0.7, 3.8),
  y = c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2, 1.4, 0.1, 2.6, 0.9, 3.3, 0.3, 2.4, 1.9),
  f = c("a", "a", "b", "a", "b", "c", "b", "c", "a", "c", "b", "c", "a", "c", "b", "a"),
  w = c(1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 0.4, 1.8, 0.6, 1.2, 0.8, 1.4)
)
for (k in 1:6) {
  set.seed(20260922 + k)
  n_k <- 20
  tie_cases[[paste0("seeded_", k)]] <- list(
    x = round(runif(n_k, 0, 6), 1),
    y = round(runif(n_k, 0, 4), 1),
    f = sample(letters[1:(3 + k %% 2)], n_k, replace = TRUE),
    w = if (k %% 3 == 0) NULL else round(runif(n_k, 0.3, 2), 2)
  )
}
# Many levels, so that a candidate side holds several levels: seven unordered
# levels, and six ordered levels (an ordered split takes a prefix of them).
for (k in 7:8) {
  set.seed(20260922 + k)
  n_k <- 40
  tie_cases[[paste0("seeded_", k)]] <- list(
    x = round(runif(n_k, 0, 6), 2),
    y = round(runif(n_k, 0, 4), 2),
    f = sample(letters[1:(14 - k)], n_k, replace = TRUE),
    w = round(runif(n_k, 0.3, 2), 3),
    ordered = k == 8
  )
}
lines <- character()
for (case in names(tie_cases)) {
  dat <- tie_cases[[case]]
  Dk <- as.matrix(dist(cbind(dat$x, dat$y)))
  f <- if (isTRUE(dat$ordered)) factor(dat$f, ordered = TRUE) else dat$f
  tr <- disstree(Dk ~ f, data = data.frame(f = f), weights = dat$w, R = 1,
    min.size = 0.1, max.depth = 3, weight.permutation = "diss")
  lines <- c(lines,
    paste(case, "D", num(t(Dk)), sep = "\t"),
    paste(case, "f", txt(dat$f), sep = "\t"),
    paste(case, "ordered", if (isTRUE(dat$ordered)) "TRUE" else "FALSE", sep = "\t"),
    paste(case, "weights", if (is.null(dat$w)) "none" else num(dat$w), sep = "\t"),
    split_lines(case, tr))
}
write_fixture(
  "disstree_split_arithmetic.txt",
  "disstree of D on one factor f with R = 1, min.size = 0.1, max.depth = 3 and diss weight permutation",
  lines
)

## --------------------------------------------------------------------------
## S4-05: a numeric threshold is printed with format(x, digits =
## getOption("digits") - 2) (disstree.R:63-70), which keeps fixed notation
## unless it is wider than scientific notation; disstree.assign then evaluates
## the printed rule.
## Consumed by tests/test_disstree_threshold_rules.py.
## --------------------------------------------------------------------------
stopifnot(getOption("digits") == 7)
lines <- paste("pos", "values", num(pos), sep = "\t")
for (base in c(123456, 1234567.5, 99999.5, 0.000123456, 12.345678, -2.5, 1e-10, 31415926.5)) {
  inc <- c(rep(base, 6), rep(base + abs(base), 6))
  tr <- disstree(Dpos ~ inc, data = data.frame(inc = inc), R = 1, min.size = 0, max.depth = 2)
  rules <- disstree.get.rules(tr)
  profile <- data.frame(inc = c(base, base + abs(base) / 1000, base + abs(base)))
  lines <- c(lines,
    paste(sprintf("%.17g", base), "breaks", sprintf("%.17g", tr$root$split$breaks), sep = "\t"),
    paste(sprintf("%.17g", base), "rules", paste(rules, collapse = " || "), sep = "\t"),
    paste(sprintf("%.17g", base), "profile", num(profile$inc), sep = "\t"),
    paste(sprintf("%.17g", base), "assign", paste(disstree.assign(rules, profile), collapse = " "), sep = "\t"))
}
write_fixture(
  "disstree_threshold_rules.txt",
  "disstree of the pos distances on one numeric predictor (base x6, base + |base| x6); rules and disstree.assign",
  lines
)

## --------------------------------------------------------------------------
## S4-06: dissassoc with weights and weight.permutation = "none" matches no
## permutation branch (dissassocweighted.R:342-359), so ret$perms is NULL and R
## stops; the weighted "diss" run on the same data is recorded for contrast.
## Consumed by tests/test_dissassoc_permutation_mechanics.py.
## --------------------------------------------------------------------------
refusal <- tryCatch({
  dissassoc(D, g, weights = w, R = 1, weight.permutation = "none")
  "no error"
}, error = function(e) paste("error", gsub("\\s+", " ", trimws(conditionMessage(e)))))
res <- dissassoc(D, g, weights = w, R = 1, weight.permutation = "diss")
write_fixture(
  "dissassoc_weights_none.txt",
  "dissassoc on the S4-01 data with weights and weight.permutation = none (R stops), and with diss",
  c(paste("none", refusal, sep = "\t"), paste("diss_t0", num(res$stat[, 1]), sep = "\t"))
)

## --------------------------------------------------------------------------
## S4-08: dissrep on an all-zero matrix sets nb to the total weight
## (dissrep.R:16, 178: diss0 rule) and divides plainly, so Q and Quality are
## NaN when DC = 0 (dissrep.R:188, 200).
## Consumed by tests/test_dissrep_zero_matrix.py.
## --------------------------------------------------------------------------
zero_cases <- list(
  unweighted = list(n = 5, weights = NULL, nrep = NULL, dmax = NULL),
  weighted_nrep = list(n = 4, weights = c(1, 2, 0.5, 1.5), nrep = 2, dmax = NULL),
  explicit_dmax = list(n = 4, weights = NULL, nrep = NULL, dmax = 1)
)
lines <- character()
for (case in names(zero_cases)) {
  z <- zero_cases[[case]]
  r <- suppressMessages(dissrep(matrix(0, z$n, z$n), criterion = "density", weights = z$weights,
    nrep = z$nrep, dmax = z$dmax))
  st <- attr(r, "Statistics")
  lines <- c(lines,
    paste(case, "n", z$n, sep = "\t"),
    paste(case, "weights", if (is.null(z$weights)) "none" else num(z$weights), sep = "\t"),
    paste(case, "nrep", if (is.null(z$nrep)) "none" else z$nrep, sep = "\t"),
    paste(case, "dmax", if (is.null(z$dmax)) "none" else num(z$dmax), sep = "\t"),
    paste(case, "representatives", paste(as.vector(r), collapse = " "), sep = "\t"),
    paste(case, "statistics", num(as.matrix(st)), sep = "\t"),
    paste(case, "quality", num(attr(r, "Quality")), sep = "\t"))
}
write_fixture(
  "dissrep_zero_matrix.txt",
  "dissrep of all-zero matrices; statistics column-major (na, na(%), nb, nb(%), SD, MD, DC, V, Q)",
  lines
)

## --------------------------------------------------------------------------
## S4-10: seqrep.grp drops cases whose group is missing: TraMineR's group()
## keeps NA, and xtabs(~grp), levels(factor(grp)) and which(grp == level)
## all leave NA out (TraMineRextras seqrep-grp.R:15-45).
## Consumed by tests/test_seqrep_grp_missing_group.py.
## --------------------------------------------------------------------------
suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
rep_rows <- rbind(c("A", "A", "B", "B"), c("A", "B", "B", "B"), c("A", "A", "A", "B"),
  c("B", "B", "C", "C"), c("B", "C", "C", "C"), c("A", "A", "B", "C"),
  c("C", "C", "C", "C"), c("A", "B", "C", "C"))
rep_seq <- suppressMessages(seqdef(rep_rows, alphabet = c("A", "B", "C")))
rep_diss <- suppressMessages(seqdist(rep_seq, method = "LCS"))
lines <- c(
  paste("input", "rows", apply(rep_rows, 1, paste, collapse = ""), sep = "\t"),
  paste("input", "diss", num(t(rep_diss)), sep = "\t"))
group_cases <- list(
  character = c("x", "x", "x", "y", "y", NA, "y", "x"),
  numeric = c(1, 1, 1, 2, 2, NA, 2, 1)
)
for (case in names(group_cases)) {
  grp <- group_cases[[case]]
  st <- suppressMessages(seqrep.grp(rep_seq, group = grp, diss = rep_diss, ret = "stat",
    criterion = "density", coverage = 0.25, weighted = TRUE))
  lines <- c(lines, paste(case, "group", txt(grp), sep = "\t"))
  for (level in names(st)) {
    lines <- c(lines, paste(case, paste0("stat_", trimws(level)), num(as.matrix(st[[level]])), sep = "\t"))
  }
}
write_fixture(
  "seqrep_grp_missing_group.txt",
  "seqrep.grp with LCS distances and a missing group value; statistics column-major per level",
  lines
)

## --------------------------------------------------------------------------
## S4-11: dissassoc, disscenter and dissvar.grp take the group levels from
## factor(), which keeps a factor's declared level order and drops unused
## levels (dissassocweighted.R:294, disscenter.R:47, dissvar-grp.R:22).
## Consumed by tests/test_group_level_order.py.
## --------------------------------------------------------------------------
xo <- c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3)
yo <- c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2)
Do <- as.matrix(dist(cbind(xo, yo)))
wo <- c(1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9)
go <- factor(c("low", "low", "high", "mid", "high", "mid", "low", "high"),
  levels = c("low", "none", "mid", "high"))
ra <- dissassoc(Do, go, weights = wo, R = 1, weight.permutation = "diss")
ac <- disscenter(Do, group = go, allcenter = TRUE, weights = wo)
md <- disscenter(Do, group = go, medoids.index = "first", weights = wo)
dv <- dissvar.grp(Do, group = go)
write_fixture(
  "group_level_order.txt",
  "group levels low < none < mid < high, none unused; results of dissassoc with diss weights, disscenter and dissvar.grp",
  c(paste("D", num(t(Do)), sep = "\t"),
    paste("weights", num(wo), sep = "\t"),
    paste("group", txt(as.character(go)), sep = "\t"),
    paste("categories", paste(levels(go), collapse = " "), sep = "\t"),
    paste("dissassoc_rows", paste(rownames(ra$groups), collapse = " "), sep = "\t"),
    paste("dissassoc_groups", num(as.matrix(ra$groups)), sep = "\t"),
    paste("allcenter_cols", paste(names(ac), collapse = " "), sep = "\t"),
    paste("allcenter", num(as.matrix(ac)), sep = "\t"),
    paste("medoid_names", paste(names(md), collapse = " "), sep = "\t"),
    paste("medoids", num(unlist(md)), sep = "\t"),
    paste("dissvar_grp_names", paste(names(dv), collapse = " "), sep = "\t"),
    paste("dissvar_grp", num(dv), sep = "\t"))
)

## --------------------------------------------------------------------------
## S4-12: with a group whose weighted size is below 1, WeightedCluster's ASW is
## NA (with a warning); dissmergegroups then stops at its first comparison when
## the merge loop runs (dissmergegroups.R:124-165), and otherwise returns the
## grouping with quality NA.
## Consumed by tests/test_dissmergegroups_small_weights.py.
## --------------------------------------------------------------------------
suppressMessages(library(WeightedCluster))
stopifnot(as.character(packageVersion("WeightedCluster")) == "2.0")
pm <- c(0.0, 0.2, 0.5, 3.0, 3.1, 3.6, 7.0, 7.4, 7.5, 12.0, 12.3, 12.9)
Dm <- as.matrix(dist(pm))
gm <- rep(1:4, each = 3)
wm <- c(0.1, 0.2, 0.3, 0.2, 0.25, 0.3, 0.5, 0.4, 0.3, 0.2, 0.2, 0.2)
merge_case <- function(min_group) {
  seen <- character()
  result <- tryCatch(
    withCallingHandlers(
      dissmergegroups(Dm, gm, weights = wm, min.group = min_group, small = 0.05, silent = TRUE),
      warning = function(w) {
        seen <<- c(seen, gsub("\\s+", " ", trimws(conditionMessage(w))))
        invokeRestart("muffleWarning")
      }),
    error = function(e) e)
  warned <- if (length(seen)) paste(unique(seen), collapse = " | ") else "no warning"
  if (inherits(result, "error")) {
    paste("error", gsub("\\s+", " ", trimws(conditionMessage(result))), sep = "\t")
  } else {
    paste("value", paste(as.vector(result), collapse = " "), num(attr(result, "qual")), warned, sep = "\t")
  }
}
write_fixture(
  "dissmergegroups_small_weights.txt",
  "dissmergegroups with four groups of weighted size below 1; merge loop running or not",
  c(paste("positions", num(pm), sep = "\t"),
    paste("group", paste(gm, collapse = " "), sep = "\t"),
    paste("weights", num(wm), sep = "\t"),
    paste("loop_runs", merge_case(2), sep = "\t"),
    paste("no_loop", merge_case(4), sep = "\t"))
)

## --------------------------------------------------------------------------
## S4-09 (dissdomassoc part): each domain matrix and jointdiss is vectorised
## with as.dist (dissdomassoc.R:55, 63), i.e. its lower triangle column by
## column; an asymmetric matrix is read from that triangle only.
## Consumed by tests/test_dissdomassoc_asymmetric.py.
## --------------------------------------------------------------------------
xa <- c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3)
ya <- c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2)
Sa <- as.matrix(dist(cbind(xa, ya)))
Aa <- Sa
Aa[lower.tri(Aa)] <- Aa[lower.tri(Aa)] * (1 + 0.1 * seq_len(sum(lower.tri(Aa))) %% 7)
Aa[upper.tri(Aa)] <- Aa[upper.tri(Aa)] * 1.5
Ja <- Sa + t(Aa) / 2
wa <- c(1.5, 0.5, 2.0, 1.0, 0.7, 1.3, 2.2, 0.9)
da <- dissdomassoc(list(Aa, Sa), jointdiss = Ja, what = c("pearson", "spearman", "R2"), weights = wa)
write_fixture(
  "dissdomassoc_asymmetric.txt",
  "dissdomassoc of an asymmetric and a symmetric domain with an asymmetric joint matrix; matrices column-major",
  c(paste("A", num(t(Aa)), sep = "\t"),
    paste("S", num(t(Sa)), sep = "\t"),
    paste("J", num(t(Ja)), sep = "\t"),
    paste("weights", num(wa), sep = "\t"),
    paste("Pearson", num(da$Pearson), sep = "\t"),
    paste("Spearman", num(da$Spearman), sep = "\t"),
    paste("Pearson.Rsquare", num(da$Pearson.Rsquare), sep = "\t"),
    paste("Spearman.Rsquare", num(da$Spearman.Rsquare), sep = "\t"))
)

cat("AUDIT20260922_S4_GOLDENS_DONE\n")
