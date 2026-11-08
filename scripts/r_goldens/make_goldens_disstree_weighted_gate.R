# Golden generation: disstree / seqtree options that no test of this project ran
# (supervision repository, docs/TRAMINEX_COMPONENT_GAPS.md, 2026-11-02): the weighted
# branches of the split permutation test, and squared, first and the rounding of weights.
#
# Permutation tests draw with sample(n) inside TraMineR.permutation
# (TraMineR-permutation.R:16). R's and NumPy's random streams differ, so base::sample is
# replaced, for the duration of one call, by a fixed script of deterministic arithmetic
# (no RNG), as in make_goldens_audit20260922_s4.R. A tree asks for permutations of several
# sizes (a node's cases, or its replicated cases, then the final adjustment), so the script
# answers any size, and each answer is written into the fixture in the order R asked: the
# test replays those rows and so also pins the order and the size of R's requests.
#
# Values are written with sprintf("%.17g"); a refusal is written with R's message.
# Consumed by tests/test_r_goldens_disstree_weighted_gate.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

num <- function(x) paste(sprintf("%.17g", as.vector(x)), collapse = " ")

## A fixed permutation of 1..n for request number k: Fisher-Yates driven by a Lehmer step,
## exact in doubles because every product stays below 2^53.
script_permutation <- function(n, k) {
  p <- seq_len(n)
  s <- k * 7919 + n
  if (n > 1L) {
    for (i in n:2L) {
      s <- (s * 48271) %% 2147483647
      j <- (s %% i) + 1
      swap <- p[i]
      p[i] <- p[j]
      p[j] <- swap
    }
  }
  as.integer(p)
}

scripted <- new.env()
scripted_sample <- function(x, size, replace = FALSE, prob = NULL) {
  stopifnot(length(x) == 1L, missing(size), !replace, is.null(prob))
  scripted$k <- scripted$k + 1L
  p <- script_permutation(x, scripted$k)
  scripted$rows[[scripted$k]] <- paste(p, collapse = " ")
  p
}
with_scripted_sample <- function(expr) {
  real_sample <- base::sample
  unlockBinding("sample", baseenv())
  assign("sample", scripted_sample, envir = baseenv())
  on.exit({
    assign("sample", real_sample, envir = baseenv())
    lockBinding("sample", baseenv())
  })
  scripted$k <- 0L
  scripted$rows <- list()
  force(expr)
}

## Two-dimensional positions, so that distances have no ties. Three real weights end in .5,
## where R's round() goes to the even number. The fourth set also rounds two weights to zero,
## which can leave a candidate group with no weight: R then stops on a NaN comparison.
x <- c(0.0, 0.3, 0.9, 1.4, 2.2, 2.5, 3.1, 3.3, 4.0, 4.6, 5.2, 6.0, 1.1, 2.9)
y <- c(1.0, 0.2, 1.7, 0.4, 2.9, 1.1, 0.5, 2.2, 1.4, 0.1, 2.6, 0.9, 3.3, 0.3)
D <- as.matrix(dist(cbind(x, y)))
dimnames(D) <- NULL
real <- c(1.5, 1.4, 2.5, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 3.5, 1.8, 0.6, 1.2)
integer <- c(2, 1, 3, 1, 1, 2, 2, 1, 1, 2, 1, 2, 1, 1)
rounds_to_zero <- c(1.5, 0.5, 2.5, 1.0, 0.7, 1.3, 2.2, 0.9, 1.1, 1.6, 0.4, 1.8, 0.6, 1.2)
predictors <- data.frame(
  phase = factor(c("a", "a", "b", "a", "b", "c", "b", "c", "a", "c", "b", "c", "a", "c")),
  cohort = c(3.5, 1.2, 4.1, 2.2, 0.4, 5.3, 3.9, 1.7, 2.8, 4.6, 0.9, 5.8, 1.1, 3.2),
  noise = factor(c("u", "v", "u", "v", "v", "u", "u", "v", "u", "v", "u", "u", "v", "v"))
)
weight_sets <- list(real = real, integer = integer, rounds_to_zero = rounds_to_zero, none = NULL)

states <- t(vapply(1:14, function(r) c("A", "B", "C")[((r + 1:6) %% 3) + 1], character(6)))

flatten_tree <- function(tree) {
  rows <- character(0)
  rec <- function(node) {
    info <- node$info
    sp <- node$split
    tail <- if (is.null(sp)) rep(NA_real_, 9) else c(
      sp$varindex, sp$info$SCres, sp$info$lpop, sp$info$rpop, sp$info$lvar, sp$info$rvar, sp$info$R2,
      if (is.null(sp$info$pval)) NA_real_ else sp$info$pval,
      if (is.null(sp$breaks)) NA_real_ else sp$breaks
    )
    rows <<- c(rows, num(c(node$id, info$depth, info$splitschedule, info$n, info$vardis, info$medoid,
                           !is.null(sp), tail)))
    if (!is.null(node$kids)) {
      rec(node$kids[[1]])
      rec(node$kids[[2]])
    }
  }
  rec(tree$root)
  rows
}

case_lines <- function(name, tree, weights_name) {
  stat <- tree$info$adjustment$stat
  c(
    paste(name, "weights", weights_name, sep = "\t"),
    if (length(scripted$rows)) paste(name, "perm", unlist(scripted$rows), sep = "\t"),
    paste(name, "node", flatten_tree(tree), sep = "\t"),
    paste(name, "fitted", num(tree$fitted[, 1]), sep = "\t"),
    paste(name, "adjust_t0", num(stat[, "t0"]), sep = "\t"),
    paste(name, "adjust_p", num(stat[, "p.value"]), sep = "\t")
  )
}

## Every argument is named at the call, so that tests/test_golden_default_coverage.py, which reads
## these calls, sees what they pass. It cannot see through R's "...".
R <- 30L
tree_case <- function(name, weights_name, R, pval = 1, weight.permutation = "replicate", squared = FALSE,
                      first = NULL, min.size = 0.05, max.depth = 3, formula = D ~ phase + cohort) {
  w <- weight_sets[[weights_name]]
  tree <- suppressMessages(with_scripted_sample(
    disstree(formula, data = predictors, weights = w, min.size = min.size, max.depth = max.depth, R = R,
             pval = pval, weight.permutation = weight.permutation, squared = squared, first = first)
  ))
  case_lines(name, tree, weights_name)
}
seqtree_case <- function(name, weights_name, R, pval = 1, weight.permutation = "replicate", squared = FALSE,
                         first = NULL) {
  w <- weight_sets[[weights_name]]
  sequences <- suppressMessages(seqdef(states, alphabet = c("A", "B", "C"), weights = w))
  tree <- suppressMessages(with_scripted_sample(
    seqtree(sequences ~ phase + cohort, data = predictors, diss = D, min.size = 0.05, max.depth = 3, R = R,
            pval = pval, weight.permutation = weight.permutation, squared = squared, first = first)
  ))
  case_lines(name, tree, weights_name)
}
refusal <- function(name, weights = NULL, weight.permutation = "replicate", first = NULL) {
  message_text <- tryCatch({
    suppressMessages(disstree(D ~ phase + cohort, data = predictors, weights = weights, min.size = 0.05,
                              max.depth = 3, R = 1, weight.permutation = weight.permutation, first = first))
    stop("R did not refuse ", name)
  }, error = function(e) gsub("[\t\n]+", " ", conditionMessage(e)))
  paste(name, "refusal", message_text, sep = "\t")
}

lines <- c(
  paste("input", "D", num(t(D)), sep = "\t"),
  paste("input", "real", num(real), sep = "\t"),
  paste("input", "integer", num(integer), sep = "\t"),
  paste("input", "rounds_to_zero", num(rounds_to_zero), sep = "\t"),
  paste("input", "phase", paste(predictors$phase, collapse = " "), sep = "\t"),
  paste("input", "cohort", num(predictors$cohort), sep = "\t"),
  paste("input", "noise", paste(predictors$noise, collapse = " "), sep = "\t"),
  paste("input", "states", paste(apply(states, 1, paste, collapse = ""), collapse = " "), sep = "\t"),
  ## Every branch of the split test, with the tree's whole path behind it.
  tree_case("gate_diss", "real", R = R, pval = 1, weight.permutation = "diss"),
  tree_case("gate_group", "real", R = R, pval = 1, weight.permutation = "group"),
  tree_case("gate_replicate", "integer", R = R, pval = 1, weight.permutation = "replicate"),
  tree_case("gate_rounded", "real", R = R, pval = 1, weight.permutation = "rounded-replicate"),
  tree_case("gate_none_weighted", "real", R = R, pval = 1, weight.permutation = "none"),
  ## A weak association and more permutations: most p-values above are near 1/R, where a wrong
  ## statistic hardly moves them; here they are mid-range, so a count over R is sensitive to it.
  tree_case("weak_replicate", "integer", R = 100, pval = 1, weight.permutation = "replicate", formula = D ~ noise, max.depth = 2),
  tree_case("weak_rounded", "real", R = 100, pval = 1, weight.permutation = "rounded-replicate", formula = D ~ noise, max.depth = 2),
  ## A threshold that gate_replicate's third split fails, so the p-value decides the tree's shape.
  tree_case("gate_rejects", "integer", R = R, pval = 0.05, weight.permutation = "replicate"),
  ## squared, first and the rounding of weights, deterministic with R = 1.
  tree_case("squared_unweighted", "none", R = 1, squared = TRUE),
  tree_case("squared_weighted", "real", R = 1, squared = TRUE, weight.permutation = "diss"),
  tree_case("first_phase", "none", R = 1, first = "phase"),
  tree_case("first_cohort", "real", R = 1, first = "cohort", weight.permutation = "diss"),
  tree_case("first_no_split", "none", R = 1, first = "cohort", min.size = 8),
  tree_case("first_gate", "real", R = R, pval = 1, first = "cohort", weight.permutation = "group"),
  tree_case("rounded_deterministic", "real", R = 1, weight.permutation = "rounded-replicate"),
  tree_case("rounded_to_zero", "rounds_to_zero", R = 1, weight.permutation = "rounded-replicate", formula = D ~ phase),
  ## The same options through seqtree.
  seqtree_case("seqtree_first_diss", "real", R = R, pval = 1, first = "phase", weight.permutation = "diss"),
  seqtree_case("seqtree_squared", "integer", R = 1, squared = TRUE, weight.permutation = "replicate"),
  ## What R refuses.
  refusal("refuse_replicate_real", weights = real, weight.permutation = "replicate"),
  refusal("refuse_first_absent", first = "absent"),
  refusal("refuse_rounded_to_zero", weights = rounds_to_zero, weight.permutation = "rounded-replicate")
)
writeLines(
  c("# disstree and seqtree with the weight permutations, squared and first; R's sample replayed from the perm rows", lines),
  file.path(out_dir, "disstree_weighted_gate.txt")
)
