# Golden generation: disstree / seqtree (Phase 2.4 / DUPLICATES A8).
#
# Anchors traminex's R-faithful weighted regression tree on a dissimilarity
# matrix to TraMineR 2.2-14 `disstree` / `seqtree`. R=1 disables the
# (stochastic) permutation gate so the grown tree is deterministic: node ids,
# weighted sizes, discrepancies, medoids, split statistics, split schedule,
# leaf memberships and the adjustment t0 statistics are all pinned.
#
# Inputs come from make_goldens_phase2_inputs.R (shared Phase 2 fixture);
# the two predictors below are deterministic.
# Requires R + TraMineR 2.2-14; see scripts/r_goldens/README.md.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D
i <- 1:inp$n

pred <- data.frame(
  A = factor(c("a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c")),
  X = ((i * 7) %% 15) + 0.1 * i
)

flatten_tree <- function(tree) {
  rows <- list()
  rec <- function(node) {
    info <- node$info
    sp <- node$split
    if (!is.null(sp)) {
      brk <- if (is.null(sp$breaks)) NaN else as.numeric(sp$breaks)
      tail <- c(
        sp$varindex, sp$info$SCres, sp$info$lpop, sp$info$rpop,
        sp$info$lvar, sp$info$rvar, sp$info$R2, brk
      )
      has_split <- 1
    } else {
      tail <- rep(NaN, 8)
      has_split <- 0
    }
    rows[[length(rows) + 1]] <<- c(
      node$id, info$depth, info$splitschedule, info$n, info$vardis,
      info$medoid, has_split, tail
    )
    if (!is.null(node$kids)) {
      rec(node$kids[[1]])
      rec(node$kids[[2]])
    }
  }
  rec(tree$root)
  do.call(rbind, rows)
}

node_header <- paste0(
  "# columns: id, depth, splitschedule, n, vardis, medoid(1-based), has_split,",
  " varindex(1-based), SCres, lpop, rpop, lvar, rvar, R2, breaks (NaN when absent)\n"
)

## --- unweighted disstree ---------------------------------------------------
tr_unw <- disstree(D ~ A + X, data = pred, R = 1, min.size = 2, max.depth = 4)
write_phase2_csv(
  "disstree_unw_nodes", flatten_tree(tr_unw),
  paste0("# R: disstree(D ~ A + X, data=pred, R=1, min.size=2, max.depth=4), preorder\n", node_header)
)
write_phase2_csv(
  "disstree_unw_fitted", matrix(tr_unw$fitted[, 1], ncol = 1),
  "# R: disstree unweighted -> fitted leaf node ids per case\n"
)

## --- weighted disstree (non-integer weights, 'diss' permutation) -----------
tr_w <- disstree(D ~ A + X,
  data = pred, weights = inp$wreal, R = 1,
  min.size = 2, max.depth = 4, weight.permutation = "diss"
)
write_phase2_csv(
  "disstree_wreal_nodes", flatten_tree(tr_w),
  paste0(
    "# R: disstree(D ~ A + X, data=pred, weights=wreal, R=1, min.size=2, max.depth=4,",
    " weight.permutation='diss'), preorder\n", node_header
  )
)
write_phase2_csv(
  "disstree_wreal_fitted", matrix(tr_w$fitted[, 1], ncol = 1),
  "# R: disstree weighted -> fitted leaf node ids per case\n"
)
write_phase2_csv(
  "disstree_wreal_adjust_t0", matrix(tr_w$info$adjustment$stat[, 1], ncol = 1),
  "# R: disstree weighted -> adjustment dissassoc t0 [Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett, Levene]\n"
)

## --- weighted seqtree (weights read from the sequence object) --------------
states <- matrix("", nrow = inp$n, ncol = 8)
for (r in 1:inp$n) {
  for (cc in 1:8) {
    states[r, cc] <- c("A", "B", "C")[((r + cc) %% 3) + 1]
  }
}
sts <- suppressMessages(seqdef(states, alphabet = c("A", "B", "C"), weights = inp$wint))
tr_s <- seqtree(sts ~ A + X,
  data = pred, diss = D, R = 1, min.size = 2, max.depth = 4,
  weighted = TRUE, weight.permutation = "replicate"
)
write_phase2_csv(
  "seqtree_wint_nodes", flatten_tree(tr_s),
  paste0(
    "# R: seqtree(sts ~ A + X, data=pred, diss=D, R=1, min.size=2, max.depth=4,",
    " weighted=TRUE, weight.permutation='replicate'), weights=wint from seqdef, preorder\n",
    node_header
  )
)
write_phase2_csv(
  "seqtree_wint_fitted", matrix(tr_s$fitted[, 1], ncol = 1),
  "# R: seqtree weighted -> fitted leaf node ids per case\n"
)

cat("DISSTREE_GOLDENS_DONE\n")
