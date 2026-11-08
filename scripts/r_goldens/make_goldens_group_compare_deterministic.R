# Golden generation: the DETERMINISTIC components of dissCompare / seqCompare /
# seqLRT / seqBIC, on a corpus where the options actually change the answer.
#
# Cycle-7 W-02. These four are stochastic at their default s = 100, which is why
# earlier cycles pinned them only lightly. Two deterministic surfaces remain,
# and both are directly comparable to R at full precision:
#
#   * s = 0 disables sampling entirely (dissCompare.R:198, seqCompare.R:198):
#     no set.seed, no permute, no sample. Every returned statistic is exact.
#   * seqxcomp (seqCompare.R:312) is a pure function of (r1, r2, diss, weights).
#
# WHY THIS GENERATOR EXISTS ALONGSIDE make_goldens_group_compare_r3.R. That one
# mines a CONSTANT dissimilarity matrix (all off-diagonal distances equal 1) so
# that its sampling cases are independent of R-versus-NumPy random bytes. The
# side effect is that its deterministic cases cannot discriminate the options
# they name: its s0_squared_false and s0_squared_true rows are byte-identical.
# This generator uses a varied matrix from real sequences, so squared, weighted,
# stat, inc.df and set each move the numbers.
#
# Every value is written with format(digits = 17); nothing is read off a print.
# Requires live R + TraMineR 2.2-14 / TraMineRextras 0.6.9.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(packageVersion("TraMineR") == "2.2.14",
          packageVersion("TraMineRextras") == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fmt <- function(x) format(x, digits = 17, scientific = FALSE)

## --------------------------------------------------------------------- ##
## Corpus. Fixed seed, and the codes are written out so the Python side   ##
## reads the same matrix rather than reproducing R's RNG.                 ##
## --------------------------------------------------------------------- ##
set.seed(20260830)
n_seq <- 18L
seq_len_ <- 8L
alphabet <- c("A", "B", "C", "D")
codes <- matrix(sample(seq_along(alphabet), n_seq * seq_len_, replace = TRUE),
                nrow = n_seq, ncol = seq_len_)
states <- matrix(alphabet[codes], nrow = n_seq, ncol = seq_len_)
weights <- round(seq(0.4, 3.1, length.out = n_seq), 3)
group <- rep(c("g1", "g2"), length.out = n_seq)
group[c(3L, 4L, 11L)] <- "g2"           # unbalanced: 7 vs 11
sets <- c(rep("early", 10L), rep("late", 8L))

sdata <- suppressMessages(
  seqdef(as.data.frame(states, stringsAsFactors = FALSE),
         alphabet = alphabet, weights = weights)
)

# 0-based codes so Python can rebuild the identical STSList.
write.csv(as.data.frame(codes - 1L),
          file.path(out_dir, "group_compare_det_codes.csv"), row.names = FALSE)
writeLines(alphabet, file.path(out_dir, "group_compare_det_alphabet.txt"))
writeLines(fmt(weights), file.path(out_dir, "group_compare_det_weights.txt"))
writeLines(group, file.path(out_dir, "group_compare_det_group.txt"))
writeLines(sets, file.path(out_dir, "group_compare_det_sets.txt"))

# A varied dissimilarity matrix -- the whole point of this generator.
D <- suppressMessages(seqdist(sdata, method = "OM", sm = "TRATE", indel = 1.0))
stopifnot(length(unique(as.vector(D))) > 10)   # guard against a constant matrix
write.csv(as.data.frame(D), file.path(out_dir, "group_compare_det_diss.csv"),
          row.names = FALSE)

## --------------------------------------------------------------------- ##
## 1. seqxcomp -- the worker, plus the intermediates it is built from.    ##
## --------------------------------------------------------------------- ##
# Fixed index sets: no sampling is involved, so this is exact.
r1 <- which(group == "g1")
r2 <- which(group == "g2")
writeLines(as.character(r1 - 1L), file.path(out_dir, "group_compare_det_r1.txt"))
writeLines(as.character(r2 - 1L), file.path(out_dir, "group_compare_det_r2.txt"))

worker_rows <- character(0)
for (sq in c(FALSE, TRUE)) {
  for (lrtpow in c(1L, 2L)) {
    for (wmode in c("none", "global", "by.group")) {
      weighted <- wmode != "none"
      w <- if (weighted) weights else NULL
      value <- TraMineRextras:::seqxcomp(
        r1, r2, D, w, is.LRT = TRUE, is.BIC = TRUE,
        squared = sq, weighted = weighted,
        weight.by = if (wmode == "by.group") "by.group" else "global",
        LRTpow = lrtpow
      )
      # The intermediates the statistic is built from, recomputed exactly as
      # seqxcomp does (seqCompare.R:360-366), so a divergence can be localised.
      n1 <- length(r1); n2 <- length(r2); n0 <- n1 + n2
      if (weighted) {
        w1 <- weights[r1]; w2 <- weights[r2]
        if (wmode == "by.group") {
          w1 <- n1 / sum(w1) * w1
          w2 <- n2 / sum(w2) * w2
        }
        wc <- c(w1, w2)
      } else {
        w1 <- rep(1, n1); w2 <- rep(1, n2); wc <- rep(1, n0)
      }
      dS  <- disscenter(D[c(r1, r2), c(r1, r2)], weights = wc, squared = sq)
      dS1 <- disscenter(D[r1, r1], weights = w1, squared = sq)
      dS2 <- disscenter(D[r2, r2], weights = w2, squared = sq)
      SS  <- sum(wc * dS^lrtpow)
      SS1 <- sum(w1 * dS1^lrtpow)
      SS2 <- sum(w2 * dS2^lrtpow)
      worker_rows <- c(worker_rows, sprintf(
        "%s,%d,%s,%s,%s,%s,%s,%s,%s,%s,%d",
        sq, lrtpow, wmode,
        fmt(value[1, 1]), fmt(value[1, 2]), fmt(value[1, 3]), fmt(value[1, 4]),
        fmt(SS), fmt(SS1), fmt(SS2), n0
      ))
    }
  }
}
writeLines(c("# squared,LRTpow,weight_mode,LRT,p,BIC,BF,SS,SS1,SS2,n0", worker_rows),
           file.path(out_dir, "group_compare_det_seqxcomp.csv"))

## --------------------------------------------------------------------- ##
## 2. The full deterministic (s = 0) option grid, for all four entries.   ##
## --------------------------------------------------------------------- ##
emit <- function(store, case, result) {
  if (is.null(dim(result))) {
    result <- matrix(result, nrow = 1, dimnames = list("1", names(result)))
  }
  rn <- if (is.null(rownames(result))) as.character(seq_len(nrow(result))) else rownames(result)
  for (i in seq_len(nrow(result))) {
    for (j in seq_len(ncol(result))) {
      store[[length(store) + 1L]] <- data.frame(
        case = case, row = rn[[i]], column = colnames(result)[[j]],
        value = fmt(result[i, j]), stringsAsFactors = FALSE
      )
    }
  }
  store
}

store <- list()
for (stat in c("all", "LRT", "BIC")) {
  for (sqname in c("FALSE", "TRUE", "LRTonly")) {
    sq <- if (sqname == "LRTonly") "LRTonly" else as.logical(sqname)
    for (wmode in c("FALSE", "TRUE", "by.group")) {
      wtd <- if (wmode == "by.group") "by.group" else as.logical(wmode)
      for (incdf in c(FALSE, TRUE)) {
        case <- sprintf("diss|stat=%s|squared=%s|weighted=%s|inc.df=%s",
                        stat, sqname, wmode, incdf)
        store <- emit(store, case, suppressMessages(dissCompare(
          D, group, s = 0, stat = stat, squared = sq,
          weighted = wtd, weights = weights, inc.df = incdf
        )))
      }
    }
  }
}

# Sets: two comparisons in one call, so row labels and per-set sizes are pinned.
#
# UPSTREAM DEFECT, recorded not copied. dissCompare.R:119-120 computes
#   idx.a[[i]] <- which(group[set==lev.set[i]] == lev.g[1])
# which are positions WITHIN the set, and then indexes the FULL matrix with
# them at line 179: diss[c(idx.a[[i]], idx.b[[i]]), ...]. For the first set
# level those positions coincide with absolute rows, so it is right by
# accident; for every later level R compares the wrong cases. Here the "late"
# level's rows are 11..18 but R reads rows 1,2,3,4,5,6,7,8. seqCompare.R:131
# does NOT share this -- it subsets with a boolean mask.
#
# So the buggy output is recorded under diss_sets|..., and the CORRECT per-set
# answer is derived from R separately by calling dissCompare on the explicitly
# subsetted matrix, where only one set exists and the defective path is not
# taken.
for (stat in c("all", "LRT", "BIC")) {
  case <- sprintf("diss_sets|stat=%s", stat)
  store <- emit(store, case, suppressMessages(dissCompare(
    D, group, set = sets, s = 0, stat = stat, weights = weights, inc.df = TRUE
  )))
}
for (lev in levels(factor(sets))) {
  rows <- which(sets == lev)
  for (stat in c("all", "LRT", "BIC")) {
    store <- emit(store, sprintf("diss_set_alone|set=%s|stat=%s", lev, stat),
                  suppressMessages(dissCompare(
                    D[rows, rows], group[rows], s = 0, stat = stat,
                    weights = weights[rows], inc.df = TRUE
                  )))
  }
}

# Sequence route: the distance is computed inside, over several methods.
for (method in c("LCS", "OM", "HAM")) {
  for (wmode in c("FALSE", "TRUE", "by.group")) {
    wtd <- if (wmode == "by.group") "by.group" else as.logical(wmode)
    # traminex spells R's `sm` as `subm` (docs/API_CONSISTENCY.md section 1b).
    common <- list(sdata, group = group, s = 0, weighted = wtd, method = method)
    if (method == "OM") common <- c(common, list(sm = "TRATE", indel = 1.0))
    store <- emit(store, sprintf("seqCompare|method=%s|weighted=%s", method, wmode),
                  suppressMessages(do.call(seqCompare, c(common, list(stat = "all")))))
    store <- emit(store, sprintf("seqLRT|method=%s|weighted=%s", method, wmode),
                  suppressMessages(do.call(seqLRT, common)))
    store <- emit(store, sprintf("seqBIC|method=%s|weighted=%s", method, wmode),
                  suppressMessages(do.call(seqBIC, common)))
  }
}

# Two-object and set-labelled sequence routes.
store <- emit(store, "seqCompare|two_object",
              suppressMessages(seqCompare(sdata[group == "g1", ], sdata[group == "g2", ],
                                          s = 0, method = "LCS", stat = "all")))
store <- emit(store, "seqCompare|sets",
              suppressMessages(seqCompare(sdata, group = group, set = sets, s = 0,
                                          method = "LCS", stat = "all", inc.df = TRUE)))

grid <- do.call(rbind, store)
write.csv(grid, file.path(out_dir, "group_compare_det_grid.csv"), row.names = FALSE)

## --------------------------------------------------------------------- ##
## 3. Shapes and labels of the SAMPLING route.                            ##
## These do not depend on the draws, so they are comparable even though   ##
## the values are not.                                                    ##
## --------------------------------------------------------------------- ##
shape_rows <- character(0)
shape <- function(tag, result) {
  if (is.null(dim(result))) {
    result <- matrix(result, nrow = 1, dimnames = list("1", names(result)))
  }
  rn <- if (is.null(rownames(result))) "" else paste(rownames(result), collapse = "|")
  sprintf("%s,%d,%d,%s,%s", tag, nrow(result), ncol(result),
          paste(colnames(result), collapse = "|"), rn)
}
for (bf in list(NULL, 1, 2)) {
  tag <- sprintf("sample|BFopt=%s", if (is.null(bf)) "NULL" else bf)
  shape_rows <- c(shape_rows, shape(tag, suppressMessages(
    dissCompare(D, group, s = 4, seed = 99, BFopt = bf, weights = weights))))
  shape_rows <- c(shape_rows, shape(paste0(tag, "|inc.df"), suppressMessages(
    dissCompare(D, group, s = 4, seed = 99, BFopt = bf, weights = weights,
                inc.df = TRUE))))
}
for (stat in c("all", "LRT", "BIC")) {
  shape_rows <- c(shape_rows, shape(sprintf("sample|stat=%s", stat), suppressMessages(
    dissCompare(D, group, s = 4, seed = 99, stat = stat, weights = weights))))
  shape_rows <- c(shape_rows, shape(sprintf("s0|stat=%s", stat), suppressMessages(
    dissCompare(D, group, s = 0, stat = stat, weights = weights))))
}
shape_rows <- c(shape_rows, shape("sample|sets", suppressMessages(
  dissCompare(D, group, set = sets, s = 4, seed = 99, weights = weights,
              inc.df = TRUE))))
writeLines(c("# tag,nrow,ncol,columns,rownames", shape_rows),
           file.path(out_dir, "group_compare_det_shapes.csv"))

## --------------------------------------------------------------------- ##
## 4. Group sizes per set, so a mis-split is caught before the statistic. ##
## --------------------------------------------------------------------- ##
size_rows <- character(0)
for (lev in levels(factor(sets))) {
  in_set <- sets == lev
  size_rows <- c(size_rows, sprintf("%s,%d,%d", lev,
                                    sum(in_set & group == "g1"),
                                    sum(in_set & group == "g2")))
}
size_rows <- c(size_rows, sprintf("ALL,%d,%d", sum(group == "g1"), sum(group == "g2")))
writeLines(c("# set,n_g1,n_g2", size_rows),
           file.path(out_dir, "group_compare_det_sizes.csv"))

cat("GROUP_COMPARE_DETERMINISTIC_GOLDENS_DONE\n")
cat("grid rows:", nrow(grid), " worker rows:", length(worker_rows),
    " shape rows:", length(shape_rows), "\n")
