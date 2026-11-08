# Golden generation: seqecmpgroup's chi-square arithmetic, bit for bit
# (2026-09-22 audit S5-16).
#
# seqecmpgroup.R builds xtabs(ww ~ group + presence) -- tapply(w, by, sum) --
# and runs stats::chisq.test on it. chisq.test (R 4.6.1, read from the
# installed stats namespace) computes sr <- rowSums(x), sc <- colSums(x),
# n <- sum(x), E <- outer(sr, sc)/n, one scalar YATES <- min(0.5, abs(x - E))
# for a 2 x 2 table (0 otherwise) and STATISTIC <- sum((abs(x - E) - YATES)^2/E).
# sum, rowSums and colSums accumulate in long double, so statistics that tie
# in exact arithmetic carry R's own last-bit differences (and ~1e-31 residue
# where the exact statistic is 0); seqecmpgroup orders subsequences on them.
#
# chisq_tables -- seeded group/presence/weight vectors, the xtabs cells and
#   chisq.test's statistic, expected counts and residuals.
# cmpgroup_rows -- seqecmpgroup on unweighted and weighted corpora with many
#   tied statistics: every row, in R's order.
#
# Consumed by tests/test_seqecmpgroup_chisq_arithmetic.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

num <- function(x) if (length(x)) paste(sprintf("%.17g", as.vector(x)), collapse = " ") else "none"
lines <- character()
add <- function(...) lines <<- c(lines, paste(..., sep = "\t"))

# ---------------------------------------------------------------------------
# chisq_tables
# ---------------------------------------------------------------------------
set.seed(20260924)
kept <- 0
while (kept < 300) {
  n <- sample(6:16, 1)
  nlev <- sample(c(2, 2, 3), 1)
  g <- sample(c("a", "b", "c")[seq_len(nlev)], n, replace = TRUE)
  p <- sample(0:1, n, replace = TRUE)
  w <- switch(sample(1:3, 1),
    rep(1, n),
    sample(c(0.1, 0.2, 0.7, 1.3), n, replace = TRUE),
    runif(n, 0.05, 4))
  gf <- factor(g)
  tab <- xtabs(w ~ gf + p)
  if (length(dim(tab)) != 2 || ncol(tab) != 2 || nrow(tab) != nlev) next
  if (any(rowSums(tab) == 0) || any(colSums(tab) == 0)) next
  chi <- suppressWarnings(chisq.test(tab))
  kept <- kept + 1
  add("table", kept, paste(as.integer(gf), collapse = " "), paste(p, collapse = " "), num(w),
    num(unclass(tab)), num(chi$statistic), num(chi$expected), num(chi$residuals))
}

# Crafted tables where the summation's precision matters: a cell of
# 2^-53 + 2^-65 next to 1 rounds differently when rowSums adds in long double
# (1) than when a double sum rounds once (1 + 2^-52).
crafted <- list(c(1, 2^-53 + 2^-65, 1, 1), c(2^-53 + 2^-65, 1, 1, 3), c(1, 1, 2^-53 + 2^-65, 1))
for (w in crafted) {
  gf <- factor(c("a", "a", "b", "b"))
  p <- c(0, 1, 0, 1)
  tab <- xtabs(w ~ gf + p)
  chi <- suppressWarnings(chisq.test(tab))
  kept <- kept + 1
  add("table", kept, paste(as.integer(gf), collapse = " "), paste(p, collapse = " "), num(w),
    num(unclass(tab)), num(chi$statistic), num(chi$expected), num(chi$residuals))
}

# r_sum -- sum, rowSums and colSums of matrices whose total depends on the
# order and precision of the additions (column-major, long double).
sums <- list(
  matrix(c(2^-64, 2^-64, 2^-64, 1, 2^-53, 0), nrow = 3),
  matrix(c(1, 2^-64, 2^-53, 2^-64, 0, 2^-64), nrow = 3),
  matrix(c(1, 2^-53 + 2^-65), nrow = 1),
  matrix(c(0.1, 0.2, 0.7, 1.3, 1e-16, 3), nrow = 2)
)
for (k in seq_along(sums)) {
  m <- sums[[k]]
  add("r_sum", k, nrow(m), num(m), num(sum(m)), num(rowSums(m)), num(colSums(m)))
}

# ---------------------------------------------------------------------------
# cmpgroup_rows
# ---------------------------------------------------------------------------
set.seed(20260925)
corpus <- data.frame(id = integer(), timestamp = numeric(), event = character())
for (i in 1:24) {
  k <- sample(2:5, 1)
  corpus <- rbind(corpus, data.frame(id = i, timestamp = sort(sample(0:9, k)),
    event = sample(c("a", "b", "c"), k, replace = TRUE)))
}
for (r in seq_len(nrow(corpus))) {
  add("corpus", corpus$id[r], sprintf("%.17g", corpus$timestamp[r]), corpus$event[r])
}
e <- seqecreate(id = corpus$id, timestamp = corpus$timestamp, event = corpus$event)
weights <- sample(c(0.1, 0.2, 0.7, 1.3), length(e), replace = TRUE)
add("weights", num(weights))
i <- seq_len(length(e))
groupings <- list(two = ifelse(i %% 2 == 0, "x", "y"), three = c("x", "y", "z")[1 + i %% 3])
for (weighted in c(FALSE, TRUE)) {
  if (weighted) seqeweight(e) <- weights
  fs <- seqefsub(e, min.support = 2, max.k = 2)
  for (tag in names(groupings)) {
    res <- seqecmpgroup(fs, group = groupings[[tag]])
    case <- paste(tag, if (weighted) "weighted" else "unweighted", sep = "_")
    sub <- as.character(res$subseq)
    for (r in seq_along(sub)) add("row", case, r, sub[r], num(unlist(res$data[r, ])))
  }
}

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c(
  "# kind\tfields (tab separated); numbers are sprintf('%.17g'), matrices column-major",
  lines
), file.path(out_dir, "audit20260922_s5b_chisq_arithmetic.txt"))
cat("AUDIT20260922_S5B_GOLDENS_DONE\n")
