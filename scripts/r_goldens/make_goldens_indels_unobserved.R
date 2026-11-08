stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for INDELS costs with an alphabet state that no sequence contains
# (audit finding TX-AUD-20261019-002).
#
# seqcost's INDELS indel is 1/frequency, and a state of frequency 0 gets the stand-in
# .Machine$double.xmax (seqcost.R). That cost never enters a distance between the data's
# sequences, but it enters the normalization bound -- OM's min(max sm, 2 * max indel)
# (OMdistance.cpp:23-33), HAM's and DHD's sum of the largest cost per position
# (seqdist.R:644-649) -- and OMloc's and OMstran's costs, so R gives 0 (or NaN) for every
# normalized distance. traminex corrects this (a DIVERGENCE row of seqdist): such a state
# changes no distance. The expected values are R's own results with that state left out of
# the alphabet (alphabet A, B); in this data A and B both occur at every position, so the
# restricted alphabet has no stand-in cost at all. R's values with the state (alphabet A, B, C)
# are recorded too, to keep the defect visible.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_indels_unobserved.R

suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

m <- rbind(c("A", "A", "B", "B"), c("A", "B", "A", "B"), c("B", "B", "B", "A"), c("A", "A", "A", "A"))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
calls <- list(
  list(method = "OM", norms = c("none", "auto", "maxlength", "gmean", "maxdist", "YujianBo"), args = list()),
  list(method = "HAM", norms = c("none", "auto", "maxlength", "gmean", "maxdist", "YujianBo"), args = list()),
  list(method = "DHD", norms = c("none", "auto", "maxlength", "gmean", "maxdist", "YujianBo"), args = list()),
  list(method = "OMloc", norms = c("none", "auto"), args = list()),
  list(method = "OMslen", norms = c("none", "auto"), args = list()),
  list(method = "OMspell", norms = c("none", "auto"), args = list()),
  list(method = "TWED", norms = c("none", "auto"), args = list(h = 0.5, nu = 0.5)),
  list(method = "OMstran", norms = c("none", "auto"), args = list(otto = 0.5))
)
rows <- list()
for (alphabet_case in list(list(name = "observed", alphabet = c("A", "B")),
                           list(name = "with_unobserved", alphabet = c("A", "B", "C")))) {
  x <- quiet(TraMineR::seqdef(m, alphabet = alphabet_case$alphabet))
  for (call in calls) {
    for (norm in call$norms) {
      args <- c(list(seqdata = x, method = call$method, sm = "INDELS", norm = norm), call$args)
      d <- tryCatch(as.matrix(quiet(do.call(TraMineR::seqdist, args))), error = function(e) conditionMessage(e))
      if (is.character(d)) {
        rows[[length(rows) + 1L]] <- data.frame(alphabet = alphabet_case$name, method = call$method,
                                                norm = norm, i = NA, j = NA, value = paste("ERROR:", d))
        next
      }
      idx <- which(upper.tri(d), arr.ind = TRUE)
      idx <- idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
      rows[[length(rows) + 1L]] <- data.frame(alphabet = alphabet_case$name, method = call$method,
                                              norm = norm, i = idx[, 1], j = idx[, 2],
                                              value = sprintf("%.17g", d[idx]))
    }
  }
}
# A state the external reference contains is compared, so its cost enters the distances:
# R's raw OM distances to a reference made of C alone (C absent from the data).
x <- quiet(TraMineR::seqdef(m, alphabet = c("A", "B", "C")))
ref <- quiet(TraMineR::seqdef(rbind(c("C", "C", "C", "C")), alphabet = c("A", "B", "C")))
d <- as.numeric(quiet(TraMineR::seqdist(x, method = "OM", sm = "INDELS", refseq = ref)))
rows[[length(rows) + 1L]] <- data.frame(alphabet = "reference_contains_it", method = "OM", norm = "none",
                                        i = seq_along(d), j = 0L, value = sprintf("%.17g", d))
# An external reference sequence and two index sets reach the distances by their own routes
# (batch-7 review I-2): OM, HAM and DHD to the reference A-B-B-A, and DHD between the sets
# {1, 2} and {3, 4}, for both alphabets. The reference contains no C either.
for (alphabet_case in list(list(name = "observed", alphabet = c("A", "B")),
                           list(name = "with_unobserved", alphabet = c("A", "B", "C")))) {
  x <- quiet(TraMineR::seqdef(m, alphabet = alphabet_case$alphabet))
  ref <- quiet(TraMineR::seqdef(rbind(c("A", "B", "B", "A")), alphabet = alphabet_case$alphabet))
  for (method in c("OM", "HAM", "DHD")) {
    for (norm in c("none", "auto", "maxdist")) {
      d <- as.numeric(quiet(TraMineR::seqdist(x, method = method, sm = "INDELS", norm = norm, refseq = ref)))
      rows[[length(rows) + 1L]] <- data.frame(alphabet = paste0("external_reference_", alphabet_case$name),
                                              method = method, norm = norm, i = seq_along(d), j = 0L,
                                              value = sprintf("%.17g", d))
    }
  }
  for (norm in c("none", "maxdist")) {
    d <- as.matrix(quiet(TraMineR::seqdist(x, method = "DHD", sm = "INDELS", norm = norm,
                                           refseq = list(1:2, 3:4))))
    rows[[length(rows) + 1L]] <- data.frame(alphabet = paste0("index_sets_", alphabet_case$name),
                                            method = "DHD", norm = norm, i = rep(1:2, times = 2),
                                            j = rep(3:4, each = 2), value = sprintf("%.17g", as.numeric(d)))
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "indels_unobserved.csv"), row.names = FALSE, quote = TRUE)
