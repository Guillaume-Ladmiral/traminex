# Live-R goldens for the orientation of the one-dimensional MDS axis (repair of
# TX-REP-20261017-006).
#
# seqrf, dissrf and seqplot's rf sort sequences by cmdscale(diss, k = 1) (wcmdscale
# with weights; dissrf.R:86-91). The axis has no sign of its own: R takes its
# eigenvector from LAPACK, and the sign follows the numerical library (on these 14
# sequences, R with OpenBLAS 0.3.26 gives the negation of SciPy's LAPACK 3.9.0).
# The groups, medoids and R2 then depend on the build. traminex orients the axis
# by a rule: the coordinate of largest absolute value is positive; coordinates
# within a relative 1e-10 of that largest value are tied, and the first of them
# decides; missing coordinates of zero-weight cases are ignored. The tolerance
# matters on symmetric data (line5, line7 and mirror6 below, Hamming distances),
# where two coordinates are equal in exact arithmetic and differ only in rounding,
# so the largest one is the build's choice (batch-5 review I-3). The rows record
# R's raw axis, the oriented axis, and R's seqrf given the oriented axis as sortv:
# the values the rule defines.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_mds_orientation.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AAABBC", "AABBCC", "ABBBCC", "AAAABB", "BBBCCC", "ABCABC", "AACCBB",
             "CCBBAA", "ABBBBB", "AAACCC", "BCCCCC", "AAAAAC", "BBAACC", "CABBBA")
rows_of <- function(s) as.data.frame(do.call(rbind, strsplit(s, "")), stringsAsFactors = FALSE)
weights <- c(1, 2, 1, 1, 0.5, 1, 1, 2, 1, 1, 1, 1, 0, 1)
plain <- suppressMessages(TraMineR::seqdef(rows_of(strings), alphabet = c("A", "B", "C")))
weighted <- suppressMessages(TraMineR::seqdef(rows_of(strings), alphabet = c("A", "B", "C"), weights = weights))
D <- suppressMessages(TraMineR::seqdist(plain, method = "LCS"))

oriented <- function(axis) {
  # the first coordinate within a relative 1e-10 of the largest |value|; NA ignored
  largest <- which(abs(axis) >= max(abs(axis), na.rm = TRUE) * (1 - 1e-10))[1]
  if (axis[largest] < 0) -axis else axis
}
symmetric <- list(
  line5 = c("AAAA", "AAAB", "AABB", "ABBB", "BBBB"),
  line7 = c("AAAAAA", "AAAAAB", "AAAABB", "AAABBB", "AABBBB", "ABBBBB", "BBBBBB"),
  mirror6 = c("AAAB", "BAAA", "AABB", "BBAA", "ABBB", "BBBA")
)
symmetric_seqs <- lapply(symmetric, function(s) suppressMessages(TraMineR::seqdef(rows_of(s), alphabet = c("A", "B"))))
symmetric_diss <- lapply(symmetric_seqs, function(s) suppressMessages(TraMineR::seqdist(s, method = "HAM")))
rows <- list()
add <- function(case, quantity, values) {
  rows[[length(rows) + 1L]] <<- data.frame(case = case, quantity = quantity, position = seq_along(values),
                                           value = format(values, digits = 17), stringsAsFactors = FALSE)
}
axes <- c(list(unweighted = cmdscale(D, k = 1)[, 1],
               weighted = vegan::wcmdscale(D, k = 1, w = weights)[, 1]),
          lapply(symmetric_diss, function(d) cmdscale(d, k = 1)[, 1]))
for (name in names(axes)) {
  add(name, "r_axis", axes[[name]])
  add(name, "oriented_axis", oriented(axes[[name]]))
}
runs <- list(
  list(case = "unweighted_first_k3", seqs = plain, axis = "unweighted", k = 3, grp = "first"),
  list(case = "unweighted_first_k4", seqs = plain, axis = "unweighted", k = 4, grp = "first"),
  list(case = "unweighted_prop_k3", seqs = plain, axis = "unweighted", k = 3, grp = "prop"),
  list(case = "unweighted_prop_k4", seqs = plain, axis = "unweighted", k = 4, grp = "prop"),
  list(case = "weighted_prop_k3", seqs = weighted, axis = "weighted", k = 3, grp = "prop")
)
for (name in names(symmetric)) {
  runs[[length(runs) + 1L]] <- list(case = paste0(name, "_first_k2"), seqs = symmetric_seqs[[name]],
                                    diss = symmetric_diss[[name]], axis = name, k = 2, grp = "first")
}
for (run in runs) {
  diss <- if (is.null(run$diss)) D else run$diss
  rf <- suppressMessages(TraMineR::seqrf(run$seqs, diss = diss, k = run$k, sortv = oriented(axes[[run$axis]]),
                                         grp.meth = run$grp))$rf
  add(run$case, "medoids", rf$medoids)
  add(run$case, "heights", rf$heights)
  add(run$case, "R2", rf$R2)
  add(run$case, "Fstat", rf$Fstat)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "mds_orientation.csv"), row.names = FALSE, quote = TRUE)
