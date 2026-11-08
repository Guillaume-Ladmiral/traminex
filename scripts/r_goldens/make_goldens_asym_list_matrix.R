stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for OM with an asymmetric substitution matrix on sequences of different
# lengths (audit finding TX-AUD-20261019-003). R applies sm[i, j] with i the state of the
# first sequence (OMdistance) and warns "Non-symmetric substitution cost matrix!".
# traminex accepted the same matrix as nested lists but applied it the other way round on
# this route and without the warning; the Python test passes R's matrix as an array, as
# nested lists and as tuples.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_asym_list_matrix.R

suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- data.frame(a = c("A", "B", "A", "C", "B"), b = c("B", "A", "B", "C", "A"),
                c = c("C", NA, "C", "B", NA), e = c(NA, NA, NA, "A", NA))
s <- suppressMessages(TraMineR::seqdef(d, alphabet = c("A", "B", "C"), right = "DEL"))
sm <- matrix(c(0, 1, 5, 3, 0, 1, 1, 4, 0), 3, 3, byrow = TRUE)

rows <- list()
for (case in list(list(name = "indel_scalar", indel = 2), list(name = "indel_vector", indel = c(1, 2, 1.5)))) {
  warned <- character()
  D <- withCallingHandlers(
    suppressMessages(TraMineR::seqdist(s, method = "OM", sm = sm, indel = case$indel)),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  D <- as.matrix(D)
  idx <- which(row(D) >= 1, arr.ind = TRUE)
  rows[[length(rows) + 1L]] <- data.frame(case = case$name, i = idx[, 1], j = idx[, 2],
                                          distance = sprintf("%.17g", D[idx]),
                                          warning = paste(warned, collapse = " | "))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "asym_list_matrix.csv"), row.names = FALSE, quote = TRUE)
