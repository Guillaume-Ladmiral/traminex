# Live-R goldens for degenerate OMstran and CHI2/EUCLID inputs (TX-AUD-20261019-080).
#
# R's OMstran drops columns without drop = FALSE (seqdist-OMstran.R:55-58): on one sequence the matrix becomes a
# vector (a 4 x 4 result, or "incorrect number of dimensions"), and it cannot keep zero columns (length-1 data
# without add.column, length-2 data with previous). transindel = "prob" needs seqtrate, which fails on length-1
# data. R's CHI2 does not cap its first overlapping episode c(1, 1 + step/2) at the sequence length
# (seqdist-CHI2.R:25-27), while it caps the others with pmin (:33). Each row records R's result (its dimensions
# and values, column by column) or its error; the step = 1 rows anchor the corrected CHI2/EUCLID values.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_degenerate_omstran.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

one <- suppressMessages(TraMineR::seqdef(rbind(c("A", "A", "B", "B", "C")), alphabet = c("A", "B", "C")))
len1 <- suppressMessages(TraMineR::seqdef(rbind("A", "B", "A", "C", "B"), alphabet = c("A", "B", "C")))
two <- suppressMessages(TraMineR::seqdef(rbind(c("A", "B"), c("B", "A"), c("A", "A")), alphabet = c("A", "B")))
cases <- list(
  one_constant_noprev_noadd = quote(TraMineR::seqdist(one, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "constant", previous = FALSE, add.column = FALSE)),
  one_constant_prev_noadd = quote(TraMineR::seqdist(one, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "constant", previous = TRUE, add.column = FALSE)),
  len1_constant_noprev_noadd = quote(TraMineR::seqdist(len1, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "constant", previous = FALSE, add.column = FALSE)),
  len1_constant_prev_noadd = quote(TraMineR::seqdist(len1, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "constant", previous = TRUE, add.column = FALSE)),
  len1_prob_noprev_add = quote(TraMineR::seqdist(len1, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "prob", previous = FALSE, add.column = TRUE)),
  two_constant_prev_noadd = quote(TraMineR::seqdist(two, method = "OMstran", sm = "CONSTANT", otto = 0.3, transindel = "constant", previous = TRUE, add.column = FALSE)),
  len1_chi2_step2_overlap = quote(TraMineR::seqdist(len1, method = "CHI2", step = 2, overlap = TRUE)),
  len1_euclid_step2_overlap = quote(TraMineR::seqdist(len1, method = "EUCLID", step = 2, overlap = TRUE)),
  len1_chi2_step1 = quote(TraMineR::seqdist(len1, method = "CHI2", step = 1)),
  len1_euclid_step1 = quote(TraMineR::seqdist(len1, method = "EUCLID", step = 1)),
  len1_chi2_step1_auto = quote(TraMineR::seqdist(len1, method = "CHI2", step = 1, norm = "auto")),
  len1_euclid_step1_auto = quote(TraMineR::seqdist(len1, method = "EUCLID", step = 1, norm = "auto"))
)
rows <- list()
for (name in names(cases)) {
  err <- ""
  dims <- ""
  value <- ""
  ok <- tryCatch({
    d <- as.matrix(suppressMessages(suppressWarnings(eval(cases[[name]]))))
    dims <- paste(dim(d), collapse = "x")
    value <- paste(format(as.vector(d), digits = 15), collapse = " ")
    TRUE
  }, error = function(e) { err <<- conditionMessage(e); FALSE })
  rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, error = err, dim = dims, value = value,
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "degenerate_omstran.csv"), row.names = FALSE, quote = TRUE)
