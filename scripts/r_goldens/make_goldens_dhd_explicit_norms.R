# Live-R goldens for DHD normalisation and explicit-cost checks (repair of audit
# findings TX-AUD-20261010-012 and -023, the low-level pairwise DHD routes).
#
# R checks an explicit DHD array with checkcost (TraMineR-checkcost.R:40-75: size,
# negative costs, NaN, diagonal) and normalises HAM and DHD by the sum over
# positions of the largest cost (seqdist.R:645-649, 771-776); norm "auto" is
# "maxlength" for DHD (seqdist.R:317-319).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_dhd_explicit_norms.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("ABBCDA", "BBCCDD", "DACBBA", "AADDCC", "CBADCB")
data <- as.data.frame(do.call(rbind, lapply(strings, function(s) strsplit(s, "")[[1]])), stringsAsFactors = FALSE)
s5 <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C", "D")))
sm5 <- abs(outer(1:4, 1:4, "-")) * 1.0
arr <- array(0, dim = c(4, 4, 6))
for (t in 1:6) arr[, , t] <- sm5 * t

norms <- c("none", "maxlength", "gmean", "maxdist", "YujianBo", "auto")
cases <- list()
for (nm in norms) {
  cases[[paste0("explicit_", nm)]] <- list(method = "DHD", sm = arr, norm = nm)
  cases[[paste0("default_", nm)]] <- list(method = "DHD", norm = nm)
}
cases[["ham_sm5_maxlength"]] <- list(method = "HAM", sm = sm5, norm = "maxlength")

rows <- list()
for (case_name in names(cases)) {
  value <- as.matrix(suppressMessages(do.call(TraMineR::seqdist, c(list(seqdata = s5), cases[[case_name]]))))
  rows[[length(rows) + 1L]] <- data.frame(
    case = case_name,
    row = rep(seq_len(nrow(value)), times = ncol(value)),
    col = rep(seq_len(ncol(value)), each = nrow(value)),
    value = sprintf("%.15g", as.vector(value)),
    stringsAsFactors = FALSE
  )
}
write.csv(do.call(rbind, rows), file.path(out_dir, "dhd_explicit_norms.csv"), row.names = FALSE, quote = TRUE)

negative <- arr; negative[1, 2, 1] <- -1; negative[2, 1, 1] <- -1
not_a_number <- arr; not_a_number[1, 2, 1] <- NaN; not_a_number[2, 1, 1] <- NaN
diagonal <- arr; diagonal[1, 1, 1] <- 1
invalid <- list(negative = negative, nan = not_a_number, diagonal = diagonal, undersized = arr[1:3, 1:3, ])
messages <- vapply(names(invalid), function(kind) {
  result <- tryCatch({
    suppressMessages(TraMineR::seqdist(s5, method = "DHD", sm = invalid[[kind]]))
    "no error"
  }, error = function(e) conditionMessage(e))
  paste0(kind, "|", trimws(gsub("\n", " ", result)))
}, character(1))
writeLines(messages, file.path(out_dir, "dhd_explicit_norms_errors.txt"))
