# Live-R record of seqdist's response to NaN method parameters (repair of audit
# finding TX-AUD-20261010-013).
#
# R refuses a NaN nu and h (TWED), expcost (OMspell, OMloc), context (OMloc),
# h (OMslen) and otto (OMstran), only through its internal error "missing value
# where TRUE/FALSE needed" (seqdist.R:213-214, 218-223, 235-236, 242-245,
# 277-283). It computes with tpow = NaN (NaN distances) and nu = Inf (infinite
# distances); those two results are recorded as values.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdist_nan_parameters.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AABBC", "ABBCC", "CCBAA", "BBCCA")
data <- as.data.frame(do.call(rbind, lapply(strings, function(s) strsplit(s, "")[[1]])), stringsAsFactors = FALSE)
s4 <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C")))

cases <- list(
  twed_nu = list(method = "TWED", sm = "CONSTANT", indel = 1, nu = NaN, h = 0.5),
  twed_h = list(method = "TWED", sm = "CONSTANT", indel = 1, nu = 0.5, h = NaN),
  omspell_expcost = list(method = "OMspell", sm = "CONSTANT", indel = 1, expcost = NaN),
  omloc_expcost = list(method = "OMloc", sm = "CONSTANT", expcost = NaN),
  omloc_context = list(method = "OMloc", sm = "CONSTANT", expcost = 0.25, context = NaN),
  omslen_h = list(method = "OMslen", sm = "CONSTANT", indel = 1, h = NaN),
  omstran_otto = list(method = "OMstran", sm = "CONSTANT", indel = 1, otto = NaN),
  omspell_tpow = list(method = "OMspell", sm = "CONSTANT", indel = 1, tpow = NaN),
  twed_nu_inf = list(method = "TWED", sm = "CONSTANT", indel = 1, nu = Inf, h = 0.5)
)
lines <- vapply(names(cases), function(case_name) {
  result <- tryCatch({
    value <- suppressMessages(do.call(TraMineR::seqdist, c(list(seqdata = s4), cases[[case_name]])))
    paste0("values|", paste(format(as.vector(as.matrix(value)[1, ]), digits = 15), collapse = " "))
  }, error = function(e) paste0("error|", trimws(gsub("\n", " ", conditionMessage(e)))))
  paste0(case_name, "|", result)
}, character(1))
writeLines(lines, file.path(out_dir, "seqdist_nan_parameters.txt"))
