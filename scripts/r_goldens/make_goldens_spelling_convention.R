# Live-R goldens for the spellings traminex accepts beside R's (owner decision of 2026-10-22, option B of triage R9-g;
# TX-AUD-20261019-061 and TX-AUD-20261010-074).
#
# traminex's seqdist reads method, sm and norm in any letter case and takes the norm aliases NULL (None), "length"
# and "max"; its seqMD reads method and norm the same way (its norm goes through seqdist; batch-10 follow-up review
# I-1). Each row is one such spelling: whether R accepts it ("ok"), R's
# error if not, and R's value for the canonical spelling ("canonical", the whole matrix column by column), which
# traminex must give for the non-canonical one.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_spelling_convention.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

x <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "A", "C", "B"), b = c("A", "B", "C", "B"),
                                                  c = c("B", "B", "B", "C"), d = c("C", "C", "A", "A")),
                                       alphabet = c("A", "B", "C")))
cases <- list(
  method_om = list(quote(TraMineR::seqdist(x, method = "om", sm = "CONSTANT")),
                   quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT"))),
  method_ham = list(quote(TraMineR::seqdist(x, method = "ham")), quote(TraMineR::seqdist(x, method = "HAM"))),
  method_lcs = list(quote(TraMineR::seqdist(x, method = "lcs")), quote(TraMineR::seqdist(x, method = "LCS"))),
  sm_trate = list(quote(TraMineR::seqdist(x, method = "OM", sm = "trate")),
                  quote(TraMineR::seqdist(x, method = "OM", sm = "TRATE"))),
  sm_constant = list(quote(TraMineR::seqdist(x, method = "OM", sm = "constant")),
                     quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT"))),
  norm_upper_case = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "MAXLENGTH")),
                         quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "maxlength"))),
  norm_length = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "length")),
                     quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "maxlength"))),
  norm_max = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "max")),
                  quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "maxdist"))),
  norm_null = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = NULL)),
                   quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "none"))),
  norm_yujianbo = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "yujianbo")),
                       quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "YujianBo"))),
  seqmd_method_om = list(quote(TraMineR::seqMD(list(x, x), method = "om", sm = list("CONSTANT", "CONSTANT"), what = "diss")),
                         quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), what = "diss"))),
  # Review m-7: more method, sm and norm spellings.
  method_omspell = list(quote(TraMineR::seqdist(x, method = "omspell", sm = "CONSTANT")),
                        quote(TraMineR::seqdist(x, method = "OMspell", sm = "CONSTANT"))),
  method_twed = list(quote(TraMineR::seqdist(x, method = "twed", sm = "CONSTANT", nu = 0.5)),
                     quote(TraMineR::seqdist(x, method = "TWED", sm = "CONSTANT", nu = 0.5))),
  method_mixed_case = list(quote(TraMineR::seqdist(x, method = "Om", sm = "CONSTANT")),
                           quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT"))),
  sm_indels = list(quote(TraMineR::seqdist(x, method = "OM", sm = "indels")),
                   quote(TraMineR::seqdist(x, method = "OM", sm = "INDELS"))),
  sm_indelslog = list(quote(TraMineR::seqdist(x, method = "OM", sm = "indelslog")),
                      quote(TraMineR::seqdist(x, method = "OM", sm = "INDELSLOG"))),
  norm_gmean_upper = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "GMEAN")),
                          quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "gmean"))),
  norm_auto_upper = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "AUTO")),
                         quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "auto"))),
  norm_none_upper = list(quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "NONE")),
                         quote(TraMineR::seqdist(x, method = "OM", sm = "CONSTANT", norm = "none"))),
  # Review I-1: seqMD's norm.
  seqmd_norm_upper_case = list(quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "GMEAN", what = "diss")),
                               quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "gmean", what = "diss"))),
  seqmd_norm_length = list(quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "length", what = "diss")),
                           quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "maxlength", what = "diss"))),
  seqmd_norm_max = list(quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "max", what = "diss")),
                        quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "maxdist", what = "diss"))),
  seqmd_norm_null = list(quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = NULL, what = "diss")),
                         quote(TraMineR::seqMD(list(x, x), method = "OM", sm = list("CONSTANT", "CONSTANT"), norm = "none", what = "diss")))
)
rows <- list()
for (name in names(cases)) {
  err <- ""
  ok <- tryCatch({ suppressMessages(suppressWarnings(eval(cases[[name]][[1]]))); TRUE },
                 error = function(e) { err <<- conditionMessage(e); FALSE })
  d <- as.matrix(suppressMessages(suppressWarnings(eval(cases[[name]][[2]]))))
  rows[[length(rows) + 1L]] <- data.frame(case = name, ok = ok, error = err,
                                          canonical = paste(format(as.vector(d), digits = 17), collapse = " "),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "spelling_convention.csv"), row.names = FALSE, quote = TRUE)
