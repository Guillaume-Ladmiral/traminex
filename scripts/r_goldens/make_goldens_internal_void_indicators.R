# Live-R goldens for per-sequence indicators on sequences with an internal void
# (record of audit finding TX-AUD-20261010-055).
#
# seqrecode(otherwise = void) can leave void cells inside a sequence. R's
# indicators read the first seqlength cells of each row (seqdss.R:12,37;
# seqdur.R:11,29-31), so after an internal void they count void cells and drop
# the states observed at the end. traminex skips void wherever it is. The same
# functions are recorded on the void-free sequences (twin), where R is right, so
# the definition's values come from R too.
#
# seqformat's SPS and DSS forms set the nr and void codes alike to missing
# (seqformat.R:120-123,197-199) and rebuild them with seqprep, so an internal
# void is written as a missing spell. traminex's seqformat follows R there; its
# string writers to_SPS and to_DSS leave void out, as the twin's forms do.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_internal_void_indicators.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

base <- suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(c("ABCCAB", "ABAACC"), ""))),
                                          alphabet = c("A", "B", "C")))
recoded <- suppressMessages(TraMineR::seqrecode(base, recodes = list(p = "A", n = "B"),
                                                otherwise = attr(base, "void")))
values <- function(expr) {
  result <- tryCatch(suppressWarnings(suppressMessages(eval(expr))), error = function(e) e)
  if (inherits(result, "error")) return(paste("ERROR:", conditionMessage(result)))
  if (inherits(result, "stslist")) return(as.character(TraMineR::seqconc(result)))
  cells <- as.vector(as.matrix(result))
  if (is.character(cells)) return(cells)
  ifelse(is.na(cells), "NA", sprintf("%.17g", cells))
}
twin <- suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(c("pnpn", "pnpp"), ""))),
                                          alphabet = c("p", "n")))
indic <- c("dlgth", "recu", "trans", "transp", "nsubs", "volat", "cplx", "turb", "turbn", "prec", "insec")
calls <- list(
  seqconc = quote(TraMineR::seqconc(obj)),
  seqlength = quote(TraMineR::seqlength(obj)),
  seqdss = quote(TraMineR::seqdss(obj)),
  seqdur = quote(TraMineR::seqdur(obj)),
  seqtransn = quote(TraMineR::seqtransn(obj)),
  seqST = quote(TraMineR::seqST(obj)),
  seqici = quote(TraMineR::seqici(obj)),
  seqivolatility = quote(TraMineR::seqivolatility(obj)),
  seqlogp = quote(TraMineR::seqlogp(obj)),
  seqidegrad = quote(TraMineR::seqidegrad(obj)),
  seqprecarity = quote(TraMineR::seqprecarity(obj)),
  seqinsecurity = quote(TraMineR::seqinsecurity(obj)),
  seqsubsn = quote(TraMineR::seqsubsn(obj)),
  seqformat_SPS = quote(TraMineR::seqformat(obj, from = "STS", to = "SPS", compress = TRUE)),
  seqformat_DSS = quote(TraMineR::seqformat(obj, from = "STS", to = "DSS", compress = TRUE))
)
rows <- list()
for (object in c("recoded", "twin")) {
  obj <- if (object == "recoded") recoded else twin
  for (name in names(calls)) {
    result <- values(calls[[name]])
    rows[[length(rows) + 1L]] <- data.frame(object = object, function_name = name, position = seq_along(result),
                                            value = result, stringsAsFactors = FALSE)
  }
  table <- suppressMessages(TraMineR::seqindic(obj, indic = indic))
  for (column in colnames(table)) {
    rows[[length(rows) + 1L]] <- data.frame(object = object, function_name = paste0("seqindic:", column),
                                            position = seq_len(nrow(table)),
                                            value = sprintf("%.17g", table[[column]]), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "internal_void_indicators.csv"), row.names = FALSE, quote = TRUE)
