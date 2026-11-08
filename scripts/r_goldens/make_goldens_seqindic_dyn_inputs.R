# Live-R goldens for the indic values seqindic.dyn accepts (repair of audit
# finding TX-AUD-20261010-051).
#
# seqindic.dyn stops for an indicator group name, uses only the first value of
# a vector indic with a message (seqindic-dyn.R:10-22), and fails inside R for
# "visitp". Its result is a matrix named after the window-end time labels and
# the sequence row names (seqindic-dyn.R:96-97).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqindic_dyn_inputs.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strings <- c("AABBCCAABB", "ABCABCABCA", "AAAAABBBBB", "CCBBAACCBB", "ABBBCCCAAA", "BBBBBBCCCC")
seqs <- suppressMessages(TraMineR::seqdef(as.data.frame(do.call(rbind, strsplit(strings, ""))),
                                          alphabet = c("A", "B", "C"), id = paste0("s", 1:6)))
attempt <- function(indic) {
  warned <- character()
  value <- withCallingHandlers(
    tryCatch(TraMineRextras::seqindic.dyn(seqs, indic = indic), error = function(e) e),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") },
    # TraMineR's msg.warn reports through message(), not warning().
    message = function(m) { warned <<- c(warned, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })
  list(value = value, warned = warned)
}
rows <- list()
add <- function(case, kind, text) rows[[length(rows) + 1L]] <<- data.frame(case = case, kind = kind, text = text,
                                                                            stringsAsFactors = FALSE)
for (indic in c("basic", "all", "complexity", "visitp")) {
  result <- attempt(indic)
  add(indic, "error", if (inherits(result$value, "error")) conditionMessage(result$value) else "<none>")
}
single <- attempt("cplx")$value
vector <- attempt(c("cplx", "trans"))
add("cplx,trans", "warning", paste(vector$warned, collapse = " | "))
add("cplx,trans", "equals_first", as.character(isTRUE(all.equal(unclass(vector$value), unclass(single)))))
add("cplx", "colnames", paste(colnames(single), collapse = "|"))
add("cplx", "rownames", paste(rownames(single), collapse = "|"))
write.csv(do.call(rbind, rows), file.path(out_dir, "seqindic_dyn_inputs.csv"), row.names = FALSE, quote = TRUE)
