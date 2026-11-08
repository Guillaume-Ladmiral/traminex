# Live-R goldens for seqMD details (TX-AUD-20261010-074).
#
# The combined sequence object: seqMD.R:197 calls seqdef(newseqdata, cnames=md.cnames), so it takes seqdef's
# default markers ("*", "%") and start (1) while copying the channels' column names. The deprecated what= values
# print R's msg.warn text (seqMD.R:17-24). cweight = c(0, 0): with link = "mean" R divides by sum(cweight) = 0
# (seqMD.R:346-348), giving NaN costs for what = "cost" and a refusal for what = "diss"; with link = "sum" the costs
# are 0. Each row records a value, a message or an error.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqmd_details.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

e1 <- suppressMessages(TraMineR::seqdef(rbind(c("A", NA, "B", "B"), c("A", "B", "B", NA)), nr = "?", void = "~",
                                        right = NA, start = 5))
e2 <- suppressMessages(TraMineR::seqdef(rbind(c("X", "X", "Y", "Y"), c("X", "Y", "Y", "Y")), nr = "?", void = "~",
                                        start = 5))
md <- suppressMessages(TraMineR::seqMD(list(e1, e2), what = "MDseq", with.missing = TRUE))
d1 <- suppressMessages(TraMineR::seqdef(rbind(c("A", "A", "B", "B"), c("A", "B", "B", "B"), c("B", "B", "A", "A"))))
d2 <- suppressMessages(TraMineR::seqdef(rbind(c("X", "X", "Y", "Y"), c("X", "Y", "Y", "Y"), c("Y", "Y", "X", "X"))))
rows <- list(
  data.frame(case = "mdseq_nr", value = attr(md, "nr")),
  data.frame(case = "mdseq_void", value = attr(md, "void")),
  data.frame(case = "mdseq_start", value = as.character(attr(md, "start"))),
  data.frame(case = "mdseq_colnames", value = paste(colnames(md), collapse = " "))
)
for (what in c("sm", "seqmc")) {
  msgs <- character(0)
  withCallingHandlers(TraMineR::seqMD(list(d1, d2), sm = list("CONSTANT", "CONSTANT"), what = what),
    message = function(m) { msgs <<- c(msgs, sub("\n$", "", conditionMessage(m))); invokeRestart("muffleMessage") })
  rows[[length(rows) + 1L]] <- data.frame(case = paste0("deprecated_", what), value = grep("deprecated", msgs, value = TRUE)[1])
}
for (lk in c("mean", "sum")) {
  for (what in c("cost", "diss")) {
    value <- tryCatch({
      r <- suppressMessages(TraMineR::seqMD(list(d1, d2), method = if (what == "diss") "OM" else NULL,
                                            sm = list("CONSTANT", "CONSTANT"), cweight = c(0, 0), link = lk, what = what))
      if (what == "cost") paste(c(format(as.vector(r)), "indel", format(attr(r, "indel"))), collapse = " ")
      else paste(format(as.vector(r)), collapse = " ")
    }, error = function(e) paste("ERROR", conditionMessage(e)))
    rows[[length(rows) + 1L]] <- data.frame(case = paste0("cweight0_", lk, "_", what), value = value)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqmd_details.csv"), row.names = FALSE, quote = TRUE)
