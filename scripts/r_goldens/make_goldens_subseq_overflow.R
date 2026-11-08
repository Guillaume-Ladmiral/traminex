# Live-R goldens for NMS, NMSMST and SVRspell where R's double sums overflow (TX-AUD-20261019-029,
# TX-AUD-20261011-047).
#
# R computes distances between distinct sequences only, so identical sequences are at distance 0. Between
# different sequences it forms Ival + Jval - 2 Aval from weighted subsequence counts (NMSdistance.cpp) in doubles,
# and once a count or a duration power overflows the result is NaN (Inf - Inf) or Inf although the distance
# itself is finite. Each case records R's whole matrix, column by column: one-state sequences of lengths 515,
# 515 and 514; kweights of 1e305 on one-state sequences of lengths 12, 12 and 11 beside two short sequences; and
# single spells of lengths 50, 50, 40 and 50 beside two short sequences, with tpow 250 (NMSMST, SVRspell) and 150
# (SVRspell).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_subseq_overflow.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pad <- function(v, width) c(v, rep(NA, width - length(v)))
build <- function(rows) {
  width <- max(lengths(rows))
  suppressMessages(seqdef(do.call(rbind, lapply(rows, pad, width)), alphabet = c("A", "B"), right = "DEL"))
}
onestate <- build(list(rep("A", 515), rep("A", 515), rep("A", 514)))
weighted <- build(list(rep("A", 12), rep("A", 12), rep("A", 11), "B", c("B", "A")))
spells <- build(list(rep("A", 50), rep("A", 50), rep("A", 40), rep("B", 50), c("A", "B"), c("B", "A")))
cases <- list(
  nms_onestate = quote(seqdist(onestate, method = "NMS")),
  nms_kweights = quote(seqdist(weighted, method = "NMS", kweights = rep(1e305, 12))),
  nmsmst_tpow_250 = quote(seqdist(spells, method = "NMSMST", tpow = 250)),
  svrspell_tpow_150 = quote(seqdist(spells, method = "SVRspell", tpow = 150)),
  svrspell_tpow_250 = quote(seqdist(spells, method = "SVRspell", tpow = 250))
)
rows <- list()
for (name in names(cases)) {
  d <- suppressMessages(suppressWarnings(eval(cases[[name]])))
  rows[[length(rows) + 1L]] <- data.frame(case = name, value = paste(format(as.vector(d), digits = 17), collapse = " "))
}
write.csv(do.call(rbind, rows), file.path(out_dir, "subseq_overflow.csv"), row.names = FALSE, quote = TRUE)
