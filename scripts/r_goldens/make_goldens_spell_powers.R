# Live-R goldens for R's power of spell lengths (TX-AUD-20261019-030).
#
# seqdist raises spell lengths to a power with R's `^`: OMslen's `dur.mat ^ (-1 * h)` (seqdist.R:619) and the
# `seqdur(...) ^ tpow` of OMspell, NMSMST and SVRspell (seqdist.R:626). The "power" rows hold R's value of d^p for
# lengths 1-300 and seven exponents; the "dist" rows hold OMslen, OMspell, NMSMST and SVRspell distances on twenty
# run-heavy sequences, column by column, and NMS's on the same data (the subsequence family's weighted sums,
# NMSdistance.cpp, TX-REP-20261022-001). Every value is written in hexadecimal (sprintf "%a"), so a test can
# compare bits.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_spell_powers.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- c(
  "AABCCDDDDDBBBAAAAAACCCCAABBBCC", "CBBBDDDDDDCCCCCCCCCCCCCCDDDDBB", "ADDDBBBBBBBAAAAAAACCCCCDDDBBBB",
  "AABBBBBBBBBBCCCCCCBDBBBDBBBBBB", "CCCCBAAABBBBBBAAAABBBBBBBBBBAA", "BBBBBAAAAAABBBBBBBCAAAABBBBBBA",
  "AAAAAAAAAAABBAAAAAAAAAAAAAAAAA", "AACCCCCAAAAAAAABDDDDDDDDAAAAAB", "AAAAAAAAAAADDDDDAAAAAACCCAAAAA",
  "BBBBCCCCCCDDDDAAAAAAABBBBBBBBB", "BBBBBBBBCCCCCCCDDCAAADDDDDDBBB", "BBDDDDDDAAAAAACCCCCCCCCCCCCBBD",
  "CCCCCAAAAAAAAAAAADDDDDDDDDDDDD", "DAAAAAAAAAAAAAAAAAAAAAAAAAAABD", "AAAABAAAAAAABBBBBBBCAAAABBDDDD",
  "AAAAADDDDDAAAAACCCBBBBCCCCCCCC", "DDDDDDBBBBBDDDDDDCCBDDDDDBBBBB", "AAAAADDDDDDDDDDDBADDDDDDDDAAAA",
  "BBBBBBBBBBDDDDCDDDDDDBBBBBBBDD", "AABBBBBBBBAABBAAAAAAAAAAAAAACC"
)
x <- suppressMessages(seqdef(do.call(rbind, strsplit(rows, "")), alphabet = c("A", "B", "C", "D")))
out <- list()
lengths <- as.double(1:300)
for (p in c(-2, -1.3, -0.5, 0.5, 1.3, 2, 2.5)) {
  out[[length(out) + 1L]] <- data.frame(kind = "power", key = format(p), index = 1:300, value = sprintf("%a", lengths ^ p))
}
cases <- list(
  OMslen_h_0.5 = quote(seqdist(x, method = "OMslen", sm = "CONSTANT", indel = 1, h = 0.5)),
  OMslen_h_1.3 = quote(seqdist(x, method = "OMslen", sm = "CONSTANT", indel = 1, h = 1.3)),
  OMspell_tpow_1.3 = quote(seqdist(x, method = "OMspell", sm = "CONSTANT", indel = 1, expcost = 0.5, tpow = 1.3)),
  NMSMST_tpow_1.3 = quote(seqdist(x, method = "NMSMST", tpow = 1.3)),
  SVRspell_tpow_1.3 = quote(seqdist(x, method = "SVRspell", tpow = 1.3)),
  NMS = quote(seqdist(x, method = "NMS"))
)
for (name in names(cases)) {
  d <- suppressMessages(eval(cases[[name]]))
  out[[length(out) + 1L]] <- data.frame(kind = "dist", key = name, index = seq_along(as.vector(d)), value = sprintf("%a", as.vector(d)))
}
write.csv(do.call(rbind, out), file.path(out_dir, "spell_powers.csv"), row.names = FALSE, quote = TRUE)
