# Live-R goldens for the alphabet that HSPELL_to_STS data gets (repair of audit
# finding TX-AUD-20261010-031).
#
# TraMineRextras::HSPELL_to_STS returns a character data frame; seqdef then
# sorts its alphabet as seqstatl does: numbers by value when every state is an R
# number, R's collation otherwise (seqdef.R, seqstatl.R:21-30). ICU is pinned to
# the root locale.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_hspell_alphabet_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  numeric = c("1", "10", "11", "12", "2", "3", "9"),
  text = c("_x", "a", "A", "b", "B", "c", "C"),
  accents = c("e", "é", "f", "E", "ê", "d", "D")
)
rows <- list()
for (name in names(cases)) {
  spells <- data.frame(b1 = rep(1, 7), e1 = rep(2, 7), s1 = cases[[name]], stringsAsFactors = FALSE)
  states <- suppressMessages(TraMineRextras::HSPELL_to_STS(spells, begin = "b1", end = "e1", status = "s1"))
  alph <- TraMineR::alphabet(suppressMessages(TraMineR::seqdef(states)))
  rows[[length(rows) + 1L]] <- data.frame(case = name, position = seq_along(alph), state = alph,
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "hspell_alphabet_order.csv"), row.names = FALSE, quote = TRUE)
