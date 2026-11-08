# Live-R goldens for seqdiff on ragged / missing data (RP3/RP4 audit F2).
#
# Exercises seqdiff's with.missing path on a 12x6 state-sequence fixture that mixes
# interior missing (nr) with trailing void (ragged). With with.missing=TRUE, R keeps
# nr as a state and DELETES void before the within-window LCS; traminex now matches
# (the void code is passed through to the distance layer). The earlier seqdiff golden
# used a no-missing fixture, so the void path was unexercised.
#
# Requires R 4.6.2 + TraMineR 2.2-14.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_phase2_inputs.R")  # write_phase2_csv()
suppressMessages(library(TraMineR))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqdiff void goldens require TraMineR 2.2-14")
}

# 12x6, alphabet 1:4. Two interior nr ("*" not at the right edge) and trailing void
# ("*" runs at the right edge -> void via right="DEL") on the second half of the rows.
M <- rbind(
  c("1", "2", "3", "4", "1", "2"),
  c("2", "1", "*", "3", "4", "1"),  # interior nr
  c("3", "3", "1", "2", "2", "3"),
  c("4", "1", "2", "*", "2", "1"),  # interior nr
  c("1", "1", "2", "3", "4", "2"),
  c("2", "3", "4", "1", "1", "3"),
  c("3", "2", "1", "4", "3", "*"),  # trailing void (1)
  c("4", "4", "1", "1", "1", "*"),  # trailing void (1)
  c("1", "3", "2", "4", "3", "*"),  # trailing void (1)
  c("2", "2", "3", "3", "4", "*"),  # trailing void (1)
  c("3", "1", "2", "2", "1", "*"),  # trailing void (1)
  c("4", "2", "1", "3", "2", "*")   # trailing void (1)
)
group <- rep(c(1, 2), length.out = nrow(M))
w <- 0.5 + 0.1 * (1:nrow(M))
sd <- suppressMessages(TraMineR::seqdef(M, alphabet = 1:4, missing = "*", weights = w))

gen <- function(tag, weighted) {
  r <- suppressMessages(TraMineR::seqdiff(
    sd, group, seqdist.args = list(method = "LCS", norm = "none"),
    with.missing = TRUE, weighted = weighted
  ))
  write_phase2_csv(
    paste0("seqdiff_void_", tag, "_stat"), r$stat,
    "# seqdiff(with.missing=TRUE) stat: PseudoF,PseudoFbf,PseudoR2,Bartlett,Levene (rows=positions)\n"
  )
  write_phase2_csv(
    paste0("seqdiff_void_", tag, "_disc"), r$discrepancy,
    "# seqdiff(with.missing=TRUE) discrepancy: group levels + Total (rows=positions)\n"
  )
}

gen("uw", FALSE)
gen("w", TRUE)

cat("SEQDIFF_VOID_GOLDENS_DONE\n")
