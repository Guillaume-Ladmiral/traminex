# Live-R goldens for seqpolyads' p-value of U when U is negative (record of
# audit finding TX-AUD-20261010-052).
#
# R computes U.tp <- 2*pt(mean.U/se, T-1, lower.tail=F) (seqpolyads.R:184), the
# two-sided p-value only when U is positive: a negative U gives U.tp above 1.
# The dyads below are more distant than random pairs, so every U is negative.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqpolyads_negative_u.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

first <- suppressMessages(TraMineR::seqdef(data.frame(a = rep(c("A", "B"), 3), b = rep(c("A", "B"), 3),
                                                      c = rep(c("A", "B"), 3)), alphabet = c("A", "B")))
second <- suppressMessages(TraMineR::seqdef(data.frame(a = rep(c("B", "A"), 3), b = rep(c("B", "A"), 3),
                                                       c = rep(c("B", "A"), 3)), alphabet = c("A", "B")))
result <- suppressMessages(TraMineRextras::seqpolyads(list(first, second), T = 200, s = 7))
out <- data.frame(dyad = seq_along(result$U), U = sprintf("%.17g", as.vector(result$U)),
                  U.tp = sprintf("%.17g", as.vector(result$U.tp)), stringsAsFactors = FALSE)
write.csv(out, file.path(out_dir, "seqpolyads_negative_u.csv"), row.names = FALSE, quote = TRUE)
