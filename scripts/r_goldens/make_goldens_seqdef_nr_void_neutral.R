# Live-R goldens for two seqdef inputs traminex refuses (record of audit finding
# TX-AUD-20261010-035).
#
# nr equal to void: R accepts it and then handles missing cells as void, so a
# missing observation disappears from the sequence. "NEUTRAL" left/gaps/right:
# R keeps a plain NA cell, neither the nr nor the void code (seqdef.R).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqdef_nr_void_neutral.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- data.frame(a = c("A", NA), b = c("B", "A"), c = c(NA, "B"))
gap <- data.frame(a = "A", b = NA, c = "B")
cases <- list(
  nr_equals_void = list(data = d, args = list(nr = "%")),
  left_neutral = list(data = d, args = list(left = "NEUTRAL")),
  right_neutral = list(data = d, args = list(right = "NEUTRAL")),
  gaps_neutral = list(data = gap, args = list(gaps = "NEUTRAL"))
)
rows <- list()
for (name in names(cases)) {
  s <- suppressMessages(do.call(TraMineR::seqdef, c(list(cases[[name]]$data), cases[[name]]$args)))
  rows[[length(rows) + 1L]] <- data.frame(case = name, sequence = seq_len(nrow(s)),
                                          seqconc = as.character(TraMineR::seqconc(s)),
                                          seqlength = as.integer(TraMineR::seqlength(s)),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqdef_nr_void_neutral.csv"), row.names = FALSE, quote = TRUE)
