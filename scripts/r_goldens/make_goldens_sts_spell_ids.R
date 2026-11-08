# Live-R goldens for the ids of the SPELL form (repair of audit finding
# TX-AUD-20261010-017).
#
# R's STS_to_SPELL labels each spell with the sequence object's row names
# (rownames(seqdata, do.NULL = FALSE, prefix = ""), seqformat-STS_to_SPELL.R:31),
# which seqdef(id = ) sets. STS_to_SPELL is defined but not exported, so it is
# called through ::: as well as through the public seqformat(to = "SPELL").
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_sts_spell_ids.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  p1_p2 = list(data = data.frame(a = c("A", "A"), b = c("A", "B"), c = c("B", "B")), id = c("p1", "p2")),
  z_x_y = list(data = data.frame(a = c("A", "B", "B"), b = c("B", "B", "A")), id = c("z", "x", "y")),
  # Trailing missing cells kept as nonresponse (right = NA; TX-AUD-20261019-014): the internal
  # STS_to_SPELL walks each row to seqlength(), which counts them, so it keeps the spell
  # (q2, 2-4, *); seqformat re-runs seqdef(right = "DEL") first and drops it.
  trailing_nr = list(data = data.frame(a = c("A", "B"), b = c("A", NA), c = c(NA, NA), d = c("B", NA)),
                     id = c("q1", "q2"), right = NA)
)
rows <- list()
for (name in names(cases)) {
  right <- if (is.null(cases[[name]]$right)) "DEL" else cases[[name]]$right
  s <- suppressMessages(TraMineR::seqdef(cases[[name]]$data, id = cases[[name]]$id, right = right))
  routes <- c("STS_to_SPELL", "seqformat", if (!is.null(cases[[name]]$right)) "STS_to_SPELL_without_missing")
  for (route in routes) {
    spell <- if (route == "STS_to_SPELL") TraMineR:::STS_to_SPELL(s) else
      if (route == "STS_to_SPELL_without_missing") TraMineR:::STS_to_SPELL(s, with.missing = FALSE) else
      suppressMessages(TraMineR::seqformat(s, from = "STS", to = "SPELL"))
    rows[[length(rows) + 1L]] <- data.frame(case = name, route = route, position = seq_len(nrow(spell)),
                                            id = as.character(spell$id), begin = spell$begin, end = spell$end,
                                            states = as.character(spell$states), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "sts_spell_ids.csv"), row.names = FALSE, quote = TRUE)
