# Live-R goldens for the colours and labels seqrecode gives the new states
# (repair of audit finding TX-AUD-20261011-024).
#
# seqrecode.R:88-104: a recoded state takes the colour of its first source
# state, the otherwise state the colour of the first original state no recode
# uses, a retained state keeps its own; a source that is the nr code gives
# missing.color. labels and cpal are plain vectors, one entry per new state.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqrecode_palette.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seqs <- suppressMessages(TraMineR::seqdef(
  data.frame(a = c("A", "A", "B", "C"), b = c("A", "B", "B", "C"), c = c("B", "B", "C", "C"),
             d = c("B", NA, NA, "D"), e = c("C", "C", NA, "D"), f = c("C", "D", NA, "D")),
  alphabet = c("A", "B", "C", "D")))
cases <- list(
  merge_ab = list(recodes = list(X = c("A", "B"))),
  merge_ab_otherwise = list(recodes = list(X = c("A", "B")), otherwise = "Z"),
  two_recodes = list(recodes = list(X = "C", Y = "D")),
  missing_source = list(recodes = list(M = c("*"), X = c("A", "B"))),
  vector_labels = list(recodes = list(X = c("A", "B")), labels = c("x", "c", "d")),
  vector_cpal = list(recodes = list(X = c("A", "B")), cpal = c("red", "green", "blue")),
  # A recode named like the nr code: R drops it from the new alphabet but still
  # picks each source by position (recodes[[i]]), so a first such recode shifts
  # the others' colours (seqrecode.R:58-62,90-92).
  nr_named_last = list(recodes = list(X = c("A", "B"), "*" = "C")),
  nr_named_first = list(recodes = list("*" = "C", X = c("A", "B"))),
  # A recode whose source is the void code takes white, with R's warning.
  void_source = list(recodes = list(X = c("A", "B"), Y = "%"))
)
rows <- list()
warnings_seen <- list()
for (name in names(cases)) {
  args <- c(list(seqs), cases[[name]])
  seen <- character()
  recoded <- withCallingHandlers(suppressMessages(do.call(TraMineR::seqrecode, args)),
    warning = function(w) { seen <<- c(seen, conditionMessage(w)); invokeRestart("muffleWarning") })
  warnings_seen[[name]] <- paste(seen, collapse = " | ")
  alph <- TraMineR::alphabet(recoded)
  rows[[length(rows) + 1L]] <- data.frame(case = name, position = seq_along(alph), state = alph,
                                          colour = attr(recoded, "cpal"), label = attr(recoded, "labels"),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "seqrecode_palette.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(case = names(warnings_seen), warning = unlist(warnings_seen), stringsAsFactors = FALSE),
          file.path(out_dir, "seqrecode_palette_warnings.csv"), row.names = FALSE, quote = TRUE)
