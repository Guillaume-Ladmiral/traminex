# Live-R goldens for order() over a sortv vector (repair of audit finding
# TX-AUD-20261010-060).
#
# plot.stslist sorts the rows with order(sortv) (plot.stslist.R:68-75): text in
# the session's collation (ICU, pinned to the root locale here), numbers by
# value, NA last.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_sortv_character_order.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  text_case_punctuation = c("b", "A", "a", "B", "_x", "c"),
  text_missing_accent = c("b", NA, "a", "é", "f"),
  numeric_missing = c(3, NA, 1),
  # Numbers, infinities and logicals with NA (batch-2 review I-1): order()
  # puts NA last, after Inf.
  numeric_missing_five = c(5, 1, NA, 4, 2),
  numeric_infinite_missing = c(Inf, NA, -Inf, 2),
  logical_missing = c(TRUE, NA, FALSE),
  # A factor whose levels are not in text order (TX-AUD-20261019-019): order() follows the
  # level codes (z before a), NA last.
  factor_levels_z_a = factor(c("a", "z", "a", "z", NA), levels = c("z", "a"))
)
rows <- list()
for (name in names(cases)) {
  values <- cases[[name]]
  rows[[length(rows) + 1L]] <- data.frame(case = name, kind = "value", position = seq_along(values),
                                          value = ifelse(is.na(values), "<NA>", as.character(values)),
                                          stringsAsFactors = FALSE)
  rows[[length(rows) + 1L]] <- data.frame(case = name, kind = "order", position = seq_along(values),
                                          value = as.character(order(values)), stringsAsFactors = FALSE)
  if (is.factor(values)) {
    rows[[length(rows) + 1L]] <- data.frame(case = name, kind = "levels", position = seq_along(levels(values)),
                                            value = levels(values), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "sortv_character_order.csv"), row.names = FALSE, quote = TRUE)
