# Live-R goldens for seqplot's group panels (repair of audit finding
# TX-AUD-20261010-044).
#
# seqplot draws one panel per level of its groups (seqplot.R:163-205) and
# titles them with the levels when main = "auto", with "<main> - <level>" for
# one main, and with main itself when it has one title per level. It builds the
# groups with its internal group() (TraMineR-group.R:15-22), which pastes a
# space before each level and re-factors, so every title starts with a space
# and numeric groups are ordered as text. traminex corrects that quirk as it does
# for seqHtplot (TX-REP-20261017-002): the levels are factor()'s. The column
# r_title is what R's seqplot draws (title() calls, traced); factor_title is R's
# title rule applied to levels(factor(group)).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqplot_group_titles.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
pdf(NULL)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seqs <- suppressMessages(TraMineR::seqdef(
  data.frame(t(sapply(c("AABBC", "ABBBC", "BBAAC", "AABBB", "BBBAA"), function(s) strsplit(s, "")[[1]]))),
  alphabet = c("A", "B", "C")))
text_group <- c("b", "a", "b", "a", "b")
cases <- list(
  list(case = "d_auto", type = "d", group = text_group, main = "auto"),
  list(case = "i_auto", type = "i", group = text_group, main = "auto"),
  list(case = "f_auto", type = "f", group = text_group, main = "auto"),
  list(case = "mt_auto", type = "mt", group = text_group, main = "auto"),
  list(case = "ms_auto", type = "ms", group = text_group, main = "auto"),
  list(case = "d_one_main", type = "d", group = text_group, main = "X"),
  list(case = "d_main_per_level", type = "d", group = text_group, main = c("P", "Q")),
  list(case = "d_numeric", type = "d", group = c(10, 2, 10, 2, 10), main = "auto"),
  list(case = "d_missing", type = "d", group = c("b", NA, "b", "a", "b"), main = "auto"),
  list(case = "d_one_sequence", type = "d", group = c("b", "a", "b", "b", "b"), main = "auto")
)
titles <- character(0)
trace("title", tracer = quote(if (!is.null(main)) titles <<- c(titles, main)), print = FALSE,
      where = asNamespace("TraMineR"))
rows <- list()
for (cs in cases) {
  titles <- character(0)
  TraMineR::seqplot(seqs, type = cs$type, group = cs$group, main = cs$main)
  levs <- levels(factor(cs$group))
  rule <- if (identical(cs$main, "auto")) levs else if (length(cs$main) == 1) paste(cs$main, "-", levs) else cs$main
  stopifnot(length(titles) == length(rule))
  rows[[length(rows) + 1L]] <- data.frame(case = cs$case, type = cs$type, position = seq_along(titles),
                                          r_title = titles, factor_level = levs, factor_title = rule,
                                          stringsAsFactors = FALSE)
}
untrace("title", where = asNamespace("TraMineR"))
write.csv(do.call(rbind, rows), file.path(out_dir, "seqplot_group_titles.csv"), row.names = FALSE, quote = TRUE)
