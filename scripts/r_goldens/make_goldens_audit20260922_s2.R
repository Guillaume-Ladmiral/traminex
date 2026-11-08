# Golden generation: construction and format findings of the 2026-09-22 audit
# (slice S2).
#
# One record per line: case, field, value (tab-separated). A state matrix is
# written row by row, rows separated by ";" and cells by "|"; a vector is
# written with its elements separated by "|".
#
# - seqdef_na_defaults / tse_na_transition (S2-05): seqdef's default missing
#   rules (left = NA, gaps = NA, right = "DEL") turn leading and internal NA
#   cells into the nonresponse marker and trailing ones into void; the NA
#   cells TSE_to_STS leaves after an NA transition are decided the same way.
# - fce_tied_events / fce_character_times / tse_tied_events /
#   toPersonPeriod_id_order (S2-07):
#   order() on a character key sorts by the session collation (ICU under a
#   UTF-8 locale, code points under LC_COLLATE=C), so tied events and
#   person ids follow collation, not code points. The collation is pinned to
#   C.UTF-8, the locale traminex's r_character_order reproduces.
# - seqMD_case_order / seqMD_case_cost (S2-08): the combined alphabet of
#   seqMD is seqdef's default, sorted in the same collation.
# - seqstart_literal_star (S2-09): seqstart fills with NA, distinct from a
#   literal "*" state.
# - seqMD_width_fill / seqMD_width_fill_reversed (S2-11): a shorter channel is
#   padded with void, which fill.with.miss does not turn into nonresponse.
#
# Consumed by tests/test_long_builder_missing_cells.py (S2-05),
# tests/test_event_tie_collation.py (S2-07), tests/test_seqmd_combined_states.py
# (S2-08, S2-11) and tests/test_marker_text_state_collision.py (S2-09).
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
stopifnot(Sys.getlocale("LC_COLLATE") == "C.UTF-8", capabilities("ICU"))

lines <- character()
record <- function(case, field, value) {
  value <- as.character(value)
  if (length(value) == 0) value <- "none"
  lines <<- c(lines, paste(case, field, paste(value, collapse = "|"), sep = "\t"))
}
record_cells <- function(case, s) {
  m <- as.matrix(s)
  record(case, "cells", paste(apply(m, 1, paste, collapse = "|"), collapse = ";"))
}
record_seq <- function(case, s) {
  record(case, "alphabet", alphabet(s))
  record_cells(case, s)
}

# S2-05: NA cells under seqdef's defaults (leading, internal, trailing, whole row).
wide <- data.frame(
  a = c("A", NA, "A", NA),
  b = c(NA, "A", "A", NA),
  c = c("B", "B", NA, NA)
)
record_seq("seqdef_na_defaults", suppressMessages(seqdef(wide)))

# S2-05: TSE_to_STS leaves NA after an NA transition; seqdef decides it.
stm <- matrix(c("S1", NA, "S2", "S2"), nrow = 2,
  dimnames = list(c("None", "S1"), c("e1", "e2")))
tse <- data.frame(id = c(1, 1, 2), time = c(1, 2, 1), event = c("e1", "e1", "e2"))
record_seq("tse_na_transition", suppressMessages(seqdef(TSE_to_STS(tse, stm = stm, tmax = 4))))

# S2-07: events tied in time, ordered by FCE_to_TSE.R:40 order(timessel, eventsel).
fce <- data.frame(pid = c(1, 2), Cev = c(3, 1), bev = c(3, 2), aev = c(5, 1), `_ev` = c(3, 1),
  check.names = FALSE)
fce_tse <- FCE_to_TSE(fce, id = "pid", cols = c("Cev", "bev", "aev", "_ev"))
record("fce_tied_events", "id", fce_tse$id)
record("fce_tied_events", "time", fce_tse$time)
record("fce_tied_events", "event", fce_tse$event)
# Character times are ordered by the same collation.
fce_chr <- data.frame(pid = 1, e1 = "b", e2 = "B", e3 = "a")
fce_chr_tse <- FCE_to_TSE(fce_chr, id = "pid", cols = c("e1", "e2", "e3"))
record("fce_character_times", "time", fce_chr_tse$time)
record("fce_character_times", "event", fce_chr_tse$event)

# S2-07: TSE_TO_STS.R:35 order(mtime, mevent) applies tied events in collation
# order; the explicit stm makes the order visible in the final state.
tie_states <- c("None", "S_b", "S_C", "S_bC", "S_b_then_C")
tie_stm <- matrix(c("S_b", "S_bC", "S_b", "S_b_then_C", "S_bC", "S_C", "S_bC", "S_bC",
  "S_b_then_C", "S_b_then_C"), ncol = 2, byrow = TRUE, dimnames = list(tie_states, c("b", "C")))
tie_tse <- data.frame(id = c(1, 1), time = c(2, 2), event = c("C", "b"))
record_cells("tse_tied_events", TSE_to_STS(tie_tse, stm = tie_stm, tmax = 4))

# S2-07: toPersonPeriod.R:15 order(pp$id, pp$timestamp, pp$state) on
# character row names.
pp_seq <- suppressMessages(seqdef(data.frame(a = c("A", "B", "A", "B"), b = c("B", "B", "A", "A")),
  id = c("b", "C", "a", "_z")))
record("toPersonPeriod_id_order", "id", as.character(toPersonPeriod(pp_seq)$id))

# S2-08: seqMD.R:197 builds the combined object with seqdef's default alphabet,
# i.e. seqstatl's sort() in the session collation; the cost matrix rows and
# columns follow it (written column-major).
md1 <- suppressMessages(seqdef(data.frame(a = c("a", "B", "c"), b = c("B", "c", "a")),
  alphabet = c("a", "B", "c")))
md2 <- suppressMessages(seqdef(data.frame(a = c("x", "x", "Y"), b = c("Y", "x", "x")),
  alphabet = c("x", "Y")))
record_seq("seqMD_case_order", suppressMessages(seqMD(list(md1, md2))))
md_cost <- suppressMessages(seqMD(list(md1, md2), sm = "CONSTANT", cweight = c(1, 3),
  what = "cost"))
record("seqMD_case_cost", "alphabet", attr(md_cost, "alphabet"))
record("seqMD_case_cost", "sm", sprintf("%.17g", as.vector(unclass(md_cost)[, ])))

# S2-09: seqstart fills with NA (seqstart.R: missing = NA), a true missing
# value that stays distinct from a literal "*" state (NA written as "NA").
star <- matrix(c("*", "*", "A", "A", "A", "A", "A", "*"), nrow = 2, byrow = TRUE)
star_start <- suppressMessages(seqstart(star, data.start = c(1, 2), new.start = 1))
star_start[is.na(star_start)] <- "NA"
record_cells("seqstart_literal_star", star_start)

# S2-11: beyond a shorter channel's width seqMD.R:169-172 pads with void, and
# fill.with.miss (lines 175-178) touches only the channel's own columns.
w1 <- suppressMessages(seqdef(data.frame(a = c("A", "B"), b = c(NA, "A"), c = c("B", "B")),
  alphabet = c("A", "B"), right = NA))
w2 <- suppressMessages(seqdef(data.frame(a = c("x", "y"), b = c("y", "y"), c = c("x", "y"),
  d = c("y", "x"), e = c("x", NA)), alphabet = c("x", "y")))
record_seq("seqMD_width_fill", suppressWarnings(suppressMessages(seqMD(list(w1, w2)))))
record_seq("seqMD_width_fill_reversed", suppressWarnings(suppressMessages(seqMD(list(w2, w1)))))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tfield\tvalue", lines), file.path(out_dir, "audit20260922_s2_cases.txt"))
cat("AUDIT20260922_S2_GOLDENS_DONE\n")
