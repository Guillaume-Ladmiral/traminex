# Golden generation: the event and group-comparison surfaces on a WIDE corpus.
#
# Cycle-7 W-04. Cycle 6 adjudicated these surfaces against 5 sequences / 13
# events and two 8x6 matrices, labelled SAMPLE. That was enough to find three
# defects; it is not enough to conclude the rest is clean. This generator uses
# a corpus an order of magnitude larger and deliberately more awkward:
#
#   * 38 event sequences instead of 5, of uneven length (1 to 14 events) --
#     two of the 40 draw no events and so never reach the TSE table at all,
#     which is itself the shape a real extract produces;
#   * 8 event types instead of 3;
#   * non-unit, non-uniform weights;
#   * simultaneous events (repeated timestamps) in a third of the sequences;
#   * an unbalanced two-group split and an unbalanced three-group
#     split, so seqecmpgroup meets a table that is not square;
#   * a separate 30 x 10 state-sequence corpus with missing data for the
#     group-comparison side, again unbalanced and weighted.
#
# Every mined quantity is recorded: the subsequence labels, Support and Count,
# the full seqeapplysub matrix under all five count methods, the seqerules
# statistics, the seqecmpgroup table, seqentrans, and seqedist over every
# interval x norm pair.
#
# seqerules is one of the unexported internals traminex promotes to a public
# name (docs/API_PARITY_LEDGER.md); it has no exported counterpart, so
# TraMineR:::seqerules is called deliberately here -- that is the only way to
# obtain its values from the oracle.
#
# Values are written with %.15g. Requires live R + TraMineR 2.2-14.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(packageVersion("TraMineR") == "2.2.14",
          packageVersion("TraMineRextras") == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fmt <- function(x) sprintf("%.15g", x)

## --------------------------------------------------------------------- ##
## 1. The wide event corpus.                                             ##
## --------------------------------------------------------------------- ##
set.seed(20260830)
n_seq <- 40L
event_types <- paste0("E", 1:8)

rows <- NULL
for (i in seq_len(n_seq)) {
  n_ev <- if (i %in% c(7L, 23L)) 0L else sample(1:14, 1L)
  if (n_ev == 0L) next
  times <- sort(sample(1:20, n_ev, replace = TRUE))   # replace=TRUE -> ties
  if (i %% 3L != 0L) times <- sort(unique(times))      # only every third keeps ties
  n_ev <- length(times)
  rows <- rbind(rows, data.frame(
    id = sprintf("s%02d", i),
    timestamp = as.numeric(times),
    event = sample(event_types, n_ev, replace = TRUE),
    stringsAsFactors = FALSE
  ))
}
write.csv(rows, file.path(out_dir, "events_wide_tse.csv"), row.names = FALSE)

eseq <- seqecreate(id = rows$id, timestamp = rows$timestamp, event = rows$event)
present_ids <- unique(rows$id)
weights <- round(seq(0.3, 2.7, length.out = length(present_ids)), 3)
seqeweight(eseq) <- weights
writeLines(fmt(weights), file.path(out_dir, "events_wide_weights.txt"))
writeLines(levels(eseq), file.path(out_dir, "events_wide_levels.txt"))
writeLines(as.character(eseq), file.path(out_dir, "events_wide_aschar.txt"))

# Unbalanced groups over the sequences that actually exist.
n_present <- length(present_ids)
group2 <- ifelse(seq_len(n_present) %% 3L == 0L, "late", "early")
group3 <- c("a", "b", "c")[1L + (seq_len(n_present) * 7L) %% 3L]
writeLines(group2, file.path(out_dir, "events_wide_group2.txt"))
writeLines(group3, file.path(out_dir, "events_wide_group3.txt"))

store <- list()
put <- function(key, value) store[[key]] <<- value

## Mining, over several support thresholds and count methods.
for (ms in c(4, 8, 15)) {
  mined <- seqefsub(eseq, min.support = ms)
  put(sprintf("fsub|ms=%g|subseq", ms), paste(as.character(mined$subseq), collapse = ";"))
  put(sprintf("fsub|ms=%g|support", ms), paste(fmt(mined$data$Support), collapse = ";"))
  put(sprintf("fsub|ms=%g|count", ms), paste(fmt(mined$data$Count), collapse = ";"))
  put(sprintf("fsub|ms=%g|n", ms), as.character(length(mined$subseq)))
}
for (ms in c(4, 8)) {
  mined <- seqefsub(eseq, min.support = ms, weighted = FALSE)
  put(sprintf("fsub_unweighted|ms=%g|subseq", ms),
      paste(as.character(mined$subseq), collapse = ";"))
  put(sprintf("fsub_unweighted|ms=%g|count", ms), paste(fmt(mined$data$Count), collapse = ";"))
}
for (mk in c(1, 2, 3)) {
  mined <- seqefsub(eseq, min.support = 6, max.k = mk)
  put(sprintf("fsub|maxk=%d|subseq", mk), paste(as.character(mined$subseq), collapse = ";"))
  put(sprintf("fsub|maxk=%d|count", mk), paste(fmt(mined$data$Count), collapse = ";"))
}
for (cm in c("COBJ", "CDIST_O", "CMINWIN", "CDIST")) {
  cst <- seqeconstraint(count.method = cm)
  mined <- seqefsub(eseq, min.support = 6, constraint = cst)
  put(sprintf("fsub|cm=%s|subseq", cm), paste(as.character(mined$subseq), collapse = ";"))
  put(sprintf("fsub|cm=%s|count", cm), paste(fmt(mined$data$Count), collapse = ";"))
}
cst_win <- seqeconstraint(count.method = "CWIN", window.size = 6)
mined_win <- seqefsub(eseq, min.support = 6, constraint = cst_win)
put("fsub|cm=CWIN|subseq", paste(as.character(mined_win$subseq), collapse = ";"))
put("fsub|cm=CWIN|count", paste(fmt(mined_win$data$Count), collapse = ";"))

# Constrained mining: max.gap, window.size, age bounds.
for (spec in list(
  list(tag = "maxgap3", cst = seqeconstraint(max.gap = 3)),
  list(tag = "window8", cst = seqeconstraint(window.size = 8)),
  list(tag = "agemax10", cst = seqeconstraint(age.max = 10)),
  list(tag = "agemin5", cst = seqeconstraint(age.min = 5)),
  list(tag = "gap2_win5", cst = seqeconstraint(max.gap = 2, window.size = 5))
)) {
  mined <- seqefsub(eseq, min.support = 5, constraint = spec$cst)
  put(sprintf("fsub|%s|subseq", spec$tag), paste(as.character(mined$subseq), collapse = ";"))
  put(sprintf("fsub|%s|count", spec$tag), paste(fmt(mined$data$Count), collapse = ";"))
}

## The count matrix, under every method the wrapper maps.
base_mined <- seqefsub(eseq, min.support = 8)
put("base|subseq", paste(as.character(base_mined$subseq), collapse = ";"))
for (method in c("count", "presence", "age", "COBJ", "CDIST_O", "CMINWIN", "CDIST")) {
  matrix_value <- seqeapplysub(base_mined, method = method)
  put(sprintf("applysub|%s", method),
      paste(fmt(as.vector(as.matrix(matrix_value))), collapse = ";"))
  put(sprintf("applysub|%s|dim", method),
      paste(dim(matrix_value), collapse = "x"))
}
matrix_win <- seqeapplysub(base_mined, constraint = seqeconstraint(count.method = "CWIN",
                                                                  window.size = 6))
put("applysub|CWIN", paste(fmt(as.vector(as.matrix(matrix_win))), collapse = ";"))

## Association rules (unexported internal; see the header note).
rules <- TraMineR:::seqerules(base_mined)
put("rules|dim", paste(dim(rules), collapse = "x"))
put("rules|Rules", paste(as.character(rules$Rules), collapse = ";"))
for (column in c("Support", "Conf", "Lift", "Standardlift", "JMeasure", "ImplicStat",
                 "p.value", "p.valueB1", "p.valueB2")) {
  put(sprintf("rules|%s", column), paste(fmt(rules[[column]]), collapse = ";"))
}
for (sortv in c("support", "lift", "jmeasure")) {
  sorted_rules <- TraMineR:::seqerules(base_mined, sortv = sortv, decreasing = TRUE)
  put(sprintf("rules_sorted|%s", sortv),
      paste(as.character(sorted_rules$Rules), collapse = ";"))
}

## Discriminating subsequences, two- and three-group, weighted and not.
for (spec in list(
  list(tag = "g2_w", group = group2, weighted = TRUE, method = "chisq"),
  list(tag = "g2_uw", group = group2, weighted = FALSE, method = "chisq"),
  list(tag = "g3_w", group = group3, weighted = TRUE, method = "chisq"),
  list(tag = "g2_bonf", group = group2, weighted = TRUE, method = "bonferroni")
)) {
  compared <- seqecmpgroup(base_mined, group = spec$group, weighted = spec$weighted,
                           method = spec$method)
  put(sprintf("cmpgroup|%s|subseq", spec$tag),
      paste(as.character(compared$subseq), collapse = ";"))
  put(sprintf("cmpgroup|%s|columns", spec$tag),
      paste(colnames(compared$data), collapse = ";"))
  put(sprintf("cmpgroup|%s|values", spec$tag),
      paste(fmt(as.vector(as.matrix(compared$data))), collapse = ";"))
}

## seqentrans, seqecontain, and the accessors.
transitions <- seqentrans(base_mined)
put("entrans|ntrans", paste(fmt(transitions$data$ntrans), collapse = ";"))
put("entrans|nevent", paste(fmt(transitions$data$nevent), collapse = ";"))
put("contain|any_E1E2", paste(as.integer(seqecontain(eseq, c("E1", "E2"))), collapse = ";"))
put("contain|all_E1E2", paste(as.integer(seqecontain(eseq, c("E1", "E2"), unknown.exclude = TRUE)),
                              collapse = ";"))
put("accessor|seqeid", paste(fmt(seqeid(eseq)), collapse = ";"))
put("accessor|seqeweight", paste(fmt(seqeweight(eseq)), collapse = ";"))
put("accessor|seqelength", paste(fmt(seqelength(eseq)), collapse = ";"))

## Event-sequence edit distance over every interval x norm pair.
idcost <- seq_along(levels(eseq)) / 4
writeLines(fmt(idcost), file.path(out_dir, "events_wide_idcost.txt"))
for (interval in c("No", "previous", "next")) {
  for (norm in c("none", "max", "YujianBo")) {
    distance <- TraMineRextras::seqedist(eseq, idcost = idcost, vparam = 0.5,
                                         interval = interval, norm = norm)
    put(sprintf("edist|%s|%s", interval, norm),
        paste(fmt(as.vector(distance)), collapse = ";"))
  }
}

events_table <- data.frame(key = names(store), value = unlist(store, use.names = FALSE),
                           stringsAsFactors = FALSE)
write.csv(events_table, file.path(out_dir, "events_wide_values.csv"), row.names = FALSE)

## --------------------------------------------------------------------- ##
## 2. A wider group-comparison corpus, with missing data.                 ##
## --------------------------------------------------------------------- ##
set.seed(20260831)
gc_n <- 30L
gc_len <- 10L
gc_alphabet <- c("A", "B", "C", "D", "E")
gc_codes <- matrix(sample(seq_along(gc_alphabet), gc_n * gc_len, replace = TRUE),
                   nrow = gc_n, ncol = gc_len)
gc_states <- matrix(gc_alphabet[gc_codes], nrow = gc_n, ncol = gc_len)
# Ragged right edge: the last few positions of a third of the rows are missing.
for (i in seq_len(gc_n)) {
  if (i %% 3L == 0L) gc_states[i, (gc_len - (i %% 4L)):gc_len] <- NA
}
gc_weights <- round(seq(0.25, 3.0, length.out = gc_n), 3)
gc_group <- ifelse(seq_len(gc_n) %% 4L == 0L, "minority", "majority")   # 7 / 23
gc_set <- c(rep("first", 16L), rep("second", 14L))

gc_data <- suppressMessages(
  seqdef(as.data.frame(gc_states, stringsAsFactors = FALSE),
         alphabet = gc_alphabet, weights = gc_weights)
)
write.csv(as.data.frame(ifelse(is.na(gc_states), "<NA>", gc_states)),
          file.path(out_dir, "events_wide_gc_states.csv"), row.names = FALSE)
writeLines(gc_alphabet, file.path(out_dir, "events_wide_gc_alphabet.txt"))
writeLines(fmt(gc_weights), file.path(out_dir, "events_wide_gc_weights.txt"))
writeLines(gc_group, file.path(out_dir, "events_wide_gc_group.txt"))
writeLines(gc_set, file.path(out_dir, "events_wide_gc_set.txt"))

gc_store <- list()
gput <- function(key, value) gc_store[[key]] <<- value

gc_diss <- suppressMessages(seqdist(gc_data, method = "LCS", with.missing = TRUE))
write.csv(as.data.frame(gc_diss), file.path(out_dir, "events_wide_gc_diss.csv"),
          row.names = FALSE)

for (stat in c("all", "LRT", "BIC")) {
  for (wmode in c("FALSE", "TRUE", "by.group")) {
    weighted <- if (wmode == "by.group") "by.group" else as.logical(wmode)
    result <- suppressMessages(dissCompare(
      gc_diss, gc_group, s = 0, stat = stat, weighted = weighted,
      weights = gc_weights, inc.df = TRUE
    ))
    gput(sprintf("diss|stat=%s|weighted=%s|columns", stat, wmode),
         paste(colnames(result), collapse = ";"))
    gput(sprintf("diss|stat=%s|weighted=%s|values", stat, wmode),
         paste(fmt(as.vector(as.matrix(result))), collapse = ";"))
  }
}
# Per-set, run on each subsetted matrix so the upstream set= index defect
# (recorded in W-02) is not baked in here.
for (level in c("first", "second")) {
  keep <- which(gc_set == level)
  result <- suppressMessages(dissCompare(
    gc_diss[keep, keep], gc_group[keep], s = 0, stat = "all",
    weights = gc_weights[keep], inc.df = TRUE
  ))
  gput(sprintf("diss_set_alone|%s|values", level),
       paste(fmt(as.vector(as.matrix(result))), collapse = ";"))
}
# HAM needs equal lengths; the ragged right edge here deletes into voids, so
# OM is used as the second method instead.
for (method in c("LCS", "OM")) {
  for (wmode in c("FALSE", "TRUE")) {
    weighted <- as.logical(wmode)
    gput(sprintf("seqCompare|%s|weighted=%s", method, wmode),
         paste(fmt(as.vector(as.matrix(suppressMessages(seqCompare(
           gc_data, group = gc_group, s = 0, stat = "all", method = method,
           weighted = weighted, with.missing = TRUE,
           sm = if (method == "OM") "CONSTANT" else NULL,
           indel = if (method == "OM") 1.0 else "auto"
         ))))), collapse = ";"))
  }
}

gc_table <- data.frame(key = names(gc_store), value = unlist(gc_store, use.names = FALSE),
                       stringsAsFactors = FALSE)
write.csv(gc_table, file.path(out_dir, "events_wide_gc_values.csv"), row.names = FALSE)

cat("EVENTS_WIDE_CORPUS_GOLDENS_DONE\n")
cat("event sequences:", length(eseq), " event rows:", nrow(rows),
    " event keys:", nrow(events_table), " group-comparison keys:", nrow(gc_table), "\n")
