# Golden generation: deep and degenerate inputs to TraMineR's event engine.
#
# Cycle-7 W-01. The vendored C++ miner is reached from Python through a
# hand-written port (traminex/events/_engine.py). Ordinary input was already
# golden-pinned; this generator pins the two shapes nobody had fed it:
#
#   1. DEEP -- sequences carrying long runs of SIMULTANEOUS events (zero gaps).
#      That run is what drives SequenceEventNode::getNextWithGap (eventseq.h:100),
#      the one genuinely self-recursive routine the port reproduces, to its
#      deepest nesting. It is reached only from tmrseqedist with interval="next"
#      (seqedist.cpp:171,186). C recurses on the machine stack; the Python port
#      must reach the same answers at depths where CPython's own limit would
#      otherwise stop it.
#
#   2. DEGENERATE -- one event, all events simultaneous, unsorted timestamps,
#      duplicated events, negative timestamps, and boundary constraints.
#
# Cases where R itself aborts (mining that keeps nothing makes TraMineR fail in
# C with "attempt access index 0/0 in VECTOR_ELT") are recorded here as that
# error text, NOT as a value: traminex deliberately returns an empty result
# instead, and the focused test pins the divergence rather than copying it.
#
# Values are written with format(digits = 17) / %.15g -- never read off a print.
# Requires live R + TraMineR 2.2-14 / TraMineRextras 0.6.9.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(packageVersion("TraMineR") == "2.2.14",
          packageVersion("TraMineRextras") == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

record <- function(store, key, expr) {
  value <- tryCatch(
    {
      out <- force(expr)
      if (is.numeric(out)) paste(format(as.vector(out), digits = 17), collapse = ";")
      else if (is.character(out)) paste(out, collapse = ";")
      else paste(as.character(out), collapse = ";")
    },
    error = function(e) paste0("ERROR: ", conditionMessage(e))
  )
  store[[key]] <- value
  store
}

## --------------------------------------------------------------------- ##
## 1. Deep simultaneous runs -> seqedist over every interval x norm pair. ##
## --------------------------------------------------------------------- ##
# Three sequences with DIFFERENT simultaneous blocks so the distances are not
# degenerate: a and b share a 1200-event block but close at different times,
# c is short. RUN_LEN exceeds CPython's default recursion limit of 1000.
RUN_LEN <- 1200L

deep_tse <- rbind(
  data.frame(id = "a", timestamp = 1,
             event = paste0("E", seq_len(RUN_LEN) %% 3L), stringsAsFactors = FALSE),
  data.frame(id = "a", timestamp = 5, event = "F", stringsAsFactors = FALSE),
  data.frame(id = "b", timestamp = 1,
             event = paste0("E", seq_len(RUN_LEN) %% 2L), stringsAsFactors = FALSE),
  data.frame(id = "b", timestamp = 7, event = "F", stringsAsFactors = FALSE),
  data.frame(id = "c", timestamp = c(1, 2, 3),
             event = c("E0", "E1", "F"), stringsAsFactors = FALSE)
)
deep <- seqecreate(id = deep_tse$id, timestamp = deep_tse$timestamp,
                   event = deep_tse$event)
writeLines(levels(deep), file.path(out_dir, "event_edges_deep_levels.txt"))
writeLines(sprintf("%d", RUN_LEN), file.path(out_dir, "event_edges_deep_runlen.txt"))

deep_idcost <- seq_along(levels(deep)) / 2  # distinct, non-unit indel costs
writeLines(sprintf("%.15g", deep_idcost),
           file.path(out_dir, "event_edges_deep_idcost.txt"))

deep_rows <- character(0)
for (iv in c("No", "previous", "next")) {
  for (nm in c("none", "max", "YujianBo")) {
    d <- TraMineRextras::seqedist(deep, idcost = deep_idcost, vparam = 0.25,
                                  interval = iv, norm = nm)
    deep_rows <- c(deep_rows, sprintf("%s,%s,%s", iv, nm,
                                      paste(format(as.vector(d), digits = 17),
                                            collapse = ";")))
  }
}
writeLines(c("# interval,norm,column-major distance matrix", deep_rows),
           file.path(out_dir, "event_edges_deep_seqedist.csv"))

## --------------------------------------------------------------------- ##
## 2. Degenerate corpora and boundary constraints.                        ##
## --------------------------------------------------------------------- ##
store <- list()

mk <- function(id, ts, ev) seqecreate(id = id, timestamp = ts, event = ev)

base <- mk(c("a", "a", "a", "b", "b", "c"), c(1, 2, 3, 1, 2, 1),
           c("A", "B", "C", "A", "C", "B"))
store <- record(store, "base|as.character", as.character(base))
store <- record(store, "base|levels", levels(base))
store <- record(store, "base|fsub_subseq",
                as.character(seqefsub(base, min.support = 1)$subseq))
store <- record(store, "base|fsub_count", seqefsub(base, min.support = 1)$data$Count)
store <- record(store, "base|fsub_support", seqefsub(base, min.support = 1)$data$Support)

one <- mk("a", 1, "A")
store <- record(store, "one_event|as.character", as.character(one))
store <- record(store, "one_event|fsub_subseq",
                as.character(seqefsub(one, min.support = 1)$subseq))
store <- record(store, "one_event|fsub_count", seqefsub(one, min.support = 1)$data$Count)
store <- record(store, "one_event|seqedist",
                TraMineRextras::seqedist(one, idcost = 1, vparam = 1))

same <- mk(rep(c("a", "b"), each = 3), rep(1, 6), rep(c("A", "B", "C"), 2))
store <- record(store, "same_ts|as.character", as.character(same))
store <- record(store, "same_ts|fsub_subseq",
                as.character(seqefsub(same, min.support = 1)$subseq))
store <- record(store, "same_ts|fsub_count", seqefsub(same, min.support = 1)$data$Count)
store <- record(store, "same_ts|fsub_support", seqefsub(same, min.support = 1)$data$Support)
store <- record(store, "same_ts|seqedist_next",
                TraMineRextras::seqedist(same, idcost = c(1, 1, 1), vparam = 1,
                                         interval = "next"))

unsorted <- mk(rep("a", 3), c(3, 1, 2), c("C", "A", "B"))
store <- record(store, "unsorted|as.character", as.character(unsorted))
store <- record(store, "unsorted|fsub_subseq",
                as.character(seqefsub(unsorted, min.support = 1)$subseq))

dup <- mk(rep("a", 3), c(1, 1, 2), c("A", "A", "B"))
store <- record(store, "dup_events|as.character", as.character(dup))
store <- record(store, "dup_events|fsub_subseq",
                as.character(seqefsub(dup, min.support = 1)$subseq))
store <- record(store, "dup_events|fsub_count", seqefsub(dup, min.support = 1)$data$Count)

neg <- mk(rep("a", 3), c(-2, -1, 3), c("A", "B", "C"))
store <- record(store, "neg_ts|as.character", as.character(neg))
store <- record(store, "neg_ts|fsub_subseq",
                as.character(seqefsub(neg, min.support = 1)$subseq))
store <- record(store, "neg_ts|fsub_count", seqefsub(neg, min.support = 1)$data$Count)
store <- record(store, "neg_ts|seqedist",
                TraMineRextras::seqedist(neg, idcost = c(1, 1, 1), vparam = 1))

# Boundary constraints. The zero-width cases still mine length-1 subsequences;
# the excluding cases are where R aborts in C.
store <- record(store, "cst|maxgap_zero",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(max.gap = 0))$subseq))
store <- record(store, "cst|windowsize_zero",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(window.size = 0))$subseq))
store <- record(store, "cst|cwin_windowsize_zero",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(count.method = "CWIN",
                                                                  window.size = 0))$subseq))
store <- record(store, "cst|cwin_windowsize_zero_count",
                seqefsub(base, min.support = 1,
                         constraint = seqeconstraint(count.method = "CWIN",
                                                     window.size = 0))$data$Count)
store <- record(store, "cst|agemax_zero",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(age.max = 0))$subseq))
store <- record(store, "cst|excludes_everything",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(age.min = 100,
                                                                  age.max = 200))$subseq))
store <- record(store, "cst|cwin_no_windowsize",
                as.character(seqefsub(base, min.support = 1,
                                      constraint = seqeconstraint(count.method = "CWIN"))$subseq))
store <- record(store, "cst|maxgap_gt_window", seqeconstraint(max.gap = 5, window.size = 2))
store <- record(store, "cst|agemin_gt_agemax", seqeconstraint(age.min = 5, age.max = 2))

# min.support / pmin.support / max.k boundaries.
store <- record(store, "ms|zero",
                as.character(seqefsub(base, min.support = 0)$subseq))
store <- record(store, "ms|negative",
                as.character(seqefsub(base, min.support = -1)$subseq))
store <- record(store, "ms|huge",
                as.character(seqefsub(base, min.support = 1000)$subseq))
store <- record(store, "ms|equal_n",
                as.character(seqefsub(base, min.support = 3)$subseq))
store <- record(store, "pms|zero",
                as.character(seqefsub(base, pmin.support = 0)$subseq))
store <- record(store, "pms|one",
                as.character(seqefsub(base, pmin.support = 1)$subseq))
store <- record(store, "pms|gt_one",
                as.character(seqefsub(base, pmin.support = 1.5)$subseq))
store <- record(store, "ms|none", as.character(seqefsub(base)$subseq))
store <- record(store, "maxk|zero",
                as.character(seqefsub(base, min.support = 1, max.k = 0)$subseq))
store <- record(store, "maxk|one",
                as.character(seqefsub(base, min.support = 1, max.k = 1)$subseq))

# seqedist degenerate arguments.
store <- record(store, "edist|zero_idcost",
                TraMineRextras::seqedist(base, idcost = c(0, 0, 0), vparam = 1))
store <- record(store, "edist|neg_idcost",
                TraMineRextras::seqedist(base, idcost = c(-1, -1, -1), vparam = 1))
store <- record(store, "edist|neg_vparam",
                TraMineRextras::seqedist(base, idcost = c(1, 1, 1), vparam = -1))
store <- record(store, "edist|wrong_idcost_len",
                TraMineRextras::seqedist(base, idcost = c(1, 1), vparam = 1))

# seqecmpgroup degenerate groupings -- R aborts on both; traminex refuses.
store <- record(store, "cmpgroup|single_group",
                dim(seqecmpgroup(seqefsub(base, min.support = 1),
                                 group = c("g", "g", "g"))$data))
store <- record(store, "cmpgroup|wrong_length_group",
                dim(seqecmpgroup(seqefsub(base, min.support = 1),
                                 group = c("g", "h"))$data))
zero_w <- mk(c("a", "a", "a", "b", "b", "c"), c(1, 2, 3, 1, 2, 1),
             c("A", "B", "C", "A", "C", "B"))
seqeweight(zero_w) <- c(0, 1, 1)
store <- record(store, "cmpgroup|zero_weight_level",
                dim(seqecmpgroup(seqefsub(zero_w, min.support = 1),
                                 group = c("g", "h", "h"))$data))
store <- record(store, "cmpgroup|two_groups_shape",
                dim(seqecmpgroup(seqefsub(base, min.support = 1),
                                 group = c("g", "g", "h"))$data))

edge <- data.frame(key = names(store), value = unlist(store, use.names = FALSE),
                   stringsAsFactors = FALSE)
write.csv(edge, file.path(out_dir, "event_edges_degenerate.csv"), row.names = FALSE)

cat("EVENT_ENGINE_EDGES_GOLDENS_DONE\n")
cat("deep rows:", length(deep_rows), " degenerate rows:", nrow(edge), "\n")
