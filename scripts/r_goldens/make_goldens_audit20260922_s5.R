# Golden generation: event-sequence behaviours from the 2026-09-22 audit (S5).
#
# fractional_event_times -- TraMineR builds each event sequence by inserting
# its events one at a time, in their input order, into a chain of gaps
# (tmrsequence.cpp:113-127, eventseq.cpp:21-48). With fractional times those
# repeated subtractions round differently from the differences of the sorted
# times, and the ages R re-adds from the gaps drive the "age" count method,
# the CWIN window edges, the CDIST equality tests and seqeage (audit S5-01).
# Consumed by tests/test_event_insertion_order_gaps.py.
#
# cmpgroup_factor -- seqecmpgroup's factor(group) semantics (audit S5-02)
# and its pvalue.limit 1:0 case (S5-07).
# Consumed by tests/test_seqecmpgroup_factor_groups.py.
#
# event_objects -- behaviours of the event objects: slicing (audit S5-03),
# seqecontain on a subsequence list (S5-06), gap formatting (S5-08), the
# count-method-6 warning (S5-09), scalar length/weight setters (S5-10).
# Consumed by tests/test_seqelist_slicing.py,
# tests/test_seqecontain_subsequence_list.py, tests/test_event_gap_format.py,
# tests/test_seqefsub_count_method_six.py and
# tests/test_seqelist_scalar_setters.py.
#
# samm_typology -- seqsammeha's names for a numeric typology (audit S5-12).
# Consumed by tests/test_seqsammeha_numeric_typology.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

num <- function(x) if (length(x)) paste(sprintf("%.17g", as.vector(x)), collapse = " ") else "none"
out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# fractional_event_times
# ---------------------------------------------------------------------------
lines <- character()
add <- function(...) lines <<- c(lines, paste(..., sep = "\t"))

# Two hand-made sequences (the audit's minimal reproductions) and a seeded
# corpus whose times sit on a 0.1 grid and arrive in shuffled order.
set.seed(20260922)
corpus <- data.frame(id = integer(), timestamp = numeric(), event = character())
for (i in 1:14) {
  k <- sample(4:7, 1)
  times <- round(runif(k, 0, 6), 1)
  events <- sample(c("a", "b", "c", "d"), k, replace = TRUE)
  if (i %% 3 == 0) {
    times <- c(times, times[1])
    events <- c(events, "d")
  }
  corpus <- rbind(corpus, data.frame(id = i, timestamp = times, event = events))
}
corpora <- list(
  minimal_age = data.frame(id = c(1, 1, 1), timestamp = c(0.1, 0.2, 2.0),
    event = c("a", "b", "c")),
  minimal_overlap = data.frame(id = c(1, 1, 1, 1, 2, 2, 2, 2),
    timestamp = c(0, 0.1, 0.2, 0.6, 0, 0.1, 1.0, 1.2),
    event = c("a", "a", "b", "b", "a", "b", "a", "b")),
  shuffled = corpus
)
build <- function(d) seqecreate(id = d$id, timestamp = d$timestamp, event = d$event)

for (name in names(corpora)) {
  d <- corpora[[name]]
  for (r in seq_len(nrow(d))) {
    add("corpus", name, d$id[r], sprintf("%.17g", d$timestamp[r]), d$event[r])
  }
  e <- build(d)
  add("levels", name, paste(levels(e), collapse = " "))
  aschar <- as.character(e)
  for (i in seq_along(aschar)) add("aschar", name, i, aschar[i])
  for (code in seq_along(levels(e))) {
    add("seqeage", name, code, num(TraMineR:::seqeage(e, code)))
  }
}

# Named subsequences counted by seqeapplysub under each method.
applied <- list(
  list("minimal_age", c("(c)", "(a)-(c)"), "age", NULL),
  list("minimal_overlap", "(a)-(b)", "CDIST", NULL),
  list("minimal_overlap", "(a)-(b)", NULL, seqeconstraint(window.size = 1, count.method = "CWIN")),
  list("shuffled", c("(a)", "(b)", "(c)", "(d)", "(a)-(b)", "(b)-(c)", "(c)-(d)", "(d)-(a)"), "age", NULL),
  list("shuffled", c("(a)-(b)", "(b)-(c)", "(c)-(d)", "(d)-(a)", "(a)-(a)"), "CDIST", NULL),
  list("shuffled", c("(a)-(b)", "(b)-(c)", "(c)-(d)", "(d)-(a)", "(a)-(a)"), NULL,
    seqeconstraint(window.size = 1.1, count.method = "CWIN")),
  list("shuffled", c("(a)-(b)", "(b)-(c)", "(c)-(d)", "(d)-(a)"), NULL,
    seqeconstraint(max.gap = 0.3, count.method = "CDIST_O"))
)
for (k in seq_along(applied)) {
  spec <- applied[[k]]
  e <- build(corpora[[spec[[1]]]])
  fs <- seqefsub(e, str.subseq = spec[[2]])
  m <- if (is.null(spec[[4]])) seqeapplysub(fs, method = spec[[3]]) else
    seqeapplysub(fs, constraint = spec[[4]])
  add("apply", k, spec[[1]], paste(spec[[2]], collapse = " "), num(m))
}

# Frequent-subsequence mining on the shuffled corpus under constraints whose
# window, gap and age limits fall on the 0.1 grid.
mined <- list(
  cdist = seqeconstraint(count.method = "CDIST"),
  cwin = seqeconstraint(window.size = 1.1, count.method = "CWIN"),
  maxgap = seqeconstraint(max.gap = 0.7, count.method = "CDIST_O"),
  agemin = seqeconstraint(age.min = 0.3, age.max = 4.7, count.method = "CDIST_O"),
  cminwin = seqeconstraint(window.size = 2.3, count.method = "CMINWIN")
)
e <- build(corpora$shuffled)
for (tag in names(mined)) {
  fs <- seqefsub(e, min.support = 3, constraint = mined[[tag]], max.k = 3)
  sub <- as.character(fs$subseq)
  add("mined", tag, length(sub))
  for (r in seq_along(sub)) {
    add("mined_row", tag, r, sub[r], sprintf("%.17g", fs$data$Support[r]),
      sprintf("%.17g", fs$data$Count[r]))
  }
}

# The event edit distance walks the same gap chain.
for (interval in c("No", "previous", "next")) {
  d <- seqedist(e, idcost = rep(1, length(levels(e))), vparam = 0.3, interval = interval,
    norm = "none")
  add("seqedist", interval, num(d))
}

writeLines(c(
  "# kind\tfields (tab separated); numbers are sprintf('%.17g'), matrices column-major",
  lines
), file.path(out_dir, "audit20260922_s5_fractional_event_times.txt"))

# ---------------------------------------------------------------------------
# cmpgroup_factor -- seqecmpgroup.R:62 turns the grouping into factor(group):
# numbers sort as numbers and print as R prints them, a factor keeps its level
# order and drops unused levels, character levels follow the collation, and an
# NA group is left out of each xtabs table while still counting in
# length(group) (audit S5-02).
# ---------------------------------------------------------------------------
lines <- character()
set.seed(20260923)
corpus <- data.frame(id = integer(), timestamp = numeric(), event = character())
for (i in 1:18) {
  k <- sample(3:6, 1)
  corpus <- rbind(corpus, data.frame(id = i, timestamp = sort(sample(0:9, k)),
    event = sample(c("a", "b", "c", "d"), k, replace = TRUE)))
}
for (r in seq_len(nrow(corpus))) {
  add("corpus", corpus$id[r], sprintf("%.17g", corpus$timestamp[r]), corpus$event[r])
}
e <- seqecreate(id = corpus$id, timestamp = corpus$timestamp, event = corpus$event)
fs <- seqefsub(e, min.support = 3, max.k = 2)
i <- seq_len(length(e))
with_na <- ifelse(i %% 3 == 0, "x", "y")
with_na[c(2, 5, 11)] <- NA
groups <- list(
  numeric = ifelse(i %% 3 == 0, 10, 2),
  character = c("b", "B", "a")[1 + i %% 3],
  factor = factor(ifelse(i %% 2 == 0, "x", "y"), levels = c("z", "y", "x")),
  missing = with_na,
  logical = i %% 2 == 0,
  fractional = ifelse(i %% 2 == 0, 1.5, 1e5)
)
for (tag in names(groups)) {
  res <- seqecmpgroup(fs, group = groups[[tag]])
  add("labels", tag, paste(res$labels, collapse = " "))
  add("columns", tag, paste(colnames(res$data), collapse = " "))
  sub <- as.character(res$subseq)
  for (r in seq_along(sub)) add("row", tag, r, sub[r], num(unlist(res$data[r, ])))
}
# When no p-value passes pvalue.limit, seqecmpgroup.R:88 keeps 1:0 = c(1, 0)
# rows, i.e. still returns the top subsequence (audit S5-07; traminex returns
# none, a recorded divergence).
res <- seqecmpgroup(fs, group = groups$numeric, pvalue.limit = 1e-300)
add("pvalue_limit", "1e-300", nrow(res$data), paste(as.character(res$subseq), collapse = " "),
  num(res$data$p.value))
writeLines(c(
  "# kind\tfields (tab separated); numbers are sprintf('%.17g')",
  lines
), file.path(out_dir, "audit20260922_s5_cmpgroup_factor.txt"))

# ---------------------------------------------------------------------------
# event_objects -- behaviours of the event-sequence objects themselves.
# ---------------------------------------------------------------------------
lines <- character()
tse <- function(x) {
  d <- TraMineR:::seqe2TSE(x)
  if (!nrow(d)) return("none")
  paste(paste(d$id, sprintf("%.17g", d$timestamp), as.character(d$event), sep = ","),
    collapse = " ")
}

# slice (audit S5-03): "[.seqelist" keeps the selected sequences, and with
# them exactly their events -- also when a sequence has no event (id 2 below
# holds only its end event) or when as.integer merges fractional ids.
with_end <- seqecreate(id = c(1, 1, 2, 3, 3), timestamp = c(0, 1, 5, 0, 2),
  event = c("a", "b", "END", "a", "c"), end.event = "END")
merged <- seqecreate(id = c(1.2, 1.7, 2.0), timestamp = c(0, 1, 0), event = c("a", "b", "c"))
slices <- list(
  list("with_end", with_end, 2:3), list("with_end", with_end, 1:2),
  list("with_end", with_end, 3), list("with_end", with_end, 2),
  list("merged", merged, 2), list("merged", merged, 1)
)
for (s in slices) {
  sub <- s[[2]][s[[3]]]
  add("slice", s[[1]], paste(s[[3]], collapse = " "), paste(seqeid(sub), collapse = " "),
    num(seqelength(sub)), tse(sub))
}

# contain (audit S5-06): seqecontain.R:6 also takes a subsequence list and
# tests each subsequence.
e <- seqecreate(id = c(1, 1, 2, 2, 2), timestamp = c(0, 5, 0, 1, 3),
  event = c("a", "b", "a", "c", "b"))
fs <- seqefsub(e, min.support = 1)
add("contain_subseq", paste(as.character(fs$subseq), collapse = " "))
contains <- list(list("c", FALSE), list(c("a", "b"), TRUE), list(c("c", "zzz"), FALSE),
  list(c("a", "c"), TRUE))
for (case in contains) {
  add("contain", paste(case[[1]], collapse = " "), case[[2]],
    paste(seqecontain(fs, case[[1]], unknown.exclude = case[[2]]), collapse = " "))
}

# gap_format (audit S5-08): Sequence::sprint renders every gap, and the
# trailing gap to the end of observation, with R's format() (eventseq.cpp:90),
# so 100000 prints as 1e+05 and 0.0001 as 1e-04; seqeapplysub names its rows
# with the same strings.
e <- seqecreate(id = c(1, 1, 2, 2, 2, 3, 3, 4, 4, 5, 5),
  timestamp = c(0, 100000, 0, 0.0001, 3, 0, 123456, 0, 1 / 3, 0, 1234567.8),
  event = c("a", "b", "a", "c", "b", "a", "b", "a", "b", "a", "b"))
add("gap_format", "created", paste(as.character(e), collapse = " "))
seqelength(e) <- c(200000, 5, 123460, 1, 2e6)
add("gap_format", "with_length", paste(as.character(e), collapse = " "))
fs <- seqefsub(e, str.subseq = "(a)-(b)")
add("gap_format", "applysub_rows", paste(rownames(seqeapplysub(fs, method = "presence")),
  collapse = " "))

# count_method_6 (audit S5-09): seqefsub.R:26-28 warns that count method 6 is
# for internal use, then mines anyway.
e <- seqecreate(id = c(1, 1, 2, 2, 2, 3, 3), timestamp = c(0, 2, 0, 1, 3, 1, 4),
  event = c("a", "b", "a", "c", "b", "b", "c"))
seen <- character()
fs <- withCallingHandlers(
  seqefsub(e, min.support = 1, constraint = seqeconstraint(count.method = 6)),
  warning = function(w) {
    seen <<- c(seen, conditionMessage(w))
    invokeRestart("muffleWarning")
  }
)
add("count_method_6", "warning", paste(seen, collapse = " | "))
add("count_method_6", "subseq", paste(as.character(fs$subseq), collapse = " "))
add("count_method_6", "count", num(fs$data$Count))

# scalar_setters (audit S5-10): seqelength<- and seqeweight<- compare
# length(value) with length(eseq), so a single number sets a one-sequence
# list and is refused for a longer one.
one <- seqecreate(id = c(1, 1), timestamp = c(0, 2), event = c("a", "b"))
seqelength(one) <- 5
seqeweight(one) <- 2.5
add("scalar_setters", "one", num(seqelength(one)), num(seqeweight(one)), as.character(one))
three <- seqecreate(id = c(1, 2, 3), timestamp = c(0, 0, 0), event = c("a", "b", "a"))
for (setter in c("length", "weight")) {
  err <- tryCatch({
    if (setter == "length") seqelength(three) <- 5 else seqeweight(three) <- 2.5
    "none"
  }, error = function(e) trimws(conditionMessage(e)))
  add("scalar_setters", setter, err)
}

writeLines(c(
  "# kind\tfields (tab separated); numbers are sprintf('%.17g')",
  lines
), file.path(out_dir, "audit20260922_s5_event_objects.txt"))

# ---------------------------------------------------------------------------
# samm_typology (audit S5-12) -- seqsammeha stores as.character(typology) and
# names the indicator columns from its values, so a numeric typology 1, 2
# gives SAMM1/SAMM2 and levels None 1 2. R fills those columns by position
# (ret[, tt]), so 2, 1 swaps the indicators and 10, 20 fails; that case is
# recorded as a divergence.
# ---------------------------------------------------------------------------
lines <- character()
m <- rbind(c("A", "A", "B", "B", "C", "C", "A"), c("B", "B", "A", "A", "A", "C", "C"),
  c("C", "A", "B", "B", "B", "C", "C"), c("A", "B", "C", "A", "B", "C", "A"),
  c("B", "C", "C", "C", "A", "A", "A"))
s <- suppressMessages(seqdef(m, alphabet = c("A", "B", "C")))
sm <- seqsamm(s, sublength = 3)
typologies <- list(whole = c(1, 2, 1, 2), fractional = c(1.5, 2.5, 1.5, 2.5),
  reversed = c(2, 1, 2, 1))
for (tag in names(typologies)) {
  eha <- seqsammeha(sm, spell = "A", typology = typologies[[tag]])
  cols <- setdiff(grep("^SAMM", colnames(eha), value = TRUE), "SAMMtypology")
  add("eha", tag, paste(cols, collapse = " "), paste(levels(eha$SAMMtypology), collapse = " "),
    paste(as.character(eha$SAMMtypology), collapse = " "))
  for (cc in cols) add("indicator", tag, cc, num(eha[[cc]]))
}
writeLines(c(
  "# kind\tfields (tab separated); numbers are sprintf('%.17g')",
  lines
), file.path(out_dir, "audit20260922_s5_samm_typology.txt"))

cat("AUDIT20260922_S5_GOLDENS_DONE\n")
