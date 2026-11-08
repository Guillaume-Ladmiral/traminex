# Live-R goldens for batch 11's event and discrete-rule items (TX-AUD-20261019-043, -044, -083;
# TX-AUD-20261010-076, -077; TX-AUD-20261011-057).
#
# Each row is case, field, value. The seqerulesdisc cases carry their simulated input (field tsef, rows
# IDPERS,A,AST,B,BST joined by '|', as TraMineRextras' seqerulesdisc help simulates it) and R's rule labels (field
# labels, joined by ' ;; '), whose numbers R prints with paste(round(...)) (seqerulesdisc.R:89).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_event_arguments.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR)); suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14", as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)

# --- seqerulesdisc labels (043) ---
erd <- function(case, seed, N, prob) {
  set.seed(seed)
  ids <- 1:N
  AST <- rbinom(N, 1, 0.5); A <- ifelse(AST == 1, sample(16:24, N, TRUE), 30)
  BST <- integer(N); B <- rep(30, N)
  for (i in 1:N) {
    if (AST[i] == 1 && runif(1) < prob) { BST[i] <- 1; B[i] <- A[i] + sample(1:3, 1) }
    else if (runif(1) < 0.25) { BST[i] <- 1; B[i] <- sample(16:28, 1) }
  }
  tsef <- data.frame(IDPERS = ids, A = A, AST = AST, B = B, BST = BST)
  agemin <- data.frame(IDPERS = ids, AGE = 15); agemax <- data.frame(IDPERS = ids, AGE = 30)
  dd <- suppressMessages(createdatadiscrete(ids, tsef, c("A", "B"), agemin, agemax))
  tse <- rbind(data.frame(id = ids[AST == 1], time = A[AST == 1], event = "A"), data.frame(id = ids[BST == 1], time = B[BST == 1], event = "B"))
  tse <- tse[order(tse$id, tse$time), ]
  seqe <- seqecreate(id = tse$id, timestamp = tse$time, event = tse$event)
  fs <- seqefsub(seqe, pmin.support = 0.01)
  r <- suppressMessages(seqerulesdisc(fs, dd, tsef, pvalue = 1, adjust = FALSE))
  add(case, "tsef", paste(apply(tsef, 1, paste, collapse = ","), collapse = "|"))
  add(case, "labels", paste(as.character(r$Hazards), collapse = " ;; "))
}
erd("erd_seed11_n80", 11, 80, 0.45)
erd("erd_seed3_n150", 3, 150, 0.6)
erd("erd_seed7_n60", 7, 60, 0.3)

# --- seqerulesdisc labels with large hazard ratios (043 follow-up, batch-11 review M-5), both links (fields labels for
# cloglog, labels_logit) and R's warnings. erd_large_hr: 40 people with A (36 have B a year later) and 200 without A
# (one has B at 25); a well-conditioned fit, HR about 1.5e3. erd_separated: the batch-11 reproduction's perfectly
# separated data (B two years after A and only after A); HR about 2.05e10, which R prints with 13 significant digits.
erd_fixed <- function(case, tsef) {
  ids <- tsef$IDPERS
  dd <- suppressMessages(createdatadiscrete(ids, tsef, c("A", "B"), data.frame(IDPERS = ids, AGE = 15), data.frame(IDPERS = ids, AGE = 30)))
  tse <- rbind(data.frame(id = ids[tsef$AST == 1], time = tsef$A[tsef$AST == 1], event = "A"),
               data.frame(id = ids[tsef$BST == 1], time = tsef$B[tsef$BST == 1], event = "B"))
  tse <- tse[order(tse$id, tse$time), ]
  fs <- seqefsub(seqecreate(id = tse$id, timestamp = tse$time, event = tse$event), pmin.support = 0.01)
  add(case, "tsef", paste(apply(tsef, 1, paste, collapse = ","), collapse = "|"))
  for (lk in c("cloglog", "logit")) {
    warns <- character(0)
    r <- withCallingHandlers(suppressMessages(seqerulesdisc(fs, dd, tsef, pvalue = 1, adjust = FALSE, link = lk)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    suffix <- if (lk == "cloglog") "" else "_logit"
    add(case, paste0("labels", suffix), paste(as.character(r$Hazards), collapse = " ;; "))
    add(case, paste0("warnings", suffix), paste(warns, collapse = " | "))
  }
}
ids <- 1:240
A <- ifelse(ids <= 40, 18 + (ids %% 6), 30)
erd_fixed("erd_large_hr", data.frame(IDPERS = ids, A = A, AST = as.integer(ids <= 40),
                                     B = ifelse(ids <= 36, A + 1, ifelse(ids == 41, 25, 30)), BST = as.integer(ids <= 36 | ids == 41)))
ids <- 1:60
A <- ifelse(ids <= 30, 18 + (ids %% 5), 30)
erd_fixed("erd_separated", data.frame(IDPERS = ids, A = A, AST = as.integer(ids <= 30), B = ifelse(ids <= 30, A + 2, 30),
                                      BST = as.integer(ids <= 30)))

# --- seqecontain with numeric event values (044): events 1 and 2; R matches factor(event.list, levels = dict)
# (seqecontain.R:10), so 1.0 and "1" match the event 1. Values are R's TRUE/FALSE joined by spaces.
ec <- seqecreate(id = c(1, 1, 2), timestamp = c(0, 1, 0), event = c(1, 2, 2))
contain <- function(case, value, exclude) add(case, "contain", paste(seqecontain(ec, value, unknown.exclude = exclude), collapse = " "))
contain("contain_double_1", 1.0, FALSE)
contain("contain_integer_1", 1L, FALSE)
contain("contain_text_1", "1", FALSE)
contain("contain_double_2_5", 2.5, FALSE)
contain("contain_double_1_and_2", c(1.0, 2.0), TRUE)
contain("contain_double_2_exclude", 2.0, TRUE)
contain("contain_double_1_and_unknown_exclude", c(1.0, 3.0), TRUE)

# --- seqerules sortv (083): R orders the rule text with order(), R's collation (seqerules.R:77-78); a (lower-case)
# sorts before B. Values are the Rules column joined by ' ;; '.
rt <- data.frame(id = c(1, 1, 1, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 6, 6, 6),
                 time = c(0, 1, 2, 0, 2, 0, 1, 3, 0, 1, 1, 2, 3, 0, 2, 3),
                 event = c("a", "B", "c", "a", "c", "B", "c", "a", "a", "B", "c", "a", "B", "B", "a", "c"))
rs <- seqefsub(seqecreate(id = rt$id, timestamp = rt$time, event = rt$event), pmin.support = 0.3)
rules <- function(case, ...) add(case, "rules", paste(as.character(TraMineR:::seqerules(rs, ...)$Rules), collapse = " ;; "))
rules("rules_sort_rule", sortv = "rule")
rules("rules_sort_rule_decreasing", sortv = "rule", decreasing = TRUE)
rules("rules_sort_support_decreasing", sortv = "support", decreasing = TRUE)
rules("rules_unsorted")

# --- seqeconstraint count.method (077): R stops unless count.method is in 1:6 or a known name (seqeconstraint.R:18-22),
# so 2.5, 0 and 7 are refused and TRUE counts as 1. Values are the stored count method as a number, R's message, or
# NULL: the text "3" passes the check but its name lookup stores nothing, and seqefsub then fails on it.
cm <- function(case, value) add(case, "count_method", tryCatch({
  stored <- seqeconstraint(count.method = value)$count.method
  if (is.null(stored)) "NULL" else format(as.numeric(stored))
}, error = conditionMessage))
cm("cm_double_2_5", 2.5)
cm("cm_double_0", 0)
cm("cm_double_7", 7)
cm("cm_double_6", 6)
cm("cm_double_2", 2)
cm("cm_logical_true", TRUE)
cm("cm_name_cwin", "CWIN")
cm("cm_text_3", "3")
fsub_error <- function(case, constraint) add(case, "seqefsub_error", tryCatch({
  seqefsub(ec, pmin.support = 0.1, constraint = constraint); "none"
}, error = conditionMessage))
# R's message for count method 3 (CWIN) names 'CMINWIN' (seqefsub.R), a slip: count method 4 is CMINWIN.
fsub_error("seqefsub_cwin_without_window", seqeconstraint(count.method = "CWIN"))
fsub_error("seqefsub_text_3", seqeconstraint(count.method = "3"))

# --- seqeapplysub method (057): after the names, R uses a method only when method %in% 1:5 (seqeapplysub.R:28), a
# match by value, so 2.5, 0 and 6 leave the constraint's count method and the text "2" counts as 2 (%in% compares
# text with text). Values are the count matrix, rows joined by '|', cells by ','; the c2 cases count under
# count.method 2.
ap <- seqecreate(id = c(1, 1, 1, 2, 2, 3), timestamp = c(0, 1, 2, 0, 3, 1), event = c("A", "B", "A", "A", "B", "B"))
apsub <- seqefsub(ap, str.subseq = c("(A)", "(A)-(B)"))
apply_case <- function(case, method, constraint = NULL) {
  m <- seqeapplysub(apsub, method = method, constraint = constraint)
  add(case, "applysub", paste(apply(m, 1, paste, collapse = ","), collapse = "|"))
}
apply_case("apply_method_null", NULL)
apply_case("apply_method_double_2_5", 2.5)
apply_case("apply_method_double_2", 2)
apply_case("apply_method_double_0", 0)
apply_case("apply_method_double_6", 6)
apply_case("apply_method_logical_true", TRUE)
apply_case("apply_method_text_2", "2")
apply_case("apply_method_text_2_0", "2.0")
apply_case("apply_method_name_count", "count")
apply_case("apply_c2_method_null", NULL, seqeconstraint(count.method = 2))
apply_case("apply_c2_method_double_1_5", 1.5, seqeconstraint(count.method = 2))
apply_case("apply_c2_method_double_1", 1, seqeconstraint(count.method = 2))

# --- seqelist indexing (058): "[.seqelist" takes any index vector (seqe.R:22-32): positions in any order, a logical
# mask, or none. Sequence 2 holds only its end event. Values are the subset's ids, weights and events (seqe2TSE rows
# id,timestamp,event), each joined by spaces; R's positions are 1-based.
el <- seqecreate(id = c(1, 1, 2, 3, 3), timestamp = c(0, 1, 5, 0, 2), event = c("a", "b", "END", "a", "c"),
                 end.event = "END")
seqeweight(el) <- c(1, 2.5, 4)
subset_case <- function(case, index) {
  sub <- el[index]
  # seqe2TSE stops on an empty list (its factor() has no labels), so an empty subset has no events to list.
  d <- if (length(sub)) TraMineR:::seqe2TSE(sub) else data.frame()
  add(case, "length", length(sub))
  add(case, "ids", paste(seqeid(sub), collapse = " "))
  add(case, "weights", paste(seqeweight(sub), collapse = " "))
  add(case, "events", if (nrow(d)) paste(paste(d$id, d$timestamp, as.character(d$event), sep = ","), collapse = " ") else "none")
}
subset_case("subset_positions_1_3", c(1, 3))
subset_case("subset_positions_3_1", c(3, 1))
subset_case("subset_mask_true_false_true", c(TRUE, FALSE, TRUE))
subset_case("subset_positions_none", integer(0))

# --- seqefsub and the caller's sequences (076). R stops on a min.support that no subsequence reaches. With
# weighted = FALSE R sets the caller's weights to 1 (seqefsub.R:12-15) and restores them only after mining
# (:84-86), so the str.subseq route, which returns first (:56), leaves them at 1. R's event sequences are shared
# objects, so writing a slice's weights changes the source's. Weights joined by spaces; data rows Support,Count by '|'.
fw <- function() {
  s <- seqecreate(id = c(1, 1, 2, 2, 2, 3, 4, 4), timestamp = c(0, 3, 1, 2, 5, 4, 0, 2),
                  event = c("A", "B", "A", "C", "B", "A", "B", "C"))
  seqeweight(s) <- c(1, 2.5, 1, 1)
  s
}
s <- fw()
add("fsub_min_support_100", "seqefsub_error", tryCatch({ seqefsub(s, min.support = 100); "none" }, error = conditionMessage))
s <- fw()
fs <- seqefsub(s, str.subseq = c("(A)-(B)", "(B)"), weighted = FALSE)
add("fsub_str_unweighted", "caller_weights", paste(seqeweight(s), collapse = " "))
add("fsub_str_unweighted", "data", paste(fs$data$Support, fs$data$Count, sep = ",", collapse = "|"))
s <- fw()
fs <- seqefsub(s, min.support = 1, weighted = FALSE)
add("fsub_min_support_unweighted", "caller_weights", paste(seqeweight(s), collapse = " "))
s <- fw()
sl <- s[1:2]
seqeweight(sl) <- c(9, 9)
add("slice_weights_written", "source_weights", paste(seqeweight(s), collapse = " "))

# --- paste(round(x, d)) at decimal ties (043 follow-up, batch-11 review M-4): R's round() (nmath/fround.c) picks the
# nearer of the two candidates x10 = floor or ceiling of x * 10^d, then / 10^d, and the even one on an exact tie, which
# is not Python's round. Ties at d = 2 (HR labels) and d = 4 (p labels), and edge values; x in %.17g.
tie2 <- (0:999 + 0.5) / 100
tie4 <- (0:999 + 0.5) / 10000
edges <- c(6.165, 5e-05, 0.00015, 0.00025, 2.675, 0.285, 1.005, 1e5, 123456.789, 1e15, 2.0529e10, 4.9999e-5, 1e-4)
round_case <- function(case, x, d) {
  add(case, "x", paste(sprintf("%.17g", x), collapse = " "))
  add(case, "rounded", paste(paste(round(x, digits = d)), collapse = " "))
}
round_case("round_ties_d2", tie2, 2)
round_case("round_ties_d4", tie4, 4)
round_case("round_edges_d2", edges, 2)
round_case("round_edges_d4", edges, 4)

write.csv(do.call(rbind, rows), file.path(out_dir, "event_arguments.csv"), row.names = FALSE, quote = TRUE)
