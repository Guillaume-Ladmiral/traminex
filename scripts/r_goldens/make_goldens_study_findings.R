# Golden generation: what R does in the cases of the 2026-10-25 study's findings F-01, F-04 and
# F-05, and the one-sequence sortv case (supervision repository, docs/FIXES_2026-11-07.md).
#
# - F-04: subsequence distances where every pair of a slice has an empty side (NMS raised
#   IndexError in traminex; R returns these matrices).
# - F-01: seqdef(states=) with new names that are also old ones. R renames one state at a time
#   (seqdef.R:78-81) and merges distinct states; traminex maps them all at once. Rows with names
#   disjoint from the alphabet are recorded too, where both agree.
# - F-05: R's stop on weights with a zero total, in dissrep and seqrep.
# - plot.stslist with one sequence draws with any sortv string (plot.stslist.R:62-76).
#
# Values are written with sprintf("%.17g"). Consumed by tests/test_r_goldens_study_findings.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

num <- function(x) paste(sprintf("%.17g", as.vector(x)), collapse = " ")
rows <- character(0)
add <- function(case, key, value) {
  rows <<- c(rows, paste(shQuote(case, "cmd"), shQuote(key, "cmd"), shQuote(value, "cmd"), sep = ","))
}
refusal <- function(expr) {
  tryCatch({
    force(expr)
    stop("R did not refuse")
  }, error = function(e) conditionMessage(e))
}

## F-04: empty sequences beside one-cell and two-cell ones.
three <- suppressMessages(seqdef(data.frame(a = c(NA, "B", NA)), right = "DEL"))
add("nms_empty_three", "distances", num(suppressMessages(TraMineR::seqdist(three, method = "NMS"))))
add("nms_soft_empty_three", "distances",
    num(suppressMessages(TraMineR::seqdist(three, method = "NMS", prox = matrix(1, 1, 1)))))
four <- suppressMessages(seqdef(data.frame(a = c(NA, "B", NA, "B"), b = c(NA, NA, NA, "C")), right = "DEL"))
for (method in c("NMS", "NMSMST", "SVRspell")) {
  add(paste0(tolower(method), "_empty_four"), "distances",
      num(suppressMessages(TraMineR::seqdist(four, method = method))))
}

## F-01: R's states= renaming.
states_rows <- function(case, data, ...) {
  s <- suppressMessages(seqdef(data, ...))
  add(case, "alphabet", paste(alphabet(s), collapse = " "))
  for (i in seq_len(nrow(s))) add(case, paste0("row", i), paste(as.character(unlist(s[i, ])), collapse = " "))
}
states_rows("states_shift", rbind(c("1", "2", "3"), c("3", "2", "1")), states = c("2", "3", "4"))
states_rows("states_swap", rbind(c("A", "B"), c("B", "A")), alphabet = c("A", "B"), states = c("B", "A"))
states_rows("states_disjoint", rbind(c("1", "2", "3"), c("3", "2", "1")), states = c("x", "y", "z"))

## F-05: weights with a zero total.
zero <- suppressMessages(seqdef(rbind(c("A", "B"), c("B", "A"), c("A", "A")), weights = c(0, 0, 0)))
zero_diss <- suppressMessages(TraMineR::seqdist(zero, method = "HAM"))
for (criterion in c("density", "freq", "dist")) {
  add(paste0("seqrep_zero_weights_", criterion), "refusal",
      refusal(suppressMessages(TraMineR::seqrep(zero, diss = zero_diss, criterion = criterion))))
}
add("dissrep_zero_weights", "refusal",
    refusal(suppressMessages(TraMineR::dissrep(matrix(0, 3, 3), weights = c(0, 0, 0)))))

## plot.stslist with one sequence: a string that names no method is a sort vector of length 1.
one <- suppressMessages(seqdef(rbind(c("A", "B", "B")), alphabet = c("A", "B")))
pdf(NULL)
drawn <- tryCatch({
  plot(one, sortv = "not-a-method")
  "drawn"
}, error = function(e) conditionMessage(e))
invisible(dev.off())
add("sortv_one_sequence_unknown_string", "outcome", drawn)

writeLines(c("# what R does in the cases of study findings F-01, F-04, F-05 and one-sequence sortv",
             "\"case\",\"key\",\"value\"", rows),
           file.path(out_dir, "study_findings.csv"))
