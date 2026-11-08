# Golden generation: small paths of what the other libraries use that no test of this project
# ran (supervision repository, docs/TRAMINEX_COMPONENT_GAPS.md, 2026-11-02): representatives of
# a single object, and what plot.stslist refuses as sortv.
#
# Values are written with sprintf("%.17g"), so R's NaN is "NaN".
# Consumed by tests/test_r_goldens_consumed_items.py.
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

## One object: nothing to choose between, the path where traminex's near-tie re-summing stops
## early (traminex/diss_p1.py, _exact_near_ties).
add_representatives <- function(case, result) {
  add(case, "representatives", num(as.vector(result)))
  add(case, "scores", num(attr(result, "Scores")))
  add(case, "statistics", num(t(as.matrix(attr(result, "Statistics")))))
}
for (criterion in c("density", "freq", "dist")) {
  add_representatives(paste0("dissrep_single_", criterion),
                      suppressMessages(dissrep(matrix(0, 1, 1), criterion = criterion)))
}
one <- suppressMessages(seqdef(rbind(c("A", "B", "B")), alphabet = c("A", "B")))
single <- suppressMessages(seqrep(one, diss = matrix(0, 1, 1)))
add("seqrep_single", "statistics", num(t(as.matrix(attr(single, "Statistics")))))

## sortv: plot.stslist.R:62-76.
three <- suppressMessages(seqdef(rbind(c("A", "B", "B"), c("B", "B", "A"), c("A", "A", "A")), alphabet = c("A", "B")))
pdf(NULL)
refusal <- function(sortv) {
  tryCatch({
    plot(three, sortv = sortv)
    stop("R did not refuse")
  }, error = function(e) conditionMessage(e))
}
add("sortv_unknown_string", "refusal", refusal("not-a-method"))
add("sortv_wrong_length", "refusal", refusal(c(1, 2)))
invisible(dev.off())
add("sortv_integers", "order", num(order(c(3L, 1L, 2L))))

writeLines(c("# representatives of a single object, and plot.stslist's sortv", "\"case\",\"key\",\"value\"", rows),
           file.path(out_dir, "consumed_items.csv"))
