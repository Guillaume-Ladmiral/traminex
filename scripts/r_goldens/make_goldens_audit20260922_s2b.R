# Golden generation: seqMD CAT costs when a combined state has a void
# component (2026-09-22 audit, finding S2-19).
#
# seqMD.R builds each CAT cost by matching every component of a combined
# state in its channel's cost alphabet. A void component ("%", from a
# shorter channel or a channel without missing values) has no entry there,
# so match() gives NA and every cost pairing that state with another one is
# NA; the diagonal keeps its initial 0 and a scalar indel is unaffected.
# what = "cost" returns that matrix; what = "diss" stops with "NA indel
# and/or substitution CAT costs, cannot compute MD distances!".
#
# One record per line: case, field, value (tab-separated). Matrices are
# written column-major, cells separated by "|"; NA is written "NA".
#
# Consumed by tests/test_seqmd_void_component_costs.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(Sys.getlocale("LC_COLLATE") == "C.UTF-8", capabilities("ICU"))

lines <- character()
record <- function(case, field, value) {
  value <- as.character(value)
  if (length(value) == 0) value <- "none"
  lines <<- c(lines, paste(case, field, paste(value, collapse = "|"), sep = "\t"))
}
number <- function(x) ifelse(is.na(x), "NA", sprintf("%.17g", x))
record_cost <- function(case, cost) {
  record(case, "alphabet", attr(cost, "alphabet"))
  record(case, "dim", dim(cost))
  record(case, "sm", number(as.vector(unclass(cost))))
  record(case, "indel", number(attr(cost, "indel")))
}
record_error <- function(case, expr) {
  message <- tryCatch({
    suppressWarnings(suppressMessages(expr))
    "none"
  }, error = function(e) trimws(conditionMessage(e)))
  record(case, "error", message)
}
quiet <- function(x) suppressWarnings(suppressMessages(x))

# Unequal widths, the shorter channel with a missing value.
w1 <- quiet(seqdef(data.frame(a = c("A", "B"), b = c(NA, "A"), c = c("B", "B")),
  alphabet = c("A", "B"), right = NA))
w2 <- quiet(seqdef(data.frame(a = c("x", "y"), b = c("y", "y"), c = c("x", "y"),
  d = c("y", "x"), e = c("x", NA)), alphabet = c("x", "y")))
record_cost("unequal_constant", quiet(seqMD(list(w1, w2), sm = "CONSTANT", what = "cost")))
record_error("unequal_diss_om", seqMD(list(w1, w2), method = "OM", sm = "CONSTANT",
  what = "diss"))

# Equal widths, no missing value, a trailing void in one channel.
e1 <- quiet(seqdef(data.frame(a = c("A", "B"), b = c("A", "B"), c = c("B", NA))))
e2 <- quiet(seqdef(data.frame(a = c("x", "y"), b = c("y", "y"), c = c("x", "y"))))
record_cost("void_constant", quiet(seqMD(list(e1, e2), sm = "CONSTANT", what = "cost")))
record_cost("void_mean_cweight", quiet(seqMD(list(e1, e2), sm = "CONSTANT", what = "cost",
  link = "mean", cweight = c(2, 0.5))))
record_cost("void_trate_time_varying", quiet(seqMD(list(e1, e2),
  sm = list(seqcost(e1, "TRATE", time.varying = TRUE)$sm,
    seqcost(e2, "TRATE", time.varying = TRUE)$sm),
  what = "cost")))
record_error("void_diss_ham", seqMD(list(e1, e2), method = "HAM", what = "diss"))
record_error("void_diss_dhd", seqMD(list(e1, e2), method = "DHD", what = "diss"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
writeLines(c("# case\tfield\tvalue", lines), file.path(out_dir, "audit20260922_s2b_cases.txt"))
cat("AUDIT20260922_S2B_GOLDENS_DONE\n")
