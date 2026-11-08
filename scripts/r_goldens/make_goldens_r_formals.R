# Golden generation: every formal argument default of every exported TraMineR
# 2.2-14 / TraMineRextras 0.6.9 function.
#
# Cycle-7 W-07. A golden that always passes an argument explicitly cannot test
# that argument's default -- the pattern that hid the seqmaintokens bug. Once a
# default is never exercised by any generator, the only remaining check is
# whether traminex's default agrees with R's, and that has to come from R
# rather than from a table someone typed. This records R's formals so the sweep
# in tests/test_golden_default_coverage.py can make that comparison offline.
#
# Deparsed source text, not evaluated values: a default like
# `1-2*expcost` has no value until the call, and its text is what the contract
# says.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(packageVersion("TraMineR") == "2.2.14",
          packageVersion("TraMineRextras") == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

records <- list()
for (pkg in c("TraMineR", "TraMineRextras")) {
  exported <- sort(getNamespaceExports(asNamespace(pkg)))
  for (name in exported) {
    fun <- tryCatch(get(name, envir = asNamespace(pkg)), error = function(e) NULL)
    if (!is.function(fun)) next
    args <- tryCatch(formals(fun), error = function(e) NULL)
    if (is.null(args) || length(args) == 0L) next
    for (arg in names(args)) {
      default <- args[[arg]]
      text <- if (missing(default) || identical(default, quote(expr = ))) {
        "<no default>"
      } else {
        paste(deparse(default), collapse = " ")
      }
      records[[length(records) + 1L]] <- data.frame(
        package = pkg, fn = name, arg = arg, default = text,
        stringsAsFactors = FALSE
      )
    }
  }
}

formals_table <- do.call(rbind, records)
formals_table <- formals_table[order(formals_table$package, formals_table$fn,
                                     formals_table$arg), ]
write.csv(formals_table, file.path(out_dir, "r_formals.csv"), row.names = FALSE)

cat("R_FORMALS_GOLDENS_DONE\n")
cat("functions:", length(unique(formals_table$fn)),
    " formals:", nrow(formals_table), "\n")
