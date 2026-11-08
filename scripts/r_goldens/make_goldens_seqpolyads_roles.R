stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for how seqpolyads compares its roles (audit finding TX-AUD-20261019-004).
#
# seqpolyads requires the roles to have the same size and identical
# alphabet(seqdata, with.missing) (seqpolyads.R:52-57); alphabet() appends the nr code
# whenever with.missing is TRUE (alphabet.R:9-11), so a role without a missing cell
# still matches one with it. Default seqdef roles whose states first appear in
# different orders share R's sorted alphabet and are computed. Two roles with the
# same states in different explicit orders are refused (identical() is false).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_seqpolyads_roles.R

suppressMessages({
  library(TraMineR)
  library(TraMineRextras)
})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

r1 <- rbind(c("a", "b", "c"), c("b", "b", "a"), c("c", "a", "a"))
r2 <- rbind(c("b", "a", "c"), c("a", "b", "b"), c("c", "c", "a"))
r3 <- rbind(c("a", "b", "c"), c("b", NA, "a"), c("c", "a", "a"))
quiet <- function(expr) suppressWarnings(suppressMessages(expr))
role <- function(m, ...) quiet(TraMineR::seqdef(m, ...))

rows <- list()
record <- function(case, expr) {
  result <- tryCatch(quiet(expr), error = function(e) e)
  if (inherits(result, "error")) {
    rows[[length(rows) + 1L]] <<- data.frame(case = case, kind = "error", position = 1L,
                                             value = conditionMessage(result))
  } else {
    d <- as.numeric(result$observed.dist)
    rows[[length(rows) + 1L]] <<- data.frame(case = case, kind = "observed.dist",
                                             position = seq_along(d), value = sprintf("%.17g", d))
  }
}

record("default_roles_first_appearance_differs",
       TraMineRextras::seqpolyads(list(role(r1), role(r2)), T = 10, method = "HAM"))
record("with_missing_true_one_role_missing",
       TraMineRextras::seqpolyads(list(role(r1), role(r3)), T = 10, method = "HAM", with.missing = TRUE))
record("with_missing_false_one_role_missing",
       TraMineRextras::seqpolyads(list(role(r1), role(r3)), T = 10, method = "HAM"))
record("explicit_orders_differ",
       TraMineRextras::seqpolyads(list(role(r1, alphabet = c("a", "b", "c")),
                                       role(r2, alphabet = c("b", "a", "c"))), T = 10, method = "HAM"))

write.csv(do.call(rbind, rows), file.path(out_dir, "seqpolyads_roles.csv"), row.names = FALSE, quote = TRUE)
