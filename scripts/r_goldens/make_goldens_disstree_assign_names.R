# Live-R goldens for disstree.assign with R-style covariate names (repair of
# audit finding TX-AUD-20261010-006).
#
# disstree.assign assigns each covariate to a variable and evaluates every rule
# with eval(parse(text = rule)) (disstree.R:143-176), so syntactic names with
# dots (N.Eastern) and backquoted names (`my var`) work. For a non-syntactic
# covariate, R's own disstree.get.rules writes the bare name (my var <= 5) and
# records it backquoted in attr(rules, "covariates"), so disstree.assign cannot
# consume those rules: its errors are recorded in
# disstree_assign_names_errors.txt. The `leaf` rows of the fitted case are the
# rule index of each case's terminal node (disstreeleaf levels are the rule
# order), the value an assignment must reproduce.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_disstree_assign_names.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
add <- function(case, kind, values) {
  rows[[length(rows) + 1L]] <<- data.frame(
    case = case, kind = kind, position = seq_along(values),
    value = as.character(values), stringsAsFactors = FALSE
  )
}

dotted_rules <- c("N.Eastern %in% c('a')", "N.Eastern %in% c('b')")
add("manual_dotted", "rule", dotted_rules)
add("manual_dotted", "assign", TraMineR::disstree.assign(
  dotted_rules, profile = data.frame(N.Eastern = c("a", "b", "a")), covar = "N.Eastern"
))

quoted_rules <- c("`my var` %in% c('x')", "`my var` %in% c('y')")
add("manual_backquoted", "rule", quoted_rules)
add("manual_backquoted", "assign", TraMineR::disstree.assign(
  quoted_rules, profile = data.frame(`my var` = c("x", "y", "x"), check.names = FALSE), covar = "my var"
))

comparison_rules <- c("N.x <= 1.5", "N.x > 1.5")
add("manual_comparison", "rule", comparison_rules)
add("manual_comparison", "assign", TraMineR::disstree.assign(
  comparison_rules, profile = data.frame(N.x = c(1, 2, 3)), covar = "N.x"
))

# A fitted tree on a dotted and a non-syntactic covariate.
seqs <- c("AABBC", "ABBCC", "CCBAA", "BBCCA", "AAABB", "CCCBB", "AABBB", "CCCBA")
data <- as.data.frame(do.call(rbind, lapply(seqs, function(s) strsplit(s, "")[[1]])), stringsAsFactors = FALSE)
sts <- suppressMessages(TraMineR::seqdef(data, alphabet = c("A", "B", "C")))
diss <- suppressMessages(TraMineR::seqdist(sts, method = "LCS"))
pred <- data.frame(N.Eastern = c("a", "a", "b", "b", "a", "b", "a", "b"),
                   `my var` = c(1, 1, 5, 5, 2, 6, 2, 6), check.names = FALSE)
set.seed(1)
tree <- TraMineR::disstree(diss ~ N.Eastern + `my var`, data = pred, R = 1, pval = 1, min.size = 0.1)
fitted_rules <- TraMineR::disstree.get.rules(tree)
add("fitted_nonsyntactic", "rule", as.vector(fitted_rules))
add("fitted_nonsyntactic", "leaf", as.integer(TraMineR::disstreeleaf(tree, label = TRUE)))

error_text <- function(expr) tryCatch({ expr; "no error" }, error = function(e) conditionMessage(e))
writeLines(c(
  "# TraMineR 2.2-14 disstree.assign on its own rules for a non-syntactic covariate",
  paste0("covar_attribute|", error_text(TraMineR::disstree.assign(fitted_rules, profile = pred))),
  paste0("covar_bare|", gsub("\n", " ", error_text(
    TraMineR::disstree.assign(fitted_rules, profile = pred, covar = c("N.Eastern", "my var"))
  )))
), file.path(out_dir, "disstree_assign_names_errors.txt"), useBytes = TRUE)

write.csv(do.call(rbind, rows), file.path(out_dir, "disstree_assign_names.csv"),
          row.names = FALSE, quote = TRUE)
