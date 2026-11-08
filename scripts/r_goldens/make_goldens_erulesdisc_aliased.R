# Live-R goldens for seqerulesdisc when a rule's predictor is constant (repair of
# audit finding TX-AUD-20261011-010).
#
# When the predictor event of a rule is constant in the data of its model, glm
# aliases its coefficient (NA) and summary(model)$coefficient[2, ] is the next
# row, AGE's, which seqerulesdisc reports under the predictor's name
# (seqerulesdisc.R:81-87). Without that misreading the coefficient has no
# estimate, its p-value is missing, and seqerulesdisc's own filter leaves such a
# rule out (!is.na(anyp), :96-98). The generator records R's rules for a control
# fixture and for the same fixture with event A at the minimum age for everyone,
# and, for each two-event rule, whether R's glm aliases the predictor
# coefficient (the same glm call seqerulesdisc makes).
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_erulesdisc_aliased.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages({library(TraMineR); library(TraMineRextras)})
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14",
          as.character(packageVersion("TraMineRextras")) == "0.6.9")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ids <- 1:48
fixture <- function(constant_a) data.frame(IDPERS = ids,
  A = if (constant_a) rep(18, 48) else 19 + ids %% 4,
  AST = if (constant_a) rep(1L, 48) else as.integer(ids %% 4 != 0),
  B = 21 + ids %% 4, BST = as.integer(ids %% 5 != 0), C = 23 + ids %% 3, CST = as.integer(ids %% 3 != 0))
rows <- list()
for (case in c("control", "constant_a")) {
  d <- fixture(case == "constant_a")
  disc <- suppressMessages(TraMineRextras::createdatadiscrete(ids, d, c("A", "B", "C"),
                                                              data.frame(IDPERS = ids, AGE = 18),
                                                              data.frame(IDPERS = ids, AGE = 26)))
  events <- do.call(rbind, lapply(ids, function(i) do.call(rbind, lapply(c("A", "B", "C"), function(e)
    if (d[i, paste0(e, "ST")] == 1) data.frame(id = i, timestamp = d[i, e], event = e) else NULL))))
  eseq <- suppressMessages(TraMineR::seqecreate(events))
  subseq <- suppressMessages(TraMineR::seqefsub(eseq, str.subseq = c("(A)-(B)", "(B)-(C)", "(A)-(C)")))
  rules <- suppressMessages(TraMineRextras::seqerulesdisc(subseq, disc, d, pvalue = 1, adjust = FALSE))
  for (rule in as.character(rules$Hazards)) {
    rows[[length(rows) + 1L]] <- data.frame(case = case, kind = "r_rule", key = rule, value = "",
                                            stringsAsFactors = FALSE)
  }
  for (pair in list(c("A", "B"), c("B", "C"), c("A", "C"))) {
    data.cur <- disc$data[[pair[2]]]
    fit <- suppressWarnings(glm(data.cur[, 3] ~ data.cur[, pair[1]] + data.cur[, "AGE"],
                                family = binomial(link = "cloglog")))
    rows[[length(rows) + 1L]] <- data.frame(case = case, kind = "aliased", key = paste(pair, collapse = "=>"),
                                            value = as.character(is.na(coef(fit)[2])), stringsAsFactors = FALSE)
  }
}
write.csv(do.call(rbind, rows), file.path(out_dir, "erulesdisc_aliased.csv"), row.names = FALSE, quote = TRUE)
