# Live-R goldens for batch 12's dissimilarity-analysis, tree and representative items (rulings R12-a to R12-r).
#
# Each row is case, field, value. Vectors are joined by spaces (numbers in %.17g, NA as NA); text vectors by ' ;; '.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_diss_trees_reps.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR)); suppressMessages(library(TraMineRextras))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14", as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.na(x), "NA", sprintf("%.17g", x)), collapse = " ")
txt <- function(x) paste(x, collapse = " ;; ")

# --- group levels as R's factor() names and groups them (R12-a: 046, 11-004). R's levels are
# unique(as.character(sorted values)): numbers as R prints them, equal to 15 significant digits one level. D: 8 points
# on a curve; each group case records dissassoc's rownames and n, disscenter's allcenter column names and medoid
# names, and dissvar.grp's names and values.
i8 <- 1:8
D8 <- sqrt(outer(i8^1.5, i8^1.5, "-")^2 + outer(i8 %% 3, i8 %% 3, "-")^2)
add("d8", "diss", num(as.vector(D8)))
groups <- list(
  g_float = c(1, 1, 1, 2, 2, 2, 10, 10),
  g_decimal = c(9, 9, 9.5, 9.5, 10, 10, 100, 100),
  g_scientific = c(1e5, 1e5, 1e5, 2, 2, 2, 3, 3),
  g_merged = c(0.3, 0.3, 0.1 + 0.2, 0.1 + 0.2, 0.7, 0.7, 0.7, 0.7),
  g_factor_numeric = factor(c(1, 1, 1, 2, 2, 2, 2, 2), levels = c(2, 1))
)
for (case in names(groups)) {
  g <- groups[[case]]
  # A factor's values are its level text read back as numbers, with its declared levels.
  add(case, "group", num(if (is.factor(g)) as.numeric(as.character(g)) else g))
  if (is.factor(g)) add(case, "levels", txt(levels(g)))
  a <- suppressMessages(dissassoc(D8, group = g, R = 1))
  add(case, "dissassoc_rows", txt(rownames(a$groups)))
  add(case, "dissassoc_n", num(a$groups$n))
  add(case, "allcenter_cols", txt(colnames(disscenter(D8, group = g, allcenter = TRUE))))
  m <- disscenter(D8, group = g, medoids.index = "first")
  add(case, "medoid_names", txt(names(m)))
  add(case, "medoids", num(m))
  v <- dissvar.grp(D8, group = g)
  add(case, "dissvar_grp_names", txt(names(v)))
  add(case, "dissvar_grp", num(v))
}
# A whole number beyond R's integer range (+-2147483647) is a double in R, so two 64-bit integers equal as doubles are
# one level, named as R writes the double (batch-12 review M-5).
g_big <- c(rep(123456789012345678, 4), rep(123456789012345679, 2), rep(5, 2))
a_big <- suppressMessages(dissassoc(D8, group = g_big, R = 1))
add("g_int64", "levels", txt(levels(factor(g_big))))
add("g_int64", "dissassoc_rows", txt(rownames(a_big$groups)))
add("g_int64", "dissassoc_n", num(a_big$groups$n))
v <- dissvar.grp(D8)
add("dissvar_grp_default", "names", txt(names(v)))
add("dissvar_grp_default", "value", num(v))

# dissassoc with group weight totals that are negative or zero (R12-f: 10-079). R computes every group whose total is
# not zero, negative totals included (a zero total gives discrepancy 0; dissassocweighted.R:313-330), and
# WeightedResidualsZ's weighted.mean is sum((x * w)[w != 0]) / sum(w), without a guard. The audit's S39 groups and
# weights on a 9-case matrix, then a group whose weights sum to 0.
i9 <- 1:9
D9 <- sqrt(outer(i9^1.5, i9^1.5, "-")^2 + outer(i9 %% 4, i9 %% 4, "-")^2)
add("d9", "diss", num(as.vector(D9)))
g9 <- rep(1:3, each = 3)
add("d9", "group", num(g9))
w9 <- list(assoc_neg_total = c(1, 1, 1, -1, -2, 0.5, 1, 1, 1), assoc_zero_total = c(1, 1, 1, 1, -0.5, -0.5, 1, 1, 1))
for (case in names(w9)) {
  a <- suppressWarnings(suppressMessages(dissassoc(D9, g9, weights = w9[[case]], R = 1, weight.permutation = "diss")))
  add(case, "weights", num(w9[[case]]))
  add(case, "groups_n", num(a$groups$n))
  add(case, "groups_disc", num(a$groups$discrepancy))
  add(case, "anova_ss", num(a$anova.table$SS))
  add(case, "anova_mse", num(a$anova.table$MSE))
  add(case, "stat", num(a$stat[, 1]))
}

# Weighted inertia in R's summation order (R12-j: 10-083). tmrWeightedInertiaDist (weightedinertia.cpp:17-34) adds
# w[j] * w[i] * D[j, i] for i outer and j > i inner (positions in the index list) to one running double, and sums the
# weights the same way; dissvar divides twice by that total (var = TRUE). Values spanning twelve orders of magnitude
# make this order differ from a pairwise sum. The subset is in a scrambled order.
i64 <- 1:64
x64 <- i64^1.37 * c(1, 1e-3, 1e3, 1e-6)[i64 %% 4 + 1]
D64 <- abs(outer(x64, x64, "-"))
w64 <- 1 + (i64 * 7) %% 11 / 13
add("d64", "diss", num(as.vector(D64)))
add("d64", "weights", num(w64))
inertia64 <- function(ind) .Call(TraMineR:::C_tmrWeightedInertiaDist, D64, 64L, 0L, as.integer(ind), as.double(w64), 0L)
sub64 <- c(40, 3, 17, 64, 1, 22, 9, 58, 33, 12)
add("inertia_order", "all", num(inertia64(i64)))
add("inertia_order", "subset", num(sub64))
add("inertia_order", "subset_value", num(inertia64(sub64)))
add("inertia_order", "dissvar", num(dissvar(D64, weights = w64)))

# disscenter's weighted mean (batch-12 review M-4): R subtracts weighted.mean(dc, w) / 2, whose sums are R's long-double
# sum(); on the 64-case input with three groups a double sum differs from it in the last bits.
add("disscenter_64_3", "group", num(rep(1:3, length.out = 64)))
add("disscenter_64_3", "values", num(disscenter(D64, group = rep(1:3, length.out = 64), weights = w64)))

# dissindic (R12-q: 11-066), TraMineRextras: a zero weight is computed (its own residual is 0/0, NaN; the other cases
# take it as R does), and a character group stops R (dissindic.R subsets data[seqok, ] on a one-column data frame,
# which drops to a bare vector that model.frame refuses for text) while the same group as a factor is computed.
indic_group <- c("a", "a", "b", "b", "b", "c", "c", "c")
indic_weights <- c(1, 0, 2, 1, 0.5, 1, 1, 3)
add("indic", "group", txt(indic_group))
add("indic", "weights", num(indic_weights))
indic_zero <- suppressMessages(dissindic(D8, group = factor(indic_group), weights = indic_weights))
add("indic_zero_weight", "marginality", num(indic_zero$marginality))
add("indic_zero_weight", "gain", num(indic_zero$gain))
indic_char <- try(suppressMessages(dissindic(D8, group = indic_group)), silent = TRUE)
add("indic_character", "result", if (inherits(indic_char, "try-error")) paste("ERROR", conditionMessage(attr(indic_char, "condition"))) else "computed")
indic_factor <- suppressMessages(dissindic(D8, group = factor(indic_group)))
add("indic_factor", "marginality", num(indic_factor$marginality))
add("indic_factor", "gain", num(indic_factor$gain))

# dissindic with a group whose weights are all zero (batch-12 review I-1). By the definition the zero-weight cases take
# no part in the centring or the projection, so the other cases' indicators are dissindic's on the subset without
# them. R gives that when the zero group is not the reference level; when it is the first level, R's qr.Q keeps a
# column for the aliased intercept and the projection gains a direction outside the design's column space.
for (case in c("indic_zero_group_ref", "indic_zero_group_other")) {
  w <- if (case == "indic_zero_group_ref") c(0, 0, 2, 1, 0.5, 1, 1, 3) else c(1, 1, 0, 0, 0, 1, 1, 3)
  full <- suppressMessages(dissindic(D8, group = factor(indic_group), weights = w))
  keep <- w > 0
  sub <- suppressMessages(dissindic(D8[keep, keep], group = factor(indic_group[keep]), weights = w[keep]))
  add(case, "weights", num(w))
  add(case, "marginality", num(full$marginality))
  add(case, "gain", num(full$gain))
  add(case, "subset_marginality", num(sub$marginality))
  add(case, "subset_gain", num(sub$gain))
}

# dissCompare's set rows and seqrep.grp's names (16 cases: two groups a/b in every set). seqrep.grp names its levels
# through TraMineR:::group, which pastes "" and the level, so they carry a leading space (seqrep-grp.R:16).
i16 <- 1:16
D16 <- sqrt(outer(i16^1.3, i16^1.3, "-")^2 + outer(i16 %% 4, i16 %% 4, "-")^2)
add("d16", "diss", num(as.vector(D16)))
cmp_group <- rep(c("a", "b"), 8)
for (case in c("set_float", "set_merged")) {
  s <- if (case == "set_float") rep(c(1, 2), each = 8) else rep(c(0.3, 0.1 + 0.2), each = 8)
  add(case, "set", num(s))
  r <- suppressMessages(dissCompare(D16, group = cmp_group, set = s, s = 0))
  add(case, "dissCompare_rows", txt(rownames(r)))
}
seq16 <- suppressMessages(seqdef(t(sapply(i16, function(k) c("A", "B", "C")[c(k %% 3, (k + 1) %% 3, k %% 2) + 1])),
                                 alphabet = c("A", "B", "C")))
add("seq16", "cells", txt(apply(as.matrix(seq16), 1, paste, collapse = "")))
for (case in c("rep_float", "rep_text")) {
  g <- if (case == "rep_float") rep(c(1, 2), each = 8) else rep(c("x", "y"), each = 8)
  r <- suppressMessages(seqrep.grp(seq16, group = g, mdis = D16))
  add(case, "seqrep_grp_names", txt(names(r)))
}

# --- tree predictor levels (R12-p: 11-055). DTNGroupFactorBinary factorises the predictor (disstree.R:404), so 0.3 and
# 0.1 + 0.2 are one level and the tree does not split; a numeric split's break is as.numeric of the level text
# (disstree.R:544), so a predictor at 1/3 splits at 0.333333333333333.
y20 <- c(rep(0, 10), rep(5, 10)) + seq(0, 0.019, length.out = 20)
D20 <- as.matrix(dist(y20))
add("d20", "diss", num(as.vector(D20)))
x_merged <- c(rep(0.3, 10), rep(0.1 + 0.2, 10))
add("tree_merged", "x", num(x_merged))
t_merged <- disstree(D20 ~ x, data = data.frame(x = x_merged), R = 1, min.size = 2)
add("tree_merged", "split", is.null(t_merged$root$split) == FALSE)
x_third <- c(rep(1 / 3, 10), rep(2 / 3, 10))
add("tree_third", "x", num(x_third))
t_third <- disstree(D20 ~ x, data = data.frame(x = x_third), R = 1, min.size = 2)
add("tree_third", "breaks", num(t_third$root$split$breaks))
add("tree_third", "rules", txt(disstree.get.rules(t_third)))
# A root-only tree (R12-b: 047). With max.depth = 1 there is no split: disstreeleaf(label = TRUE) is NA for every case
# (the factor has no levels), disstree.get.rules turns the empty levels into one NA rule through its 1:length loop
# (disstree.R:117-137), and disstree.assign evaluates that NA rule to NA for every profile (disstree.R:166-173).
t_root <- disstree(D20 ~ x, data = data.frame(x = x_third), R = 1, min.size = 2, max.depth = 1)
add("tree_root", "leaf_labels_na", paste(is.na(disstreeleaf(t_root, label = TRUE)), collapse = " "))
r_root <- disstree.get.rules(t_root)
add("tree_root", "rules_na", paste(is.na(r_root), collapse = " "))
add("tree_root", "rules_covariates", txt(attr(r_root, "covariates")))
add("tree_root", "assign", num(disstree.assign(r_root, data.frame(x = x_third[c(1, 11, 20)]))))
# Rule terms are R expressions (R12-e: 10-065): disstree.assign parses each rule (disstree.R:166-167), so a quoted level
# holding ' & ' stays one string. A level with a quote makes rules R cannot parse (R's assign stops); the same level
# written with an R escape or in double quotes assigns; a doubled quote '' is not an R string and stops.
D6 <- as.matrix(dist(c(0, 0.1, 0.2, 5, 5.1, 5.2)))
add("d6", "diss", num(as.vector(D6)))
assign_or_error <- function(rules, profile) {
  a <- try(disstree.assign(rules, profile = profile), silent = TRUE)
  if (inherits(a, "try-error")) "ERROR" else num(a)
}
for (case in c("rule_amp", "rule_quote")) {
  g <- rep(c(if (case == "rule_amp") "u & v" else "o'k & z", "w"), each = 3)
  add(case, "g", txt(g))
  r <- disstree.get.rules(disstree(D6 ~ g, data = data.frame(g = g), R = 1, min.size = 2))
  add(case, "rules", txt(r))
  add(case, "assign", assign_or_error(r, data.frame(g = g)))
}
g_quote <- rep(c("o'k & z", "w"), each = 3)
hand <- list(rule_escaped = "g %in% c('o\\'k & z')", rule_double = "g %in% c(\"o'k & z\")",
             rule_doubled = "g %in% c('o''k & z')")
for (case in names(hand)) {
  r <- structure(c(hand[[case]], "g %in% c('w')"), covariates = "g")
  add(case, "rules", txt(r))
  add(case, "assign", assign_or_error(r, data.frame(g = g_quote)))
}
# Medoid ties in R's order (R12-k: 11-023). tmrWeightedInertiaContrib (weightedinertia.cpp:68-101) adds each case's
# terms in index order, so identical sequences tie exactly and which.min takes the first. The weighted mvad tree (OM,
# constant costs; every node's medoid) and dissrf with k = 10 on an index-order sortv, weighted and not (R's groups list
# cases in sort order, dissrf.R:141, and on the default MDS axis the groups are mirrored: TX-REP-20261017-006). mvad is
# stored as alphabet positions, one string per case; the predictors keep their factor levels; three sums check the
# distance matrix the test computes.
data(mvad)
alph_mvad <- sort(unique(as.vector(as.matrix(mvad[, 17:86]))))
add("mvad", "alphabet", txt(alph_mvad))
add("mvad", "codes", paste(apply(as.matrix(mvad[, 17:86]), 1, function(r) paste(match(r, alph_mvad), collapse = "")),
                           collapse = " "))
add("mvad", "weights", num(mvad$weight))
mvad_pred <- c("male", "funemp", "gcse5eq", "Grammar", "fmpr", "livboth")
for (v in mvad_pred) {
  add("mvad", v, txt(as.character(mvad[[v]])))
  add("mvad", paste0(v, "_levels"), txt(levels(mvad[[v]])))
}
seq_mvad <- suppressMessages(seqdef(mvad[, 17:86], alphabet = alph_mvad, weights = mvad$weight))
D_mvad <- suppressMessages(seqdist(seq_mvad, method = "OM", sm = "CONSTANT", indel = 1))
add("mvad", "diss_sums", num(c(sum(D_mvad), sum(D_mvad^2), sum(D_mvad * outer(1:712, 1:712)))))
t_mvad <- disstree(D_mvad ~ male + funemp + gcse5eq + Grammar + fmpr + livboth, data = mvad[mvad_pred],
                   weights = mvad$weight, R = 1, weight.permutation = "diss", min.size = 0.05, max.depth = 5)
node_names <- character(0); node_medoids <- integer(0)
walk_nodes <- function(node, name) {
  node_names <<- c(node_names, name); node_medoids <<- c(node_medoids, node$info$medoid)
  if (!is.null(node$split)) { walk_nodes(node$kids[[1]], paste0(name, "L")); walk_nodes(node$kids[[2]], paste0(name, "R")) }
}
walk_nodes(t_mvad$root, "N")
add("tree_medoids", "names", txt(node_names))
add("tree_medoids", "medoids", num(node_medoids))
for (case in c("rf_weighted", "rf_unweighted")) {
  rf <- suppressMessages(dissrf(D_mvad, k = 10, weights = if (case == "rf_weighted") mvad$weight, sortv = 1:712))
  add(case, "medoids", num(rf$medoids))
  add(case, "R2", num(rf$R2))
}

# dissrf's default MDS axis does not depend on the scale (R12-o: 11-052): R's cmdscale has no tie test. mvad cases 1 to
# 80, LCS, k = 8, at scales 1, 1e-4 and 1e-10 (two sums check the distance matrix the test computes).
D80 <- suppressMessages(seqdist(seq_mvad[1:80, ], method = "LCS"))
add("rf_scale", "diss_sums", num(c(sum(D80), sum(D80^2))))
scales80 <- c(1, 1e-4, 1e-10)
add("rf_scale", "scales", num(scales80))
add("rf_scale", "R2", num(sapply(scales80, function(s) suppressMessages(dissrf(D80 * s, k = 8))$R2)))

# dissrf orders a text or factor sortv with order(sortv) (R12-r: TX-REP-20261021-001; dissrf.R:94): text in R's
# collation, a factor by its level codes, NA last. The mvad 1:80 LCS matrix, k = 4; each sortv also with two NA.
sortv_text <- rep(c("c", "a", "B", "b"), 20)
sortv_na <- c(5, 30)
add("rf_sortv", "text", txt(sortv_text))
add("rf_sortv", "levels", txt(c("c", "b", "B", "a")))
add("rf_sortv", "na", num(sortv_na))
for (case in c("rf_sortv_text", "rf_sortv_factor", "rf_sortv_text_na", "rf_sortv_factor_na")) {
  s <- sortv_text
  if (grepl("_na$", case)) s[sortv_na] <- NA
  if (grepl("factor", case)) s <- factor(s, levels = c("c", "b", "B", "a"))
  rf <- suppressMessages(dissrf(D80, k = 4, sortv = s))
  add(case, "medoids", num(rf$medoids))
  add(case, "R2", num(rf$R2))
}

# The controls of dissrf's former near-tie tolerance: R takes the exact first minimum. [0, s, 2s] at four scales (the
# middle case is the medoid at every scale), and a 4-case matrix whose d(1, 2) is 1 + 4e-15 (the third case's
# contribution, 3, is the minimum; the second's is 3 + 4e-15).
base3 <- abs(outer(0:2, 0:2, "-"))
scales3 <- c(1, 1e-10, 1e-14, 1e-20)
add("rf_near_tie", "scales", num(scales3))
add("rf_near_tie", "medoids", num(sapply(scales3, function(s)
  suppressMessages(dissrf(base3 * s, k = 1, sortv = NULL, weights = rep(1, 3)))$medoids)))
perturbed4 <- matrix(c(0, 1, 1, 2, 1, 0, 1, 1, 1, 1, 0, 1, 2, 1, 1, 0), 4)
perturbed4[1, 2] <- perturbed4[2, 1] <- 1 + 4e-15
add("medoid_last_ulp", "diss", num(as.vector(perturbed4)))
add("medoid_last_ulp", "medoid", num(disscenter(perturbed4, medoids.index = "first", weights = rep(1, 4))))

# seqtree passes seqdist.args as given (R12-g: 10-080; seqtree.R:3, 15-16): the default list(method = "LCS",
# norm = "auto") builds the tree, and list() calls seqdist without a method, so R stops.
seqtree_cells <- c("AABBC", "ABBCC", "CCBAA", "BBCCA", "AAABB", "CCCBB")
seqtree_g <- c("u", "u", "v", "v", "u", "v")
add("seqtree_args", "cells", txt(seqtree_cells))
add("seqtree_args", "g", txt(seqtree_g))
seqtree_seq <- suppressMessages(seqdef(t(sapply(seqtree_cells, function(x) strsplit(x, "")[[1]])), alphabet = c("A", "B", "C")))
seqtree_default <- suppressWarnings(suppressMessages(seqtree(seqtree_seq ~ g, data = data.frame(g = seqtree_g), R = 10)))
add("seqtree_args", "default_root_n", num(seqtree_default$root$info$n))
seqtree_empty <- try(suppressWarnings(suppressMessages(seqtree(seqtree_seq ~ g, data = data.frame(g = seqtree_g), R = 10,
                                                               seqdist.args = list()))), silent = TRUE)
add("seqtree_args", "empty", if (inherits(seqtree_empty, "try-error")) paste("ERROR", conditionMessage(attr(seqtree_empty, "condition"))) else "tree")

# Escapes in a rule's levels (batch-12 review M-3): disstree.assign parses each rule, so a level is read with R's string
# escapes (\n, \t, octal, \x, \u, \U ...) and an unrecognised escape, a missing hex digit or a nul stops R. Each case
# records the literal as written in the rule, the code points R reads (or ERROR) and disstree.assign's result on a
# profile holding that text and "w".
escape_literals <- c("a\\nb", "a\\tb", "a\\\\b", "\\x41", "\\101", "\\u00e9", "\\u{e9}", "\\U0001F600", "\\U{1F600}",
                     "a\\\"b", "a\\`b", "a\\ b", "\\ab", "\\x4g", "\\1234", "a\\'b", "a\\qb", "a\\0b", "\\x", "\\8",
                     "\\777", "\\x00", "\\u", "\\u{e9")
for (k in seq_along(escape_literals)) {
  lit <- escape_literals[k]
  case <- paste0("rule_escape_", k)
  add(case, "literal", lit)
  value <- try(eval(parse(text = paste0("'", lit, "'"))), silent = TRUE)
  rules <- structure(c(paste0("g %in% c('", lit, "')"), "g %in% c('w')"), covariates = "g")
  if (inherits(value, "try-error")) {
    add(case, "codes", "ERROR")
    add(case, "assign", assign_or_error(rules, data.frame(g = c("x", "w"))))
  } else {
    add(case, "codes", num(utf8ToInt(value)))
    add(case, "assign", assign_or_error(rules, data.frame(g = c(value, "w"))))
  }
}

# Plot panels and legends follow factor(group) (seqedplot.R:53-66, seqplot-tentrop.R:10): one level for 0.3 and 0.1 + 0.2.
plot_group <- c(1, 1, 0.3, 0.1 + 0.2)
add("plot_levels", "group", num(plot_group))
add("plot_levels", "levels", txt(levels(factor(plot_group))))

# dissmergegroups numbers the groups with as.integer(factor(group)) (dissmergegroups.R:105), so 0.3 and 0.1 + 0.2 are one
# group from the start.
merge_group <- rep(c(0.3, 0.1 + 0.2, 0.7, 0.9), each = 4)
add("merge_merged", "group", num(merge_group))
add("merge_merged", "membership", num(suppressMessages(dissmergegroups(D16, group = merge_group, silent = TRUE))))

write.csv(do.call(rbind, rows), file.path(out_dir, "diss_trees_reps.csv"), row.names = FALSE, quote = TRUE)
