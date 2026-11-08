stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
# Live-R goldens for batch 18's ledger items (TX-AUD-20261011-005, rulings R18-j): what TraMineR returns for the inputs
# the tests use, where the ledger rows said EXACT although traminex departs from R, and for gower_matrix's weight
# totals of zero and below, where traminex now returns R's values.
#
# Each row is case, field, value. Vectors are joined by spaces (numbers in %.17g, NaN as NaN, NA as NA); text vectors
# by ' ;; '.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_ledger_items.R
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.nan(x), "NaN", ifelse(is.na(x), "NA", sprintf("%.17g", x))), collapse = " ")
txt <- function(x) paste(ifelse(is.na(x), "NA", x), collapse = " ;; ")
err <- function(expr) tryCatch({suppressWarnings(suppressMessages(expr)); "computed"},
                               error = function(e) paste("ERROR:", conditionMessage(e)))
quiet <- function(expr) suppressMessages(suppressWarnings(expr))

cells <- rbind(c("A", "A", "B", "C"), c("A", "B", "C", "C"), c("B", "C", "C", "A"), c("C", "C", "A", "B"))
alphabet <- c("A", "B", "C")

# --- is.stslist (TraMineR-is_helpers.R:11-19). It requires the weights' names to equal the row names; seqdef with id=
# names the rows by id but the weights 1..n, so R rejects an object its own seqdef built. Recorded: the names and the
# results with id= and weights=, weights= alone and id= alone.
s_idw <- quiet(seqdef(cells, alphabet = alphabet, id = c("p1", "p2", "p3", "p4"), weights = c(1, 2, 3, 4)))
add("is_stslist", "cells", txt(apply(cells, 1, paste, collapse = "-")))
add("is_stslist", "row_names", txt(rownames(s_idw)))
add("is_stslist", "weight_names", txt(names(attr(s_idw, "weights"))))
add("is_stslist", "id_and_weights", as.character(is.stslist(s_idw)))
add("is_stslist", "weights_only", as.character(is.stslist(quiet(seqdef(cells, alphabet = alphabet, weights = c(1, 2, 3, 4))))))
add("is_stslist", "id_only", as.character(is.stslist(quiet(seqdef(cells, alphabet = alphabet, id = c("p1", "p2", "p3", "p4"))))))

# --- seqLLCS and seqLLCP with an interior void (seqLLCS.R, seqLLCP.R). Both pass every cell to C with seqlength() as
# the length, which assumes the voids trail; seqrecode can put one inside. Recorded: the recoded cells and R's 4 x 4
# matrices of seqLLCS(s[i, ], s[j, ]) and seqLLCP(s[i, ], s[j, ]), row by row.
s_iv <- quiet(seqrecode(quiet(seqdef(cells, alphabet = alphabet)), recodes = list("%" = "B")))
add("llc_interior_void", "cells", txt(apply(as.matrix(s_iv), 1, paste, collapse = "-")))
add("llc_interior_void", "seqlength", num(as.vector(seqlength(s_iv))))
pairs <- function(s, fun) as.vector(t(outer(1:4, 1:4, Vectorize(function(i, j) fun(s[i, ], s[j, ])))))
add("llc_interior_void", "seqLLCS", num(pairs(s_iv, seqLLCS)))
add("llc_interior_void", "seqLLCP", num(pairs(s_iv, seqLLCP)))

# The same with a missing cell beside each interior void (batch-18 review M-3): seqdef makes NA the missing state "*",
# which both R's C code and traminex compare as a state. Recorded: the input, the recoded cells, seqlength() and the
# two matrices, row by row.
cells_nr <- rbind(c("A", NA, "B", "C"), c("B", "A", NA, "C"), c(NA, "B", "C", "A"), c("C", "B", NA, "A"))
s_nr <- quiet(seqrecode(quiet(seqdef(cells_nr, alphabet = alphabet)), recodes = list("%" = "B")))
add("llc_void_beside_missing", "input", txt(apply(cells_nr, 1, paste, collapse = "-")))
add("llc_void_beside_missing", "cells", txt(apply(as.matrix(s_nr), 1, paste, collapse = "-")))
add("llc_void_beside_missing", "seqlength", num(as.vector(seqlength(s_nr))))
add("llc_void_beside_missing", "seqLLCS", num(pairs(s_nr, seqLLCS)))
add("llc_void_beside_missing", "seqLLCP", num(pairs(s_nr, seqLLCP)))

# --- seqecontain on a single event sequence (seqecontain.R:6-7): R stops; traminex returns one value (owner decision
# 2026-10-24: kept and recorded). Recorded: the events, R's error and R's values on the sequence list.
es <- seqecreate(id = c(1, 1, 2), timestamp = c(0, 1, 0), event = c("a", "b", "b"))
add("seqecontain_single", "events", txt(c("1:0:a", "1:1:b", "2:0:b")))
add("seqecontain_single", "single", err(seqecontain(es[[1]], "a")))
add("seqecontain_single", "list", txt(as.character(seqecontain(es, "a"))))

# --- gower_matrix with weight totals of zero and below (dissmfacw.R:11-22). s = weights / sum(weights): a zero total
# gives NaN, a negative total finite values (the shares still sum to 1). Recorded: D and the two matrices, row by row.
D <- matrix(c(0, 1, 2, 1, 0, 1, 2, 1, 0), 3, 3)
add("gower_totals", "D", num(as.vector(t(D))))
add("gower_totals", "zero_total", num(as.vector(t(TraMineR:::gower_matrix(D, weights = c(1, 1, -2))))))
add("gower_totals", "negative_total", num(as.vector(t(TraMineR:::gower_matrix(D, weights = c(-1, -1, -2))))))

write.csv(do.call(rbind, rows), file.path(out_dir, "ledger_items.csv"), row.names = FALSE, quote = TRUE)
