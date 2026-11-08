# Live-R goldens for seqdomassoc's windowed cross-tables (repair of audit
# finding TX-AUD-20261011-009).
#
# With wrange, seqdomassoc hands seqxtabs.win the column-major unlist of each
# domain (all cases at time 1, then time 2, ...) and windows it in blocks of
# ncol/ncases as if each block were one case's sequence (seqdomassoc.R:170,210,
# 277-301), so it counts co-occurrences across different people, ignores the
# weights and uses the alphabet only even with with.missing = TRUE. The
# definition counts, for each person, the state at position p in the first domain
# against the states at p + wrange[1] .. p + wrange[2] in the second. Its values
# come from R too: seqdomassoc's position-wise route (right for unweighted data)
# on sequences that list, for each person, every (p, window position) pair.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_domassoc_window.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_md_inputs.R"))

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The (position, window position) pairs of a sequence of length L.
window_pairs <- function(L, wrange) {
  do.call(rbind, lapply(seq_len(L), function(p) {
    w <- max(1, p + wrange[1]):min(L, p + wrange[2])
    cbind(p = rep(p, length(w)), w = w)
  }))
}
expand <- function(seqs, columns, alphabet) {
  cells <- as.matrix(seqs)
  cells[cells %in% c(attr(seqs, "nr"), attr(seqs, "void"))] <- NA
  suppressMessages(TraMineR::seqdef(as.data.frame(cells[, columns, drop = FALSE]), alphabet = alphabet,
                                    right = NA))
}
rows <- list()
record <- function(case, pair, kind, result) {
  stats <- unclass(result)
  for (name in colnames(stats)) {
    rows[[length(rows) + 1L]] <<- data.frame(case = case, pair = pair, kind = kind, section = "stat",
                                             name = name, value = sprintf("%.17g", stats[1, name]),
                                             stringsAsFactors = FALSE)
  }
  table <- attr(result, "cross.tables")[[1]]
  for (i in seq_len(nrow(table))) for (j in seq_len(ncol(table))) {
    rows[[length(rows) + 1L]] <<- data.frame(case = case, pair = pair, kind = kind, section = "table",
                                             name = paste(rownames(table)[i], colnames(table)[j], sep = "|"),
                                             value = sprintf("%.17g", table[i, j]), stringsAsFactors = FALSE)
  }
}
run <- function(case, domains, alphabets, wrange, with.missing = FALSE) {
  L <- ncol(domains[[1]])
  pairs <- window_pairs(L, wrange)
  for (a in 1:(length(domains) - 1)) for (b in (a + 1):length(domains)) {
    pair <- paste0("Dom", a, "_with_Dom", b)
    windowed <- suppressMessages(TraMineR::seqdomassoc(list(domains[[a]], domains[[b]]), wrange = wrange,
                                                       cross.table = TRUE, with.missing = with.missing,
                                                       dnames = c(paste0("Dom", a), paste0("Dom", b))))
    record(case, pair, "r_window", windowed)
    expanded <- list(expand(domains[[a]], pairs[, "p"], alphabets[[a]]),
                     expand(domains[[b]], pairs[, "w"], alphabets[[b]]))
    definition <- suppressMessages(TraMineR::seqdomassoc(expanded, cross.table = TRUE,
                                                         with.missing = with.missing,
                                                         dnames = c(paste0("Dom", a), paste0("Dom", b))))
    record(case, pair, "definition", definition)
  }
}

# The reproduced two-person case.
d1 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B"), b = c("A", "B"), c = c("A", "B")),
                                        alphabet = c("A", "B")))
d2 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("X", "Y"), b = c("X", "Y"), c = c("X", "Y")),
                                        alphabet = c("X", "Y")))
run("tiny", list(d1, d2), list(c("A", "B"), c("X", "Y")), c(0L, 1L))
# The domains of make_goldens_seqdomassoc.R's overall_wrange_cross case.
inp <- md_inputs()
run("md_wrange", inp$s, rep(list(as.character(1:3)), 3), c(-1L, 1L))
# Both domains with missing cells (leading or internal, so missing, not void),
# with.missing = TRUE.
m1 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("A", "B", NA), b = c("A", NA, "B"), c = c("B", "B", "A")),
                                        alphabet = c("A", "B")))
m2 <- suppressMessages(TraMineR::seqdef(data.frame(a = c("X", NA, "Y"), b = c("X", "Y", NA), c = c("Y", "X", "Y")),
                                        alphabet = c("X", "Y")))
stopifnot(!any(as.matrix(m1) == attr(m1, "void")), !any(as.matrix(m2) == attr(m2, "void")))
run("missing", list(m1, m2), list(c("A", "B"), c("X", "Y")), c(0L, 1L), with.missing = TRUE)
write.csv(do.call(rbind, rows), file.path(out_dir, "domassoc_window.csv"), row.names = FALSE, quote = TRUE)
