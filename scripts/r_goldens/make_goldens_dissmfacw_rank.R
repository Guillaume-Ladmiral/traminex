# Live-R goldens for dissmfacw on a rank-deficient design (repair of audit
# finding TX-AUD-20261010-040).
#
# dissmfacw builds every hat matrix from qr.Q(qr(X)) with all ncol(X) columns
# (dissmfacw.R:67-69) and counts degrees of freedom in columns, so a design whose
# rank is lower than its column count gains one arbitrary direction. In mvad the
# five region dummies sum to one in every row, so with the intercept the design
# below has rank 6 for 7 columns. The same model without the redundant column
# (Western), and a numeric predictor with a collinear copy (weight, 2 * weight)
# against the predictor alone, are recorded beside it: on a full-rank design R is
# right, and the column space, hence the definition's values, is the same. A
# constant numeric predictor adds no rank to the intercept at all.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_dissmfacw_rank.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data(mvad)
m <- mvad[1:40, ]
regions <- c("Belfast", "N.Eastern", "Southern", "S.Eastern", "Western")
stopifnot(all(rowSums(m[, regions] == "yes") == 1))
alph <- sort(unique(as.vector(as.matrix(mvad[, 15:86]))))
seqs <- suppressMessages(TraMineR::seqdef(m[, 15:86], alphabet = alph))
D <- suppressMessages(TraMineR::seqdist(seqs, method = "LCS"))
m$weight2 <- 2 * m$weight
m$constant <- 1

fits <- list(
  regions_full = D ~ male + Belfast + N.Eastern + Southern + S.Eastern + Western,
  regions_reduced = D ~ male + Belfast + N.Eastern + Southern + S.Eastern,
  numeric_collinear = D ~ weight + weight2,
  numeric_alone = D ~ weight,
  constant = D ~ constant
)
rows <- list()
for (name in names(fits)) {
  set.seed(1)
  r <- TraMineR::dissmfacw(fits[[name]], data = m, R = 1)
  rows[[length(rows) + 1L]] <- data.frame(case = name, Variable = as.character(r$mfac$Variable),
                                          PseudoF = sprintf("%.17g", r$mfac$PseudoF),
                                          PseudoR2 = sprintf("%.17g", r$mfac$PseudoR2),
                                          stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "dissmfacw_rank.csv"), row.names = FALSE, quote = TRUE)
write.csv(data.frame(m[, c("male", regions, "weight")]), file.path(out_dir, "dissmfacw_rank_data.csv"),
          row.names = FALSE, quote = TRUE)
write.table(format(D, digits = 17), file.path(out_dir, "dissmfacw_rank_diss.csv"), sep = ",",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
