# Live-R goldens for R's rank rule (TX-REP-20261019-002, batch-11 ruling R11-m).
#
# R's rank is LINPACK dqrdc2's (R's src/appl/dqrdc2.f): a column whose remaining norm falls below tol times its own
# original norm is moved to the end. qr() uses tol = 1e-7, and dissmfacw builds its projections with qr() of the
# weighted design (dissmfacw.R:66-70). glm.fit calls Cdqrls with tol = min(1e-7, epsilon / 1000) = 1e-11 on the design
# weighted by the working weights at every iteration and aliases the columns past the rank (NA coefficients).
#
# Each row is case, field, value. Matrices are rows joined by '|', cells by ',' (%.17g); vectors are joined by spaces,
# with NA for an aliased coefficient.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_rank_and_convergence.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
rows <- list()
add <- function(case, field, value) rows[[length(rows) + 1L]] <<- data.frame(case = case, field = field, value = value)
num <- function(x) paste(ifelse(is.na(x), "NA", sprintf("%.17g", x)), collapse = " ")
mat <- function(X) paste(apply(X, 1, function(r) paste(sprintf("%.17g", r), collapse = ",")), collapse = "|")

# --- qr() rank and pivot (tol 1e-7) on designs of three columns near rank deficiency or of very different scales,
# and on designs that move a middle column, hold a zero or copied column, or have more columns than rows.
set.seed(20261022)
n <- 30
x1 <- rnorm(n); z <- rnorm(n); w <- rnorm(n)
designs <- list()
for (k in 3:10) designs[[sprintf("near_collinear_1e-%d", k)]] <- cbind(1, x1, x1 + 10^-k * z)
for (k in c(3, 6, 9, 12)) {
  designs[[sprintf("small_scale_1e-%d", k)]] <- cbind(1, x1, 10^-k * w)
  designs[[sprintf("large_scale_1e%d", k)]] <- cbind(1, 10^k * x1, w)
}
designs[["dependent_middle"]] <- cbind(1, x1, 2 * x1 + 1e-9 * z, w)
designs[["copied_column"]] <- cbind(1, x1, 3 * x1)
designs[["zero_column"]] <- cbind(1, 0, x1)
designs[["wide_3x5"]] <- matrix(rnorm(15), 3, 5)
for (case in names(designs)) {
  X <- designs[[case]]
  q <- qr(X)
  add(case, "design", mat(X))
  add(case, "qr_rank", q$rank)
  add(case, "qr_pivot", paste(q$pivot, collapse = " "))
}

# --- dissmfacw on a design with a column scaled by 1e9, unweighted and weighted: R's qr() keeps all three columns.
set.seed(20261023)
pts <- matrix(rnorm(2 * n), n)
mfac_d <- as.matrix(dist(pts))
mfac_data <- data.frame(a = 1e9 * rnorm(n), b = rnorm(n))
mfac_w <- runif(n, 0.5, 2)
add("mfac_scaled", "diss", mat(mfac_d))
add("mfac_scaled", "a", num(mfac_data$a))
add("mfac_scaled", "b", num(mfac_data$b))
add("mfac_scaled", "weights", num(mfac_w))
for (weighted in c(FALSE, TRUE)) {
  case <- if (weighted) "mfac_scaled_weighted" else "mfac_scaled"
  r <- dissmfacw(mfac_d ~ a + b, data = mfac_data, R = 1, weights = if (weighted) mfac_w else NULL)
  add(case, "variables", paste(r$mfac$Variable, collapse = " "))
  add(case, "pseudo_f", num(r$mfac$PseudoF))
  add(case, "pseudo_r2", num(r$mfac$PseudoR2))
}

# --- glm aliasing (Cdqrls tol 1e-11 on the weighted design): the third column of cbind(1, x1, x1 + eps * z) is
# aliased at eps = 1e-12 and kept at eps = 1e-10. The response is the one the reproduction used.
set.seed(1)
y <- rbinom(n, 1, 0.5)
add("glm_response", "y", paste(y, collapse = " "))
for (k in c(10, 12)) {
  X <- cbind(1, x1, x1 + 10^-k * z)
  add(sprintf("glm_near_collinear_1e-%d", k), "design", mat(X))
  for (link in c("logit", "cloglog")) {
    case <- sprintf("glm_near_collinear_1e-%d_%s", k, link)
    nc_warns <- character(0)
    fit <- withCallingHandlers(glm(y ~ X - 1, family = binomial(link = link)),
                               warning = function(w) { nc_warns <<- c(nc_warns, conditionMessage(w)); invokeRestart("muffleWarning") })
    cf <- unname(coef(fit))
    p <- rep(NA_real_, length(cf))
    p[!is.na(cf)] <- summary(fit)$coefficients[, 4]
    add(case, "coefficients", num(cf))
    add(case, "pvalues", num(p))
    # R's warnings (TX-REP-20261022-005): none on these designs.
    add(case, "warnings", paste(nc_warns, collapse = " | "))
  }
}

# --- glm without convergence (TX-AUD-20261019-085): on these separated data glm.fit stops after maxit = 25
# iterations, warns and returns the last iterate. The data are the reproduction's (seed 1, n = 21).
set.seed(1)
n_sep <- 20 + 1 %% 40
x_sep <- rnorm(n_sep)
y_sep <- as.integer(x_sep + rnorm(n_sep, sd = 0.05 * (1 %% 3)) > 0)
add("glm_separated", "x", num(x_sep))
add("glm_separated", "y", paste(y_sep, collapse = " "))
for (link in c("logit", "cloglog")) {
  case <- sprintf("glm_separated_%s", link)
  warns <- character(0)
  fit <- withCallingHandlers(glm(y_sep ~ x_sep, family = binomial(link = link)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  add(case, "converged", fit$converged)
  add(case, "iterations", fit$iter)
  add(case, "warnings", paste(warns, collapse = " | "))
  add(case, "coefficients", num(unname(coef(fit))))
}

# --- glm on separated data beyond |eta| = 30 (TX-REP-20261022-004): R's logit linkinv and mu.eta (C code) give
# mu = eps/(1+eps) or 1/(1+eps) and mu.eta = eps outside [-30, 30], which decides the working weights there. Two of
# the batch-11 review's data sets: gap (converges) and rand30 (does not), y = (x > 0), both links.
glm_fit_case <- function(case, x, y, link) {
  warns <- character(0)
  fit <- withCallingHandlers(glm(y ~ x, family = binomial(link = link)),
                             warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") })
  add(case, "converged", fit$converged)
  add(case, "iterations", fit$iter)
  add(case, "warnings", paste(warns, collapse = " | "))
  add(case, "coefficients", num(unname(coef(fit))))
}
sep_sets <- list(gap = c(-2, -1.5, -1, 1, 1.5, 2, 5, 9))
set.seed(5)
sep_sets$rand30 <- sort(rnorm(30))
for (nm in names(sep_sets)) {
  x_set <- sep_sets[[nm]]
  y_set <- as.integer(x_set > 0)
  add(sprintf("glm_%s", nm), "x", num(x_set))
  add(sprintf("glm_%s", nm), "y", paste(y_set, collapse = " "))
  for (link in c("logit", "cloglog")) glm_fit_case(sprintf("glm_%s_%s", nm, link), x_set, y_set, link)
}

write.csv(do.call(rbind, rows), file.path(out_dir, "rank_and_convergence.csv"), row.names = FALSE, quote = TRUE)
