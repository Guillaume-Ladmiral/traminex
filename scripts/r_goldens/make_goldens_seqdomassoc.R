# Live-R goldens for seqdomassoc (RP4, TraMineR seqdomassoc, rep.method="overall").
#
# Cross-domain LRT / Cramer's V on the shared 3-channel multidomain fixture.
#
# UNWEIGHTED: straight from TraMineR seqdomassoc (faithful).
# WEIGHTED: traminex intentionally DIVERGES from R's seqdomassoc(weighted=TRUE).
#   R builds the weighted cross-table as xtabs(rep(weights,each=L) ~ unlist(d1) + unlist(d2)),
#   but unlist(stslist) is COLUMN-major (position-blocked) while the weight-rep is CASE-major,
#   so R assigns each cell the weight of the WRONG case (a latent R bug -- see the RP3/RP4 audit
#   report). traminex weights each case's co-occurrences by THAT case's weight. The weighted
#   golden below is therefore R's own loglin / Cramer's V applied to the CORRECTLY (per-case)
#   weighted cross-table -- a live-R computation of the sensible quantity, independent of
#   traminex's own closed-form loglin. See docs/API_PARITY.md.
#
# Requires R 4.6.2 + TraMineR 2.2-14.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_md_inputs.R")  # md_inputs(), write_phase2_csv()
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

inp <- md_inputs()
slist <- inp$s
smat <- inp$smat  # 0-based code matrices, one per domain
w <- inp$w
ndom <- length(slist)
L <- inp$L

# Unweighted: faithful seqdomassoc.
r_uw <- suppressMessages(seqdomassoc(slist, assoc = c("LRT", "V"), rep.method = "overall"))
write_phase2_csv(
  "seqdomassoc_uw", matrix(as.numeric(r_uw), nrow = nrow(r_uw)),
  "# seqdomassoc(weighted=FALSE) cols: df, LRT, p(LRT), v, p(v); rows = domain pairs\n"
)

# Weighted (sensible per-case alignment): R loglin / Cramer V on the correctly-weighted table.
sensible_row <- function(ca, cb) {
  # ca, cb: n x L 0-based code matrices. Per-case weighting: weight each case's
  # co-occurrences by that case's weight (row-major flatten aligns case-major weights).
  a <- as.vector(t(ca))
  b <- as.vector(t(cb))
  wl <- rep(w, each = L)
  xt <- xtabs(wl ~ a + b)  # no missing in this fixture -> all states present
  start <- rep(1, length(as.vector(xt)))
  start[as.vector(xt) == 0] <- 0
  nzeros <- sum(1 - start)
  ll <- suppressWarnings(loglin(xt, margin = list(1, 2), start = start, print = FALSE))
  df <- ll$df
  if (nzeros < df) df <- df - nzeros
  n <- sum(xt)
  v <- sqrt(ll$pearson / (n * (min(dim(xt)) - 1)))
  c(df, ll$lrt, pchisq(ll$lrt, df, lower.tail = FALSE), v, pchisq(ll$pearson, df, lower.tail = FALSE))
}

rows <- list()
for (d1 in 1:(ndom - 1)) {
  for (d2 in (d1 + 1):ndom) {
    rows[[length(rows) + 1]] <- sensible_row(smat[[d1]], smat[[d2]])
  }
}
r_w <- do.call(rbind, rows)
write_phase2_csv(
  "seqdomassoc_w", r_w,
  "# seqdomassoc weighted (traminex per-case alignment; R loglin on correctly-weighted table)\n"
)

cat("SEQDOMASSOC_GOLDENS_DONE\n")

## R2 complete rep.method contract -------------------------------------------
slist_w <- lapply(
  smat,
  function(m) suppressMessages(seqdef(m + 1L, alphabet = 1:3, weights = w))
)
dlist <- inp$d

r2_rows <- list()
r2_add <- function(case, contract, section, pair = "", row = 0L, col = 0L,
                   value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
    case = case,
    contract = contract,
    section = section,
    pair = pair,
    row = row,
    col = col,
    value = as.numeric(value),
    text = as.character(text),
    stringsAsFactors = FALSE
  )
}

r2_record <- function(case, contract, out, section = "result") {
  for (j in seq_len(ncol(out)))
    r2_add(case, contract, "column", col = j, text = colnames(out)[[j]])
  for (i in seq_len(nrow(out))) {
    r2_add(case, contract, "row_label", row = i, text = rownames(out)[[i]])
    for (j in seq_len(ncol(out)))
      r2_add(
        case, contract, section, pair = rownames(out)[[i]], row = i,
        col = j, value = out[i, j], text = colnames(out)[[j]]
      )
  }
  tabs <- attr(out, "cross.tables")
  if (!is.null(tabs)) {
    for (tabname in names(tabs)) {
      tab <- tabs[[tabname]]
      for (i in seq_len(nrow(tab)))
        for (j in seq_len(ncol(tab)))
          r2_add(
            case, contract, "cross_table", pair = tabname, row = i,
            col = j, value = tab[i, j],
            text = paste(rownames(tab)[[i]], colnames(tab)[[j]], sep = "|")
          )
    }
  }
}

r2_run <- function(case, contract = "R_EXACT", ...) {
  conditions <- character()
  out <- withCallingHandlers(
    seqdomassoc(...),
    message = function(m) {
      conditions <<- c(conditions, paste0("MESSAGE:", conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(wrn) {
      conditions <<- c(conditions, paste0("WARNING:", conditionMessage(wrn)))
      invokeRestart("muffleWarning")
    }
  )
  r2_record(case, contract, out)
  if (length(conditions))
    for (i in seq_along(conditions))
      r2_add(case, contract, "condition", row = i, text = conditions[[i]])
  invisible(out)
}

## Compute the selected Python divergence using R's loglin on correctly
## case-aligned weights.  This is independent of the Python implementation.
r2_assoc_from_matrices <- function(mats, case_weights, assoc = c("LRT", "V"),
                                   struct.zero = TRUE, p.value = TRUE,
                                   dnames = paste0("Dom", seq_along(mats))) {
  assoc <- toupper(assoc)
  result <- matrix(NaN, nrow = length(mats) * (length(mats) - 1L) / 2L, ncol = 5L)
  colnames(result) <- c("df", "LRT", "p(LRT)", "v", "p(v)")
  tabnames <- character()
  k <- 0L
  dn.split <- strsplit(dnames, " x ")
  for (d1 in seq_len(length(mats) - 1L)) {
    for (d2 in (d1 + 1L):length(mats)) {
      if (!any(dn.split[[d1]] %in% dn.split[[d2]])) {
        k <- k + 1L
        a <- as.vector(t(mats[[d1]]))
        b <- as.vector(t(mats[[d2]]))
        ww <- rep(case_weights, each = ncol(mats[[d1]]))
        xt <- xtabs(
          ww ~ factor(a, levels = 1:3) + factor(b, levels = 1:3)
        )
        dimnames(xt) <- list(as.character(1:3), as.character(1:3))
        start <- rep(1, length(as.vector(xt)))
        if (struct.zero) start[as.vector(xt) == 0] <- 0
        nzeros <- sum(1 - start)
        ll <- suppressWarnings(loglin(xt, margin = list(1, 2), start = start, print = FALSE))
        df <- ll$df
        if (nzeros < df) df <- df - nzeros
        result[k, "df"] <- df
        if ("LRT" %in% assoc) {
          result[k, "LRT"] <- ll$lrt
          if (p.value) result[k, "p(LRT)"] <- pchisq(ll$lrt, df, lower.tail = FALSE)
        }
        if ("V" %in% assoc) {
          result[k, "v"] <- sqrt(ll$pearson / (sum(xt) * (min(dim(xt)) - 1L)))
          if (p.value) result[k, "p(v)"] <- pchisq(ll$pearson, df, lower.tail = FALSE)
        }
        tabnames[[k]] <- paste(dnames[[d1]], dnames[[d2]], sep = "_with_")
      }
    }
  }
  keep <- which(!apply(result, 2, function(x) all(is.nan(x))))
  result <- result[seq_len(k), keep, drop = FALSE]
  rownames(result) <- tabnames
  result
}

r2_rep_matrices <- function(domains, distances, args, weights) {
  lapply(seq_along(domains), function(d) {
    callargs <- c(
      list(seqdata = domains[[d]], diss = distances[[d]], weights = weights),
      args
    )
    reps <- suppressMessages(do.call(seqrep, callargs))
    dr <- attr(reps, "Distances")
    ids <- colnames(dr)[apply(dr, 1, which.min)]
    ids <- unlist(lapply(ids, function(id) which(rownames(domains[[1]]) == id)))
    as.matrix(domains[[d]])[ids, , drop = FALSE]
  })
}

r2_run(
  "overall_uw_default_cross", seqdata.dom = slist, weighted = FALSE,
  cross.table = TRUE, dnames = c("Home", "Work", "Care")
)
r2_run(
  "overall_lrt_no_p_structfalse", seqdata.dom = slist, weighted = FALSE,
  assoc = "LRT", p.value = FALSE, struct.zero = FALSE
)
r2_run(
  "overall_v_no_p_custom_overlap", seqdata.dom = slist, weighted = FALSE,
  assoc = "V", p.value = FALSE,
  dnames = c("Life x Work", "Life", "Health")
)
r2_run(
  "overall_wrange_cross", seqdata.dom = slist, weighted = TRUE,
  wrange = c(-1L, 1L), cross.table = TRUE
)

## Missing-state routes (overall only; representative routes receive explicit
## dissimilarities and are covered separately).
missing_mats <- lapply(smat, function(m) m + 1L)
missing_mats[[1]][1, 1] <- NA
missing_mats[[2]][2, 3] <- NA
missing_mats[[3]][3, 5] <- NA
missing_domains <- lapply(
  missing_mats,
  function(m) suppressMessages(seqdef(m, alphabet = 1:3, nr = "*"))
)
r2_run(
  "overall_missing_drop", seqdata.dom = missing_domains,
  weighted = FALSE, with.missing = FALSE
)
r2_run(
  "overall_missing_keep", seqdata.dom = missing_domains,
  weighted = FALSE, with.missing = TRUE
)

rep_density_args <- list(coverage = .6, pradius = .2)
rep_dist_args <- list(criterion = "dist", nrep = 3L, pradius = .1)
r2_run(
  "rep_uw_density", seqdata.dom = slist, diss.dom = dlist,
  rep.method = "rep", weighted = FALSE, seqrep.args = rep_density_args,
  cross.table = TRUE
)
r2_run(
  "rep_uw_dist", seqdata.dom = slist, diss.dom = dlist,
  rep.method = "rep", weighted = FALSE, seqrep.args = rep_dist_args,
  assoc = "V"
)
r2_run(
  "eq_uw_prop", seqdata.dom = slist, diss.dom = dlist,
  rep.method = "eq.group", weighted = FALSE,
  seqrf.args = list(k = 3L, sortv = seq_len(inp$n)), cross.table = TRUE
)
r2_run(
  "eq_uw_first", seqdata.dom = slist, diss.dom = dlist,
  rep.method = "eq.group", weighted = FALSE,
  seqrf.args = list(k = 4L, sortv = seq_len(inp$n), grp.meth = "first")
)
r2_run(
  "eq_weighted_prop", seqdata.dom = slist_w, diss.dom = dlist,
  rep.method = "eq.group", weighted = TRUE,
  seqrf.args = list(k = 3L, sortv = seq_len(inp$n))
)

## Record upstream R's recycled-weight values as evidence, then the selected
## correctly aligned values that Python must return.
upstream_overall_w <- suppressMessages(seqdomassoc(
  slist_w, rep.method = "overall", weighted = TRUE
))
r2_record(
  "overall_weighted_correct", "SELECTED_DIVERGENCE",
  upstream_overall_w, section = "upstream_result"
)
selected_overall_w <- r2_assoc_from_matrices(
  lapply(smat, function(m) m + 1L), w
)
r2_record(
  "overall_weighted_correct", "SELECTED_DIVERGENCE", selected_overall_w
)

upstream_rep_w <- suppressMessages(seqdomassoc(
  slist_w, diss.dom = dlist, rep.method = "rep", weighted = TRUE,
  seqrep.args = rep_density_args
))
r2_record(
  "rep_weighted_correct", "SELECTED_DIVERGENCE",
  upstream_rep_w, section = "upstream_result"
)
rep_w_mats <- r2_rep_matrices(slist_w, dlist, rep_density_args, w)
selected_rep_w <- r2_assoc_from_matrices(rep_w_mats, w)
r2_record(
  "rep_weighted_correct", "SELECTED_DIVERGENCE", selected_rep_w
)

utils::write.csv(
  do.call(rbind, r2_rows),
  file.path("tests", "r_goldens", "seqdomassoc_r2_cases.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

r2_error <- function(expr) {
  tryCatch(
    {
      withCallingHandlers(
        suppressMessages(force(expr)),
        warning = function(wrn) invokeRestart("muffleWarning")
      )
      "VALUE"
    },
    error = function(e) paste0(
      "ERROR:", sub("[[:space:]]+$", "", conditionMessage(e))
    )
  )
}
short_domain <- slist[[2]][-inp$n, ]
narrow_domain <- slist[[2]][, -inp$L]
writeLines(
  c(
    paste0("not_list=", r2_error(seqdomassoc(slist[[1]]))),
    paste0("one_domain=", r2_error(seqdomassoc(list(slist[[1]])))),
    paste0("rep_missing_diss=", r2_error(seqdomassoc(slist, rep.method = "rep"))),
    paste0("eq_diss_wrong_length=", r2_error(seqdomassoc(
      slist, diss.dom = dlist[1:2], rep.method = "eq.group"
    ))),
    paste0("bad_assoc=", r2_error(seqdomassoc(slist, assoc = "tau"))),
    paste0("bad_rep_method=", r2_error(seqdomassoc(
      slist, diss.dom = dlist, rep.method = "bogus"
    ))),
    paste0("reversed_wrange=", r2_error(seqdomassoc(
      slist, wrange = c(1L, -1L)
    ))),
    paste0("incompatible_ncases=", r2_error(seqdomassoc(
      list(slist[[1]], short_domain), weighted = FALSE
    ))),
    paste0("incompatible_width=", r2_error(seqdomassoc(
      list(slist[[1]], narrow_domain), weighted = FALSE
    ))),
    paste0("seqrep_bad_args=", r2_error(seqdomassoc(
      slist, diss.dom = dlist, rep.method = "rep",
      seqrep.args = list(criterion = "bogus")
    )))
  ),
  file.path("tests", "r_goldens", "seqdomassoc_r2_errors.txt")
)

cat("SEQDOMASSOC_R2_GOLDENS_DONE\n")
