# Live-R goldens for dissdomassoc (RP4 T14).
#
# Weighted inter-domain Pearson/Spearman correlations + R2 + Cronbach on the
# shared 3-channel multidomain fixture.
#
# Requires R 4.6.2 + TraMineR 2.2-14; see scripts/r_goldens/README.md.

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source("scripts/r_goldens/make_goldens_md_inputs.R")  # md_inputs(), write_phase2_csv()
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

inp <- md_inputs()
dlist <- inp$d
w <- inp$w

gen <- function(tag, weights) {
  r <- suppressMessages(dissdomassoc(
    dlist, what = c("pearson", "spearman", "R2", "cronbach"), weights = weights
  ))
  write_phase2_csv(paste0("dissdomassoc_", tag, "_pearson"), r$Pearson, "# Pearson corr\n")
  write_phase2_csv(paste0("dissdomassoc_", tag, "_ppearson"), r$p.pearson, "# Pearson p\n")
  write_phase2_csv(
    paste0("dissdomassoc_", tag, "_r2"), matrix(r$Pearson.Rsquare, nrow = 1), "# Pearson R2\n"
  )
  write_phase2_csv(paste0("dissdomassoc_", tag, "_spearman"), r$Spearman, "# Spearman corr\n")
  write_phase2_csv(
    paste0("dissdomassoc_", tag, "_cronbach"), matrix(r$Cronbach, nrow = 1), "# Cronbach alpha\n"
  )
}

gen("uw", NULL)
gen("w", w)

cat("DISSDOMASSOC_GOLDENS_DONE\n")

## R2 weighted-midrank and complete route contract ---------------------------
r2_rows <- list()
r2_add <- function(case, section, domain = "", row = 0L, col = 0L,
                   value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
    case = case,
    section = section,
    domain = domain,
    row = row,
    col = col,
    value = as.numeric(value),
    text = as.character(text),
    stringsAsFactors = FALSE
  )
}

matrix_from_lower <- function(values, n = 5L) {
  out <- matrix(0, nrow = n, ncol = n)
  out[lower.tri(out)] <- values
  out <- out + t(out)
  out
}

unique_dists <- list(
  matrix_from_lower(1:10),
  matrix_from_lower(c(10, 2, 8, 4, 6, 1, 9, 3, 7, 5)),
  matrix_from_lower(c(3, 9, 1, 8, 2, 10, 4, 7, 5, 6))
)
tie_dists <- list(
  matrix_from_lower(c(1, 1, 2, 2, 3, 3, 4, 5, 5, 5)),
  matrix_from_lower(c(2, 1, 2, 1, 3, 4, 3, 4, 4, 4)),
  matrix_from_lower(c(1, 2, 1, 3, 2, 3, 5, 5, 5, 5))
)
equal_dists <- lapply(1:3, function(i) matrix_from_lower(rep(i, 10)))

r2_record_value <- function(case, section, value) {
  if (is.matrix(value)) {
    for (i in seq_len(nrow(value)))
      for (j in seq_len(ncol(value)))
        r2_add(
          case, section, row = i, col = j, value = value[i, j],
          text = paste(rownames(value)[[i]], colnames(value)[[j]], sep = "|")
        )
  } else if (is.list(value) && !is.atomic(value)) {
    stop("unexpected nested result")
  } else {
    value <- unlist(value)
    for (i in seq_along(value))
      r2_add(
        case, section, row = 1L, col = i, value = value[[i]],
        text = gsub(",", ";", names(value)[[i]], fixed = TRUE)
      )
  }
}

r2_run <- function(case, distances, weights = NULL, what = c("pearson", "spearman"),
                   w.rank = FALSE, jointdiss = NULL, dnames = NULL) {
  if (is.null(dnames)) dnames <- paste0("D", seq_along(distances))
  n <- nrow(distances[[1]])
  pair_weights <- as.numeric(as.dist(
    (if (is.null(weights)) rep(1, n) else weights) %*%
      t(if (is.null(weights)) rep(1, n) else weights)
  ))
  all_distances <- distances
  all_names <- dnames
  if (!is.null(jointdiss)) {
    all_distances <- c(all_distances, list(jointdiss))
    all_names <- c(all_names, "joint")
  }
  for (d in seq_along(all_distances)) {
    values <- as.numeric(as.dist(all_distances[[d]]))
    ordinary <- rank(values)
    weighted_rank <- TraMineR:::w.rank(values, w = pair_weights)
    weighted_rank[is.na(weighted_rank)] <- ordinary[is.na(weighted_rank)]
    for (i in seq_along(values)) {
      r2_add(case, "distance", all_names[[d]], i, value = values[[i]])
      r2_add(case, "ordinary_rank", all_names[[d]], i, value = ordinary[[i]])
      r2_add(case, "weighted_rank", all_names[[d]], i, value = weighted_rank[[i]])
    }
  }

  conditions <- character()
  out <- tryCatch(
    withCallingHandlers(
      dissdomassoc(
        distances, jointdiss = jointdiss, what = what, dnames = dnames,
        weights = weights, w.rank = w.rank
      ),
      message = function(m) {
        conditions <<- c(conditions, paste0(
          "MESSAGE:", sub("[[:space:]]+$", "", conditionMessage(m))
        ))
        invokeRestart("muffleMessage")
      },
      warning = function(wrn) {
        conditions <<- c(conditions, paste0(
          "WARNING:", sub("[[:space:]]+$", "", conditionMessage(wrn))
        ))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    r2_add(case, "result_error", text = conditionMessage(out))
  } else {
    for (name in names(out)) {
      if (name == "dnames") {
        for (i in seq_along(out[[name]]))
          r2_add(case, "dnames", col = i, text = out[[name]][[i]])
      } else {
        r2_record_value(case, name, out[[name]])
      }
    }
  }
  if (length(conditions))
    for (i in seq_along(conditions))
      r2_add(case, "condition", row = i, text = conditions[[i]])
}

r2_run(
  "base_unweighted_all", dlist, what = "all",
  dnames = c("Home", "Work", "Care")
)
r2_run(
  "base_fractional_weighted_rank_joint", dlist, weights = w,
  what = c("pearson", "spearman", "R2"), w.rank = TRUE,
  jointdiss = Reduce(`+`, dlist), dnames = c("Home", "Work", "Care")
)
r2_run(
  "base_integer_weighted_rank_spearman",
  dlist, weights = rep(c(1, 2, 3), length.out = inp$n),
  what = "spearman", w.rank = TRUE
)
r2_run(
  "unique_fractional_weighted_rank", unique_dists,
  weights = c(.5, 1, 1.5, 2, 3), what = c("spearman", "R2"), w.rank = TRUE
)
r2_run(
  "ties_fractional_weighted_rank", tie_dists,
  weights = c(.5, 1, 1.5, 2, 3), what = c("pearson", "spearman"), w.rank = TRUE
)
r2_run(
  "ties_fractional_plain_rank", tie_dists,
  weights = c(.5, 1, 1.5, 2, 3), what = c("pearson", "spearman"), w.rank = FALSE
)
r2_run(
  "ties_zero_weights", tie_dists,
  weights = c(0, 1, 0, 2, 3), what = "spearman", w.rank = TRUE
)
r2_run(
  "all_equal_rank", equal_dists,
  weights = c(.5, 1, 1.5, 2, 3), what = "spearman", w.rank = TRUE
)
r2_run("route_pearson", unique_dists, what = "pearson")
r2_run("route_r2", unique_dists, what = "R2")
r2_run("route_cronbach", unique_dists, what = "cronbach")
r2_run("route_cron_subsets", unique_dists, what = "cron.subsets")

utils::write.csv(
  do.call(rbind, r2_rows),
  file.path("tests", "r_goldens", "dissdomassoc_r2_rank.csv"),
  row.names = FALSE,
  na = "NaN",
  quote = FALSE
)

r2_error <- function(expr) {
  tryCatch(
    {
      suppressMessages(force(expr))
      "VALUE"
    },
    error = function(e) paste0(
      "ERROR:", sub("[[:space:]]+$", "", conditionMessage(e))
    )
  )
}
writeLines(
  c(
    paste0("invalid_Ppearson=", r2_error(dissdomassoc(dlist, what = "Ppearson"))),
    paste0("not_list=", r2_error(dissdomassoc(dlist[[1]]))),
    paste0("one_domain=", r2_error(dissdomassoc(dlist[1]))),
    paste0("bad_what=", r2_error(dissdomassoc(dlist, what = "kendall"))),
    paste0("weights_short=", r2_error(dissdomassoc(dlist, weights = w[-1L]))),
    paste0("joint_nonnumeric=", r2_error(dissdomassoc(dlist, jointdiss = "x"))),
    paste0("nonnumeric_domain=", r2_error(dissdomassoc(list(dlist[[1]], "x"))))
  ),
  file.path("tests", "r_goldens", "dissdomassoc_r2_errors.txt")
)

cat("DISSDOMASSOC_R2_GOLDENS_DONE\n")
