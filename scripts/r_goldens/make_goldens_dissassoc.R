# Golden generation: dissassoc (Phase 2.2 / DUPLICATES A6).
#
# Anchors traminex's R-faithful weighted `dissassoc` (pseudo ANOVA of a
# dissimilarity matrix) to TraMineR 2.2-14 `dissassoc(diss, group, weights=,
# R=, weight.permutation=, squared=)`.
#
# R=1 is used so the observed statistics (t0) are deterministic: the five
# test values [Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett, Levene], the
# pseudo-ANOVA table and the per-group discrepancy table are pinned; p-values
# are permutation-based and stochastic, so they are checked structurally on
# the Python side only (see MIGRATION_PLAN.md Phase 2.2).
#
# Inputs come from make_goldens_phase2_inputs.R (shared Phase 2 fixture).
# Requires R + TraMineR 2.2-14; see scripts/r_goldens/README.md.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D
grp <- inp$grp

dump_dissassoc <- function(tag, res, rcall) {
  write_phase2_csv(
    paste0("dissassoc_", tag, "_anova"), as.matrix(res$anova.table),
    paste0("# R: ", rcall, " -> anova.table [SS, df, MSE] rows Exp/Res/Total\n")
  )
  write_phase2_csv(
    paste0("dissassoc_", tag, "_groups"), as.matrix(res$groups),
    paste0("# R: ", rcall, " -> groups [n, discrepancy] rows g1,g2,g3,Total\n")
  )
  write_phase2_csv(
    paste0("dissassoc_", tag, "_t0"), matrix(res$stat[, 1], ncol = 1),
    paste0(
      "# R: ", rcall,
      " -> stat t0 rows [Pseudo F, Pseudo Fbf, Pseudo R2, Bartlett, Levene]\n"
    )
  )
}

res_unw <- dissassoc(D, grp, R = 1)
dump_dissassoc("unw", res_unw, "dissassoc(D, grp, R=1)")

res_wreal <- dissassoc(D, grp, weights = inp$wreal, R = 1, weight.permutation = "diss")
dump_dissassoc(
  "wreal", res_wreal,
  "dissassoc(D, grp, weights=wreal, R=1, weight.permutation='diss')"
)

res_wint <- dissassoc(D, grp, weights = inp$wint, R = 1, weight.permutation = "replicate")
dump_dissassoc(
  "wint_rep", res_wint,
  "dissassoc(D, grp, weights=wint, R=1, weight.permutation='replicate')"
)

res_wreal_sq <- dissassoc(D, grp,
  weights = inp$wreal, R = 1,
  weight.permutation = "diss", squared = TRUE
)
write_phase2_csv(
  "dissassoc_wreal_sq_t0", matrix(res_wreal_sq$stat[, 1], ncol = 1),
  "# R: dissassoc(D, grp, weights=wreal, R=1, weight.permutation='diss', squared=TRUE) -> stat t0\n"
)

cat("DISSASSOC_GOLDENS_DONE\n")

## R2: complete weight.permutation contract ----------------------------------
r2_rows <- list()
r2_add <- function(case, method, weights_kind, R, squared, group_kind,
                   contract, section, row = 0L, col = 0L,
                   value = NaN, text = "") {
  r2_rows[[length(r2_rows) + 1L]] <<- data.frame(
    case = case,
    method = method,
    weights_kind = weights_kind,
    R = R,
    squared = squared,
    group_kind = group_kind,
    contract = contract,
    section = section,
    row = row,
    col = col,
    value = as.numeric(value),
    text = as.character(text),
    stringsAsFactors = FALSE
  )
}

r2_run <- function(case, method, weights_kind, R, squared = FALSE,
                   group_kind = "three", seed = 8101L) {
  ww <- switch(
    weights_kind,
    none = NULL,
    integer = inp$wint,
    fractional = inp$wreal
  )
  gg <- if (group_kind == "two") rep(c(1L, 2L), length.out = nrow(D)) else grp
  set.seed(seed)
  conditions <- character()
  out <- withCallingHandlers(
    dissassoc(
      D, gg, weights = ww, R = R,
      weight.permutation = method, squared = squared
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
  )
  effective_method <- out$weight.permutation
  r2_add(
    case, method, weights_kind, R, squared, group_kind,
    "EXACT", "attribute", text = paste0("effective_method=", effective_method)
  )
  for (i in seq_len(nrow(out$groups)))
    for (j in seq_len(ncol(out$groups)))
      r2_add(
        case, method, weights_kind, R, squared, group_kind,
        "EXACT", "groups", i, j, out$groups[i, j],
        paste(rownames(out$groups)[[i]], colnames(out$groups)[[j]], sep = "|")
      )
  for (i in seq_len(nrow(out$anova.table)))
    for (j in seq_len(ncol(out$anova.table)))
      r2_add(
        case, method, weights_kind, R, squared, group_kind,
        "EXACT", "anova", i, j, out$anova.table[i, j],
        paste(rownames(out$anova.table)[[i]], colnames(out$anova.table)[[j]], sep = "|")
      )
  for (i in seq_len(nrow(out$stat))) {
    r2_add(
      case, method, weights_kind, R, squared, group_kind,
      "EXACT", "t0", i, 1L, out$stat[i, 1], rownames(out$stat)[[i]]
    )
    r2_add(
      case, method, weights_kind, R, squared, group_kind,
      if (R > 1L) "R_RNG_ONLY" else "EXACT",
      "pvalue", i, 1L, out$stat[i, 2], rownames(out$stat)[[i]]
    )
  }
  if (R > 1L) {
    for (j in seq_len(ncol(out$perms$t))) {
      values <- out$perms$t[, j]
      summary_values <- c(
        min = min(values), max = max(values), mean = mean(values),
        sum = sum(values), finite = sum(is.finite(values))
      )
      for (i in seq_along(summary_values))
        r2_add(
          case, method, weights_kind, R, squared, group_kind,
          "R_RNG_ONLY", "permutation_summary", j, i,
          summary_values[[i]], names(summary_values)[[i]]
        )
    }
  }
  if (length(conditions))
    for (i in seq_along(conditions))
      r2_add(
        case, method, weights_kind, R, squared, group_kind,
        "EXACT", "condition", i, text = conditions[[i]]
      )
}

r2_run("uw_none_r31", "none", "none", 31L, seed = 8201L)
r2_run("uw_requested_permute_r0", "permute", "none", 0L, seed = 8202L)
r2_run(
  "uw_requested_random_sampling_r1", "random-sampling", "none", 1L,
  group_kind = "two", seed = 8203L
)
r2_run("wint_diss_r0", "diss", "integer", 0L, seed = 8204L)
r2_run("wreal_diss_r31_sq", "diss", "fractional", 31L, squared = TRUE, seed = 8205L)
r2_run("wint_replicate_r1", "replicate", "integer", 1L, seed = 8206L)
r2_run("wint_replicate_r31", "replicate", "integer", 31L, seed = 8207L)
r2_run(
  "wreal_rounded_replicate_r1", "rounded-replicate", "fractional", 1L,
  seed = 8208L
)
r2_run(
  "wreal_rounded_replicate_r31_two", "rounded-replicate", "fractional", 31L,
  group_kind = "two", seed = 8209L
)
r2_run("wreal_group_r1", "group", "fractional", 1L, seed = 8210L)
r2_run("wreal_group_r31", "group", "fractional", 31L, seed = 8211L)
r2_run("wint_permute_r1", "permute", "integer", 1L, seed = 8212L)
r2_run("wint_permute_r31", "permute", "integer", 31L, seed = 8213L)
r2_run(
  "wreal_permute_r31_sq_two", "permute", "fractional", 31L,
  squared = TRUE, group_kind = "two", seed = 8214L
)
r2_run("wint_random_sampling_r1", "random-sampling", "integer", 1L, seed = 8215L)
r2_run("wint_random_sampling_r31", "random-sampling", "integer", 31L, seed = 8216L)
r2_run(
  "wreal_random_sampling_r31_sq_two", "random-sampling", "fractional", 31L,
  squared = TRUE, group_kind = "two", seed = 8217L
)

utils::write.csv(
  do.call(rbind, r2_rows),
  file.path("tests", "r_goldens", "dissassoc_r2_cases.csv"),
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
    paste0("invalid_method=", r2_error(dissassoc(
      D, grp, weights = inp$wreal, R = 1L, weight.permutation = "bogus"
    ))),
    paste0("fractional_replicate=", r2_error(dissassoc(
      D, grp, weights = inp$wreal, R = 1L, weight.permutation = "replicate"
    ))),
    paste0("weights_short=", r2_error(dissassoc(
      D, grp, weights = inp$wreal[-1L], R = 1L, weight.permutation = "diss"
    ))),
    paste0("group_short=", r2_error(dissassoc(
      D, grp[-1L], weights = inp$wreal, R = 1L, weight.permutation = "diss"
    ))),
    paste0("negative_R=", r2_error(dissassoc(
      D, grp, weights = inp$wreal, R = -1L, weight.permutation = "diss"
    ))),
    paste0("negative_sampling_weight=", r2_error(dissassoc(
      D, grp, weights = replace(inp$wreal, 1L, -1), R = 3L,
      weight.permutation = "random-sampling"
    ))),
    paste0("zero_sampling_weights=", r2_error(dissassoc(
      D, grp, weights = rep(0, nrow(D)), R = 3L,
      weight.permutation = "random-sampling"
    )))
  ),
  file.path("tests", "r_goldens", "dissassoc_r2_errors.txt")
)

cat("DISSASSOC_R2_GOLDENS_DONE\n")
