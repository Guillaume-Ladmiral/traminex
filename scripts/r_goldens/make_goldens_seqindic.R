# Golden generation: seqindic (T11) -- composite indicator table.
#
# Pins traminex.seqindic (re-targeted from a one-hot matrix to the R indicator
# table) to TraMineR 2.2-14 seqindic(indic, with.missing). Pins the default set
# {visited,trans,entr,cplx,turb2n} across A/B/C, plus a broad basic+diversity+
# complexity set on A. Column order/labels follow R.
#
# Run from the repository root with the oracle library (scripts/r_goldens/README.md):
#   Rscript scripts/r_goldens/make_goldens_seqindic.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_corestats.R"))

if (!identical(as.character(packageVersion("TraMineR")), "2.2.14")) {
  stop("seqindic goldens require TraMineR 2.2-14")
}

DEF <- c("visited", "trans", "entr", "cplx", "turb2n")
BROAD <- c(
  "lgth", "nonm", "dlgth", "visited", "visitp", "recu",
  "trans", "transp", "entr", "cplx", "turb", "turbn", "turb2", "turb2n"
)

emit <- function(sd, tag, wm, indic) {
  t <- suppressMessages(seqindic(sd, indic = indic, with.missing = wm))
  write_num_matrix(
    paste0("seqindic_", tag), as.matrix(t),
    sprintf("# R: seqindic(%s, indic=c(%s), with.missing=%s)\n", tag, paste(indic, collapse = ","), wm)
  )
}

csv_quote <- function(value) {
  paste0('"', gsub('"', '""', as.character(value), fixed = TRUE), '"')
}

format_numeric <- function(value) {
  if (is.nan(value)) return("NaN")
  if (is.na(value)) return("NA")
  if (is.infinite(value)) return(if (value > 0) "Inf" else "-Inf")
  trimws(formatC(value, format = "g", digits = 15))
}

write_seqindic_cases <- function(filename, cases) {
  lines <- "case,row,column,value"
  for (case in names(cases)) {
    values <- as.data.frame(cases[[case]])
    if (ncol(values) == 0L) {
      for (row in seq_len(nrow(values))) {
        lines <- c(lines, paste(csv_quote(case), row, csv_quote(""), "", sep = ","))
      }
      next
    }
    for (row in seq_len(nrow(values))) {
      for (column in seq_len(ncol(values))) {
        lines <- c(
          lines,
          paste(
            csv_quote(case),
            row,
            csv_quote(names(values)[[column]]),
            format_numeric(values[[column]][[row]]),
            sep = ","
          )
        )
      }
    }
  }
  writeLines(
    lines,
    con = file.path(out_dir, filename),
    sep = "\n",
    useBytes = TRUE
  )
}

duration_cases <- function(sd, with_missing) {
  list(
    together_scrambled_w05 = suppressMessages(
      seqindic(
        sd,
        indic = c("volat", "dustd2", "meand2", "nsubs", "dustd", "meand"),
        with.missing = with_missing,
        w = .5
      )
    ),
    meand = suppressMessages(
      seqindic(sd, indic = "meand", with.missing = with_missing)
    ),
    meand2 = suppressMessages(
      seqindic(sd, indic = "meand2", with.missing = with_missing)
    ),
    dustd = suppressMessages(
      seqindic(sd, indic = "dustd", with.missing = with_missing)
    ),
    dustd2 = suppressMessages(
      seqindic(sd, indic = "dustd2", with.missing = with_missing)
    ),
    nsubs = suppressMessages(
      seqindic(sd, indic = "nsubs", with.missing = with_missing)
    ),
    volat_w0 = suppressMessages(
      seqindic(sd, indic = "volat", with.missing = with_missing, w = 0)
    ),
    volat_w05 = suppressMessages(
      seqindic(sd, indic = "volat", with.missing = with_missing, w = .5)
    ),
    volat_w1 = suppressMessages(
      seqindic(sd, indic = "volat", with.missing = with_missing, w = 1)
    )
  )
}

capture_error_line <- function(case, expression) {
  condition <- tryCatch(
    {
      suppressMessages(eval.parent(substitute(expression)))
      NULL
    },
    error = function(err) err
  )
  if (is.null(condition)) stop("expected R error for case: ", case)
  message <- trimws(gsub("[\r\n]+", " ", conditionMessage(condition)))
  paste("error", case, message, sep = "|")
}

capture_value_warnings <- function(case, expression) {
  warning_lines <- character()
  value <- withCallingHandlers(
    eval.parent(substitute(expression)),
    message = function(condition) {
      message <- trimws(gsub("[\r\n]+", " ", conditionMessage(condition)))
      if (grepl("\\[!!\\]", message)) {
        warning_lines <<- c(warning_lines, paste("warning", case, message, sep = "|"))
      }
      invokeRestart("muffleMessage")
    },
    warning = function(condition) {
      message <- trimws(gsub("[\r\n]+", " ", conditionMessage(condition)))
      warning_lines <<- c(warning_lines, paste("warning", case, message, sep = "|"))
      invokeRestart("muffleWarning")
    }
  )
  list(value = value, lines = warning_lines)
}

write_condition_lines <- function(filename, lines) {
  writeLines(
    c("kind|case|message", lines),
    con = file.path(out_dir, filename),
    sep = "\n",
    useBytes = TRUE
  )
}

fx <- corestats_fixtures()

binary_all <- c("ppos", "nvolat", "vpos", "integr")
binary_args_a_false <- list(
  pos.states = c("A", "B"),
  neg.states = c("C", "D"),
  dss = FALSE,
  pow = 2,
  w = .25
)
binary_args_a_true <- binary_args_a_false
binary_args_a_true$dss <- TRUE
binary_args_b <- list(
  pos.states = c("A", "*"),
  neg.states = c("B", "C"),
  dss = FALSE,
  pow = 2,
  w = .75
)

binary_cases <- list(
  all_dss_false = suppressMessages(
    seqindic(fx$A, indic = binary_all, ipos.args = binary_args_a_false)
  ),
  all_dss_true = suppressMessages(
    seqindic(fx$A, indic = binary_all, ipos.args = binary_args_a_true)
  ),
  ppos_forces_nondss = suppressMessages(
    seqindic(
      fx$A,
      indic = "ppos",
      ipos.args = list(pos.states = c("A", "B"), neg.states = c("C", "D"), dss = TRUE)
    )
  ),
  ppos_unclassified_nan = suppressMessages(
    seqindic(
      fx$A,
      indic = "ppos",
      ipos.args = list(pos.states = "A", neg.states = "B")
    )
  ),
  ppos_phantom_missing_state = suppressMessages(
    seqindic(
      fx$A,
      indic = "ppos",
      with.missing = TRUE,
      ipos.args = list(pos.states = "*", neg.states = "A")
    )
  ),
  nvolat_forces_dss = suppressMessages(
    seqindic(
      fx$A,
      indic = "nvolat",
      ipos.args = list(pos.states = c("A", "B"), neg.states = c("C", "D"), dss = FALSE)
    )
  ),
  vpos_dss_false_w0 = suppressMessages(
    seqindic(
      fx$A,
      indic = "vpos",
      ipos.args = list(
        pos.states = c("A", "B"), neg.states = c("C", "D"), dss = FALSE, w = 0
      )
    )
  ),
  vpos_dss_true_w1 = suppressMessages(
    seqindic(
      fx$A,
      indic = "vpos",
      ipos.args = list(
        pos.states = c("A", "B"), neg.states = c("C", "D"), dss = TRUE, w = 1
      )
    )
  ),
  integr_dss_false_pow2 = suppressMessages(
    seqindic(
      fx$A,
      indic = "integr",
      ipos.args = list(
        pos.states = c("A", "B"), neg.states = c("C", "D"), dss = FALSE, pow = 2
      )
    )
  ),
  missing_all = suppressMessages(
    seqindic(
      fx$B,
      indic = binary_all,
      with.missing = TRUE,
      ipos.args = binary_args_b
    )
  )
)

binary_overwrite <- capture_value_warnings(
  "overwritten_binary_args",
  seqindic(
    fx$A,
    indic = binary_all,
    with.missing = FALSE,
    ipos.args = list(
      seqdata = fx$B,
      with.missing = TRUE,
      index = "integr",
      pos.states = c("A", "B"),
      neg.states = c("C", "D"),
      dss = FALSE
    )
  )
)
binary_cases[["overwritten_binary_args"]] <- binary_overwrite$value

binary_conditions <- c(
  binary_overwrite$lines,
  capture_error_line("empty_ipos_args", seqindic(fx$A, indic = "ppos")),
  capture_error_line("invalid_indicator", seqindic(fx$A, indic = "not-an-indicator")),
  capture_error_line(
    "invalid_pos_state",
    seqindic(fx$A, indic = "ppos", ipos.args = list(pos.states = "Z"))
  ),
  capture_error_line(
    "duplicate_pos_state",
    seqindic(fx$A, indic = "ppos", ipos.args = list(pos.states = c("A", "A")))
  )
)

ranked_all <- c("degrad", "bad", "prec", "insec")
ranked_reverse_args <- list(
  state.order = c("D", "C", "B", "A"),
  method = "RANK",
  pow = 2,
  spow = 1,
  otto = .35,
  a = .8,
  b = 1.1,
  unused.argument = "filtered"
)
ranked_missing_args <- list(
  state.order = c("*", "C", "B", "A"),
  method = "RANK",
  pow = 2
)

ranked_default <- capture_value_warnings(
  "empty_prec_args",
  seqindic(fx$A, indic = ranked_all)
)
ranked_overwrite <- capture_value_warnings(
  "overwritten_ranked_args",
  seqindic(
    fx$A,
    indic = ranked_all,
    with.missing = FALSE,
    prec.args = list(seqdata = fx$B, with.missing = TRUE, method = "RANK")
  )
)
ranked_ipos_type <- capture_value_warnings(
  "ipos_type_source_bug",
  seqindic(
    fx$A,
    indic = "prec",
    ipos.args = list(type = 2),
    prec.args = list(method = "TRATEDSS")
  )
)
ranked_prec_type <- capture_value_warnings(
  "prec_type_silently_overwritten",
  seqindic(fx$A, indic = "prec", prec.args = list(type = 2))
)
if (length(ranked_prec_type$lines) != 0L) {
  stop("TraMineR source behavior changed: prec.args$type unexpectedly warned")
}

ranked_cases <- list(
  all_default = ranked_default$value,
  all_explicit_reverse = suppressMessages(
    seqindic(fx$A, indic = ranked_all, prec.args = ranked_reverse_args)
  ),
  all_partial_order = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      prec.args = list(state.order = c("B", "C"))
    )
  ),
  all_state_equiv = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      prec.args = list(state.equiv = list(c("B", "C")))
    )
  ),
  all_partial_equiv = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      prec.args = list(
        state.order = c("A", "D"),
        state.equiv = list(c("D", "B"))
      )
    )
  ),
  all_explicit_stprec = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      prec.args = list(stprec = c(0, .5, .5, -1))
    )
  ),
  all_stprec_state_equiv = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      prec.args = list(
        stprec = c(0, .2, .8, 1),
        state.equiv = list(c("B", "C"))
      )
    )
  ),
  phantom_missing_all = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      with.missing = TRUE,
      prec.args = list(unused.argument = 42)
    )
  ),
  phantom_missing_equiv = suppressMessages(
    seqindic(
      fx$A,
      indic = ranked_all,
      with.missing = TRUE,
      prec.args = list(state.equiv = list(c("*", "A")))
    )
  ),
  degrad_default = suppressMessages(seqindic(fx$A, indic = "degrad")),
  degrad_single_state_order = suppressMessages(
    seqindic(
      fx$A,
      indic = "degrad",
      prec.args = list(state.order = "A")
    )
  ),
  degrad_constant_stprec = suppressMessages(
    seqindic(
      fx$A,
      indic = "degrad",
      prec.args = list(stprec = rep(1, 4))
    )
  ),
  bad_default = suppressMessages(seqindic(fx$A, indic = "bad")),
  bad_single_state_order = suppressMessages(
    seqindic(
      fx$A,
      indic = "bad",
      prec.args = list(state.order = "A")
    )
  ),
  bad_constant_stprec = suppressMessages(
    seqindic(
      fx$A,
      indic = "bad",
      prec.args = list(stprec = rep(1, 4))
    )
  ),
  prec_default = suppressMessages(seqindic(fx$A, indic = "prec")),
  insec_default_spow0 = suppressMessages(seqindic(fx$A, indic = "insec")),
  insec_explicit_spow1 = suppressMessages(
    seqindic(fx$A, indic = "insec", prec.args = list(pow = 2, spow = 1))
  ),
  prec_explicit_correction = suppressMessages(
    seqindic(
      fx$A,
      indic = "prec",
      prec.args = list(correction = c(-.1, .2, 0, .5))
    )
  ),
  insec_explicit_correction = suppressMessages(
    seqindic(
      fx$A,
      indic = "insec",
      prec.args = list(correction = c(-.1, .2, 0, .5))
    )
  ),
  private_degrad_options = suppressMessages(
    seqindic(
      fx$A,
      indic = c("prec", "insec"),
      prec.args = list(
        penalized = "POS",
        method = "TRATEDSS",
        weight.type = "INV",
        border.effect = 20,
        pow = 2
      )
    )
  ),
  missing_all = suppressMessages(
    seqindic(
      fx$B,
      indic = ranked_all,
      with.missing = TRUE,
      prec.args = ranked_missing_args
    )
  ),
  filtered_unknown = suppressMessages(
    seqindic(fx$A, indic = ranked_all, prec.args = list(unused.argument = 42))
  ),
  overwritten_ranked_args = ranked_overwrite$value,
  prec_type_silently_overwritten = ranked_prec_type$value,
  ipos_type_source_bug = ranked_ipos_type$value
)

ranked_conditions <- c(
  ranked_default$lines,
  ranked_overwrite$lines,
  ranked_ipos_type$lines,
  capture_error_line(
    "prec_spell_integr_duplicate",
    seqindic(fx$A, indic = "prec", prec.args = list(spell.integr = FALSE))
  ),
  capture_error_line(
    "insec_spell_integr_duplicate",
    seqindic(fx$A, indic = "insec", prec.args = list(spell.integr = TRUE))
  ),
  capture_error_line(
    "prec_single_state_order",
    seqindic(fx$A, indic = "prec", prec.args = list(state.order = "A"))
  ),
  capture_error_line(
    "insec_single_state_order",
    seqindic(fx$A, indic = "insec", prec.args = list(state.order = "A"))
  ),
  capture_error_line(
    "prec_constant_stprec",
    seqindic(fx$A, indic = "prec", prec.args = list(stprec = rep(1, 4)))
  ),
  capture_error_line(
    "insec_constant_stprec",
    seqindic(fx$A, indic = "insec", prec.args = list(stprec = rep(1, 4)))
  )
)

group_prec_args <- list(method = "RANK")
group_cases <- list(
  basic = suppressMessages(seqindic(fx$A, indic = "basic")),
  diversity = suppressMessages(seqindic(fx$A, indic = "diversity")),
  diversity_plus_visited = suppressMessages(
    seqindic(fx$A, indic = c("diversity", "visited"))
  ),
  visited_visitp_phantom_missing = suppressMessages(
    seqindic(
      fx$A,
      indic = c("visited", "visitp"),
      with.missing = TRUE
    )
  ),
  complexity = suppressMessages(seqindic(fx$A, indic = "complexity")),
  binary = suppressMessages(
    seqindic(fx$A, indic = "binary", ipos.args = binary_args_a_false)
  ),
  ranked = suppressMessages(
    seqindic(fx$A, indic = "ranked", prec.args = group_prec_args)
  ),
  all_only = suppressMessages(seqindic(fx$A, indic = "all")),
  aliases = suppressMessages(
    seqindic(
      fx$A,
      indic = c("prec2", "inpos", "visit"),
      ipos.args = binary_args_a_false,
      prec.args = group_prec_args
    )
  ),
  mixed_scrambled_duplicates = suppressMessages(
    seqindic(
      fx$A,
      indic = c("volat", "meand", "trans", "meand", "cplx", "dustd2")
    )
  ),
  direct_visitp = suppressMessages(seqindic(fx$A, indic = "visitp"))
)

emit(fx$A, "A", FALSE, DEF)
emit(fx$B, "B_wm0", FALSE, DEF)
emit(fx$B, "B_wm1", TRUE, DEF)
emit(fx$C, "C_wm0", FALSE, DEF)
emit(fx$C, "C_wm1", TRUE, DEF)
emit(fx$A, "A_broad", FALSE, BROAD)
write_seqindic_cases("seqindic_A_duration.csv", duration_cases(fx$A, FALSE))
write_seqindic_cases("seqindic_C_duration_wm1.csv", duration_cases(fx$C, TRUE))
write_seqindic_cases("seqindic_binary.csv", binary_cases)
write_condition_lines("seqindic_binary_errors.txt", binary_conditions)
write_seqindic_cases("seqindic_ranked.csv", ranked_cases)
write_seqindic_cases("seqindic_groups.csv", group_cases)
write_condition_lines("seqindic_ranked_warnings.txt", ranked_conditions)
# Row names of the result (TX-AUD-20261010-093): R names the rows by rownames(seqdata), the sequence
# ids, "1".."n" when seqdef gets none; a subset keeps its rows' names.
ids_data <- rbind(c("A", "A", "B"), c("A", "B", "B"), c("B", "B", "A"))
rowname_cases <- list(
  character_ids = suppressMessages(seqdef(ids_data, alphabet = c("A", "B"), id = c("p1", "p2", "p3"))),
  default_ids = suppressMessages(seqdef(ids_data, alphabet = c("A", "B"))),
  numeric_ids = suppressMessages(seqdef(ids_data, alphabet = c("A", "B"), id = c(10, 20, 30)))
)
rowname_cases$reordered_subset <- rowname_cases$character_ids[c(3, 1), ]
rowname_lines <- "case,row,rowname"
for (case in names(rowname_cases)) {
  names_out <- rownames(suppressMessages(seqindic(rowname_cases[[case]], indic = "visited")))
  rowname_lines <- c(rowname_lines, paste(csv_quote(case), seq_along(names_out), csv_quote(names_out), sep = ","))
}
writeLines(rowname_lines, con = file.path(out_dir, "seqindic_rownames.csv"), sep = "\n", useBytes = TRUE)
cat("SEQINDIC_GOLDENS_DONE\n")
