# Golden generation for TraMineRextras 0.6.9 transversal-entropy plots.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw <- data.frame(
  t1 = c("A", "A", "B", "C", "A", "B"),
  t2 = c("A", "B", "B", "C", NA, "A"),
  t3 = c("B", "B", "A", "C", "B", "A"),
  t4 = c("B", "C", "A", "C", NA, "A"),
  t5 = c("C", "C", "A", "C", NA, NA),
  stringsAsFactors = FALSE
)
seqdata <- suppressMessages(
  seqdef(
    raw,
    alphabet = c("A", "B", "C"),
    right = "DEL",
    missing = NA,
    weights = c(1, 2, 1, 3, 2, 1),
    start = 15,
    xtstep = 2,
    tick.last = TRUE
  )
)
group <- factor(
  c("beta", "alpha", "beta", "gamma", NA, "alpha"),
  levels = c("unused", "gamma", "beta", "alpha")
)

entropy_records <- list()
record_index <- 1L
record_series <- function(case_name, curve_name, result) {
  entropy_records[[record_index]] <<- data.frame(
    case = case_name,
    curve = curve_name,
    position_one = seq_along(result$Entropy),
    entropy = result$Entropy,
    valid = result$ValidStates,
    stringsAsFactors = FALSE
  )
  record_index <<- record_index + 1L
}

factored <- factor(group)
grouped <- by(seqdata, factored, seqstatd)
for (curve_name in names(grouped)) {
  record_series("grouped", curve_name, grouped[[curve_name]])
}

short <- seqdata[, 1:3]
multiple <- list(full = seqdata, short = short, first_three = seqdata[1:3, ])
for (curve_name in names(multiple)) {
  record_series("multiple", curve_name, seqstatd(multiple[[curve_name]]))
}

write.csv(
  do.call(rbind, entropy_records),
  file.path(out_dir, "tentrop_r4_series.csv"),
  row.names = FALSE,
  quote = TRUE
)

all_entropy <- unlist(lapply(grouped, function(value) value$Entropy))
ylim <- c(floor(10 * min(all_entropy)), ceiling(10 * max(all_entropy))) / 10
group_colors <- colorspace::qualitative_hcl(length(grouped), palette = "Dark 3")
list_colors <- RColorBrewer::brewer.pal(9, "Set1")[seq_along(multiple)]
metadata <- data.frame(
  case = c(
    rep("grouped", 8),
    rep("multiple", 6),
    rep("contract", 2)
  ),
  key = c(
    "curve_names", "curve_count", "colors", "line_types", "line_widths",
    "ylim", "xtlab", "tick_positions",
    "curve_names", "curve_count", "colors", "npos", "xtlab", "unequal_widths",
    "auto_title", "group_factor_levels"
  ),
  value = c(
    paste(names(grouped), collapse = "\036"),
    length(grouped),
    paste(group_colors, collapse = "\036"),
    paste(c("solid", "dashed", "dotted")[seq_along(grouped)], collapse = "\036"),
    paste(rep(3.5, length(grouped)), collapse = "\036"),
    paste(ylim, collapse = "\036"),
    paste(names(seqdata), collapse = "\036"),
    paste(c(1, 3, 5), collapse = "\036"),
    paste(names(multiple), collapse = "\036"),
    length(multiple),
    paste(list_colors, collapse = "\036"),
    max(vapply(multiple, ncol, integer(1))),
    paste(paste0("t", seq_len(max(vapply(multiple, ncol, integer(1))))), collapse = "\036"),
    "TRUE",
    "Tranversal Entropies",
    paste(levels(factored), collapse = "\036")
  ),
  stringsAsFactors = FALSE
)
write.csv(
  metadata,
  file.path(out_dir, "tentrop_r4_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

capture_condition <- function(expr) {
  warning_message <- NULL
  error_message <- NULL
  withCallingHandlers(
    tryCatch(force(expr), error = function(condition) {
      error_message <<- conditionMessage(condition)
      NULL
    }),
    warning = function(condition) {
      warning_message <<- conditionMessage(condition)
      invokeRestart("muffleWarning")
    }
  )
  if (!is.null(error_message)) return(paste("error", error_message, sep = "\036"))
  if (!is.null(warning_message)) return(paste("warning", warning_message, sep = "\036"))
  "none\036"
}

pdf(NULL)
conditions <- data.frame(
  case = c("group_length_mismatch", "multiple_fewer_than_two", "multiple_unequal_widths"),
  condition = c(
    capture_condition(seqplot.tentrop(seqdata, group = c("a", "b"))),
    capture_condition(seqplot.tentrop.m(list(seqdata))),
    capture_condition(seqplot.tentrop.m(list(seqdata, short), with.legend = FALSE))
  ),
  stringsAsFactors = FALSE
)
dev.off()
write.csv(
  conditions,
  file.path(out_dir, "tentrop_r4_conditions.csv"),
  row.names = FALSE,
  quote = TRUE
)
