# Golden generation: RP6b seqimplic (implicative statistics of group membership).
# Anchors traminex's seqimplic to TraMineRextras 0.6.9. Run in the live-R env:
#   Rscript scripts/r_goldens/make_goldens_seqimplic.R
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
stopifnot(as.character(packageVersion("TraMineRextras")) == "0.6.9")
out_dir <- file.path("tests", "r_goldens")

write_num_matrix <- function(name, m, header) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  cat(header, file = con)
  txt <- apply(as.matrix(m), 1, function(r) paste(formatC(r, format = "g", digits = 15), collapse = ","))
  writeLines(txt, con)
  close(con)
}

A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))
sdAw <- seqdef(A, alphabet = c("A", "B", "C", "D"), weights = c(1.5, 2, 0.5, 3))
B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)
grp <- c(1, 1, 2, 2)

# $indices: array (levels x alphabet x positions); as.vector is column-major
# (level fastest, then alphabet, then position) -- matches numpy flatten(order="F").
imp_A <- suppressWarnings(seqimplic(sdA, group = grp))
write_num_matrix("seqimplic_A", as.vector(imp_A$indices), "# R: seqimplic(sdA, c(1,1,2,2)); dims 2x4x6, F-order\n")

imp_Aw <- suppressWarnings(seqimplic(sdAw, group = grp))
write_num_matrix("seqimplic_A_w", as.vector(imp_Aw$indices), "# R: seqimplic(sdAw weighted, c(1,1,2,2)); dims 2x4x6\n")

imp_B <- suppressWarnings(seqimplic(sdB, group = grp, with.missing = TRUE))
write_num_matrix("seqimplic_B_wm", as.vector(imp_B$indices), "# R: seqimplic(sdB, c(1,1,2,2), wm=T); dims 2x4x6 (A,B,C,nr)\n")

custom_time <- A
colnames(custom_time) <- paste0("wave-", LETTERS[seq_len(ncol(custom_time))])
sd_custom_time <- seqdef(custom_time, alphabet = c("A", "B", "C", "D"))
imp_custom_time <- suppressWarnings(seqimplic(sd_custom_time, group = grp))

# R4 plot.seqimplic transformation and metadata. These artifacts intentionally
# leave the numerical RP6b arrays above unchanged.
plot_records <- list()
plot_index <- 1L
for (case_name in names(list(A = imp_A, B_missing = imp_B))) {
  object <- list(A = imp_A, B_missing = imp_B)[[case_name]]
  plotted <- -object$indices
  plotted[plotted < 0] <- NA_real_
  for (level_index in seq_along(object$levels)) {
    for (state_index in seq_along(object$alphabet)) {
      for (position_index in seq_len(dim(object$indices)[3])) {
        plot_records[[plot_index]] <- data.frame(
          case = case_name,
          level_index_zero = level_index - 1L,
          level = object$levels[level_index],
          state_index_zero = state_index - 1L,
          state = object$alphabet[state_index],
          position_index_zero = position_index - 1L,
          position_label = dimnames(object$indices)[[3]][position_index],
          stored = object$indices[level_index, state_index, position_index],
          plotted = plotted[level_index, state_index, position_index],
          stringsAsFactors = FALSE
        )
        plot_index <- plot_index + 1L
      }
    }
  }
}
write.csv(
  do.call(rbind, plot_records),
  file.path(out_dir, "seqimplic_plot_r4_series.csv"),
  row.names = FALSE,
  quote = TRUE,
  na = "NA"
)

metadata <- data.frame(
  case = c(
    rep("A", 10),
    rep("B_missing", 5),
    "custom_time_labels",
    rep("contract", 4)
  ),
  key = c(
    "type", "levels", "alphabet", "labels", "cpal", "xtstep", "tick.last",
    "positions", "default_ylim", "qnorm_0.95",
    "levels", "alphabet", "labels", "cpal", "positions",
    "positions",
    "stored_sign", "negative_plotted_values", "intensity_ylim", "intensity_threshold"
  ),
  value = c(
    imp_A$type,
    paste(imp_A$levels, collapse = "\036"),
    paste(imp_A$alphabet, collapse = "\036"),
    paste(imp_A$labels, collapse = "\036"),
    paste(imp_A$cpal, collapse = "\036"),
    imp_A$xtstep,
    imp_A$tick.last,
    paste(dimnames(imp_A$indices)[[3]], collapse = "\036"),
    paste(c(0, max(-imp_A$indices, na.rm = TRUE)), collapse = "\036"),
    qnorm(0.95),
    paste(imp_B$levels, collapse = "\036"),
    paste(imp_B$alphabet, collapse = "\036"),
    paste(imp_B$labels, collapse = "\036"),
    paste(imp_B$cpal, collapse = "\036"),
    paste(dimnames(imp_B$indices)[[3]], collapse = "\036"),
    paste(dimnames(imp_custom_time$indices)[[3]], collapse = "\036"),
    "plot negates x$indices on a local copy",
    "NA",
    "0\0361",
    "conf"
  ),
  stringsAsFactors = FALSE
)
for (entry in list(c("layout_right_0.2", 0.2), c("layout_right_0.4", 0.4))) {
  proportion <- as.numeric(entry[[2]])
  r_layout <- TraMineR:::TraMineR.setlayout(
    length(imp_A$levels), NA, NA, "right", "all", proportion
  )
  metadata <- rbind(
    metadata,
    data.frame(
      case = entry[[1]],
      key = c("plot_fraction", "legend_fraction"),
      value = c(
        sum(r_layout$widths[-length(r_layout$widths)]),
        r_layout$widths[length(r_layout$widths)]
      ),
      stringsAsFactors = FALSE
    )
  )
}
write.csv(
  metadata,
  file.path(out_dir, "seqimplic_plot_r4_metadata.csv"),
  row.names = FALSE,
  quote = TRUE
)

conditions <- data.frame(
  case = c("default_levels", "duplicate_levels", "unknown_levels", "confidence_bounds"),
  value = c(
    paste(imp_A$levels, collapse = "\036"),
    "membership selection plots each declared level at most once",
    "unknown levels select no panel",
    "confidence values are probabilities passed to qnorm for indice objects"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  conditions,
  file.path(out_dir, "seqimplic_plot_r4_conditions.csv"),
  row.names = FALSE,
  quote = TRUE
)

cat("SEQIMPLIC_GOLDENS_DONE\n")
