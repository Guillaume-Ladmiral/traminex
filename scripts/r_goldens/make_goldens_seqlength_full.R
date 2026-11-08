# Golden generation for the full TraMineR::seqlength contract.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_lengths <- function(values, name, call) {
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  on.exit(close(con), add = TRUE)
  cat("# R: ", call, "\n", sep = "", file = con)
  cat("# rownames: ", paste(rownames(values), collapse = ","), "\n", sep = "", file = con)
  utils::write.table(
    values,
    con,
    sep = ",",
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
}

# B contains interior missing states but no structural void. right=NA keeps
# trailing missing cells as missing, rather than applying seqdef's DEL policy.
B <- rbind(
  c("A", "A", "*", "B", "B", "C"),
  c("*", "A", "A", "A", "*", "*"),
  c("A", "*", "A", "B", "B", "B"),
  c("C", "C", "C", "C", "C", "C")
)
rownames(B) <- paste0("B", seq_len(nrow(B)))
sdB <- seqdef(B, alphabet = c("A", "B", "C"), missing = "*", right = NA)

# C is ragged. TraMineR's default right="DEL" converts trailing missing cells
# to structural void while retaining C3's interior missing state.
C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
rownames(C) <- paste0("C", seq_len(nrow(C)))
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")

B_wm0 <- seqlength(sdB, with.missing = FALSE)
B_wm1 <- seqlength(sdB, with.missing = TRUE)
C_wm0 <- seqlength(sdC, with.missing = FALSE)
C_wm1 <- seqlength(sdC, with.missing = TRUE)

stopifnot(
  identical(as.integer(B_wm0[, 1]), c(5L, 3L, 5L, 6L)),
  identical(as.integer(B_wm1[, 1]), c(6L, 6L, 6L, 6L)),
  identical(as.integer(C_wm0[, 1]), c(6L, 3L, 3L, 2L)),
  identical(as.integer(C_wm1[, 1]), c(6L, 3L, 4L, 2L))
)

write_lengths(B_wm0, "seqlength_B_wm0", "seqlength(sdB, with.missing=FALSE)")
write_lengths(B_wm1, "seqlength_B_wm1", "seqlength(sdB, with.missing=TRUE)")
write_lengths(C_wm0, "seqlength_C_wm0", "seqlength(sdC, with.missing=FALSE)")
write_lengths(C_wm1, "seqlength_C_wm1", "seqlength(sdC, with.missing=TRUE)")

cat("SEQLENGTH_FULL_GOLDENS_DONE\n")
