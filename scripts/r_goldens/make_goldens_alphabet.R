# Golden generation for TraMineR::alphabet against TraMineR 2.2-14.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_chars <- function(values, name) {
  writeLines(as.character(values), file.path(out_dir, name), useBytes = TRUE)
}

# A has no missing cells. with.missing=TRUE must still append the declared
# missing marker, because alphabet() reads metadata rather than observed cells.
A <- seqdef(
  data.frame(t1 = c("A", "B"), t2 = c("B", "A"), stringsAsFactors = FALSE),
  alphabet = c("A", "B")
)

# B contains one interior missing state and one right-side structural void.
# Structural void is never part of alphabet(), while nr is appended on request.
B <- seqdef(
  data.frame(
    t1 = c("A", "B"),
    t2 = c(NA, "A"),
    t3 = c("B", NA),
    stringsAsFactors = FALSE
  ),
  alphabet = c("A", "B")
)

# Event levels come from levels(seqelist); with.missing is ignored for events.
E <- seqecreate(data.frame(
  id = c(1, 1, 2, 2),
  timestamp = c(0, 2, 0, 1),
  event = c("B", "A", "C", "B"),
  stringsAsFactors = FALSE
))

write_chars(alphabet(A, with.missing = FALSE), "alphabet_A_wm0.txt")
write_chars(alphabet(A, with.missing = TRUE), "alphabet_A_wm1.txt")
write_chars(alphabet(B, with.missing = TRUE), "alphabet_B_wm1.txt")
write_chars(alphabet(E, with.missing = TRUE), "alphabet_E_events.txt")

stopifnot(
  identical(alphabet(A, with.missing = TRUE), c(attr(A, "alphabet"), attr(A, "nr"))),
  identical(alphabet(B, with.missing = TRUE), c(attr(B, "alphabet"), attr(B, "nr"))),
  identical(alphabet(E, with.missing = FALSE), alphabet(E, with.missing = TRUE))
)

cat("ALPHABET_GOLDENS_DONE\n")
cat("A wm0:", paste(alphabet(A), collapse = ","), "\n")
cat("A wm1:", paste(alphabet(A, with.missing = TRUE), collapse = ","), "\n")
cat("B wm1:", paste(alphabet(B, with.missing = TRUE), collapse = ","), "\n")
cat("E:", paste(alphabet(E), collapse = ","), "\n")
