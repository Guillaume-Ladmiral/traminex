# Golden generation: what R does in the cases of part 4 of the consumed-component gaps
# (supervision repository, docs/TRAMINEX_COMPONENT_GAPS.md; 2026-11-08).
#
# - seqformat to SPS with a right= state, and with complex or logical cells (as.character).
# - as.character of complex numbers, and strsplit with an empty separator.
# - R's as.numeric on nan/hex spellings, which traminex's r_string_to_double reproduces.
# - seqdef on complex cells: R keeps only the real parts (an R defect; traminex keeps the states).
# - seqdist OM with sm = "TRATE" whose largest cost is 0: the full matrix and refseq index sets.
# - seqdist OM with sm = "INDELS" against an external reference with a void cell.
#
# Values are written with sprintf("%.17g"). Consumed by tests/test_consumed_gaps_part4_construction.py
# and tests/test_consumed_gaps_part4_distances.py.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

num <- function(x) paste(sprintf("%.17g", as.vector(x)), collapse = " ")
rows <- character(0)
add <- function(case, key, value) {
  rows <<- c(rows, paste(shQuote(case, "cmd"), shQuote(key, "cmd"), shQuote(value, "cmd"), sep = ","))
}
spells <- function(sps) {
  # One string per sequence, its spells joined by "-", as traminex returns SPS.
  m <- as.matrix(sps)
  paste(apply(m, 1, function(row) paste(row[!is.na(row)], collapse = "-")), collapse = " | ")
}

## seqformat to SPS.
right_x <- data.frame(t1 = c("A", "B"), t2 = c("B", "B"), t3 = c("C", NA))
add("sps_right_x", "sps", spells(suppressMessages(seqformat(right_x, from = "STS", to = "SPS", right = "X"))))
complex_cells <- data.frame(a = c(1 + 1i, 2 + 0i), b = c(1 + 1i, 1 + 1i))
add("sps_complex", "sps", spells(suppressMessages(seqformat(complex_cells, from = "STS", to = "SPS"))))
logical_cells <- data.frame(a = c(TRUE, FALSE), b = c(TRUE, TRUE))
add("sps_logical", "sps", spells(suppressMessages(seqformat(logical_cells, from = "STS", to = "SPS"))))

## as.character of complex numbers and strsplit with an empty separator.
z <- c(1 + 1i, -1 - 2.5i, 0.1 + 0i, 1e10 + 1i, 0 - 1i, complex(real = 1, imaginary = -0),
       complex(real = Inf, imaginary = -Inf), complex(real = NaN, imaginary = 1), complex(real = 1 / 3, imaginary = -2 / 3))
add("complex_text", "text", paste(as.character(z), collapse = " | "))
add("strsplit_empty", "parts", paste(strsplit("ABC", "")[[1]], collapse = " "))

## as.numeric on R's nan and hexadecimal spellings.
spellings <- c("nan", "NaN", " nan ", "nanx", "0x1G", "0x10", "0x1p5000", "-0x1p5000")
add("as_numeric", "values", num(suppressWarnings(as.numeric(spellings))))

## seqdef on complex cells.
complex_states <- suppressMessages(seqdef(data.frame(a = c(2 + 0i, 1 + 1i, 1 - 1i), b = c(1 + 1i, 1 + 1i, 2 + 0i))))
add("seqdef_complex", "alphabet", paste(alphabet(complex_states), collapse = " "))

## TRATE costs whose largest cost is 0: two states that always alternate.
alternating <- suppressMessages(seqdef(data.frame(a = c("A", "B", "A", "B"), b = c("B", "A", "B", "A"),
                                                  c = c("A", "B", "A", "B"), d = c("B", "A", "B", "A"))))
add("trate_zero_full", "distances", num(suppressMessages(seqdist(alternating, method = "OM", sm = "TRATE"))))
add("trate_zero_refsets", "distances",
    num(suppressMessages(seqdist(alternating, method = "OM", sm = "TRATE", refseq = list(c(1, 2), c(3, 4))))))

## INDELS costs against an external reference with a void cell.
base <- suppressMessages(seqdef(data.frame(a = c("A", "B", "A"), b = c("B", "B", NA), c = c("C", NA, NA)),
                                alphabet = c("A", "B", "C")))
external <- suppressMessages(seqdef(data.frame(a = "A", b = "B", c = NA), alphabet = c("A", "B", "C")))
add("indels_external_void", "distances", num(suppressMessages(seqdist(base, method = "OM", sm = "INDELS", refseq = external))))

writeLines(
  c("# what R does in the cases of part 4 of the consumed-component gaps", '"case","key","value"', rows),
  file.path(out_dir, "gap_part4.csv")
)
