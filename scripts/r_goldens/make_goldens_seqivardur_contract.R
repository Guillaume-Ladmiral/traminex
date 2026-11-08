# Golden generation for the full TraMineR::seqivardur result contract.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_contract <- function(result, name, call) {
  values <- rbind(
    variance = as.vector(result),
    vmax = as.vector(attr(result, "vmax")),
    meand = as.vector(attr(result, "meand"))
  )
  con <- file(file.path(out_dir, paste0(name, ".csv")), open = "wt")
  on.exit(close(con), add = TRUE)
  cat("# R: ", call, "\n", sep = "", file = con)
  cat("# rows: variance,vmax,meand\n", file = con)
  text <- apply(
    values,
    1,
    function(row) paste(formatC(row, format = "g", digits = 15), collapse = ",")
  )
  writeLines(text, con)
}

A <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "B", "C", "A"),
  c("D", "D", "D", "D", "D", "D"),
  c("A", "B", "A", "B", "A", "B")
)
sdA <- seqdef(A, alphabet = c("A", "B", "C", "D"))

C <- rbind(
  c("A", "A", "B", "C", "C", "C"),
  c("A", "B", "B", "*", "*", "*"),
  c("A", "*", "A", "B", "*", "*"),
  c("C", "C", "*", "*", "*", "*")
)
sdC <- seqdef(C, alphabet = c("A", "B", "C"), missing = "*")

write_contract(
  seqivardur(sdA, type = 1, with.missing = FALSE),
  "seqivardur_contract_A_type1",
  "seqivardur(sdA, type=1, with.missing=FALSE)"
)
write_contract(
  seqivardur(sdA, type = 2, with.missing = FALSE),
  "seqivardur_contract_A_type2",
  "seqivardur(sdA, type=2, with.missing=FALSE)"
)
write_contract(
  seqivardur(sdC, type = 1, with.missing = TRUE),
  "seqivardur_contract_C_wm1",
  "seqivardur(sdC, type=1, with.missing=TRUE)"
)

# A row that is all structural void is supported inside a mixed object, but an
# object made entirely from that row is refused by TraMineR's seqdss path.
V <- rbind(
  c("A", "B", "A", "B"),
  c(NA, NA, NA, NA)
)
sdV <- seqdef(V, alphabet = c("A", "B"))
mixed_type1 <- seqivardur(sdV, type = 1)
mixed_type2 <- seqivardur(sdV, type = 2)
all_void <- sdV[2, , drop = FALSE]
oracle_error <- tryCatch(seqivardur(all_void), error = identity)
stopifnot(inherits(oracle_error, "error"))
scalar_text <- function(value) trimws(formatC(value, format = "g", digits = 15))
error_lines <- c(
  paste0("class=", paste(class(oracle_error), collapse = ",")),
  paste0("message=", conditionMessage(oracle_error)),
  paste0("mixed_type1_variance=", scalar_text(mixed_type1[2])),
  paste0("mixed_type1_vmax=", scalar_text(attr(mixed_type1, "vmax")[2])),
  paste0("mixed_type1_meand=", scalar_text(attr(mixed_type1, "meand")[2])),
  paste0("mixed_type2_variance=", scalar_text(mixed_type2[2])),
  paste0("mixed_type2_vmax=", scalar_text(attr(mixed_type2, "vmax")[2])),
  paste0("mixed_type2_meand=", scalar_text(attr(mixed_type2, "meand")[2]))
)
writeLines(error_lines, file.path(out_dir, "seqivardur_contract_errors.txt"))

cat("SEQIVARDUR_CONTRACT_GOLDENS_DONE\n")
