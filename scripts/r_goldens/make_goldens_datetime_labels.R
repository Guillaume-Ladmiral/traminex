# Live-R goldens for time labels made from dates and date-times (repair of
# audit finding TX-AUD-20261010-018, seqdef_from_long's time axis).
#
# A long table with Date or POSIXct times reaches TraMineR as a wide table whose
# column names are as.character of the times, and seqdef keeps them as the
# sequence column names. as.character gives the date alone for a Date and for a
# date-time at midnight, element by element, and the clock time otherwise
# (fractional seconds without trailing zeros), never the time zone.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_datetime_labels.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cases <- list(
  dates = as.Date(c("2001-01-01", "2002-01-01", "2005-01-01")),
  midnight_utc = as.POSIXct(c("2001-01-01 00:00:00", "2002-01-01 00:00:00", "2005-01-01 00:00:00"), tz = "UTC"),
  mixed_utc = as.POSIXct(c("2001-01-01 00:00:00", "2002-01-01 12:30:00", "2005-01-01 12:30:00.5"), tz = "UTC"),
  midnight_paris = as.POSIXct(c("2001-01-01 00:00:00", "2002-01-01 00:00:00", "2005-07-01 00:00:00"),
                              tz = "Europe/Paris")
)
rows <- list()
for (name in names(cases)) {
  times <- cases[[name]]
  wide <- as.data.frame(matrix(c("A", "B", "A"), nrow = 1, dimnames = list(NULL, as.character(times))))
  seqs <- suppressMessages(TraMineR::seqdef(wide))
  zone <- if (inherits(times, "POSIXct")) attr(times, "tzone") else ""
  rows[[length(rows) + 1L]] <- data.frame(case = name, kind = class(times)[1], tz = zone,
                                          position = seq_along(times),
                                          time = format(times, "%Y-%m-%d %H:%M:%OS6"),
                                          label = colnames(seqs), stringsAsFactors = FALSE)
}
write.csv(do.call(rbind, rows), file.path(out_dir, "datetime_labels.csv"), row.names = FALSE, quote = TRUE)
