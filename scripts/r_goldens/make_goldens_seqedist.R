# Golden generation: TraMineRextras::seqedist -- the public event-sequence
# distance wrapper. Public fixtures cover every interval and normalization
# string. Internal TraMineR engine fixtures separately preserve Python's
# numeric-code extension.
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
suppressMessages(library(TraMineRextras))
stopifnot(
  as.character(packageVersion("TraMineR")) == "2.2.14",
  as.character(packageVersion("TraMineRextras")) == "0.6.9"
)

out_dir <- file.path("tests", "r_goldens")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tse <- data.frame(
  id    = c("s1", "s1", "s1", "s2", "s2", "s3", "s3", "s3", "s4", "s4"),
  time  = c(0,    2,    2,    0,    3,    0,    1,    4,    0,    2),
  event = c("A",  "B",  "C",  "A",  "C",  "B",  "A",  "C",  "A",  "B"),
  stringsAsFactors = FALSE
)
e <- seqecreate(tse)
idcost <- c(1, 2, 3)  # A, B, C
V <- 0.5

dump_mat <- function(m, tag, call) {
  con <- file(file.path(out_dir, paste0("seqedist_", tag, ".csv")), "w")
  cat("# R: ", call, "\n", sep = "", file = con)
  cat(sprintf("# %d x %d distance matrix\n", nrow(m), ncol(m)), file = con)
  for (i in seq_len(nrow(m))) {
    cat(paste(sprintf("%.15g", m[i, ]), collapse = ","), "\n", file = con, sep = "")
  }
  close(con)
}

dump_mat(
  TraMineRextras::seqedist(e, idcost = idcost, vparam = V),
  "public_default",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V)"
)
dump_mat(
  TraMineRextras::seqedist(
    e, idcost = idcost, vparam = V, interval = "No", norm = "none"
  ),
  "public_no_none",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V, interval='No', norm='none')"
)
dump_mat(
  TraMineRextras::seqedist(
    e, idcost = idcost, vparam = V, interval = "No", norm = "max"
  ),
  "public_no_max",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V, interval='No', norm='max')"
)
dump_mat(
  TraMineRextras::seqedist(
    e, idcost = idcost, vparam = V, interval = "No", norm = "YujianBo"
  ),
  "public_no_yujianbo",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V, interval='No', norm='YujianBo')"
)
dump_mat(
  TraMineRextras::seqedist(
    e, idcost = idcost, vparam = V, interval = "previous", norm = "none"
  ),
  "public_previous_none",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V, interval='previous', norm='none')"
)
dump_mat(
  TraMineRextras::seqedist(
    e, idcost = idcost, vparam = V, interval = "next", norm = "none"
  ),
  "public_next_none",
  "TraMineRextras::seqedist(e, idcost=idcost, vparam=V, interval='next', norm='none')"
)

# Internal engine evidence for Python's documented numeric-code extension.
dump_mat(
  TraMineR:::seqedist(e, idcost = idcost, vparam = V, interval = 1, norm = 0),
  "internal_i1_n0",
  "TraMineR:::seqedist(e, idcost=idcost, vparam=V, interval=1, norm=0)"
)
dump_mat(
  TraMineR:::seqedist(e, idcost = idcost, vparam = V, interval = 1, norm = 1),
  "internal_i1_n1",
  "TraMineR:::seqedist(e, idcost=idcost, vparam=V, interval=1, norm=1)"
)
dump_mat(
  TraMineR:::seqedist(e, idcost = idcost, vparam = V, interval = 1, norm = 2),
  "internal_i1_n2",
  "TraMineR:::seqedist(e, idcost=idcost, vparam=V, interval=1, norm=2)"
)

cat("SEQEDIST_GOLDENS_DONE\n")
