# Live-R goldens for the text R writes into a discrepancy-tree node (repair of
# audit finding TX-AUD-20261010-014).
#
# disstree2dot draws a node image with image.fun and writes the node text into
# that image's title (disstree2dot.R:365-383): main is the branch label ("Root",
# "<= x", "> x", "[ a,b ]") and the node line "n: .. s2: ..", sub is the split
# line "Split: <variable> R2:..". The DOT label of such a node is " ". traminex
# draws no node image unless asked, so it writes this text as the DOT label.
# The title() calls are recorded by tracing graphics::title while disstree2dot
# runs with an image function that draws an empty plot (pdf device). Same tree
# as make_goldens_disstree_dot.R.
#
# Run from the repository root:
#   R_LIBS=<oracle library> Rscript scripts/r_goldens/make_goldens_disstree_node_labels.R

# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
icuSetCollate(locale = "root")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D
i <- 1:inp$n
pred <- data.frame(
  A = factor(c("a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c")),
  X = ((i * 7) %% 15) + 0.1 * i
)
tr <- TraMineR::disstree(D ~ A + X, data = pred, R = 1, min.size = 2, max.depth = 4)

out_dir <- normalizePath(file.path("tests", "r_goldens"), mustWork = TRUE)
titles <- list()
trace("title", tracer = quote(titles[[length(titles) + 1L]] <<- list(main = main, sub = sub)),
      print = FALSE, where = asNamespace("TraMineR"))
td <- tempfile("nodelabels")
dir.create(td)
oldwd <- setwd(td)
TraMineR::disstree2dot(tr, "tree", image.fun = function(...) plot.new(), device = "pdf",
                       image.format = "pdf")
setwd(oldwd)
untrace("title", where = asNamespace("TraMineR"))
dot <- readLines(file.path(td, "tree.dot"))
names <- sub('^"([^"]+)"\\[shape=box.*', "\\1", grep("shape=box", dot, value = TRUE))
names <- names[!startsWith(names, "node_")]
stopifnot(length(names) == length(titles))
rows <- data.frame(position = seq_along(titles), node_name = names,
                   main = vapply(titles, function(t) if (is.null(t$main)) "" else t$main, ""),
                   sub = vapply(titles, function(t) if (is.null(t$sub)) "" else t$sub, ""),
                   stringsAsFactors = FALSE)
write.csv(rows, file.path(out_dir, "disstree_node_labels.csv"), row.names = FALSE, quote = TRUE)
