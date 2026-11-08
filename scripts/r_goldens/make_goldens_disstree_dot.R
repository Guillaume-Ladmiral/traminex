# Golden generation: faithful disstree2dot DOT output (RP2 / T20, R2 image route).
#
# Anchors traminex's R-faithful disstree2dot to TraMineR 2.2-14 `disstree2dot`.
# Uses the text/topology path (image.fun=NULL, use.title=FALSE) so node labels
# carry the n / s2 / Split / R2 text and no plot device is needed. The R2 cases
# additionally freeze the semantic per-node image protocol. Raster bytes and R
# device/title-placement details are intentionally not goldened: the callback
# records its subsets while the generated DOT and call manifest record filenames,
# labels, title mode, and depth behavior. Same deterministic tree as
# make_goldens_disstree.R (R=1 -> no permutation gate). The filename "tree" is
# the DOT node-name prefix, so it must match the Python test's filename.
#
# Requires R + TraMineR 2.2-14; see memory "live-traminer-r-env".
# Fixed locale and message language: byte-identical fixtures from any shell (TX-AUD-20261010-019).
stopifnot(nzchar(Sys.setlocale("LC_ALL", "C.UTF-8")), capabilities("ICU")); invisible(Sys.setLanguage("en")); icuSetCollate(locale = "root")
suppressMessages(library(TraMineR))
# Oracle pin: a generator must refuse to run against the wrong TraMineR,
# otherwise a regeneration can silently record the wrong provenance.
stopifnot(as.character(packageVersion("TraMineR")) == "2.2.14")
source(file.path("scripts", "r_goldens", "make_goldens_phase2_inputs.R"))

inp <- phase2_inputs()
D <- inp$D
i <- 1:inp$n
pred <- data.frame(
  A = factor(c("a", "a", "b", "b", "a", "b", "b", "b", "c", "c", "c", "c", "c", "a", "c")),
  X = ((i * 7) %% 15) + 0.1 * i
)
tr <- disstree(D ~ A + X, data = pred, R = 1, min.size = 2, max.depth = 4)

out_dir <- normalizePath(file.path("tests", "r_goldens"), mustWork = TRUE)
td <- tempdir()
oldwd <- getwd()
setwd(td)
disstree2dot(tr, "tree", image.fun = NULL, use.title = FALSE)
file.copy(file.path(td, "tree.dot"),
  file.path(out_dir, "disstree_dot_plain.dot"), overwrite = TRUE)
disstree2dot(tr, "tree", image.fun = NULL, use.title = FALSE, show.depth = TRUE)
file.copy(file.path(td, "tree.dot"),
  file.path(out_dir, "disstree_dot_depth.dot"), overwrite = TRUE)

# R2 portable image protocol -------------------------------------------------

formd <- function(x) format(x, digits = 3)

split_label <- function(sp, side) {
  if (!is.null(sp$breaks)) {
    out <- c(paste("<=", formd(sp$breaks)), paste(">", formd(sp$breaks)))
  } else {
    out <- c(
      paste("[", paste(sp$labels[sp$index == 1], collapse = ","), "]"),
      paste("[", paste(sp$labels[sp$index == 2], collapse = ","), "]")
    )
  }
  if (!is.null(sp$naGroup)) {
    out[sp$naGroup] <- paste(out[sp$naGroup], "with NA")
  }
  out[side]
}

node_manifest <- list()
walk_nodes <- function(node, nodename = "tree_none", branch = "Root") {
  node_text <- paste0("n: ", formd(node$info$n), " s2: ", formd(node$info$vardis))
  split_text <- if (is.null(node$split)) {
    ""
  } else {
    paste0(
      "Split: ", colnames(tr$data)[node$split$varindex],
      " R2:", formd(node$split$info$R2)
    )
  }
  node_manifest[[length(node_manifest) + 1L]] <<- data.frame(
    node_name = nodename,
    is_leaf = is.null(node$split),
    node_text = node_text,
    split_text = split_text,
    branch_label = branch,
    stringsAsFactors = FALSE
  )
  if (!is.null(node$split)) {
    walk_nodes(
      node$kids[[1]], paste0(nodename, "_left"), split_label(node$split, 1)
    )
    walk_nodes(
      node$kids[[2]], paste0(nodename, "_right"), split_label(node$split, 2)
    )
  }
}
walk_nodes(tr$root)
node_manifest <- do.call(rbind, node_manifest)

calls <- list()
run_image_case <- function(case_name, image_data, only_leaf, image_format,
                           use_title, show_depth, with_references = FALSE) {
  selected <- node_manifest[!only_leaf | node_manifest$is_leaf, , drop = FALSE]
  call_no <- 0L
  image_callback <- function(x) {
    call_no <<- call_no + 1L
    current <- selected[call_no, , drop = FALSE]
    values <- if (is.data.frame(x)) x$row_id else as.vector(x)
    calls[[length(calls) + 1L]] <<- data.frame(
      case = case_name,
      call_order = call_no,
      node_name = current$node_name,
      image_filename = paste0(current$node_name, ".", image_format),
      subset_kind = if (is.data.frame(x)) "data.frame" else "vector",
      subset_values = paste(values, collapse = ";"),
      only_leaf = only_leaf,
      image_format = image_format,
      use_title = use_title,
      show_depth = show_depth,
      node_text = current$node_text,
      split_text = current$split_text,
      branch_label = current$branch_label,
      stringsAsFactors = FALSE
    )
    plot.new()
  }

  extra <- if (with_references) {
    list(
      legend.text = "<B>R2 legend</B>",
      image.legend = "legend.png",
      image.quality = "quality.png"
    )
  } else {
    list()
  }
  do.call(disstree2dot, c(
    list(
      tree = tr,
      filename = "tree",
      image.fun = image_callback,
      image.data = image_data,
      only.leaf = only_leaf,
      device = "png",
      image.format = image_format,
      use.title = use_title,
      show.depth = show_depth
    ),
    extra
  ))
}

explicit_data <- data.frame(
  row_id = sprintf("row-%02d", seq_len(inp$n)),
  value = seq_len(inp$n) * 10,
  stringsAsFactors = FALSE
)

# The first case is also the structural image DOT golden. The other cases add
# the leaf/default/title matrix to the callback manifest.
run_image_case(
  "explicit_all", explicit_data, FALSE, "png", FALSE, TRUE,
  with_references = TRUE
)
file.copy(file.path(td, "tree.dot"),
  file.path(out_dir, "disstree_dot_r2_images.dot"), overwrite = TRUE)
run_image_case("explicit_leaf", explicit_data, TRUE, "png", TRUE, FALSE)
run_image_case("default_all", NULL, FALSE, "jpg", TRUE, TRUE)
run_image_case("default_leaf", NULL, TRUE, "jpg", FALSE, FALSE)

call_manifest <- do.call(rbind, calls)
write.csv(
  call_manifest,
  file.path(out_dir, "disstree_dot_r2_image_calls.csv"),
  row.names = FALSE,
  quote = TRUE
)
setwd(oldwd)

cat("=== plain ===\n")
cat(readLines(file.path(out_dir, "disstree_dot_plain.dot")), sep = "\n")
cat("\n=== depth ===\n")
cat(readLines(file.path(out_dir, "disstree_dot_depth.dot")), sep = "\n")
cat("\n=== R2 images ===\n")
cat(readLines(file.path(out_dir, "disstree_dot_r2_images.dot")), sep = "\n")
cat("\n=== R2 image calls ===\n")
print(call_manifest, row.names = FALSE)
cat("\nDISSTREE_DOT_GOLDENS_DONE\n")
