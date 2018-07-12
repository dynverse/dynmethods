library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(monocle)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds("/input/data.rds")
params <- jsonlite::read_json("/input/params.json")

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  counts,
  groups_n = NULL,
  reduction_method = "DDRTree",
  max_components = 2,
  norm_method = "vstExprs",
  auto_param_selection = TRUE
) {
  requireNamespace("monocle")
  requireNamespace("BiocGenerics")
  requireNamespace("igraph")
  requireNamespace("Biobase")
  requireNamespace("Matrix")

  # just in case
  if (is.factor(norm_method)) norm_method <- as.character(norm_method)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # load in the new dataset
  pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
  fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
  cds <- monocle::newCellDataSet(t(counts), pd, fd)

  # estimate size factors and dispersions
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds)

  # reduce the dimensionality
  cds <- monocle::reduceDimension(
    cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    auto_param_selection = auto_param_selection
  )

  # order the cells
  cds <- monocle::orderCells(cds, num_paths = groups_n)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract the igraph and which cells are on the trajectory
  gr <-
    if (reduction_method == "DDRTree") {
      cds@auxOrderingData[[reduction_method]]$pr_graph_cell_proj_tree
    } else if (reduction_method == "ICA") {
      cds@auxOrderingData[[reduction_method]]$cell_ordering_tree
    }
  to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

  # convert to milestone representation
  cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)

  if ("weight" %in% colnames(cell_graph)) {
    cell_graph <- cell_graph %>% rename(length = weight)
  } else {
    cell_graph <- cell_graph %>% mutate(length = 1)
  }

  cell_graph <- cell_graph %>% select(from, to, length, directed)

  # retrieve data for visualisation
  if (exists("postprocess_monocle_cds")) {
    plot_data <- postprocess_monocle_cds(cds)
  } else {
    plot_data <- NULL
  }


  # wrap output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep,
    plot_data = plot_data,
    reduction_method = reduction_method
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, "/output/output.rds")