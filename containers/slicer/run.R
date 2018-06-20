library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(SLICER)
library(lle)
library(igraph)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  start_id,
  features_id = NULL,
  end_id = NULL,
  kmin = 10,
  m = 2,
  verbose = FALSE
) {
  requireNamespace("SLICER")
  requireNamespace("lle")
  requireNamespace("igraph")

  start_cell <- sample(start_id, 1)

  if (!is.null(features_id)) {
    # use 'neighbourhood variance' to identify genes that vary smoothly
    features_id <- SLICER::select_genes(expression)
  }
  expr_filt <- expression[, features_id]

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # try catch for sink
  tryCatch({
    # stop output if not verbose
    if (!verbose) {
      sink("/dev/null")
    }

    # determine k for knn
    k <- SLICER::select_k(expr_filt, kmin = kmin)

    # perform local linear embedding
    traj_lle <- lle::lle(expr_filt, m = m, k = k)$Y
    rownames(traj_lle) <- rownames(expr_filt)
    colnames(traj_lle) <- paste0("comp_", seq_len(ncol(traj_lle)))

  }, finally = {
    # resume output if not verbose
    if (!verbose) {
      sink()
    }
  })

  # get LLE KNN graph
  traj_graph <- SLICER::conn_knn_graph(traj_lle, k = k)

  # find extreme cells
  if (is.null(end_id)) {
    ends <- SLICER::find_extreme_cells(traj_graph, traj_lle, do_plot = FALSE)
  } else {
    ends <- match(c(start_cell, end_id), rownames(expression))
  }

  # order cells
  start <- which(rownames(expr_filt) == start_cell)
  cells_ordered <- SLICER::cell_order(traj_graph, start)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get shortest paths to start and all other nodes
  shortest_paths <- igraph::shortest_paths(traj_graph, start)
  gr <- lapply(shortest_paths$vpath, function(path) {
    P <- rbind(path[-length(path)], path[-1]) %>% as.vector
    igraph::E(traj_graph, P = P)
  })
  subgr <- igraph::subgraph.edges(traj_graph, eids = unique(unlist(gr)))

  # prepare sample graph simplification
  cell_graph <- igraph::as_data_frame(subgr, "edges") %>%
    dplyr::select(from, to, length = weight) %>%
    mutate(
      from = rownames(expr_filt)[from],
      to = rownames(expr_filt)[to],
      directed = TRUE
    )
  sh_p_to_ends <- igraph::shortest_paths(subgr, start, ends)
  nodes_to_keep <- unique(sh_p_to_ends$vpath %>% unlist)
  to_keep <- setNames(igraph::V(traj_graph) %in% nodes_to_keep, rownames(expr_filt))

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expr_filt)
  ) %>% add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep,
    traj_graph = subgr,
    start = start,
    ends = ends,
    is_kept = to_keep
  ) %>% add_dimred(
    dimred = traj_lle
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')