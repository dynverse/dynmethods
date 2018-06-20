library(dynwrap)
library(jsonlite)
library(readr)
library(dplyr)
library(purrr)

library(sincell)

#   ____________________________________________________________________________
#   Load data                                                               ####

data <- read_rds('/input/data.rds')
params <- jsonlite::read_json('/input/params.json')

#   ____________________________________________________________________________
#   Infer trajectory                                                        ####

run_fun <- function(
  expression,
  distance_method = "euclidean",
  dimred_method = "none",
  clust.method = "max.distance",
  mutual = TRUE,
  max.distance = 0,
  k = 3,
  shortest.rank.percent = 10,
  graph.algorithm = "MST",
  graph.using.cells.clustering = FALSE,
  k_imc = 3,
  pct_leaf_node_cutoff = 0.5
) {
  requireNamespace("sincell")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # initialise sincell object
  SO <- sincell::sc_InitializingSincellObject(t(expression))

  # calculate distances
  SO <- SO %>% sincell::sc_distanceObj(
    method = distance_method
  )

  # perform dimred, if necessary
  if (dimred_method != "none") {
    SO <- SO %>% sincell::sc_DimensionalityReductionObj(
      method = dimred_method
    )
  }

  # cluster cells
  SO <- SO %>% sincell::sc_clusterObj(
    clust.method = clust.method,
    mutual = mutual,
    max.distance = max.distance,
    shortest.rank.percent = shortest.rank.percent,
    k = k
  )

  # build graph
  SO <- SO %>% sincell::sc_GraphBuilderObj(
    graph.algorithm = graph.algorithm,
    graph.using.cells.clustering = graph.using.cells.clustering,
    k = k_imc
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract cell hierarchy
  cell_graph <- SO$cellstateHierarchy %>%
    igraph::as_data_frame() %>%
    rename(length = weight) %>%
    mutate(directed = FALSE)

  # Leaf nodes are iteratively removed until the percentage of leaf nodes
  # is below the given cutoff. Removed nodes are projected to their closest
  # neighbour.
  # This is to constrain the number of milestones being created.
  gr <- SO$cellstateHierarchy
  deg <- igraph::degree(gr)
  prev_deg <- deg * 0
  while (length(deg) > 10 && mean(deg <= 1) > pct_leaf_node_cutoff && any(deg != prev_deg)) {
    del_v <- names(which(deg == 1))
    cat("Removing ", length(del_v), " vertices with degree 1\n", sep = "")
    gr <- igraph::delete_vertices(gr, del_v)
    prev_deg <- deg[names(igraph::V(gr))]
    deg <- igraph::degree(gr)
  }
  to_keep <- setNames(rownames(expression) %in% names(igraph::V(gr)), rownames(expression))

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

args <- params[intersect(names(params), names(formals(run_fun)))]

model <- do.call(run_fun, c(args, data))

#   ____________________________________________________________________________
#   Save output                                                             ####

write_rds(model, '/output/output.rds')