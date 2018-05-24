#' Inferring trajectories with Sincell
#'
#' @inherit ti_angle description
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_sincell <- create_ti_method(
  name = "Sincell",
  short_name = "sincell",
  package_loaded = c(),
  package_required = c("sincell"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distance_method", default = "euclidean", values = c("euclidean", "cosine", "pearson", "spearman", "L1", "MI")),
    makeDiscreteParam(id = "dimred_method", default = "none", values = c("none", "PCA", "ICA", "tSNE", "classical-MDS", "nonmetric-MDS")),
    makeDiscreteParam(id = "cluster_method", default = "max.distance", values = c("max.distance", "percent", "knn", "k-medoids", "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
    makeLogicalParam(id = "mutual", default = TRUE), #, requires = quote(cluster_method == "knn")),
    makeNumericParam(id = "max_distance", default = 0, lower = 0, upper = 5), #, requires = quote(cluster_method == "max.distance")),
    makeIntegerParam(id = "k", default = 3L, lower=1L, upper=99L), #, requires = quote(cluster_method == "knn")),
    makeNumericParam(id = "shortest_rank_percent", default = 10, lower=0, upper=100), #, requires = quote(cluster_method == "percent")),
    makeDiscreteParam(id = "graph_method", default = "MST", values = c("MST", "SST", "IMC")),
    makeLogicalParam(id = "graph_using_cells_clustering", default = FALSE), #, requires = quote(graph_method == "MST")),
    makeIntegerParam(id = "k_imc", default = 3L, lower=1L, upper=99L), #, requires = quote(graph_method == "IMC")),
    makeNumericParam(id = "pct_leaf_node_cutoff", default = .5, lower = .01, upper = .8) #, requires = quote(!graph_using_cells_clustering)),
  ),
  run_fun = "run_sincell",
  plot_fun = "plot_sincell"
)

run_sincell <- function(
  expression,
  distance_method = "spearman",
  dimred_method = "none",
  cluster_method = "max.distance",
  mutual = TRUE,
  max_distance = 0,
  k = 3L,
  shortest_rank_percent = 10,
  graph_method = "MST",
  graph_using_cells_clustering = FALSE,
  k_imc = 3L,
  pct_leaf_node_cutoff = .5
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
    clust.method = cluster_method,
    mutual = mutual,
    max.distance = max_distance,
    shortest.rank.percent = shortest_rank_percent,
    k = k
  )

  # build graph
  SO <- SO %>% sincell::sc_GraphBuilderObj(
    graph.algorithm = graph_method,
    graph.using.cells.clustering = graph_using_cells_clustering,
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

#' @importFrom magrittr set_colnames
plot_sincell<- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  prediction$milestone_network %>%
    igraph::graph_from_data_frame(directed = F) %>%
    tidygraph::as_tbl_graph() %>%
    ggraph::ggraph() +
    ggraph::geom_node_point() +
    ggraph::geom_edge_link() +
    ggraph::theme_graph()
}
