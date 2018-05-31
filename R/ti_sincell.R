#' Inferring trajectories with Sincell
#'
#' @inherit ti_angle description
#'
#' @param distance_method Distance method to be used. The available distances are the Euclidean distance (method="euclidean"), Manhattan distance (also called L1 distance, method="L1"), cosine distance (method="cosine") , distance based on Pearson (method="pearson") or Spearman (method="spearman") correlation coefficients, and distance based on Mutual Information (method="MI"). Intervals used to assess Mutual Information are indicated in the parameter “bins”.
#' @param dimred_method Dimensionality reduction algorithm to be used. Options are: Principal Component Analysis (method="PCA"), Independent Component Analysis (method="ICA"; using fastICA() function in fastICA package), t-Distributed Stochastic Neighbor Embedding (method="tSNE"; using Rtsne() function in Rtsne package with parameters tsne.perplexity=1 and tsne.theta=0.25), classical Multidimensional Scaling (method="classical-MDS"; using the cmdscale() function) and non-metric Multidimensional Scaling (method="nonmetric-MDS";using the isoMDS() function in MASS package). if method="PCA" is chosen, the proportion of variance explained by each of the principal axes is plotted. We note that Sincell makes use of the Rtsne implementation of the Barnes-Hut algorithm, which approximates the likelihood. The user should be aware that this is a less accurate version of t-SNE than e.g. the one used as basis of viSNE (Amir,E.D. et al. 2013, Nat Biotechnol 31, 545–552).
#' @inheritParams sincell::sc_clusterObj
#' @inheritParams sincell::sc_GraphBuilderObj
#' @param k_imc If IMC algorithm is selected, the number of nearest neighbors used in the underlying K-Mutual Nearest Neighbour (K-MNN) algorithm is set to k.
#' @param pct_leaf_node_cutoff Leaf nodes are iteratively removed until the percentage of leaf nodes is below the given cutoff. Removed nodes are projected to their closest neighbour. This is to constrain the number of milestones being created.
#'
#' @seealso [sincell::sc_distanceObj()], [sincell::sc_DimensionalityReductionObj()], [sincell::sc_clusterObj()]
#'
#' @export
ti_sincell <- create_ti_method(
  name = "Sincell",
  short_name = "sincell",
  package_loaded = c(),
  package_required = c("sincell"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distance_method", default = "euclidean", values = c("euclidean", "cosine", "pearson", "spearman", "L1", "MI")),
    makeDiscreteParam(id = "dimred_method", default = "none", values = c("none", "PCA", "ICA", "tSNE", "classical-MDS", "nonmetric-MDS")),
    makeDiscreteParam(id = "clust.method", default = "max.distance", values = c("max.distance", "percent", "knn", "k-medoids", "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
    makeLogicalParam(id = "mutual", default = TRUE), #, requires = quote(cluster_method == "knn")),
    makeNumericParam(id = "max.distance", default = 0, lower = 0, upper = 5), #, requires = quote(cluster_method == "max.distance")),
    makeIntegerParam(id = "k", default = 3L, lower=1L, upper=99L), #, requires = quote(cluster_method == "knn")),
    makeNumericParam(id = "shortest.rank.percent", default = 10, lower=0, upper=100), #, requires = quote(cluster_method == "percent")),
    makeDiscreteParam(id = "graph.algorithm", default = "MST", values = c("MST", "SST", "IMC")),
    makeLogicalParam(id = "graph.using.cells.clustering", default = FALSE), #, requires = quote(graph_method == "MST")),
    makeIntegerParam(id = "k_imc", default = 3L, lower=1L, upper=99L), #, requires = quote(graph_method == "IMC")),
    makeNumericParam(id = "pct_leaf_node_cutoff", default = .5, lower = .01, upper = .8) #, requires = quote(!graph_using_cells_clustering)),
  ),
  run_fun = "dynmethods::run_sincell",
  plot_fun = "dynmethods::plot_sincell"
)

run_sincell <- function(
  expression,
  distance_method = "spearman",
  dimred_method = "none",
  clust.method = "max.distance",
  mutual = TRUE,
  max.distance = 0,
  k = 3L,
  shortest.rank.percent = 10,
  graph.algorithm = "MST",
  graph.using.cells.clustering = FALSE,
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
