#' Description for sincell
#' @export
description_sincell <- function() create_description(
  name = "sincell",
  short_name = "sincell",
  package_loaded = c(),
  package_required = c("sincell"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "distance_method", default = "euclidean", values = c("euclidean", "cosine", "pearson", "spearman", "L1", "MI")),
    makeDiscreteParam(id = "dimred_method", default = "none", values = c("none", "PCA", "ICA", "tSNE", "classical-MDS", "nonmetric-MDS")),
    makeDiscreteParam(id = "cluster_method", default = "max.distance", values = c("max.distance", "percent", "knn", "k-medoids", "ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid")),
    makeLogicalParam(id="mutual", default = TRUE, requires = quote(cluster_method == "knn")),
    makeNumericParam(id="max.distance", default = 0, requires = quote(cluster_method == "max.distance"), lower=0, upper=5),
    makeIntegerParam(id="k", default = 3L, requires = quote(cluster_method == "knn"), lower=1L, upper=99L),
    makeNumericParam(id="shortest.rank.percent", default = 10, requires = quote(cluster_method == "percent"), lower=0, upper=100),
    makeDiscreteParam(id = "graph_method", default = "MST", values = c("MST", "SST", "IMC")),
    makeLogicalParam(id = "graph.using.cells.clustering", default = FALSE, requires = quote(graph_method == "MST")),
    makeIntegerParam(id = "k_imc", default = 3L, requires = quote(graph_method == "IMC"), lower=1L, upper=99L)
    # makeIntegerParam(id = "ndim", lower = 2L, default = 3L, upper = 20L),
    # makeIntegerParam(id = "k", lower = 0L, default = 4L, upper = 20L),
    # makeNumericParam(id = "thresh", lower = -5, upper = 5, default = -3, trafo = function(x) 10^x),
    # makeIntegerParam(id = "maxit", lower = 0L, upper = 50L, default = 10L),
    # makeNumericParam(id = "stretch", lower = 0, upper = 5, default = 0),
    # makeDiscreteParam(id = "smoother", default = "smooth.spline", values = c("smooth.spline", "lowess", "periodic.lowess"))
  ),
  properties = c(),
  run_fun = run_sincell,
  plot_fun = plot_sincell
)

run_sincell <- function(
  expression,
  distance_method = "spearman",
  dimred_method = "none",
  cluster_method = "max.distance",
  mutual = TRUE,
  max.distance = 0,
  k = 3L,
  shortest.rank.percent = 10,
  graph_method = "MST",
  graph.using.cells.clustering = FALSE,
  k_imc = 3L
  ) {
  requireNamespace("sincell")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # initialise sincell object
  SO <- sincell::sc_InitializingSincellObject(t(expression))

  # calculate distances
  SO <- sincell::sc_distanceObj(SO, distance_method)

  # perform dimred, if necessary
  if (dimred_method != "none") {
    SO <- sincell::sc_DimensionalityReductionObj(SO, dimred_method)
  }

  cell2celldist <- SO$cell2celldist

  # cluster cells
  SO <- sincell::sc_clusterObj(SO, cluster_method, mutual, max.distance, shortest.rank.percent, k)

  # build graph
  SO <- sincell::sc_GraphBuilderObj(SO, graph_method, graph.using.cells.clustering, k_imc)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract cell hierarchy
  edges <- SO$cellstateHierarchy %>%
    igraph::as_data_frame() %>%
    rename(length = weight) %>%
    mutate(directed = FALSE)

  # keep all cells
  to_keep <- setNames(rep(TRUE, nrow(expression)), rownames(expression))

  # simplify sample graph
  out <- dynutils::simplify_sample_graph(edges, to_keep, is_directed=FALSE)

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # remove extra data in SO for visualisation
  SO$expressionmatrix <- NULL
  SO$cell2celldist <- NULL

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    SO = SO
  ) %>% attach_timings_attribute(tl)
}

#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr set_colnames
plot_sincell<- function(prediction) {
  requireNamespace("sincell")
  requireNamespace("ggraph")

  prediction$SO$cellstateHierarchy %>%
    tidygraph::as_tbl_graph(SO$cellstateHierarchy) %>%
    ggraph::ggraph() +
    ggraph::geom_node_point() +
    ggraph::geom_edge_link() +
    ggraph::theme_graph()
}
