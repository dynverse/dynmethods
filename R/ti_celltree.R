
#' Inferring trajectories with \code{celltree}
#'
#' Arguments passed to this function will be used as default parameters for the method.
#'
#' @param num_topics_lower The lower bound of topics to be fitted in the model.
#' @param num_topics_upper The upper bound of topics to be fitted in the model.
#' @param num_topics The number of topics to fit in the model.
#' @param tot_iter Number of iterations of the LDA inference.
#' @param tolerance Tolerance values of the LDA inference.
#' @param sd_filter Standard-deviation threshold below which genes should be removed from the data.
#' @param absolute_width Distance threshold below which a cell vertex is considered to be attached to a backbone vertex (see paper for more details).
#'   By default, this threshold is computed dynamically, based on the distance distribution for each branch.
#' @param width_scale_factor A scaling factor for the dynamically-computed distance threshold (ignored if absolute_width is provided).
#'   Higher values will result in less branches in the backbone tree, while lower values might lead to a large number of backbone branches.
#' @param outlier_distance_factor Proportion of vertices, out of the total number of vertices divided by the total number of branches,
#'   that can be left at the end of the backbone tree-building algorithm.
#' @param rooting_method Method used to root the backbone tree. Must be one of: ‘null’, ‘longest.path’, ‘center.start.group’ or ‘average.start.group’.
#' ‘longest.path' picks one end of the longest shortest-path between two vertices.
#' 'center.start.group’ picks the vertex in the starting group with lowest mean-square-distance to the others.
#' ‘average.start.group’ creates a new artificial vertex, as the average of all cells in the starting group.
#' ‘null’ picks the best method based on the type of grouping and start group information available.
#'
#' @rdname celltree
#'
#' @include wrapper_create_description.R
abstract_celltree_description <- function(
  method
) {
  method_value <- c(maptpx = "maptpx", gibbs = "Gibbs", vem = "VEM")[[method]]

  common_params <- list(
    makeDiscreteParam(id = "method", values = method_value, default = method_value),
    makeNumericParam(id = "sd_filter", lower = log(.01), upper = log(5.0), default = log(.5), special.vals = list(FALSE), trafo = exp),
    makeDiscreteParam(id = "absolute_width", values = 0, default = 0, tunable = FALSE),
    makeNumericParam(id = "width_scale_factor", lower = log(.1), default = log(1.5), upper = log(100), trafo = exp),
    makeNumericParam(id = "outlier_tolerance_factor", lower = log(.0001), default = log(.1), upper = log(1000), trafo = exp),
    makeDiscreteParam(id = "rooting_method", values = c("longest.path", "center.start.group", "average.start.group", "null"), default = "null")
  )

  par_set <- switch(
    method,
    maptpx = makeParamSet(
      params = c(common_params, list(
        makeIntegerParam(id = "num_topics_lower", lower = 2L, upper = 15L, default = 2L),
        makeIntegerParam(id = "num_topics_upper", lower = 2L, upper = 15L, default = 15L),
        makeNumericParam(id = "tot_iter", lower = log(1e4), upper = log(1e7), default = log(1e6), trafo = function(x) round(exp(x))),
        makeNumericParam(id = "tolerance", lower = log(.001), upper = log(.5), default = log(.05), trafo = exp)
      )),
      forbidden = quote(num_topics_lower > num_topics_upper)
    ),
    gibbs = makeParamSet(
      params = c(common_params, list(
        makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
        makeNumericParam(id = "tot_iter", lower = log(50), upper = log(500), default = log(200), trafo = function(x) round(exp(x))),
        makeNumericParam(id = "tolerance", lower = log(1e-7), upper = log(1e-3), default = log(1e-5), trafo = exp)
      ))
    ),
    vem = makeParamSet(
      params = c(common_params, list(
        makeIntegerParam(id = "num_topics", lower = 2L, default = 4L, upper = 15L),
        makeNumericParam(id = "tot_iter", lower = log(1e4), upper = log(1e7), default = log(1e6), trafo = function(x) round(exp(x))),
        makeNumericParam(id = "tolerance", lower = log(1e-7), upper = log(1e-3), default = log(1e-5), trafo = exp)
      ))
    )
  )

  create_description(
    name = pritt("cellTree with {method}"),
    short_name = pritt("ct{method}"),
    package_loaded = c(),
    package_required = c("cellTree"),
    par_set = par_set,
    properties = c(),
    run_fun = "run_celltree",
    plot_fun = "plot_celltree"
  )
}


#' @rdname celltree
#' @export
description_ctmaptpx <- abstract_celltree_description("maptpx")

#' @rdname celltree
#' @export
description_ctgibbs <- abstract_celltree_description("gibbs")

#' @rdname celltree
#' @export
description_ctvem <- abstract_celltree_description("vem")

#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
run_celltree <- function(
  # transcriptomics data
  expression,

  # prior information
  start_cells = NULL,
  grouping_assignment = NULL,

  # parameters
  method,
  num_topics_lower = NULL,
  num_topics_upper = NULL,
  num_topics = NULL,
  sd_filter,
  tot_iter,
  tolerance,
  absolute_width,
  width_scale_factor,
  outlier_tolerance_factor,
  rooting_method
) {
  requireNamespace("cellTree")

  start_cell <-
    if (!is.null(start_cells)) {
      sample(start_cells, 1)
    } else {
      NULL
    }

  if (rooting_method == "null") {
    rooting_method <- NULL
  }

  if (is.null(num_topics)) {
    num_topics <- seq(num_topics_lower, num_topics_upper)
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # infer the LDA model
  lda_out <- cellTree::compute.lda(
    t(expression) + min(expression) + 1,
    k.topics = num_topics,
    method = method,
    log.scale = FALSE,
    sd.filter = sd_filter,
    tot.iter = tot_iter,
    tol = tolerance
  )

  # put the parameters for the backbones in a list,
  # for adding optional grouping_assignment and (if grouping is given) start group
  backbone_params <- list(
    lda.results = lda_out,
    absolute.width = absolute_width,
    width.scale.factor = width_scale_factor,
    outlier.tolerance.factor = outlier_tolerance_factor,
    rooting.method = rooting_method,
    only.mst = FALSE,
    merge.sequential.backbone = FALSE
  )

  # if these parameters are available, add them to the list
  if(!is.null(grouping_assignment)) {
    backbone_params$grouping <- grouping_assignment %>% slice(match(cell_id, rownames(expression))) %>% pull(group_id)
    if(!is.null(start_cell)) {
      backbone_params$start.group.label <- grouping_assignment %>% filter(cell_id == start_cell) %>% pull(group_id)
    }
  }

  # construct the backbone tree
  mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # simplify sample graph to just its backbone
  cell_graph <- igraph::as_data_frame(mst_tree, "edges") %>%
    select(from, to, length = weight) %>%
    mutate(
      from = rownames(expression)[from],
      to = rownames(expression)[to],
      directed = FALSE
    )
  to_keep <- igraph::V(mst_tree)$is.backbone %>%
    setNames(rownames(expression))

  # extract data for visualisations
  tree <- cellTree:::.compute.tree.layout(mst_tree, ratio = 1)
  vertices <- igraph::as_data_frame(tree, "vertices") %>% as_data_frame()
  edges <- igraph::as_data_frame(tree, "edges") %>% as_data_frame()

  # wrap output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep,
    is_directed = FALSE,
    plot_vertices = vertices,
    plot_edges = edges
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom grDevices rainbow
plot_celltree <- function(prediction) {
  requireNamespace("ggforce")

  # Based on cellTree::ct.plot.topics(prediction$mst_tree)
  vertices <- prediction$plot_vertices
  edges <- prediction$plot_edges

  # calculate pie sizes
  pie_df <- map_df(seq_len(nrow(vertices)), function(i) {
    pieval <- vertices$pie[[i]]
    data.frame(
      vertices[i,] %>% select(-pie),
      topic = paste0("Topic ", seq_along(pieval)),
      stringsAsFactors = FALSE
    ) %>% mutate(
      topic = factor(topic, levels = topic),
      value = pieval,
      arc = value * 2 * pi / sum(value),
      end = cumsum(arc),
      start = end - arc
    )
  })

  # obtain edge positioning
  edges_df <- data.frame(
    edges,
    from = vertices[edges$from,c("x","y")],
    to = vertices[edges$to,c("x","y")]
  )

  # get color scheme
  num_topics <- length(vertices$pie[[1]])
  ann_cols <- setNames(grDevices::rainbow(num_topics), paste0("Topic ", seq_len(num_topics)))

  # make pie graph plot
  g <- ggplot() +
    geom_segment(aes(x = from.x, xend = to.x, y = from.y, yend = to.y), edges_df) +
    ggforce::geom_arc_bar(aes(x0 = x, y0 = y, r0 = 0, r = size*2,
                              start = start, end = end, fill = topic, group = cell.name), data = pie_df, colour = NA) +
    scale_fill_manual(values = ann_cols) +
    scale_size_identity() +
    labs(fill = "Topic") +
    theme(legend.position = c(0.9, 0.125)) +
    guides(fill = guide_legend(ncol = ceiling(num_topics / 8)))
  process_dynplot(g, prediction$id)
}
