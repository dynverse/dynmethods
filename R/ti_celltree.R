abstract_celltree_description <- function(method) {
  method_value <- c(maptpx = "maptpx", gibbs = "Gibbs", vem = "VEM")[[method]]

  common_params <- list(
    method = list(
      type = "discrete",
      default = method_value,
      values = method_value
    ),
    sd_filter = list(
      type = "numeric",
      lower = 0.01,
      upper = 5,
      default = 0.5
    ),
    absolute_width = list(
      type = "numeric",
      values = 0,
      default = 0,
      tunable = FALSE
    ),
    width_scale_factor = list(
      type = "numeric",
      lower = 0.1,
      default = 1.5,
      upper = 100
    ),
    outlier_tolerance_factor = list(
      type = "numeric",
      lower = 0.0001,
      default = 0.1,
      upper = 1000,
      distribution = "exponential",
      rate = 1
    ),
    rooting_method = list(
      type = "discrete",
      values = c("longest.path", "center.start.group", "average.start.group", "null"),
      default = "null"
    )
  )

  parameters <- switch(
    method,
    maptpx = c(common_params, list(
        num_topics_lower = list(
          type = "integer",
          lower = 2,
          upper = 15,
          default = 2
        ),
        num_topics_upper = list(
          type = "integer",
          lower = 2,
          upper = 15,
          default = 15
        ),
        tot_iter = list(
          type = "numeric",
          lower = 1e4,
          upper = 1e7,
          default = 1e6
        ),
        tolerance = list(
          type = "numeric",
          lower = 0.001,
          upper = 0.5,
          default = 0.05
        ),
        forbidden = "num_topics_lower > num_topics_upper"
      )
    ),
    gibbs = parameters <- c(common_params, list(
        num_topics = list(
          type = "integer",
          lower = 2,
          default = 4,
          upper = 15
        ),
        tot_iter = list(
          type = "numeric",
          lower = 50,
          upper = 500,
          default = 200
        ),
        tolerance = list(
          type = "numeric",
          lower = 1e-7,
          upper = 1e-3,
          default = 1e-5
        )
      )
    ),
    vem = parameters <- c(common_params, list(
      num_topics = list(
        type = "integer",
        lower = 2,
        default = 4,
        upper = 15
      ),
      tot_iter = list(
        type = "numeric",
        lower = 1e4,
        upper = 1e7,
        default = 1e6
      ),
      tolerance = list(
        type = "numeric",
        lower = 1e-7,
        upper = 1e-3,
        default = 1e-5
      )
    ))
  )

  create_ti_method(
    name = pritt("cellTree with {method}"),
    short_name = pritt("celltree_{method}"),
    implementation_id = "celltree",
    package_loaded = c(),
    package_required = c("cellTree"),
    parameters = parameters,
    run_fun = "dynmethods::run_celltree",
    plot_fun = "dynmethods::plot_celltree",
    apt_dependencies = "libgsl-dev",
    doi = "10.1186/s12859-016-1175-6",
    trajectory_types = c("linear", "bifurcation", "convergence", "multifurcation", "binary_tree", "tree"),
    topology_inference = "free",
    type = "algorithm",
    license = "Artistic-2.0",
    authors = list(
      list(
        given = "David",
        family = "duVerle",
        email = "dave@cb.k.u-tokyo.ac.jp",
        role = "aut"
      ),
      list(
        given = "Koji",
        family = "Tsuda",
        email = "tsuda@k.u-tokyo.ac.jp",
        role = "aut"
      )
    ),
    publication_date = "2016-08-13",
    version = "1.10.0"
  )
}

#' Inferring trajectories with cellTree
#'
#' @inherit ti_angle description
#'
#' @param method LDA inference method to use. Can be any unique prefix of ‘maptpx’, ‘Gibbs’ or ‘VEM’ (defaults to ‘maptpx’)
#' @param num_topics_lower The lower bound of topics to be fitted in the model.
#' @param num_topics_upper The upper bound of topics to be fitted in the model.
#' @param tot_iter Number of iterations of the LDA inference.
#' @param tolerance Tolerance values of the LDA inference.
#' @param sd_filter Standard-deviation threshold below which genes should be removed from the data.
#' @param absolute_width Distance threshold below which a cell vertex is considered to be attached to a backbone vertex (see paper for more details).
#'   By default, this threshold is computed dynamically, based on the distance distribution for each branch.
#' @param width_scale_factor A scaling factor for the dynamically-computed distance threshold (ignored if absolute_width is provided).
#'   Higher values will result in less branches in the backbone tree, while lower values might lead to a large number of backbone branches.
#' @param outlier_tolerance_factor Proportion of vertices, out of the total number of vertices divided by the total number of branches,
#'   that can be left at the end of the backbone tree-building algorithm.
#' @param rooting_method Method used to root the backbone tree. Must be one of: ‘null’, ‘longest.path’, ‘center.start.group’ or ‘average.start.group’.
#' ‘longest.path' picks one end of the longest shortest-path between two vertices.
#' 'center.start.group’ picks the vertex in the starting group with lowest mean-square-distance to the others.
#' ‘average.start.group’ creates a new artificial vertex, as the average of all cells in the starting group.
#' ‘null’ picks the best method based on the type of grouping and start group information available.
#'
#' @export
ti_celltree_maptpx <- abstract_celltree_description("maptpx")

#' @inheritParams ti_celltree_maptpx
#' @param num_topics The number of topics to fit in the model.
#' @export
ti_celltree_gibbs <- abstract_celltree_description("gibbs")

#' @inheritParams ti_celltree_maptpx
#' @param num_topics The number of topics to fit in the model.
#' @export
ti_celltree_vem <- abstract_celltree_description("vem")

#' @importFrom igraph degree distances get.vertex.attribute induced_subgraph
run_celltree <- function(
  # transcriptomics data
  expression,

  # prior information
  start_id = NULL,
  groups_id = NULL,

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
    if (!is.null(start_id)) {
      sample(start_id, 1)
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
  # for adding optional groups_id and (if grouping is given) start group
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
  if(!is.null(groups_id)) {
    backbone_params$grouping <- groups_id %>% slice(match(cell_id, rownames(expression))) %>% pull(group_id)
    if(!is.null(start_cell)) {
      backbone_params$start.group.label <- groups_id %>% filter(cell_id == start_cell) %>% pull(group_id)
    }
  }

  # construct the backbone tree
  mst_tree <- do.call(cellTree::compute.backbone.tree, backbone_params)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # simplify sample graph to just its backbone
  cell_graph <- igraph::as_data_frame(mst_tree, "edges") %>%
    dplyr::select(from, to, length = weight) %>%
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
      vertices[i,] %>% dplyr::select(-pie),
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
