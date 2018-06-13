#' Inferring trajectories with Growing Neural Gas
#'
#' @inherit ti_angle description
#'
#' @param dimred A character vector specifying which dimensionality reduction method to use.
#'   See [dyndimred::dimred] for the list of available dimensionality reduction methods.
#' @inheritParams dyndimred::dimred
#' @inheritParams GNG::gng
#' @param apply_mst If true, an MST post-processing of the GNG is performed.
#'
#' @export
ti_gng <- create_ti_method(
  name = "Growing Neural Gas",
  short_name = "gng",
  package_loaded = c(),
  package_required = c("GNG", "igraph", "dyndimred"),
  parameters = list(
    dimred = list(
      type = "discrete",
      default = "pca",
      values = c("pca", "mds", "tsne", "ica", "lle", "mds_sammon", "mds_isomds", "mds_smacof", "umap"),
      description = "A character vector specifying which dimensionality reduction method to use.\nSee \\link[dyndimred:dimred]{dyndimred::dimred} for the list of available dimensionality reduction methods."
    ),
    ndim = list(
      type = "integer",
      default = 5L,
      upper = 10L,
      lower = 2L,
      description = "The number of dimensions"
    ),
    max_iter = list(
      type = "numeric",
      default = 13.8155105579643,
      upper = 18.4206807439524,
      lower = 4.60517018598809,
      description = "The max number of iterations"
    ),
    max_nodes = list(
      type = "integer",
      default = 8L,
      upper = 30L,
      lower = 2L,
      description = "The maximum number of nodes"
    ),
    apply_mst = list(
      type = "logical",
      default = TRUE,
      values = c("TRUE", "FALSE"),
      description = "If true, an MST post-processing of the GNG is performed.")
  ),
  run_fun = "dynmethods::run_gng",
  plot_fun = "dynmethods::plot_gng",
  apt_dependencies = "libudunits2-dev"
)

#' @importFrom stats dist
run_gng <- function(
  expression,
  dimred,
  ndim,
  max_iter,
  max_nodes,
  apply_mst
) {
  requireNamespace("GNG")
  requireNamespace("igraph")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # perform dimensionality reduction
  space <- dyndimred::dimred(expression, method = dimred, ndim = ndim)

  # calculate GNG
  gng_out <- GNG::gng(
    space,
    max_iter = max_iter,
    max_nodes = max_nodes,
    assign_cluster = FALSE
  )
  node_dist <- stats::dist(gng_out$node_space) %>% as.matrix

  # transform to milestone network
  node_names <- gng_out$nodes %>% mutate(name = as.character(name))
  milestone_network <- gng_out$edges %>%
    select(from = i, to = j) %>%
    mutate(
      length = node_dist[cbind(from, to)],
      directed = FALSE
    ) %>%
    select(from, to, length, directed)

  # apply MST, if required
  if (apply_mst) {
    gr <- igraph::graph_from_data_frame(milestone_network, directed = F, vertices = node_names$name)
    milestone_network <- igraph::minimum.spanning.tree(gr, weights = igraph::E(gr)$length) %>% igraph::as_data_frame()
  }

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_dimred_projection(
    milestone_ids = rownames(gng_out$node_space),
    milestone_network = milestone_network,
    dimred_milestones = gng_out$node_space,
    dimred = space
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom viridis scale_colour_viridis
plot_gng <- function(prediction) {
  dimred <- data.frame(prediction$dimred) %>% rownames_to_column("cell_id") %>%
    left_join(dynwrap::get_cell_grouping(prediction$milestone_percentages), by = "cell_id")
  dimred_milestones <- data.frame(prediction$dimred_milestones) %>% rownames_to_column("milestone_id")
  segments <- data.frame(prediction$dimred_trajectory_segments)

  g <- ggplot() +
    geom_segment(aes(x = from_comp_1, xend = to_comp_1, y = from_comp_2, yend = to_comp_2), segments, colour = "darkgray") +
    geom_point(aes(comp_1, comp_2, colour = group_id), dimred) +
    geom_point(aes(comp_1, comp_2, colour = milestone_id), dimred_milestones, size = 10, shape = 8) +
    theme(legend.position = "none") +
    labs(x = "comp_1", y = "comp_2")

  if (length(prediction$milestone_ids) <= 9) {
    g <- g + scale_colour_brewer(palette = "Dark2")
  }
  process_dynplot(g, prediction$id)
}
