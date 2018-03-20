#' Description for gng
#' @export
description_gng <- function() create_description(
  name = "Growing Neural Gas",
  short_name = "gng",
  package_loaded = c(),
  package_required = c("GNG", "igraph"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "dimred", default = "pca", values = names(list_dimred_methods())),
    makeIntegerParam(id = "ndim", default = 5L, lower = 2L, upper = 10L),
    makeNumericParam(id = "max_iter", lower = log(1e2), default = log(1e6), upper = log(1e8), trafo = function(x) round(exp(x))),
    makeIntegerParam(id = "max_nodes", default = 8L, lower = 2L, upper = 30L),
    makeLogicalParam(id = "apply_mst", default = TRUE)
  ),
  properties = c(),
  run_fun = run_gng,
  plot_fun = plot_gng
)

#' @importFrom stats dist
run_gng <- function(
  expression,
  dimred = "pca",
  ndim = 5,
  max_iter = 1e6,
  max_nodes = 8,
  apply_mst = TRUE
) {
  requireNamespace("GNG")
  requireNamespace("igraph")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  space <- list_dimred_methods()[[dimred]](expression, ndim)

  gng_out <- GNG::gng(
    space,
    max_iter = max_iter,
    max_nodes = max_nodes,
    assign_cluster = FALSE
  )

  node_dist <- stats::dist(gng_out$node_space) %>% as.matrix

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  node_names <- gng_out$nodes %>% mutate(name = as.character(name))
  milestone_network <- gng_out$edges %>%
    left_join(node_names %>% select(i = index, from = name), by = "i") %>%
    left_join(node_names %>% select(j = index, to = name), by = "j") %>%
    mutate(
      length = node_dist[cbind(from, to)],
      directed = FALSE
    ) %>%
    select(from, to, length, directed)

  if (apply_mst) {
    gr <- igraph::graph_from_data_frame(milestone_network, directed = F, vertices = node_names$name)
    milestone_network <- igraph::minimum.spanning.tree(gr, weights = igraph::E(gr)$length) %>% igraph::as_data_frame()
  }

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_dimred_projection_to_wrapper(
    milestone_network = milestone_network,
    dimred_milestones = gng_out$node_space,
    dimred_cells = space
  ) %>% add_timings_to_wrapper(
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
    geom_segment(aes(x = from_Comp1, xend = to_Comp1, y = from_Comp2, yend = to_Comp2), segments, colour = "darkgray") +
    geom_point(aes(Comp1, Comp2, colour = group_id), dimred) +
    geom_point(aes(Comp1, Comp2, colour = milestone_id), dimred_milestones, size = 10, shape = 8) +
    theme(legend.position = "none") +
    labs(x = "Comp1", y = "Comp2")

  if (length(prediction$milestone_ids) <= 9) {
    g <- g + scale_colour_brewer(palette = "Dark2")
  }
  process_dynplot(g, prediction$id)
}
