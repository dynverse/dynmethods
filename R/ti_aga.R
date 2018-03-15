#' Description for aga
#' @export
description_aga <- function() create_description(
  name = "AGA",
  short_name = "aga",
  package_loaded = c(),
  package_required = c("aga"),
  par_set = makeParamSet(
    makeNumericParam(id = "n_neighbours", lower = 1, default = 30, upper = 100),
    makeNumericParam(id = "n_pcs", lower = 0, default = 50, upper = 100),
    makeNumericParam(id = "n_dcs", lower = 2, default = 10, upper = 50)

  ),
  properties = c(),
  run_fun = run_aga,
  plot_fun = plot_aga
)

## TODO: handle start cells
run_aga <- function(
  counts,
  grouping_assignment=NULL,
  start_cells=NULL,
  n_neighbours = 30,
  n_pcs = 50,
  n_dcs = 10,
  verbose=F,
  num_cores=1
) {
  requireNamespace("aga")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  aga_args <- as.list(environment()) %>% {.[intersect(methods::formalArgs(aga::aga), names(.))]}
  aga_out <- do.call(aga::aga, aga_args)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  if(is.null(start_cells)) {
    milestone_percentages <- tibble(
      cell_id = aga_out$obs$cell_id,
      milestone_id = aga_out$obs$group_id,
      percentage = 1
    )
    milestone_network <- aga_out$adj %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(aga_adjacency_tree_confidence > 0) %>%
      mutate(length = 1) %>%
      mutate(directed=TRUE) %>%
      select(from, to, length, directed)
    divergence_regions <- milestone_network %>%
      group_by(from) %>%
      filter(n() > 1) %>%
      mutate(divergence_id = from) %>%
      gather(fromto, milestone_id, from, to) %>%
      mutate(is_start = fromto == "from") %>%
      select(divergence_id, milestone_id, is_start) %>%
      distinct()

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, milestone_percentages$milestone_id))
    prediction <- wrap_prediction_model(
      cell_ids = rownames(counts)
    ) %>%
      add_trajectory_to_wrapper(
        milestone_ids = milestone_ids,
        milestone_network = milestone_network,
        milestone_percentages = milestone_percentages,
        divergence_regions = divergence_regions,
        aga_out = aga_out
      )
  } else {
    stop("Not supported yet, have to combine pseudotimes (located in obs dataframe) with network structure. Probably will have to convert the graph to its line graph and put the cells on that by scaling the pseudotime for each branch")
  }

  prediction
}


plot_aga <- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  milestone_graph <- prediction$aga_out$adj %>% tidygraph::as_tbl_graph()
  layout <- ggraph::create_layout()
  milestone_graph %>%
    ggraph::ggraph("tree") +
      ggraph::geom_node_point(aes(color = name), size=5) +
      geom_text(aes(x, y, label=name)) +
      geom_
      ggraph::theme_graph() +
      theme(legend.position = "none")
}
