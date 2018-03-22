#' Description for aga
#' @export
description_aga <- function() abstract_aga_description("aga")

#' Description for agapt
#' @export
description_agapt <- function() abstract_aga_description("agapt")

abstract_aga_description <- function(method) {
  par_set <- makeParamSet(
    makeIntegerParam(id = "n_neighbours", lower = 1, default = 30, upper = 100),
    makeIntegerParam(id = "n_pcs", lower = 0, default = 50, upper = 100),
    makeIntegerParam(id = "n_dcs", lower = 2, default = 10, upper = 50),
    makeNumericParam(id = "resolution", lower = 0.1, default = 1, upper = 10),
    makeLogicalParam(id = "tree_based_confidence", default = TRUE)
  )

  run_fun <- run_aga
  if (method == "agapt") {
    formals(run_fun)$start_cells <- formals(run_fun)$expression
  }
  name <- switch(
    method,
    aga = "AGA",
    agapt = "AGA pseudotime"
  )

  create_description(
    name = name,
    short_name = method,
    package_loaded = c(),
    package_required = c("aga", "igraph"),
    par_set = par_set,
    properties = c(),
    run_fun = run_fun,
    plot_fun = plot_aga
  )
}

run_aga <- function(
  expression,
  start_cells = NULL,
  grouping_assignment = NULL,
  n_neighbours = 30,
  n_pcs = 50,
  n_dcs = 10,
  resolution = 1,
  tree_based_confidence = TRUE,
  verbose = FALSE,
  num_cores = 1
) {
  requireNamespace("aga")
  requireNamespace("igraph")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # sample one start cell, if any are given
  if (!is.null(start_cells)) {
    start_cell <- sample(start_cells, 1)
  } else {
    start_cell <- NULL
  }

  aga_out <- aga::aga(
    counts = expression,
    start_cell = start_cell,
    grouping_assignment = grouping_assignment,
    n_neighbours = n_neighbours,
    n_pcs = n_pcs,
    n_dcs = n_dcs,
    resolution = resolution,
    tree_based_confidence = tree_based_confidence,
    verbose = verbose,
    num_cores = num_cores
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # Reformat everything into milestone_network and milestone_percentages
  cell_ids <- rownames(expression)

  if(is.null(start_cells)) {
    milestone_percentages <- tibble(
      cell_id = aga_out$obs$cell_id,
      milestone_id = aga_out$obs$group_id,
      percentage = 1
    )
    milestone_network <- aga_out$adj %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(aga_adjacency_tree_confidence > 0) %>%
      mutate(
        length = 1,
        directed = TRUE
      ) %>%
      select(from, to, length, directed)

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, milestone_percentages$milestone_id))
  } else {
    # create network between branches
    branch_network <- aga_out$adj %>%
      mutate_at(vars(from, to), as.character) %>%
      filter(aga_adjacency_tree_confidence > 0) %>%
      select(from, to)

    # determine order of branches, based on location of root cell
    branch_graph <- igraph::graph_from_data_frame(branch_network)
    branch_order <- igraph::dfs(
      branch_graph,
      aga_out$obs %>%
        filter(cell_id == start_cells[[1]]) %>%
        pull(group_id)
    )$order %>%
      names()

    # now flip order of branch network if from branch is later than to branch
    branch_network <- branch_network %>%
      mutate(from_original = from, to_original = to) %>%
      mutate(flip = map2(from, to, ~diff(match(c(.x, .y), branch_order)) < 0)) %>%
      mutate(
        from = ifelse(flip, to_original, from_original),
        to = ifelse(flip, from_original, to_original)
      ) %>%
      select(from, to)

    # now create milestone network by giving each branch an edge, and adding a zero-length edge between each branch
    branch_ids <- unique(c(branch_network$from, branch_network$to, aga_out$obs$louvain_groups))

    milestone_network <- bind_rows(
      tibble(
        from = paste0(branch_ids, "_from"),
        to = paste0(branch_ids, "_to"),
        length = 1,
        directed = TRUE
      ),
      branch_network %>% mutate(from = paste0(from, "_to"), to = paste0(to, "_from"), length=0, directed=TRUE)
    )

    progressions <- aga_out$obs %>%
      mutate(from = paste0(group_id, "_from"), to = paste0(group_id, "_to")) %>%
      group_by(group_id) %>%
      mutate(percentage = dynutils::scale_minmax(aga_pseudotime)) %>%
      ungroup()

    progressions <- progressions %>%
      select(cell_id, from, to, percentage)

    milestone_ids <- unique(c(milestone_network$from, milestone_network$to, progressions$from, progressions$to))

    milestone_percentages <- dynwrap::convert_progressions_to_milestone_percentages(
      cell_ids,
      milestone_ids,
      milestone_network,
      progressions
    )
  }

  divergence_regions <- milestone_network %>%
    group_by(from) %>%
    filter(n() > 1) %>%
    mutate(divergence_id = from) %>%
    gather(fromto, milestone_id, from, to) %>%
    mutate(is_start = fromto == "from") %>%
    select(divergence_id, milestone_id, is_start) %>%
    distinct()

  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_trajectory_to_wrapper(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    milestone_percentages = milestone_percentages,
    divergence_regions = divergence_regions,
    aga_out = aga_out
  ) %>% add_timings_to_wrapper(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}


plot_aga <- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  milestone_graph <- prediction$aga_out$adj %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(edges) %>%
    filter(aga_adjacency_full_attachedness > 0) %>%
    # remove duplicates, retain the edge with higher weight
    arrange(-aga_adjacency_tree_confidence, -aga_adjacency_full_confidence) %>%
    mutate(edge_id = map2_chr(from, to, ~paste0(sort(c(.x, .y)), collapse="#"))) %>%
    group_by(edge_id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # add edge type
    mutate(
      edge_type = ifelse(aga_adjacency_tree_confidence > 0, "tree", ifelse(aga_adjacency_full_confidence > 0, "full", "attached"))
    )
  layout <- ggraph::create_layout(
    milestone_graph,
    "fr",
    weights = milestone_graph %>% tidygraph::activate(edges) %>% pull(aga_adjacency_tree_confidence) %>% {.+1}
  )
  layout %>%
    ggraph::ggraph() +
    ggraph::geom_edge_link(aes(edge_linetype = edge_type, alpha=edge_type)) +
    ggraph::geom_edge_link(aes(x=x+(xend-x)/2, y=y+(yend-y)/2, xend = x+(xend-x)/1.999, yend=y+(yend-y)/1.999), arrow=arrow(length=unit(0.1, "inches"))) +
    ggraph::scale_edge_linetype_manual(values=c(tree="solid", full="longdash", attached="dashed")) +
    ggraph::scale_edge_alpha_manual(values=c(tree=1, full=1, attached=0.5)) +
    ggraph::geom_node_point(aes(color = name), size=10) +
    geom_text(aes(x, y, label=name), size=8) +
    ggraph::theme_graph() +
    theme(legend.position = "none")
}
