abstract_paga_description <- function(method) {
  par_set <- makeParamSet(
    makeIntegerParam(id = "n_neighbours", lower = 1, default = 30, upper = 100),
    makeIntegerParam(id = "n_comps", lower = 0, default = 50, upper = 100),
    makeNumericParam(id = "resolution", lower = 0.1, default = 2.5, upper = 10)
  )

  run_fun <- switch(
    method,
    paga = "dynmethods::run_paga",
    agapt = "dynmethods::run_agapt"
  )

  name <- switch(
    method,
    paga = "PAGA",
    agapt = "AGA pseudotime"
  )

  create_ti_method(
    name = name,
    short_name = method,
    package_loaded = c(),
    package_required = c("paga", "igraph"),
    par_set = par_set,
    run_fun = run_fun,
    plot_fun = "dynmethods::plot_paga"
  )
}

#' Inferring trajectories with PAGA
#'
#' @inherit ti_angle description
#'
#' @param n_neighbours Number of neighbours for knn
#' @param n_comps Number of principal components
#' @param resolution Resolution of louvain clustering, which determines the granularity of the clustering. Higher values will result in more clusters.
#'
#' @rdname paga
#' @export
ti_paga <- abstract_paga_description("paga")

# #' Description for agapt
# #' @export
# ti_agapt <- abstract_aga_description("agapt")

run_paga <- function(
  expression,
  grouping_assignment = NULL,
  n_neighbours = 30,
  n_comps = 50,
  resolution = 2.5,
  num_cores = 1
) {
  requireNamespace("paga")
  requireNamespace("igraph")

  set_cores(num_cores)

  # load matcher
  use_virtualenv(file.path(find.package("paga"), "venv"))
  pymatcher <- import("scanpy")

  sc <- import("scanpy.api", as="sc")
  anndata <- import("anndata")

  # preprocess grouping
  if (!is.null(grouping_assignment)) {
    obs <- grouping_assignment %>% mutate(louvain=group_id) %>% select(louvain) %>% r_to_py()
    obs$louvain <- obs$louvain$astype("category")

    adata <- anndata$AnnData(expression, obs)
  } else {
    adata <- anndata$AnnData(expression)
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # run kNN
  sc$pp$neighbors(adata, n_neighbors=as.integer(n_neighbours))

  # run PCA
  sc$tl$pca(adata, n_comps=as.integer(n_comps))

  # run Louvain (if grouping is not supplied as prior information)
  if (is.null(grouping_assignment)) {
    sc$tl$louvain(adata, resolution=resolution)
  }

  # run PAGA
  sc$tl$paga(adata)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # get data
  py$adata <- adata

  grouping <- py_eval("adata.obs.astype({'louvain':str})") %>%  # avoid conversion error
    py_to_r() %>%
    pull(louvain) %>%
    set_names(rownames(expression))
  milestone_ids <- py_eval("adata.obs.louvain.cat.categories.astype('str').tolist()")

  adj_ids <- c("confidence", "confidence_tree", "connectivities")
  adj <- map(adj_ids, function(adj_id) {
    py_eval(glue::glue("adata.uns['paga']['{adj_id}'].todense()")) %>%
      reshape2::melt(varnames=c("from", "to"), value.name=adj_id) %>%
      mutate(from = milestone_ids[as.integer(from)], to = milestone_ids[as.integer(to)])
  }) %>% bind_cols() %>% select(from, to, !!adj_ids)

  # We use the `confidence_tree` to construct the milestone network.
  milestone_network <- adj %>%
    filter(confidence_tree > 0) %>%
    mutate(
      length = 1,
      directed = FALSE
    ) %>%
    select(from, to, length, directed)

  # Wrap the output
  prediction <- wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_grouping(
    group_ids = milestone_ids,
    grouping = grouping
  ) %>% add_cluster_graph(
    milestone_network = milestone_network,
    adj = adj
  ) %>% add_timings(
    tl %>% add_timing_checkpoint("method_afterpostproc")
  )
  prediction
}

# agapt is still in here, in case we still want to add this
#
# run_agapt <- function(
#   expression,
#   start_cells,
#   grouping_assignment = NULL,
#   n_neighbours = 30,
#   n_pcs = 50,
#   n_dcs = 10,
#   resolution = 1,
#   tree_based_confidence = TRUE,
#   verbose = FALSE,
#   num_cores = 1
# ) {
#   requireNamespace("aga")
#   requireNamespace("igraph")
#
#   # TIMING: done with preproc
#   tl <- add_timing_checkpoint(NULL, "method_afterpreproc")
#
#   # sample one start cell
#   start_cell <- sample(start_cells, 1)
#
#   aga_out <- aga::aga(
#     expression = expression,
#     start_cell = start_cell,
#     grouping_assignment = grouping_assignment,
#     n_neighbours = n_neighbours,
#     n_pcs = n_pcs,
#     n_dcs = n_dcs,
#     resolution = resolution,
#     tree_based_confidence = tree_based_confidence,
#     verbose = verbose,
#     num_cores = num_cores
#   )
#   aga_out$obs <- aga_out$obs %>% mutate(group_id = paste0("B", group_id))
#   aga_out$adj <- aga_out$adj %>% mutate(from = paste0("B", from), to = paste0("B", to))
#
#   branch_ids <- unique(c(aga_out$adj$from, aga_out$adj$to, aga_out$obs$group_id))
#
#   # TIMING: done with method
#   tl <- tl %>% add_timing_checkpoint("method_aftermethod")
#
#   # create network between branches
#   branch_network <- aga_out$adj %>%
#     mutate_at(vars(from, to), as.character) %>%
#     filter(aga_adjacency_tree_confidence > 0) %>%
#     select(from, to)
#
#   # determine order of branches, based on location of root cell
#   branch_graph <- igraph::graph_from_data_frame(branch_network)
#   branch_order <-
#     branch_graph %>%
#     igraph::dfs(aga_out$obs %>% filter(cell_id == start_cell) %>% pull(group_id)) %>%
#     .$order %>%
#     names()
#
#   # flip order of branch network if from branch is later than to branch
#   branch_network <- branch_network %>%
#     mutate(from_original = from, to_original = to) %>%
#     mutate(flip = map2(from, to, ~diff(match(c(.x, .y), branch_order)) < 0)) %>%
#     mutate(
#       from = ifelse(flip, to_original, from_original),
#       to = ifelse(flip, from_original, to_original)
#     ) %>%
#     select(from, to)
#
#   # now create milestone network by giving each branch an edge, and adding a zero-length edge between each branch
#   milestone_network <- bind_rows(
#     tibble(
#       from = paste0(branch_ids, "_from"),
#       to = paste0(branch_ids, "_to"),
#       length = 1,
#       directed = TRUE
#     ),
#     mutate(
#       branch_network,
#       from = paste0(from, "_to"),
#       to = paste0(to, "_from"),
#       length = 0,
#       directed = TRUE
#     )
#   )
#
#   # Place each cell along an edge of the milestone network:
#   milestone_ids <- unlist(map(branch_ids, ~ paste0(., c("_from", "_to"))))
#
#   progressions <- aga_out$obs %>%
#     mutate(from = paste0(group_id, "_from"), to = paste0(group_id, "_to")) %>%
#     group_by(group_id) %>%
#     mutate(percentage = dynutils::scale_minmax(aga_pseudotime)) %>%
#     ungroup() %>%
#     select(cell_id, from, to, percentage)
#
#   # Create and return the predicted trajectory
#   wrap_prediction_model(
#     cell_ids = rownames(expression)
#   ) %>% add_trajectory(
#     milestone_ids = milestone_ids,
#     milestone_network = milestone_network,
#     progressions = progressions,
#     divergence_regions = NULL,
#     aga_out = aga_out
#   ) %>% add_timings(
#     tl %>% add_timing_checkpoint("method_afterpostproc")
#   )
# }


plot_paga <- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  milestone_graph <- prediction$adj %>%
    tidygraph::as_tbl_graph() %>%
    tidygraph::activate(edges) %>%
    filter(confidence > 0) %>%
    # remove duplicates, retain the edge with higher weight
    arrange(-confidence_tree, -confidence) %>%
    mutate(edge_id = map2_chr(from, to, ~paste0(sort(c(.x, .y)), collapse="#"))) %>%
    group_by(edge_id) %>%
    filter(row_number() == 1) %>%
    ungroup() %>%
    # add edge type
    mutate(
      edge_type = ifelse(confidence_tree > 0, "tree", ifelse(confidence > 1e-6, "full", "attached"))
    )

  layout <- ggraph::create_layout(
    milestone_graph,
    "fr",
    weights = milestone_graph %>% tidygraph::activate(edges) %>% pull(confidence_tree) %>% {.+1}
  )
  aga_plot <- layout %>%
    ggraph::ggraph() +
    ggraph::geom_edge_link(aes(edge_linetype = edge_type, alpha=edge_type, edge_width=edge_type)) +
    ggraph::geom_edge_link(aes(x=x+(xend-x)/2, y=y+(yend-y)/2, xend = x+(xend-x)/1.999, yend=y+(yend-y)/1.999)) +
    ggraph::scale_edge_linetype_manual(values=c(tree="solid", full="longdash", attached="dashed")) +
    ggraph::scale_edge_alpha_manual(values=c(tree=1, full=0.5, attached=0.5)) +
    ggraph::scale_edge_width_manual(values=c(tree=2, full=1, attached=1)) +
    ggraph::geom_node_point(aes(color = name), size=10) +
    geom_text(aes(x, y, label=name), size=8) +
    theme(legend.position = "none")

  aga_plot %>% process_dynplot(prediction$id)
}
