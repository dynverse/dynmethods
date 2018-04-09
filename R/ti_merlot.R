description_merlot <- function(short_name) {
  create_description(
    name = "MERLot",
    short_name = "merlot",
    package_loaded = c("merlot", "destiny"),
    package_required = c("destiny"),
    par_set = makeParamSet(
      makeLogicalParam("density_norm", TRUE),
      makeIntegerParam("n_components", 2, 20, default=20),
      makeIntegerParam("n_components_to_use", 2, 20, default=3),
      makeIntegerParam("NumberOfNodes", 2, 1000, default=100)
    ),
    properties = c(),
    run_fun = run_merlot,
    plot_fun = plot_merlot
  )
}

run_merlot <- function(
  expression,
  start_cell_ids = NULL,
  density_norm = TRUE,
  n_components = 20,
  n_components_to_use = 3,
  NumberOfNodes = 100

) {
  requireNamespace("destiny")
  requireNamespace("merlot")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # From inst/examples/ExampleGuo2010.R
  # Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
  DatasetDM <- destiny::DiffusionMap(expression, density_norm = density_norm, verbose = F, n_eigs = n_components)

  # The first 3 diffusion map components will be used for this example
  CellCoordinates=DatasetDM@eigenvectors[,1:ncomp]

  # We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
  ScaffoldTree=CalculateScaffoldTree(CellCoordinates = CellCoordinates)

  # Set the number of nodes to be used to build the Principal Elastic Tree.
  # NumberOfNodes=100

  # We calculate the elastic principal tree using the scaffold tree for its initialization
  ElasticTree= CalculateElasticTree(ScaffoldTree = ScaffoldTree, N_yk = NumberOfNodes)
  # plot_elastic_tree(ElasticTree)

  # Embedd the principal elastic tree into the gene expression space from which it was calculated.
  EmbeddedTree= GenesSpaceEmbedding(ExpressionMatrix = expression, ElasticTree = ElasticTree)

  # Calculate Pseudotimes for the nodes in the Tree in the full gene expression space.
  # T0=3 means that the Endpoint number 3 in the Endpoints list corresponds to the zygote fate and is used as initial pseudotime t0
  # Any given cell can be used as t0 by specifying its index using the parameter C0=cell_index
  if (is.null(start_cell_ids)) {
    Pseudotimes=CalculatePseudotimes(EmbeddedTree, T0=1)
  } else {
    Pseudotimes=CalculatePseudotimes(EmbeddedTree, C0=which(rownames(expression) == start_cell_ids[[1]]))
  }

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # first add both the milestone network (without lengths) and progressions
  milestone_network <- ElasticTree$Connectivity %>%
    as.data.frame() %>%
    purrr::set_names(c("from", "to")) %>%
    mutate(edge_id = row_number())
  progressions <- tibble(
    cell_id = rownames(expression),
    edge_id = ElasticTree$Cells2Branches,
    pseudotime = Pseudotimes$Proyected_Times_Cells
  ) %>%
    left_join(milestone_network, "edge_id")

  # now calcualte milestone network lengths
  milestone_network <- left_join(
    milestone_network,
    progressions %>%
      group_by(edge_id) %>%
      summarise(length = max(pseudotime) - min(pseudotime)),
    "edge_id"
  ) %>%
    mutate(length = ifelse(is.na(length), mean(length, na.rm=T), length))

  # now calculate percentages of progression
  progressions <- progressions %>%
    group_by(edge_id) %>%
    mutate(percentage = (pseudotime - min(pseudotime))/(max(pseudotime) - min(pseudotime))) %>%
    select(cell_id, from, to, percentage)

  # wrap output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>%
    add_trajectory_to_wrapper(
      unique(c(milestone_network$from, milestone_network$to)),
      milestone_network,
      NULL,
      progressions = progressions,
      ElasticTree
    ) %>%
    add_timings_to_wrapper(
      timings = tl %>% add_timing_checkpoint("method_afterpostproc")
    )
}

plot_merlot <- function(prediction) {
  requireNamespace("ggraph")
  requireNamespace("tidygraph")

  space_elastic_tree <- prediction$ElasticTree$Nodes %>%
    as.data.frame() %>%
    {colnames(.) <- paste0("Comp", seq_len(ncol(.))); .}
  colnames(space_elastic_tree)[1:2] <- c("x", "y")
  network_elastic_tree <- prediction$ElasticTree$Branches %>% map(function(branches) {
    tibble(from = lead(as.character(branches), 1), to = as.character(branches)) %>%
      drop_na()
  }) %>% bind_rows()
  graph_elastic_tree <- tidygraph::tbl_graph(space_elastic_tree, network_elastic_tree %>% mutate_all(as.numeric))

  space_cells <- prediction$ElasticTree$CellCoords %>%
    as.data.frame() %>%
    {colnames(.) <- paste0("Comp", seq_len(ncol(.))); .} %>%
    mutate(cell_id = prediction$cell_ids) %>%
    mutate(branch_id = factor(prediction$ElasticTree$Cells2Branches))

  plot <- ggraph::ggraph(graph_elastic_tree, layout="manual", node.positions=space_elastic_tree) +
    ggraph::geom_node_point() +
    ggraph::geom_edge_fan() +
    geom_point(aes(Comp1, Comp2, color=branch_id), space_cells)
  process_dynplot(plot, "")
}
