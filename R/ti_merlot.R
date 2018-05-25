#' Inferring trajectories with MERLoT
#'
#' @inherit ti_angle description
#'
#' @inheritParams ti_dpt
#' @param n_components_to_use Which components to use in downstream analysis
#' @inheritParams merlot::CalculateElasticTree
#' @param FixEndpoints Documentation not provided by authors
#' @param increaseFactor_mu factor by which the mu will be increased for the embedding
#' @param increaseFactor_lambda factor by which the mu will be increased for the embedding
#'
#' @export
#'
#' @include wrapper_create_ti_method.R
ti_merlot <- create_ti_method(
  name = "MERLoT",
  short_name = "merlot",
  package_loaded = c("merlot", "destiny"),
  package_required = c("destiny"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "sigma", default = "local", values = c("local", "global")),
    makeDiscreteParam(id = "distance", default = "euclidean", values = c("euclidean", "cosine", "rankcor")),
    makeIntegerParam(id = "ndim", lower = 2L, upper = 20L, default = 20L),
    makeLogicalParam(id = "density_norm", default = TRUE),
    makeIntegerParam(id = "n_local_lower", lower = 2L, upper = 20L, default = 5L),
    makeIntegerParam(id = "n_local_upper", lower = 2L, upper = 20L, default = 7L),
    makeNumericParam(id = "w_width", lower = -4, upper = 0, default = log(.1), trafo = exp),
    makeIntegerParam(id = "n_components_to_use", lower=2, upper=20, default=3),
    makeIntegerParam(id = "N_yk", lower=2, upper=1000, default=100),
    makeNumericParam(id = "lambda_0", lower = -15, upper = -4, default = -10, trafo = exp),
    makeNumericParam(id = "mu_0", lower = 0.0005, upper = 0.005, default = 0.0025),
    makeNumericParam(id = "increaseFactor_mu", lower=2, upper = 50, default = 20),
    makeNumericParam(id = "increaseFactor_lambda", lower=2, upper = 50, default = 20),
    makeLogicalParam(id = "FixEndpoints", default = F),
    forbidden = quote(n_local_lower > n_local_upper)
  ),
  run_fun = "run_merlot",
  plot_fun = "plot_merlot"
)

run_merlot <- function(
  expression,
  start_cell_ids = NULL,
  n_end_states = NULL,
  sigma,
  distance,
  ndim,
  density_norm,
  n_local_lower,
  n_local_upper,
  w_width,
  n_components_to_use,
  N_yk,
  lambda_0,
  mu_0,
  FixEndpoints,
  increaseFactor_mu,
  increaseFactor_lambda
) {
  requireNamespace("destiny")
  requireNamespace("merlot")

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # create n_local vector
  n_local <- seq(n_local_lower, n_local_upper, by = 1)

  #### Example fromrom inst/examples/ExampleGuo2010.R
  if(!is.null(n_end_states)) {
    n_components_to_use <- n_end_states - 1
  }
  n_components <- max(n_components_to_use, n_components) # always make sure that enough components are extracted, even if the provided n_components is too low

  # Embed Cells into their manifold, in this case we use Diffusion Maps as calculated by Destiny
  DatasetDM <- destiny::DiffusionMap(
    data = expression,
    sigma = sigma,
    distance = distance,
    n_eigs = ndim,
    density_norm = density_norm,
    n_local = n_local,
    verbose = F
  )

  # Extract dimensionality reduction
  CellCoordinates <- DatasetDM@eigenvectors[,seq_len(n_components_to_use)]

  # We calculate the scaffold tree using the first 3 diffusion components from the diffusion map
  ScaffoldTree <- merlot::CalculateScaffoldTree(
    CellCoordinates = CellCoordinates,
    NEndpoints = n_end_states
  )

  # Set the number of nodes to be used to build the Principal Elastic Tree.
  # This is now a parameter of the method

  # We calculate the elastic principal tree using the scaffold tree for its initialization
  ElasticTree <- merlot::CalculateElasticTree(
    ScaffoldTree = ScaffoldTree,
    N_yk = N_yk,
    lambda_0 = lambda_0,
    mu_0 = mu_0,
    FixEndpoints = FixEndpoints
  )

  # Embedd the principal elastic tree into the gene expression space from which it was calculated.
  EmbeddedTree <- merlot::GenesSpaceEmbedding(
    ExpressionMatrix = expression,
    ElasticTree = ElasticTree,
    lambda_0 = lambda_0,
    mu_0 = mu_0,
    increaseFactor_mu = increaseFactor_mu,
    increaseFactor_lambda = increaseFactor_lambda
  )

  # Calculate Pseudotimes for the nodes in the Tree in the full gene expression space.
  # T0=3 means that the Endpoint number 3 in the Endpoints list corresponds to the zygote fate and is used as initial pseudotime t0
  # Any given cell can be used as t0 by specifying its index using the parameter C0=cell_index
  if (is.null(start_cell_ids)) {
    Pseudotimes <- merlot::CalculatePseudotimes(EmbeddedTree, T0=1)
  } else {
    Pseudotimes <- merlot::CalculatePseudotimes(EmbeddedTree, C0=which(rownames(expression) == start_cell_ids[[1]]))
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

  # now calculate milestone network lengths
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
  ) %>% add_trajectory(
    unique(c(milestone_network$from, milestone_network$to)),
    milestone_network,
    NULL,
    progressions = progressions,
    ElasticTree
  ) %>% add_timings(
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
