#' Description for monocle DDRTree
#' @export
description_monocle2_ddrtree <- function() abstract_monocle_description("DDRTree")

#' Description for monocle ICA
#' @export
description_monocle1_ica <- function() abstract_monocle_description("ICA")

# These reduction methods are not implemented yet.
#
# #' Description for monocle SimplePPT
# #' @export
# description_monocle2_simpleppt <- function() abstract_monocle_description("SimplePPT")
#
# #' Description for monocle L1-graph
# #' @export
# description_monocle2_l1graph <- function() abstract_monocle_description("L1-graph")
#
# #' Description for monocle SGL-tree
# #' @export
# description_monocle2_sgltree <- function() abstract_monocle_description("SGL-tree")

abstract_monocle_description <- function(reduction_method) {
  par_set <- switch(
    reduction_method,
    DDRTree = makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = reduction_method, default = reduction_method),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
      makeLogicalParam(id = "auto_param_selection", default = TRUE)
    ),
    ICA = makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = reduction_method, default = reduction_method),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none"))
    )
  )

  short_name <- c(
    "DDRTree" = "mnclDDR",
    "ICA" = "mnclICA",
    "tSNE" = "mncltSNE",
    "SimplePPT" = "mnclSPPT",
    "L1-graph" = "mnclL1gr",
    "SGL-tree" = "mnclSGLT"
  )
  create_description(
    name = pritt("monocle with {reduction_method}"),
    short_name = short_name[reduction_method],
    package_loaded = c("monocle"),
    package_required = c("BiocGenerics", "igraph", "Biobase"),
    par_set = par_set,
    properties = c(),
    run_fun = run_monocle,
    plot_fun = plot_monocle
  )
}

run_monocle <- function(counts,
                        n_end_states = NULL,
                        reduction_method,
                        max_components = 2,
                        norm_method = "vstExprs",
                        auto_param_selection = TRUE) {
  requireNamespace("monocle")
  requireNamespace("BiocGenerics")
  requireNamespace("igraph")
  requireNamespace("Biobase")
  requireNamespace("Matrix")

  # just in case
  if (is.factor(norm_method)) norm_method <- as.character(norm_method)

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # load in the new dataset
  pd <- Biobase::AnnotatedDataFrame(data.frame(row.names = rownames(counts)))
  fd <- Biobase::AnnotatedDataFrame(data.frame(row.names = colnames(counts), gene_short_name = colnames(counts)))
  cds <- monocle::newCellDataSet(t(counts), pd, fd)

  # estimate size factors and dispersions
  cds <- BiocGenerics::estimateSizeFactors(cds)
  cds <- BiocGenerics::estimateDispersions(cds)

  # reduce the dimensionality
  cds <- monocle::reduceDimension(
    cds,
    max_components = max_components,
    reduction_method = reduction_method,
    norm_method = norm_method,
    auto_param_selection = auto_param_selection
  )

  # order the cells
  cds <- monocle::orderCells(cds, num_paths = n_end_states)

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract the igraph and which cells are on the trajectory
  gr <-
    if (reduction_method == "DDRTree") {
      cds@auxOrderingData[[reduction_method]]$pr_graph_cell_proj_tree
    } else if (reduction_method == "ICA") {
      cds@auxOrderingData[[reduction_method]]$cell_ordering_tree
    }
  to_keep <- setNames(rep(TRUE, nrow(counts)), rownames(counts))

  # convert to milestone representation
  edges <- igraph::as_data_frame(gr, "edges")
  if ("weight" %in% edges) {
    edges <- edges %>% rename(length = weight)
  } else {
    edges <- edges %>% mutate(length = 1)
  }
  out <- dynutils::simplify_sample_graph(
    edges =  edges %>% mutate(directed = FALSE),
    to_keep = to_keep,
    is_directed = FALSE
  )

  # retrieve data for visualisation
  plot_data <- postprocess_monocle_cds(cds)

  # TIMING: after postproc
  tl <- tl %>% add_timing_checkpoint("method_afterpostproc")

  # wrap output
  wrap_prediction_model(
    cell_ids = rownames(counts),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    plot_data = plot_data,
    reduction_method = reduction_method
  ) %>% attach_timings_attribute(tl)
}

plot_monocle <- function(prediction) {
  requireNamespace("monocle")
  # Based on monocle::plot_cell_trajectory(cds)
  reduction_method <- prediction$reduction_method
  plot_data <- prediction$plot_data

  g <- ggplot() +
    geom_point(aes_string(x="data_dim_1", y="data_dim_2", color = "State"), plot_data$data_df, size=1.5, na.rm = TRUE) +
    geom_segment(aes_string(x="source_prin_graph_dim_1", y="source_prin_graph_dim_2", xend="target_prin_graph_dim_1", yend="target_prin_graph_dim_2"), size=0.75, linetype="solid", na.rm=TRUE, plot_data$edge_df) +
    theme(legend.position = c(.9, .15))

  if (prediction$reduction_method == "DDRTree") {
    branch_point_df <- plot_data$branch_point_df

    g <- g +
      geom_point(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2"),
                 size=5, na.rm=TRUE, branch_point_df) +
      geom_text(aes_string(x="prin_graph_dim_1", y="prin_graph_dim_2", label="branch_point_idx"),
                size=4, color="white", na.rm=TRUE, branch_point_df)
  }

  process_dynplot(g, prediction$id)
}

postprocess_monocle_cds <- function(cds) {
  requireNamespace("igraph")
  requireNamespace("monocle")
  requireNamespace("Biobase")

  # adapted from monocle::plot_cell_trajectory(cds)
  lib_info_with_pseudo <- Biobase::pData(cds)
  sample_state <- Biobase::pData(cds)$State

  if (cds@dim_reduce_type == "ICA"){
    reduced_dim_coords <- monocle::reducedDimS(cds)
  } else if (cds@dim_reduce_type %in% c("simplePPT", "DDRTree") ){
    reduced_dim_coords <- monocle::reducedDimK(cds)
  }
  ica_space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))

  edge_df <- cds %>%
    monocle::minSpanningTree() %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(ica_space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(ica_space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

  data_df <- t(monocle::reducedDimS(cds)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = 1, data_dim_2 = 2) %>%
    rownames_to_column("sample_name") %>%
    mutate(sample_state) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

  out <- lst(
    ica_space_df,
    data_df,
    edge_df
  )

  if (cds@dim_reduce_type == 'DDRTree'){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    out$branch_point_df <- ica_space_df %>%
      slice(match(mst_branch_nodes, sample_name)) %>%
      mutate(branch_point_idx = seq_len(n()))
  }

  out
}
