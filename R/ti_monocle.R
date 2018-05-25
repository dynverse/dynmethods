abstract_monocle_description <- function(short_name) {
  reduction_method <- c(
    "mnclddr" = "DDRTree",
    "mnclica" = "ICA"
  )[short_name] %>% setNames(NULL)

  par_set <- switch(
    short_name,
    mnclddr = makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = reduction_method, default = reduction_method),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none")),
      makeLogicalParam(id = "auto_param_selection", default = TRUE)
    ),
    mnclica = makeParamSet(
      makeDiscreteParam(id = "reduction_method", values = reduction_method, default = reduction_method),
      makeIntegerParam(id = "max_components", lower = 2L, default = 2L, upper = 20L),
      makeDiscreteParam(id = "norm_method", default = "vstExprs", values = c("vstExprs", "log", "none"))
    )
  )

  create_ti_method(
    name = pritt("Monocle {reduction_method}"),
    short_name = short_name,
    package_loaded = c("monocle"),
    package_required = c("BiocGenerics", "igraph", "Biobase"),
    par_set = par_set,
    run_fun = "run_monocle",
    plot_fun = "plot_monocle"
  )
}

#' Inferring trajectories with Monocle
#'
#' @inherit ti_angle description
#'
#' @inheritParams monocle::reduceDimension
#' @inheritParams monocle::orderCells
#'
#' @seealso [monocle::reduceDimension()], [monocle::orderCells()]
#'
#' @rdname monocle
#'
#' @include wrapper_create_ti_method.R
#'
#' @export
ti_mnclddr <- abstract_monocle_description("mnclddr")

#' @rdname monocle
#' @export
ti_mnclica <- abstract_monocle_description("mnclica")


run_monocle <- function(
  counts,
  n_branches = NULL,
  reduction_method,
  max_components = 2,
  norm_method = "vstExprs",
  auto_param_selection = TRUE
) {
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
  cds <- monocle::orderCells(cds, num_paths = n_branches)

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
  cell_graph <- igraph::as_data_frame(gr, "edges") %>% mutate(directed = FALSE)

  if ("weight" %in% colnames(cell_graph)) {
    cell_graph <- cell_graph %>% rename(length = weight)
  } else {
    cell_graph <- cell_graph %>% mutate(length = 1)
  }

  cell_graph <- cell_graph %>% select(from, to, length, directed)

  # retrieve data for visualisation
  plot_data <- postprocess_monocle_cds(cds)

  # wrap output
  wrap_prediction_model(
    cell_ids = rownames(counts)
  ) %>% add_cell_graph(
    cell_graph = cell_graph,
    to_keep = to_keep,
    plot_data = plot_data,
    reduction_method = reduction_method
  ) %>% add_timings(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
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
  space_df <- Matrix::t(reduced_dim_coords) %>%
    as.data.frame() %>%
    select_(prin_graph_dim_1 = 1, prin_graph_dim_2 = 2) %>%
    mutate(sample_name = rownames(.), sample_state = rownames(.))

  edge_df <- cds %>%
    monocle::minSpanningTree() %>%
    igraph::as_data_frame() %>%
    select_(source = "from", target = "to") %>%
    left_join(space_df %>% select_(source="sample_name", source_prin_graph_dim_1="prin_graph_dim_1", source_prin_graph_dim_2="prin_graph_dim_2"), by = "source") %>%
    left_join(space_df %>% select_(target="sample_name", target_prin_graph_dim_1="prin_graph_dim_1", target_prin_graph_dim_2="prin_graph_dim_2"), by = "target")

  data_df <- t(monocle::reducedDimS(cds)) %>%
    as.data.frame() %>%
    select_(data_dim_1 = 1, data_dim_2 = 2) %>%
    rownames_to_column("sample_name") %>%
    mutate(sample_state) %>%
    left_join(lib_info_with_pseudo %>% rownames_to_column("sample_name"), by = "sample_name")

  out <- lst(
    space_df,
    data_df,
    edge_df
  )

  if (cds@dim_reduce_type == "DDRTree"){
    mst_branch_nodes <- cds@auxOrderingData[[cds@dim_reduce_type]]$branch_points
    out$branch_point_df <- space_df %>%
      slice(match(mst_branch_nodes, sample_name)) %>%
      mutate(branch_point_idx = seq_len(n()))
  }

  out
}
