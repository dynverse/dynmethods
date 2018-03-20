#' Description for SLICE
#' @export
description_slice <- function() create_description(
  name = "SLICE",
  short_name = "slice",
  package_loaded = c(),
  package_required = c("SLICE", "igraph"),
  par_set = makeParamSet(
    makeDiscreteParam(id = "lm.method", default = "clustering", values = c("clustering", "graph")),
    makeDiscreteParam(id = "model.type", default = "tree", values = c("tree", "graph")),
    makeDiscreteParam(id = "ss.method", default = "all", values = c("all", "top", "pcst")),
    makeNumericParam(id = "ss.threshold", default=0.25, lower=0, upper=1),
    makeDiscreteParam(id = "community.method", default = "louvain", values = c("fast_greedy", "edge_betweenness", "label_prop", "leading_eigen", "louvain", "spinglass", "walktrap", "auto")),
    makeDiscreteParam(id = "cluster.method", default = "kmeans", values = c("kmeans", "pam")),
    makeDiscreteParam(id = "k", default = 0, values = c(0, 3:20)),
    makeIntegerParam(id = "k.max", lower = 3L, upper = 20L, default = 10L),
    makeIntegerParam(id = "B", lower = 3L, upper = 500L, default = 100L),
    makeDiscreteParam(id = "k.opt.method", default = "firstmax", values = c("firstmax", "globalmax", "Tibs2001SEmax", "firstSEmax", "globalSEmax"))
  ),
  properties = c(),
  run_fun = run_slice,
  plot_fun = plot_slice
)

run_slice <- function(
  expression,
  grouping_assignment = NULL,
  marker_feature_ids = NULL,
  lm.method = "clustering",
  model.type = "tree",
  ss.method = "all",
  ss.threshold = 0.25,
  community.method = "louvain",
  cluster.method = "kmeans",
  k = 0,
  k.max = 10,
  B = 100,
  k.opt.method = "firstmax"
) {
  requireNamespace("SLICE")
  requireNamespace("igraph")

  # if k is 0, set to NULL
  if (k == 0) {
    k <- NULL
  }

  # if grouping_assignment is not given, fill it with 1's
  if(!is.null(grouping_assignment)) {
    cellidentity <- grouping_assignment %>%
      slice(match(rownames(expression), cell_id)) %>%
      pull(group_id) %>%
      factor()
  } else {
    cellidentity <- factor(rep(1, nrow(expression)))
  }

  # TIMING: done with preproc
  tl <- add_timing_checkpoint(NULL, "method_afterpreproc")

  # wrap data
  sc <- SLICE::construct(
    exprmatrix = as.data.frame(t(expression)),
    cellidentity = cellidentity
  )

  num_genes <- ncol(expression)
  km <- matrix(
    runif(num_genes * num_genes),
    ncol = num_genes,
    dimnames = list(colnames(expression), colnames(expression))
  )

  # calculate the entropy of individual cells
  sc <- SLICE::getEntropy(sc, km = km)

  # reduce expression space
  sc <- SLICE::getRDS(
    sc,
    method = "pca",
    num_dim = 2,
    log.base = 2,
    do.center = TRUE,
    do.scale = FALSE,
    use.cor = TRUE,
    min.var = 0,
    min.cells = 0,
    genes.use = marker_feature_ids
  )

  # infer entropy-directed cell lineage model
  sc <- SLICE::getLineageModel(
    sc,
    lm.method = lm.method,
    model.type = model.type,
    ss.method = ss.method,
    ss.threshold = ss.threshold,
    community.method = community.method,
    cluster.method = cluster.method,
    k = k,
    k.max = k.max,
    B = B,
    k.opt.method = k.opt.method,
    do.plot = FALSE
  )

  # TIMING: done with method
  tl <- tl %>% add_timing_checkpoint("method_aftermethod")

  # extract the milestone ids
  lin_model <- sc@model$lineageModel
  milestone_ids <- names(igraph::V(lin_model))

  # extract the milestone network
  milestone_network <- lin_model %>%
    igraph::as_data_frame() %>%
    rename(length = weight) %>%
    mutate(directed = TRUE)

  # extract the pseudotimes
  pseudotimes <- map_df(seq_len(nrow(milestone_network)), function(i) {
    from <- milestone_network[i, 1]
    to <- milestone_network[i, 2]
    sc_tmp <- SLICE::getTrajectories(
      sc,
      method = "pc",
      start = match(from, milestone_ids),
      end = match(to, milestone_ids),
      do.plot = FALSE,
      do.trim = FALSE
    )
    sc_tmp@transitions[[1]]$i.pseudotime %>%
      rownames_to_column("cell_id") %>%
      mutate(from = from, to = to) %>%
      select(cell_id, from, to, percentage = ptime)
  })

  # check whether the state of a cell is in the network's
  # from or to, and get the earliest timepoint
  progressions <- sc@model$cells.df %>%
    rownames_to_column("cell_id") %>%
    slice(match(rownames(expression), cell_id)) %>%
    mutate(state = paste0("slice.ss.", slice.state)) %>%
    select(cell_id, state) %>%
    right_join(pseudotimes, by = "cell_id") %>%
    filter((state == from) | (state == to)) %>%
    group_by(cell_id) %>%
    arrange(percentage) %>%
    slice(1) %>%
    select(-state) %>%
    ungroup()

  # collect data for visualisation
  cells.df <- sc@model$cells.df
  edge.df <- igraph::get.edgelist(lin_model) %>%
    as.data.frame() %>%
    mutate(
      ix = match(V1, rownames(cells.df)),
      iy = match(V2, rownames(cells.df)),
      src.x = cells.df$x[ix],
      src.y = cells.df$y[ix],
      dst.x = cells.df$x[iy],
      dst.y = cells.df$y[iy]
    )

  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_trajectory_to_wrapper(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_dimred_to_wrapper(
    dimred = cells.df %>% subset(slice.realcell == 1) %>% as.matrix %>% .[,c("x", "y")],
    dimred_milestones = cells.df %>% subset(slice.realcell == 0) %>% as.matrix %>% .[,c("x", "y")],
    dimred_trajectory_segments = edge.df[,c("src.x", "src.y", "dst.x", "dst.y")] %>% as.matrix %>% set_colnames(c("from_x", "from_y", "to_x", "to_y")),
    cells.df = cells.df, # todo: remove this and modify plot function
    edge.df = edge.df
  ) %>% add_timings_to_wrapper(
    timings = tl %>% add_timing_checkpoint("method_afterpostproc")
  )
}

#' @importFrom grid arrow
plot_slice <- function(prediction) {
  requireNamespace("igraph")

  # Adapted from SLICE code
  cells.df <- prediction$cells.df
  edge.df <- prediction$edge.df

  # split up the "cells"
  milestones <- cells.df %>% subset(slice.realcell == 0)
  cells_stable <- cells.df %>% subset(slice.realcell==1 & slice.stablestate != "NA")
  cells_unstable <- cells.df %>% subset(slice.realcell==1 & slice.stablestate == "NA")

  # make plot
  g <- ggplot(mapping = aes(x, y, size = entropy)) +
    geom_point(aes(col = slice.state), cells_stable) +
    geom_point(aes(col = slice.state), cells_unstable) +
    geom_segment(aes(x = src.x, y = src.y, xend = dst.x, yend = dst.y), edge.df,
                 size = 2, linetype = "solid", col = "black",
                 alpha = 0.6, arrow = grid::arrow(), na.rm = TRUE) +
    geom_point(data = milestones, col = "black") +
    theme(legend.position = c(.92, .2)) +
    labs(colour = "Group", size = "Entropy")
  process_dynplot(g, prediction$id)
}
