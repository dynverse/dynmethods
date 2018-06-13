#' Inferring trajectories with SLICE
#'
#' @inherit ti_angle description
#'
#' @inheritParams SLICE::getEntropy
#' @inheritParams SLICE::getLineageModel
#'
#' @export
ti_slice <- create_ti_method(
  name = "SLICE",
  short_name = "slice",
  package_loaded = c(),
  package_required = c("SLICE", "igraph"),
  parameters = list(
    lm.method = list(
      type = "discrete",
      default = "clustering",
      values = c("clustering", "graph"),
      description = "Select \"clustering\" based or \"graph\" based method to infer lineage model"),
    model.type = list(
      type = "discrete",
      default = "tree",
      values = c("tree", "graph"),
      description = "The type of models that will be infered: \"tree\" - directed minimum spanning tree based, \"graph\" - directed graph based"),
    ss.method = list(
      type = "discrete",
      default = "all",
      values = c("all", "top", "pcst"
      ),
      description = "The method for defining core cell set for stable state detection: \nall - all the cells in a cluster constitute the core cell set; \ntop - cells with scEntropy lower than the ss.threshold quantile of all the values in a cluster constitute the core cell set; \npcst - cells with scEntropy lower than the ss.threshold quantile of all the values in a cluster constitute the prize nodes, linear prize-collecting steiner tree algorithm is used to approximate an optimal subnetwork, the cells in the subnetwork constitute the core cell set. Stable states are defined as the centroids of the core cell sets."),

    ss.threshold = list(
      type = "numeric",
      default = 0.25,
      upper = 1,
      lower = 0,
      description = "The threshold used when ss.method is \"top\" or \"pcst\". Default: 0.25."),
    community.method = list(
      type = "discrete",
      default = "louvain",
      values = c("fast_greedy", "edge_betweenness", "label_prop", "leading_eigen", "louvain", "spinglass", "walktrap", "auto"),
      description = "The method for network community detection. \nMost of the community detection methods implemented in the igraph package are supported, \nincluding \"fast_greedy\", \"edge_betweenness\", \"label_prop\", \"leading_eigen\",\"louvain\",\"spinglass\", \"walktrap\". \nIf this parameter is set to \"auto\", the algorithm will perform all the community detection methods and select the one that generates the communities with best modularity. \nOnly take effect when lm.method is \"graph\""),

    cluster.method = list(
      type = "discrete",
      default = "kmeans",
      values = c("kmeans", "pam"),
      description = "Use \"kmeans\" or \"pam\" to divide cells into clusters. Only take effect when lm.method is \"clustering\""),
    k = list(
      type = "discrete",
      default = 0,
      values = c(0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20),
      description = "The number of cell clusters. If NULL, Gap statistic will be used to determine an optimal k."),
    k.max = list(

      type = "integer",
      default = 10L,
      upper = 20L,
      lower = 3L,
      description = "The \"k.max\" parameter of cluster::clusGap(); used when k is NULL."),
    B = list(
      type = "integer",
      default = 100L,
      upper = 500L,
      lower = 3L,
      description = "The \"B\" parameter of cluster::clusGap(); used when k is NULL"),
    k.opt.method = list(
      type = "discrete",
      default = "firstmax",
      values = c("firstmax", "globalmax", "Tibs2001SEmax", "firstSEmax", "globalSEmax"),
      description = "The \"method\" parameter of cluster::maxSE(); used when k is NULL")
  ),
  run_fun = "dynmethods::run_slice",
  plot_fun = "dynmethods::plot_slice"
)

run_slice <- function(
  expression,
  groups_id = NULL,
  features_id = NULL,
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

  # if groups_id is not given, fill it with 1's
  if(!is.null(groups_id)) {
    cellidentity <- groups_id %>%
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
    genes.use = features_id
  )

  # infer entropy-directed cell lineage model
  sc <- SLICE::getLineageModel(
    sc,
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

  # extract the pseudotime
  pseudotime <- map_df(seq_len(nrow(milestone_network)), function(i) {
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
      tibble::rownames_to_column("cell_id") %>%
      mutate(from = from, to = to) %>%
      select(cell_id, from, to, percentage = ptime)
  })

  # check whether the state of a cell is in the network's
  # from or to, and get the earliest timepoint
  progressions <- sc@model$cells.df %>%
    tibble::rownames_to_column("cell_id") %>%
    slice(match(rownames(expression), cell_id)) %>%
    mutate(state = paste0("slice.ss.", slice.state)) %>%
    select(cell_id, state) %>%
    right_join(pseudotime, by = "cell_id") %>%
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

  dimred <- cells.df %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(slice.realcell == 1) %>%
    tibble::column_to_rownames("cell_id") %>%
    select(x, y) %>%
    as.matrix()

  dimred_milestones <- cells.df %>%
    tibble::rownames_to_column("cell_id") %>%
    filter(slice.realcell == 0) %>%
    tibble::column_to_rownames("cell_id") %>%
    select(x, y) %>%
    as.matrix()



  # return output
  wrap_prediction_model(
    cell_ids = rownames(expression)
  ) %>% add_trajectory(
    milestone_ids = milestone_ids,
    milestone_network = milestone_network,
    progressions = progressions,
    divergence_regions = NULL
  ) %>% add_dimred(
    dimred = dimred,
    dimred_milestones = dimred_milestones,
    dimred_trajectory_segments = edge.df[,c("src.x", "src.y", "dst.x", "dst.y")] %>%
      mutate_all(as.numeric) %>%
      as.matrix %>%
      magrittr::set_colnames(c("from_comp_1", "from_comp_2", "to_comp_1", "to_comp_2")),
    cells.df = cells.df, # todo: remove this and modify plot function
    edge.df = edge.df
  ) %>% add_timings(
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
