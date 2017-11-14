#' Description for SLICER
#' @export
description_slicer <- function() create_description(
  name = "SLICER",
  short_name = "SLICER",
  package_loaded = c(),
  package_required = c("SLICER", "lle", "igraph"),
  par_set = makeParamSet(
    makeIntegerParam(id = "kmin", lower = 2L, upper = 20L, default = 10L),
    makeIntegerParam(id = "m", lower = 2L, upper = 20L, default = 2L),
    makeNumericParam(id = "min_branch_len", lower = 0.5, upper = 20, default = 5),
    makeNumericParam(id = "min_representative_percentage", lower = 0.5, upper = 1, default = 0.8),
    makeNumericParam(id = "max_same_milestone_distance", lower = 0.1, upper = 10, default = 0.1)
  ),
  properties = c(),
  run_fun = run_slicer,
  plot_fun = plot_slicer
)

run_slicer <- function(counts,
                       start_cells,
                       end_cells = NULL,
                       kmin = 10,
                       m = 2,
                       min_branch_len = 5,
                       min_representative_percentage = 0.8,
                       max_same_milestone_distance = 0.1,
                       verbose = FALSE) {
  requireNamespace("SLICER")
  requireNamespace("lle")
  requireNamespace("igraph")

  start_cell <- sample(start_cells, 1)

  # log transform expresison
  expr <- log2(counts + 1)

  # use 'neighbourhood variance' to identify genes that vary smoothly
  genes <- SLICER::select_genes(expr)
  expr_filt <- expr[, genes]

  # stop output if not verbose
  if (!verbose) {
    sink("/dev/null")
  }

  # determine k for knn
  k <- SLICER::select_k(expr_filt, kmin = kmin)

  # perform local linear embedding
  traj_lle <- lle::lle(expr_filt, m = m, k = k)$Y
  rownames(traj_lle) <- rownames(expr_filt)
  colnames(traj_lle) <- paste0("Comp", seq_len(ncol(traj_lle)))

  # resume output if not verbose
  if (!verbose) {
    sink()
  }

  # get LLE KNN graph
  traj_graph <- SLICER::conn_knn_graph(traj_lle, k = k)

  # find extreme cells
  if (is.null(end_cells)) {
    ends <- SLICER::find_extreme_cells(traj_graph, traj_lle, do_plot = FALSE)
  } else {
    ends <- match(c(start_cell, end_cells), rownames(counts))
  }

  # order cells
  start <- which(rownames(expr_filt) == start_cell)
  cells_ordered <- SLICER::cell_order(traj_graph, start)

  # get shortest paths to start and all other nodes
  shortest_paths <- igraph::shortest_paths(traj_graph, start)
  edges <- lapply(shortest_paths$vpath, function(path) {
    P <- rbind(path[-length(path)], path[-1]) %>% as.vector
    igraph::E(traj_graph, P = P)
  })
  subgr <- igraph::subgraph.edges(traj_graph, eids = unique(unlist(edges)))

  # prepare sample graph simplification
  simp_edges <- igraph::as_data_frame(subgr, "edges") %>%
    select(from, to, length = weight) %>%
    mutate(
      from = rownames(expr_filt)[from],
      to = rownames(expr_filt)[to],
      directed = TRUE
    )
  sh_p_to_ends <- igraph::shortest_paths(subgr, start, ends)
  nodes_to_keep <- unique(sh_p_to_ends$vpath %>% unlist)
  to_keep <- setNames(igraph::V(traj_graph) %in% nodes_to_keep, rownames(expr_filt))

  # simplify graph by projecting cells to the closest point on a trajectory between
  # the start node and one of the end nodes
  out <- simplify_sample_graph(simp_edges, to_keep, is_directed = FALSE)

  # return output
  wrap_ti_prediction(
    ti_type = "multifurcating",
    id = "SLICER",
    cell_ids = rownames(expr_filt),
    milestone_ids = out$milestone_ids,
    milestone_network = out$milestone_network,
    progressions = out$progressions,
    dimred_samples = traj_lle %>% as.data.frame() %>% rownames_to_column("cell_id"),
    traj_graph = subgr,
    start = start,
    ends = ends,
    to_keep = to_keep
  )
}

#' @importFrom grDevices colorRampPalette
plot_slicer <- function(prediction) {
  requireNamespace("SLICER")
  requireNamespace("igraph")

  # based on SLICER::graph_process_distance(traj_graph, dimred_samples[,c("Comp1", "Comp2")], start)

  # calculate the geodesic distances between samples
  dimred_samples <- prediction$dimred_samples
  geodesic_dists <- SLICER::process_distance(prediction$traj_graph, prediction$start)[1,] %>% dynutils::scale_minmax()

  # get colour scale
  plotclr <- grDevices::colorRampPalette(c("black", "red", "yellow"), space="rgb")(50)

  # construct cell df
  cell_df <- dimred_samples %>% mutate(dist = geodesic_dists)

  # construct edge_df
  edge_df <- prediction$traj_graph %>%
    igraph::as_data_frame("edges") %>%
    do(with(., data.frame(row.names = NULL, from = dimred_samples[from,], to = dimred_samples[to,], stringsAsFactors = F))) %>%
    mutate(edge_kept = prediction$to_keep[from.cell_id] & prediction$to_keep[to.cell_id])

  # make plot
  aes_segm <- aes(x = from.Comp1, xend = to.Comp1, y = from.Comp2, yend = to.Comp2)
  g <- ggplot() +
    geom_segment(aes_segm, edge_df %>% filter(!edge_kept), colour = "gray", size = .35) +
    geom_segment(aes_segm, edge_df %>% filter(edge_kept), size = .5) +
    geom_point(aes(Comp1, Comp2, colour = dist), cell_df) +
    scale_colour_gradientn(colours = plotclr) +
    theme(legend.position = c(.92, .12))

  process_dynplot(g, prediction$id)
}
